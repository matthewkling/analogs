// src/fft_grid.cpp
//
// Binned-convolution (FFT) engine for kernel-weighted aggregation.
//
// Every aggregate statistic the exact engine computes for a focal f is built
// from sums of the form
//
//     S(f) = sum_j  m_j * K(z_f - z_j)
//
// where z are the active axes (environment and/or geography), m_j is a
// per-pool-point "mass" (e.g. pw_j, pw_j * y_j, pw_j * x_j * y_j), and K is a
// product of per-family kernels. One such sum is a "field". This engine
// evaluates each field for all focals at once:
//
//   1. Bin:       linear binning of each pool point's mass onto the 2^D nodes
//                 of a regular grid surrounding it.
//   2. Convolve:  for each axis group (one per kernel family), convolve every
//                 slice of the grid along that group's axes with the family's
//                 kernel, via real FFTs on a zero-padded buffer. Because the
//                 kernel is a product across families, convolving group by
//                 group equals convolving with the full product kernel.
//   3. Read off:  multilinear interpolation of the convolved grid at each
//                 focal (and, for leave-one-out, subtraction of the focal's
//                 own binned contribution).
//
// fft_finalize_cpp() then turns per-focal field values into the exact
// engine's output columns, reusing focal_stats.h and ridge_solve.h.
//
// Grid geometry, kernel arrays, padded FFT sizes, and the field definitions
// are all decided on the R side (R/fft.R); this file only does arithmetic.

#define POCKETFFT_NO_MULTITHREADING
#include "pocketfft_hdronly.h"

#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <functional>
#include "codes.h"
#include "se_code.h"
#include "focal_stats.h"
#include "ridge_solve.h"

using namespace Rcpp;

namespace analogs {
namespace fftgrid {

typedef std::complex<double> cplx;

// Regular grid over D axes. Row-major node layout: axis D-1 varies fastest.
struct Grid {
      int D = 0;
      std::vector<double> origin, delta;
      std::vector<int> n;          // nodes per axis
      std::vector<size_t> stride;  // row-major strides (in nodes)
      size_t N = 1;                // total nodes

      Grid(const NumericVector& origin_in, const NumericVector& delta_in,
           const IntegerVector& n_in) {
            D = n_in.size();
            origin.assign(origin_in.begin(), origin_in.end());
            delta.assign(delta_in.begin(), delta_in.end());
            n.assign(n_in.begin(), n_in.end());
            stride.assign(D, 1);
            for (int k = D - 2; k >= 0; --k) stride[k] = stride[k + 1] * n[k + 1];
            N = 1;
            for (int k = 0; k < D; ++k) N *= static_cast<size_t>(n[k]);
      }
};

// The 2^D grid nodes surrounding a point, with multilinear weights (summing
// to 1). Corners with zero weight are dropped, so a point exactly on a node
// (or on the last node of an axis) never indexes past the grid. Per-axis
// node indices are kept for the leave-one-out self-term.
struct Corners {
      int count = 0;
      size_t idx[32];    // linear node index
      double w[32];      // multilinear weight
      int node[32][5];   // per-axis node index (D <= 5)
};

// Point coordinates are read from row `i` of an R matrix with `nrow` rows.
// Returns false (with no corners) for a point outside the grid: the grid is
// trimmed to where pool mass can reach focals, so such a pool point
// contributes nothing and such a focal has no analogs.
inline bool point_corners(const Grid& g, const double* mat, size_t nrow,
                          size_t i, Corners& c) {
      int base[5];
      double frac[5];
      c.count = 0;
      for (int k = 0; k < g.D; ++k) {
            double u = (mat[i + k * nrow] - g.origin[k]) / g.delta[k];
            const double umax = static_cast<double>(g.n[k] - 1);
            if (!(u >= -1e-9 && u <= umax + 1e-9)) return false;  // outside (or NaN)
            if (u < 0.0) u = 0.0;
            if (u > umax) u = umax;
            int b = static_cast<int>(std::floor(u));
            if (b >= g.n[k] - 1) b = g.n[k] - 1;
            base[k] = b;
            frac[k] = u - b;
      }
      c.count = 0;
      const int n_corner = 1 << g.D;
      for (int bits = 0; bits < n_corner; ++bits) {
            double w = 1.0;
            size_t lin = 0;
            for (int k = 0; k < g.D; ++k) {
                  const bool hi = (bits >> k) & 1;
                  const double wk = hi ? frac[k] : (1.0 - frac[k]);
                  if (wk == 0.0) { w = 0.0; break; }
                  w *= wk;
                  const int node = base[k] + (hi ? 1 : 0);
                  c.node[c.count][k] = node;
                  lin += static_cast<size_t>(node) * g.stride[k];
            }
            if (w == 0.0) continue;
            c.idx[c.count] = lin;
            c.w[c.count] = w;
            ++c.count;
      }
      return true;
}

// One axis group (a kernel family's contiguous block of grid axes).
struct Group {
      int first = 0;                 // first grid axis in the group
      int d = 0;                     // number of axes in the group
      std::vector<int> n;            // unpadded nodes per group axis
      std::vector<int> m;            // padded FFT size per group axis
      std::vector<int> h;            // kernel half-width (nodes) per group axis
      size_t n_slice_nodes = 1;      // prod(n)
      size_t m_real = 1;             // prod(m)
      size_t m_cplx = 1;             // prod(m) with last axis m/2+1
      size_t n_slices = 1;           // grid nodes / n_slice_nodes
      size_t outer_stride = 0;       // grid stride of the axis before the group
      size_t inner_size = 1;         // product of n over axes after the group
      std::vector<size_t> pad_index; // slice node (row-major over n) -> padded index
      pocketfft::shape_t shape;
      pocketfft::stride_t stride_r, stride_c;
      pocketfft::shape_t axes;

      // Grid node index of slice s's first node; slice node t (row-major
      // over n) is at slice_base(s) + t * inner_size.
      inline size_t slice_base(size_t s) const {
            const size_t outer = s / inner_size;
            const size_t inner = s % inner_size;
            return outer * n_slice_nodes * inner_size + inner;
      }
};

// A kernel's per-group arrays (from R) and their spectra.
struct KernelSpec {
      std::vector<std::vector<double>> arr;   // per group: column-major, dims 2h+1
      std::vector<std::vector<cplx>> spec;    // per group: r2c of padded kernel
      double kmax = 0.0;                      // max |K| of the full product kernel
};

// Value of a group's kernel array at an integer offset (0 outside support).
inline double kernel_at(const Group& gr, const std::vector<double>& arr,
                        const int* off) {
      size_t lin = 0, mult = 1;
      for (int a = 0; a < gr.d; ++a) {
            const int o = off[a];
            if (o < -gr.h[a] || o > gr.h[a]) return 0.0;
            lin += static_cast<size_t>(o + gr.h[a]) * mult;
            mult *= static_cast<size_t>(2 * gr.h[a] + 1);
      }
      return arr[lin];
}

// Convolve one slice in place: copy into the zero-padded buffer, multiply
// spectra, transform back, copy the unpadded region out. Empty slices (common:
// most climate nodes hold no pool mass) are left as zeros without any FFT.
inline void convolve_slice(const Group& gr, const std::vector<cplx>& kspec,
                           double* grid, size_t s,
                           std::vector<double>& pad, std::vector<cplx>& freq) {
      const size_t base = gr.slice_base(s);
      const size_t step = gr.inner_size;
      bool empty = true;
      for (size_t t = 0; t < gr.n_slice_nodes; ++t) {
            if (grid[base + t * step] != 0.0) { empty = false; break; }
      }
      if (empty) return;

      std::fill(pad.begin(), pad.end(), 0.0);
      for (size_t t = 0; t < gr.n_slice_nodes; ++t) {
            pad[gr.pad_index[t]] = grid[base + t * step];
      }
      pocketfft::r2c(gr.shape, gr.stride_r, gr.stride_c, gr.axes, pocketfft::FORWARD,
                     pad.data(), freq.data(), 1.0);
      for (size_t q = 0; q < gr.m_cplx; ++q) freq[q] *= kspec[q];
      pocketfft::c2r(gr.shape, gr.stride_c, gr.stride_r, gr.axes, pocketfft::BACKWARD,
                     freq.data(), pad.data(), 1.0 / static_cast<double>(gr.m_real));
      for (size_t t = 0; t < gr.n_slice_nodes; ++t) {
            grid[base + t * step] = pad[gr.pad_index[t]];
      }
}

// Runs f(begin, end) over parallel ranges.
template <class F>
struct RangeWorker : public RcppParallel::Worker {
      F& f;
      explicit RangeWorker(F& f_) : f(f_) {}
      void operator()(std::size_t begin, std::size_t end) { f(begin, end); }
};

struct SliceWorker : public RcppParallel::Worker {
      const Group& gr;
      const std::vector<cplx>& kspec;
      double* grid;
      SliceWorker(const Group& gr_, const std::vector<cplx>& kspec_, double* grid_)
            : gr(gr_), kspec(kspec_), grid(grid_) {}
      void operator()(std::size_t begin, std::size_t end) {
            std::vector<double> pad(gr.m_real);
            std::vector<cplx> freq(gr.m_cplx);
            for (std::size_t s = begin; s < end; ++s) {
                  convolve_slice(gr, kspec, grid, s, pad, freq);
            }
      }
};

// Shared, read-only state for evaluating fields.
struct Engine {
      Grid grid;
      std::vector<Group> groups;
      std::vector<KernelSpec> kernels;
      const double* pool;  size_t n_pool;
      const double* focal; size_t n_focal;
      const double* mass;  // n_pool x F, column-major
      std::vector<int> field_kernel;
      bool exclude_self;
      std::vector<double> self_term;  // n_pool x n_kernels (column-major), if exclude_self
      double* out;                    // n_focal x F, column-major

      Engine(const NumericVector& origin, const NumericVector& delta,
             const IntegerVector& n_nodes)
            : grid(origin, delta, n_nodes) {}

      // Bin, convolve and read off field f, using `buf` (length grid.N) as the
      // grid. With `parallel_inner`, slices are convolved in parallel.
      void run_field(int f, std::vector<double>& buf, bool parallel_inner) const {
            std::fill(buf.begin(), buf.end(), 0.0);
            const double* mf = mass + static_cast<size_t>(f) * n_pool;
            const KernelSpec& ks = kernels[field_kernel[f]];

            // 1. Bin
            Corners c;
            double total_abs = 0.0;
            for (size_t j = 0; j < n_pool; ++j) {
                  const double mj = mf[j];
                  if (mj == 0.0) continue;
                  total_abs += std::fabs(mj);
                  point_corners(grid, pool, n_pool, j, c);
                  for (int q = 0; q < c.count; ++q) buf[c.idx[q]] += mj * c.w[q];
            }
            double* of = out + static_cast<size_t>(f) * n_focal;
            if (total_abs == 0.0) {
                  std::fill(of, of + n_focal, 0.0);
                  return;
            }

            // 2. Convolve, group by group
            for (size_t gi = 0; gi < groups.size(); ++gi) {
                  const Group& gr = groups[gi];
                  if (parallel_inner && gr.n_slices > 1) {
                        SliceWorker w(gr, ks.spec[gi], buf.data());
                        RcppParallel::parallelFor(0, gr.n_slices, w);
                  } else {
                        std::vector<double> pad(gr.m_real);
                        std::vector<cplx> freq(gr.m_cplx);
                        for (size_t s = 0; s < gr.n_slices; ++s) {
                              convolve_slice(gr, ks.spec[gi], buf.data(), s, pad, freq);
                        }
                  }
            }

            // 3. Read off at focals; values below FFT round-off are zeroed so
            // an empty neighborhood reads as exactly empty.
            const double tol = 1e-11 * total_abs * ks.kmax;
            const int kid = field_kernel[f];
            auto read_off = [&](std::size_t begin, std::size_t end) {
                  Corners cc;
                  for (std::size_t i = begin; i < end; ++i) {
                        point_corners(grid, focal, n_focal, i, cc);
                        double v = 0.0;
                        for (int q = 0; q < cc.count; ++q) v += buf[cc.idx[q]] * cc.w[q];
                        if (exclude_self) v -= mf[i] * self_term[kid * n_pool + i];
                        of[i] = (std::fabs(v) < tol) ? 0.0 : v;
                  }
            };
            if (parallel_inner) {
                  RangeWorker<decltype(read_off)> w(read_off);
                  RcppParallel::parallelFor(0, n_focal, w, 4096);
            } else {
                  read_off(0, n_focal);
            }
      }
};

struct FieldWorker : public RcppParallel::Worker {
      const Engine& eng;
      FieldWorker(const Engine& eng_) : eng(eng_) {}
      void operator()(std::size_t begin, std::size_t end) {
            std::vector<double> buf(eng.grid.N);
            for (std::size_t f = begin; f < end; ++f) {
                  eng.run_field(static_cast<int>(f), buf, false);
            }
      }
};

// Leave-one-out self-term: the share of point j's own mass that binning
// followed by interpolation returns to j's location,
//   sum_{c, c'} w_c w_c' K(node_c - node_c'),
// with K the full product kernel at integer node offsets.
inline double self_term(const Engine& eng, const KernelSpec& ks, const Corners& c) {
      double s = 0.0;
      int off[5];
      for (int a = 0; a < c.count; ++a) {
            for (int b = 0; b < c.count; ++b) {
                  double k = 1.0;
                  for (size_t gi = 0; gi < eng.groups.size() && k != 0.0; ++gi) {
                        const Group& gr = eng.groups[gi];
                        for (int q = 0; q < gr.d; ++q) {
                              const int ax = gr.first + q;
                              off[q] = c.node[a][ax] - c.node[b][ax];
                        }
                        k *= kernel_at(gr, ks.arr[gi], off);
                  }
                  s += c.w[a] * c.w[b] * k;
            }
      }
      return s;
}

} // namespace fftgrid
} // namespace analogs


//' Evaluate kernel-convolved fields at focal points (FFT engine)
//'
//' @param pool_pts n_pool x D matrix of pool coordinates on the grid axes.
//' @param focal_pts n_focal x D matrix of focal coordinates on the grid axes.
//' @param origin,delta,n_nodes Per-axis grid origin, node spacing, node count.
//' @param group_sizes Number of axes in each kernel family's axis group, in
//'   axis order (groups are contiguous and cover all D axes).
//' @param pad_sizes Padded FFT length per axis.
//' @param kernels List over kernels; each a list over groups of numeric
//'   arrays (column-major, dims 2h+1 per group axis), the kernel's values at
//'   integer node offsets. All kernels share each group's dims.
//' @param masses n_pool x F matrix of field masses.
//' @param field_kernel Integer (0-based) kernel index per field.
//' @param exclude_self Subtract each focal's own binned contribution (focal i
//'   is pool point i).
//' @return n_focal x F matrix of field values.
//' @noRd
// [[Rcpp::export(.fft_fields_cpp)]]
NumericMatrix fft_fields_cpp(NumericMatrix pool_pts, NumericMatrix focal_pts,
                             NumericVector origin, NumericVector delta,
                             IntegerVector n_nodes, IntegerVector group_sizes,
                             IntegerVector pad_sizes, List kernels,
                             NumericMatrix masses, IntegerVector field_kernel,
                             bool exclude_self) {
      using namespace analogs::fftgrid;

      const int D = n_nodes.size();
      if (D > 5) stop("Internal error: FFT engine supports at most 5 grid axes.");
      if (pool_pts.ncol() != D || focal_pts.ncol() != D) {
            stop("Internal error: point matrices must have one column per grid axis.");
      }
      const int F = masses.ncol();
      if (masses.nrow() != pool_pts.nrow() || field_kernel.size() != F) {
            stop("Internal error: masses / field_kernel do not match the pool.");
      }
      if (exclude_self && focal_pts.nrow() != pool_pts.nrow()) {
            stop("Internal error: exclude_self requires focal points identical to pool points.");
      }

      Engine eng(origin, delta, n_nodes);
      eng.pool = REAL(pool_pts);   eng.n_pool = pool_pts.nrow();
      eng.focal = REAL(focal_pts); eng.n_focal = focal_pts.nrow();
      eng.mass = REAL(masses);
      eng.field_kernel.assign(field_kernel.begin(), field_kernel.end());
      eng.exclude_self = exclude_self;

      // Axis groups
      int first = 0;
      for (int g = 0; g < group_sizes.size(); ++g) {
            Group gr;
            gr.first = first;
            gr.d = group_sizes[g];
            gr.n.assign(n_nodes.begin() + first, n_nodes.begin() + first + gr.d);
            gr.m.assign(pad_sizes.begin() + first, pad_sizes.begin() + first + gr.d);
            gr.n_slice_nodes = 1; gr.m_real = 1;
            for (int a = 0; a < gr.d; ++a) {
                  gr.n_slice_nodes *= gr.n[a];
                  gr.m_real *= gr.m[a];
            }
            gr.m_cplx = (gr.m_real / gr.m[gr.d - 1]) * (gr.m[gr.d - 1] / 2 + 1);
            gr.n_slices = eng.grid.N / gr.n_slice_nodes;
            gr.inner_size = 1;
            for (int k = first + gr.d; k < D; ++k) gr.inner_size *= n_nodes[k];

            // pocketfft geometry (byte strides, row-major)
            gr.shape.assign(gr.m.begin(), gr.m.end());
            gr.stride_r.assign(gr.d, 0);
            gr.stride_c.assign(gr.d, 0);
            gr.axes.clear();
            ptrdiff_t sr = sizeof(double), sc = sizeof(cplx);
            for (int a = gr.d - 1; a >= 0; --a) {
                  gr.stride_r[a] = sr;
                  gr.stride_c[a] = sc;
                  sr *= gr.m[a];
                  sc *= (a == gr.d - 1) ? (gr.m[a] / 2 + 1) : gr.m[a];
            }
            for (int a = 0; a < gr.d; ++a) gr.axes.push_back(a);

            // Slice node (row-major over n) -> padded index (row-major over m)
            gr.pad_index.resize(gr.n_slice_nodes);
            std::vector<int> idx(gr.d, 0);
            for (size_t t = 0; t < gr.n_slice_nodes; ++t) {
                  size_t p = 0;
                  for (int a = 0; a < gr.d; ++a) p = p * gr.m[a] + idx[a];
                  gr.pad_index[t] = p;
                  for (int a = gr.d - 1; a >= 0; --a) {
                        if (++idx[a] < gr.n[a]) break;
                        idx[a] = 0;
                  }
            }
            eng.groups.push_back(gr);
            first += gr.d;
      }
      if (first != D) stop("Internal error: axis groups do not cover all grid axes.");

      // Kernels: per-group arrays, half-widths, and spectra
      const int n_kernels = kernels.size();
      eng.kernels.resize(n_kernels);
      for (int kk = 0; kk < n_kernels; ++kk) {
            List kl = kernels[kk];
            if (kl.size() != static_cast<int>(eng.groups.size())) {
                  stop("Internal error: kernel must have one array per axis group.");
            }
            KernelSpec& ks = eng.kernels[kk];
            ks.kmax = 1.0;
            for (size_t gi = 0; gi < eng.groups.size(); ++gi) {
                  Group& gr = eng.groups[gi];
                  NumericVector a = kl[gi];
                  IntegerVector dims = a.attr("dim");
                  if (dims.size() != gr.d) stop("Internal error: kernel array rank mismatch.");
                  std::vector<int> h(gr.d);
                  for (int q = 0; q < gr.d; ++q) h[q] = (dims[q] - 1) / 2;
                  if (kk == 0) {
                        gr.h = h;
                        for (int q = 0; q < gr.d; ++q) {
                              if (gr.m[q] < gr.n[q] + gr.h[q]) {
                                    stop("Internal error: padded FFT size too small for kernel support.");
                              }
                        }
                  } else if (h != gr.h) {
                        stop("Internal error: kernels must share each group's support.");
                  }
                  ks.arr.emplace_back(a.begin(), a.end());

                  double gmax = 0.0;
                  for (double v : ks.arr.back()) gmax = std::max(gmax, std::fabs(v));
                  ks.kmax *= gmax;

                  // Padded kernel with offset o stored at (o mod m), then r2c.
                  std::vector<double> kp(gr.m_real, 0.0);
                  const size_t n_arr = ks.arr.back().size();
                  std::vector<int> o(gr.d);
                  for (size_t lin = 0; lin < n_arr; ++lin) {
                        size_t rem = lin;
                        for (int q = 0; q < gr.d; ++q) {      // column-major decode
                              const int len = 2 * gr.h[q] + 1;
                              o[q] = static_cast<int>(rem % len) - gr.h[q];
                              rem /= len;
                        }
                        size_t p = 0;
                        for (int q = 0; q < gr.d; ++q) {      // row-major padded index
                              const int pos = (o[q] % gr.m[q] + gr.m[q]) % gr.m[q];
                              p = p * gr.m[q] + pos;
                        }
                        kp[p] = ks.arr.back()[lin];
                  }
                  std::vector<cplx> spec(gr.m_cplx);
                  pocketfft::r2c(gr.shape, gr.stride_r, gr.stride_c, gr.axes,
                                 pocketfft::FORWARD, kp.data(), spec.data(), 1.0);
                  ks.spec.push_back(std::move(spec));
            }
      }

      // Leave-one-out self-terms (per pool point, per kernel)
      if (exclude_self) {
            eng.self_term.assign(eng.n_pool * n_kernels, 0.0);
            Corners c;
            for (size_t j = 0; j < eng.n_pool; ++j) {
                  point_corners(eng.grid, eng.pool, eng.n_pool, j, c);
                  for (int kk = 0; kk < n_kernels; ++kk) {
                        eng.self_term[kk * eng.n_pool + j] = self_term(eng, eng.kernels[kk], c);
                  }
            }
      }

      NumericMatrix out(eng.n_focal, F);
      eng.out = REAL(out);

      // Parallelize over fields when every group is a single slice (one
      // family active); otherwise over slices within each field.
      bool single_slice = true;
      for (const Group& gr : eng.groups) single_slice = single_slice && (gr.n_slices == 1);
      try {
            if (single_slice) {
                  FieldWorker w(eng);
                  RcppParallel::parallelFor(0, static_cast<std::size_t>(F), w, 1);
            } else {
                  std::vector<double> buf(eng.grid.N);
                  for (int f = 0; f < F; ++f) eng.run_field(f, buf, true);
            }
      } catch (std::exception& e) {
            stop(std::string("FFT engine error: ") + e.what());
      }
      return out;
}


//' Turn per-focal field values into aggregate output columns (FFT engine)
//'
//' @param field_vals n_focal x F matrix from .fft_fields_cpp().
//' @param roles Named list of 0-based field indices (-1 = absent) giving
//'   which field feeds each accumulator slot; see .fft_build_fields().
//' @param aggregate_codes,se_code,n_vars,n_covs,n_classes_per_var,lambda As
//'   for query_analog_index_cpp().
//' @param has_values Whether y was supplied.
//' @param values_have_na Whether any y value is NA (per-variable denominators).
//' @param user_weight_mean Mean of positive user weights (1 if none); the
//'   ridge information measure is sum_w / user_weight_mean.
//' @return n_focal x n_cols matrix in the exact engine's column layout.
//' @noRd
// [[Rcpp::export(.fft_finalize_cpp)]]
NumericMatrix fft_finalize_cpp(NumericMatrix field_vals, List roles,
                               IntegerVector aggregate_codes, int se_code,
                               int n_vars, int n_covs,
                               IntegerVector n_classes_per_var,
                               double lambda, bool has_values,
                               bool values_have_na,
                               double user_weight_mean) {
      using namespace analogs;

      std::vector<AggregateCode> acodes(aggregate_codes.size());
      for (int s = 0; s < aggregate_codes.size(); ++s) {
            acodes[s] = static_cast<AggregateCode>(aggregate_codes[s]);
      }
      const SeCode sc = static_cast<SeCode>(se_code);
      const std::vector<int> ncls(n_classes_per_var.begin(), n_classes_per_var.end());
      const StatLayout L(acodes, sc, has_values, n_vars, n_covs, ncls, values_have_na);
      const int n_cols = L.n_cols();
      const size_t n_focal = field_vals.nrow();

      auto role = [&](const char* name) {
            IntegerVector r = roles[name];
            return std::vector<int>(r.begin(), r.end());
      };
      const int r_count = role("count")[0];
      const int r_sw = role("sw")[0];
      const int r_sw2 = role("sw2")[0];
      const std::vector<int> r_obs_count = role("obs_count");
      const std::vector<int> r_obs_sw = role("obs_sw");
      const std::vector<int> r_obs_sw2 = role("obs_sw2");
      const std::vector<int> r_value = role("value");
      const std::vector<int> r_wm_wy2 = role("wm_wy2");
      const std::vector<int> r_wm_w2y = role("wm_w2y");
      const std::vector<int> r_wm_w2y2 = role("wm_w2y2");
      const std::vector<int> r_tab = role("tab");
      const std::vector<int> r_xtwx = role("xtwx");
      const std::vector<int> r_xtwy = role("xtwy");
      const std::vector<int> r_reg_wy2 = role("reg_wy2");
      const std::vector<int> r_own = role("own");
      const std::vector<int> r_own_xtwx = role("own_xtwx");
      const std::vector<int> r_own_sw = role("own_sw");
      const std::vector<int> r_own_sw2 = role("own_sw2");

      const double* fv = REAL(field_vals);
      std::vector<double> agg(n_focal * static_cast<size_t>(n_cols), NA_REAL);
      const bool do_reg = L.has_regression && has_values && n_covs >= 0;
      const int dim = n_covs + 1;

      struct FinalizeWorker : public RcppParallel::Worker {
            std::function<void(size_t)> body;
            void operator()(std::size_t begin, std::size_t end) {
                  for (std::size_t i = begin; i < end; ++i) body(i);
            }
      } worker;

      worker.body = [&](size_t i) {
            auto get = [&](int r) { return (r < 0) ? 0.0 : fv[i + static_cast<size_t>(r) * n_focal]; };

            FocalAccum a;
            a.reset(L, has_values);
            a.count = get(r_count);
            a.sum_weights = get(r_sw);
            a.sum_weights_squared = get(r_sw2);
            for (int idx = 0; idx < L.n_regular; ++idx) {
                  const AggregateCode code = L.acodes[L.regular_pos[idx]];
                  if (code == AggregateCode::SUM_WEIGHTS || code == AggregateCode::MEAN_WEIGHTS) {
                        a.regular_accum[idx] = a.sum_weights;
                  }
            }
            for (size_t k = 0; k < a.value_accum.size(); ++k) a.value_accum[k] = get(r_value[k]);
            for (size_t v = 0; v < a.wm_sum_wy2.size(); ++v) a.wm_sum_wy2[v] = get(r_wm_wy2[v]);
            for (size_t v = 0; v < a.wm_sum_w2y.size(); ++v) a.wm_sum_w2y[v] = get(r_wm_w2y[v]);
            for (size_t v = 0; v < a.wm_sum_w2y2.size(); ++v) a.wm_sum_w2y2[v] = get(r_wm_w2y2[v]);
            a.wm_sum_w2 = a.sum_weights_squared;
            for (size_t b = 0; b < a.tabulate_accum.size(); ++b) a.tabulate_accum[b] = get(r_tab[b]);
            auto role_at = [](const std::vector<int>& r, size_t v) { return r.empty() ? -1 : r[v]; };
            for (size_t v = 0; v < a.obs_count.size(); ++v) {
                  a.obs_count[v]  = get(role_at(r_obs_count, v));
                  a.obs_sum_w[v]  = get(role_at(r_obs_sw, v));
                  a.obs_sum_w2[v] = get(role_at(r_obs_sw2, v));
            }

            double* row = agg.data() + i * static_cast<size_t>(n_cols);
            const int col_idx = write_nonregression_stats(L, a, has_values, row);

            if (!do_reg) return;
            if (!(a.sum_weights > 0.0)) {
                  write_regression_stats(L, nullptr, nullptr, row, col_idx);
                  return;
            }

            RidgeMoments m;
            m.reset(n_covs, n_vars);
            for (int q = 0; q < dim * dim; ++q) m.XtWX[q] = get(r_xtwx[q]);
            m.sum_w = a.sum_weights;
            m.sum_w_sq = a.sum_weights_squared;
            for (int q = 0; q < n_vars * dim; ++q) m.XtWy[q] = get(r_xtwy[q]);
            for (int v = 0; v < n_vars; ++v) {
                  const int slot_in = r_own.empty() ? -1 : r_own[v];
                  if (slot_in < 0) continue;
                  const int s = m.add_own(v);
                  for (int q = 0; q < dim * dim; ++q) {
                        m.own_XtWX[s * dim * dim + q] = get(r_own_xtwx[slot_in * dim * dim + q]);
                  }
                  m.own_sum_w[s] = get(r_own_sw[slot_in]);
                  m.own_sum_w_sq[s] = r_own_sw2.empty() ? 0.0 : get(r_own_sw2[slot_in]);
            }
            std::vector<double> sum_wy2(n_vars, 0.0);
            for (int v = 0; v < static_cast<int>(r_reg_wy2.size()); ++v) sum_wy2[v] = get(r_reg_wy2[v]);

            std::vector<double> coeffs(n_vars * dim);
            std::vector<double> se(L.reg_se ? n_vars * dim : 0);
            solve_ridge_from_moments(m, sum_wy2.data(), lambda,
                                     a.sum_weights / user_weight_mean, sc,
                                     coeffs.data(), L.reg_se ? se.data() : nullptr);
            write_regression_stats(L, coeffs.data(), L.reg_se ? se.data() : nullptr,
                                   row, col_idx);
      };

      RcppParallel::parallelFor(0, n_focal, worker);

      NumericMatrix res(n_focal, n_cols);
      for (size_t i = 0; i < n_focal; ++i) {
            for (int s = 0; s < n_cols; ++s) {
                  res(i, s) = agg[i * n_cols + s];
            }
      }
      return res;
}
