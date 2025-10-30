#include <Rcpp.h>
#include "types.hpp"
#include "metrics.hpp"
#include "constraints.hpp"
#include "lattice_index.hpp"

using namespace Rcpp;
using analogs::size_tu;
using analogs::index_t;

// ---- helpers -------------------------------------------------------------
static inline analogs::MatrixView as_view_mat(NumericMatrix M) {
      analogs::MatrixView v;
      v.data = REAL(M);
      v.nrow = (size_tu)M.nrow();
      v.ncol = (size_tu)M.ncol();
      return v;
}

static inline bool is_null(SEXP x) {
      return x == R_NilValue;
}

static inline analogs::GeoSpace parse_geo_flag(SEXP geo_flag) {
      if (is_null(geo_flag)) return analogs::GeoSpace::LonLat;
      if (TYPEOF(geo_flag) == STRSXP && Rf_length(geo_flag) == 1) {
            SEXP s = STRING_ELT(geo_flag, 0);
            if (s == NA_STRING) return analogs::GeoSpace::LonLat;
            std::string g = CHAR(s);
            if (g=="lonlat"||g=="LonLat"||g=="LL") return analogs::GeoSpace::LonLat;
            if (g=="projected"||g=="Planar"||g=="XY") return analogs::GeoSpace::Planar;
      }
      if (TYPEOF(geo_flag) == INTSXP || TYPEOF(geo_flag) == LGLSXP) {
            int v = as<int>(geo_flag);
            return (v==1)? analogs::GeoSpace::Planar : analogs::GeoSpace::LonLat;
      }
      return analogs::GeoSpace::LonLat;
}

// ---- core impl -----------------------------------------------------------
static SEXP find_analogs_core_impl(
            SEXP focal_clim_, SEXP ref_clim_,
            SEXP focal_geo_, SEXP ref_geo_,
            SEXP mode_, SEXP k_,
            SEXP climate_band_, SEXP radius_km_,
            SEXP geo_flag_, SEXP compact_bins_) {

      // Parse matrices
      NumericMatrix focal_clim(focal_clim_);
      NumericMatrix ref_clim(ref_clim_);
      const size_tu n_focal = focal_clim.nrow();
      const size_tu n_ref = ref_clim.nrow();
      const size_tu n_vars = ref_clim.ncol();

      if (focal_clim.ncol() != (int)n_vars) {
            stop("focal_clim and ref_clim must have same number of columns (variables).");
      }

      // Parse geographic data (optional)
      NumericMatrix focal_geo, ref_geo;
      bool have_geo = !is_null(focal_geo_) && !is_null(ref_geo_);
      if (have_geo) {
            focal_geo = NumericMatrix(focal_geo_);
            ref_geo = NumericMatrix(ref_geo_);
            if (focal_geo.ncol()!=2 || ref_geo.ncol()!=2) {
                  stop("focal_geo and ref_geo must have 2 columns (lon,lat or x,y).");
            }
            if (focal_geo.nrow()!=(int)n_focal || ref_geo.nrow()!=(int)n_ref) {
                  stop("focal_geo/ref_geo row counts must match focal_clim/ref_clim.");
            }
      }

      // Parse climate constraint: scalar radius OR per-var band
      analogs::ClimateConstraint clim;
      clim.nvars = n_vars;
      clim.use_scalar = true;
      clim.radius = R_PosInf;
      std::vector<double> band_vec;

      if (!is_null(climate_band_)) {
            NumericVector v(climate_band_);
            if ((size_tu)v.size() == 1) {
                  clim.use_scalar = true;
                  clim.radius = v[0];
            } else if ((size_tu)v.size() == n_vars) {
                  band_vec.assign(v.begin(), v.end());
                  clim.use_scalar = false;
                  clim.max_abs_diff = band_vec.data();
            } else {
                  stop("`max_clim` must be length 1 (scalar radius) or length ncol(*_clim).");
            }
      }

      // Parse radius_km (not fully implemented yet)
      if (!is_null(radius_km_)) {
            (void)as<double>(radius_km_);  // parse but don't use yet
      }

      // Parse k parameter
      int k = 0;
      if (!is_null(k_)) k = as<int>(k_);

      // Build lattice index with quantile binning
      // Note: we're using default target_occupancy of 20
      analogs::MatrixView ref_clim_v = as_view_mat(ref_clim);
      analogs::MatrixView ref_geo_v;
      if (have_geo) ref_geo_v = as_view_mat(ref_geo);

      analogs::LatticeIndex index(ref_clim_v, ref_geo_v, n_vars);
      index.build();

      // Get diagnostics for reporting
      size_tu total_bins;
      double avg_occ, min_occ, max_occ;
      index.get_diagnostics(total_bins, avg_occ, min_occ, max_occ);

      // Process each focal point
      List out(n_focal);
      std::vector<index_t> ids;
      ids.reserve(n_ref);
      std::vector<double> focal_buf(n_vars);

      for (size_tu i = 0; i < n_focal; ++i) {
            // Extract focal climate values
            for (size_tu kk = 0; kk < n_vars; ++kk) {
                  focal_buf[kk] = focal_clim(i, kk);
            }

            // Search for analogs
            if (k > 0) {
                  index.best_first_search(focal_buf.data(), clim, ids, (index_t)k);
            } else {
                  index.hard_filter(focal_buf.data(), clim, ids);
            }

            // Convert to 1-based indices for R
            IntegerVector idx(ids.size());
            for (size_tu j = 0; j < ids.size(); ++j) {
                  idx[j] = ids[j] + 1;
            }
            out[i] = idx;
      }

      // Add attributes for diagnostics
      out.attr("n_ref") = (int)n_ref;
      out.attr("n_vars") = (int)n_vars;
      out.attr("total_bins") = (int)total_bins;
      out.attr("avg_bin_occupancy") = avg_occ;
      out.attr("min_bin_occupancy") = min_occ;
      out.attr("max_bin_occupancy") = max_occ;
      out.attr("geo") = "lonlat";  // default for now
      out.attr("binning_method") = "quantile";

      return out;
}

// [[Rcpp::export]]
SEXP find_analogs_core(
            SEXP focal_clim_, SEXP ref_clim_,
            SEXP focal_geo_, SEXP ref_geo_,
            SEXP mode_, SEXP k_,
            SEXP climate_band_, SEXP radius_km_,
            SEXP geo_flag_, SEXP compact_bins_) {
      return find_analogs_core_impl(
            focal_clim_, ref_clim_,
            focal_geo_, ref_geo_,
            mode_, k_,
            climate_band_, radius_km_,
            geo_flag_, compact_bins_);
}
