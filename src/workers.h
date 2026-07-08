#pragma once

#include <Rcpp.h>
#include <RcppParallel.h>
#include "types.h"
#include "lattice.h"
#include "geometry.h"
#include "climate.h"
#include "weights.h"
#include "mahalanobis.h"
#include "se_code.h"

#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_map>

using namespace RcppParallel;

namespace analogs {

// -------------------------------------------------------------------------
// Worker for pair-returning modes: knn_clim, knn_geog, all
//
// Writes into out_indices[i] a vector<int> of 1-based ref indices.
// Writes into the three parallel weight vectors:
//   out_sample_weights[i] - bin-level downsampling weights (>=1.0 typically)
//   out_area_weights[i]   - per-point cell-area weights (mean-1 normalized)
//   out_user_weights[i]   - per-point user-provided weights
// All three weight vectors are populated for every analog. The R-side
// decides which to emit as columns based on which mechanisms are in use.
// -------------------------------------------------------------------------
struct PairWorker : public Worker {
      const double* focal_ptr;      // focal matrix (original coords + clim)
      const double* ref_ptr;        // reference matrix (original coords + clim)
      const double* ref_latt_ptr;   // matrix used by lattice (may be ref_ptr or ECEF+clim)
      int n_focal;
      int n_ref;
      int n_clim;
      int stride_f;
      int stride_r;                 // stride for ref_ptr (n_ref)
      int stride_latt_r;            // stride for ref_latt_ptr

      bool use_lattice;
      bool use_geog_filter;
      bool use_haversine;
      bool use_scalar_clim;
      bool use_pervar_clim;
      bool use_ecef;                // true when lonlat + geog filtering → use ECEF

      double max_clim_scalar;
      double max_geog;
      double max_geog_chord;        // chord distance threshold for ECEF mode
      std::vector<double> max_clim_pervar;

      SelectCode scode;     // Selection strategy

      int k;                        // k for kNN modes (>=1), ignored for ALL
      Lattice* lattice_ptr;         // may be nullptr if !use_lattice

      double R_earth;               // Earth radius (km)

      // Mahalanobis distance support
      bool use_mahalanobis;
      const double* x_cov_ptr;
      int x_cov_stride;
      int n_cov_components;

      // Per-point pool weights (length n_ref, indexed by 0-based pool index).
      // When the corresponding has_* flag is false the pointer may be nullptr
      // and the worker treats the weight as 1.0.
      bool has_area_weight;
      const double* area_weight_ptr;
      bool has_user_weight;
      const double* user_weight_ptr;

      std::vector< std::vector<int> >&    out_indices;
      std::vector< std::vector<double> >& out_sample_weights;
      std::vector< std::vector<double> >& out_area_weights;
      std::vector< std::vector<double> >& out_user_weights;

      // Thread-local storage for reusable allocations
      struct ThreadLocalStorage {
            std::vector<double> fgeo_vec;
            std::vector<double> fclim_vec;
            std::vector<double> q_geo;
            std::vector<double> q_clim;
            std::vector<index_t> cand;
            std::vector<double> cand_weights;

            ThreadLocalStorage() {
                  fgeo_vec.reserve(3);
                  fclim_vec.reserve(16);
                  q_geo.reserve(3);
                  q_clim.reserve(16);
                  cand.reserve(256);
                  cand_weights.reserve(256);
            }
      };

      static thread_local ThreadLocalStorage tls;  // Defined in workers.cpp

      PairWorker(const Rcpp::NumericMatrix& focal_mm,
                 const Rcpp::NumericMatrix& ref_mm,
                 const double* ref_latt_ptr_,
                 int stride_latt_r_,
                 bool use_lattice_,
                 bool use_geog_filter_,
                 bool use_haversine_,
                 bool use_scalar_clim_,
                 bool use_pervar_clim_,
                 double max_clim_scalar_,
                 double max_geog_,
                 double max_geog_chord_,
                 const std::vector<double>& max_clim_pervar_,
                 SelectCode scode_,
                 int k_,
                 Lattice* lattice_ptr_,
                 bool use_ecef_,
                 double R_earth_,
                 bool use_mahalanobis_,
                 const double* x_cov_ptr_,
                 int x_cov_stride_,
                 int n_cov_components_,
                 bool has_area_weight_,
                 const double* area_weight_ptr_,
                 bool has_user_weight_,
                 const double* user_weight_ptr_,
                 std::vector< std::vector<int> >&    out_indices_,
                 std::vector< std::vector<double> >& out_sample_weights_,
                 std::vector< std::vector<double> >& out_area_weights_,
                 std::vector< std::vector<double> >& out_user_weights_)
            : focal_ptr(REAL(focal_mm)),
              ref_ptr(REAL(ref_mm)),
              ref_latt_ptr(ref_latt_ptr_),
              n_focal(focal_mm.nrow()),
              n_ref(ref_mm.nrow()),
              n_clim(focal_mm.ncol() - 2),
              stride_f(focal_mm.nrow()),
              stride_r(ref_mm.nrow()),
              stride_latt_r(stride_latt_r_),
              use_lattice(use_lattice_),
              use_geog_filter(use_geog_filter_),
              use_haversine(use_haversine_),
              use_scalar_clim(use_scalar_clim_),
              use_pervar_clim(use_pervar_clim_),
              use_ecef(use_ecef_),
              max_clim_scalar(max_clim_scalar_),
              max_geog(max_geog_),
              max_geog_chord(max_geog_chord_),
              max_clim_pervar(max_clim_pervar_),
              scode(scode_),
              k(k_),
              lattice_ptr(lattice_ptr_),
              R_earth(R_earth_),
              use_mahalanobis(use_mahalanobis_),
              x_cov_ptr(x_cov_ptr_),
              x_cov_stride(x_cov_stride_),
              n_cov_components(n_cov_components_),
              has_area_weight(has_area_weight_),
              area_weight_ptr(area_weight_ptr_),
              has_user_weight(has_user_weight_),
              user_weight_ptr(user_weight_ptr_),
              out_indices(out_indices_),
              out_sample_weights(out_sample_weights_),
              out_area_weights(out_area_weights_),
              out_user_weights(out_user_weights_)
      {}

      void operator()(std::size_t begin, std::size_t end);
};


// -------------------------------------------------------------------------
// Worker for aggregate modes: COUNT / SUM / MEAN / REGRESSION (potentially multiple)
// Writes into agg[i * n_stats + s] the scalar aggregate for focal i, stat s.
//
// Combined weight semantics:
//   combined_weight = dist_weight * sample_weight * area_weight * user_weight
// Each weight is 1.0 when its corresponding mechanism is inactive.
//
// COUNT preserves its existing semantics: count += static_cast<int>(sample_weight),
// adjusting for downsampling but unaffected by area/user weights. Users wanting
// area- or user-weighted counts should request stat = "sum_weights" with
// kernel = "uniform".
// -------------------------------------------------------------------------
struct AggWorker : public Worker {
      const double* focal_ptr;
      const double* ref_ptr;
      const double* ref_latt_ptr;   // ECEF+clim matrix if use_ecef, else ref_ptr
      int n_focal;
      int n_ref;
      int n_clim;
      int stride_f;
      int stride_r;
      int stride_latt_r;

      bool use_lattice;
      bool use_geog_filter;
      bool use_haversine;
      bool use_scalar_clim;
      bool use_pervar_clim;
      bool use_ecef;

      double max_clim_scalar;
      double max_geog;
      double max_geog_chord;
      std::vector<double> max_clim_pervar;

      SelectCode scode;
      std::vector<AggregateCode> acodes;
      FamilyKernel clim_kernel;
      FamilyKernel geog_kernel;
      double clim_wparam;
      double geog_wparam;

      Lattice* lattice_ptr;
      double R_earth;

      // Mahalanobis distance support
      bool use_mahalanobis;
      const double* x_cov_ptr;
      int x_cov_stride;
      int n_cov_components;

      // Values support for user-provided aggregation variables
      bool has_values;
      const double* values_ptr;
      int values_stride;
      int n_vars;

      // Covariates support for regression stat
      bool has_covariates;
      const double* covariates_ptr;
      int covariates_stride;
      int n_covs;
      double lambda;

      // Per-point pool weights (length n_ref, indexed by 0-based pool index).
      // When the has_* flag is false the pointer may be nullptr and the
      // worker treats the weight as 1.0.
      bool has_area_weight;
      const double* area_weight_ptr;
      bool has_user_weight;
      const double* user_weight_ptr;

      // SE variant (applies to weighted_mean and regression)
      SeCode se_code;

      // TABULATE: per-variable count of unique classes (length n_vars when
      // tabulate is requested, length 0 otherwise). Used to size and index
      // the per-focal class-vote accumulator.
      std::vector<int> n_classes_per_var;

      // Self-exclusion: when true, skip the pool row with the same index as
      // the focal. Only valid when x and pool are identical (enforced at the
      // R level). Used for leave-one-out cross-validation.
      bool exclude_self;

      std::vector<double>& agg; // flat output: size = n_focal * n_stats
      int n_stats;              // number of stats to compute (total columns in output)

      // Thread-local storage for reusable allocations
      struct ThreadLocalStorage {
            std::vector<double> q_geo;
            std::vector<double> q_clim;
            std::vector<index_t> cand;
            std::vector<double> cand_weights;

            ThreadLocalStorage() {
                  q_geo.reserve(3);
                  q_clim.reserve(16);
                  cand.reserve(256);
                  cand_weights.reserve(256);
            }
      };

      static thread_local ThreadLocalStorage tls;  // Defined in workers.cpp

      AggWorker(const Rcpp::NumericMatrix& focal_mm,
                const Rcpp::NumericMatrix& ref_mm,
                const double* ref_latt_ptr_,
                int stride_latt_r_,
                bool use_lattice_,
                bool use_geog_filter_,
                bool use_haversine_,
                bool use_scalar_clim_,
                bool use_pervar_clim_,
                double max_clim_scalar_,
                double max_geog_,
                double max_geog_chord_,
                const std::vector<double>& max_clim_pervar_,
                SelectCode scode_,
                const std::vector<AggregateCode>& acodes_,
                FamilyKernel clim_kernel_,
                FamilyKernel geog_kernel_,
                double clim_wparam_,
                double geog_wparam_,
                Lattice* lattice_ptr_,
                bool use_ecef_,
                double R_earth_,
                bool use_mahalanobis_,
                const double* x_cov_ptr_,
                int x_cov_stride_,
                int n_cov_components_,
                bool has_values_,
                const double* values_ptr_,
                int values_stride_,
                int n_vars_,
                bool has_covariates_,
                const double* covariates_ptr_,
                int covariates_stride_,
                int n_covs_,
                double lambda_,
                bool has_area_weight_,
                const double* area_weight_ptr_,
                bool has_user_weight_,
                const double* user_weight_ptr_,
                SeCode se_code_,
                const std::vector<int>& n_classes_per_var_,
                bool exclude_self_,
                std::vector<double>& agg_,
                int n_stats_)
            : focal_ptr(REAL(focal_mm)),
              ref_ptr(REAL(ref_mm)),
              ref_latt_ptr(ref_latt_ptr_),
              n_focal(focal_mm.nrow()),
              n_ref(ref_mm.nrow()),
              n_clim(focal_mm.ncol() - 2),
              stride_f(focal_mm.nrow()),
              stride_r(ref_mm.nrow()),
              stride_latt_r(stride_latt_r_),
              use_lattice(use_lattice_),
              use_geog_filter(use_geog_filter_),
              use_haversine(use_haversine_),
              use_scalar_clim(use_scalar_clim_),
              use_pervar_clim(use_pervar_clim_),
              use_ecef(use_ecef_),
              max_clim_scalar(max_clim_scalar_),
              max_geog(max_geog_),
              max_geog_chord(max_geog_chord_),
              max_clim_pervar(max_clim_pervar_),
              scode(scode_),
              acodes(acodes_),
              clim_kernel(clim_kernel_),
              geog_kernel(geog_kernel_),
              clim_wparam(clim_wparam_),
              geog_wparam(geog_wparam_),
              lattice_ptr(lattice_ptr_),
              R_earth(R_earth_),
              use_mahalanobis(use_mahalanobis_),
              x_cov_ptr(x_cov_ptr_),
              x_cov_stride(x_cov_stride_),
              n_cov_components(n_cov_components_),
              has_values(has_values_),
              values_ptr(values_ptr_),
              values_stride(values_stride_),
              n_vars(n_vars_),
              has_covariates(has_covariates_),
              covariates_ptr(covariates_ptr_),
              covariates_stride(covariates_stride_),
              n_covs(n_covs_),
              lambda(lambda_),
              has_area_weight(has_area_weight_),
              area_weight_ptr(area_weight_ptr_),
              has_user_weight(has_user_weight_),
              user_weight_ptr(user_weight_ptr_),
              se_code(se_code_),
              n_classes_per_var(n_classes_per_var_),
              exclude_self(exclude_self_),
              agg(agg_),
              n_stats(n_stats_)
      {}

      void operator()(std::size_t begin, std::size_t end);
};

} // namespace analogs
