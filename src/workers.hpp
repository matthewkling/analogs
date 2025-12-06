#pragma once

#include <Rcpp.h>
#include <RcppParallel.h>
#include "types.hpp"
#include "lattice.hpp"
#include "geometry.hpp"
#include "climate.hpp"
#include "weights.hpp"
#include "mahalanobis.hpp"

#include <vector>
#include <queue>
#include <algorithm>

using namespace RcppParallel;

namespace analogs {

// -------------------------------------------------------------------------
// Worker for pair-returning modes: knn_clim, knn_geog, all
// Writes into out_indices[i] a vector<int> of 1-based ref indices.
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
      bool use_mahalanobis;                           // true if x_cov provided
      const double* x_cov_ptr;                        // pointer to x_cov matrix (n_focal × n_cov_components)
      int x_cov_stride;                               // stride for x_cov
      int n_cov_components;                           // n_clim * (n_clim + 1) / 2

      std::vector< std::vector<int> >& out_indices;

      // Thread-local storage for reusable allocations
      struct ThreadLocalStorage {
            std::vector<double> fgeo_vec;
            std::vector<double> fclim_vec;
            std::vector<double> q_geo;
            std::vector<double> q_clim;
            std::vector<index_t> cand;

            ThreadLocalStorage() {
                  // Pre-allocate with reasonable sizes
                  fgeo_vec.reserve(3);   // max 3 for ECEF
                  fclim_vec.reserve(16); // typical climate dims
                  q_geo.reserve(3);
                  q_clim.reserve(16);
                  cand.reserve(256);     // typical candidate count
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
                 std::vector< std::vector<int> >& out_indices_)
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
              out_indices(out_indices_)
      {}

      void operator()(std::size_t begin, std::size_t end);
};


// -------------------------------------------------------------------------
// Worker for aggregate modes: COUNT / SUM / MEAN (potentially multiple)
// Writes into agg[i * n_stats + s] the scalar aggregate for focal i, stat s.
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
      bool use_ecef;                // use ECEF for geog filtering

      double max_clim_scalar;
      double max_geog;
      double max_geog_chord;        // chord threshold for ECEF
      std::vector<double> max_clim_pervar;

      SelectCode scode;             // Selection strategy
      std::vector<AggregateCode> acodes;  // Which aggregations to perform
      WeightCode wcode;             // Weight function (for sum/mean)
      double weight_param1;         // Pre-computed weight parameter 1
      double weight_param2;         // Pre-computed weight parameter 2

      Lattice* lattice_ptr;
      double R_earth;

      // Mahalanobis distance support
      bool use_mahalanobis;                           // true if x_cov provided
      const double* x_cov_ptr;                        // pointer to x_cov matrix
      int x_cov_stride;                               // stride for x_cov
      int n_cov_components;                           // n_clim * (n_clim + 1) / 2

      std::vector<double>& agg; // flat output: size = n_focal * n_stats
      int n_stats;              // number of stats to compute

      // Thread-local storage for reusable allocations
      struct ThreadLocalStorage {
            std::vector<double> q_geo;
            std::vector<double> q_clim;
            std::vector<index_t> cand;

            ThreadLocalStorage() {
                  q_geo.reserve(3);
                  q_clim.reserve(16);
                  cand.reserve(256);
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
                WeightCode wcode_,
                double weight_param1_,
                double weight_param2_,
                Lattice* lattice_ptr_,
                bool use_ecef_,
                double R_earth_,
                bool use_mahalanobis_,
                const double* x_cov_ptr_,
                int x_cov_stride_,
                int n_cov_components_,
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
              wcode(wcode_),
              weight_param1(weight_param1_),
              weight_param2(weight_param2_),
              lattice_ptr(lattice_ptr_),
              R_earth(R_earth_),
              use_mahalanobis(use_mahalanobis_),
              x_cov_ptr(x_cov_ptr_),
              x_cov_stride(x_cov_stride_),
              n_cov_components(n_cov_components_),
              agg(agg_),
              n_stats(n_stats_)
      {}

      void operator()(std::size_t begin, std::size_t end);
};

} // namespace analogs
