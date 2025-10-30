#include <Rcpp.h>
#include "metrics.hpp"
#include "types.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
double analogs_euclid_cpp(NumericVector a, NumericVector b) {
      if (a.size() != b.size()) stop("Length mismatch");
      analogs::Euclidean m;
      return m.dist(a.begin(), b.begin(), static_cast<analogs::size_tu>(a.size()));
}

// [[Rcpp::export]]
double analogs_haversine_cpp(NumericVector lonlat_a, NumericVector lonlat_b) {
      if (lonlat_a.size() != 2 || lonlat_b.size() != 2) {
            stop("Expect vectors of length 2: c(lon, lat) in degrees");
      }
      analogs::Haversine2D h;
      return h.dist(lonlat_a.begin(), lonlat_b.begin(), 2);
}
