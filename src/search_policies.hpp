#pragma once
#include "types.hpp"
#include "metrics.hpp"
#include <vector>
#include <utility>
#include <algorithm>

namespace analogs {

template <class MetricT>
struct RankByGeo {
      const MetricT& geo_metric;
      const MatrixView ref_geo;
      RankByGeo(const MetricT& m, MatrixView g) : geo_metric(m), ref_geo(g) {}

      void rank(const double* focal_geo, std::vector<index_t>& ids) const {
            std::vector<std::pair<double,index_t>> tmp;
            tmp.reserve(ids.size());
            for (auto id : ids) {
                  tmp.emplace_back(geo_metric.dist(focal_geo, &ref_geo(id,0), 2), id);
            }
            std::sort(tmp.begin(), tmp.end(),
                      [](auto& a, auto& b){ return a.first < b.first; });
            for (size_t i=0;i<ids.size();++i) ids[i]=tmp[i].second;
      }
};

template <class MetricT>
struct RankByClim {
      const MetricT& clim_metric;
      const MatrixView ref_clim; size_tu nvars;
      RankByClim(const MetricT& m, MatrixView c, size_tu nv) : clim_metric(m), ref_clim(c), nvars(nv) {}

      void rank(const double* focal_clim, std::vector<index_t>& ids) const {
            std::vector<std::pair<double,index_t>> tmp;
            tmp.reserve(ids.size());
            for (auto id : ids) {
                  tmp.emplace_back(clim_metric.dist(focal_clim, &ref_clim(id,0), nvars), id);
            }
            std::sort(tmp.begin(), tmp.end(),
                      [](auto& a, auto& b){ return a.first < b.first; });
            for (size_t i=0;i<ids.size();++i) ids[i]=tmp[i].second;
      }
};

} // namespace analogs
