#ifndef SPANNER_WASSERSTEIN_LOWER_BOUND_H
#define SPANNER_WASSERSTEIN_LOWER_BOUND_H

#include "basic_defs_ws.h"
#include "dnn/geometry/euclidean-fixed.h"
#include "dnn/local/kd-tree.h"

namespace hera {
    namespace ws {

        // match every bidder to the nearest neigbour to get lower bound for Wasserstein distance
        template<class Real, class PointContainer>
        Real
        lower_bound_(const PointContainer& bidders, const PointContainer& items, const AuctionParams <Real>& params)
        {

            const Real q = params.wasserstein_power;

            const bool is_bottleneck = (q == std::numeric_limits<Real>::infinity() or q == hera::get_infinity());

            Real result { 0.0 };

            using DnnPoint          = dnn::Point<2, Real>;
            using DnnTraits         = dnn::PointTraits<DnnPoint>;

            std::vector<DnnPoint> dnn_points;
            std::vector<DnnPoint*> dnn_point_handles;

            IdxType true_idx { 0 };

            dnn_points.clear();
            dnn_points.reserve(items.size());
            // store normal items in kd-tree
            for (const DiagramPoint<Real>& g : items) {
                if (g.is_normal()) {
                    // index of items is id of dnn-point
                    DnnPoint p(true_idx);
                    p[0] = g.getRealX();
                    p[1] = g.getRealY();
                    dnn_points.push_back(p);
                }
                true_idx++;
            }

            for (size_t i = 0; i < dnn_points.size(); ++i) {
                dnn_point_handles.push_back(&dnn_points[i]);
            }

            DnnTraits traits;
            traits.internal_p = params.internal_p;
            dnn::KDTree<DnnTraits> kdtree { traits, dnn_point_handles, is_bottleneck ? 1.0 : q };

            for (const DiagramPoint<Real>& bidder : bidders) {
                if (bidder.is_diagonal())
                    continue;
                DnnPoint bidder_dnn;
                bidder_dnn[0] = bidder.getRealX();
                bidder_dnn[1] = bidder.getRealY();
                auto best_item = kdtree.findK(bidder_dnn, 1);
                if (is_bottleneck)
                    result = std::max(result, std::min(best_item[0].d, std::pow(bidder.persistence_lp(params.internal_p), q)));
                else
                    result += std::min(best_item[0].d, std::pow(bidder.persistence_lp(params.internal_p), q));
            }

            if (not is_bottleneck)
                result = std::pow(result, 1.0 / q );

            return result;
        }

        template<class Real, class PointContainer>
        Real lower_bound(const PointContainer& bidders, const PointContainer& items, const AuctionParams <Real>& params)
        {
            return std::max(lower_bound_<Real, PointContainer>(bidders, items, params),
                            lower_bound_<Real, PointContainer>(items, bidders, params));
        }

    }
}

#endif //SPANNER_WASSERSTEIN_LOWER_BOUND_H
