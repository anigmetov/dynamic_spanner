#ifndef WASSERSTEIN_SPACE_POINT_H
#define WASSERSTEIN_SPACE_POINT_H

#include <iostream>
#include <utility>
#include <vector>
#include <limits>
#include <algorithm>
#include <math.h>
#include <random>
#include <ostream>

#include <cassert>

#include "spdlog/spdlog.h"

namespace spd = spdlog;

constexpr int INVALID_ID = -1;



template<class Real>
struct Pair_of_points_info {

      Real distance;
      bool distance_requested;
      bool exact_distance_used;
      Real upper_bound;
      Real lower_bound;

      Pair_of_points_info() { }

      Pair_of_points_info(Real distance)
              :distance(distance)
      {
          distance_requested = false;
          exact_distance_used = false;
          upper_bound = std::numeric_limits<Real>::max();
          lower_bound = 0.0;


      }

    };


template<class R>
std::ostream& operator<<(std::ostream& os, const Pair_of_points_info<R>& x)
{
  os << "(d = " << x.distance << ", exact_distance_used = " << x.exact_distance_used;
  os << ", low = " << x.lower_bound << ", upp = " << x.upper_bound << ")";
  return os;
}

template<class Real_ = double>
class DynamicSpanner {

public:

    using Real = Real_;

//protected:


public:

    using Matrix = std::vector<std::vector<Pair_of_points_info<Real>>>;

    using MatrixReal = typename std::vector<std::vector<Real>>;

    using VertexDescriptor = size_t;

    size_t m_num_points;
    Matrix m_matrix;
    MatrixReal m_distance_timings_matrix;
    double m_min_dist;
    double m_max_dist;
    std::mt19937 m_twister;
    double m_epsilon { 0.0 };
    double m_time_elapsed { 0.0 };
    std::string m_strategy;
    double m_hera_factor {  1.0 / (1.01) };

public:

    DynamicSpanner(const MatrixReal& distance_matrix, const MatrixReal& distance_timings_matrix)
            :
            m_distance_timings_matrix(distance_timings_matrix),
            m_min_dist(std::numeric_limits<double>::max()),
            m_max_dist(0.0),
            //m_twister(std::random_device()())
            m_twister(1)
    {
        m_num_points = distance_matrix.size();
        m_matrix.resize(m_num_points);

        for (size_t i = 0; i < m_num_points; i++) {
            m_matrix[i].resize(m_num_points);

            for (size_t j = 0; j < m_num_points; j++) {
                m_matrix[i][j] = Pair_of_points_info<Real>(distance_matrix[i][j]);
            }
        }

        // On the diagonal, bounds are exact
        for (size_t i = 0; i < m_num_points; i++) {
            m_matrix[i][i].lower_bound = m_matrix[i][i].upper_bound = 0.0;
            m_matrix[i][i].exact_distance_used = true; // not sure that is ever queried
        }

        for (int i = 0; i < m_num_points; ++i) {
            for (int j = i + 1; j < m_num_points; ++j) {
                m_max_dist = std::max(m_max_dist, distance_matrix[i][j]);
                m_min_dist = std::min(m_min_dist, distance_matrix[i][j]);
            }
        }
    }

    void make_exact_distance_keeper()
    {
        for(int i = 0; i < m_num_points; ++i) {
            for(int j = 0; j < m_num_points; ++j) {
                m_matrix[i][j].lower_bound = m_matrix[i][j].upper_bound = m_matrix[i][j].distance;
                m_matrix[i][j].distance_requested = false;
                m_matrix[i][j].exact_distance_used = false;
             }
        }
    }

    std::vector<std::pair<double, double>> get_buckets()
    {
        int n_buckets = ceil(log2(m_max_dist / m_min_dist));
        std::vector<std::pair<double, double>> result;
        result.reserve(n_buckets);
        double scaling_factor = m_min_dist;
        for (int k = 0; k < n_buckets; ++k) {
            result.emplace_back(scaling_factor * std::pow(2.0, k), 4.0 * scaling_factor * std::pow(2.0, k));
        }
        return result;
    }

    int get_bucket_index(const std::vector<std::pair<double, double>>& buckets, double value)
    {
        std::vector<int> cand_list;
        for (int i = 0; i < buckets.size(); ++i) {
            auto bucket = buckets[i];
            if (value >= buckets[i].first and value <= buckets[i].second) {
                cand_list.push_back(i);
            }
        }
        assert(not cand_list.empty());
        if (cand_list.size() == 1)
            return cand_list[0];
        else {
            std::random_shuffle(cand_list.begin(), cand_list.end());
            return cand_list[0];
        }
    }

    decltype(auto) weakly_sorted_edges()
    {
        auto buckets = get_buckets();
        int n_buckets = buckets.size();
        using EdgeInfo = std::tuple<double, int, int>;
        using EdgeVectorVec = std::vector<EdgeInfo>;
        std::vector<EdgeVectorVec> result_1(n_buckets, EdgeVectorVec());
        for (int i = 0; i < m_num_points; ++i) {
            for (int j = i + 1; j < m_num_points; ++j) {
                double d = m_matrix[i][j].distance;
                EdgeInfo ei{d, i, j};
                int bucket_idx = get_bucket_index(buckets, d);
                result_1[bucket_idx].emplace_back(d, i, j);
            }
        }
        for (auto& v : result_1)
            std::random_shuffle(v.begin(), v.end());
        EdgeVectorVec result;
        for (auto&& v : result_1) {
            for (auto&& x : v) {
                result.push_back(x);
            }
        }
        return result;
    }

    Real get_distance_no_cache(VertexDescriptor i, VertexDescriptor j) const
    {
        return m_matrix[i][j].distance;

    }

    // TODO; use to avoid case distinctions?
    struct Comparator {

      bool m_strict;

      Comparator(bool strict)
              :m_strict(strict) { }

      bool operator()(Real x, Real y)
      {
          if (m_strict) {
              return x < y;
          }
          return x <= y;
      }
    };

    bool is_distance_less(VertexDescriptor i, VertexDescriptor j,
            Real value, bool strict)
    {

        //std::cout << "Call less than " << i <<  " " << j << ", with dist " << m_matrix[i][j].distance << " val=" << value << " bounds: {" << m_matrix[i][j].lower_bound << "," << m_matrix[i][j].upper_bound << "]" << std::endl;

        if (i == j) {
            if (strict) {
                return 0.0 < value;
            }
            else {
                return 0.0 <= value;
            }
        }

        Pair_of_points_info<Real>& info = m_matrix[i][j];
        Pair_of_points_info<Real>& info_mirror = m_matrix[j][i];

        info.distance_requested = true;
        info_mirror.distance_requested = true;

        if (strict && info.upper_bound < value) {
            return true;
        }

        if (not strict && info.upper_bound <= value) {
            return true;
        }

        if (strict && info.lower_bound >= value) {
            return false;
        }

        if (not strict && info.lower_bound > value) {
            return false;
        }

        Real dist = get_distance(i, j);

        if (strict) {
            return dist < value;
        }
        else {
            return dist <= value;
        }

    }

    void update_bounds_using_distance(VertexDescriptor i, VertexDescriptor j, Real dist_ij = -1.0)
    {

        if (dist_ij == -1.0)
            dist_ij = m_matrix[i][j].distance;

//#pragma omp parallel
        {
//#pragma omp for schedule(dynamic)
            for (size_t x = 0; x < m_num_points; x++) {
                for (size_t y = x + 1; y < m_num_points; y++) {
                    Pair_of_points_info<Real>& info_xy = m_matrix[x][y];
                    Pair_of_points_info<Real>& info_yx = m_matrix[y][x];
                    if (info_xy.exact_distance_used) {
                        continue;
                    }

                    Real new_upper_bound_1;
                    Real new_upper_bound_1_1 = m_matrix[x][i].upper_bound;
                    Real new_upper_bound_1_2 = m_matrix[j][y].upper_bound;
                    if (new_upper_bound_1_1 == std::numeric_limits<Real>::max() or new_upper_bound_1_2 == std::numeric_limits<Real>::max()) {
                        new_upper_bound_1 = std::numeric_limits<Real>::max();
                    }
                    else {
                        new_upper_bound_1 = new_upper_bound_1_1 + dist_ij + new_upper_bound_1_2;
                    }
                    Real new_upper_bound_2;
                    Real new_upper_bound_2_1 = m_matrix[x][j].upper_bound;
                    Real new_upper_bound_2_2 = m_matrix[i][y].upper_bound;
                    if (new_upper_bound_2_1 == std::numeric_limits<Real>::max() or new_upper_bound_2_2 == std::numeric_limits<Real>::max()) {
                        new_upper_bound_2 = std::numeric_limits<Real>::max();
                    }
                    else {
                        new_upper_bound_2 = new_upper_bound_2_1 + dist_ij + new_upper_bound_2_2;
                    }
                    Real new_upper_bound = std::min(new_upper_bound_1, new_upper_bound_2);

                    //if (new_upper_bound < m_hera_factor * info_xy.distance) {
                    //    std::cerr << "ERROR HERE: " << "i = " << i << ", j = " << j << ", dist_ij = " << dist_ij << ", x = " << x << ", y = " << y << ", " << info_xy << ", new_ub_1  = " << new_upper_bound_1 << ", new ub2 = " << new_upper_bound_2 << std::endl;
                    //    throw std::runtime_error("bad upper bound");
                    //}

                    info_xy.upper_bound = std::min(info_xy.upper_bound, new_upper_bound);
                    info_yx.upper_bound = info_xy.upper_bound;
                }
            }
        }

//#pragma omp for schedule(dynamic)
        for (size_t x = 0; x < m_num_points; x++) {
            for (size_t y = x + 1; y < m_num_points; y++) {
                Pair_of_points_info<Real>& info_xy = m_matrix[x][y];
                Pair_of_points_info<Real>& info_yx = m_matrix[y][x];
                if (info_xy.exact_distance_used) {
                    continue;
                }

                Real new_lower_bound_1;
                Real new_lower_bound_1_1 = m_matrix[x][i].upper_bound;
                Real new_lower_bound_1_2 = m_matrix[j][y].upper_bound;
                //std::cout << x << " " << y << "Upper bounds xi " <<  new_lower_bound_1_1 << " jy " << new_lower_bound_1_2 << std::endl;
                if (new_lower_bound_1_1 == std::numeric_limits<Real>::max() or new_lower_bound_1_2 == std::numeric_limits<Real>::max()) {
                    new_lower_bound_1 = 0.0;
                }
                else {
                    new_lower_bound_1 = m_hera_factor *  dist_ij - new_lower_bound_1_1 - new_lower_bound_1_2;
                }
                Real new_lower_bound_2;
                Real new_lower_bound_2_1 = m_matrix[x][j].upper_bound;
                Real new_lower_bound_2_2 = m_matrix[i][y].upper_bound;
                //if (x == 52 and y == 98) std::cout << "Upper bounds xj " <<  new_lower_bound_2_1 << " iy " << new_lower_bound_2_2 << std::endl;
                if (new_lower_bound_2_1 == std::numeric_limits<Real>::max() or new_lower_bound_2_2 == std::numeric_limits<Real>::max()) {
                    new_lower_bound_2 = 0.0;
                }
                else {
                    new_lower_bound_2 = m_hera_factor * dist_ij - new_lower_bound_2_1 - new_lower_bound_2_2;
                }
                Real new_lower_bound = std::max(new_lower_bound_1, new_lower_bound_2);

                //if (new_lower_bound > info_xy.distance) {
                //    std::cerr << "ERROR HERE: " << "x = " << x << ", y = " << y << ", " << info_xy << ", new_lb_1  = " << new_lower_bound_1 << ", new_lb2 = " << new_lower_bound_2 << std::endl;
                //    std::cerr << "ERROR HERE: " << "i = " << i << ", j = " << j << ", dist_ij = " << dist_ij <<", with factor: "
                //        << m_hera_factor * dist_ij <<
                //        ", Upper bounds xj " <<  new_lower_bound_2_1 << " iy " << new_lower_bound_2_2 << std::endl;
                //    throw std::runtime_error("bad lower bound");
                //}

                info_xy.lower_bound = std::max(info_xy.lower_bound, new_lower_bound);
                info_yx.lower_bound = info_xy.lower_bound;
            }
        }

#if 1
//#pragma omp parallel
        {
//#pragma omp for schedule(dynamic)
            for (size_t x = 0; x < m_num_points; x++) {
                for (size_t y = x + 1; y < m_num_points; y++) {
                    Pair_of_points_info<Real>& info_xy = m_matrix[x][y];
                    Pair_of_points_info<Real>& info_yx = m_matrix[y][x];
                    if (info_xy.exact_distance_used) {
                        continue;
                    }

                    Real new_lower_bound_1;
                    Real new_lower_bound_1_1 = m_matrix[x][i].upper_bound;
                    Real new_lower_bound_1_2 = m_matrix[j][y].lower_bound;
                    //std::cout << x << " " << y << "Upper bounds xi " <<  new_lower_bound_1_1 << " jy " << new_lower_bound_1_2 << std::endl;
                    if (new_lower_bound_1_1 == std::numeric_limits<Real>::max() or new_lower_bound_1_2 == 0.0) {
                        new_lower_bound_1 = 0.0;
                    }
                    else {
                        new_lower_bound_1 = new_lower_bound_1_2 - dist_ij - new_lower_bound_1_1;
                    }
                    Real new_lower_bound_2;
                    Real new_lower_bound_2_1 = m_matrix[x][j].upper_bound;
                    Real new_lower_bound_2_2 = m_matrix[i][y].lower_bound;
                    //std::cout << "Upper bounds xj " <<  new_lower_bound_2_1 << " iy " << new_lower_bound_2_2 << std::endl;
                    if (new_lower_bound_2_1 == std::numeric_limits<Real>::max() or new_lower_bound_2_2 == 0.0) {
                        new_lower_bound_2 = 0.0;
                    }
                    else {
                        new_lower_bound_2 = new_lower_bound_2_2 - dist_ij - new_lower_bound_2_1;
                    }
                    Real new_lower_bound_3;
                    Real new_lower_bound_3_1 = m_matrix[j][y].upper_bound;
                    Real new_lower_bound_3_2 = m_matrix[x][i].lower_bound;
                    //std::cout << x << " " << y << "Upper bounds xi " <<  new_lower_bound_1_1 << " jy " << new_lower_bound_1_2 << std::endl;
                    if (new_lower_bound_3_1 == std::numeric_limits<Real>::max() or new_lower_bound_3_2 == 0.0) {
                        new_lower_bound_3 = 0.0;
                    }
                    else {
                        new_lower_bound_3 = new_lower_bound_3_2 - dist_ij - new_lower_bound_3_1;
                    }
                    Real new_lower_bound_4;
                    Real new_lower_bound_4_1 = m_matrix[i][y].upper_bound;
                    Real new_lower_bound_4_2 = m_matrix[x][j].lower_bound;
                    //std::cout << "Upper bounds xj " <<  new_lower_bound_2_1 << " iy " << new_lower_bound_2_2 << std::endl;
                    if (new_lower_bound_4_1 == std::numeric_limits<Real>::max() or new_lower_bound_4_2 == 0.0) {
                        new_lower_bound_4 = 0.0;
                    }
                    else {
                        new_lower_bound_4 = new_lower_bound_4_2 - dist_ij - new_lower_bound_4_1;
                    }

                    Real new_lower_bound = std::max(std::max(new_lower_bound_1, new_lower_bound_2),
                            std::max(new_lower_bound_3, new_lower_bound_4));


                    //if (new_lower_bound > info_xy.distance) {
                    //    std::cerr << "ERROR HERE: " << "x = " << x << ", y = " << y << ", " << info_xy << ", new_lb_1  = " << new_lower_bound_1 << ", new_lb2 = " << new_lower_bound_2 << ", new lb 3 = " << new_lower_bound_3 << ", new_lb_4 = " << new_lower_bound_4 << std::endl;
                    //    std::cerr << "info_xi = " << m_matrix[x][i] << std::endl;
                    //    std::cerr << "info_jy = " << m_matrix[j][y] << std::endl;
                    //    throw std::runtime_error("bad lower bound");
                    //}

                    info_xy.lower_bound = std::max(info_xy.lower_bound, new_lower_bound);
                    info_yx.lower_bound = info_xy.lower_bound;
                }

            }
        }

#endif

#if DEBUG
        for(size_t x=0;x<m_num_points;x++) {
            for(size_t y=0;y<m_num_points;y++) {
                Pair_of_points_info<Real>& info_xy = m_matrix[x][y];
                Pair_of_points_info<Real>& info_yx = m_matrix[y][x];
                assert(info_xy.upper_bound>=m_hera_factor * info_xy.distance);
                assert(info_xy.lower_bound<=info_xy.distance);
                assert(info_xy.upper_bound==info_yx.upper_bound);
                assert(info_xy.lower_bound==info_yx.lower_bound);
            }
        }
#endif

    }

    bool is_distance_greater(VertexDescriptor i, VertexDescriptor j, Real value, bool strict)
    {
        return not is_distance_less(i, j, value, not strict);
    }

    Real get_distance(VertexDescriptor i, VertexDescriptor j)
    {
        //console->debug("get_distance called for {} {}", i, j);
        Pair_of_points_info<Real>& info = m_matrix[i][j];
        Pair_of_points_info<Real>& info_mirror = m_matrix[j][i];
        info.distance_requested = true;
        info_mirror.distance_requested = true;
        info.exact_distance_used = true;
        info_mirror.exact_distance_used = true;
        info.upper_bound = info.distance;
        info.lower_bound = m_hera_factor * info.distance;
        info_mirror.upper_bound =  info.distance;
        info_mirror.lower_bound = m_hera_factor * info.distance;
        update_bounds_using_distance(i, j);
        m_time_elapsed += m_distance_timings_matrix[i][j];
        return info.distance;
    }

    Real get_distance_no_cache(VertexDescriptor i, VertexDescriptor j)
    {
        return m_matrix[i][j].distance;
    }

    size_t get_num_points() const
    {
        return m_num_points;
    }

    size_t get_number_of_computed_distances() const
    {
        size_t result = 0;
        for (size_t i = 0; i < m_num_points; i++) {
            for (size_t j = i + 1; j < m_num_points; j++) {
                if (m_matrix[i][j].exact_distance_used) {
                    result++;
                }
            }
        }
        return result;
    }

    double get_microseconds_to_compute_distances() const
    {
        return m_time_elapsed;
    }

    size_t get_number_of_requested_distances() const
    {
        size_t result = 0;
        for (size_t i = 0; i < m_num_points; i++) {
            for (size_t j = i + 1; j < m_num_points; j++) {
                if (m_matrix[i][j].distance_requested) {
                    result++;
                }
            }
        }
        return result;
    }

    double get_fraction_of_computed_distances() const
    {

        double num = (double) get_number_of_computed_distances();
        double denom = m_num_points * (m_num_points - 1) / 2.0;

        return num / denom;

    }

    double get_fraction_of_requested_distances() const
    {

        double num = (double) get_number_of_requested_distances();
        double denom = m_num_points * (m_num_points - 1) / 2.0;

        return num / denom;

    }

    std::string get_statistics() const
    {
        double n_pairs = m_num_points * (m_num_points - 1) / 2;
        double distances_computed = get_number_of_computed_distances();
        double sparseness_1 = distances_computed / n_pairs;
        double sparseness_2 = distances_computed / m_num_points;
        return fmt::format("{};{};{};{};{};{}", m_strategy, m_epsilon, m_num_points, distances_computed, sparseness_1, sparseness_2);
    }

    bool find_random_bad_ratio(int& i, int& j, bool connect_first, bool lower_bound_first, double epsilon)
    {
        std::vector<std::pair<size_t, size_t>> candidates;

        double result = 0.0;
        bool is_disconnected = false;
        bool is_lower_bound_zero = false;
        for (size_t i = 0; i < m_num_points; i++) {
            for (size_t j = i + 1; j < m_num_points; j++) {
                if (m_matrix[i][j].exact_distance_used)
                    continue;
                Real low = m_matrix[i][j].lower_bound;
                Real upp = m_matrix[i][j].upper_bound;
                if (upp == std::numeric_limits<Real>::max()) {
                    if (connect_first and not is_disconnected) {
                        candidates.clear();
                        is_disconnected = true;
                    }
                    candidates.emplace_back(i, j);
                }
                else if (not is_disconnected and low == 0.0) {
                    if (lower_bound_first and not is_lower_bound_zero) {
                        candidates.clear();
                        is_lower_bound_zero = true;
                    }
                    candidates.emplace_back(i, j);
                }
                else if (not is_disconnected and
                        not is_lower_bound_zero and
                        upp / low > 1.0 + epsilon) {
                    candidates.emplace_back(i, j);
                }
            }
        }

        if (candidates.empty())
            return false;

        std::uniform_int_distribution<int> dis(0, candidates.size()-1);
        int idx = dis(m_twister);
        i = candidates[idx].first;
        j = candidates[idx].second;
        return true;
    }

    void check_ratio(double _eps)
    {
        int n_used = 0;
        for (size_t i = 0; i<m_num_points; i++) {
            for (size_t j = i+1; j<m_num_points; j++) {
                Real true_dist = m_matrix[i][j].distance;
                n_used += m_matrix[i][j].distance_requested;
                Real low = m_matrix[i][j].lower_bound;
                Real upp = m_matrix[i][j].upper_bound;
                Real approx_dist = upp;
                //if ( fabs(true_dist - approx_dist) / true_dist > _eps) {
                //    std::cerr << "ERROR -1 : " << i << ", " << j << " , low = " << low << ", true = " << true_dist << ", upp = " << upp << std::endl;
                //    throw std::runtime_error("Bad spanner - bad ratio");
                //}

                //if (low >true_dist || upp < true_dist) {
                //    std::cerr << "ERROR - 2: " << i << ", " << j << " , low = " << low << ", true = " << true_dist << ", upp = " << upp << std::endl;
                //    throw std::runtime_error("Bad spanner - bad bounds");
                //}
            }
        }
        std::cout << "CHECKED SPANNER, n_used = " << n_used << std::endl;
    }

    void find_worst_ratio(size_t& x, size_t& y, double& ratio)
    {
        std::vector<std::pair<size_t, size_t>> candidates;

        double result = 0.0;
        for (size_t i = 0; i < m_num_points; i++) {
            for (size_t j = i + 1; j < m_num_points; j++) {
                Real low = m_matrix[i][j].lower_bound;
                Real upp = m_matrix[i][j].upper_bound;
                if (low == 0.0 or upp == std::numeric_limits<Real>::max()) {
                    if (result != std::numeric_limits<Real>::max()) {
                        candidates.clear();
                    }
                    result = std::numeric_limits<Real>::max();
                    candidates.push_back(std::make_pair(i, j));
                }
                else {
                    double curr = upp / low;
                    if (curr > result) {
                        candidates.clear();
                    }
                    if (curr >= result) {
                        candidates.push_back(std::make_pair(i, j));
                        result = curr;
                    }
                }
            }
        }
        std::random_shuffle(candidates.begin(), candidates.end());
        x = candidates.begin()->first;
        y = candidates.begin()->second;
        ratio = result;
//        std::cout << "Worst ratio at " << x << ", " << y << ": " << m_matrix[x][y].lower_bound << " " << m_matrix[x][y].upper_bound << std::endl;
    }

    double find_worst_ratio()
    {
        double result = 0.0;
        for (size_t i = 0; i < m_num_points; i++) {
            for (size_t j = i + 1; j < m_num_points; j++) {
                Real low = m_matrix[i][j].lower_bound;
                Real upp = m_matrix[i][j].upper_bound;
                if (low == 0.0 or upp == std::numeric_limits<Real>::max()) {
                    return std::numeric_limits<Real>::max();
                }
                result = std::max(result, upp / low);
            }
        }
        return result;
    }

    void construct_blind_greedy_eps_spanner(double eps)
    {
        m_epsilon = eps;
        m_strategy = "blind-greedy";
        size_t i, j;
        double ratio;

        find_worst_ratio(i, j, ratio);

        while (ratio > (1 + eps)) {
            get_distance(i, j);
            find_worst_ratio(i, j, ratio);
        }
    }

    void construct_blind_random_eps_spanner(double eps)
    {
        m_epsilon = eps;
        m_strategy = "blind-random";

        size_t i, j;
        double ratio;

        ratio = find_worst_ratio();

        std::vector<std::pair<size_t, size_t>> vec;

        for (size_t i = 0; i < m_num_points; i++) {
            for (size_t j = i + 1; j < m_num_points; j++) {
                vec.emplace_back(i, j);
            }
        }

        std::random_shuffle(vec.begin(), vec.end());
        auto vec_iter = vec.begin();

        while (ratio > (1 + eps)) {
            std::tie(i, j) = *vec_iter++;
            get_distance(i, j);
            //std::cout << "adding i = " << i << ", j = " << j << std::endl;
            ratio = find_worst_ratio();
        }
    }

    void construct_blind_quasi_sorted_greedy_eps_spanner(double eps)
    {
        m_epsilon = eps;
        m_strategy = "blind-quasi-sorted-greedy";

        update_bounds_from_approximate_distances(2.0);
        auto quasi_sorted_edges = weakly_sorted_edges();
        int i, j;
        int edge_idx = 0;
        double dist;
        double ratio = find_worst_ratio();

        while (ratio > (1 + eps)) {
            std::tie(dist, i, j) = quasi_sorted_edges[edge_idx++];
            get_distance(i, j);
            ratio = find_worst_ratio();
        }
    }

    void construct_blind_quasi_sorted_shaker_eps_spanner(double eps)
    {
        m_epsilon = eps;
        m_strategy = "blind-quasi-sorted-shaker";

        update_bounds_from_approximate_distances(2.0);
        auto quasi_sorted_edges = weakly_sorted_edges();
        int i, j;
        int edge_idx_low = 0;
        int edge_idx_high = quasi_sorted_edges.size() - 1;
        double dist;
        double ratio = find_worst_ratio();

        bool go_up = false;

        while (ratio > (1 + eps)) {
            if (go_up)
                std::tie(dist, i, j) = quasi_sorted_edges[edge_idx_low++];
            else
                std::tie(dist, i, j) = quasi_sorted_edges[edge_idx_high--];
//            std::cout << "quasi_shaker: dist = " << dist << ", go up = " << go_up << ", edge_idx_low = " << edge_idx_low << ", edge_idx_high = " << edge_idx_high << std::endl;
            go_up = !go_up;
            get_distance(i, j);
            ratio = find_worst_ratio();
        }
    }

    void construct_blind_random_ratio_eps_spanner(bool connect_first, bool lower_bound_first, double eps)
    {
        m_epsilon = eps;
        if (connect_first and lower_bound_first)
            m_strategy = "blind-random-bad-ratio-connect-first-lower-bound-first";
        else if (connect_first)
            m_strategy = "blind-random-bad-ratio-connect-first";
        else if (lower_bound_first)
            m_strategy = "blind-random-bad-ratio-lower-bound-first";
        else
            m_strategy = "blind-random-bad-ratio";

        int i, j;
        while (true) {
            find_random_bad_ratio(i, j, connect_first, lower_bound_first, eps);
            if (find_worst_ratio() <= 1.0 + eps)
                break;
            get_distance(i, j);
        }
    }

    void update_bounds_from_approximate_distances(const double approximation_factor)
    {

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(1.0, approximation_factor);


//#pragma omp parallel
        {
//#pragma omp for schedule(dynamic)
            for (int i = 0; i < m_num_points; ++i) {
                for (int j = i + 1; j < m_num_points; ++j) {
                    auto& info_ij = m_matrix[i][j];
                    auto& info_ji = m_matrix[j][i];
                    double approx_value = dis(gen) * info_ij.distance;
                    info_ij.upper_bound = approx_value;
                    info_ij.lower_bound = (1.0 / approximation_factor) * approx_value;
                    info_ji.upper_bound = approx_value;
                    info_ji.lower_bound = (1.0 / approximation_factor) * approx_value;
                }
            }
        }

//        for (int i = 0; i < m_num_points; ++i) {
//            for (int j = i + 1; j < m_num_points; ++j) {
//                update_bounds_using_distance(i, j);
//            }
//        }
    }

    struct Value_with_index {
      double v;
      size_t i;
      size_t j;

      Value_with_index(double d, size_t i, size_t j)
              :v(d), i(i), j(j) { }
    };

    struct Comp_value_with_index {

      bool operator()(Value_with_index& a, Value_with_index& b)
      {
          return a.v < b.v;
      }

      Comp_value_with_index() { }
    };

    void construct_greedy_eps_spanner(double eps)
    {
        m_epsilon = eps;
        m_strategy = "greedy-non-blind";

        std::vector<Value_with_index> vec;

        for (size_t i = 0; i < m_num_points; i++) {
            for (size_t j = i + 1; j < m_num_points; j++) {
                vec.push_back(Value_with_index(m_matrix[i][j].distance, i, j));
            }
        }

        std::sort(vec.begin(), vec.end(), Comp_value_with_index());
        //std::reverse(vec.begin(),vec.end());

        //std::random_shuffle(vec.begin(),vec.end());

        for (size_t k = 0; k < vec.size(); k++) {
            Pair_of_points_info<Real>& info = m_matrix[vec[k].i][vec[k].j];

            //std::cout << "New pair " << vec[k].i << " " << vec[k].j << " with dist " << info.distance << " bounds: " << info.lower_bound << " " << info.upper_bound << " ratio: " << info.upper_bound / info.distance << std::endl;

            if (info.upper_bound / info.distance > (1 + eps)) {
                get_distance(vec[k].i, vec[k].j);
            }
        }
    }

    void print_ratios()
    {
        for (size_t i = 0; i < m_num_points; i++) {
            for (size_t j = 0; j < m_num_points; j++) {
                if (i == j) {
                    std::cout << "0 ";
                }
                else {
                    Real low = m_matrix[i][j].lower_bound;
                    Real upp = m_matrix[i][j].upper_bound;
                    if (low == 0.0 or upp == std::numeric_limits<Real>::max()) {
                        std::cout << "inf ";
                    }
                    else {
                        std::cout << upp / low << " ";
                    }
                }
            }
            std::cout << std::endl;
        }
    }

};


#endif // WASSERSTEIN_SPACE_POINT_H

