#ifndef NEAREST_NEIGHBOR_SEARCH_H
#define NEAREST_NEIGHBOR_SEARCH_H

#include <iostream>
#include <utility>
#include <vector>
#include <limits>
#include <algorithm>

#include "dynamic_spanner.h"

template<class Real_ = double>
class Nearest_neighbor_search {

public:

    using Real = Real_;

protected:

    struct Point_info {

      Real distance;
      bool exact_distance_used;
      Real upper_bound;
      Real lower_bound;

      Point_info() { }

      Point_info(Real distance)
              :distance(distance)
      {
          exact_distance_used = false;
          upper_bound = std::numeric_limits<Real>::max();
          lower_bound = 0.0;
      }
    };

public:

    using Vector = std::vector<Point_info>;

    using VectorReal = typename std::vector<Real>;

    using VertexDescriptor = size_t;

    size_t m_num_points;
    Vector m_vector;

    DynamicSpanner<Real>* m_dspanner;

    const Real MAX = std::numeric_limits<Real>::max();

public:

    Nearest_neighbor_search(const VectorReal& distance_vector,
            DynamicSpanner<Real>* ptr_dspanner)
            :m_dspanner(ptr_dspanner)
    {

        m_num_points = m_dspanner->get_num_points();
        m_vector.resize(m_num_points);

        for (size_t i = 0; i < m_num_points; i++) {
            m_vector[i] = Point_info(distance_vector[i]);
        }
    }

    VertexDescriptor get_pos_of_smallest_upper_bound()
    {
        Real best_bound = MAX;
        VertexDescriptor best_index = 0;
        for (size_t i = 0; i < m_num_points; i++) {
            if (m_vector[i].upper_bound < best_bound) {
                best_index = i;
                best_bound = m_vector[i].upper_bound;
            }
        }
        return best_index;
    }

    VertexDescriptor find_nearest_neighbor()
    {

        while (true) {
            VertexDescriptor curr = get_pos_of_smallest_upper_bound();
            std::vector<VertexDescriptor> candidates_smaller_than_curr;

            for (size_t i = 0; i < m_num_points; i++) {
                if (i == curr) {
                    continue;
                }
                if (m_vector[i].lower_bound < m_vector[curr].upper_bound) {
                    // i could still be closer than curr
                    candidates_smaller_than_curr.push_back(i);
                }
            }

            if (candidates_smaller_than_curr.empty()) {
                return curr;
            }

            if (not m_vector[curr].exact_distance_used) {
                candidates_smaller_than_curr.push_back(curr);
            }
            std::random_shuffle(candidates_smaller_than_curr.begin(),
                    candidates_smaller_than_curr.end());
            get_distance_to(*candidates_smaller_than_curr.begin());

        }
    }

    void get_distance_to(VertexDescriptor i)
    {

        Point_info& info = m_vector[i];
        info.exact_distance_used = true;
        info.upper_bound = info.lower_bound = info.distance;
        update_bounds_using_distance(i);
    }

    void update_bounds_using_distance(VertexDescriptor i)
    {

        Real dist_i = m_vector[i].distance;

        for (size_t x = 0; x < m_num_points; x++) {
            Point_info& info_x = m_vector[x];

            if (info_x.exact_distance_used) {
                continue;
            }

            Real new_upper_bound;
            Real new_upper_bound_1 = m_dspanner->m_matrix[x][i].upper_bound;

            if (new_upper_bound_1 == MAX) {
                new_upper_bound = MAX;
            }
            else {
                new_upper_bound = new_upper_bound_1 + dist_i;
            }

            info_x.upper_bound = std::min(info_x.upper_bound, new_upper_bound);

        }

        for (size_t x = 0; x < m_num_points; x++) {
            Point_info& info_x = m_vector[x];

            if (info_x.exact_distance_used) {
                continue;
            }

            Real new_lower_bound;
            Real new_lower_bound_1 = m_dspanner->m_matrix[x][i].upper_bound;

            if (new_lower_bound_1 == MAX) {
                new_lower_bound = 0.0;
            }
            else {
                new_lower_bound = dist_i - new_lower_bound_1;
            }

            info_x.lower_bound = std::max(info_x.lower_bound, new_lower_bound);
        }
    }

    size_t get_num_points() const
    {
        return m_num_points;
    }

    size_t get_number_of_computed_distances() const
    {
        size_t result = 0;
        for (size_t i = 0; i < m_num_points; i++) {
            if (m_vector[i].exact_distance_used) {
                result++;
            }
        }
        return result;
    }

    double get_fraction_of_computed_distances() const
    {

        double num = (double) get_number_of_computed_distances();
        double denom = m_num_points;

        return num / denom;

    }

    void print_bounds()
    {
        for (size_t i = 0; i < m_num_points; i++) {

            Real low = m_vector[i].lower_bound;
            Real upp = m_vector[i].upper_bound;
            std::cout << "[ " << low << " : ";
            if (upp == MAX) {
                std::cout << "inf";
            }
            else {
                std::cout << upp;
            }
            std::cout << " ] ";
        }
        std::cout << std::endl;
    }

};

#endif // WASSERSTEIN_SPACE_POINT_H

