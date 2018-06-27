#ifndef SPANNER_WASSERSTEIN_EXPERIMENT_HELPER_H
#define SPANNER_WASSERSTEIN_EXPERIMENT_HELPER_H

#include <vector>
#include <ostream>

#include "wasserstein/wasserstein_space_point.h"

extern std::mt19937_64 twister;

namespace wasser_spanner {

    using DynamicSpannerR = DynamicSpanner<double>;
    using WassersteinSpacePointR = DynamicSpannerR::WassersteinSpacePoint;
    using DiagramPointR = DynamicSpannerR::DiagramPointR;
    using DiagramR = DynamicSpannerR::DiagramR;

    using DgmVec = std::vector<WassersteinSpacePointR>;
    using AuctionParamsR = hera::AuctionParams<double>;
    using MatrixR= DynamicSpannerR::MatrixR;

    struct RandomDiagParams
    {
        enum class Strategy { NORMAL, CLUSTERED };
        friend std::ostream& operator<<(std::ostream& os, const RandomDiagParams& params);

        Strategy strategy { Strategy::CLUSTERED };
        int random_generator_seed { 1 };
        int n_diagrams { 10 }; // for NORMAL only
        int min_diagram_size { 10 };
        int max_diagram_size { 1000 };
        double std_dev { 100.0 };
        double lower_bound_hor { 0.0 };
        double upper_bound_hor { 100.0 };
        double lower_bound_vert { 0.0 };
        double upper_bound_vert { 100.0 };

        // for clustering
        int n_seeds { 30 };
        double survive_prob { 1.0 };
        int n_hierarchy_levels { 3 };
        int branching_factor{ 3 };
        double std_dev_decay { 0.1 };

        AuctionParamsR auction_params;
    };


    void generate_random_diag(DgmVec& dgms, int n_diag, int n_points, int seed = 1);

    // check triangle inequality for all possible combinations of three m_points
    // (note that a priori approximate distance computation function is not symmetric)
    bool is_triangle_ineq_satisfied(int a_idx, int b_idx, int c_idx, DynamicSpannerR& ds);

    bool is_triangle_ineq_satisfied(DgmVec& dgms, DynamicSpannerR& ds);

    double get_expansion_constant(const std::vector<std::vector<double>>& dist_matrix, double base = 2.0);

    std::vector<std::vector<double>> get_dist_matrix(DgmVec& dgms, DynamicSpannerR& ds);

    void read_mcgill_dgms(std::vector<WassersteinSpacePointR>& dgms, char* dir_name, int dim);

    void setup_random_test(const RandomDiagParams& par, DgmVec& dgms, MatrixR& distance_matrix, double& max_distance,
                           double& min_distance, const std::string& dir_name);

    void create_dir_if_not_exists(const std::string& dirName);

    void create_dist_matr_from_dir(const std::string& dir_name, const AuctionParamsR& ap, MatrixR& distance_matrix,
                                   double& max_distance, double& min_distance);

    void load_test_from_dir(DgmVec& dgms, MatrixR& distance_matrix, double& max_distance, double& min_distance,
                            const std::string& dir_name);

    void subsample_test(DgmVec& dgms, MatrixR& distance_matrix, double& max_distance, double& min_distance, double fraction);


    void setup_test_from_dist_matr(MatrixR& distance_matrix, double& max_distance, double& min_distance,
                                   const std::string& dir_name);

}

#endif //SPANNER_WASSERSTEIN_EXPERIMENT_HELPER_H
