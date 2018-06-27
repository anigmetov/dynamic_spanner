#include <iostream>

//#define LOG_AUCTION

#include "spdlog/spdlog.h"

#include "wasserstein/wasserstein_space_point.h"
#include "cover_tree/cover_tree.h"
#include "cover_tree/wspd.h"
#include "experiment_helper.h"

std::mt19937_64 twister;

namespace spd = spdlog;
using namespace wasser_spanner;

int main(int argc, char** argv)
{
    auto console = spd::stdout_color_mt("console");
    console->set_level(spd::level::info);

    DgmVec dgms;
    MatrixR dist_matrix;
    RandomDiagParams par;
    double max_dist = -1.0;
    double min_dist = std::numeric_limits<double>::max();

//    par.auction_params.wasserstein_power = atof(argv[2]);
//    create_dist_matr_from_dir(argv[1], par.auction_params, dist_matrix, max_dist, min_dist);
//    return 0;

    if (false) {

        if (argc < 10) {
            std::cout << "Usage: " << argv[0]
                      << " dir_name wasser_power number_of_diagrams min_diagram_size max_diagram_size max_x std_dev std_dev_decay WSPD_epsilon [seed]"
                      << std::endl;
            return 0;
        }

        int arg_idx = 1;
        std::string dir_name = argv[arg_idx++];
        create_dir_if_not_exists(dir_name);

        par.auction_params.wasserstein_power = atof(argv[arg_idx++]);
        par.n_diagrams = atoi(argv[arg_idx++]);
        par.min_diagram_size = atoi(argv[arg_idx++]);
        par.max_diagram_size = atoi(argv[arg_idx++]);
        par.upper_bound_hor = atof(argv[arg_idx++]);
        par.std_dev = atof(argv[arg_idx++]);
        par.std_dev_decay = atof(argv[arg_idx++]);
        double eps = atof(argv[arg_idx++]);
        par.random_generator_seed = (argc > arg_idx) ? atoi(argv[arg_idx++]) : 1;

        par.strategy = RandomDiagParams::Strategy::CLUSTERED;

        twister.seed(par.random_generator_seed);

    }

    if (argc < 5) {
        std::cout << "Usage: " << argv[0]
                  << " dir_name wasser_power WSPD_epsilon subsample_prob"
                  << std::endl;
        return 0;
    }

    int arg_idx = 1;
    std::string dir_name = argv[arg_idx++];

    par.auction_params.wasserstein_power = atof(argv[arg_idx++]);
    double eps = atof(argv[arg_idx++]);
    double subsample_factor = atof(argv[arg_idx++]);

    std::string log_fname = dir_name + "/experiment_log.txt";
    auto file_log = spd::basic_logger_mt("experiment_logger", log_fname);
    file_log->set_level(spd::level::info);
    file_log->set_pattern("[%H:%M:%S.%e] %v");
    file_log->flush_on(spd::level::info);
    file_log->info("Started running. Parameters: {}", par);


    // setup_random_test(par, dgms, dist_matrix, max_dist, min_dist, dir_name);
    load_test_from_dir(dgms, dist_matrix, max_dist, min_dist, dir_name);

    twister.seed(1);
    subsample_test(dgms, dist_matrix, max_dist, min_dist, subsample_factor);

    par.n_diagrams = dist_matrix.size();

    int n_dist_pairs = dgms.size() * (dgms.size() - 1) / 2;

    console->info("dist_matrix size = {}, dgms.size = {}", dist_matrix.size(), dgms.size());

    file_log->info("Distances calculated. max_dist / min_dist =  {} / {} = {} ",
            max_dist, min_dist, max_dist / min_dist);


    CoverTree ct(max_dist, dgms, dist_matrix, par.auction_params);
    // ct.print_tree(std::cout);
    DynamicSpannerR& ds = ct.m_dspanner;

    console->info("Cover tree built. Requested distance : {}, computed distances : {}, cached calls: {}, total distances: {}, geom_lb_useful = {}",
                  ds.m_requested_distances.size() / 2, ds.m_distance_cache.size() / 2, ds.m_cached_calls, n_dist_pairs, ds.m_geom_lower_bound_useful);
    file_log->info("Cover tree built. Requested distance = {}, computed distances =  {}, cached calls = {}, total distances = {}, geom_lb_usefule = {}",
                   ds.m_requested_distances.size() / 2, ds.m_distance_cache.size() / 2, ds.m_cached_calls, n_dist_pairs, ds.m_geom_lower_bound_useful);

    WspdNode::dspanner = &ct.m_dspanner;
    WSPD wspd(ct, eps);

    console->info(
            "WSPD built, size = {}, requested distance : {}, computed distances : {}, cached calls: {}, total distances: {}",
            wspd.size(), ds.m_requested_distances.size() / 2, ds.m_distance_cache.size() / 2,
            ds.m_cached_calls, n_dist_pairs);
    file_log->info(
            "WSPD built, size = {}, requested distance : {}, computed distances : {}, cached calls: {}, total distances: {}",
            wspd.size(), ds.m_requested_distances.size() / 2, ds.m_distance_cache.size() / 2,
            ds.m_cached_calls, n_dist_pairs);

    file_log->info("Started estimating spanner quality");
    console->info("Started estimating spanner quality");
    double max_rel_error = wspd.get_max_relative_error();
    file_log->info("Max. relative error: {}", max_rel_error);
    console->info("Max. relative error: {}", max_rel_error);

    file_log->info("Started estimating dynamic spanner quality");
    std::vector<double> deltas { 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0 };
    auto spanner_estimate = ct.m_dspanner.get_bounds_quality(deltas);
    for(const auto x : spanner_estimate) {
        file_log->info("Spanner estimates with relative error at most {}: {}", x.first, x.second);
        console->info("Spanner estimates with relative error at most {}: {}", x.first, x.second);
    }

    console->info("Started calculating expansion constant");
    file_log->info("Started calculating expansion constant");
    double exp_const = get_expansion_constant(dist_matrix);
    console->info("expansion constant = {}", exp_const);
    file_log->info("expansion constant = {}", exp_const);



    // validation of data structures

    bool is_valid_tree = ct.is_valid_tree();
    console->info("Is valid tree: {}", is_valid_tree);
    file_log->info("Is valid tree: {}", is_valid_tree);

    bool is_valid_wspd = wspd.is_valid();

    console->info("WSPD is valid: {}", is_valid_wspd);
    file_log->info("WSPD is valid: {}", is_valid_wspd);

    return 0;
}
