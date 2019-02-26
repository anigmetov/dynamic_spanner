#include <iostream>

//#define LOG_AUCTION

#include "spdlog/spdlog.h"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>

#include "cover_tree/cover_tree.h"
#include "cover_tree/wspd.h"
#include "dynamic_spanner.h"
#include "nearest_neighbor_search.h"

#include <armadillo>

//using namespace arma;

//std::mt19937_64 twister;

namespace spd = spdlog;

using namespace wasser_spanner;



int get_num_points(const std::string& fname, const char delimiter = ' ')
{
    std::ifstream matr_file(fname);
    if (not matr_file.good()) {
        std::cerr << "Cannot read matrix from file " << fname << std::endl;
        return -1;
    }

    std::string s;
    std::getline(matr_file, s);

    int result = 0;
    bool in_delim = false;
    for(int i = 0; i < s.length(); ++i) {
        if (not in_delim and s[i] == delimiter) {
            result++;
            in_delim = true;
        } else if (in_delim and s[i] != delimiter) {
            in_delim = false;
        }
    }

    matr_file.close();
    return result;
}

int get_dimension(std::string fname, const char delim = ' ')
{
    return get_num_points(fname, delim);
}

bool read_points_file(MatrixR& distance_matrix, double& max_distance, double& min_distance, std::string fname)
{
    std::ifstream points_file(fname);
    if (not points_file.good()) {
        std::cerr << "Cannot read matrix from file " << fname << std::endl;
        return false;
    }

    std::vector<std::vector<double>> points;
    int dim = get_dimension(fname);

    std::string s;
    while(std::getline(points_file, s)) {
        points.emplace_back();
        std::stringstream ss(s);
        double a;
        for (int i = 0; i < dim; ++i) {
            ss >> a;
            points.back().push_back(a);
        }
    }

    int n_points = points.size();
    min_distance = std::numeric_limits<double>::max();
    max_distance = -1.0;

    distance_matrix = MatrixR(n_points, std::vector<double>(n_points, 0.0));

    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < i; ++j) {
            double d = 0.0;
            for (int k = 0; k < dim; ++k) {
                d += (points[i][k] - points[j][k]) * (points[i][k] - points[j][k]);
            }
            distance_matrix[i][j] = distance_matrix[j][i] = sqrt(d);

            max_distance = std::max(distance_matrix[i][j], max_distance);
            if (i != j) {
                min_distance = std::min(distance_matrix[i][j], min_distance);
            }
        }
    }
}

bool read_distance_matrix(MatrixR& distance_matrix, double& max_distance, double& min_distance, const std::string& fname)
{
    std::ifstream matr_file(fname);
    if (not matr_file.good()) {
        std::cerr << "Cannot read matrix from file " << fname << std::endl;
        return false;
    }

    min_distance = std::numeric_limits<double>::max();
    max_distance = -1.0;
    int n_points = get_num_points(fname);

    distance_matrix = MatrixR(n_points, std::vector<double>(n_points, 0.0));

    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < n_points; ++j) {
            matr_file >> distance_matrix[i][j];
            max_distance = std::max(distance_matrix[i][j], max_distance);
            if (i != j) {
                min_distance = std::min(distance_matrix[i][j], min_distance);
            }
        }
    }

    assert(0.0 < min_distance and min_distance <= max_distance);
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = i; j < n_points; ++j) {
            if (i == j)
                assert(distance_matrix[i][j] == 0.0);
            else
                assert(distance_matrix[i][j] == distance_matrix[j][i]);
        }
    }

    return true;
}

bool read_distance_matrix_mcgill(MatrixR& distance_matrix, double& max_distance, double& min_distance, const std::string& fname, std::vector<int>::size_type n_samples)
{
    std::ifstream matr_file(fname);
    if (not matr_file.good()) {
        std::cerr << "Cannot read matrix from file " << fname << std::endl;
        return false;
    }

    int n_points;

    matr_file >> n_points;

    distance_matrix = MatrixR(n_points, std::vector<double>(n_points, 0.0));

    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < n_points; ++j) {
            matr_file >> distance_matrix[i][j];
        }
    }

    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = i; j < n_points; ++j) {
            if (i == j)
                assert(distance_matrix[i][j] == 0.0);
            else
                assert(distance_matrix[i][j] == distance_matrix[j][i]);
        }
    }

    if (n_samples == 0)
        n_samples = n_points;

    std::vector<size_t> indices(n_points);
    std::iota(indices.begin(), indices.end(), 0);
    std::random_device rd;
    std::mt19937_64 twister(rd());
    std::shuffle(indices.begin(), indices.end(), twister);

    indices.resize(n_samples);

    min_distance = std::numeric_limits<double>::max();
    max_distance = -1.0;

    assert(0.0 < min_distance and min_distance <= max_distance);

    MatrixR distance_matrix_final {n_samples, std::vector<double>(n_samples, 0.0)};
    for(size_t i = 0; i < n_samples; ++i) {
        for(size_t j = 0; j < n_samples; ++j) {
            distance_matrix_final[i][j] = distance_matrix[indices[i]][indices[j]];
            max_distance = std::max(distance_matrix_final[i][j], max_distance);
            if (i != j) {
                min_distance = std::min(distance_matrix_final[i][j], min_distance);
            }
        }
    }

    distance_matrix = distance_matrix_final;

    return true;
}


bool read_distance_matrix_and_queries(MatrixR& distance_matrix, MatrixR& queries, double& max_distance, double& min_distance, const std::string& fname)
{
    using Real = std::remove_reference<decltype(min_distance)>::type;
    std::ifstream matr_file(fname);
    if (not matr_file.good()) {
        std::cerr << "Cannot read matrix from file " << fname << std::endl;
        return false;
    }

    min_distance = std::numeric_limits<Real>::max();
    max_distance = -1.0;
    size_t n_diagrams;
    matr_file >> n_diagrams;

    distance_matrix = MatrixR(n_diagrams, std::vector<Real>(n_diagrams, 0.0));

    for (size_t i = 0; i < n_diagrams; ++i) {
        for (size_t j = 0; j < n_diagrams; ++j) {
            matr_file >> distance_matrix[i][j];
            max_distance = std::max(distance_matrix[i][j], max_distance);
            if (i != j) {
                min_distance = std::min(distance_matrix[i][j], min_distance);
            }
        }
    }

    assert(0.0 < min_distance and min_distance <= max_distance);
    for (size_t i = 0; i < n_diagrams; ++i) {
        for (size_t j = i; j < n_diagrams; ++j) {
            if (i == j)
                assert(distance_matrix[i][j] == 0.0);
            else
                assert(distance_matrix[i][j] == distance_matrix[j][i]);
        }
    }

    size_t no_queries;

    matr_file >> no_queries;

    queries = MatrixR(no_queries, std::vector<Real>(n_diagrams, 0.0));

    for (size_t i = 0; i < no_queries; i++) {
        for (size_t j = 0; j < n_diagrams; j++) {
            matr_file >> queries[i][j];
        }
    }

    return true;
}

double get_expansion_constant(const MatrixR& dist_matrix, double base)
{
    assert(base > 1.0);
    double result = base;
    for (auto v : dist_matrix) {
        std::sort(v.begin(), v.end());
        assert(v[0] == 0.0);
        for (size_t i = 1; i < v.size(); ++i) {
            if (i + 1 < v.size() and v[i] == v[i + 1]) {
                continue;
            }
            double r = v[i];
            size_t n_points_in_ring = 0;
            for (size_t j = i + 1; j < v.size() and v[j] < base * r; ++j)
                ++n_points_in_ring;
            result = std::max(result, (i + n_points_in_ring + 1.0) / (i + 1.0));
        }
    }
    return result;
}

//arma::mat generate_uniform_points(int dim, int n_points)
//{
//    return arma::randu(n_points, dim);
//}

//MatrixR get_distance_matrix(const arma::mat& points)
//{
//    int n_points = points.n_rows;
////    arma::mat dist_matr = arma::zeros(n_points, n_points);
//    MatrixR result (n_points, std::vector<double>(n_points, 0.0));
//    for(int i = 0; i < n_points; ++i) {
//        for (int j = i + 1; j < n_points; ++j) {
//            double d = arma::norm(points.row(i) - points.row(j), 2);
//            result[i][j] = d;
//            result[j][i] = d;
////            dist_matr(i, j) = d;
////            dist_matr(j, i) = d;
//        }
//    }
//    return result;
//}

int main(int argc, char** argv)
{
    auto console = spd::stdout_color_mt("console");
    console->set_level(spd::level::info);

#if 0
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " DIM NPOINTS EPSILON" << std::endl;
        return 0;
    }
    int dim = std::atoi(argv[1]);
    int n_points = std::atoi(argv[2]);
    double eps = std::atof(argv[3]);
    if (dim < 1) {
        std::cerr << "Bad dimension" << std::endl;
        return 1;
    }

    if (n_points < 1) {
        std::cerr << "Bad number of points" << std::endl;
        return 1;
    }

    if (eps <= 0.0) {
        std::cerr << "Bad epsilon" << std::endl;
        return 1;
    }
    auto points = generate_uniform_points(dim, n_points);
    auto dist_matrix = get_distance_matrix(points);
    DynamicSpannerR spanner(dist_matrix);
    spanner.construct_blind_random_eps_spanner(eps);
    std::cout << "sparseness: " << spanner.get_fraction_of_computed_distances() << std::endl;
#else
    MatrixR dist_matrix, queries;
    double max_dist = -1.0;
    double min_dist = std::numeric_limits<double>::max();

    if (argc < 4) {
        std::cout << "Usage: " << argv[0]
                  << " input_name epsilon is_point [n_samples] [logname]"
                  << std::endl;
        return 0;
    }

    int arg_idx = 1;
    std::string dist_name = argv[arg_idx++];

    double eps = atof(argv[arg_idx++]);
    console->info("Reading from file {}, epsilon = {}", argv[1], argv[2]);


    std::string is_input_dm = argv[arg_idx++];
    bool is_input_dist_matrix = is_input_dm == "y" or is_input_dm == "Y" or is_input_dm == "d" or is_input_dm == "D";

    std::vector<int>::size_type n_samples = (argc > arg_idx) ? atoi(argv[arg_idx++]) : 0;

    std::string log_name  = (argc > arg_idx) ? argv[arg_idx++] : "experiment_log.txt";
    console->info("log_name = {}, epsilon = {}, n_samples = {}", log_name, eps, n_samples);
    auto exp_logger = spd::basic_logger_st("experiment_logger", log_name);
    exp_logger->set_pattern("%v");

    if (is_input_dist_matrix) {
        console->info("Reading distance matrix");
        read_distance_matrix_mcgill(dist_matrix, max_dist, min_dist, dist_name, n_samples);
    } else {
        console->info("Reading points");
        read_points_file(dist_matrix, max_dist, min_dist, dist_name);
    }

    size_t n = dist_matrix.size();

    console->info("dist_matrix size = {}, spread = {} / {} = {}", dist_matrix.size(), max_dist, min_dist, max_dist / min_dist);

    DynamicSpannerR spanner(dist_matrix);
    DynamicSpannerR ct_spanner(dist_matrix);
    std::cout << "Spanner initialized" << std::endl;

#if 0

    DynamicSpannerR copy(spanner);
//
//    copy.construct_blind_greedy_eps_spanner(eps);
//
//    console->info("eps Spanner built from scratch. Requested distances: {}, computed distances: {}\n", copy.get_fraction_of_requested_distances(),
//            copy.get_fraction_of_computed_distances());


    ct_spanner.make_exact_distance_keeper();

    CoverTree ct(max_dist, ct_spanner);

    // TODO very ugly
    DynamicSpannerR* ds = &ct.m_dspanner;

    console->info("Cover tree built. Requested distance: {}, fraction = {} / {} = {}",
            ds->get_number_of_requested_distances(),
            ds->get_number_of_requested_distances(),
            ds->m_num_points * (ds->m_num_points  -1 )/ 2,
            ds->get_fraction_of_requested_distances());

//    std::cout << "Cover tree built. Requested distance : " << ds->get_fraction_of_requested_distances() << ", " << ds->get_number_of_requested_distances() << " out of " <<  << std::endl;
//              << ", computed distances : " << ds->get_fraction_of_computed_distances() << std::endl << std::endl;

    copy = *ds;

//    copy.construct_blind_greedy_eps_spanner(eps);
//
//    std::cout << "eps Spanner built from cover tree. Requested distance : " << copy.get_fraction_of_requested_distances()
//              << ", computed distances : " << copy.get_fraction_of_computed_distances() << std::endl << std::endl;

    WspdNode::dspanner = &ct.m_dspanner;
    WSPD wspd(ct, eps);

    console->info("WSPD built. Requested distance: {}, fraction = {}", ds->get_number_of_requested_distances(), ds->get_fraction_of_requested_distances());
//    std::cout << "WSPD built. Requested distance : " << ds->get_fraction_of_requested_distances() << ", " << ds->get_number_of_requested_distances() << std::endl;
//              << ", computed distances : " << ds->get_fraction_of_computed_distances() << std::endl << std::endl;

//    copy = *ds;
//
//    copy.construct_blind_greedy_eps_spanner(eps);
//
//    std::cout << "eps Spanner built from wspd. Requested distance : " << copy.get_fraction_of_requested_distances()
//              << ", computed distances : " << copy.get_fraction_of_computed_distances() << std::endl << std::endl;

    //ds->print_ratios();
    wspd.make_spanner1();
    console->info("WSPD spanner built. Requested distances = {}, fraction = {}, #eddges in spanner = {}", ds->get_number_of_requested_distances(), ds->get_fraction_of_computed_distances(), wspd.m_num_spanner_edges);
//    std::cout << "eps-Spanner built with WSPD method. Requested distance : " << ds->get_fraction_of_requested_distances() << ", " << ds->get_number_of_requested_distances() << std::endl;
//              << ", computed distances : " << ds->get_fraction_of_computed_distances() << std::endl << std::endl;

    std::cout << "Requested distance / # points " << (double) (ds->get_number_of_requested_distances())  / ds->m_num_points << std::endl;

    exp_logger->info("{};wspd;{};{};{};{}", dist_name, eps, ct.m_num_points, wspd.m_num_spanner_edges, ds->get_number_of_requested_distances());
//    console->log("ratio: {}",
    //ds->print_ratios();
    //std::cout << std::endl << std::endl;

    /*
    for(size_t i=0;i<n;i++) {
      for(size_t j=0;j<n;j++) {
    if(i==j) {
      std::cout << "1 ";
    } else {
      std::cout << wspd.get_spanner_distance(i,j)/ds->m_matrix[i][j].distance << " ";
    }
      }
      std::cout << std::endl;
    }
    */

    //double exp_const = get_expansion_constant(dist_matrix);
    //std::cout << "expansion_constant: " << exp_const << std::endl;

    bool is_valid_tree = ct.is_valid_tree();
    std::cout << "Is valid tree:" << is_valid_tree << std::endl;

//    bool is_valid_wspd = wspd.is_valid();
//
//    std::cout << "WSPD is valid: " << is_valid_wspd << std::endl;

#else

    {
        DynamicSpannerR greedy_spanner(dist_matrix);
        greedy_spanner.construct_greedy_eps_spanner(eps);
        console->info("{};{}", dist_name, greedy_spanner.get_statistics());
        exp_logger->info("{};{}", dist_name, greedy_spanner.get_statistics());
        console->flush();
        exp_logger->flush();
    }

    {
        DynamicSpannerR quasi_greedy_spanner(dist_matrix);
        quasi_greedy_spanner.construct_blind_quasi_sorted_greedy_eps_spanner(eps);
        console->info("{};{}", dist_name, quasi_greedy_spanner.get_statistics());
        exp_logger->info("{};{}", dist_name, quasi_greedy_spanner.get_statistics());
    }

    {
        DynamicSpannerR quasi_shaker_spanner(dist_matrix);
        quasi_shaker_spanner.construct_blind_quasi_sorted_shaker_eps_spanner(eps);
        exp_logger->info("{};{}", dist_name, quasi_shaker_spanner.get_statistics());
    }

    {
        DynamicSpannerR blind_random_spanner_1_1(dist_matrix);
        blind_random_spanner_1_1.construct_blind_random_ratio_eps_spanner(true, true, eps);
        console->info("{};{}", dist_name, blind_random_spanner_1_1.get_statistics());
        exp_logger->info("{};{}", dist_name, blind_random_spanner_1_1.get_statistics());
        console->flush();
        exp_logger->flush();
    }

    {
        DynamicSpannerR blind_random_spanner_0_0(dist_matrix);
        blind_random_spanner_0_0.construct_blind_random_ratio_eps_spanner(false, false, eps);
        console->info("{};{}", dist_name, blind_random_spanner_0_0.get_statistics());
        exp_logger->info("{};{}", dist_name, blind_random_spanner_0_0.get_statistics());
        console->flush();
        exp_logger->flush();
    }

    {
        DynamicSpannerR blind_random_spanner_1_0(dist_matrix);
        blind_random_spanner_1_0.construct_blind_random_ratio_eps_spanner(true, false, eps);
        console->info("{};{}", dist_name, blind_random_spanner_1_0.get_statistics());
        exp_logger->info("{};{}", dist_name, blind_random_spanner_1_0.get_statistics());
    }

    {
        DynamicSpannerR blind_random_spanner_0_1(dist_matrix);
        blind_random_spanner_0_1.construct_blind_random_ratio_eps_spanner(false, true, eps);
        console->info("{};{}", dist_name, blind_random_spanner_0_1.get_statistics());
        exp_logger->info("{};{}", dist_name, blind_random_spanner_0_1.get_statistics());
    }

    {
        DynamicSpannerR blind_greedy_spanner(dist_matrix);
        blind_greedy_spanner.construct_blind_greedy_eps_spanner(eps);
        console->info("{};{}", dist_name, blind_greedy_spanner.get_statistics());
        exp_logger->info("{};{}", dist_name, blind_greedy_spanner.get_statistics());
        exp_logger->flush();
    }

    {
        DynamicSpannerR blind_random_spanner(dist_matrix);
        blind_random_spanner.construct_blind_random_eps_spanner(eps);
        console->info("{};{}", dist_name, blind_random_spanner.get_statistics());
        exp_logger->info("{};{}", dist_name, blind_random_spanner.get_statistics());
        exp_logger->flush();
    }

//    for (size_t i = 0; i < queries.size(); i++) {
//        std::vector<double>& query = queries[i];
//        Nearest_neighbor_search<double> nns(query, &spanner);
//        size_t result = nns.find_nearest_neighbor();
//        std::cout << "Query " << i << ", answer: " << result << " - Computed distances: "
//                  << nns.get_fraction_of_computed_distances() << std::endl;
//        std::cout << "(correct answer is " << std::distance(query.begin(), std::min_element(query.begin(), query.end()))
//                  << ")" << std::endl << std::endl;
//    }
#endif

#endif
    return 0;
}
