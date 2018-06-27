#include <random>
#include <fstream>
#include <iostream>
#include <dirent.h>
#include <sys/stat.h>

#include "boost/filesystem/path.hpp"
#include "boost/filesystem/operations.hpp"

#include "experiment_helper.h"

namespace fs = boost::filesystem;

namespace wasser_spanner {


//bool get_random_diagram_UnifShear(DiagramR& diag,
//                                 const int diagram_size,
//                                 const double lower_bound_hor,
//                                 const double upper_bound_hor,
//                                 const double lower_bound_vert,
//                                 const double upper_bound_vert)
//{
//    if (diagram_size <= 0) {
//        std::cerr << "Error! diagram_size is not positive: " << diagram_size << std::endl;
//        return false;
//    }
//    static std::uniform_real_distribution<double> unif_hor(lower_bound_hor, upper_bound_hor);
//    static std::default_random_engine rand_engineHor;
//    static std::uniform_real_distribution<double> unifVert(lower_bound_vert, upper_bound_vert);
//    static std::default_random_engine rand_engineVert;
//    diag.clear();
//    diag.reserve(diagram_size);
//
//    for(int k = 0; k < diagram_size; ++k) {
//        double rand_x = unif_hor(rand_engineHor);
//        double rand_y = rand_x + unifVert(rand_engineVert);
//        diag.push_back(Point(rand_x, rand_y));
//    }
//    return true;
//}


//bool getRandomDiagramUnif(DiagramR& diag,
//                            const int diagram_size,
//                            const double lower_bound_hor,
//                            const double upper_bound_hor,
//                            const double lower_bound_vert,
//                            const double upper_bound_vert)
//{
//    if (diagram_size <= 0) {
//        std::cerr << "Error! diagram_size is not positive: " << diagram_size << std::endl;
//        return false;
//    }
//    static std::default_random_engine rand_engine;
//    static std::uniform_real_distribution<double> unif_hor(lower_bound_hor, upper_bound_hor);
//    static std::uniform_real_distribution<double> unifVert(lower_bound_vert, upper_bound_vert);
//    //static std::default_random_engine rand_engineVert;
//    diag.clear();
//    diag.reserve(diagram_size);
//
//
//    size_t pointsGenerated { 0 };
//
//    while(pointsGenerated < diagram_size) {
//        double rand_x = unif_hor(rand_engine);
//        double rand_y = unifVert(rand_engine);
//        if ( fabs(rand_x - rand_y) < 0.000001 * ( upper_bound_hor - lower_bound_hor) ) {
//            continue;
//        }
//        pointsGenerated++;
//        if ( rand_x > rand_y ) {
//            double t = rand_x;
//            rand_x = rand_y;
//            rand_y = t;
//        }
//        diag.push_back(Point(rand_x, rand_y));
//    }
//    return true;
//}
//
//bool writeRandomDiagramUnif(const std::string diagram_name,
//                                 const int diagram_size,
//                                 const double lower_bound_hor,
//                                 const double upper_bound_hor,
//                                 const double lower_bound_vert,
//                                 const double upper_bound_vert)
//{
//    DiagramR diag;
//    return getRandomDiagramUnif(diag, diagram_size, lower_bound_hor, upper_bound_hor, lower_bound_vert, upper_bound_vert) and write_diagram(diag, diagram_name);
//}
//
//


    // create random persistence diagram with points concentrated
    // near diagonal and rare persistence features
    bool write_random_diagram_real(const std::string diagram_name,
                                   const int diagram_size,
                                   const double lower_bound_hor,
                                   const double upper_bound_hor,
                                   const double vert_mean,
                                   const double vert_std_dev)
    {
        if (diagram_size <= 0) {
            std::cerr << "Error! diagram_size is not positive: " << diagram_size << std::endl;
            return false;
        }

        static std::uniform_real_distribution<double> unif_hor(lower_bound_hor, upper_bound_hor);
        static std::normal_distribution<double> gauss_vert(vert_mean, vert_std_dev);

        std::ofstream diagram_file(diagram_name);
        if (!diagram_file.good()) {
            std::cerr << "Cannot write to file " << diagram_name << std::endl;
            return false;
        }
        diagram_file.precision(15);

        for (int k = 0; k < diagram_size; ++k) {
            double rand_x = unif_hor(twister);
            double rand_y = rand_x + abs(gauss_vert(twister));
            if (rand_y - rand_x < 0.0001) {
                rand_y += 0.0001 * (upper_bound_hor - lower_bound_hor);
            }
            diagram_file << rand_x << " " << rand_y << std::endl;
        }
        diagram_file.close();
        return true;
    }

    bool write_diagram(const DiagramR& diag, std::string fname)
    {
        std::ofstream dgm_file(fname);
        if (!dgm_file.good()) {
            std::cerr << "Cannot write to file " << fname << std::endl;
            return false;
        }

        for (const auto& p : diag) {
            dgm_file << p.first << " " << p.second << std::endl;
        }
        dgm_file.close();
        return true;
    }

    bool get_random_diagram_normal(DiagramR& diag,
                                   const int diagram_size,
                                   const double lower_bound_hor,
                                   const double upper_bound_hor,
                                   const double std_dev)
    {
        if (diagram_size <= 0) {
            std::cerr << "Error! diagram_size is not positive: " << diagram_size << std::endl;
            return false;
        }

        static std::uniform_real_distribution<double> unif_hor(lower_bound_hor, upper_bound_hor);
        static std::normal_distribution<double> gauss_vert(0.0, std_dev);
        diag.clear();
        diag.reserve(diagram_size);

        for (int k = 0; k < diagram_size; ++k) {
            double a = unif_hor(twister);
            double b = abs(gauss_vert(twister));
            double rand_x = a - 0.5 * b;
            double rand_y = a + 0.5 * b;
            if (rand_y - rand_x < 0.01) {
                rand_y += 0.0001 * (upper_bound_hor - lower_bound_hor);
            }
            diag.emplace_back(rand_x, rand_y);
        }
        return true;
    }


    bool get_seed_diagram(DiagramR& diag,
                          const int diagram_size,
                          const double lower_bound_hor,
                          const double upper_bound_hor,
                          const double std_dev)
    {
        get_random_diagram_normal(diag, diagram_size, lower_bound_hor, upper_bound_hor, std_dev);
    }

    bool perturb_diagram(const DiagramR& diag_in,
                         DiagramR& diag_out,
                         const double std_dev,
                         const double surviving_prob)
    {
        diag_out.clear();
        diag_out.reserve(diag_in.size());
        std::normal_distribution<double> gauss_noise(0.0, std_dev);
        std::uniform_real_distribution<double> survive_distr(0.0, 1.0);
        for (const auto& pt : diag_in) {
            if (survive_distr(twister) <= surviving_prob) {
                double x = pt.first + gauss_noise(twister);
                double y = pt.second + gauss_noise(twister);
                diag_out.emplace_back(x, y);
            }
        }
    }


    bool write_random_diagram_normal(const std::string fname,
                                     const int diagram_size,
                                     const double lower_bound_hor,
                                     const double upper_bound_hor,
                                     const double std_dev)
    {
        DiagramR diag;
        return get_random_diagram_normal(diag, diagram_size, lower_bound_hor, upper_bound_hor, std_dev) and
               write_diagram(diag, fname);
    }


    bool dir_exists(const std::string& dir_name)
    {
        struct stat st = { 0 };
        return (stat(dir_name.c_str(), &st) != -1);
    }

    void create_dir_if_not_exists(const std::string& dirName)
    {
        if (not dir_exists(dirName)) {
            mkdir(dirName.c_str(), 0777);
        }
    }


    void generate_random_diag(DgmVec& dgms, int n_diag, int n_points, int seed)
    {

        std::uniform_real_distribution<double> distr(0.0, 1000.0);

        for (int i = 0; i < n_diag; ++i) {
            WassersteinSpacePointR dgm;
            for (int j = 0; j < n_points; ++j) {
                double x = distr(twister);
                double y = 9000.0 + x + distr(twister);
                dgm.dgm.emplace_back(x, y);
            }
            dgms.push_back(dgm);
        }
    }

// check triangle inequality for all possible combinations of three m_points
// (note that a priori approximate distance computation function is not symmetric)
    bool is_triangle_ineq_satisfied(int a_idx, int b_idx, int c_idx, DynamicSpannerR& ds)
    {
        std::vector<double> d_ab_v { ds.get_distance(a_idx, b_idx), ds.get_distance(b_idx, a_idx) };
        std::vector<double> d_ac_v { ds.get_distance(a_idx, c_idx), ds.get_distance(c_idx, a_idx), };
        std::vector<double> d_bc_v { ds.get_distance(b_idx, c_idx), ds.get_distance(c_idx, b_idx), };
        for (const auto d_ab : d_ab_v) {
            for (const auto d_ac : d_ac_v) {
                for (const auto d_bc : d_bc_v) {
                    if (d_ab > d_ac + d_bc or
                        d_ac > d_ab + d_bc or
                        d_bc > d_ab + d_ac) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool is_triangle_ineq_satisfied(DgmVec& dgms, DynamicSpannerR& ds)
    {
        int n = dgms.size();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = j + 1; k < n; ++k) {
                    std::cout << "checking triple " << i << ", " << j << ", " << ", " << k << std::endl;
                    if (not is_triangle_ineq_satisfied(i, j, k, ds)) {
                        return false;
                    }
                }
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

    void read_mcgill_dgms(std::vector<WassersteinSpacePointR>& dgms, char* dir_name, int dim)
    {
        std::string sdir_name(dir_name);
        std::string sdim = std::to_string(dim) + ".txt";
        DIR* dir;
        struct dirent* ent;
        if ((dir = opendir(dir_name))) {
            int num_subdirs = 0;
            while ((ent = readdir(dir))) {
                std::string sub_dir_name = sdir_name + ent->d_name;
                if (sub_dir_name.find("Ply") != std::string::npos) {
                    num_subdirs++;
                    if (num_subdirs > 6) {
                        break;
                    }
                    std::cout << sub_dir_name << '\n';
                    DIR* sub_dir;
                    struct dirent* sub_ent;
                    if ((sub_dir = opendir(sub_dir_name.c_str()))) {
                        while ((sub_ent = readdir(sub_dir))) {
                            std::string file_name = sub_dir_name + "/" + sub_ent->d_name;
                            if (file_name.find(sdim) != std::string::npos) {
                                dgms.emplace_back(file_name, dgms.size());
                                std::cout << file_name << '\n';
                            }
                        }
                    }
                }
            }
        }
    }

    void save_diagrams(const std::string& dir_name, const DgmVec& dgms)
    {
    }

    // save square matrix to plain text file
    // first line is matrix size
    // every entry of matrix is on separate line
    bool save_distance_matrix(const std::string& fname, const MatrixR& distance_matrix)
    {
        std::ofstream matr_file(fname);
        if (!matr_file.good()) {
            std::cerr << "Cannot write to file " << fname << std::endl;
            return false;
        }

        matr_file << distance_matrix.size() << std::endl;
        matr_file.precision(18);

        for (const auto& v : distance_matrix) {
            for (auto d : v) {
                matr_file << d << std::endl;
            }
        }
        matr_file.close();
        return true;
    }

    bool save_parameters(const std::string& fname, const RandomDiagParams& par)
    {
        std::ofstream param_file(fname);
        if (not param_file.good()) {
            std::cerr << "Cannot write parameters to " << fname;
            return false;
        }
        param_file << "seed: " << par.random_generator_seed << std::endl;
        param_file << "n_diagrams: " << par.n_diagrams << std::endl;
        param_file << "min_diagram_size: " << par.min_diagram_size << std::endl;
        param_file << "max_diagram_size: " << par.max_diagram_size << std::endl;
        param_file << "std_dev: " << par.std_dev << std::endl;
        param_file << "lower_bound_hor: " << par.lower_bound_hor << std::endl;
        param_file << "upper_bound_hor: " << par.upper_bound_hor << std::endl;
        param_file << "lower_bound_vert: " << par.lower_bound_vert << std::endl;
        param_file << "upper_bound_vert: " << par.upper_bound_vert << std::endl;
        param_file << std::endl;
        param_file << "Auction.wasserstein_power: " << par.auction_params.wasserstein_power << std::endl;
        param_file << "Auction.internal_p: " << par.auction_params.internal_p << std::endl;
        param_file << "Auction.delta: " << par.auction_params.delta << std::endl;
        param_file << "Auction.epsilon_common_ratio: " << par.auction_params.epsilon_common_ratio << std::endl;
        param_file << "Auction.initial_epsilon: " << par.auction_params.initial_epsilon << std::endl;
        param_file.close();
        return true;
    }

    void
    get_distance_matrix(const AuctionParamsR& ap, const DgmVec& dgms, MatrixR& distance_matrix, double& max_distance,
                        double& min_distance)
    {
        size_t n_diagrams = dgms.size();
        distance_matrix = MatrixR(n_diagrams, std::vector<double>(n_diagrams, 0.0));

        max_distance = -1.0;
        min_distance = std::numeric_limits<double>::max();
        for (size_t i = 0; i < n_diagrams; ++i) {
            for (size_t j = i + 1; j < n_diagrams; ++j) {
                DiagramR dgm_i = dgms[i].dgm;
                DiagramR dgm_j = dgms[j].dgm;
                double d;
                if (ap.wasserstein_power != std::numeric_limits<double>::infinity() and
                    ap.wasserstein_power != hera::get_infinity()) {
                        d = hera::wasserstein_dist(dgm_i, dgm_j, ap);
                } else {
                        d = hera::bottleneckDistApprox(dgm_i, dgm_j, ap.delta);
                }
                distance_matrix[j][i] = distance_matrix[i][j] = d;
                max_distance = std::max(max_distance, d);
                min_distance = std::min(min_distance, d);
            }
        }
    }

    bool
    read_distance_matrix(MatrixR& distance_matrix, double& max_distance, double& min_distance, const std::string& fname)
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
    }


    void read_diagrams_from_dir(DgmVec& dgms, const std::string& dir_name)
    {
        fs::path dir(dir_name);
        for (int diagram_idx = 0;; ++diagram_idx) {
            std::string fname = "diagram_" + std::to_string(diagram_idx) + ".txt";
            fs::path file(fname);
            fs::path full_path = dir / file;
            if (not fs::exists(full_path)) {
                break;
            }
            dgms.emplace_back(full_path.c_str(), diagram_idx);
        }
    }

    void
    load_test_from_dir(DgmVec& dgms, MatrixR& distance_matrix, double& max_distance, double& min_distance,
                       const std::string& dir_name)
    {
        std::string dist_matr_fname = dir_name + "/distance_matrix.txt";
        read_diagrams_from_dir(dgms, dir_name);
        read_distance_matrix(distance_matrix, max_distance, min_distance, dist_matr_fname);
    }

    void get_random_dgm_vec_low_dim(const RandomDiagParams& par, DgmVec& out_dgms, const std::vector<DiagramR>& seeds,
                                    const std::string& dir_name)
    {
        int hidden_dim = 5;

        // choose random directions
        std::vector<double> displ_angles(hidden_dim, 0.0);
        std::uniform_real_distribution<double> angle_distr(M_PI / 4.0, 5.0 * M_PI / 4.00);
        for(auto& x : displ_angles) {
            x = angle_distr(twister);
        }

        // choose random displacement vectors
        std::vector<std::pair<double, double>> displ_vectors;

        for(int i = 0; i <hidden_dim; ++i) {
            displ_vectors.emplace_back(cos(displ_angles[i]), sin(displ_angles[i]));
        }

        std::uniform_int_distribution<size_t> seed_idx_distr(0, seeds.size());
        for(int i = 0; i  <par.n_diagrams; ++i) {
            out_dgms.emplace_back();
            size_t seed_idx = seed_idx_distr(twister);
            std::vector<double> lambdas;
            for(int j = 0; j < seeds[seed_idx].size(); ++j) {
                double xs = seeds[seed_idx][j].first;
                double ys = seeds[seed_idx][j].second;
                for(int k = 0; k < hidden_dim; ++k) {
                    xs += lambdas[k] * displ_vectors[k].first;
                    ys += lambdas[k] * displ_vectors[k].second;
                }
                out_dgms[i].dgm.emplace_back(xs, ys);
            }
        }
    }

    void get_random_dgm_vec_normal(const RandomDiagParams& par, DgmVec& dgms, const std::string& dir_name)
    {
        dgms.clear();
        dgms.reserve(par.n_diagrams);
        std::uniform_int_distribution<int> diagram_size_distribution(par.min_diagram_size, par.max_diagram_size);

        // init seeds
        for (int i = 0; i < par.n_diagrams; ++i) {
            DiagramR dgm;
            dgm.clear();
            int dgm_size = diagram_size_distribution(twister);
            dgm.reserve(dgm_size);
            get_random_diagram_normal(dgm, dgm_size, par.lower_bound_hor, par.upper_bound_hor, par.std_dev);
            std::string fname = dir_name + "/diagram_" + std::to_string(i) + ".txt";
            write_diagram(dgm, fname);
            dgms.emplace_back(dgm, i);
        }

    }

    void get_random_dgm_vec_clustered(const RandomDiagParams& par, DgmVec& dgms, const std::string& dir_name)
    {
        dgms.clear();

        std::vector<std::vector<DiagramR>> result;

        std::uniform_int_distribution<int> diagram_size_distribution(par.min_diagram_size, par.max_diagram_size);

        // init seeds
        std::vector<DiagramR> seed_dgms;
        seed_dgms.reserve(par.n_seeds);
        for (int i = 0; i < par.n_seeds; ++i) {
            DiagramR seed_dgm;
            seed_dgm.clear();
            int dgm_size = diagram_size_distribution(twister);
            seed_dgm.reserve(dgm_size);
            get_random_diagram_normal(seed_dgm, dgm_size, par.lower_bound_hor, par.upper_bound_hor, par.std_dev);
            seed_dgms.push_back(seed_dgm);
        }

        result.push_back(seed_dgms);

        double curr_std_dev = par.std_dev;
        for (int level = 1; level < par.n_hierarchy_levels; ++level) {
            curr_std_dev *= par.std_dev_decay;
            result.push_back(std::vector<DiagramR>());
            for (const auto& seed_dgm : result[level - 1]) {
                for (int j = 0; j < par.branching_factor; ++j) {
                    DiagramR dgm;
                    perturb_diagram(seed_dgm, dgm, curr_std_dev, par.survive_prob);
                    result[level].push_back(dgm);
                }
            }
        }

        int idx = 0;
        for (const auto& hier_level : result) {
            for (const DiagramR& dgm : hier_level) {
                dgms.emplace_back(dgm, idx);
                std::string fname = dir_name + "/diagram_" + std::to_string(idx) + ".txt";
                write_diagram(dgm, fname);
                idx++;
            }
        }
    }


//    void get_random_dgm_vec_clustered(const RandomDiagParams& par, DgmVec& dgms, const std::string& dir_name)
//    {
//        dgms.clear();
//        dgms.reserve(par.n_diagrams);
//        std::uniform_int_distribution<int> diagram_size_distribution(par.min_diagram_size, par.max_diagram_size);
//
//        // init seeds
//        std::vector<DiagramR> seed_dgms;
//        seed_dgms.reserve(par.n_seeds);
//        for (int i = 0; i < par.n_diagrams; ++i) {
//            DiagramR seed_dgm;
//            seed_dgm.clear();
//            int dgm_size = diagram_size_distribution(twister);
//            seed_dgm.reserve(dgm_size);
//            get_random_diagram_normal(seed_dgm, dgm_size, par.lower_bound_hor, par.upper_bound_hor, par.std_dev);
//            seed_dgms.push_back(seed_dgm);
//        }
//
//        int idx = 0;
//        for(const auto& seed_dgm : seed_dgms) {
//            for(int i = 0; i < par.n_diagrams_pro_seed; ++i) {
//                DiagramR dgm;
//                dgm.clear();
//                int dgm_size = diagram_size_distribution(twister);
//                dgm.reserve(dgm_size);
//                perturb_diagram(seed_dgm, dgm, par.perturbation_std_dev, par.survive_prob);
//                std::string fname = dir_name + "/diagram_" + std::to_string(idx) + ".txt";
//                write_diagram(dgm, fname);
//                dgms.emplace_back(dgm, idx);
//                idx++;
//            }
//        }
//    }

    void get_random_dgm_vec(const RandomDiagParams& par, DgmVec& dgms, const std::string& dir_name)
    {
        switch (par.strategy) {
            case RandomDiagParams::Strategy::NORMAL :
                get_random_dgm_vec_normal(par, dgms, dir_name);
                break;
            case RandomDiagParams::Strategy::CLUSTERED :
                get_random_dgm_vec_clustered(par, dgms, dir_name);
                break;
        }
    }

    void
    setup_random_test(const RandomDiagParams& par, DgmVec& dgms, MatrixR& distance_matrix, double& max_distance,
                      double& min_distance, const std::string& dir_name)
    {
        create_dir_if_not_exists(dir_name);
        get_random_dgm_vec(par, dgms, dir_name);
        get_distance_matrix(par.auction_params, dgms, distance_matrix, max_distance, min_distance);
        std::string dist_matr_fname = dir_name + "/distance_matrix.txt";
        save_distance_matrix(dist_matr_fname, distance_matrix);
        std::string param_fname = dir_name + "/parameters.txt";
        save_parameters(param_fname, par);
    }

    void setup_test_from_dist_matr(MatrixR& distance_matrix, double& max_distance, double& min_distance,
                                   const std::string& dir_name)
    {
        const std::string dist_matr_fname = dir_name + "/distance_matrix.txt";
        read_distance_matrix(distance_matrix, max_distance, min_distance, dist_matr_fname);
    }

    std::ostream& operator<<(std::ostream& os, const RandomDiagParams& params)
    {
        switch (params.strategy) {
            case RandomDiagParams::Strategy::NORMAL :
                os << "Normal distributrion";
                os << "seed: " << params.random_generator_seed << " n_diagrams: " << params.n_diagrams
                   << " diagram_size: "
                   << params.min_diagram_size << " std_dev: " << params.std_dev << " lower_bound_hor: "
                   << params.lower_bound_hor
                   << " upper_bound_hor: " << params.upper_bound_hor << " lower_bound_vert: " << params.lower_bound_vert
                   << " upper_bound_vert: " << params.upper_bound_vert;
                break;
            case RandomDiagParams::Strategy::CLUSTERED :
                os << "Clustered distribution";
                os << "seed: " << params.random_generator_seed << " min_diagram_size: " << params.min_diagram_size
                   << " max_diagram_size: "
                   << params.max_diagram_size << " std_dev: " << params.std_dev << " lower_bound_hor: "
                   << params.lower_bound_hor
                   << " upper_bound_hor: " << params.upper_bound_hor << " lower_bound_vert: " << params.lower_bound_vert
                   << " upper_bound_vert: " << params.upper_bound_vert
                   << " n_seeds: " << params.n_seeds << " branching factor: " << params.branching_factor
                   << " survive_prob: " << params.survive_prob << " std-dev decay: " << params.std_dev_decay;
                break;
        };
        os << " auction_params: " << params.auction_params;
        return os;
    }

    void create_dist_matr_from_dir(const std::string& dir_name, const AuctionParamsR& ap, MatrixR& distance_matrix,
                                   double& max_distance, double& min_distance)
    {
        DgmVec dgms;
        read_diagrams_from_dir(dgms, dir_name);
        get_distance_matrix(ap, dgms, distance_matrix, max_distance, min_distance);
        std::string dist_matr_fname = dir_name + "/distance_matrix.txt";
        save_distance_matrix(dist_matr_fname, distance_matrix);
    }



    void subsample_test(DgmVec& dgms, MatrixR& distance_matrix, double& max_distance, double& min_distance, double fraction)
    {
        std::uniform_real_distribution<double> prob_distr(0.0, 1.0);
        std::unordered_set<size_t> survived_indices;
        for(size_t i = 0; i < dgms.size(); ++i) {
            double p = prob_distr(twister);
            if (p < fraction) {
                survived_indices.insert(i);
            }
        }
        DgmVec new_dgms;
        size_t new_dgms_size = survived_indices.size();
        size_t k = 0;
        for(size_t i = 0; i < dgms.size(); ++i) {
            if (survived_indices.count(i) == 1) {
                new_dgms.emplace_back(dgms[i]);
                new_dgms.back().id = k++;
            }
        }
        size_t i_1 = 0;
        MatrixR new_dist_matr = MatrixR(new_dgms_size, std::vector<double>(new_dgms_size, 0.0));
        max_distance = 0.0;
        min_distance = std::numeric_limits<double>::max();
        for(size_t i = 0; i < dgms.size(); ++i) {
            if (survived_indices.count(i) ==  0)
                continue;
            size_t j_1 = 0;
            for (size_t j = 0; j < dgms.size(); ++j) {
                if (survived_indices.count(j) == 1) {
                    double d = distance_matrix[i][j];
                    new_dist_matr[i_1][j_1] = d;
                    max_distance = std::max(max_distance, d);
                    if (i != j)
                        min_distance = std::min(min_distance, d);
                    ++j_1;
                };
            }
            ++i_1;
        }
        dgms = new_dgms;
        distance_matrix = new_dist_matr;
    }

} // namespace wasser_spanner
