#include <iostream>

//#define LOG_AUCTION

//#include "spdlog/spdlog.h"

#include <fstream>
#include <iostream>

#include "cover_tree/cover_tree.h"
#include "cover_tree/wspd.h"
#include "dynamic_spanner.h"

//std::mt19937_64 twister;

//namespace spd = spdlog;
using namespace wasser_spanner;

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


int main(int argc, char** argv)
{
  //auto console = spd::stdout_color_mt("console");
  //console->set_level(spd::level::info);

    MatrixR dist_matrix;
    double max_dist = -1.0;
    double min_dist = std::numeric_limits<double>::max();

    if (argc < 3) {
        std::cout << "Usage: " << argv[0]
                  << " input_name epsilon"
                  << std::endl;
        return 0;
    }

    int arg_idx = 1;
    std::string dist_name = argv[arg_idx++];

    double eps = atof(argv[arg_idx++]);

    std::cout << "eps=" << eps << std::endl;

    read_distance_matrix(dist_matrix, max_dist, min_dist, dist_name);

    size_t n = dist_matrix.size();
    
    std::cout << "Distance Marix of size " << n << " computed; max distance = " << max_dist << std::endl;

    //console->info("dist_matrix size = {}", dist_matrix.size());
    
    /*
    file_log->info("Distances calculated. max_dist / min_dist =  {} / {} = {} ",
            max_dist, min_dist, max_dist / min_dist);
    */

    DynamicSpannerR spanner(dist_matrix);

    std::cout << "Spanner initialized" << std::endl << std::endl;

#if 1

    DynamicSpannerR copy(spanner);

    copy.construct_greedy_eps_spanner(eps);

    std::cout << "eps Spanner built from scratch. Requested distance : " << copy.get_fraction_of_requested_distances()
	      << ", computed distances : " << copy.get_fraction_of_computed_distances() << std::endl << std::endl;
    
    CoverTree ct(max_dist, spanner);

    //std::cout << "cover tree constructed" << std::endl;
    // ct.print_tree(std::cout);

    // TODO very ugly
    DynamicSpannerR* ds = &ct.m_dspanner;

    std::cout << "Cover tree built. Requested distance : " << ds->get_fraction_of_requested_distances()
	      << ", computed distances : " << ds->get_fraction_of_computed_distances() << std::endl << std::endl;



    copy = *ds;

    copy.construct_greedy_eps_spanner(eps);

    std::cout << "eps Spanner built from cover tree. Requested distance : " << copy.get_fraction_of_requested_distances()
	      << ", computed distances : " << copy.get_fraction_of_computed_distances() << std::endl << std::endl;


    WspdNode::dspanner = &ct.m_dspanner;
    WSPD wspd(ct, eps);

    std::cout << "WSPD built. Requested distance : " << ds->get_fraction_of_requested_distances()
	      << ", computed distances : " << ds->get_fraction_of_computed_distances() << std::endl << std::endl;

    copy = *ds;

    copy.construct_greedy_eps_spanner(eps);

    std::cout << "eps Spanner built from wspd. Requested distance : " << copy.get_fraction_of_requested_distances()
	      << ", computed distances : " << copy.get_fraction_of_computed_distances() << std::endl << std::endl;

    

    //ds->print_ratios();
    wspd.make_spanner();
    std::cout << "eps-Spanner built with WSPD method. Requested distance : " << ds->get_fraction_of_requested_distances()
	      << ", computed distances : " << ds->get_fraction_of_computed_distances() << std::endl << std::endl;

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

    bool is_valid_wspd = wspd.is_valid();

    std::cout << "WSPD is valid: " << is_valid_wspd << std::endl;

#else

    spanner.construct_greedy_eps_spanner(eps);

    std::cout << "Spanner built. Requested distance : " << spanner.get_fraction_of_requested_distances()
	      << ", computed distances : " << spanner.get_fraction_of_computed_distances() << std::endl;

    std::cout << std::endl << std::endl;
    spanner.print_ratios();

    

#endif

    return 0;
}
