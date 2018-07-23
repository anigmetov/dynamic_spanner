#ifndef WASSERSTEIN_SPACE_POINT_H
#define WASSERSTEIN_SPACE_POINT_H

#include <iostream>
#include <utility>
#include <vector>
#include <limits>
#include <algorithm>

#include <cassert>

//#include "spdlog/spdlog.h"

//namespace spd = spdlog;

constexpr int INVALID_ID = -1;

template<class Real_ = double>
class DynamicSpanner
{

  public:
  
  using Real = Real_;

  protected:

  struct Pair_of_points_info {
    
    Real distance;
    bool distance_requested;
    bool exact_distance_used;
    Real upper_bound;
    Real lower_bound;

    Pair_of_points_info() {}

    Pair_of_points_info(Real distance) 
      : distance(distance) {
      distance_requested = false;
      exact_distance_used = false;
      upper_bound = std::numeric_limits<Real>::max();
      lower_bound = 0.0;
    }
  };

public:

  using Matrix = std::vector<std::vector<Pair_of_points_info>>;

  using MatrixReal = typename std::vector<std::vector<Real>>;
  
  using VertexDescriptor = size_t;

  size_t m_num_points;
  Matrix m_matrix;

public:

  DynamicSpanner(const MatrixReal& distance_matrix) {
    
    m_num_points = distance_matrix.size();
    m_matrix.resize(m_num_points);
    
    for(size_t i=0;i<m_num_points;i++) {
      m_matrix[i].resize(m_num_points);
      
      for(size_t j=0;j<m_num_points;j++) {
	m_matrix[i][j] = Pair_of_points_info(distance_matrix[i][j]);
      }
    }

    // On the diagonal, bounds are exact
    for(size_t i=0;i<m_num_points;i++) {
      m_matrix[i][i].lower_bound = m_matrix[i][i].upper_bound = 0.0;
      m_matrix[i][i].exact_distance_used=true; // not sure that is ever queried
    }
  }
  
  
  Real get_distance_no_cache(VertexDescriptor i, VertexDescriptor j) const
  {
    return m_matrix[i][j].distance;
    
  }

  // TODO; use to avoid case distinctions?
  struct Comparator {

    bool m_strict;

    Comparator(bool strict) : m_strict(strict) {}

    bool operator() (Real x, Real y) {
      if(m_strict) {
	return x<y;
      }
      return x<=y;
    }
  };

  bool is_distance_less(VertexDescriptor i, VertexDescriptor j,
                          Real value, bool strict)
  {

    //std::cout << "Call less than " << i <<  " " << j << ", with dist " << m_matrix[i][j].distance << " val=" << value << " bounds: {" << m_matrix[i][j].lower_bound << "," << m_matrix[i][j].upper_bound << "]" << std::endl;

    if (i == j) {
      if (strict) {
	return 0.0 < value;
      } else {
	return 0.0 <= value;
      }
    }

    Pair_of_points_info& info = m_matrix[i][j];
    Pair_of_points_info& info_mirror = m_matrix[j][i];
        
    info.distance_requested=true;
    info_mirror.distance_requested=true;    

    if(strict && info.upper_bound < value) {
      return true;
    }

    if(not strict && info.upper_bound <= value) {
      return true;
    }

    if(strict && info.lower_bound >= value) {
      return false;
    }

    if(not strict && info.lower_bound > value) {
      return false;
    }

    Real dist = get_distance(i,j);
    
    if(strict) {
      return dist<value;
    } else {
      return dist<=value;
    }
    
  }

  void update_bounds_using_distance(VertexDescriptor i,
				   VertexDescriptor j) {

    Real dist_ij = m_matrix[i][j].distance;

    for(size_t x=0;x<m_num_points;x++) {
      for(size_t y=x+1;y<m_num_points;y++) {
	Pair_of_points_info& info_xy = m_matrix[x][y];
	Pair_of_points_info& info_yx = m_matrix[y][x];
	if(info_xy.exact_distance_used) {
	  continue;
	}

	Real new_upper_bound_1;
	Real new_upper_bound_1_1 = m_matrix[x][i].upper_bound;
	Real new_upper_bound_1_2 = m_matrix[j][y].upper_bound;
	if(new_upper_bound_1_1==std::numeric_limits<Real>::max()
	   or 
	   new_upper_bound_1_2==std::numeric_limits<Real>::max()) {
	  new_upper_bound_1 = std::numeric_limits<Real>::max();
	} else {
	  new_upper_bound_1 = new_upper_bound_1_1 + dist_ij+new_upper_bound_1_2;
	}
	Real new_upper_bound_2;
	Real new_upper_bound_2_1 = m_matrix[x][j].upper_bound;
	Real new_upper_bound_2_2 = m_matrix[i][y].upper_bound;
	if(new_upper_bound_2_1==std::numeric_limits<Real>::max()
	   or 
	   new_upper_bound_2_2==std::numeric_limits<Real>::max()) {
	  new_upper_bound_2 = std::numeric_limits<Real>::max();
	} else {
	  new_upper_bound_2 = new_upper_bound_2_1 + dist_ij+new_upper_bound_2_2;
	}
	Real new_upper_bound = std::min(new_upper_bound_1,new_upper_bound_2);
		
	info_xy.upper_bound = std::min(info_xy.upper_bound, new_upper_bound);
	info_yx.upper_bound = info_xy.upper_bound;
      }
    }

    for(size_t x=0;x<m_num_points;x++) {
      for(size_t y=x+1;y<m_num_points;y++) {
	Pair_of_points_info& info_xy = m_matrix[x][y];
	Pair_of_points_info& info_yx = m_matrix[y][x];
	if(info_xy.exact_distance_used) {
	  continue;
	}

	Real new_lower_bound_1;
	Real new_lower_bound_1_1 = m_matrix[x][i].upper_bound;
	Real new_lower_bound_1_2 = m_matrix[j][y].upper_bound;
	//std::cout << x << " " << y << "Upper bounds xi " <<  new_lower_bound_1_1 << " jy " << new_lower_bound_1_2 << std::endl;
	if(new_lower_bound_1_1==std::numeric_limits<Real>::max()
	   or 
	   new_lower_bound_1_2==std::numeric_limits<Real>::max()) {
	  new_lower_bound_1 = 0.0;
	} else {
	  new_lower_bound_1 = dist_ij-new_lower_bound_1_1-new_lower_bound_1_2;
	}
	Real new_lower_bound_2;
	Real new_lower_bound_2_1 = m_matrix[x][j].upper_bound;
	Real new_lower_bound_2_2 = m_matrix[i][y].upper_bound;
	//std::cout << "Upper bounds xj " <<  new_lower_bound_2_1 << " iy " << new_lower_bound_2_2 << std::endl;
	if(new_lower_bound_2_1==std::numeric_limits<Real>::max()
	   or 
	   new_lower_bound_2_2==std::numeric_limits<Real>::max()) {
	  new_lower_bound_2 = 0.0;
	} else {
	  new_lower_bound_2 = dist_ij-new_lower_bound_2_1-new_lower_bound_2_2;
	}
	Real new_lower_bound = std::max(new_lower_bound_1,new_lower_bound_2);
		
	info_xy.lower_bound = std::max(info_xy.lower_bound, new_lower_bound);
	info_yx.lower_bound = info_xy.lower_bound;
      }
    }

    #if 1
    {
      for(size_t x=0;x<m_num_points;x++) {
	for(size_t y=x+1;y<m_num_points;y++) {
	  Pair_of_points_info& info_xy = m_matrix[x][y];
	  Pair_of_points_info& info_yx = m_matrix[y][x];
	  if(info_xy.exact_distance_used) {
	    continue;
	  }
	  
	  Real new_lower_bound_1;
	  Real new_lower_bound_1_1 = m_matrix[x][i].upper_bound;
	  Real new_lower_bound_1_2 = m_matrix[j][y].lower_bound;
	  //std::cout << x << " " << y << "Upper bounds xi " <<  new_lower_bound_1_1 << " jy " << new_lower_bound_1_2 << std::endl;
	  if(new_lower_bound_1_1==std::numeric_limits<Real>::max()
	     or 
	     new_lower_bound_1_2==0.0) {
	    new_lower_bound_1 = 0.0;
	  } else {
	    new_lower_bound_1 = new_lower_bound_1_2-dist_ij-new_lower_bound_1_1;
	  }
	  Real new_lower_bound_2;
	  Real new_lower_bound_2_1 = m_matrix[x][j].upper_bound;
	  Real new_lower_bound_2_2 = m_matrix[i][y].lower_bound;
	  //std::cout << "Upper bounds xj " <<  new_lower_bound_2_1 << " iy " << new_lower_bound_2_2 << std::endl;
	  if(new_lower_bound_2_1==std::numeric_limits<Real>::max()
	     or 
	     new_lower_bound_2_2== 0.0) {
	    new_lower_bound_2 = 0.0;
	  } else {
	    new_lower_bound_2 = new_lower_bound_2_2-dist_ij-new_lower_bound_2_1;
	  }
	  Real new_lower_bound_3;
	  Real new_lower_bound_3_1 = m_matrix[j][y].upper_bound;
	  Real new_lower_bound_3_2 = m_matrix[x][i].lower_bound;
	  //std::cout << x << " " << y << "Upper bounds xi " <<  new_lower_bound_1_1 << " jy " << new_lower_bound_1_2 << std::endl;
	  if(new_lower_bound_3_1==std::numeric_limits<Real>::max()
	     or 
	     new_lower_bound_3_2==0.0) {
	    new_lower_bound_3 = 0.0;
	  } else {
	    new_lower_bound_3 = new_lower_bound_3_2-dist_ij-new_lower_bound_3_1;
	  }
	  Real new_lower_bound_4;
	  Real new_lower_bound_4_1 = m_matrix[i][y].upper_bound;
	  Real new_lower_bound_4_2 = m_matrix[x][j].lower_bound;
	  //std::cout << "Upper bounds xj " <<  new_lower_bound_2_1 << " iy " << new_lower_bound_2_2 << std::endl;
	  if(new_lower_bound_4_1==std::numeric_limits<Real>::max()
	     or 
	     new_lower_bound_4_2== 0.0) {
	    new_lower_bound_4 = 0.0;
	  } else {
	    new_lower_bound_4 = new_lower_bound_4_2-dist_ij-new_lower_bound_4_1;
	  }

	  Real new_lower_bound = std::max(std::max(new_lower_bound_1,new_lower_bound_2),
					  std::max(new_lower_bound_3,new_lower_bound_4));
	  
	  info_xy.lower_bound = std::max(info_xy.lower_bound, new_lower_bound);
	  info_yx.lower_bound = info_xy.lower_bound;
	}
      }
    }
     
    #endif
   
#if DEBUG
    for(size_t x=0;x<m_num_points;x++) {
	for(size_t y=0;y<m_num_points;y++) {
	  Pair_of_points_info& info_xy = m_matrix[x][y];
	  Pair_of_points_info& info_yx = m_matrix[y][x];
	  assert(info_xy.upper_bound>=info_xy.distance);
	  assert(info_xy.lower_bound<=info_xy.distance);
	  assert(info_xy.upper_bound==info_yx.upper_bound);
	  assert(info_xy.lower_bound==info_yx.lower_bound);
	}
    }
#endif

  }

	   

  bool is_distance_greater(VertexDescriptor i, VertexDescriptor j,
			   Real value, bool strict)
  {
    return not is_distance_less(i, j, value, not strict);
  }
  

  Real get_distance(VertexDescriptor i, VertexDescriptor j) {
    //console->debug("get_distance called for {} {}", i, j);
    Pair_of_points_info& info = m_matrix[i][j];
    Pair_of_points_info& info_mirror = m_matrix[j][i];
    info.exact_distance_used = true;
    info_mirror.exact_distance_used = true;
    info.upper_bound=info.lower_bound=info.distance;
    info_mirror.upper_bound=info_mirror.lower_bound=info.distance;
    update_bounds_using_distance(i,j);
    return info.distance;
  }

  Real get_distance_no_cache(VertexDescriptor i, VertexDescriptor j) {
    return m_matrix[i][j].distance;
  }

  size_t get_num_points() const {
    return m_num_points;
  }


  size_t get_number_of_computed_distances() const {
    size_t result=0;
    for(size_t i=0;i<m_num_points;i++) {
      for(size_t j=i+1;j<m_num_points;j++) {
	if(m_matrix[i][j].exact_distance_used) {
	  result++;
	}
      }
    }
    return result;
  }

  size_t get_number_of_requested_distances() const {
    size_t result=0;
    for(size_t i=0;i<m_num_points;i++) {
      for(size_t j=i+1;j<m_num_points;j++) {
	if(m_matrix[i][j].distance_requested) {
	  result++;
	}
      }
    }
    return result;
  }

  double get_fraction_of_computed_distances() const {
    
    double num=(double)get_number_of_computed_distances();
    double denom = m_num_points*(m_num_points-1)/2.0;

    return num/denom;

  }

  double get_fraction_of_requested_distances() const {
    
    double num=(double)get_number_of_requested_distances();
    double denom = m_num_points*(m_num_points-1)/2.0;

    return num/denom;

  }

  void find_worst_ratio(size_t& x, size_t& y, double& ratio) {
    
    std::vector<std::pair<size_t,size_t> > candidates;

    double result=0.0;
    for(size_t i=0;i<m_num_points;i++) {
      for(size_t j=i+1;j<m_num_points;j++) {
	Real low = m_matrix[i][j].lower_bound;
	Real upp = m_matrix[i][j].upper_bound;
	if(low==0.0 or upp==std::numeric_limits<Real>::max()) {
	  if(result != std::numeric_limits<Real>::max()) {
	    candidates.clear();
	  }
	  result = std::numeric_limits<Real>::max();
	  candidates.push_back(std::make_pair(i,j));
	    
	} else {
	  double curr = upp/low;
	  if(curr>result) {
	    candidates.clear();
	  }
	  if(curr>=result) {
	    candidates.push_back(std::make_pair(i,j));
	    result = curr;
	  }
	}
      }
    }
    std::random_shuffle(candidates.begin(),candidates.end());
    x = candidates.begin()->first;
    y = candidates.begin()->second;
    ratio = result;
    //std::cout << "Worst ratio at " << x << ", " << y << ": " << m_matrix[x][y].lower_bound << " " << m_matrix[x][y].upper_bound << std::endl;
  }

  void construct_blind_greedy_eps_spanner(double eps) {
    size_t i,j;
    double ratio;

    find_worst_ratio(i,j,ratio);

    while(ratio>(1+eps)) {
      get_distance(i,j);
      find_worst_ratio(i,j,ratio);
    }

  }
  
  struct Value_with_index {
    double v;
    size_t i;
    size_t j;
    Value_with_index(double d, size_t i, size_t j) : v(d), i(i), j(j) {}
  };

  struct Comp_value_with_index {
    
    bool operator() (Value_with_index& a, Value_with_index& b) {
      return a.v<b.v;
    }
    Comp_value_with_index() {}
  };

  void construct_greedy_eps_spanner(double eps) {
    
    std::vector<Value_with_index> vec;
    
    for(size_t i=0;i<m_num_points;i++) {
      for(size_t j=i+1;j<m_num_points;j++) {
	vec.push_back(Value_with_index(m_matrix[i][j].distance,i,j));
      }
    }
    
    std::sort(vec.begin(),vec.end(),Comp_value_with_index());
    //std::reverse(vec.begin(),vec.end());

    //std::random_shuffle(vec.begin(),vec.end());

    for(size_t k=0;k<vec.size();k++) {
      Pair_of_points_info& info = m_matrix[vec[k].i][vec[k].j];
      
      //std::cout << "New pair " << vec[k].i << " " << vec[k].j << " with dist " << info.distance << " bounds: " << info.lower_bound << " " << info.upper_bound << " ratio: " << info.upper_bound / info.distance << std::endl;

      if(info.upper_bound/info.distance>(1+eps)) {
	get_distance(vec[k].i,vec[k].j);
      }
    }
  }
    


  void print_ratios() {
    for(size_t i=0;i<m_num_points;i++) {
      for(size_t j=0;j<m_num_points;j++) {
	if(i==j) {
	  std::cout <<"0 ";
	} else {
	  Real low = m_matrix[i][j].lower_bound;
	  Real upp = m_matrix[i][j].upper_bound;
	  if(low==0.0 or upp==std::numeric_limits<Real>::max()) {
	    std::cout << "inf ";
	  } else {
	    std::cout << upp/low << " ";
	  }
	}
      }
      std::cout << std::endl;
    }
  }


};

#endif // WASSERSTEIN_SPACE_POINT_H

