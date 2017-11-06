#ifndef _kvegasH_
#define _kvegasH_

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <ctime>
#include <sstream>


//#include <random> // or other random generator... This is for quick initial testing,
// is experimental, and requires -std=c++11 compiler option, which may lead to trouble


#define def_alpha 1.5
#define def_grid_dim 50


typedef int (*k_integrand_t)(int ndim, const double x[],
			     double& f, void *userdata, double weight, int iter);

typedef std::vector<std::vector<double> > grid_t;


class vegas_integral{

private:

  double alpha;
  int int_dim, grid_dim, last_neval, last_iter;

  grid_t grid, dd, di;

  double si,si2,swgt,schi;
  double last_ti, last_tsi, last_sdi;

public:

  vegas_integral();
  vegas_integral(int dim);

  int integrate(k_integrand_t func, int neval, int verbose, void* userdata);
  void refine_grid();
  void iteration_output(std::ostream & out1 = std::cout);

  void reset();
  void reset_resize(int int_dim_new, int grid_dim_new = def_grid_dim);
  void reset_cumulatives();
  void set_alpha(double);
  void copy(vegas_integral&);
  int copy_from_file(std::istream & infile);
  
  void combine_partial_iters(int nparts, vegas_integral* iter_parts);

  int niter;
  
// read only
  int r_dim();
  int r_gdim();
  double r_grid(int,int);
  double r_si();
  double r_si2();
  double r_swgt();
  double r_schi();
  int r_neval();
  int r_iter();
  void info_dump(std::ostream & out1 = std::cout);
  
  
};



// Random machinery

//  extern std::default_random_engine generator;
//  extern std::uniform_real_distribution<double> distribution;

  void vegas_random_init(int);
  void vegas_random(int, double*);


  extern "C" {
    void randa_(int*,double*);
  }

  extern "C" {
    void randin_(int*);
  }

#endif
