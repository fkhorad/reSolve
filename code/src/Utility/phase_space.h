#ifndef _phase_spaceH_
#define _phase_spaceH_

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#include "constants.h"
#include "lorentz.h"

// Phase Space (PS) point class + some generic routines using it.
// For the purposes of this code, "PS point" means "a collection of 4-momenta"

// This is the most ill-developed and prone to possible changes part of the code
// for now, expecially the routines that define the PS in particular cases (2->2,
// 2->3, ...)


class PSpoint{

public:

// Set
  void set_products();
  void set_dim(int);
  int set_mom(int, double*);
  int set_mom(int, four_momentum);

// Get
  int dim();
  four_momentum mom(int);
  double ss(int,int);
  std::complex<double> sa(int,int);
  std::complex<double> sb(int,int);


private:

  std::vector<four_momentum> momenta;
  std::vector<std::vector<double> > DotProducts;
  std::vector<std::vector<std::complex<double> > > SpinProductsA;
  std::vector<std::vector<std::complex<double> > > SpinProductsB;

// For Debugging
  void PS_checker(int);

};


// Various basic PS generation routines: NOT memeber functions of the PS class!

  double set_PS_twobody(const double* rands, const double* masses, PSpoint&);
  double set_PS_twobody(const double* rands, const double* masses, double, double, PSpoint&);
  double set_PS_threebody(const double* rands, const double* masses, PSpoint&);
  double set_PS_threebody(const double* rands, const double* masses, double, double, PSpoint&);
  void set_PS_fromfile(const char* filename, PSpoint&);
  double nonQCD2to2PS(const double* x, double en, PSpoint&);
  double nonQCD2to3PS(const double* x, double en, PSpoint&);
  double nonQCD2to3PS(const double* x, double en, double, double, PSpoint&);
  double QCD_0_2to3PS(const double* x, double en, PSpoint&);
  double QCD_0_2to2PS(const double* x, double en, PSpoint&);
  double QCD_1_2to3PS(const double* x, double en, double, double, PSpoint&);


#endif
