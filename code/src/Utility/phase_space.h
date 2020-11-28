#ifndef _phase_spaceH_
#define _phase_spaceH_

#include <vector>
#include <complex>

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
  void add_mom(four_momentum);

// Get
  int dim() const;
  four_momentum mom(int) const;
  double ss(int,int) const;
  std::complex<double> sa(int,int) const;
  std::complex<double> sb(int,int) const;


private:

  std::vector<four_momentum> momenta;
  std::vector<std::vector<double> > DotProducts;
  std::vector<std::vector<std::complex<double> > > SpinProductsA;
  std::vector<std::vector<std::complex<double> > > SpinProductsB;

// For Debugging
  void PS_checker(int);

};

void set_PS_fromfile(const char* filename, PSpoint&);
/* void set_PS_fromfile_jets(const char* filename, PSpoint&); */

////////////////////////////////////////////////////////////////

void set_PS_twobody(double costheta, double phi, const four_momentum P_in, double m1, double m2, PSpoint& PS);
//
void set_PS_threebody(double cos1, double phi1, double mQ, double cos2, double phi2, const four_momentum P_in, double m1, double m2, double m3, PSpoint& PS);


#endif
