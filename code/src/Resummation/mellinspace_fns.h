#ifndef MELLINSPACE_FNS_H
#define MELLINSPACE_FNS_H

// Routine to do Inverse Mellin transform integrals for N-dependent resummation functions which are process independent

//Set up the contour and weights - inverse mellin transform is done over a contour in the complex plane with a postive and negative branch from -infinity to +infinity

#include <complex>
#include <vector>
#include <array>

#include "constants.h"


class contour_info {
 public:
  std::vector<double> weights;
  std::vector<std::complex<double> > Np, Nm;
  std::complex<double> phicomplexp, phicomplexm;
};

class resuDep_N {
 public:
  std::complex<double> qqi, qgf, gqi, ggi, ggf, ns1mi, ns1pi, ns1f,
    qq1f, qg1f, gq1i, gq1f, gg1i, gg1f;
  std::complex<double> C1qg, C1gq, C1qq, C1gg;
  std::complex<double> C2qg, C2NSqqb, C2Sqqb, C2NSqq;
  std::complex<double> gamma1qq, gamma1qg, gamma1gq, gamma1gg;
  std::complex<double> gamma2qq, gamma2qqV, gamma2qqbV, gamma2qqS, gamma2qqbS, gamma2qg, gamma2gq, gamma2gg;

};


class C2values {
 public:
  std::complex<double> xp, C2qgpos, C2NSqqbpos, C2Sqqbpos, C2NSqqpos;
  std::complex<double> xm, C2qgneg, C2NSqqbneg, C2Sqqbneg, C2NSqqneg;
};

void inv_mel_init(contour_info& contours, std::array<resuDep_N, k_constants::mellin_points>& PosBranch, std::array<resuDep_N, k_constants::mellin_points>& NegBranch, int verbosity, int Nf, double Ca, double Cf);

void C2valuescalc(C2values& b, std::complex<double> xp, std::complex<double> xm,
                  int verbosity, int Nf, double Ca, double Cf);


//////////////////////////////////////////////////////////////
// VARIOUS FUNCTIONS IN MELLIN SPACE -- ADAPTED FROM BLUEMLEIN
//////////////////////////////////////////////////////////////

std::complex<double> psi(std::complex<double> z, int verbosity);
std::complex<double> psideriv1(std::complex<double> z, int verbosity);
std::complex<double> psideriv2(std::complex<double> z, int verbosity);
std::vector<std::complex<double> > anomcalc(std::complex<double> x, int verbosity);
std::complex<double> *anomcalc( std::complex<double> ancalc[15], std::complex<double> x , int verbosity);
std::complex<double> *C2calc( std::complex<double> c2arrcalc[5], std::complex<double> x , int verbosity);
std::complex<double> psi0 (std::complex<double> x, int verbosity);
std::complex<double> psi1 (std::complex<double> x, int verbosity);
std::complex<double> psi2 (std::complex<double> x, int verbosity);
std::complex<double> psi3 (std::complex<double> x, int verbosity);
std::complex<double> bet (std::complex<double> x, int verbosity);
std::complex<double> bet1 (std::complex<double> x, int verbosity);
std::complex<double> bet2 (std::complex<double> x, int verbosity);
std::complex<double> bet3 (std::complex<double> x, int verbosity);
std::complex<double> acg1 (std::complex<double> x, int verbosity);
std::complex<double> acg1p (std::complex<double> x, int verbosity);
std::complex<double> acg1pp (std::complex<double> x, int verbosity);
std::complex<double> acg2 (std::complex<double> x, int verbosity);
std::complex<double> acg2p (std::complex<double> x, int verbosity);
std::complex<double> acg3 (std::complex<double> x, int verbosity);
std::complex<double> acg4 (std::complex<double> x, int verbosity);
std::complex<double> acg4p (std::complex<double> x, int verbosity);
std::complex<double> acg5 (std::complex<double> x, int verbosity);
std::complex<double> acg6 (std::complex<double> x, int verbosity);
std::complex<double> acg7 (std::complex<double> x, int verbosity);
std::complex<double> acg8 (std::complex<double> x, int verbosity);
std::complex<double> acg9 (std::complex<double> x, int verbosity);
std::complex<double> acg13 (std::complex<double> x, int verbosity);
std::complex<double> acg20 (std::complex<double> x, int verbosity);
std::complex<double> acg21 (std::complex<double> x, int verbosity);
std::complex<double> cbeta(std::complex<double> z1, std::complex<double> z2, int verbosity);
std::complex<double> lngam(std::complex<double> z, int verbosity);


#endif
