#ifndef HARDFNS_H
#define HARDFNS_H

#include <complex>
#include <vector>

class ResummationInfo;
class resummationIndepfns;

#include "mellinspace_fns.h"


class PSdep_variables {
 public:
  double q2, qt2, eta, mur2, muf2,
    mures2, a, b0p, x, alphas, alphaqf;
  int ifit, ifit2;
  std::complex<double> H1q, H1g, H2g;
  std::complex<double> H2q[5];
  double sigmaij[11][11];
};


class AllDep_N {
public:
//  std::complex<double> c1qg, c1gq, c1qq, c1gg, c2qgM, c2NSqqM, c2SqqbM, c2NSqqbM;
//  std::complex<double> gamma1qq, gamma1qg, gamma1gq, gamma1gg, gamma2qq, gamma2qqV, gamma2qqbV, gamma2qqS, gamma2qqbS, gamma2qg, gamma2gq, gamma2gg;
  std::vector<std::complex<double> > FP;

void preent();

};

void GetResuPars (int i, int isign, int ibeam, ResummationInfo* resuminfo,
                  PSdep_variables* resu, std::complex<double> alpq,
                  AllDep_N& ParamsforResu);


//Calculates Hard functions for inverse mellin transform argument
void Hardfns_calc (int order,
                   PSdep_variables* resu, resummationIndepfns* indep,
                   AllDep_N PR1, AllDep_N PR2, resuDep_N N1, resuDep_N N2,
                   std::complex<double> alpq,
                   std::complex<double> aexp, std::complex<double> aexpb,
                   std::complex<double> sudakq, std::complex<double> sudakg,
                   std::complex<double>& HCRNqq, std::complex<double>& HCRNgg);



#endif
