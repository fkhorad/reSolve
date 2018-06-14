#ifndef _resu_PSH_
#define _resu_PSH_

#include "lorentz.h"

struct ResummationInfo;

struct resu_PS {
  double q2, qt2, eta, mur2, muf2,
    mures2, a, b0p, x, alphas, alphaqf;
  int ifit, ifit2;
  std::complex<double> H1q, H1g, H2g;
  std::vector<std::complex<double> > H2q;
  std::vector<std::vector<double> > sigmaij;

  void set(const ResummationInfo&, double q2_in, double eta_in, double qT2_in);

};


void in_state_w_recoil(four_momentum qvec, double sqs, double alpha, four_momentum& k1, four_momentum& k2);


#endif
