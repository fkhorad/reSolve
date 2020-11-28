#ifndef _drellyan_res_psH_
#define _drellyan_res_psH_


#include "lorentz.h"

class PSpoint;
struct resu_PS;
struct ResummationInfo;
struct drellyan_input;


extern "C" {
// for DYres comparison
  void setqmass_(double*, double*);
  double alphas_ellis_(double* mu, double* amz, int* nloop);
}

//DYRes phase space method:

double drellyan_res_ps(const double x[], drellyan_input* drellyan_in, PSpoint&, resu_PS&);

#endif
