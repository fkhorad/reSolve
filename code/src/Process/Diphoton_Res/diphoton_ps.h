#ifndef _diphoton_psH_
#define _diphoton_psH_


#include "lorentz.h"

class PSpoint;
struct resu_PS;
struct ResummationInfo;


extern "C" {
  double alphas_(double *MUR);
}


double diphoton_ps(const double x[], ResummationInfo*, PSpoint&, resu_PS&);


#endif
