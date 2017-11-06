#ifndef _diphoton_psH_
#define _diphoton_psH_


#include "lorentz.h"

class PSpoint;
class PSdep_variables;
class ResummationInfo;


extern "C" {
  double alphas_(double *MUR);
}


double diphoton_ps(const double x[], double CM_energy, ResummationInfo*,
		   double& randsjacob, PSpoint&, PSdep_variables&);

double kinematics(double M2, double qT2, double eta, double thetaCM,
                  double phiCM, double sqs, four_momentum& k1,
                  four_momentum& k2, four_momentum& qq1, four_momentum& qq2,
                  int verbosity);

#endif
