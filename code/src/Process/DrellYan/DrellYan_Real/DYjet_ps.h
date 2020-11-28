#ifndef _DYjet_psH_
#define _DYjet_psH_


#include "lorentz.h"

struct resu_PS;
struct drellyan_input;
class drellyan_amplitude;


//DYRes phase space method:

double DYjet_ps(const double x[], drellyan_input* drellyan_in, drellyan_amplitude&, resu_PS&, double &ss_hat, double &etaa_hat);

void in_state_plus_jet(four_momentum qvec, double CM_energy, double kjet_plus, four_momentum& k1, four_momentum& k2, four_momentum& kjet);


#endif
