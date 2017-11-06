#ifndef _diphoton_cutsH_
#define _diphoton_cutsH_

#include "phase_space.h"
#include "lorentz.h"
#include "diphoton_input.h"

bool diph_cuts(diphoton_input* diph_in, double q2, double qt2, double eta, double mures2, PSpoint* PS_);
bool PScuts_1(diphoton_input*, PSpoint*);

#endif
