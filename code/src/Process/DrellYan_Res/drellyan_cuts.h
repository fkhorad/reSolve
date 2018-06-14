#ifndef _drellyan_cutsH_
#define _drellyan_cutsH_

#include "phase_space.h"
#include "lorentz.h"
#include "drellyan_input.h"

bool drellyan_cuts(drellyan_input* drellyan_in, double q2, double qt2, double eta, double mures2, PSpoint* PS_);
bool PSDYcuts_1(drellyan_input*, PSpoint*);

#endif
