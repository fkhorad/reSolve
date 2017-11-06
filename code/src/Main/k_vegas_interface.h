#ifndef _k_vegas_interfaceH_
#define _k_vegas_interfaceH_

#include "k_vegas.h"
class InputPars;

void k_vegas_call(InputPars& input_1, k_integrand_t integrand, void* userdata);

void k_vegas_combiner(InputPars& input_1);


#endif
