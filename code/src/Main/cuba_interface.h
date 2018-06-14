#ifndef _cuba_interfaceH_
#define _cuba_interfaceH_

#include <string>

#include "cuba.h"
struct InputPars;

// Some auxiliary routines to launch CUBA and save/collect events. I leave them mostly
// uncommented
/*
typedef int (*integrand_t_full)(const int *ndim, const cubareal x[],
    const int *ncomp, cubareal f[], void *userdata,
    const int *, const int *, const double * weight, const int * iter);
*/

void cuba_vegas_call(const InputPars&, int ndim, int ncomp, integrand_t, void* userdata);
void collect_events(std::string, int maxeval, int nstart, int nincrease,
                    std::string machine_tag, int save_events);



#endif
