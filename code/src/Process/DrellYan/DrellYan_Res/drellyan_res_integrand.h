#ifndef drellyan_res_integrand_H
#define drellyan_res_integrand_H

struct resu_PS;
struct drellyan_input;
class PSpoint;


// Function declarations
int drellyan_res_integrand(int ndim, const double x[], double& f, drellyan_input* drellyan_in, double weight, int iter);

double drellyan_res_calc(const double x[], drellyan_input* drellyan_in, resu_PS& resuvars, PSpoint& pp1, double randsjacob);

#endif
