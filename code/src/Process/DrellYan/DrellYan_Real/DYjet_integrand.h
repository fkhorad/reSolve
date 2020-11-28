#ifndef DYjet_integrand_H
#define DYjet_integrand_H

struct resu_PS;
struct drellyan_input;
class drellyan_amplitude;


// Function declarations
int DYjet_integrand(int ndim, const double x[], double& f, drellyan_input* drellyan_in, double weight, int iter);

double DYjet_calc(const double x[], drellyan_input* drellyan_in, resu_PS& resuvars, drellyan_amplitude& pp1, double randsjacob, double ss_hat, double etaa_hat);

#endif
