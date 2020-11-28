#ifndef drellyan_integrand_H
#define drellyan_integrand_H


// Function declarations
int drellyan_integrand(int ndim, const double x[], double& f, void *userdata, double weight, int iter);

int drellyan_integrand_cuba(const int* ndim, const double x[], const int *ncomp, double f[], void *userdata, const int*, const int*, const double* weight, const int* iter);

#endif
