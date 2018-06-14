#ifndef diphoton_integrand_H
#define diphoton_integrand_H

class diphoton_input;
class ResummationInfo;


// Function declarations
int diphoton_integrand(int ndim, const double x[], double& f, void *userdata, double weight, int iter);

int diphoton_integrand_cuba(const int* ndim, const double x[], const int *ncomp, double f[], void *userdata, const int*, const int*, const double* weight, const int* iter);


#endif
