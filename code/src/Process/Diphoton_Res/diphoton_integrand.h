#ifndef diphoton_integrand_H
#define diphoton_integrand_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "phase_space.h"
#include "constants.h"
#include "events_out.h"

// Generic resummation includes
#include "resu_preproc.h"
#include "inv_fourier.h"
#include "hardfns.h"

// Specific includes
#include "diphoton_input.h"
#include "diphoton_ps.h"
#include "diphoton_cuts.h"
#include "diphoton_hard.h"


class DiphRes{
public:

  diphoton_input* diph_in;
  ResummationInfo* res_1;

};


// Function declarations
int diphoton_integrand(int ndim, const double x[], double& f, void *userdata, double weight, int iter);

int diphoton_integrand_cuba(const int* ndim, const double x[], const int *ncomp,
                            double f[], void *userdata,
                            const int*, const int*, const double* weight, const int* iter);
                       
                       
#endif
