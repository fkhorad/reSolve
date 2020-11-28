//Integrands for drellyan case
#include "drellyan_integrand.h"

#include <fstream>
#include <sstream>
#include <iostream>

#include "drellyan_input.h"
#include "drellyan_res_integrand.h"
#include "DYjet_integrand.h"


int drellyan_integrand(int ndim, const double x[], double& f, void *userdata, double weight, int iter){

    drellyan_input* drellyan_in = (drellyan_input*) userdata;

// Born cases uses same PS routine as resummation
    if(drellyan_in->res_1.order == 0 || drellyan_in->res_1.resum_flag == 1 || drellyan_in->res_1.CT_flag == 1){
      drellyan_res_integrand(ndim, x, f, drellyan_in, weight, iter);
    }
    else{
      DYjet_integrand(ndim, x, f, drellyan_in, weight, iter);
    }

    return 0;

}


int drellyan_integrand_cuba(const int* ndim, const double x[], const int *ncomp, double f[], void *userdata, const int* unused1, const int* unused2, const double* weight, const int* iter){

  double fres;
  int outcode;

  outcode = drellyan_integrand(*ndim, x, fres, userdata, *weight, *iter);

  if(*ncomp > 1) outcode = 1;

  f[0] = fres;

  return outcode;

}
