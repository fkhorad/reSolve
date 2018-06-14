//Integrands for drellyan case
#include "drellyan_integrand.h"

#include <fstream>
#include <sstream>
#include <iostream>

#include "phase_space.h"
#include "events_out.h"

// Generic resummation includes
#include "inv_fourier.h"
#include "resu_PS.h"
#include "resummation_input.h"

// Specific includes
#include "drellyan_ps.h"
#include "drellyan_cuts.h"
#include "drellyan_hard.h"
#include "drellyan_input.h"


int drellyan_integrand(int ndim, const double x[], double& f, void *userdata, double weight, int iter) {

    drellyan_input* drellyan_in = (drellyan_input*) userdata;
    ResummationInfo* resuminfo = &drellyan_in->res_1;

//    double pi = k_constants::pi;

// Phase space
    PSpoint pp1;
    double randsjacob;
    resu_PS resuvars;

    randsjacob = drellyan_ps(x, drellyan_in, pp1, resuvars);


// CUTS:
    bool cutcheck = false; //i.e not cut
    cutcheck = drellyan_cuts(drellyan_in, resuvars.q2, resuvars.qt2, resuvars.eta, resuvars.mures2, &pp1);

    double tempf;
    if(cutcheck==true || randsjacob==0.){
      tempf = 0.;
    }
    else{

      int Nf = resuminfo->Nf;
      sigmaijdrellyancalc(drellyan_in, pp1, resuvars.sigmaij, resuvars.alphas, resuvars.q2);

//All hard factors are 0 for DY as in DY scheme - everything for virtual corrections is in the C factors
      double H1q = 0.0;
      double H1g = 0.0;
      double H2q[Nf];
      for (int j = 0; j < Nf; j++) {
        H2q[j] = 0.0;
      }
      double H2g = 0.0;
      resuvars.H1q = std::complex<double>(H1q,0.0);
      resuvars.H1g = std::complex<double>(H1g,0.0);
      resuvars.H2g = std::complex<double>(H2g,0.0);
      for (int i = 0; i<Nf; i++) {
        resuvars.H2q[i] = H2q[i];
      }

// Call main resummation function
      tempf = resummed(resuvars, resuminfo);
    }

    f = tempf*randsjacob;

// Print event to file

    if(drellyan_in->event_info.save_events==1 || drellyan_in->event_info.save_events==2){
      events_out(drellyan_in->event_info, iter, f, weight, std::sqrt(resuvars.mur2), resuminfo->alpha_QED, resuvars.alphas, pp1, ndim, x);
    }

    return 0;
}



int drellyan_integrand_cuba(const int* ndim, const double x[], const int *ncomp,
                            double f[], void *userdata,
                            const int* unused1, const int* unused2, const double* weight, const int* iter){

  double fres;
  int outcode;

  outcode = drellyan_integrand(*ndim, x, fres, userdata, *weight, *iter);

  if(*ncomp > 1) outcode = 1;

  f[0] = fres;

  return outcode;

}
