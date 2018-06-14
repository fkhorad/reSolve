//Integrands for diphoton case

#include "diphoton_integrand.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "phase_space.h"
#include "events_out.h"

// Generic resummation includes
#include "resummation_input.h"
#include "inv_fourier.h"
#include "resu_PS.h"

// Specific includes
#include "diphoton_input.h"
#include "diphoton_ps.h"
#include "diphoton_cuts.h"
#include "diphoton_hard.h"


// Debugging
extern int k_count00, k_count0;
int k_count00 = 0;


int diphoton_integrand(int ndim, const double x[], double& f, void* userdata, double weight, int iter) {

    diphoton_input* diph_in = (diphoton_input*) userdata;
    ResummationInfo* resuminfo = &diph_in->res_1;

// Phase space
    PSpoint pp1;
    resu_PS resuvars;

    double randsjacob = diphoton_ps(x, resuminfo, pp1, resuvars);

    
// CUTS:
    bool cutcheck = false; //i.e not cut
    cutcheck = diph_cuts(diph_in, resuvars.q2, resuvars.qt2, resuvars.eta, resuvars.mures2, &pp1);


    double tempf;
    if(cutcheck==true){
      tempf = 0.;
    }
    else{

// Calculate the process-dependent stuff for resummation
      int Nf = resuminfo->Nf;
      sigmaijdiphotoncalc(diph_in, pp1, resuvars.sigmaij, resuvars.alphas);

      // for (int i = 0; i < 11; i++) {
      // 	for (int j = 0; j < 11; j++) {
      // 	  if (resuvars.sigmaij[i][j] != 0) {
      // 	    std::cout << "sigmaij[" << i << "][" << j << " ]= " << resuvars.sigmaij[i][j] << std::endl;
      // 	  }
      // 	}
      // }
      // std::cout << "resuvars.sigmaij[6][4] = " << resuvars.sigmaij[6][4] << std::endl;

      double s = 2.*pp1.ss(0,1);
      double t = -2.*pp1.ss(0,2);
      double costheta = 1. + 2.*t/s;

      double H1q = 2*H1qdiphoton_DY(costheta, resuminfo->verbosity);
      double H1g = 0.0;
      double H2q[Nf];
      for (int j = 0; j < Nf; j++) {
        H2q[j] = 4*H2qdiphotoncalc_DY(&pp1, costheta, j+1, resuminfo->Cf, resuminfo->Nf, resuminfo->Nc, resuminfo->verbosity);
      }
      double H2g = 0.0;
      resuvars.H1q = std::complex<double>(H1q,0.0);
      resuvars.H1g = std::complex<double>(H1g,0.0);
      resuvars.H2g = std::complex<double>(H2g,0.0);
      for (int i = 0; i<Nf; i++) {
        resuvars.H2q[i] = H2q[i];
      }

// Call main resummation function

      tempf = resummed(resuvars, resuminfo)*std::sqrt(1-pow(costheta,2)); //Accounting for the sintheta jacobian factor (it's not sinthetaCM it's sintheta so have to include here not in diphoton_ps.cc)
      // std::cout << "tempf pre cuts = " << tempf << std::endl;
    }

    f = tempf*randsjacob;


    // std::cout << "q2 = " << resuvars.q2 << std::endl;
    // std::cout << "qt2 = " << resuvars.qt2 << std::endl;
    // std::cout << "eta = " << resuvars.eta << std::endl;
    
    // std::cout << "tempf = " << tempf << std::endl;
    // std::cout << "randsjacob = " << randsjacob << std::endl;
    // std::cout << "f = " << f << std::endl;
//

// Print event to file

    if(diph_in->event_info.save_events==1 || diph_in->event_info.save_events==2){
      events_out(diph_in->event_info, iter, f, weight, std::sqrt(resuvars.mur2), resuminfo->alpha_QED, resuvars.alphas, pp1, ndim, x);
    }

    return 0;
}



int diphoton_integrand_cuba(const int* ndim, const double x[], const int *ncomp, double f[], void *userdata, const int* unused1, const int* unused2, const double* weight, const int* iter){

  double fres;
  int outcode;

  outcode = diphoton_integrand(*ndim, x, fres, userdata, *weight, *iter);

  if(*ncomp > 1) outcode = 1;

  f[0] = fres;

  return outcode;

}
