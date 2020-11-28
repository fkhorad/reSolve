//Integrands for drellyan case
#include "drellyan_res_integrand.h"

#include <fstream>
#include <sstream>
#include <iostream>

#include "phase_space.h"
#include "events_out.h"

// Generic resummation includes
#include "resummation_input.h"
#include "inv_fourier.h"
#include "resu_PS.h"
#include "CT_NLO.h"

// Specific includes
#include "drellyan_res_ps.h"
#include "drellyan_res_hard.h"
#include "drellyan_input.h"
#include "drellyan_cuts.h"

#include "xsection.h"


int drellyan_res_integrand(int ndim, const double x[], double& f, drellyan_input* drellyan_in, double weight, int iter){

    ResummationInfo* resuminfo = &drellyan_in->res_1;

// Phase space
    PSpoint pp1;
    double randsjacob;
    resu_PS resuvars;

    randsjacob = drellyan_res_ps(x, drellyan_in, pp1, resuvars);
//
    f = drellyan_res_calc(x, drellyan_in, resuvars, pp1, randsjacob);

    //TCRIDGE CHANGE MOD TO Q2 TO HELP SAMPLING
       double mod = 1.;
    // double mod = resuvars.q2;
    f = mod*f;
    
    if(drellyan_in->AFB_flag == 1){

       double y[ndim];
       for(int ii=0; ii<ndim; ii++) y[ii] = x[ii];
       y[2] = 1. - x[2]; // changes costheta --> -costheta
//
       PSpoint pp2;
       resu_PS resuvars2;
       randsjacob = drellyan_res_ps(y, drellyan_in, pp2, resuvars);
       double f1 = drellyan_res_calc(y, drellyan_in, resuvars, pp2, randsjacob);
//
       f = (f1 - f)/2.;
       mod = resuvars.q2;
       f = mod*f;

       double s = 2.*pp2.ss(1,2); double t = -2.*pp2.ss(1,3);
       double costheta = 1. + 2.*t/s; // Note that theta here is the CM angle between parton 1 and

       if(costheta < 0) f = -f;
       if(resuvars.eta < 0) f = -f;

    }

    else if(drellyan_in->AFB_flag == 2){
       mod = resuvars.q2;
       f = mod*f;
    }

// Print event to file

    if(drellyan_in->event_info.save_events==1 || drellyan_in->event_info.save_events==2){
      events_out(drellyan_in->event_info, iter, f/mod, weight, std::sqrt(resuvars.mur2), resuminfo->alpha_QED, resuvars.alphas, pp1, ndim, x);
    }


    return 0;

}


double drellyan_res_calc(const double x[], drellyan_input* drellyan_in, resu_PS& resuvars, PSpoint& pp1, double randsjacob){

// CUTS:
    bool cutcheck = false; //i.e not cut
    cutcheck = drellyan_cuts(drellyan_in, resuvars.q2, resuvars.qt2, resuvars.eta, resuvars.mures2, &pp1);

    double tempf;
    if(cutcheck==true || randsjacob==0.){
      tempf = 0.;
    }
    else{

      ResummationInfo* resuminfo = &drellyan_in->res_1;

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

// Call main resummation function or the main CT function
      if(resuminfo->CT_flag == 1){
        tempf = CT_NLO(resuvars, resuminfo, x+6);
	//TCRIDGE alter by *x/q2 factor
	//tempf = CT_NLO(resuvars, resuminfo, x+6)*resuvars.x/resuvars.q2;
      }
      else{
        tempf = resummed(resuvars, resuminfo);
	// std::cout << "tempf = " << tempf << std::endl;
      }
      //Set to xsec for now to test PS vol
      // tempf =xsection(resuvars,resuminfo);

    }

    double f = tempf*randsjacob;
    // std::cout << "tempf = " << tempf << " randsjacob = " << randsjacob << " f = " << f << std::endl;
    // std::cout << "f*x/q2 = " << f*resuvars.x/resuvars.q2 << std::endl;
    //Test volume of phase space, this is 4 body PS -> should have volume E^4/(24576*pi^5)
    //double f = randsjacob;
    //double f = 1;

    return f;

}
