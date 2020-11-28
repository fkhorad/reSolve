//Integrands for drellyan case
#include "DYjet_integrand.h"

#include <fstream>
#include <sstream>
#include <iostream>

#include "phase_space.h"
#include "events_out.h"

// Generic resummation includes
#include "resummation_input.h"
#include "xsection.h"
#include "resu_PS.h"

// Specific includes
#include "DYjet_ps.h"
#include "DYjet_cuts.h"
#include "drellyan_cuts.h"
#include "DYjet_hard.h"
#include "drellyan_input.h"

int DYjet_integrand(int ndim, const double x[], double& f, drellyan_input* drellyan_in, double weight, int iter) {

    ResummationInfo* resuminfo = &drellyan_in->res_1;

// Phase space
    drellyan_amplitude amp1;
    double randsjacob;
    resu_PS resuvars;
    double ss_hat, etaa_hat;

    // std::cout << "ss_hat = " << ss_hat << std::endl;
    // std::cout << "etaa_hat = " << etaa_hat << std::endl;
    randsjacob = DYjet_ps(x, drellyan_in, amp1, resuvars, ss_hat, etaa_hat);
    // std::cout << "ss_hat = " << ss_hat << std::endl;
    // std::cout << "etaa_hat = " << etaa_hat << std::endl;
    //std::cout << "x = " << x << std::endl;
    // std::cout << "randsjacob in integrand= " << randsjacob << std::endl;
//
    f = DYjet_calc(x, drellyan_in, resuvars, amp1, randsjacob, ss_hat, etaa_hat);
    // std::cout << "f = " << f << std::endl;
    //TCRIDGE NEXT TWO LINES TO ALTER SAMPLING AS ALREADY DONE FOR AFB LATER IN THIS FILE, USED TO BE MOD = 1 HERE
    //double mod = resuvars.q2;
    //TCRIDGE need mod=1 for DYjetps validation for some reason
    double mod = 1.0;
    f = mod*f;


    //double mod = 1.; //Commented TCRIDGE FOR ALTERED SAMPLING
    if(drellyan_in->AFB_flag == 1){

       double y[ndim];
       for(int ii=0; ii<ndim; ii++) y[ii] = x[ii];
       y[3] = 1. - x[3]; // changes costheta --> -costheta
//
       drellyan_amplitude amp2;
       resu_PS resuvars2;
       randsjacob = DYjet_ps(y, drellyan_in, amp2, resuvars, ss_hat, etaa_hat);
       double f1 = DYjet_calc(y, drellyan_in, resuvars, amp2, randsjacob, ss_hat, etaa_hat);
//
       f = (f1 - f)/2.;
       mod = resuvars.q2;
       f = mod*f;
       // std::cout << "AFB1" << std::endl;

       double s = 2.*amp2.ss(1,2); double t = -2.*amp2.ss(1,3);
       double costheta = 1. + 2.*t/s; // Note that theta here is the CM angle between parton 1 and

       if(costheta < 0) f = -f;
       if(resuvars.eta < 0) f = -f;

    }

    else if(drellyan_in->AFB_flag == 2){
       mod = resuvars.q2;
       f = mod*f;
       // std::cout << "AFB2" << std::endl;
    }

// Print event to file

    if(drellyan_in->event_info.save_events==1 || drellyan_in->event_info.save_events==2){
      // std::cout << "mod = " << mod << std::endl;
      // std::cout << "f/mod = " << f/mod << std::endl;
      PSpoint pp1 = amp1;
      events_out(drellyan_in->event_info, iter, f/mod, weight, std::sqrt(resuvars.mur2), resuminfo->alpha_QED, resuvars.alphas, pp1, ndim, x);
    }

    return 0;

}


double DYjet_calc(const double x[], drellyan_input* drellyan_in, resu_PS& resuvars, drellyan_amplitude& amp1, double randsjacob, double ss_hat, double etaa_hat){

// CUTS:
    bool cutcheck = false; //i.e not cut
    cutcheck = drellyan_cuts(drellyan_in, resuvars.q2, resuvars.qt2, resuvars.eta, resuvars.mures2, &amp1);
    cutcheck = cutcheck || DYjet_cuts(drellyan_in, &amp1);

    double tempf;
    if(cutcheck==true || randsjacob==0.){
      tempf = 0.;
    }
    else{

      ResummationInfo* resuminfo = &drellyan_in->res_1;

      int Nf = resuminfo->Nf;
      amp1.sigmaij_dyjet_calc(drellyan_in, resuvars.alphas, resuvars.sigmaij);
      // for (int ii=0;ii<=10;ii++) {
      // 	for (int jj=0;jj<=10;jj++) {
      // 	  std::cout << "sigmaij[" << ii << "][" << jj << "]=" << resuvars.sigmaij[ii][jj] << std::endl;
      // 	}
      // }
      tempf = xsection_nlojet(resuvars, resuminfo, ss_hat, etaa_hat);
    }
    // double f = tempf*randsjacob;
    //TCRIDGE TRY x/q2 factor
    //double f = tempf*randsjacob*resuvars.x/resuvars.q2;
    //TCRIDGE tempf=1 for DYjet PS validation
    //tempf = 1;
    double f = tempf*randsjacob;
    // std::cout << "f = " << f << std::endl;
    // std::cout << "tempf = " << tempf << " randsjacob = " << randsjacob << std::endl;
    // std::cout << "f = " << f << std::endl;
    
    return f;

}
