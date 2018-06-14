//Invbtoqt - Does inverse Fourier transform, calls inverse_mellin

#include "inv_fourier.h"

#include <iostream>
#include <ctime>

#include "resummation_input.h"
#include "resu_preproc.h"
#include "intde2.h"
#include "inv_mellin.h"
#include "resu_PS.h"
#include "xsection.h"


#include <fstream>
#include <iomanip>

extern int k_count, k_count0, k_count00;

namespace global_fitpars{
  extern fitparams fitparamsa;
}


int k_count, k_count0 = 0;

double resummed(resu_PS& resuvars, ResummationInfo* resuminfo) {

//   k_count0++;

    double xborn, resummedval, xnorm;

    xborn = xsection(resuvars, resuminfo);

// Extra normalization for dsigma/dqT + flux factor
    xnorm = resuvars.x/resuvars.q2;

    if (resuminfo->resum_flag !=0){
// Get fit sectors information
      fitcalc(resuvars, resuminfo);
//
//Extra normalization which helps compensate for small errors in the fit
      double xborn2, xnormal;
      xborn2 = xsection2(resuvars, resuminfo);
      xnormal = xborn/xborn2;
      if (xnormal> 1.15 || xnormal < 0.85) {
        xnormal = 1.0; //Avoid problems at extreme kinematics (and possibly divisions by 0!)
      }

// Main thing:
      double resuv = invbtoqt(resuvars, resuminfo);
      resummedval = xnorm*xnormal*resuv;
      // std::cout << "resuv = " << resuv << std::endl;
      // std::cout << "xnorm = " << xnorm << std::endl;
      // std::cout << "xnormal = " << xnormal << std::endl;

/*
      double temp = 64000./3.*resuv*3.141592653589793;
      std::ofstream debborah;
      debborah.open("debborah.dat", std::ios_base::app);
      debborah << std::setprecision(17) << temp << "\n";
      debborah << "----------\n\n";
      debborah.close();
*/


    }

    else resummedval = xnorm*xborn;

    return resummedval;
}



double invbtoqt(resu_PS& resuvars, ResummationInfo* resuminfo){

  double qt = std::pow(resuvars.qt2,0.5);
  double errt = 0.0;
  double resum;


  if (qt > 0){
    intdeo_data userdata;
    userdata.resuvars = &resuvars;
    userdata.resuminfo = resuminfo;

    k_count = 0;

    intdeo(invres,0.0,qt,resuminfo->aw,&resum,&errt, &userdata);

// DEBUG
    // std::cout << "# of invres evaluations:  " << k_count << std::endl;
  }
  else {
    std::cout << "qt must be > 0!" << std::endl;
  }

  return resum;

}

double invres (double b, void* userdata){

    std::complex<double> bb;
    double invresval = 0.0;
    intdeo_data* data = (intdeo_data*) userdata;

    bb = std::complex<double>(b,0.0);

    invresval = inversemellin_resummed(bb, data->resuvars, data->resuminfo);

    // std::cout << "b = " << b << " invresval = " << invresval << std::endl;

    return invresval;
}




void fitcalc(resu_PS& resuvars, ResummationInfo* resuminfo){

    int verbosity = resuminfo->verbosity;
    double eta = resuvars.eta;

// Set 1st fit parameter IFIT: rapidity
    int IFIT;
    if (eta < -4.50001) { IFIT = 13;}
    else if (eta < -4.0001) { IFIT = 12;}
    else if (eta < -3.50001) { IFIT = 11;}
    else if (eta < -3.0001) { IFIT = 10;}
    else if (eta < -2.0001) { IFIT = 9;}
    else if (eta < -1.0001) { IFIT = 8;}
    else if (eta < 0.0001) { IFIT = 7;}
    else if (eta < 1.0001) { IFIT = 0;}
    else if (eta < 2.0001) { IFIT = 1;}
    else if (eta < 3.0001) { IFIT = 2;}
    else if (eta < 3.50001) { IFIT = 3;}
    else if (eta < 4.0001) { IFIT = 4;}
    else if (eta < 4.50001) { IFIT = 5;}
    else if (eta < 5.0001) { IFIT = 6;}
    else {
//      std::cout << "WARNING: very large eta: " << eta << std::endl;
      IFIT = 6;
    }

//Set 2nd fit parameter IFIT2: CM total final state energy
    int IFIT2 = 0;
    if (resuminfo->muF_flag == 1) {
      double ensecmultiplier = global_fitpars::fitparamsa.en_sec_multiplier;
      double qtemp = global_fitpars::fitparamsa.minmuf;
      double qq = std::pow(resuvars.q2,0.5);

      if (verbosity >= 12) {
        std::cout << "ensecmultiplier = " << ensecmultiplier << std::endl;
        std::cout << "qtemp = qmin = " << qtemp << std::endl;
      }

      qtemp = ensecmultiplier*qtemp;
      while (qtemp < qq) {
        qtemp = ensecmultiplier*qtemp;
        IFIT2 = IFIT2 + 1;
      }
    }

     if (verbosity >= 12) {
        std::cout << "IFIT = " << IFIT << std::endl;
        std::cout << "IFIT2 = " << IFIT2 << std::endl;
     }

    resuvars.ifit = IFIT;
    resuvars.ifit2 = IFIT2;

}
