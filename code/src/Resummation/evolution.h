#ifndef _evolutionH_
#define _evolutionH_

#include <complex>
#include <cmath>
#include <iostream>

#include "resu_procindep.h"


void RENO2 ( std::complex<double> XN, int I, int ISIGN, int IBEAM, int flag1, int nf,
             std::complex<double> QQI, std::complex<double> QGF, std::complex<double> GQI,
             std::complex<double> GGI, std::complex<double> GGF, std::complex<double> NS1MI,
             std::complex<double> NS1PI, std::complex<double> NS1F, std::complex<double> QQ1F,
             std::complex<double> QG1F, std::complex<double> GQ1I, std::complex<double> GQ1F,
             std::complex<double> GG1I,std::complex<double>GG1F, std::complex<double> UVI,
             std::complex<double> DVI, std::complex<double> USI, std::complex<double> DSI,
             std::complex<double> SSI, std::complex<double>GLI, std::complex<double>CHI,
             std::complex<double> BOI,
             double ALPS, double q2, double b0p, std::complex<double> ALPQ,
             double a_param, double as, double mur2, int verbosity,
            std::complex<double>* FN);


void ANOM (std::complex<double>& ANS, std::complex<double>& AM, std::complex<double>& AP, std::complex<double>& AL,
           std::complex<double>& BE, std::complex<double>& AB, std::complex<double>& RMIN, std::complex<double>& RPLUS,
           std::complex<double>& RQQ, std::complex<double>& RQG, std::complex<double>& RGQ, std::complex<double>& RGG,
           std::complex<double>& C2Q, std::complex<double>& C2G, std::complex<double>& CDYQ,
           std::complex<double>& CDYG, std::complex<double>& XN, int& FR, std::complex<double>& QQI,
           std::complex<double>& QGF, std::complex<double>& GQI, std::complex<double>& GGI, std::complex<double>& GGF,
           std::complex<double>& NS1MI, std::complex<double>& NS1PI, std::complex<double>& NS1F,
           std::complex<double>& QQ1F, std::complex<double>& QG1F, std::complex<double>& GQ1I,
           std::complex<double>& GQ1F, std::complex<double>& GG1I, std::complex<double>& GG1F,
           std::complex<double>& C2QI, std::complex<double>& C2GF, std::complex<double>& CDYQI,
           std::complex<double>& CDYGI, int verbosity, int nf);


void alphaslcalc(double q2, double b0p, double a_param, double as, double mur2,
                 std::complex<double> nq2, resummationIndepfns* resufns,
                 std::complex<double>& alphasl, std::complex<double>& aexp,
                 std::complex<double>& aexpb, int verbosity, int flag1);


#endif
