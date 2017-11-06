//Header file for sudakov.cc

#ifndef SUDAKOV_H
#define SUDAKOV_H

#include <complex>
#include <iostream>

#include "resu_procindep.h"


std::complex<double> Sc(std::complex<double> b, double qsquared, double b0p, double alphas, double muresumsquared, double a, resummationIndepfns* resufns, int order, double gqnp, double ggnp, char qorg, int verbosity);
std::complex<double> f0(std::complex<double> y, double beta0, double A1q, double A1g, char qorg, int verbosity);
std::complex<double> f1(std::complex<double> y, double q2, double a, double beta0, double beta1, double A1q, double A1g, double A2q, double A2g, double B1q, double B1g, double mur2, char qorg, int verbosity); 
std::complex<double> f2(std::complex<double> y, double q2, double a, double beta0, double beta1, double beta2, double A1q, double A1g, double A2q, double A2g, double B1q, double B1g, double A3q, double A3g, double B2q, double B2g, double C1qqn, double C1ggn, double mur2, char qorg, int verbosity); 

#endif
