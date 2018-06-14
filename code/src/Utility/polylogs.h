#ifndef _polylogsH
#define _polylogsH_

#include <complex>

double Li2fn (double xx);
double Li2(double x);
std::complex<double> Li2(double x, double y);
double Li3fn(double xx);
double Li3(double x);
double Li4(double x);

std::complex<double> Log(double x);
std::complex<double> L0(double x, double y);
std::complex<double> L1(double x, double y);
std::complex<double> L2(double x, double y);
std::complex<double> Ls1(double x, double y, double z, double w);
std::complex<double> Ls0(double x, double y, double z, double w);
std::complex<double> Ls_1(double x, double y, double z, double w);


#endif
