//Header file for the process dependent parts of the hard function calculations for the case of diphoton

#ifndef diphoton_hard_H
#define diphoton_hard_H

#include <complex>

#include "phase_space.h"
#include "lorentz.h"
#include "constants.h"
#include "diphoton_input.h"

void sigmaijdiphotoncalc(diphoton_input* diph_in, PSpoint& PS, double jacob,
                         std::vector<std::vector<double> >& sigmaij, double alphas);

double ggBoxdiphotoncalc(double costheta, int verbosity);

double H1qdiphoton_DY(double costheta, int verbosity);

double H1diphotoncalc(PSpoint* PS, double Cf, int verbosity);

double H2qdiphotoncalc_DY (PSpoint* PS, double costheta, int quark,
                           double Cf, int Nf, int Nc, int verbosity);

double H2stqqdiphotoncalc(PSpoint* PS, int quark,
                          double Cf, int Nf, int Nc, int verbosity);

double Li4(double x);
double Li3(double x);
double Li3fn(double xx);
double Li2(double x);
double Li2fn (double xx);
double Asmitad(double u, double t, double s, int verbosity);
double Bsmitad(double u,double t,double s,double muR2, int verbosity);
double D2smitad(double u,double t,double s,double muR2, int verbosity);
double E3smitad(double u,double t,double s,double muR2, int verbosity);
double G1smitad(double u,double t,double s,double muR2, int verbosity);
double Finite0x2(double t, double u, double s, double muR2, double Qsqd, double Cf, int Nf, int Nc, int verbosity);
double Finite1x1(double t, double u, double s,double muR2, double Cf, int Nc, int verbosity);


#endif
