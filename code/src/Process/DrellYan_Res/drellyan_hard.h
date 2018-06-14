//Header file for the process dependent parts of the hard function calculations for the case of diphoton

#ifndef drellyan_hard_H
#define drellyan_hard_H

#include <complex>
#include <vector>

class drellyan_input;
class PSpoint;


void sigmaijdrellyancalc(drellyan_input* drellyan_in, PSpoint& PS, std::vector<std::vector<double> >& sigmaij, double alphas, double& q2);

#endif
