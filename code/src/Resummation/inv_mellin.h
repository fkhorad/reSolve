//Header file for inv_mellin.cc which is the program file which takes the resummed form factor, multiplies by the sudakov factor and does the inverse mellin transform.

#ifndef INV_MELLIN_H
#define INV_MELLIN_H

#include <complex>
class PSdep_variables;
class ResummationInfo;


double inversemellin_resummed(std::complex<double> b,
                              PSdep_variables* resu, ResummationInfo* resuminfo);


#endif
 
