//Header file for inv_mellin.cc which is the program file which takes the resummed form factor, multiplies by the sudakov factor and does the inverse mellin transform.

#ifndef INV_MELLIN_H
#define INV_MELLIN_H

#include <complex>

struct resu_PS;
struct ResummationInfo;


double inversemellin_resummed(std::complex<double> b, resu_PS* resu, ResummationInfo* resuminfo);


#endif
