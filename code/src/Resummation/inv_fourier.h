//header file for inv_fourier.cc

#ifndef INV_FOURIER_H
#define INV_FOURIER_H

struct resu_PS;
struct ResummationInfo;

class intdeo_data {
public:
  resu_PS* resuvars;
  ResummationInfo* resuminfo;
};


double resummed(resu_PS& resuvars, ResummationInfo* resuminfo);

double invbtoqt(resu_PS& resuvars, ResummationInfo* resuminfo);

// Integrand function for intdeo --> invbtoqt
double invres (double b, void* data);

void fitcalc(resu_PS& resuvars, ResummationInfo* resuminfo);

#endif
