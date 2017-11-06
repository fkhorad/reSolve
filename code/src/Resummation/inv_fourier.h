//header file for inv_fourier.cc

#ifndef INV_FOURIER_H
#define INV_FOURIER_H

class PSdep_variables;
class ResummationInfo;


extern "C" {
  double alphas_(double*);
}

extern "C" {
  void initalphas_(int*, double*, double*, double*, double*, double*, double*);
}


class intdeo_data {
public:
  PSdep_variables* resuvars;
  ResummationInfo* resuminfo;
};


double resummed(PSdep_variables& resuvars, ResummationInfo* resuminfo);

double invbtoqt(PSdep_variables& resuvars, ResummationInfo* resuminfo);

// Integrand function for intdeo --> invbtoqt
double invres (double b, void* data);

void fitcalc(PSdep_variables& resuvars, ResummationInfo* resuminfo);

#endif
