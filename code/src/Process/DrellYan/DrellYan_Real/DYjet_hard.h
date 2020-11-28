//Header file for the process dependent parts of the hard function calculations for the case of diphoton

#ifndef DYjet_hard_H
#define DYjet_hard_H

#include <complex>
#include <vector>
#include "phase_space.h"

class drellyan_input;

class drellyan_amplitude : public PSpoint{

public:
  void sigmaij_dyjet_calc(drellyan_input* drellyan_in, double as, std::vector<std::vector<double> >& sigmaij);

private:
  void ampkinfacssq_tom(std::vector<double> & ampfacsq_tom_vec,std::vector<std::complex<double> > ampfac_tom_vec);
  void ampkinfacs(int i1, int i2, int i3, int i4, int i5, std::vector<std::complex<double> >& ampfac);
  void ampkinfacs_tom(std::vector<std::complex<double> >& ampfac_tom_vec, int proc);
  std::complex<double> ampfac_tom(int i1, int i2, int i3, int i4, int i5,int bra_ket, int gq);
  
};


#endif
