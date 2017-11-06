//To calculate the process independent functions A1, A2, B1, B2 etc that appear in resummation in sudakov form factor and in the hard factor, cf with init_resu.f of old code. Here only done N-independent parameters so far. By T Cridge

#ifndef RESU_PROCINDEP_H
#define RESU_PROCINDEP_H

#include <iostream>

#include "constants.h"


class resummationIndepfns {
 public:
  double beta0, beta1, beta2, A1g, A2g, A3g, B1g, B2g, C1ggn, A1q, A2q, A3q, B1q, B2q, C1qqn;
};

void resu_init_0(resummationIndepfns& b, int Nf, double Ca, double Cf);


#endif

