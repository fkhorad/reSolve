
#include "resu_procindep.h"


//To calculate the process independent functions A1, A2, B1, B2 etc that appear in resummation in sudakov form factor and in the hard factor, cf with init_resu.f of old code. Here only done N-independent parameters so far. By T Cridge


void resu_init_0(resummationIndepfns& b, int Nf, double Ca, double Cf) {
  
  double beta0temp = 0.0, beta1temp = 0.0, beta2temp = 0.0, A1gtemp = 0.0, A2gtemp = 0.0, A3gtemp = 0.0, B1gtemp = 0.0, B2gtemp = 0.0, C1ggntemp = 0.0, A1qtemp = 0.0, A2qtemp = 0.0, A3qtemp = 0.0, B1qtemp = 0.0, B2qtemp = 0.0, C1qqntemp = 0.0;
  
  double pi = k_constants::pi;
  double zeta3 = k_constants::zeta3;

  beta0temp = (33 - 2*Nf)/12.0;
  beta1temp = (153-19*Nf)/24.0;
  beta2temp = 2857/128.0 - 5033*Nf/1152.0 + 325*Nf*Nf/3456.0;
//Gluon coefficients
  A1gtemp = Ca;
  A2gtemp = Ca/2*((67/18.0 - pi*pi/6.0)*Ca - 5/9.0*Nf);
  A3gtemp = Ca*(13.81-2.15*Nf - Nf*Nf/108.0);
  B1gtemp = -(11*Ca-2*Nf)/6.0;
  B2gtemp = Ca*Ca*(23/24.0 + (11*pi*pi)/18.0 - 3*zeta3/2.0) + Cf*Nf/2.0 - Ca*Nf*(1/12.0 + pi*pi/9) -11/8.0*Cf*Ca;
  C1ggntemp = (pi*pi/2 + 11/2.0 + pi*pi)/2.0; // Only N-independent part
//Quark coefficients
  A1qtemp=Cf;
  A2qtemp=Cf/2.0*(67/6.0-(pi*pi)/2.0-5/9.0*Nf);
  A3qtemp=Cf*(13.81-2.15*Nf-Nf*Nf/108.0) + Cf*(Ca*(29.9259-28.0*zeta3)-8.2963*Nf/2)*2*(beta0temp*4.0)/64.0; //A3 from Becher & Neubert
  B1qtemp=-(3.0*Cf)/2.0;
  B2qtemp=Cf*Cf*(pi*pi/4-3/16.0-3*zeta3) + Ca*Cf*(11*pi*pi/36-193/48.0+3*zeta3/2.0) + Cf*Nf*(17/24.0-pi*pi/18);
  C1qqntemp=Cf/2.0*(pi*pi/2-4.0);// Only N-independent part

  b.beta0 = beta0temp;
  b.beta1 = beta1temp;
  b.beta2 = beta2temp;
  b.A1g = A1gtemp;
  b.A2g = A2gtemp;
  b.A3g = A3gtemp;
  b.B1g = B1gtemp;
  b.B2g = B2gtemp;
  b.C1ggn = C1ggntemp;
  b.A1q = A1qtemp;
  b.A2q = A2qtemp;
  b.A3q = A3qtemp;
  b.B1q = B1qtemp;
  b.B2q = B2qtemp;
  b.C1qqn = C1qqntemp;
}
