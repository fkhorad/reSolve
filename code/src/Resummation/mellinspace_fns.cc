#include "mellinspace_fns.h"

#include <iostream>


// Initialisation routine which evaluates the (process-independent) N-dependent resummation functions at all points which will be needed for the inverse Mellin transform.

// Inverse Mellin transform is done over a contour in the complex plane with a postive and negative branch from (approximate) -infinity to +infinity - contour points and weights are also set here.


void inv_mel_init(contour_info& contours, std::array<resuDep_N, k_constants::mellin_points>& PosBranch, std::array<resuDep_N, k_constants::mellin_points>& NegBranch, int verbosity, int Nf, double Ca, double Cf) {

    std::complex<double> OneComplex;
    OneComplex = std::complex<double>(1.0,0.0);

    double pi = k_constants::pi;
    int m_points = k_constants::mellin_points;

// WZ and ZS are constants used in setting weights
// A Gaussian quadrature of fixed order 8 is used here.
    const double WZ[8] = {0.101228536290376, 0.222381034453374, 0.313706645877887, 0.362683783378362, 0.362683783378362, 0.313706645877887, 0.222381034453374, 0.101228536290376};
    const double ZS[8] = {-0.960289856497536, -0.796666477413627, -0.525532409916329, -0.183434642495650, 0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536};
//Phi is related to angle contour is taken at in the complex plane
    double phi = 0, cosphi = 0, sinphi = 0;
    std::complex<double> Xn;
    std::complex<double> phicomplexp, phicomplexm;
    phicomplexp = std::complex<double>(0.0,0.0);
    phicomplexm = std::complex<double>(0.0,0.0);
    phi = pi*3./4;
    cosphi = std::cos(phi);
    sinphi = std::sin(phi);
    phicomplexp = std::complex<double>(cosphi,sinphi);
    phicomplexm = std::complex<double>(cosphi,-sinphi);
    contours.phicomplexp = phicomplexp;
    contours.phicomplexm = phicomplexm;
    double Z =0;
    std::vector<double> weights(m_points);
    std::vector<std::complex<double> > Np(m_points), Nm(m_points);
    for (int p = 0; p<m_points; p++) {
      Np[p] = std::complex<double>(0,0);
      Nm[p] = std::complex<double>(0,0);
    }
// This part contains the UNWANTED fixed constant "17" -- basically it's a repeated Gaussian
// rule. This is a "bad thing" of different scale of the fixed order of the Gaussian integration.
// To be changed later.

// up and down are integration contour parameters
    const double down[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0,
      21.0, 24.0, 27.0, 30.0, 33.0};
    double up[17];
    double sum = 0, diff = 0;
    for (int i = 0; i<16; i++) {
        up[i] = down[i+1];
    }
    up[16] = 36.0;
    int k = 0;
    for (int i = 0; i<17; i++) {
      sum = up[i]+down[i];
      diff = up[i]-down[i];
      if (verbosity >= 25) {
        std::cout << up[i] << "       " << down[i] << std::endl;
      }
      for (int j = 0; j<8; j++) {
// Initialisation for Inverse Mellin transform integral
        Z = 0.5*(sum + diff*ZS[j]);
        weights[k] = diff/2.0*WZ[j];
        Np[k] = std::complex<double>(2+cosphi*Z,sinphi*Z);
        Nm[k] = std::complex<double>(2+cosphi*Z,-sinphi*Z);
        k = k+1;
      }
    }
    if (verbosity >= 25) {
      for (int k = 0; k< m_points; k++) {
        std::cout << k << " " << Np[k] << " " << Nm[k] << std::endl;
        std::cout << k << " " << weights[k] << std::endl;
      }
    }
    contours.weights = weights;
    contours.Np = Np;
    contours.Nm = Nm;


//Set up positive branch of contour
    std::complex<double> ZeroComplex;
    ZeroComplex = std::complex<double>(0.0,0.0);
    std::vector<std::complex<double> > C1qgpos(m_points), C1gqpos(m_points), C1qqpos(m_points), C1ggpos(m_points);
    std::vector<std::complex<double> > qqip(m_points), qgfp(m_points), gqip(m_points), ggip(m_points), ggfp(m_points), ns1mip(m_points), ns1pip(m_points),
        ns1fp(m_points), qq1fp(m_points), qg1fp(m_points), gq1ip(m_points), gq1fp(m_points), gg1ip(m_points), gg1fp(m_points);

    for (int k=0; k<m_points; k++) {
      Xn = Np[k];
      if (verbosity >= 50) {
        std::cout << "k = " << k << " Xn = " << Xn << std::endl;
      }

//C1 coefficients on positive branch
      C1qgpos[k] = 1.0/((Xn+OneComplex)*(Xn+2.0*OneComplex));
      C1gqpos[k] = 4.0/(3.0*(Xn+OneComplex));
      C1qqpos[k] = 2*pi*pi/3.0-16/3.0+4.0/(3.0*Xn*(Xn+OneComplex));
      C1ggpos[k] = pi*pi/2+11/2.0+pi*pi;

      if (verbosity >= 50) {
        std::cout << "C1qgpos[" << k << "] = " << C1qgpos[k] << std::endl;
        std::cout << "C1gqpos[" << k << "] = " << C1gqpos[k] << std::endl;
        std::cout << "C1qqpos[" << k << "] = " << C1qqpos[k] << std::endl;
        std::cout << "C1ggpos[" << k << "] = " << C1ggpos[k] << std::endl;
        std::cout << "C1 coefficients along positive branch of contour: " << std::endl;
        std::cout << C1qgpos[k].real() << "   " << C1qgpos[k].imag() << "   "
          << C1gqpos[k].real() << "   " << C1gqpos[k].imag() << "   " << C1qqpos[k].real()
          << "   " << C1qqpos[k].imag() << "   " << C1ggpos[k].real() << "   "
          << C1ggpos[k].imag() << std::endl;
      }


//Anomalous dimensions on positive branch
      std::complex<double> tempp[15];
      for (int i = 0; i<15; i++) { tempp[i] = ZeroComplex;}
      *anomcalc(tempp,Xn, verbosity);
      qqip[k] = tempp[1];
      qgfp[k] = tempp[2];
      gqip[k] = tempp[3];
      ggip[k] = tempp[4];
      ggfp[k] = tempp[5];
      ns1mip[k] = tempp[6];
      ns1pip[k] = tempp[7];
      ns1fp[k] = tempp[8];
      qq1fp[k] = tempp[9];
      qg1fp[k] = tempp[10];
      gq1ip[k] = tempp[11];
      gq1fp[k] = tempp[12];
      gg1ip[k] = tempp[13];
      gg1fp[k] = tempp[14];
      if (verbosity >= 50) {
        std::cout << "qqip[" << k << "] = " << qqip[k] << std::endl;
        std::cout << "qgfp[" << k << "] = " << qgfp[k] << std::endl;
        std::cout << "gqip[" << k << "] = " << gqip[k] << std::endl;
        std::cout << "ggip[" << k << "] = " << ggip[k] << std::endl;
        std::cout << "ggfp[" << k << "] = " << ggfp[k] << std::endl;
        std::cout << "ns1mip[" << k << "] = " << ns1mip[k] << std::endl;
        std::cout << "ns1pip[" << k << "] = " << ns1pip[k] << std::endl;
        std::cout << "ns1fp[" << k << "] = " << ns1fp[k] << std::endl;
        std::cout << "qq1fp[" << k << "] = " << qq1fp[k] << std::endl;
        std::cout << "qg1fp[" << k << "] = " << qg1fp[k] << std::endl;
        std::cout << "gq1ip[" << k << "] = " << gq1ip[k] << std::endl;
        std::cout << "gq1fp[" << k << "] = " << gq1fp[k] << std::endl;
        std::cout << "gg1ip[" << k << "] = " << gg1ip[k] << std::endl;
        std::cout << "gg1fp[" << k << "] = " << gg1fp[k] << std::endl;
      }
    }


//Set up negative branch of contour
    std::vector<std::complex<double> > C1qgneg(m_points), C1gqneg(m_points), C1qqneg(m_points), C1ggneg(m_points);
    std::vector<std::complex<double> > qqim(m_points), qgfm(m_points), gqim(m_points), ggim(m_points), ggfm(m_points), ns1mim(m_points), ns1pim(m_points),
        ns1fm(m_points), qq1fm(m_points), qg1fm(m_points), gq1im(m_points), gq1fm(m_points), gg1im(m_points), gg1fm(m_points);

    for (int k=0; k<m_points; k++) {
      Xn = Nm[k];

//C1 coefficients on negative branch
      C1qgneg[k] = 1.0/((Xn+OneComplex)*(Xn+2.0*OneComplex));
      C1gqneg[k] = 4.0/(3.0*(Xn+OneComplex));
      C1qqneg[k] = 2*pi*pi/3.0-16/3.0+4.0/(3.0*Xn*(Xn+OneComplex));
      C1ggneg[k] = pi*pi/2+11/2.0+pi*pi;
      if (verbosity >= 50) {
        std::cout << k << "   " << C1qgneg[k] << "   " << C1gqneg[k]
          << "   " << C1qqneg[k] << "   " << C1ggneg[k] << std::endl;
        std::cout << C1qgneg[k].real() << "   " << C1qgneg[k].imag() << "   "
          << C1gqneg[k].real() << "   " << C1gqneg[k].imag() << "   " << C1qqneg[k].real()
          << "   " << C1qqneg[k].imag() << "   " << C1ggneg[k].real() << "   "
          << C1ggneg[k].imag() << std::endl;
      }


//Anomalous dimensions on negative branch
      std::complex<double> tempm[15];
      for (int i = 0; i<15; i++) { tempm[i] = std::complex<double>(0.0,0.0);}
      *anomcalc(tempm,Xn, verbosity);
      qqim[k] = tempm[1];
      qgfm[k] = tempm[2];
      gqim[k] = tempm[3];
      ggim[k] = tempm[4];
      ggfm[k] = tempm[5];
      ns1mim[k] = tempm[6];
      ns1pim[k] = tempm[7];
      ns1fm[k] = tempm[8];
      qq1fm[k] = tempm[9];
      qg1fm[k] = tempm[10];
      gq1im[k] = tempm[11];
      gq1fm[k] = tempm[12];
      gg1im[k] = tempm[13];
      gg1fm[k] = tempm[14];
      if (verbosity >= 50) {
        std::cout << "Xn = " << Xn << " qqim[" << k << "] = " << qqim[k] << std::endl;
        std::cout << "qqim[" << k << "] = " << qqim[k] << std::endl;
        std::cout << "qgfm[" << k << "] = " << qgfm[k] << std::endl;
        std::cout << "gqim[" << k << "] = " << gqim[k] << std::endl;
        std::cout << "ggim[" << k << "] = " << ggim[k] << std::endl;
        std::cout << "ggfm[" << k << "] = " << ggfm[k] << std::endl;
        std::cout << "ns1mim[" << k << "] = " << ns1mim[k] << std::endl;
        std::cout << "ns1pim[" << k << "] = " << ns1pim[k] << std::endl;
        std::cout << "ns1fm[" << k << "] = " << ns1fm[k] << std::endl;
        std::cout << "qq1fm[" << k << "] = " << qq1fm[k] << std::endl;
        std::cout << "qg1fm[" << k << "] = " << qg1fm[k] << std::endl;
        std::cout << "gq1im[" << k << "] = " << gq1im[k] << std::endl;
        std::cout << "gq1fm[" << k << "] = " << gq1fm[k] << std::endl;
        std::cout << "gg1im[" << k << "] = " << gg1im[k] << std::endl;
        std::cout << "gg1fm[" << k << "] = " << gg1fm[k] << std::endl;
      }
    }


//C2 coefficients
    std::vector<std::complex<double> > C2qgpos(m_points), C2NSqqbpos(m_points), C2NSqqpos(m_points), C2Sqqbpos(m_points);
    std::vector<std::complex<double> > C2qgneg(m_points), C2NSqqbneg(m_points), C2NSqqneg(m_points), C2Sqqbneg(m_points);
    for (int k = 0; k<m_points; k++) {
      C2values C2s;
      C2s.C2qgpos = std::complex<double>(0.0,0.0);
      C2s.C2NSqqbpos = std::complex<double>(0.0,0.0);
      C2s.C2NSqqpos = std::complex<double>(0.0,0.0);
      C2s.C2Sqqbpos = std::complex<double>(0.0,0.0);
      C2s.C2qgneg = std::complex<double>(0.0,0.0);
      C2s.C2NSqqbneg = std::complex<double>(0.0,0.0);
      C2s.C2NSqqneg = std::complex<double>(0.0,0.0);
      C2s.C2Sqqbneg = std::complex<double>(0.0,0.0);

// Main calc
      C2valuescalc(C2s, Np[k], Nm[k], verbosity, Nf, Ca, Cf);

// Assign positive branch
      C2qgpos[k] = C2s.C2qgpos;
      C2NSqqbpos[k] = C2s.C2NSqqbpos;
      C2NSqqpos[k] = C2s.C2NSqqpos;
      C2Sqqbpos[k] = C2s.C2Sqqbpos;

// Assign negative branch
      C2qgneg[k] = C2s.C2qgneg;
      C2NSqqbneg[k] = C2s.C2NSqqbneg;
      C2NSqqneg[k] = C2s.C2NSqqneg;
      C2Sqqbneg[k] = C2s.C2Sqqbneg;

      if (verbosity >= 50) {
        std::cout << "Np[" << k << "] = " << Np[k] << std::endl;
        std::cout << "Nm[" << k << "] = " << Nm[k] << std::endl;
        std::cout << "k = " << k << std::endl;
        std::cout << "C2qgpos[" << k << "] = " << C2qgpos[k] << std::endl;
        std::cout << "C2NSqqbpos[" << k << "] = " << C2NSqqbpos[k] << std::endl;
        std::cout << "C2NSqqpos[" << k << "] = " << C2NSqqpos[k] << std::endl;
        std::cout << "C2Sqqbpos[" << k << "] = " << C2Sqqbpos[k] << std::endl;
        std::cout << "C2qgneg[" << k << "] = " << C2qgneg[k] << std::endl;
        std::cout << "C2NSqqbneg[" << k << "] = " << C2NSqqbneg[k] << std::endl;
        std::cout << "C2NSqqneg[" << k << "] = " << C2NSqqneg[k] << std::endl;
        std::cout << "C2Sqqbneg[" << k << "] = " << C2Sqqbneg[k] << std::endl;
        std::cout << "k = " << k << std::endl;
        std::cout << "qqip[k] = " << qqip[k] << std::endl;
        std::cout << "qgfp[k] = " << qgfp[k] << std::endl;
        std::cout << "gqip[k] = " << gqip[k] << std::endl;
        std::cout << "ggip[k] = " << ggip[k] << std::endl;
        std::cout << "ggfp[k] = " << ggfp[k] << std::endl;
        std::cout << "ns1mip[k] = " << ns1mip[k] << std::endl;
        std::cout << "ns1pip[k] = " << ns1pip[k] << std::endl;
        std::cout << "ns1fp[k] = " << ns1fp[k] << std::endl;
        std::cout << "qq1fp[k] = " << qq1fp[k] << std::endl;
        std::cout << "qg1fp[k] = " << qg1fp[k] << std::endl;
        std::cout << "gq1ip[k] = " << gq1ip[k] << std::endl;
        std::cout << "gg1ip[k] = " << gg1ip[k] << std::endl;
        std::cout << "gq1fp[k] = " << gq1fp[k] << std::endl;
        std::cout << "gg1fp[k] = " << gg1fp[k] << std::endl;
        std::cout << "qqim[k] = " << qqim[k] << std::endl;
        std::cout << "qgfm[k] = " << qgfm[k] << std::endl;
        std::cout << "gqim[k] = " << gqim[k] << std::endl;
        std::cout << "ggfm[k] = " << ggfm[k] << std::endl;
        std::cout << "ns1mim[k] = " << ns1mim[k] << std::endl;
        std::cout << "ns1pim[k] = " << ns1pim[k] << std::endl;
        std::cout << "ns1fm[k] = " << ns1fm[k] << std::endl;
        std::cout << "qq1fm[k] = " << qq1fm[k] << std::endl;
        std::cout << "qg1fm[k] = " << qg1fm[k] << std::endl;
        std::cout << "gq1im[k] = " << gq1im[k] << std::endl;
        std::cout << "gq1fm[k] = " << gq1fm[k] << std::endl;
        std::cout << "gg1im[k] = " << gg1im[k] << std::endl;
        std::cout << "gg1fm[k] = " << gg1fm[k] << std::endl;
        std::cout << "C1qgpos[k] = " << C1qgpos[k] << std::endl;
        std::cout << "C1gqpos[k] = " << C1gqpos[k] << std::endl;
        std::cout << "C1qqpos[k] = " << C1qqpos[k] << std::endl;
        std::cout << "C1ggpos[k] = " << C1ggpos[k] << std::endl;
        std::cout << "C2qgpos[k] = " << C2qgpos[k] << std::endl;
        std::cout << "C2NSqqbpos[k] = " << C2NSqqbpos[k] << std::endl;
        std::cout << "C2Sqqbpos[k] = " << C2Sqqbpos[k] << std::endl;
        std::cout << "C2NSqqpos[k] = " << C2NSqqpos[k] << std::endl;
        std::cout << "C2qgneg[k] = " << C2qgneg[k] << std::endl;
        std::cout << "C2NSqqbneg[k] = " << C2NSqqbneg[k] << std::endl;
        std::cout << "C2Sqqbneg[k] = " << C2Sqqbneg[k] << std::endl;
        std::cout << "C2NSqqneg[k] = " << C2NSqqneg[k] << std::endl;
      }

    }


// Set values for output
    std::complex<double> Nfcomplex = std::complex<double>(Nf,0.0);
    for(int kk=0; kk<m_points; kk++){
      PosBranch[kk].qqi = qqip[kk];
      PosBranch[kk].qgf = qgfp[kk];
      PosBranch[kk].gqi = gqip[kk];
      PosBranch[kk].ggi = ggip[kk];
      PosBranch[kk].ggf = ggfp[kk];
      PosBranch[kk].ns1mi = ns1mip[kk];
      PosBranch[kk].ns1pi = ns1pip[kk];
      PosBranch[kk].ns1f = ns1fp[kk];
      PosBranch[kk].qq1f = qq1fp[kk];
      PosBranch[kk].qg1f = qg1fp[kk];
      PosBranch[kk].gq1i = gq1ip[kk];
      PosBranch[kk].gg1i = gg1ip[kk];
      PosBranch[kk].gq1f = gq1fp[kk];
      PosBranch[kk].gg1f = gg1fp[kk];
      PosBranch[kk].C1qg = C1qgpos[kk];
      PosBranch[kk].C1gq = C1gqpos[kk];
      PosBranch[kk].C1qq = C1qqpos[kk];
      PosBranch[kk].C1gg = C1ggpos[kk];
      PosBranch[kk].C2qg = C2qgpos[kk];
      PosBranch[kk].C2NSqqb = C2NSqqbpos[kk];
      PosBranch[kk].C2Sqqb = C2Sqqbpos[kk];
      PosBranch[kk].C2NSqq = C2NSqqpos[kk];
// Define renormalised anomalous dimensions from basic coefficients
      PosBranch[kk].gamma1qq = -(qqip[kk]/4.0);
      PosBranch[kk].gamma1qg = -(qgfp[kk]/8.0);
      PosBranch[kk].gamma1gq = -(gqip[kk]/4.0);
      PosBranch[kk].gamma1gg = -((ggip[kk] + Nfcomplex*ggfp[kk])/4.0);
      PosBranch[kk].gamma2qq = -(((ns1pip[kk] + Nfcomplex*ns1fp[kk]) + Nfcomplex*qq1fp[kk])/8.0);
      PosBranch[kk].gamma2qqV = -(ns1pip[kk] + 2.0*Nfcomplex*ns1fp[kk] + ns1mip[kk])/16.0;
      PosBranch[kk].gamma2qqbV = -(ns1pip[kk]-ns1mip[kk])/16.0;
      PosBranch[kk].gamma2qqS = -(qq1fp[kk]/16.0);
      PosBranch[kk].gamma2qqbS = -(qq1fp[kk]/16.0);
      PosBranch[kk].gamma2qg = -(qg1fp[kk]/16.0);
      PosBranch[kk].gamma2gq = -((gq1ip[kk] + Nfcomplex*gq1fp[kk])/8.0);
      PosBranch[kk].gamma2gg = -((gg1ip[kk] + Nfcomplex*gg1fp[kk])/8.0);

      NegBranch[kk].qqi = qqim[kk];
      NegBranch[kk].qgf = qgfm[kk];
      NegBranch[kk].gqi = gqim[kk];
      NegBranch[kk].ggi = ggim[kk];
      NegBranch[kk].ggf = ggfm[kk];
      NegBranch[kk].ns1mi = ns1mim[kk];
      NegBranch[kk].ns1pi = ns1pim[kk];
      NegBranch[kk].ns1f = ns1fm[kk];
      NegBranch[kk].qq1f = qq1fm[kk];
      NegBranch[kk].qg1f = qg1fm[kk];
      NegBranch[kk].gq1i = gq1im[kk];
      NegBranch[kk].gq1f = gq1fm[kk];
      NegBranch[kk].gg1i = gg1im[kk];
      NegBranch[kk].gg1f = gg1fm[kk];
      NegBranch[kk].C1qg = C1qgneg[kk];
      NegBranch[kk].C1gq = C1gqneg[kk];
      NegBranch[kk].C1qq = C1qqneg[kk];
      NegBranch[kk].C1gg = C1ggneg[kk];
      NegBranch[kk].C2qg = C2qgneg[kk];
      NegBranch[kk].C2NSqqb = C2NSqqbneg[kk];
      NegBranch[kk].C2Sqqb = C2Sqqbneg[kk];
      NegBranch[kk].C2NSqq = C2NSqqneg[kk];
// Define renormalised anomalous dimensions from basic coefficients
      NegBranch[kk].gamma1qq = -(qqim[kk]/4.0);
      NegBranch[kk].gamma1qg = -(qgfm[kk]/8.0);
      NegBranch[kk].gamma1gq = -(gqim[kk]/4.0);
      NegBranch[kk].gamma1gg = -((ggim[kk] + Nfcomplex*ggfm[kk])/4.0);
      NegBranch[kk].gamma2qq = -(((ns1pim[kk] + Nfcomplex*ns1fm[kk]) + Nfcomplex*qq1fm[kk])/8.0);
      NegBranch[kk].gamma2qqV = -(ns1pim[kk] + 2.0*Nfcomplex*ns1fm[kk] + ns1mim[kk])/16.0;
      NegBranch[kk].gamma2qqbV = -(ns1pim[kk]-ns1mim[kk])/16.0;
      NegBranch[kk].gamma2qqS = -(qq1fm[kk]/16.0);
      NegBranch[kk].gamma2qqbS = -(qq1fm[kk]/16.0);
      NegBranch[kk].gamma2qg = -(qg1fm[kk]/16.0);
      NegBranch[kk].gamma2gq = -((gq1im[kk] + Nfcomplex*gq1fm[kk])/8.0);
      NegBranch[kk].gamma2gg = -((gg1im[kk] + Nfcomplex*gg1fm[kk])/8.0);
    }

}



//C2 coefficients at each value of x along the contour as N dependent
void C2valuescalc(C2values& c, std::complex<double> xp, std::complex<double> xm, int verbosity, int Nf, double Ca, double Cf) {
    std::complex<double> x = std::complex<double>(0.0,0.0);
    x = xp;
    std::complex<double> OneComplex;
    OneComplex = std::complex<double>(1.0,0.0);

// This constant is only used here, I believe
    const double plog41f2=0.5174790616738993;

    double pi = k_constants::pi;
    double zeta2 = k_constants::zeta2;
    double zeta3 = k_constants::zeta3;
    double zeta4 = k_constants::zeta4;
    double Eulerconst = k_constants::Eulerconst;

    double ln2 = std::log(2);

    std::complex<double> XnP1, XnP2, XnP3, XnM1, XnO2P1, XnO2P2, XnO2P1O2, XnO2P3O2, XnO2;
    XnP1 = x + OneComplex;
    XnP2 = x + 2.0*OneComplex;
    XnP3 = x + 3.0*OneComplex;
    XnM1 = x - OneComplex;
    XnO2P1 = x/2.0+OneComplex;
    XnO2P2 = x/2.0+2.0*OneComplex;
    XnO2P3O2 = x/2.0+3.0*OneComplex/2.0;
    XnO2P1O2 = x/2.0+OneComplex/2.0;
    XnO2 = x/2.0;
    double log2 = std::log(2);
    std::complex<double> log2comp;
    log2comp = std::complex<double>(log2,0.0);
    std::complex<double> psi0 (std::complex<double> x, int verbosity);
    std::complex<double> psi1 (std::complex<double> x, int verbosity);
    std::complex<double> psi2 (std::complex<double> x, int verbosity);
    std::complex<double> psi3 (std::complex<double> x, int verbosity);
    std::complex<double> bet (std::complex<double> x, int verbosity);
    std::complex<double> bet1 (std::complex<double> x, int verbosity);
    std::complex<double> bet2 (std::complex<double> x, int verbosity);
    std::complex<double> bet3 (std::complex<double> x, int verbosity);
    std::complex<double> ps0np1, ps0np2, ps0np3, ps0n, ps0nf2p1, ps0nf2p3f2, ps0nf2p1f2, ps0nf2, ps0nf2p2;
    std::complex<double> ps1np1, ps1np2, ps1np3, ps1n, ps1nf2p1, ps1nf2p3f2, ps1nf2p1f2, ps1nf2, ps1nf2p2;
    std::complex<double> ps2np1, ps2np2, ps2np3, ps2n, ps2nf2p1, ps2nf2p3f2, ps2nf2p1f2, ps2nf2, ps2nf2p2;
    std::complex<double> ps3n, ps3nf2p1f2, ps3nf2;
    std::complex<double> beta0n, beta0np1, beta0np2, beta0np3, beta1n, beta2n, beta3n;
    ps0np1 = psi0(XnP1, verbosity);
    ps0np2 = psi0(XnP2, verbosity);
    ps0np3 = psi0(XnP3, verbosity);
    ps0n = psi0(x, verbosity);
    ps0nf2p1 = psi0(XnO2P1, verbosity);
    ps0nf2p3f2 = psi0(XnO2P3O2, verbosity);
    ps0nf2p1f2 = psi0(XnO2P1O2, verbosity);
    ps0nf2 = psi0(XnO2, verbosity);
    ps0nf2p2 = psi0(XnO2P2, verbosity);
    ps1np1 = psi1(XnP1, verbosity);
    ps1np2 = psi1(XnP2, verbosity);
    ps1np3 = psi1(XnP3, verbosity);
    ps1n = psi1(x, verbosity);
    ps1nf2p1 = psi1(XnO2P1, verbosity);
    ps1nf2p3f2 = psi1(XnO2P3O2, verbosity);
    ps1nf2p1f2 = psi1(XnO2P1O2, verbosity);
    ps1nf2 = psi1(XnO2, verbosity);
    ps1nf2p2 = psi1(XnO2P2, verbosity);
    ps2np1 = psi2(XnP1, verbosity);
    ps2np2 = psi2(XnP2, verbosity);
    ps2np3 = psi2(XnP3, verbosity);
    ps2n = psi2(x, verbosity);
    ps2nf2p1 = psi2(XnO2P1, verbosity);
    ps2nf2p3f2 = psi2(XnO2P3O2, verbosity);
    ps2nf2p1f2 = psi2(XnO2P1O2, verbosity);
    ps2nf2 = psi2(XnO2, verbosity);
    ps2nf2p2 = psi2(XnO2P2, verbosity);
    ps3n = psi3(x, verbosity);
    ps3nf2p1f2 = psi3(XnO2P1O2, verbosity);
    ps3nf2 = psi3(XnO2, verbosity);
    beta0n = bet(x, verbosity);
    beta0np1 = bet(XnP1, verbosity);
    beta0np2 = bet(XnP2, verbosity);
    beta0np3 = bet(XnP3, verbosity);
    beta1n = bet1(x, verbosity);
    beta2n = bet2(x, verbosity);
    beta3n = bet3(x, verbosity);
    std::complex<double> s1n, s1np1, s1np2, s1nm1, s1nf2, s1nm1f2, s1np1f2, s1nm2f2, s1np2f2;
    std::complex<double> s2n, s2np1, s2np2, s2nm1, s2nf2, s2nm1f2, s2np1f2, s2nm2f2, s2np2f2;
    std::complex<double> s3n, s3np1, s3np2, s3nm1, s3nf2, s3nm1f2, s3np1f2, s3nm2f2, s3np2f2;
    std::complex<double> s4nm1, s4nm1f2, s4nm2f2 ;
    s1n = ps0np1 + Eulerconst;
    s1np1 = ps0np2 + Eulerconst;
    s1np2 = ps0np3 + Eulerconst;
    s1nm1 = ps0n + Eulerconst;
    s1nf2 = ps0nf2p1 + Eulerconst;
    s1nm1f2 = ps0nf2p1f2 + Eulerconst;
    s1np1f2 = ps0nf2p3f2 + Eulerconst;
    s1nm2f2 = ps0nf2 + Eulerconst;
    s1np2f2 = ps0nf2p2 + Eulerconst;
    s2n = -ps1np1 + zeta2;
    s2np1 = -ps1np2 + zeta2;
    s2np2 = -ps1np3 + zeta2;
    s2nm1 = -ps1n + zeta2;
    s2nf2 = -ps1nf2p1 + zeta2;
    s2nm1f2 = -ps1nf2p1f2 + zeta2;
    s2np1f2 = -ps1nf2p3f2 + zeta2;
    s2nm2f2 = -ps1nf2 + zeta2;
    s2np2f2 = -ps1nf2p2 + zeta2;
    s3n = ps2np1/2.0 + zeta3;
    s3np1 = ps2np2/2.0 + zeta3;
    s3np2 = ps2np3/2.0 + zeta3;
    s3nm1 = ps2n/2.0 + zeta3;
    s3nf2 = ps2nf2p1/2.0 + zeta3;
    s3nm1f2 = ps2nf2p1f2/2.0 + zeta3;
    s3np1f2 = ps2nf2p3f2/2.0 + zeta3;
    s3nm2f2 = ps2nf2/2.0 + zeta3;
    s3np2f2 = ps2nf2p2/2.0 + zeta3;
    s4nm1 = -ps3n/6.0 + zeta4;
    s4nm1f2 = -ps3nf2p1f2/6.0 + zeta4;
    s4nm2f2 = -ps3nf2/6.0 + zeta4;
    std::complex<double> acg21 (std::complex<double> x, int verbosity);
    std::complex<double> acg20 (std::complex<double> x, int verbosity);
    std::complex<double> acg3 (std::complex<double> x, int verbosity);
    std::complex<double> acg6 (std::complex<double> x, int verbosity);
    std::complex<double> acg5 (std::complex<double> x, int verbosity);
    std::complex<double> acg8 (std::complex<double> x, int verbosity);
    std::complex<double> s211nm1, s31nm1, sm1nm1, sm1n, sm1np1, sm1np2, sm2nm1, sm3nm1, sm4nm1;
    std::complex<double> sm21n, sm21np1, sm21np2, sm31nm1, sm22nm1, sm211nm1;
    s211nm1 = -acg21(XnM1, verbosity)+6.0/5.0*zeta2*zeta2;
    s31nm1 = acg20(XnM1, verbosity)+zeta2*s2nm1-0.5*zeta2*zeta2;
    sm1nm1 = beta0n-ln2;
    sm1n = -beta0np1-ln2;
    sm1np1 = beta0np2-ln2;
    sm1np2 = -beta0np3 - ln2;
    sm2nm1 = -beta1n - zeta2/2.0;
    sm3nm1 = beta2n/2.0 - zeta3*0.75;
    sm4nm1 = -beta3n/6.0-zeta4*7.0/8.0;
    sm21n = acg3(x, verbosity) + zeta2*sm1n - 5.0/8.0*zeta3+zeta2*ln2;
    sm21np1 = -acg3(XnP1, verbosity) + zeta2*sm1np1 - 5.0/8.0*zeta3+zeta2*ln2;
    sm21np2 = acg3(XnP2, verbosity) + zeta2*sm1np2 - 5.0/8.0*zeta3+zeta2*ln2;
    sm31nm1 = acg6(XnM1, verbosity) + zeta2*sm2nm1 - zeta3*sm1nm1 - 3.0/5.0*zeta2*zeta2 + 2*plog41f2 + 0.75*zeta3*ln2 - zeta2/2.0*ln2*ln2 + ln2*ln2*ln2*ln2/12.0;
    sm22nm1 = acg5(XnM1, verbosity) -2.0*sm31nm1 + 2.0*zeta2*sm2nm1 + 0.075*zeta2*zeta2;
    sm211nm1 = -acg8(XnM1, verbosity) + zeta3*sm1nm1 - plog41f2 + zeta2*zeta2/8.0 + zeta3*ln2/8.0 + zeta2*ln2*ln2/4.0 - ln2*ln2*ln2*ln2/24.0;
    std::complex<double> acg13 (std::complex<double> x, int verbosity);
    std::complex<double> acg1p (std::complex<double> x, int verbosity);
    std::complex<double> acg1 (std::complex<double> x, int verbosity);
    std::complex<double> acg1pp (std::complex<double> x, int verbosity);
    std::complex<double> acg2p (std::complex<double> x, int verbosity);
    std::complex<double> acg4 (std::complex<double> x, int verbosity);
    std::complex<double> acg4p (std::complex<double> x, int verbosity);
    std::complex<double> acg7 (std::complex<double> x, int verbosity);
    std::complex<double> acg9 (std::complex<double> x, int verbosity);
    std::complex<double> acg13xnm1, acg1xnm1, acg1pxn, acg1pxnp1, acg1ppxnm1, acg2pxnm1, acg4xn, acg4xnp1, acg4pxnm1, acg5xnm1, acg6xnm1, acg7xnm1, acg9xnm1;
    acg13xnm1 = acg13(XnM1, verbosity);
    acg1xnm1 = acg1(XnM1, verbosity);
    acg1pxn = acg1p(x, verbosity);
    acg1pxnp1 = acg1p(XnP1, verbosity);
    acg1ppxnm1 = acg1pp(XnM1, verbosity);
    acg2pxnm1 = acg2p(XnM1, verbosity);
    acg4xn = acg4(x, verbosity);
    acg4xnp1 = acg4(XnP1, verbosity);
    acg4pxnm1 = acg4p(XnM1, verbosity);
    acg5xnm1 = acg5(XnM1, verbosity);
    acg6xnm1 = acg6(XnM1, verbosity);
    acg7xnm1 = acg7(XnM1, verbosity);
    acg9xnm1 = acg9(XnM1, verbosity);
    std::complex<double> C2Sqqb, C2NSqq, C2NSqqb, C2qg;
    C2Sqqb = std::complex<double>(0.0,0.0);
    C2NSqq = std::complex<double>(0.0,0.0);
    C2NSqqb = std::complex<double>(0.0,0.0);
    C2qg = std::complex<double>(0.0,0.0);
    C2Sqqb = Cf*(172.0/XnM1-54.0/(x*x*x*x) - 27.0/(x*x*x) - 126.0/(x*x) - 315.0/x - 54.0/(XnP1*XnP1*XnP1*XnP1) - 27.0/(XnP1*XnP1*XnP1) + 180.0/(XnP1*XnP1) + 279.0/XnP1 - 72.0/(XnP2*XnP2*XnP2) - 192.0/(XnP2*XnP2) - 136.0/XnP2 - (72.0*zeta2)/XnM1 + (108.0*zeta2)/x - (108.0*zeta2)/XnP1 + (72.0*zeta2)/XnP2 + (72.0*s2nm1)/XnM1 - (108.0*s2n)/x + (108.0*s2np1)/XnP1 - (72.0*s2np2)/XnP2)/432.0;
    C2NSqq = Cf*Nf*(127.0/192.0 - 1.0/(24.0*x*x*x) + 5.0/(72.0*x*x) - 37.0/(216.0*x) - 1.0/(24.0*XnP1*XnP1*XnP1) + 5.0/(72.0*XnP1*XnP1) - 19.0/(216.0*XnP1) - (2.0*zeta2)/3.0 + (7.0*zeta3)/36.0 - (7.0*s1nm1)/27.0 + (5.0*s2nm1)/36.0 -s3nm1/12.0) + Ca*Cf*(-1535.0/384.0 - 1.0/(8.0*x*x*x*x) + 11.0/(48.0*x*x*x) - 47.0/(72.0*x*x) + 50.0/(27.0*x) - 1.0/(8.0*XnP1*XnP1*XnP1*XnP1) - 1.0/(48.0*XnP1*XnP1*XnP1)-83.0/(72.0*XnP1*XnP1) + 1.0/(54.0*XnP1) + (97.0*zeta2)/24.0 - (3.0*zeta2)/(8.0*x) + (3.0*zeta2)/(8.0*XnP1) + (85.0*zeta3)/72.0 - (5.0*zeta3)/(8.0*x) - (5.0*zeta3)/(8.0*XnP1) - (29.0*zeta4)/16.0 + (101.0*s1nm1)/54.0 - (5.0*zeta3*s1nm1)/4.0 + (zeta2*s1nm1*s1nm1)/4.0 + (zeta2*s1n)/(4.0*x) - s1np1/(8.0*XnP1) + (zeta2*s1np1)/(4.0*XnP1) - (19.0*s2nm1)/18.0 + (zeta2*s2nm1)/4.0 - (s1nm1*s1nm1*s2nm1)/4.0 + s2n/(4.0*x*x) + s2n/(4.0*x) - (s1n*s2n)/(4.0*x) + s2np1/(4.0*XnP1*XnP1) - s2np1/(4.0*XnP1) - (s1np1*s2np1)/(4.0*XnP1) + s211nm1/2.0 + (11.0*s3nm1)/24.0 - (s1nm1*s3nm1)/2.0 - s3n/(4.0*x) - s3np1/(4.0*XnP1) - s4nm1/2.0) + Cf*Cf*(255.0/128.0 - 1.0/(8.0*x*x*x*x) + 3.0/(8.0*x*x) -19.0/(8.0*x) - 1.0/(8.0*XnP1*XnP1*XnP1*XnP1) + 1.0/(4.0*XnP1*XnP1*XnP1) + 2.0/(XnP1*XnP1) + 19.0/(8.0*XnP1) - (35.0*zeta2)/16.0 + (3.0*zeta2)/(2.0*x) - (3.0*zeta2)/(2.0*XnP1) - (3.0*zeta3)/2.0 - (3.0*zeta3)/(4.0*x) - (3.0*zeta3)/(4.0*XnP1) + (61.0*zeta4)/16.0 - (3.0*zeta3*s1nm1)/2.0 - (zeta2*s1nm1*s1nm1)/2.0 - s1n/(4.0*x*x) - (zeta2*s1n)/(2.0*x) + s1n*s1n/(8.0*x*x) + s1np1/(4.0*XnP1*XnP1) + s1np1/(8.0*XnP1) - (zeta2*s1np1)/(2.0*XnP1) + s1np1*s1np1/(8.0*XnP1*XnP1) + s2nm1 - (zeta2*s2nm1)/2.0 + (s1nm1*s1nm1*s2nm1)/2.0 - (3.0*s2n)/(8.0*x*x) - (3.0*s2n)/(4.0*x) + (s1n*s2n)/(2.0*x) - (3.0*s2np1)/(8.0*XnP1*XnP1) + (3.0*s2np1)/(4.0*XnP1) + (s1np1*s2np1)/(2.0*XnP1) - s211nm1/2.0 - (3.0*s3nm1)/8.0 + (3.0*s1nm1*s3nm1)/2.0 + (3.0*s3n)/(4.0*x) + (3.0*s3np1)/(4.0*XnP1) - s31nm1/2.0 + s4nm1);
    C2NSqqb = (-Ca/2.0+Cf)*Cf/24.0*(168.0 + 4.0*x*(180.0 + 3.0*x*(97.0+ x*(79.0 + x*(32.0 + x*(7.0+x)))) - 2.0*x*XnP1*XnP1*XnP1*(1.0+2.0*x)*pi*pi) - (x*x*XnP1*XnP1)*(192.0*x*x*XnP1*XnP1*acg13xnm1+16.0*x*x*XnP1*XnP1*pi*pi*acg1xnm1 + (48.0)*(ps1nf2-ps1n)+x*(96.0*XnP1*XnP1*acg1pxn - 96.0*x*XnP1*acg1pxnp1 - 48.0*x*acg1ppxnm1 - 96.0*x*x*acg1ppxnm1 - 48.0*x*x*x*acg1ppxnm1 + 96.0*x*acg2pxnm1 + 192.0*x*x*acg2pxnm1 + 96.0*x*x*x*acg2pxnm1 + 96.0*acg4xn + 192.0*x*acg4xn + 96.0*x*x*acg4xn - 96.0*x*acg4xnp1 - 96.0*x*x*acg4xnp1 + 96.0*x*acg4pxnm1 + 192.0*x*x*acg4pxnm1 + 96.0*x*x*x*acg4pxnm1 + 96.0*x*acg5xnm1 + 192.0*x*x*acg5xnm1 + 96.0*x*x*x*acg5xnm1 - 192.0*x*acg6xnm1 - 384.0*x*x*acg6xnm1 - 192.0*x*x*x*acg6xnm1 - 288.0*x*acg7xnm1 - 576.0*x*x*acg7xnm1 - 288.0*x*x*x*acg7xnm1 + 192.0*x*acg9xnm1 + 384.0*x*x*acg9xnm1 + 192.0*x*x*x*acg9xnm1 - 4.0*x*XnP1*pi*pi*ps0nf2p3f2 - 12.0*(-6.0+(x-3.0)*x)*ps1nf2 + 48.0*x*XnP1*ps1n - 12.0*x*(x+3.0)*ps1nf2p3f2 - 3.0*XnP1*(2.0+3.0*x)*ps2nf2 + 24.0*XnP1*XnP1*ps2n + 3.0*x*XnP1*ps2nf2p3f2 - x*XnP1*XnP1*ps3nf2 + 8.0*x*XnP1*XnP1*ps3n + 8.0*XnP1*XnP1*log2comp*(pi*pi+6.0*x*zeta3) - 8.0*XnP1*XnP1*ps0n*(pi*pi+6.0*x*zeta3)+4.0*XnP1*ps0nf2*((2.0+3.0*x)*pi*pi + 12.0*x*XnP1*zeta3))))/(4.0*x*x*x*x*XnP1*XnP1*XnP1*XnP1);
   C2qg= Cf*(1.0/(16.0*x*x*x*x) + 1.0/(32.0*x*x*x) - 1.0/(4.0*x*x) -13.0/(32.0*x) - 1.0/(8.0*XnP1*XnP1*XnP1*XnP1) + 3.0/(8.0*XnP1*XnP1*XnP1)-15.0/(32.0*XnP1*XnP1) + 43.0/(32.0*XnP1) + 1.0/(4.0*XnP2*XnP2*XnP2*XnP2) - 1.0/(4.0*XnP2*XnP2*XnP2) + 1.0/(4.0*XnP2*XnP2) - 5.0/(4.0*XnP2) + (3.0*zeta2)/(4.0*XnP1) - (3.0*zeta2)/(4.0*XnP2) + zeta3/x - (2.0*zeta3)/XnP1 + (2.0*zeta3)/XnP2 - s1n*s1n/(16.0*x*x) + s1n*s1n*s1n/(48.0*x) - s1np1/(4.0*XnP1*XnP1) + (3.0*s1np1)/(16.0*XnP1) + s1np1*s1np1/(8.0*XnP1*XnP1) + s1np1*s1np1/(8.0*XnP1) - s1np1*s1np1*s1np1/(24.0*XnP1) + s1np2/(4.0*XnP2*XnP2) - s1np2/(4.0*XnP2) - s1np2*s1np2/(8.0*XnP2*XnP2) - s1np2*s1np2/(8.0*XnP2) + s1np2*s1np2*s1np2/(24.0*XnP2) - s2n/(16.0*x*x) + (s1n*s2n)/(16.0*x) + s2np1/(8.0*XnP1*XnP1) - s2np1/(8.0*XnP1) - (s1np1*s2np1)/(8.0*XnP1) - s2np2/(8.0*XnP2*XnP2) + s2np2/(8.0*XnP2) + (s1np2*s2np2)/(8.0*XnP2) - s3n/(12.0*x) + s3np1/(6.0*XnP1) - s3np2/(6.0*XnP2)) + Ca*(43.0/(108.0*XnM1) - 1.0/(8.0*x*x*x*x) - 1.0/(16.0*x*x*x) - 7.0/(24.0*x*x) - 35.0/(48.0*x) - 1.0/(4.0*XnP1*XnP1*XnP1*XnP1) + 1.0/(4.0*XnP1*XnP1*XnP1) + 5.0/(12.0*XnP1*XnP1) + 43.0/(48.0*XnP1) - 11.0/(12.0*XnP2*XnP2*XnP2) - 17.0/(18.0*XnP2*XnP2) - 149.0/(216.0*XnP2) - zeta2/(6.0*XnM1) + zeta2/(4.0*x*x) + zeta2/(4.0*x) - 6.0*zeta2/(8.0*XnP1) + zeta2/(2.0*XnP2*XnP2) + 16.0*zeta2/(24.0*XnP2) - 4.0*zeta3/(16.0*x) + 4.0*zeta3/(8.0*XnP1) - 4.0*zeta3/(8.0*XnP2) - s1n*s1n*s1n/(48.0*x) - 3.0*s1np1/(16.0*XnP1) - s1np1*s1np1/(8.0*XnP1) + s1np1*s1np1*s1np1/(24.0*XnP1) + s1np2/(4.0*XnP2) + s1np2*s1np2/(8.0*XnP2) - s1np2*s1np2*s1np2/(24.0*XnP2) + s2nm1f2/(8.0*x*x) - s1n*s2nm1f2/(8.0*x) + s2nm1/(6.0*XnM1) - s2n/(4.0*x*x) - s2n/(4.0*x) + s1n*s2n/(16.0*x)- s2np1f2/(4.0*XnP1*XnP1) + 2.0*s2np1f2/(16.0*XnP1) + s2np1f2/(4.0*XnP2*XnP2) - 2.0*s2np1f2/(16.0*XnP2) + s1np1*s2np1f2/(4.0*XnP1) - s1np2*s2np1f2/(4.0*XnP2) + 5.0*s2np1/(8.0*XnP1) - s1np1*s2np1/(8.0*XnP1) -s2np2/(2.0*XnP2*XnP2) - 13.0*s2np2/(24.0*XnP2) + s1np2*s2np2/(8.0*XnP2) - 2.0*s3nm1f2/(64.0*x) -s3n/(24.0*x) + 2.0*s3np1f2/(32.0*XnP1) - 2.0*s3np1f2/(32.0*XnP2) + s3np1/(12.0*XnP1) - s3np2/(12.0*XnP2) + sm21n/(4.0*x) - sm21np1/(2.0*XnP1) + sm21np2/(2.0*XnP2));
   C2qg=4.0*C2qg;
   C2NSqqb=4.0*C2NSqqb;
   C2Sqqb=4.0*C2Sqqb;
   C2NSqq=4.0*C2NSqq;
   c.xp = xp;
   c.C2qgpos = C2qg;
   c.C2NSqqbpos = C2NSqqb;
   c.C2Sqqbpos = C2Sqqb;
   c.C2NSqqpos = C2NSqq;

   //Now negative branch
   x = xm;
   XnP1 = x + OneComplex;
   XnP2 = x + 2.0*OneComplex;
   XnP3 = x + 3.0*OneComplex;
   XnM1 = x - OneComplex;
   XnO2P1 = x/2.0+OneComplex;
   XnO2P2 = x/2.0+2.0*OneComplex;
   XnO2P3O2 = x/2.0+3.0*OneComplex/2.0;
   XnO2P1O2 = x/2.0+OneComplex/2.0;
   XnO2 = x/2.0;
   ps0np1 = psi0(XnP1, verbosity);
   ps0np2 = psi0(XnP2, verbosity);
   ps0np3 = psi0(XnP3, verbosity);
   ps0n = psi0(x, verbosity);
   ps0nf2p1 = psi0(XnO2P1, verbosity);
   ps0nf2p3f2 = psi0(XnO2P3O2, verbosity);
   ps0nf2p1f2 = psi0(XnO2P1O2, verbosity);
   ps0nf2 = psi0(XnO2, verbosity);
   ps0nf2p2 = psi0(XnO2P2, verbosity);
   ps1np1 = psi1(XnP1, verbosity);
   ps1np2 = psi1(XnP2, verbosity);
   ps1np3 = psi1(XnP3, verbosity);
   ps1n = psi1(x, verbosity);
   ps1nf2p1 = psi1(XnO2P1, verbosity);
   ps1nf2p3f2 = psi1(XnO2P3O2, verbosity);
   ps1nf2p1f2 = psi1(XnO2P1O2, verbosity);
   ps1nf2 = psi1(XnO2, verbosity);
   ps1nf2p2 = psi1(XnO2P2, verbosity);
   ps2np1 = psi2(XnP1, verbosity);
   ps2np2 = psi2(XnP2, verbosity);
   ps2np3 = psi2(XnP3, verbosity);
   ps2n = psi2(x, verbosity);
   ps2nf2p1 = psi2(XnO2P1, verbosity);
   ps2nf2p3f2 = psi2(XnO2P3O2, verbosity);
   ps2nf2p1f2 = psi2(XnO2P1O2, verbosity);
   ps2nf2 = psi2(XnO2, verbosity);
   ps2nf2p2 = psi2(XnO2P2, verbosity);
   ps3n = psi3(x, verbosity);
   ps3nf2p1f2 = psi3(XnO2P1O2, verbosity);
   ps3nf2 = psi3(XnO2, verbosity);
   beta0n = bet(x, verbosity);
   beta0np1 = bet(XnP1, verbosity);
   beta0np2 = bet(XnP2, verbosity);
   beta0np3 = bet(XnP3, verbosity);
   beta1n = bet1(x, verbosity);
   beta2n = bet2(x, verbosity);
   beta3n = bet3(x, verbosity);
   s1n = ps0np1 + Eulerconst;
   s1np1 = ps0np2 + Eulerconst;
   s1np2 = ps0np3 + Eulerconst;
   s1nm1 = ps0n + Eulerconst;
   s1nf2 = ps0nf2p1 + Eulerconst;
   s1nm1f2 = ps0nf2p1f2 + Eulerconst;
   s1np1f2 = ps0nf2p3f2 + Eulerconst;
   s1nm2f2 = ps0nf2 + Eulerconst;
   s1np2f2 = ps0nf2p2 + Eulerconst;
   s2n = -ps1np1 + zeta2;
   s2np1 = -ps1np2 + zeta2;
   s2np2 = -ps1np3 + zeta2;
   s2nm1 = -ps1n + zeta2;
   s2nf2 = -ps1nf2p1 + zeta2;
   s2nm1f2 = -ps1nf2p1f2 + zeta2;
   s2np1f2 = -ps1nf2p3f2 + zeta2;
   s2nm2f2 = -ps1nf2 + zeta2;
   s2np2f2 = -ps1nf2p2 + zeta2;
   s3n = ps2np1/2.0 + zeta3;
   s3np1 = ps2np2/2.0 + zeta3;
   s3np2 = ps2np3/2.0 + zeta3;
   s3nm1 = ps2n/2.0 + zeta3;
   s3nf2 = ps2nf2p1/2.0 + zeta3;
   s3nm1f2 = ps2nf2p1f2/2.0 + zeta3;
   s3np1f2 = ps2nf2p3f2/2.0 + zeta3;
   s3nm2f2 = ps2nf2/2.0 + zeta3;
   s3np2f2 = ps2nf2p2/2.0 + zeta3;
   s4nm1 = -ps3n/6.0 + zeta4;
   s4nm1f2 = -ps3nf2p1f2/6.0 + zeta4;
   s4nm2f2 = -ps3nf2/6.0 + zeta4;
   s211nm1 = -acg21(XnM1, verbosity)+6.0/5.0*zeta2*zeta2;
   s31nm1 = acg20(XnM1, verbosity)+zeta2*s2nm1-0.5*zeta2*zeta2;
   sm1nm1 = beta0n-ln2;
   sm1n = -beta0np1-ln2;
   sm1np1 = beta0np2-ln2;
   sm1np2 = -beta0np3 - ln2;
   sm2nm1 = -beta1n - zeta2/2.0;
   sm3nm1 = beta2n/2.0 - zeta3*0.75;
   sm4nm1 = -beta3n/6.0-zeta4*7.0/8.0;
   sm21n = acg3(x, verbosity) + zeta2*sm1n - 5.0/8.0*zeta3+zeta2*ln2;
   sm21np1 = -acg3(XnP1, verbosity) + zeta2*sm1np1 - 5.0/8.0*zeta3+zeta2*ln2;
   sm21np2 = acg3(XnP2, verbosity) + zeta2*sm1np2 - 5.0/8.0*zeta3+zeta2*ln2;
   sm31nm1 = acg6(XnM1, verbosity) + zeta2*sm2nm1 - zeta3*sm1nm1 - 3.0/5.0*zeta2*zeta2 + 2*plog41f2 + 0.75*zeta3*ln2 - zeta2/2.0*ln2*ln2 + ln2*ln2*ln2*ln2/12.0;
   sm22nm1 = acg5(XnM1, verbosity) -2.0*sm31nm1 + 2.0*zeta2*sm2nm1 + 0.075*zeta2*zeta2;
   sm211nm1 = -acg8(XnM1, verbosity) + zeta3*sm1nm1 - plog41f2 + zeta2*zeta2/8.0 + zeta3*ln2/8.0 + zeta2*ln2*ln2/4.0 - ln2*ln2*ln2*ln2/24.0;
   acg13xnm1 = acg13(XnM1, verbosity);
   acg1xnm1 = acg1(XnM1, verbosity);
   acg1pxn = acg1p(x, verbosity);
   acg1pxnp1 = acg1p(XnP1, verbosity);
   acg1ppxnm1 = acg1pp(XnM1, verbosity);
   acg2pxnm1 = acg2p(XnM1, verbosity);
   acg4xn = acg4(x, verbosity);
   acg4xnp1 = acg4(XnP1, verbosity);
   acg4pxnm1 = acg4p(XnM1, verbosity);
   acg5xnm1 = acg5(XnM1, verbosity);
   acg6xnm1 = acg6(XnM1, verbosity);
   acg7xnm1 = acg7(XnM1, verbosity);
   acg9xnm1 = acg9(XnM1, verbosity);
   C2Sqqb = std::complex<double>(0.0,0.0);
   C2NSqq = std::complex<double>(0.0,0.0);
   C2NSqqb = std::complex<double>(0.0,0.0);
   C2qg = std::complex<double>(0.0,0.0);
   C2Sqqb = Cf*(172.0/XnM1-54.0/(x*x*x*x) - 27.0/(x*x*x) - 126.0/(x*x) - 315.0/x - 54.0/(XnP1*XnP1*XnP1*XnP1) - 27.0/(XnP1*XnP1*XnP1) + 180.0/(XnP1*XnP1) + 279.0/XnP1 - 72.0/(XnP2*XnP2*XnP2) - 192.0/(XnP2*XnP2) - 136.0/XnP2 - (72.0*zeta2)/XnM1 + (108.0*zeta2)/x - (108.0*zeta2)/XnP1 + (72.0*zeta2)/XnP2 + (72.0*s2nm1)/XnM1 - (108.0*s2n)/x + (108.0*s2np1)/XnP1 - (72.0*s2np2)/XnP2)/432.0;
   C2NSqq = Cf*Nf*(127.0/192.0 - 1.0/(24.0*x*x*x) + 5.0/(72.0*x*x) - 37.0/(216.0*x) - 1.0/(24.0*XnP1*XnP1*XnP1) + 5.0/(72.0*XnP1*XnP1) - 19.0/(216.0*XnP1) - (2.0*zeta2)/3.0 + (7.0*zeta3)/36.0 - (7.0*s1nm1)/27.0 + (5.0*s2nm1)/36.0 -s3nm1/12.0) + Ca*Cf*(-1535.0/384.0 - 1.0/(8.0*x*x*x*x) + 11.0/(48.0*x*x*x) - 47.0/(72.0*x*x) + 50.0/(27.0*x) - 1.0/(8.0*XnP1*XnP1*XnP1*XnP1) - 1.0/(48.0*XnP1*XnP1*XnP1)-83.0/(72.0*XnP1*XnP1) + 1.0/(54.0*XnP1) + (97.0*zeta2)/24.0 - (3.0*zeta2)/(8.0*x) + (3.0*zeta2)/(8.0*XnP1) + (85.0*zeta3)/72.0 - (5.0*zeta3)/(8.0*x) - (5.0*zeta3)/(8.0*XnP1) - (29.0*zeta4)/16.0 + (101.0*s1nm1)/54.0 - (5.0*zeta3*s1nm1)/4.0 + (zeta2*s1nm1*s1nm1)/4.0 + (zeta2*s1n)/(4.0*x) - s1np1/(8.0*XnP1) + (zeta2*s1np1)/(4.0*XnP1) - (19.0*s2nm1)/18.0 + (zeta2*s2nm1)/4.0 - (s1nm1*s1nm1*s2nm1)/4.0 + s2n/(4.0*x*x) + s2n/(4.0*x) - (s1n*s2n)/(4.0*x) + s2np1/(4.0*XnP1*XnP1) - s2np1/(4.0*XnP1) - (s1np1*s2np1)/(4.0*XnP1) + s211nm1/2.0 + (11.0*s3nm1)/24.0 - (s1nm1*s3nm1)/2.0 - s3n/(4.0*x) - s3np1/(4.0*XnP1) - s4nm1/2.0) + Cf*Cf*(255.0/128.0 - 1.0/(8.0*x*x*x*x) + 3.0/(8.0*x*x) -19.0/(8.0*x) - 1.0/(8.0*XnP1*XnP1*XnP1*XnP1) + 1.0/(4.0*XnP1*XnP1*XnP1) + 2.0/(XnP1*XnP1) + 19.0/(8.0*XnP1) - (35.0*zeta2)/16.0 + (3.0*zeta2)/(2.0*x) - (3.0*zeta2)/(2.0*XnP1) - (3.0*zeta3)/2.0 - (3.0*zeta3)/(4.0*x) - (3.0*zeta3)/(4.0*XnP1) + (61.0*zeta4)/16.0 - (3.0*zeta3*s1nm1)/2.0 - (zeta2*s1nm1*s1nm1)/2.0 - s1n/(4.0*x*x) - (zeta2*s1n)/(2.0*x) + s1n*s1n/(8.0*x*x) + s1np1/(4.0*XnP1*XnP1) + s1np1/(8.0*XnP1) - (zeta2*s1np1)/(2.0*XnP1) + s1np1*s1np1/(8.0*XnP1*XnP1) + s2nm1 - (zeta2*s2nm1)/2.0 + (s1nm1*s1nm1*s2nm1)/2.0 - (3.0*s2n)/(8.0*x*x) - (3.0*s2n)/(4.0*x) + (s1n*s2n)/(2.0*x) - (3.0*s2np1)/(8.0*XnP1*XnP1) + (3.0*s2np1)/(4.0*XnP1) + (s1np1*s2np1)/(2.0*XnP1) - s211nm1/2.0 - (3.0*s3nm1)/8.0 + (3.0*s1nm1*s3nm1)/2.0 + (3.0*s3n)/(4.0*x) + (3.0*s3np1)/(4.0*XnP1) - s31nm1/2.0 + s4nm1);
   C2NSqqb = (-Ca/2.0+Cf)*Cf/24.0*(168.0 + 4.0*x*(180.0 + 3.0*x*(97.0+ x*(79.0 + x*(32.0 + x*(7.0+x)))) - 2.0*x*XnP1*XnP1*XnP1*(1.0+2.0*x)*pi*pi) - (x*x*XnP1*XnP1)*(192.0*x*x*XnP1*XnP1*acg13xnm1+16.0*x*x*XnP1*XnP1*pi*pi*acg1xnm1 + (48.0)*(ps1nf2-ps1n)+x*(96.0*XnP1*XnP1*acg1pxn - 96.0*x*XnP1*acg1pxnp1 - 48.0*x*acg1ppxnm1 - 96.0*x*x*acg1ppxnm1 - 48.0*x*x*x*acg1ppxnm1 + 96.0*x*acg2pxnm1 + 192.0*x*x*acg2pxnm1 + 96.0*x*x*x*acg2pxnm1 + 96.0*acg4xn + 192.0*x*acg4xn + 96.0*x*x*acg4xn - 96.0*x*acg4xnp1 - 96.0*x*x*acg4xnp1 + 96.0*x*acg4pxnm1 + 192.0*x*x*acg4pxnm1 + 96.0*x*x*x*acg4pxnm1 + 96.0*x*acg5xnm1 + 192.0*x*x*acg5xnm1 + 96.0*x*x*x*acg5xnm1 - 192.0*x*acg6xnm1 - 384.0*x*x*acg6xnm1 - 192.0*x*x*x*acg6xnm1 - 288.0*x*acg7xnm1 - 576.0*x*x*acg7xnm1 - 288.0*x*x*x*acg7xnm1 + 192.0*x*acg9xnm1 + 384.0*x*x*acg9xnm1 + 192.0*x*x*x*acg9xnm1 - 4.0*x*XnP1*pi*pi*ps0nf2p3f2 - 12.0*(-6.0+(x-3.0)*x)*ps1nf2 + 48.0*x*XnP1*ps1n - 12.0*x*(x+3.0)*ps1nf2p3f2 - 3.0*XnP1*(2.0+3.0*x)*ps2nf2 + 24.0*XnP1*XnP1*ps2n + 3.0*x*XnP1*ps2nf2p3f2 - x*XnP1*XnP1*ps3nf2 + 8.0*x*XnP1*XnP1*ps3n + 8.0*XnP1*XnP1*log2comp*(pi*pi+6.0*x*zeta3) - 8.0*XnP1*XnP1*ps0n*(pi*pi+6.0*x*zeta3)+4.0*XnP1*ps0nf2*((2.0+3.0*x)*pi*pi + 12.0*x*XnP1*zeta3))))/(4.0*x*x*x*x*XnP1*XnP1*XnP1*XnP1);
   C2qg= Cf*(1.0/(16.0*x*x*x*x) + 1.0/(32.0*x*x*x) - 1.0/(4.0*x*x) -13.0/(32.0*x) - 1.0/(8.0*XnP1*XnP1*XnP1*XnP1) + 3.0/(8.0*XnP1*XnP1*XnP1)-15.0/(32.0*XnP1*XnP1) + 43.0/(32.0*XnP1) + 1.0/(4.0*XnP2*XnP2*XnP2*XnP2) - 1.0/(4.0*XnP2*XnP2*XnP2) + 1.0/(4.0*XnP2*XnP2) - 5.0/(4.0*XnP2) + (3.0*zeta2)/(4.0*XnP1) - (3.0*zeta2)/(4.0*XnP2) + zeta3/x - (2.0*zeta3)/XnP1 + (2.0*zeta3)/XnP2 - s1n*s1n/(16.0*x*x) + s1n*s1n*s1n/(48.0*x) - s1np1/(4.0*XnP1*XnP1) + (3.0*s1np1)/(16.0*XnP1) + s1np1*s1np1/(8.0*XnP1*XnP1) + s1np1*s1np1/(8.0*XnP1) - s1np1*s1np1*s1np1/(24.0*XnP1) + s1np2/(4.0*XnP2*XnP2) - s1np2/(4.0*XnP2) - s1np2*s1np2/(8.0*XnP2*XnP2) - s1np2*s1np2/(8.0*XnP2) + s1np2*s1np2*s1np2/(24.0*XnP2) - s2n/(16.0*x*x) + (s1n*s2n)/(16.0*x) + s2np1/(8.0*XnP1*XnP1) - s2np1/(8.0*XnP1) - (s1np1*s2np1)/(8.0*XnP1) - s2np2/(8.0*XnP2*XnP2) + s2np2/(8.0*XnP2) + (s1np2*s2np2)/(8.0*XnP2) - s3n/(12.0*x) + s3np1/(6.0*XnP1) - s3np2/(6.0*XnP2)) + Ca*(43.0/(108.0*XnM1) - 1.0/(8.0*x*x*x*x) - 1.0/(16.0*x*x*x) - 7.0/(24.0*x*x) - 35.0/(48.0*x) - 1.0/(4.0*XnP1*XnP1*XnP1*XnP1) + 1.0/(4.0*XnP1*XnP1*XnP1) + 5.0/(12.0*XnP1*XnP1) + 43.0/(48.0*XnP1) - 11.0/(12.0*XnP2*XnP2*XnP2) - 17.0/(18.0*XnP2*XnP2) - 149.0/(216.0*XnP2) - zeta2/(6.0*XnM1) + zeta2/(4.0*x*x) + zeta2/(4.0*x) - 6.0*zeta2/(8.0*XnP1) + zeta2/(2.0*XnP2*XnP2) + 16.0*zeta2/(24.0*XnP2) - 4.0*zeta3/(16.0*x) + 4.0*zeta3/(8.0*XnP1) - 4.0*zeta3/(8.0*XnP2) - s1n*s1n*s1n/(48.0*x) - 3.0*s1np1/(16.0*XnP1) - s1np1*s1np1/(8.0*XnP1) + s1np1*s1np1*s1np1/(24.0*XnP1) + s1np2/(4.0*XnP2) + s1np2*s1np2/(8.0*XnP2) - s1np2*s1np2*s1np2/(24.0*XnP2) + s2nm1f2/(8.0*x*x) - s1n*s2nm1f2/(8.0*x) + s2nm1/(6.0*XnM1) - s2n/(4.0*x*x) - s2n/(4.0*x) + s1n*s2n/(16.0*x)- s2np1f2/(4.0*XnP1*XnP1) + 2.0*s2np1f2/(16.0*XnP1) + s2np1f2/(4.0*XnP2*XnP2) - 2.0*s2np1f2/(16.0*XnP2) + s1np1*s2np1f2/(4.0*XnP1) - s1np2*s2np1f2/(4.0*XnP2) + 5.0*s2np1/(8.0*XnP1) - s1np1*s2np1/(8.0*XnP1) -s2np2/(2.0*XnP2*XnP2) - 13.0*s2np2/(24.0*XnP2) + s1np2*s2np2/(8.0*XnP2) - 2.0*s3nm1f2/(64.0*x) -s3n/(24.0*x) + 2.0*s3np1f2/(32.0*XnP1) - 2.0*s3np1f2/(32.0*XnP2) + s3np1/(12.0*XnP1) - s3np2/(12.0*XnP2) + sm21n/(4.0*x) - sm21np1/(2.0*XnP1) + sm21np2/(2.0*XnP2));
   C2qg=4.0*C2qg;
   C2NSqqb=4.0*C2NSqqb;
   C2Sqqb=4.0*C2Sqqb;
   C2NSqq=4.0*C2NSqq;
   c.xm = xm;
   c.C2qgneg = C2qg;
   c.C2NSqqbneg = C2NSqqb;
   c.C2Sqqbneg = C2Sqqb;
   c.C2NSqqneg = C2NSqq;

}



//////////////////////////////////////////////////////////////
// VARIOUS FUNCTIONS IN MELLIN SPACE -- ADAPTED FROM BLUEMLEIN
//////////////////////////////////////////////////////////////

//Psi function used in determination of C2 coefficients and of anomalous dimension
std::complex<double> psi(std::complex<double> z, int verbosity) {
    std::complex<double> OneComplex;
    OneComplex = std::complex<double>(1.0,0.0);
    std::complex<double> sub, zz, rz, dz, psi;
    sub = std::complex<double>(0.0,0.0);
    zz = std::complex<double>(0.0,0.0);
    rz = std::complex<double>(0.0,0.0);
    dz = std::complex<double>(0.0,0.0);
    psi = std::complex<double>(0.0,0.0);
    zz = z;
    while (real(zz)<10.0) {
	sub = sub -1.0/zz;
	zz = zz + OneComplex;
	// std::cout << "sub = " << sub << std::endl;
	// std::cout << "zz = " << zz << std::endl;
    }
    rz = 1.0/zz;
    dz = rz*rz;
    // std::cout << "rz = " << rz << std::endl;
    // std::cout << "dz = " << dz << std::endl;
    psi = sub + std::log(zz) -rz/2.0 - dz/2520.0*(210.0*OneComplex+dz*(-21.0*OneComplex+10.0*dz));
    // std::cout << "psi = " << psi << std::endl;
    // std::cout << "psi pieces: " << sub << " " << log(zz) << " " << -rz/2.0 << " " << -dz/2520*(210.0*OneComplex+dz*(-21.0*OneComplex+10.0*dz)) << std::endl;
    return psi;
}

//Psi first derivative function used in determination of C2 coefficients and of anomalous dimension
std::complex<double> psideriv1(std::complex<double> z, int verbosity) {
    std::complex<double> OneComplex;
    OneComplex = std::complex<double>(1.0,0.0);
    std::complex<double> sub, zz, rz, dz, psifirstderiv;
    sub = std::complex<double>(0.0,0.0);
    zz = std::complex<double>(0.0,0.0);
    rz = std::complex<double>(0.0,0.0);
    dz = std::complex<double>(0.0,0.0);
    psifirstderiv = std::complex<double>(0.0,0.0);
    zz = z;
    while (real(zz)<10.0) {
	sub = sub + 1.0/(zz*zz);
	zz = zz + OneComplex;
    }
    rz = 1.0/zz;
    dz = rz*rz;
    psifirstderiv = sub + rz +dz/2.0*(OneComplex + rz/630.0*(210.0*OneComplex - dz*(42.0*OneComplex - dz*(30.0*OneComplex - 42.0*dz))));

    return psifirstderiv;
}

//Psi second derivative function used in determination of C2 coefficients and of anomalous dimension
std::complex<double> psideriv2(std::complex<double> z, int verbosity) {
    std::complex<double> OneComplex;
    OneComplex = std::complex<double>(1.0,0.0);
    std::complex<double> sub, zz, rz, dz, psisecondderiv;
    sub = std::complex<double>(0.0,0.0);
    zz = std::complex<double>(0.0,0.0);
    rz = std::complex<double>(0.0,0.0);
    dz = std::complex<double>(0.0,0.0);
    psisecondderiv = std::complex<double>(0.0,0.0);
    zz = z;
    while (real(zz)<10.0) {
	sub = sub -2.0/(zz*zz*zz);
	zz = zz + OneComplex;
    }
    rz = 1.0/zz;
    dz = rz*rz;
    psisecondderiv = sub - dz/60.0*(60.0*OneComplex + rz*(60.0*OneComplex + rz*(30.0*OneComplex - dz*(10.0*OneComplex - dz*(10.0*OneComplex - dz*(18.0*OneComplex - 50.0*dz))))));
    return psisecondderiv;
}

//Anomalous Dimensions calculations up to dependence on number of active flavours Nf
std::complex<double> *anomcalc( std::complex<double> ancalc[15], std::complex<double> x , int verbosity) {
    std::complex<double> xsquared, xP1, xP2, xM1;
    std::complex<double> OneComplex;

    double zeta2 = k_constants::zeta2;
    double zeta3 = k_constants::zeta3;

    OneComplex = std::complex<double>(1.0,0.0);
    xsquared = std::complex<double>(0.0,0.0);
    xP1 = std::complex<double>(0.0,0.0);
    xP2 = std::complex<double>(0.0,0.0);
    xM1 = std::complex<double>(0.0,0.0);
    xsquared = x*x;
    xP1 = x + OneComplex;
    xP2 = x + 2.0*OneComplex;
    xM1 = x- OneComplex;
    std::complex<double> cpsi, qqi, qgf, gqi, ggi, ggf;
    cpsi = std::complex<double>(0.0,0.0);
    qqi = std::complex<double>(0.0,0.0);
    qgf = std::complex<double>(0.0,0.0);
    gqi = std::complex<double>(0.0,0.0);
    ggi = std::complex<double>(0.0,0.0);
    ggf = std::complex<double>(0.0,0.0);
    //Leading order
    cpsi = psi(xP1, verbosity) + 0.577216;
    qqi = (8.0/3.0)*(-3.0-2.0/(x*xP1) + 4.0*cpsi);
    qgf = -4.0*(xsquared+x+2.0*OneComplex)/(x*xP1*xP2);
    gqi = (-16.0/3.0)*(xsquared+x+2.0*OneComplex)/(x*xP1*xM1);
    ggi = -22.0 - 24.0/(x*xM1) - 24.0/(xP1*xP2) + 24.0*cpsi;
    ggf = 4.0/3.0;
    //Next to Leading order
    std::complex<double> xcubed, xto4, xP1squared, xP1cubed;
    xcubed = std::complex<double>(0.0,0.0);
    xto4 = std::complex<double>(0.0,0.0);
    xP1squared = std::complex<double>(0.0,0.0);
    xP1cubed = std::complex<double>(0.0,0.0);
    xcubed = x*x*x;
    xto4 = xcubed*x;
    xP1squared = xP1*xP1;
    xP1cubed = xP1*xP1*xP1;
    //Analytic continuations of N-sums
    std::complex<double> cpsi1, spmom, slc, slv,sschlm,sstr2m,sstr3m,sschlp,sstr2p,sstr3p;
    cpsi1 = std::complex<double>(0.0,0.0);
    spmom = std::complex<double>(0.0,0.0);
    slc = std::complex<double>(0.0,0.0);
    slv = std::complex<double>(0.0,0.0);
    sschlm = std::complex<double>(0.0,0.0);
    sstr2m = std::complex<double>(0.0,0.0);
    sstr3m = std::complex<double>(0.0,0.0);
    sschlp = std::complex<double>(0.0,0.0);
    sstr2p = std::complex<double>(0.0,0.0);
    sstr3p = std::complex<double>(0.0,0.0);
    cpsi1 = zeta2 - psideriv1(xP1, verbosity);
    spmom = 1.004/xP1 -0.846/xP2 + 1.342/(x+3.0*OneComplex) - 1.532/(x+4.0*OneComplex) + 0.839/(x+5.0*OneComplex);
    slc = -5.0/8.0*zeta3;
    slv = -zeta2/2.0*(psi(xP1/2.0, verbosity)-psi(x/2.0, verbosity)) + cpsi/xsquared + spmom;
    sschlm = slc - slv;
    sstr2m = zeta2 - psideriv1(xP1/2.0, verbosity);
    sstr3m = 0.5*psideriv2(xP1/2.0, verbosity) + zeta3;
    sschlp = slc + slv;
    sstr2p = zeta2 - psideriv1(xP2/2.0, verbosity);
    sstr3p = 0.5*psideriv2(xP2/2.0, verbosity)+zeta3;
    //non-singlet pieces
    std::complex<double> ns1ma, ns1pa, ns1b, ns1c, ns1mi, ns1pi, ns1f;
    ns1ma = std::complex<double>(0.0,0.0);
    ns1pa = std::complex<double>(0.0,0.0);
    ns1c = std::complex<double>(0.0,0.0);
    ns1mi = std::complex<double>(0.0,0.0);
    ns1pi = std::complex<double>(0.0,0.0);
    ns1f = std::complex<double>(0.0,0.0);
    ns1ma = 16.0*cpsi*(2.0*x + OneComplex)/(xsquared*xP1squared) + 16.0*(2.0*cpsi -1.0/(x*xP1))*(cpsi1-sstr2m) + 64.0*sschlm + 24.0*cpsi1 -3.0*OneComplex - 8.0*sstr3m - 8.0*(3.0*xcubed + xsquared - OneComplex)/(xcubed*xP1cubed) + 16.0*(2.0*xsquared + 2.0*x + OneComplex)/(xcubed*xP1cubed);
    ns1pa = 16.0*cpsi*(2.0*x + OneComplex)/(xsquared*xP1squared) + 16.0*(2.0*cpsi -1.0/(x*xP1))*(cpsi1-sstr2p) + 64.0*sschlp + 24.0*cpsi1 -3.0*OneComplex - 8.0*sstr3p - 8.0*(3.0*xcubed + xsquared - OneComplex)/(xcubed*xP1cubed) - 16.0*(2.0*xsquared + 2.0*x + OneComplex)/(xcubed*xP1cubed);
    ns1b = cpsi*(536.0/9.0 + 8.0*(2.0*x + OneComplex)/(xsquared*xP1squared)) - (16.0*cpsi + 52.0/3.0 -8.0/(x*xP1))*cpsi1 - 43.0/6.0 - (151.*xto4 + 263.0*xcubed + 97.0*xsquared + 3.0*x +9.0*OneComplex)*4.0/(9.0*xcubed*xP1cubed);
    ns1c = -160.0/9.0*cpsi + 32.0/3.0*cpsi1 + 4.0/3.0 + 16.0*(11.0*xsquared + 5.0*x -3.0*OneComplex)/(9.0*xsquared*xP1squared);
    ns1mi = -2.0/9.0*ns1ma + 4.0*ns1b;
    ns1pi = -2.0/9.0*ns1pa + 4.0*ns1b;
    ns1f = 2.0/3.0*ns1c;
    //Singlet pieces
    std::complex<double> xto5, xto6, xto7, xto8, xto9, xM1squared, xP2squared, xP2cubed;
    xto5 = std::complex<double>(0.0,0.0);
    xto6 = std::complex<double>(0.0,0.0);
    xto7 = std::complex<double>(0.0,0.0);
    xto8 = std::complex<double>(0.0,0.0);
    xto9 = std::complex<double>(0.0,0.0);
    xM1squared = std::complex<double>(0.0,0.0);
    xP2squared = std::complex<double>(0.0,0.0);
    xP2cubed = std::complex<double>(0.0,0.0);
    xto5 = x*x*x*x*x;
    xto6 = xto5*x;
    xto7 = xto6*x;
    xto8 = xto7*x;
    xto9 = xto8*x;
    xM1squared = xM1*xM1;
    xP2squared = xP2*xP2;
    xP2cubed = xP2*xP2*xP2;
    std::complex<double> qq1f, qg1a, qg1b, qg1f, gq1a, gq1b, gq1c, gq1i, gq1f, gg1a, gg1b, gg1c, gg1i, gg1f;
    qq1f = std::complex<double>(0.0,0.0);
    qg1a = std::complex<double>(0.0,0.0);
    qg1b = std::complex<double>(0.0,0.0);
    qg1f = std::complex<double>(0.0,0.0);
    gq1a = std::complex<double>(0.0,0.0);
    gq1b = std::complex<double>(0.0,0.0);
    gq1c = std::complex<double>(0.0,0.0);
    gq1i = std::complex<double>(0.0,0.0);
    gq1f = std::complex<double>(0.0,0.0);
    gg1a = std::complex<double>(0.0,0.0);
    gg1b = std::complex<double>(0.0,0.0);
    gg1i = std::complex<double>(0.0,0.0);
    gg1f = std::complex<double>(0.0,0.0);
    qq1f = (5.0*xto5+32.0*xto4 + 49.0*xcubed + 38.0*xsquared + 28.0*x + 8.0*OneComplex)/(xM1*xcubed*xP1cubed*xP2squared)*(-32.0/3.0);
    qg1a = (-2.0*cpsi*cpsi + 2.0*cpsi1 -2.0*sstr2p)*(xsquared+x+2.0*OneComplex)/(x*xP1*xP2) + (8.0*cpsi*(2.0*x+3.0*OneComplex))/(xP1squared*xP2squared) + 2.0*(xto9+6.0*xto8 + 15.0*xto7 + 25.0*xto6 + 36.0*xto5 + 85.0*xto4 + 128.0*xcubed + 104.0*xsquared + 64.0*x + 16.0*OneComplex)/(xM1*xcubed*xP1cubed*xP2cubed);
    qg1b = (2.0*cpsi*cpsi - 2.0*cpsi1 + 5.0)*(xsquared + x + 2.0*OneComplex)/(x*xP1*xP2) - 4.0*cpsi/xsquared + (11.0*xto4 + 26.0*xcubed + 15.0*xsquared + 8.0*x + 4.0*OneComplex)/(xcubed*xP1cubed*xP2);
    qg1f = -12.0*qg1a - 16.0/3.0*qg1b;
    gq1a = (-2.0*cpsi*cpsi + 10.0*cpsi - 2.0*cpsi1)*(xsquared+x+2.0*OneComplex)/(xM1*x*xP1) - 4.0*cpsi/xP1squared - (12.0*xto6 + 30.0*xto5 + 43.0*xto4 + 28.0*xcubed - xsquared - 12.0*x -4.0*OneComplex)/(xM1*xcubed*xP1cubed);
    gq1b = (cpsi*cpsi + cpsi1 - sstr2p)*(xsquared + x + 2.0*OneComplex)/(xM1*x*xP1) - cpsi*(17.0*xto4 + 41.0*xsquared - 22.0*x - 12.0)/(3.0*xM1squared*xsquared*xP1) + (109.0*xto9 + 621.0*xto8 + 1400.0*xto7 + 1678.0*xto6 + 695.0*xto5 -1031.0*xto4 - 1304.0*xcubed - 152.0*xsquared + 432.0*x + 144.0)/(9.0*xM1squared*xcubed*xP1cubed*xP2squared);
    gq1c = (cpsi-8.0/3.0)*(xsquared + x + 2.0*OneComplex)/(xM1*x*xP1) + 1.0/xP1squared;
    gq1i = -64.0/9.0*gq1a - 32.0*gq1b;
    gq1f = -64.0/9.0*gq1c;
    gg1a = 16.0/9.0*(38.0*xto4 + 76.0*xcubed + 94.0*xsquared + 56.0*x + 12.0 )/(xM1*xsquared*xP1squared*xP2) - 160.0/9.0*cpsi + 32.0/3.0;
    gg1b = (2.0*xto6 + 4.0*xto5 + xto4 -10.0*xcubed - 5.0*xsquared -4.0*x - 4.0*OneComplex)*16.0/(xM1*xcubed*xP1cubed*xP2) + 8.0*OneComplex;
    gg1c = (2.0*xto5 + 5.0*xto4 + 8.0*xcubed + 7.0*xsquared -2.0*x -2.0*OneComplex)*64.0*cpsi/(xM1squared*xsquared*xP1squared*xP2squared) + 536.0/9.0*cpsi -64.0/3.0 + 32.0*sstr2p*(xsquared+x+OneComplex)/(xM1*x*xP1*xP2) - 16.0*cpsi*sstr2p + 32.0*sschlp - 4.0*sstr3p - 4.0*(457.0*xto9 + 2742.0*xto8 + 6040.0*xto7 + 6098.0*xto6 + 1567.0*xto5 - 2344.0*xto4 - 1632.0*xcubed + 560.0*xsquared + 1488.0*x + 576.0*OneComplex)/(9.0*xM1squared*xcubed*xP1cubed*xP2cubed);
    gg1i = 9.0*gg1c;
    gg1f = 3.0/2.0*gg1a + 2.0/3.0*gg1b;
    for (int i=0; i<15; i++) {ancalc[i] = std::complex<double> (0.0,0.0);}
    ancalc[1] = qqi; ancalc[2] = qgf; ancalc[3] = gqi; ancalc[4] = ggi; ancalc[5] = ggf; ancalc[6] = ns1mi; ancalc[7] = ns1pi, ancalc[8] = ns1f; ancalc[9] = qq1f; ancalc[10] = qg1f; ancalc[11] = gq1i; ancalc[12] = gq1f; ancalc[13] = gg1i; ancalc[14] = gg1f; ancalc[0] = x;

    return ancalc;
}


//psi0 function as needed in C2 coefficient calculations
std::complex<double> psi0 (std::complex<double> x, int verbosity) {
  std::complex<double> OneComplex;
  OneComplex = std::complex<double>(1.0,0.0);
  std::complex<double> zz,t, y, y2, t0, psi0;
  zz = x;
  t = std::complex<double>(0.0,0.0);
  y = std::complex<double>(0.0,0.0);
  y2 = std::complex<double>(0.0,0.0);
  t0 = std::complex<double>(0.0,0.0);
  double r = 0;
  r = std::abs(zz);
  while (r<=10) {
      t = t - OneComplex/zz;
      zz = zz + OneComplex;
      r = std::abs(zz);
  }
  y = OneComplex/zz;
  // y2 = (y*y).real(); ///NOTE I TAKE REAL PART HERE TO MATCH BLUEMLEIN.f FOR TESTING BUT ULTIMATELY SHOULDN'T BE REAL PART (THIS WAS A BUG IN BLUEMLEIN.f) AS IS A POLYGAMMA FUNCTION
  y2 = y*y;
  t0 = (-0.5 + (-1.0/12.0 + (1.0/12.0+(-1.0/252.0 + (1.0/240.0 + (-1.0/132.0 + (691.0/32760.0 + (-OneComplex/12.0 + OneComplex*3617.0/8160.0*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y)*y - std::log(y); //Check this they seem to have commented it out after 240.0 in blumlein.f
  psi0 = t + t0;
  return psi0;
}

//psi1 function as needed in C2 coefficient calculations
std::complex<double> psi1 (std::complex<double> x, int verbosity) {
  std::complex<double> OneComplex;
  OneComplex = std::complex<double>(1.0,0.0);
  std::complex<double> zz,t, y, y2, t0, psi1;
  zz = x;
  t = std::complex<double>(0.0,0.0);
  y = std::complex<double>(0.0,0.0);
  y2 = std::complex<double>(0.0,0.0);
  t0 = std::complex<double>(0.0,0.0);
  double r = 0;
  r = std::abs(zz);
  while (r<=10) {
    t = t + OneComplex/(zz*zz);
    zz = zz + OneComplex;
    r = std::abs(zz);
  }
  y = OneComplex/zz;
  // y2 = (y*y).real(); ///NOTE I TAKE REAL PART HERE TO MATCH BLUEMLEIN.f FOR TESTING BUT ULTIMATELY SHOULDN'T BE REAL PART (THIS WAS A BUG IN BLUEMLEIN.f) AS IS A POLYGAMMA FUNCTION
  y2 = y*y;
  t0 = (1.0 + (0.5 + (1.0/6.0 + (-1.0/30.0 + (1.0/42.0 + (-1.0/30.0 + (5.0/66.0 - 691.0/2730.0*y2)*y2)*y2)*y2)*y2)*y)*y)*y;
  psi1 = t + t0;
  return psi1;
}

//psi2 function as needed in C2 coefficient calculations
std::complex<double> psi2 (std::complex<double> x, int verbosity) {
  std::complex<double> OneComplex;
  OneComplex = std::complex<double>(1.0,0.0);
  std::complex<double> zz,t, y, y2, t0, psi2;
  zz = x;
  t = std::complex<double>(0.0,0.0);
  y = std::complex<double>(0.0,0.0);
  y2 = std::complex<double>(0.0,0.0);
  t0 = std::complex<double>(0.0,0.0);
  double r = 0;
  r = std::abs(zz);
  while (r<=10) {
    t = t - 2.0*OneComplex/(zz*zz*zz);
    zz = zz + OneComplex;
    r = std::abs(zz);
  }
  y = OneComplex/zz;
  // y2 = (y*y).real(); ///NOTE I TAKE REAL PART HERE TO MATCH BLUEMLEIN.f FOR TESTING BUT ULTIMATELY SHOULDN'T BE REAL PART (THIS WAS A BUG IN BLUEMLEIN.f) AS IS A POLYGAMMA FUNCTION
  y2 = y*y;
  t0 = (-1.0 + (-1.0 + (-0.5 + (1.0/6.0 + (-1.0/6.0 + (3.0/10.0 + (-5.0/6.0 + 691.0/210.0*y2)*y2)*y2)*y2)*y2)*y)*y)*y2;
  psi2 = t + t0;
  return psi2;
}


//psi3 function as needed in C2 coefficient calculations
std::complex<double> psi3 (std::complex<double> x, int verbosity) {
  std::complex<double> OneComplex;
  OneComplex = std::complex<double>(1.0,0.0);
  std::complex<double> zz,t, y, y2, t0, psi3;
  zz = x;
  t = std::complex<double>(0.0,0.0);
  y = std::complex<double>(0.0,0.0);
  y2 = std::complex<double>(0.0,0.0);
  t0 = std::complex<double>(0.0,0.0);
  double r = 0;
  r = std::abs(zz);
  while (r<=10) {
    t = t + 6.0*OneComplex/(zz*zz*zz*zz);
    zz = zz + OneComplex;
    r = std::abs(zz);
  }
  y = OneComplex/zz;
  // y2 = (y*y).real(); ///NOTE I TAKE REAL PART HERE TO MATCH BLUEMLEIN.f FOR TESTING BUT ULTIMATELY SHOULDN'T BE REAL PART (THIS WAS A BUG IN BLUEMLEIN.f) AS IS A POLYGAMMA FUNCTION
  y2 = y*y;
  t0 = (2.0 + (3.0 + (2.0 + (-1.0 + (4.0/3.0 + (-3.0 + (10.0 + (-691.0/15.0 + (280.0-10851.0/5.0*y2)*y2)*y2)*y2)*y2)*y2)*y2)*y)*y)*y2*y;
  psi3 = t + t0;
  return psi3;
}

//bet function as needed in C2 coefficient calculations
std::complex<double> bet (std::complex<double> x, int verbosity) {
  std::complex<double> OneComplex;
  OneComplex = std::complex<double>(1.0,0.0);
  std::complex<double> zz,z1, z2, v1, v2, bet;
  z1 = std::complex<double>(0.0,0.0);
  z2 = std::complex<double>(0.0,0.0);
  v1 = std::complex<double>(0.0,0.0);
  v2 = std::complex<double>(0.0,0.0);
  bet = std::complex<double>(0.0,0.0);
  std::complex<double> psi0 (std::complex<double> x, int verbosity);
  zz = x;
  z1 = (zz + OneComplex)/2.0;
  z2 = zz/2.0;
  v1 = psi0(z1, verbosity);
  v2 = psi0(z2, verbosity);
  bet = (v1-v2)/2.0;
  return bet;
}

//bet1 function as needed in C2 coefficient calculations
std::complex<double> bet1 (std::complex<double> x, int verbosity) {
  std::complex<double> OneComplex;
  OneComplex = std::complex<double>(1.0,0.0);
  std::complex<double> zz,z1, z2, v1, v2, bet1;
  z1 = std::complex<double>(0.0,0.0);
  z2 = std::complex<double>(0.0,0.0);
  v1 = std::complex<double>(0.0,0.0);
  v2 = std::complex<double>(0.0,0.0);
  bet1 = std::complex<double>(0.0,0.0);
  std::complex<double> psi1 (std::complex<double> x, int verbosity);
  zz = x;
  z1 = (zz + OneComplex)/2.0;
  z2 = zz/2.0;
  v1 = psi1(z1, verbosity);
  v2 = psi1(z2, verbosity);
  bet1 = (v1-v2)/4.0;
  return bet1;
}

//bet2 function as needed in C2 coefficient calculations
std::complex<double> bet2 (std::complex<double> x, int verbosity) {
  std::complex<double> OneComplex;
  OneComplex = std::complex<double>(1.0,0.0);
  std::complex<double> zz,z1, z2, v1, v2, bet2;
  z1 = std::complex<double>(0.0,0.0);
  z2 = std::complex<double>(0.0,0.0);
  v1 = std::complex<double>(0.0,0.0);
  v2 = std::complex<double>(0.0,0.0);
  bet2 = std::complex<double>(0.0,0.0);
  std::complex<double> psi2 (std::complex<double> x, int verbosity);
  zz = x;
  z1 = (zz + OneComplex)/2.0;
  z2 = zz/2.0;
  v1 = psi2(z1, verbosity);
  v2 = psi2(z2, verbosity);
  bet2 = (v1-v2)/8.0;
  return bet2;
}

//bet3 function as needed in C2 coefficient calculations
std::complex<double> bet3 (std::complex<double> x, int verbosity) {
  std::complex<double> OneComplex;
  OneComplex = std::complex<double>(1.0,0.0);
  std::complex<double> zz,z1, z2, v1, v2, bet3;
  z1 = std::complex<double>(0.0,0.0);
  z2 = std::complex<double>(0.0,0.0);
  v1 = std::complex<double>(0.0,0.0);
  v2 = std::complex<double>(0.0,0.0);
  bet3 = std::complex<double>(0.0,0.0);
  std::complex<double> psi3 (std::complex<double> x, int verbosity);
  zz = x;
  z1 = (zz + OneComplex)/2.0;
  z2 = zz/2.0;
  v1 = psi3(z1, verbosity);
  v2 = psi3(z2, verbosity);
  bet3 = (v1-v2)/16.0;
  return bet3;
}

//Mellin transform of log(1+x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg1 (std::complex<double> x, int verbosity) {
 std::complex<double> acg1;
 acg1 = std::complex<double>(0.0,0.0);
 double ak2[10];
 ak2[0] = 0.999999980543793; ak2[1] = -0.999995797779624; ak2[2] = 0.916516447393493;
 ak2[3] = -0.831229921350708; ak2[4] = 0.745873737923571; ak2[5] = -0.634523908078600;
 ak2[6] = 0.467104011423750; ak2[7] = -0.261348046799178; ak2[8] = 0.0936814286867420;
 ak2[9] = -0.0156249375012462;
 double log2 = log(2);
 std::complex<double> t;
 t = std::complex<double>(0.0,0.0);
 int k = 0;
 double m = 0;
 for (int i = 1; i<=10; i++) {
   k = i-1;
   m = i+1;
   t = t + ak2[k]/((x)+m);
 }
 acg1 = (log2*log2 - x*t)/2.0;
 return acg1;
}


//Mellin transform of log(1+x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg1p (std::complex<double> x, int verbosity) {
 std::complex<double> acgp1;
 acgp1 = std::complex<double>(0.0,0.0);
 double ak2[10];
 ak2[0] = 0.999999980543793; ak2[1] = -0.999995797779624; ak2[2] = 0.916516447393493;
 ak2[3] = -0.831229921350708; ak2[4] = 0.745873737923571; ak2[5] = -0.634523908078600;
 ak2[6] = 0.467104011423750; ak2[7] = -0.261348046799178; ak2[8] = 0.0936814286867420;
 ak2[9] = -0.0156249375012462;
 double log2 = log(2);
 std::complex<double> t, t1;
 t = std::complex<double>(0.0,0.0);
 t1 = std::complex<double>(0.0,0.0);
 int k = 0;
 double m = 0;
 for (int i = 1; i<=10; i++) {
   k = i-1;
   m = i+1;
   t = t + ak2[k]/((x)+m);
   t1 = t1 + ak2[k]/(((x)+m)*((x)+m));
 }
 acgp1 = (-t + x*t1)/2.0;
 return acgp1;
}

//Mellin transform of log(1+x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg1pp (std::complex<double> x, int verbosity) {
 std::complex<double> acgpp1;
 acgpp1 = std::complex<double>(0.0,0.0);
 double ak2[10];
 ak2[0] = 0.999999980543793; ak2[1] = -0.999995797779624; ak2[2] = 0.916516447393493;
 ak2[3] = -0.831229921350708; ak2[4] = 0.745873737923571; ak2[5] = -0.634523908078600;
 ak2[6] = 0.467104011423750; ak2[7] = -0.261348046799178; ak2[8] = 0.0936814286867420;
 ak2[9] = -0.0156249375012462;
 double log2 = log(2);
 std::complex<double> t, t1;
 t = std::complex<double>(0.0,0.0);
 t1 = std::complex<double>(0.0,0.0);
 int k = 0;
 double m = 0;
 for (int i = 1; i<=10; i++) {
   k = i-1;
   m = i+1;
   t = t + ak2[k]/(((x)+m)*((x)+m));
   t1 = t1 + ak2[k]/(((x)+m)*((x)+m)*((x)+m));
 }
 acgpp1 = (t - x*t1);
 return acgpp1;
}

//Mellin transform of log(1+x)^2/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg2 (std::complex<double> x, int verbosity) {
 std::complex<double> acg2;
 acg2 = std::complex<double>(0.0,0.0);
 double ak3[11];
 ak3[0] = 0.999999989322696; ak3[1] = -1.49999722020708; ak3[2] = 1.74988008499745; ak3[3] = -1.87296689068405;
 ak3[4] = 1.91539974617231; ak3[5] = -1.85963744001295; ak3[6] = 1.62987195424434; ak3[7] = -1.17982353224299;
 ak3[8] = 0.628710122994999; ak3[9] = -0.211307487211713; ak3[10] = 0.0328953352932140;
 double log2 = log(2);
 std::complex<double> t;
 t = std::complex<double>(0.0,0.0);
 int k = 0;
 double m = 0;
 for (int i = 2; i<=12; i++) {
   k = i-2;
   m = i+1;
   t = t + ak3[k]/((x)+m);
 }
 acg2 = (log2*log2*log2 - x*t)/3.0;
 return acg2;
}

//Mellin transform of log(1+x)^2/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg2p (std::complex<double> x, int verbosity) {
 std::complex<double> acg2p;
 acg2p = std::complex<double>(0.0,0.0);
 double ak3[11];
 ak3[0] = 0.999999989322696; ak3[1] = -1.49999722020708; ak3[2] = 1.74988008499745; ak3[3] = -1.87296689068405;
 ak3[4] = 1.91539974617231; ak3[5] = -1.85963744001295; ak3[6] = 1.62987195424434; ak3[7] = -1.17982353224299;
 ak3[8] = 0.628710122994999; ak3[9] = -0.211307487211713; ak3[10] = 0.0328953352932140;
 double log2 = log(2);
 std::complex<double> t, t1;
 t = std::complex<double>(0.0,0.0);
 t1 = std::complex<double>(0.0,0.0);
 int k = 0;
 double m = 0;
 for (int i = 2; i<=12; i++) {
   k = i-2;
   m = i+1;
   t = t + ak3[k]/((x)+m);
   t1 = t1 + ak3[k]/(((x)+m)*((x)+m));
 }
 acg2p = (-t + x*t1)/3.0;
 return acg2p;
}

//Mellin transform of Li2(x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg3 (std::complex<double> x, int verbosity) {
 std::complex<double> acg3;
 acg3 = std::complex<double>(0.0,0.0);
 double ak1[9];
 ak1[0] = 0.999999974532240; ak1[1] = -0.499995525890027; ak1[2] = 0.333203435554182; ak1[3] = -0.248529457735332;
 ak1[4] = 0.191451164493502; ak1[5] = -0.137466222203386; ak1[6] = 0.0792107405737825; ak1[7] = -0.0301109652783781;
 ak1[8] = 0.00538406198111749;
 double log2 = log(2);
 std::complex<double> t, zz, zzz;
 t = std::complex<double>(0.0,0.0);
 zz = x + std::complex<double>(1.0,0.0);
 zzz = zz + std::complex<double>(1.0,0.0);
 t = 1.004/zzz - 0.846/(zzz+std::complex<double>(1.0,0.0)) + 1.342/(x + 4.0*std::complex<double>(1.0,0.0)) - 1.532/(x + 5.0*std::complex<double>(1.0,0.0)) + 0.839/(x + 6.0*std::complex<double>(1.0,0.0));
 acg3 = t;
 return acg3;
}

//Mellin transform of Li2(-x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg4 (std::complex<double> x, int verbosity) {
 std::complex<double> acg4;
 acg4 = std::complex<double>(0.0,0.0);

 double zeta2 = k_constants::zeta2;

 double ak1[9];
 ak1[0] = 0.999999974532240; ak1[1] = -0.499995525890027; ak1[2] = 0.333203435554182; ak1[3] = -0.248529457735332;
 ak1[4] = 0.191451164493502; ak1[5] = -0.137466222203386; ak1[6] = 0.0792107405737825; ak1[7] = -0.0301109652783781;
 ak1[8] = 0.00538406198111749;
 double log2 = log(2);
 std::complex<double> t, zz, zzz, zzz1, v1, log2comp;
 std::complex<double> bet (std::complex<double> x, int verbosity);
 t = std::complex<double>(-zeta2/2.0*log2,0.0);
 zz = std::complex<double>(0.0,0.0);
 zzz = std::complex<double>(0.0,0.0);
 zzz1 = std::complex<double>(0.0,0.0);
 v1 = std::complex<double>(0.0,0.0);
 log2comp = std::complex<double>(log2,0.0);
 double m = 0;
 for (int i =1; i<=9; i++) {
 m = i;
 zz = std::complex<double>(m,0.0);
 zzz = x + zz;
 zzz1 = zzz + std::complex<double>(1.0,0.0);
 v1 = bet(zzz1, verbosity);
 t = t + ak1[i-1]*(x/zzz*zeta2/2.0 + zz/(zzz*zzz)*(log2comp-v1));
 }
 acg4 = t;
 return acg4;
}

//Mellin transform of Li2(-x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg4p (std::complex<double> x, int verbosity) {
 std::complex<double> acg4p;
 acg4p = std::complex<double>(0.0,0.0);

 double zeta2 = k_constants::zeta2;

 double ak1[9];
 ak1[0] = 0.999999974532240; ak1[1] = -0.499995525890027; ak1[2] = 0.333203435554182; ak1[3] = -0.248529457735332;
 ak1[4] = 0.191451164493502; ak1[5] = -0.137466222203386; ak1[6] = 0.0792107405737825; ak1[7] = -0.0301109652783781;
 ak1[8] = 0.00538406198111749;
 double log2 = log(2);
 std::complex<double> t, zz, zzz, zzz1, v1, v2, log2comp;
 std::complex<double> bet (std::complex<double> x, int verbosity);
 std::complex<double> bet1 (std::complex<double> x, int verbosity);
 t = std::complex<double>(0.0,0.0);
 zz = std::complex<double>(0.0,0.0);
 zzz = std::complex<double>(0.0,0.0);
 zzz1 = std::complex<double>(0.0,0.0);
 v1 = std::complex<double>(0.0,0.0);
 v2 = std::complex<double>(0.0,0.0);
 log2comp = std::complex<double>(log2,0.0);
 double m = 0;
 for (int i =1; i<=9; i++) {
 m = i;
 zz = std::complex<double>(m,0.0);
 zzz = x + zz;
 zzz1 = zzz + std::complex<double>(1.0,0.0);
 v1 = bet(zzz1, verbosity);
 v2 = bet1(zzz1, verbosity);
 t = t + ak1[i-1]*((zzz-x)/(zzz*zzz)*zeta2/2.0 - 2.0*zz/(zzz*zzz*zzz)*(log2comp-v1)-zz/(zzz*zzz)*v2);
 }
 acg4p = t;
 return acg4p;
}

//Mellin transform of log(x)*Li2(x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg5 (std::complex<double> x, int verbosity) {
 std::complex<double> acg5;
 acg5 = std::complex<double>(0.0,0.0);

 double zeta2 = k_constants::zeta2;
 double Eulerconst = k_constants::Eulerconst;

 double ak1[9];
 ak1[0] = 0.999999974532240; ak1[1] = -0.499995525890027; ak1[2] = 0.333203435554182; ak1[3] = -0.248529457735332;
 ak1[4] = 0.191451164493502; ak1[5] = -0.137466222203386; ak1[6] = 0.0792107405737825; ak1[7] = -0.0301109652783781;
 ak1[8] = 0.00538406198111749;
 double log2 = log(2);
 std::complex<double> t, zz, zzz, zzz1, psi0a, psi1a, log2comp;
 std::complex<double> psi0 (std::complex<double> x, int verbosity);
 std::complex<double> psi1 (std::complex<double> x, int verbosity);
 t = std::complex<double>(0.0,0.0);
 zz = std::complex<double>(0.0,0.0);
 zzz = std::complex<double>(0.0,0.0);
 zzz1 = std::complex<double>(0.0,0.0);
 psi0a = std::complex<double>(0.0,0.0);
 psi1a = std::complex<double>(0.0,0.0);
 log2comp = std::complex<double>(log2,0.0);
 double m = 0;
 for (int i =1; i<=9; i++) {
 m = i;
 zz = std::complex<double>(m,0.0);
 zzz = x + zz;
 zzz1 = zzz + std::complex<double>(1.0,0.0);
 psi0a = psi0(zzz1, verbosity);
 psi1a = psi1(zzz1, verbosity);
 t = t - ak1[i-1]*zz/(zzz*zzz)*(zeta2+psi1a-2.0*(psi0a+Eulerconst)/zzz);
 }
 acg5 = t;
 return acg5;
}

//Mellin transform of Li3(x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg6 (std::complex<double> x, int verbosity) {
 std::complex<double> acg6;
 acg6 = std::complex<double>(0.0,0.0);

 double zeta2 = k_constants::zeta2;
 double zeta3 = k_constants::zeta3;
 double Eulerconst = k_constants::Eulerconst;

 double ak1[9];
 ak1[0] = 0.999999974532240; ak1[1] = -0.499995525890027; ak1[2] = 0.333203435554182; ak1[3] = -0.248529457735332;
 ak1[4] = 0.191451164493502; ak1[5] = -0.137466222203386; ak1[6] = 0.0792107405737825; ak1[7] = -0.0301109652783781;
 ak1[8] = 0.00538406198111749;
 double log2 = log(2);
 std::complex<double> t, zz, zzz, zzz1, v1, log2comp;
 std::complex<double> psi0 (std::complex<double> x, int verbosity);
 std::complex<double> psi1 (std::complex<double> x, int verbosity);
 t = std::complex<double>(log2*zeta3,0.0);
 zz = std::complex<double>(0.0,0.0);
 zzz = std::complex<double>(0.0,0.0);
 zzz1 = std::complex<double>(0.0,0.0);
 v1 = std::complex<double>(0.0,0.0);
 log2comp = std::complex<double>(log2,0.0);
 double m = 0;
 for (int i =1; i<=9; i++) {
 m = i;
 zz = std::complex<double>(m,0.0);
 zzz = x + zz;
 zzz1 = zzz + std::complex<double>(1.0,0.0);
 v1 = psi0(zzz1, verbosity);
 t = t - ak1[i-1]*(x/(zzz)*zeta3 + zz/(zzz*zzz)*(zeta2 - (v1+Eulerconst)/zzz));
 }
 acg6 = t;
 return acg6;
}

//Mellin transform of Li3(-x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg7 (std::complex<double> x, int verbosity) {
 std::complex<double> acg7;
 acg7 = std::complex<double>(0.0,0.0);

 double zeta2 = k_constants::zeta2;
 double zeta3 = k_constants::zeta3;

 double ak1[9];
 ak1[0] = 0.999999974532240; ak1[1] = -0.499995525890027; ak1[2] = 0.333203435554182; ak1[3] = -0.248529457735332;
 ak1[4] = 0.191451164493502; ak1[5] = -0.137466222203386; ak1[6] = 0.0792107405737825; ak1[7] = -0.0301109652783781;
 ak1[8] = 0.00538406198111749;
 double log2 = log(2);
 std::complex<double> t, zz, zzz, zzz1, v1, log2comp;
 std::complex<double> bet (std::complex<double> x, int verbosity);
 t = std::complex<double>(-3.0*zeta3/4.0*log2,0.0);
 zz = std::complex<double>(0.0,0.0);
 zzz = std::complex<double>(0.0,0.0);
 zzz1 = std::complex<double>(0.0,0.0);
 v1 = std::complex<double>(0.0,0.0);
 log2comp = std::complex<double>(log2,0.0);
 double m = 0;
 for (int i =1; i<=9; i++) {
     m = i;
     zz = std::complex<double>(m,0.0);
     zzz = x + zz;
     zzz1 = zzz + std::complex<double>(1.0,0.0);
     v1 = bet(zzz1, verbosity);
     t = t + ak1[i-1]*(x/(zzz)*3.0*zeta3/4.0 + zz/(2.0*zzz*zzz)*zeta2 - zz/(zzz*zzz*zzz)*(log2comp-v1));
 }
 acg7 = t;
 return acg7;
}

//Mellin transform of S12(x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg8 (std::complex<double> x, int verbosity) {
 std::complex<double> acg8;
 acg8 = std::complex<double>(0.0,0.0);

 double zeta2 = k_constants::zeta2;
 double zeta3 = k_constants::zeta3;
 double Eulerconst = k_constants::Eulerconst;

 double ak1[9];
 ak1[0] = 0.999999974532240; ak1[1] = -0.499995525890027; ak1[2] = 0.333203435554182; ak1[3] = -0.248529457735332;
 ak1[4] = 0.191451164493502; ak1[5] = -0.137466222203386; ak1[6] = 0.0792107405737825; ak1[7] = -0.0301109652783781;
 ak1[8] = 0.00538406198111749;
 double log2 = log(2);
 std::complex<double> t, zz, zzz, zzz1, psi0a, psi1a, log2comp;
 std::complex<double> psi0 (std::complex<double> x, int verbosity);
 std::complex<double> psi1 (std::complex<double> x, int verbosity);
 t = std::complex<double>(0.0,0.0);
 zz = std::complex<double>(0.0,0.0);
 zzz = std::complex<double>(0.0,0.0);
 zzz1 = std::complex<double>(0.0,0.0);
 psi0a = std::complex<double>(0.0,0.0);
 psi1a = std::complex<double>(0.0,0.0);
 log2comp = std::complex<double>(log2,0.0);
 double m = 0;
 for (int i =1; i<=9; i++) {
 m = i;
 zz = std::complex<double>(m,0.0);
 zzz = x + zz;
 zzz1 = zzz + std::complex<double>(1.0,0.0);
 psi0a = psi0(zzz1, verbosity);
 psi1a = psi1(zzz1, verbosity);
 t = t - ak1[i-1]*(x*zeta3/zzz + zz/(2.0*zzz*zzz)*((psi0a+Eulerconst)*(psi0a+Eulerconst)+(zeta2 - psi1a)));
 }
 acg8 = t + log2comp*zeta3;
 return acg8;
}

//Mellin transform of S12(-x)/(1+x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg9 (std::complex<double> x, int verbosity) {
 std::complex<double> acg9;
 acg9 = std::complex<double>(0.0,0.0);

 double zeta3 = k_constants::zeta3;

 double ak1[9];
 ak1[0] = 0.999999974532240; ak1[1] = -0.499995525890027; ak1[2] = 0.333203435554182; ak1[3] = -0.248529457735332;
 ak1[4] = 0.191451164493502; ak1[5] = -0.137466222203386; ak1[6] = 0.0792107405737825; ak1[7] = -0.0301109652783781;
 ak1[8] = 0.00538406198111749;
 double ak2[10];
 ak2[0] = 0.999999980543793; ak2[1] = -0.999995797779624; ak2[2] = 0.916516447393493;
 ak2[3] = -0.831229921350708; ak2[4] = 0.745873737923571; ak2[5] = -0.634523908078600;
 ak2[6] = 0.467104011423750; ak2[7] = -0.261348046799178; ak2[8] = 0.0936814286867420;
 ak2[9] = -0.0156249375012462;
 double ak3[11];
 ak3[0] = 0.999999989322696; ak3[1] = -1.49999722020708; ak3[2] = 1.74988008499745; ak3[3] = -1.87296689068405;
 ak3[4] = 1.91539974617231; ak3[5] = -1.85963744001295; ak3[6] = 1.62987195424434; ak3[7] = -1.17982353224299;
 ak3[8] = 0.628710122994999; ak3[9] = -0.211307487211713; ak3[10] = 0.0328953352932140;
 double log2 = log(2);
 std::complex<double> t, t1, t2, znkl, zk, zz;
 t = std::complex<double>(zeta3*log2/8.0,0.0);
 t2 = std::complex<double>(0.0,0.0);
 zz = std::complex<double>(0.0,0.0);
 int l1 = 0, m = 0, n = 0, a = 0;
 for (int k = 0; k<=8; k++) {
     t1 = std::complex<double>(0.0,0.0);
     for (int l = 1; l<=10; l++) {
	 l1 = l-1;
	 m = k+l+2.0;
	 znkl = x + std::complex<double>(m,0);
	 t1 = t1 + ak2[l1]/znkl;
     }
     n = k + 1.0;
     zk = std::complex<double>(n,0.0);
     zz = x + zk;
     t = t - ak1[k]*x/zz*(zeta3/8.0-t1/2.0);
 }
 for (int p = 2; p<13; p++) {
     a = p -2;
     n = p + 1.0;
     zk = std::complex<double>(n,0.0);
     zz = x + zk;
     t2 = t2 + ak3[a]/zz;
 }
 t = t - t2/2.0;
 acg9 = t;
 return acg9;
}

//Mellin transform of log(1+x)/(1+x)*LI2(-x) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg13 (std::complex<double> x, int verbosity) {
 std::complex<double> acg13;
 acg13 = std::complex<double>(0.0,0.0);

 double zeta2 = k_constants::zeta2;

 double ak2[10];
 ak2[0] = 0.999999980543793; ak2[1] = -0.999995797779624; ak2[2] = 0.916516447393493;
 ak2[3] = -0.831229921350708; ak2[4] = 0.745873737923571; ak2[5] = -0.634523908078600;
 ak2[6] = 0.467104011423750; ak2[7] = -0.261348046799178; ak2[8] = 0.0936814286867420;
 ak2[9] = -0.0156249375012462;
 double ak3[11];
 ak3[0] = 0.999999989322696; ak3[1] = -1.49999722020708; ak3[2] = 1.74988008499745; ak3[3] = -1.87296689068405;
 ak3[4] = 1.91539974617231; ak3[5] = -1.85963744001295; ak3[6] = 1.62987195424434; ak3[7] = -1.17982353224299;
 ak3[8] = 0.628710122994999; ak3[9] = -0.211307487211713; ak3[10] = 0.0328953352932140;
 double log2 = log(2);
 std::complex<double> t0, t1, t2, znk1, v1, log2comp;
 std::complex<double> bet (std::complex<double> x, int verbosity);
 t0 = std::complex<double>(-1.0/4.0*zeta2*log2*log2,0.0);
 t1 = std::complex<double>(0.0,0.0);
 t2 = std::complex<double>(0.0,0.0);
 znk1 = std::complex<double>(0.0,0.0);
 v1 = std::complex<double>(0.0,0.0);
 log2comp = std::complex<double>(log2,0.0);
 double n = 0; int a = 0;
 for (int p = 2; p<13; p++) {
     a = p -2;
     n = p + 1.0;
     t1 = t1 + ak3[a]/(x+n);
 }
 for (int l = 1; l < 11; l++) {
     a = l-1;
     n = l + 2.0;
     znk1 = x + n;
     v1 = bet(znk1, verbosity);
     t2 = t2 + ak2[a]*x/(znk1 - std::complex<double>(1.0,0.0))*(zeta2/2.0 - (log2comp-v1)/(znk1 - std::complex<double>(1.0,0.0)));
 }
 acg13 = t0 + (t1+t2)/2.0;
 return acg13;
}

//Mellin transform of (Li3(x)-zeta3)/(x-1) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg20 (std::complex<double> x, int verbosity) {
 std::complex<double> acg20;

 double zeta2 = k_constants::zeta2;
 double zeta3 = k_constants::zeta3;
 double Eulerconst = k_constants::Eulerconst;

 double ck2[13];
 ck2[0] = 0.0, ck2[1] = 0.0; ck2[2] = 2.1253046159349207; ck2[3] = -1.0523034804772378;
 ck2[4] = 0.160179976133047; ck2[5] = -0.351982236677917; ck2[6] = 1.41033316846244; ck2[7] = -3.53343997843391;
 ck2[8] = 5.93934696819832; ck2[9] = -6.60019784804042; ck2[10] = 4.66330349413063; ck2[11] = -1.89825467489058;
 ck2[12] = 0.339772909487512;
 double ck4[13];
 ck4[0] =  0.0; ck4[1] = 2.215008697869307; ck4[2] = -0.9133677154535804; ck4[3] = 3.4783104357500143;
 ck4[4] = -2.823955592989266; ck4[5] = 0.992890266001707; ck4[6] = -1.30026190226546; ck4[7] = 3.41870577921103;
 ck4[8] = -5.76763902370864; ck4[9] = 6.45554138192407; ck4[10] = -4.59405622046138; ck4[11] = 1.88510809558304;
 ck4[12] = -0.340476080290674;
 double p24[5], p34[5], p22[4], p12[4], p14[5];
 p12[0] = zeta3 - 11.0/6.0*zeta2 + 4.0/3.0;
 p12[1] = 3.0*zeta2 - 13.0/4.0;
 p12[2] = -3.0/2.0*zeta2 + 5.0/2.0;
 p12[3] = 1.0/3.0*zeta2 - 7.0/12.0;
 p14[0] = 257.0/144.0 - 205.0/72.0*zeta2 + zeta2*zeta2;
 p14[1] = -167.0/36.0 + 25.0/6.0*zeta2;
 p14[2] = 101.0/24.0 - 23.0/12.0*zeta2;
 p14[3] = -59.0/36.0 + 13.0/18.0*zeta2;
 p14[4] = 41.0/144.0 - zeta2/8.0;
 p22[0] = -1.0;
 p22[1] = 5.0/2.0;
 p22[2] = -2.0;
 p22[3] = 0.5;
 p24[0] = -167.0/36.0 + 25.0/6.0*zeta2;
 p24[1] = 235.0/18.0 - 8.0*zeta2;
 p24[2] = -40.0/3.0 + 6.0*zeta2;
 p24[3] = 109.0/18.0 - 8.0/3.0*zeta2;
 p24[4] = -41.0/36.0 + zeta2/2.0;
 p34[0] = 35.0/12.0;
 p34[1] = -26.0/3.0;
 p34[2] = 19.0/2.0;
 p34[3] = -14.0/3.0;
 p34[4] = 11.0/12.0;
 std::complex<double> t, zn1, ps, znk, znk1, ps1, r;
 std::complex<double> psi0(std::complex<double> x, int verbosity);
 std::complex<double> psi1(std::complex<double> x, int verbosity);
 t = std::complex<double>(zeta2*zeta2/2.0,0.0);
 zn1 = std::complex<double>(0.0,0.0);
 ps = std::complex<double>(0.0,0.0);
 znk = std::complex<double>(0.0,0.0);
 znk1 = std::complex<double>(0.0,0.0);
 ps1 = std::complex<double>(0.0,0.0);
 r = std::complex<double>(0.0,0.0);
 zn1 = x + std::complex<double>(1.0,0.0);
 ps = psi0(zn1, verbosity);
 t = t - zeta3*(ps + Eulerconst);
 double l = 0;
 for (int k = 0; k<13; k++) {
     l = k;
     znk = x + std::complex<double>(l,0.0);
     znk1 = znk + std::complex<double>(1.0,0.0);
     ps = psi0(znk1, verbosity);
     r = x/(znk);
     t = t + ck2[k]*r*(ps+Eulerconst);
 }
 std::complex<double> t1;
 t1 = std::complex<double>(0.0,0.0);
 for (int k = 0; k<4; k++) {
     l = k;
     znk = x + std::complex<double>(l,0.0);
     znk1 = znk + std::complex<double>(1.0,0.0);
     ps = psi0(znk1, verbosity);
     ps1 = psi1(znk1, verbosity);
     r = x/znk;
     t1 = (ps+Eulerconst)*(ps+Eulerconst) - ps1 + zeta2;
     t = t - p22[k]*r*t1;
 }
 for (int k = 0; k<13; k++) {
     l = k;
     znk = x + std::complex<double>(l,0.0);
     znk1 = znk + std::complex<double>(1.0,0.0);
     ps = psi0(znk1, verbosity);
     t = t - ck4[k]/znk*x/2.0;
 }
 std::complex<double> t2;
 t1 = std::complex<double>(0.0,0.0);
 t2 = std::complex<double>(0.0,0.0);
 for (int k = 0; k < 5; k++) {
     l = k;
     znk = x + std::complex<double>(l,0.0);
     znk1 = znk + std::complex<double>(1.0,0.0);
     ps = psi0(znk1, verbosity);
     ps1 = psi1(znk1, verbosity);
     t1 = -(ps+Eulerconst)/znk;
     t2 = ((ps+Eulerconst)*(ps+Eulerconst) + (-ps1+zeta2))/znk;
     t = t - (p24[k]*t1+p34[k]*t2)*x/2.0;
 }
 acg20 = t;
 return acg20;
}

//Mellin transform of (s12(x)-zeta3)/(x-1) for complex argument as needed in C2 coefficient calculations
std::complex<double> acg21 (std::complex<double> x, int verbosity) {
 std::complex<double> acg21;
 acg21 =  std::complex<double> (0.0,0.0);

 double zeta2 = k_constants::zeta2;
 double zeta3 = k_constants::zeta3;
 double Eulerconst = k_constants::Eulerconst;

 double ck3[10];
 ck3[0] = 0.0; ck3[1] = 1.4236162474052558; ck3[2] = -0.08001203559240111; ck3[3] = -0.39875367195395994;
 ck3[4] = 0.339241791547134, ck3[5] = -0.0522116678353452; ck3[6] = -0.0648354706049337; ck3[7] = 0.0644165053822532;
 ck3[8] = -0.0394927322542075; ck3[9] = 0.0100879370657869;
 double p13[5], p23[5], p33[5];
 p13[0] = zeta3 - 2035.0/1728.0;
 p13[1] = 205.0/144.0;
 p13[2] = -95.0/288.0;
 p13[3] = 43.0/432.0;
 p13[4] = -1.0/64.0;
 p23[0] = 205.0/144.0;
 p23[1] = -25.0/12.0;
 p23[2] = 23.0/24.0;
 p23[3] = -13.0/36.0;
 p23[4] = 1.0/16.0;
 p33[0] = -25.0/24.0;
 p33[1] = 2.0;
 p33[2] = -1.5;
 p33[3] = 2.0/3.0;
 p33[4] = -1.0/8.0;
 std::complex<double> t, zn1, ps, ps1, ps2;
 t =  std::complex<double> (0.0,0.0);
 zn1 =  std::complex<double> (0.0,0.0);
 ps =  std::complex<double> (0.0,0.0);
 ps1 =  std::complex<double> (0.0,0.0);
 ps2 =  std::complex<double> (0.0,0.0);
 std::complex<double> psi0 (std::complex<double> x, int verbosity);
 std::complex<double> psi1 (std::complex<double> x, int verbosity);
 std::complex<double> psi2 (std::complex<double> x, int verbosity);
 zn1 = x  + std::complex<double>(1.0,0.0);
 ps = psi0(zn1, verbosity);
 ps1 = psi1(zn1, verbosity);
 ps2 = psi2(zn1, verbosity);
 t = t - zeta3*(ps+Eulerconst);
 t = t + ((ps+Eulerconst)*(ps+Eulerconst)*(ps+Eulerconst) + 3.0*(ps+Eulerconst)*(zeta2-ps1) + 2.0*(ps2/2.0 + zeta3))/(2.0*x);
 double l = 0;
 std::complex<double> znk,znk1;
 znk =  std::complex<double> (0.0,0.0);
 znk1 =  std::complex<double> (0.0,0.0);
 for (int k=0;k<10;k++) {
     l = k;
     znk = x + std::complex<double>(l,0.0);
     znk1 = znk + + std::complex<double>(1.0,0.0);
     ps = psi0(znk1, verbosity);
     t = t + ck3[k]*(ps+Eulerconst)*x/znk;
 }
 ps =  std::complex<double> (0.0,0.0);
 znk =  std::complex<double> (0.0,0.0);
 znk1 =  std::complex<double> (0.0,0.0);
 std::complex<double> t2, t3;
 t2 =  std::complex<double> (0.0,0.0);
 t3 =  std::complex<double> (0.0,0.0);
 for (int k = 0; k<5; k++) {
     l = k;
     znk = x + std::complex<double>(l,0.0);
     znk1 = znk + std::complex<double>(1.0,0.0);
     ps = psi0(znk1, verbosity);
     ps1 = psi1(znk1, verbosity);
     ps2 = psi2(znk1, verbosity);
     t2 = (ps+Eulerconst)*(ps+Eulerconst)+zeta2-ps1;
     t3 = (ps+Eulerconst)*(ps+Eulerconst)*(ps+Eulerconst)+3.0*(ps+Eulerconst)*(zeta2-ps1) + 2.0*(ps2/2.0 + zeta3);
     t = t + x/znk*(p33[k]*t3-p23[k]*t2);
 }
 acg21 = t;
 return acg21;
}


//Function used in mellin moments of pdfs
std::complex<double> cbeta(std::complex<double> z1, std::complex<double> z2, int verbosity) {
  std::complex<double> sub, zz1, zz2, lg1, lg2, lg12, cbeta;
  sub = std::complex<double>(0.0,0.0);
  zz1 = std::complex<double>(0.0,0.0);
  zz2 = std::complex<double>(0.0,0.0);
  lg1 = std::complex<double>(0.0,0.0);
  lg2 = std::complex<double>(0.0,0.0);
  lg12 = std::complex<double>(0.0,0.0);
  cbeta = std::complex<double>(0.0,0.0);
  std::complex<double> lngam(std::complex<double> z, int verbosity);
  zz1 = z1;
  while (std::real(zz1) < 15.0) {
      sub = sub + std::log((zz1+z2)/zz1);
      zz1 = zz1 + std::complex<double>(1.0,0.0);
  }
  zz2 = z2;
  while (std::real(zz2) < 15.0) {
      sub = sub + std::log((zz1+zz2)/zz2);
      zz2 = zz2 + std::complex<double>(1.0,0.0);
  }
  lg1 = lngam(zz1, verbosity);
  lg2 = lngam(zz2, verbosity);
  lg12 = lngam(zz1+zz2, verbosity);
  cbeta = std::exp(lg1+lg2-lg12+sub);
  return cbeta;
}

std::complex<double> lngam(std::complex<double> z, int verbosity) {
  std::complex<double> l;
  l = std::complex<double>(0.0,0.0);
  l = (z-std::complex<double>(0.5,0.0))*std::log(z) - z + std::complex<double>(0.91893853,0.0) + 1.0/(12.0*z)*(1.0 - 1.0/(30.0*z*z)*(1.0 - 1.0/(3.5*z*z)*(1.0 - 4.0/(3.0*z*z))));
  return l;
}
