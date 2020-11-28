
#include "DYjet_hard.h"


#include "constants.h"
#include "drellyan_input.h"
#include "diphoton_hard.h" //for Asmitad etc



void drellyan_amplitude::sigmaij_dyjet_calc(drellyan_input* drellyan_in, double as, std::vector<std::vector<double> >& sigmaij){
// Note that despite the function being named "sigmaij", this is not really a proper partonic cross-section,
// but rather the squared amplitude |M|^2 multiplied by units, fluz factor and spin and colour averages,
// that is, it does not include any PS measure (which has to be accounted for elsewhere)
//
// Overall factor (not including coupling constants)
    double q2 = 2.*ss(1,2);
    // std::cout << "q2 = 2*ss(1,2) = " << q2 << std::endl;
    double gevpb = drellyan_in->res_1.gevm2topb;
    double pi = k_constants::pi;
// Overall factor for amplitudes
    double flux = 1./2./q2;
    double spinavg = 1/4.;
    int Nc = drellyan_in->res_1.Nc;
    double colfac = (Nc*Nc-1.)/2.;
//
    double fac_overall = flux*spinavg*colfac*gevpb; // Note that colour average is channel-dependent and not included here; coupling constants are also added later.
//
// Input parameters
/* DYprocess: selection of mediator boson:
1 = W+, 2 = W-, 3 = both W+ and W- (final state charge not measured)
4 = Z only, 5 = gamma* + Z, 6 = gamma* + Z + Z' */
    int DYproc = drellyan_in->DYprocess;
    int DYnwapprox= drellyan_in->DYnarrowwidthapprox; // if 1 --> use narrow width approx on W and Z resonance
    int Nf = drellyan_in->res_1.Nf;
// Photon
    double alphaem = drellyan_in->res_1.alpha_QED;
    double ee = std::sqrt(4*pi*alphaem);
// W
    double mw = drellyan_in->mw;
    double ww = drellyan_in->ww;
    double gw = drellyan_in->gw;
// Z
    double mz = drellyan_in->mz;
    double zw = drellyan_in->zw;
    double gz = drellyan_in->gz;
    double sw2 = drellyan_in->sw2;
// Z'
    double mzp = drellyan_in->mzp; // Z prime mass
    double gzp = drellyan_in->gzp; // Z prime coupling
    double zpw = drellyan_in->zpw; // Z prime width
// Couplings to fermions:
// Photon
    double Qu = 2./3.;
    double Qd = -1./3.;
// Array-style (useful for Z' generalisation)
    double QQ[5] = {Qu, Qd, Qd, Qu, Qd};
// W
    double cLw = gw/(std::sqrt(2.));
    double Vud = drellyan_in->Vud;
    double Vus = drellyan_in->Vus;
    double Vub = drellyan_in->Vub;
    double Vcd = drellyan_in->Vcd;
    double Vcs = drellyan_in->Vcs;
    double Vcb = drellyan_in->Vcb;
// Z
    double cLzu=gz*(0.5-Qu*sw2);
    double cLzd=gz*(-0.5-Qd*sw2);
    double cRzu=-gz*Qu*sw2;
    double cRzd=-gz*Qd*sw2;
    double cLzl=gz*(-0.5+sw2);
    double cRzl=gz*sw2;
// Array-style (useful for Z' generalisation)
    double gLzQ[5] = {cLzu, cLzd, cLzd, cLzu, cLzd};
    double gRzQ[5] = {cRzu, cRzd, cRzd, cRzu, cRzd};
//
// Z' -- assuming no FCNC for now
    double cLu1 = gzp*drellyan_in->cLu1;
    double cLu2 = gzp*drellyan_in->cLu2;
    double cLd1 = gzp*drellyan_in->cLd1;
    double cLd2 = gzp*drellyan_in->cLd2;
    double cLd3 = gzp*drellyan_in->cLd3;
    double cLZpl = gzp*drellyan_in->cLl1;
    double cRu1 = gzp*drellyan_in->cRu1;
    double cRu2 = gzp*drellyan_in->cRu2;
    double cRd1 = gzp*drellyan_in->cRd1;
    double cRd2 = gzp*drellyan_in->cRd2;
    double cRd3 = gzp*drellyan_in->cRd3;
    double cRZpl = gzp*drellyan_in->cRl1;
    double gLZpQ[5] = {cLu1, cLd1, cLd2, cLu2, cLd3};
    double gRZpQ[5] = {cRu1, cRd1, cRd2, cRu2, cRd3};
//
// Propagator factors for photon, Z and Z'
    std::complex<double> II(0., 1.);
    double PropDen_y = 1./q2;
    std::complex<double> PropDen_Z = 1./(q2 - mz*mz + II*mz*zw);
    std::complex<double> PropDen_Zp = 1./(q2 - mzp*mzp + II*mzp*zpw);

////////////////////////

// Charged current: subprocesses 1-3 -- NOT IMPLEMENTED YET
//
    if(DYproc == 1 || DYproc == 2 || DYproc == 3) {
      std::cout << "W+- channel yet TO DO" << std::endl;
      exit(EXIT_FAILURE);
    }
//
// Neutral current: subprocesses 4-6
//
// Z boson only (no photon): subprocess 4

    else if(DYproc == 4 || DYproc == 5 || DYproc == 6){

      std::vector<std::complex<double> > ampkinfac_qqb(8), ampkinfac_qbq(8), ampkinfac_qg(8), ampkinfac_gq(8), ampkinfac_gqb(8), ampkinfac_qbg(8);
// See member function definition for an explanation of the notation
      this->ampkinfacs(1, 2, -3, -4, -5, ampkinfac_qqb);
      this->ampkinfacs(2, 1, -3, -4, -5, ampkinfac_qbq);
      this->ampkinfacs(1, 5, -3, -4, -2, ampkinfac_qg);
      this->ampkinfacs(5, 1, -3, -4, -2, ampkinfac_gq);
      this->ampkinfacs(5, 2, -3, -4, -1, ampkinfac_gqb);
      this->ampkinfacs(2, 5, -3, -4, -1, ampkinfac_qbg);

      std::vector<std::complex<double> > amp_dyn(Nf);
      std::complex<double> pre_amp;
      double pre_ampsq, col_avg;

      if (DYproc == 4) {
        for(unsigned int ii=0; ii<Nf; ii++){
          amp_dyn[ii] = gLzQ[ii]*cLzl*PropDen_Z;
        }
      }
// gamma* + Z -- subprocess 5
      else if (DYproc == 5) {
        for(unsigned int ii=0; ii<Nf; ii++){
          amp_dyn[ii] = -QQ[ii]*ee*ee*PropDen_y + gLzQ[ii]*cLzl*PropDen_Z;
        }
      }
//
// gamma* + Z + Z' -- subprocess 6
      else if (DYproc == 6) {
        for(unsigned int ii=0; ii<Nf; ii++){
          amp_dyn[ii] = -QQ[ii]*ee*ee*PropDen_y + gLzQ[ii]*cLzl*PropDen_Z + gLZpQ[ii]*cLZpl*PropDen_Zp;
        }
      }

// Put pieces together
      for(unsigned int ii=0; ii<Nf; ii++){
// qqb channel
        pre_ampsq = 0.;
        col_avg = 1./(Nc*Nc);
        for(unsigned int jj=0; jj<8; jj++){
          pre_amp = ampkinfac_qqb[jj]*amp_dyn[ii];
          pre_ampsq += 4.*pi*as*std::norm(pre_amp);
        }
        sigmaij[Nf+ii][Nf-ii] = fac_overall*pre_ampsq;
// qbq channel
        pre_ampsq = 0.;
        col_avg = 1./(Nc*Nc);
        for(unsigned int jj=0; jj<8; jj++){
          pre_amp = ampkinfac_qbq[jj]*amp_dyn[ii];
          pre_ampsq += 4.*pi*as*std::norm(pre_amp);
        }
        sigmaij[Nf-ii][Nf+ii] = fac_overall*pre_ampsq;
// qg channel
        pre_ampsq = 0.;
        col_avg = 1./(Nc*(Nc*Nc-1.));
        for(unsigned int jj=0; jj<8; jj++){
          pre_amp = ampkinfac_qg[jj]*amp_dyn[ii];
          pre_ampsq += 4.*pi*as*std::norm(pre_amp);
        }
        sigmaij[Nf+ii][Nf] = fac_overall*pre_ampsq;
// gq channel
        pre_ampsq = 0.;
        col_avg = 1./(Nc*(Nc*Nc-1.));
        for(unsigned int jj=0; jj<8; jj++){
          pre_amp = ampkinfac_gq[jj]*amp_dyn[ii];
          pre_ampsq += 4.*pi*as*std::norm(pre_amp);
        }
        sigmaij[Nf][Nf+ii] = fac_overall*pre_ampsq;
// gqb channel
        pre_ampsq = 0.;
        col_avg = 1./(Nc*(Nc*Nc-1.));
        for(unsigned int jj=0; jj<8; jj++){
          pre_amp = ampkinfac_gqb[jj]*amp_dyn[ii];
          pre_ampsq += 4.*pi*as*std::norm(pre_amp);
        }
        sigmaij[Nf][Nf-ii] = fac_overall*pre_ampsq;
// qbg channel
        pre_ampsq = 0.;
        col_avg = 1./(Nc*(Nc*Nc-1.));
        for(unsigned int jj=0; jj<8; jj++){
          pre_amp = ampkinfac_qbg[jj]*amp_dyn[ii];
          pre_ampsq += 4.*pi*as*std::norm(pre_amp);
        }
        sigmaij[Nf-ii][Nf] = fac_overall*pre_ampsq;
      }

    }
    else{
      std::cout << "Wrong DYprocess flag -- STOPPING" << std::endl;
      exit(EXIT_FAILURE);
    }

    //     //For now set sigmaij to just include the PS fac * the spincolavg as want to test PS vol
    // for(int ii=0; ii<2*Nf; ii++){
    //   for(int ij=0; ij<2*Nf; ij++){
    // 	// std::cout << "Nf = " << Nf << std::endl;
    // 	// std::cout << ii << " " << ij << std::endl;
    // 	sigmaij[ii][ij] = 1.0;
    // 	//sigmaij[ii][ij] = 1/(32*pi*pi);
    //   }
    // }


}



void drellyan_amplitude::ampkinfacs(int i1, int i2, int i3, int i4, int i5, std::vector<std::complex<double> >& ampfac){

  // for(int iii=1;iii<=5;iii++) {
  //   for(int jjj=1;jjj<=5;jjj++) {
  //     // std::cout << "sa(" << iii << "," << jjj < ")=" << std::real(sa(iii,jjj))<< "+i*" << std::imag(sa(iii,jjj)) << std::endl;
  //     // std::cout << "sb(" << iii << "," << jjj < ")=" << sb(iii,jjj).real() << std::endl;
  //   }
  // }
  // for(int iii=1;iii<=5;iii++) {
  //   for(int jjj=iii;jjj<=5;jjj++) {
  //     std::cout << "sa(" << iii << "," << jjj << ") = " << std::real(sa(iii,jjj)) << "+i*" << std::imag(sa(iii,jjj)) << std::endl;
  //     std::cout << "sb(" << iii << "," << jjj << ") = " << std::real(sb(iii,jjj)) << "+i*" << std::imag(sb(iii,jjj)) << std::endl;
  //     std::cout << "sa^2( " << iii << "," << jjj << ") = " << std::norm(sa(iii,jjj)) << std::endl;
  //     std::cout << "sb^2( " << iii << "," << jjj << ") = " << std::norm(sb(iii,jjj)) << std::endl;
  //   }
  // }

// helicity - qqb: -1, 1, l-l+: -1, 1, g: -1;
      ampfac[0] = (sa(i3,i5)*sb(i2,i4))/sb(i1,i5) + (sa(i1,i3)*sb(i1,i2)*sb(i2,i4))/(sb(i1,i5)*sb(i2,i5));

// helicity - qqb: -1, 1, l-l+: -1, 1, g: 1;
      ampfac[1] = -((sa(i1,i2)*sa(i1,i3)*sb(i2,i4))/(sa(i1,i5)*sa(i2,i5))) + (sa(i1,i3)*sb(i4,i5))/sa(i2,i5);

// helicity - qqb: -1, 1, l-l+: 1, -1, g: -1;
      ampfac[2] = -((sa(i4,i5)*sb(i2,i3))/sb(i1,i5)) - (sa(i1,i4)*sb(i1,i2)*sb(i2,i3))/(sb(i1,i5)*sb(i2,i5));

// helicity - qqb: -1, 1, l-l+: 1, -1, g: 1;
      ampfac[3] = (sa(i1,i2)*sa(i1,i4)*sb(i2,i3))/(sa(i1,i5)*sa(i2,i5)) - (sa(i1,i4)*sb(i3,i5))/sa(i2,i5);

// helicity - qqb: 1, -1, l-l+: -1, 1, g: -1;
      ampfac[4] = (sa(i3,i5)*sb(i1,i4))/sb(i2,i5) - (sa(i2,i3)*sb(i1,i2)*sb(i1,i4))/(sb(i1,i5)*sb(i2,i5));

// helicity - qqb: 1, -1, l-l+: -1, 1, g: 1;
      ampfac[5] = (sa(i1,i2)*sa(i2,i3)*sb(i1,i4))/(sa(i1,i5)*sa(i2,i5)) + (sa(i2,i3)*sb(i4,i5))/sa(i1,i5);

// helicity - qqb: 1, -1, l-l+: 1, -1, g: -1;
      ampfac[6] = -((sa(i4,i5)*sb(i1,i3))/sb(i2,i5)) + (sa(i2,i4)*sb(i1,i2)*sb(i1,i3))/(sb(i1,i5)*sb(i2,i5));

// helicity - qqb: 1, -1, l-l+: 1, -1, g: 1;
      ampfac[7] = -((sa(i1,i2)*sa(i2,i4)*sb(i1,i3))/(sa(i1,i5)*sa(i2,i5))) - (sa(i2,i4)*sb(i3,i5))/sa(i1,i5);

      for(unsigned int ii = 0; ii<8; ii++) ampfac[ii] *= 2*std::sqrt(2.);
      
}
