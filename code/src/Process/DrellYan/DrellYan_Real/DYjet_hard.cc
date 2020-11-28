
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
    // std::cout << "2*ss(3,4) = " << 2*ss(3,4) << std::endl;
    // std::cout << "2*ss(4,5) = " << 2*ss(4,5) << std::endl;
    // std::cout << "2*ss(3,5) = " << 2*ss(3,5) << std::endl;
    double gevpb = drellyan_in->res_1.gevm2topb;
    double pi = k_constants::pi;
    // std::cout << "gevpb = " << gevpb << std::endl;
// Overall factor for amplitudes
    double flux = 1./2./q2;
    double spinavg = 1/4.;
    int Nc = drellyan_in->res_1.Nc;
    double colfac = (Nc*Nc-1.)/2.;
    //
    double fac_overall = flux*spinavg*colfac*gevpb; // Note that colour average is channel-dependent and not included here; coupling constants are also added later.
    // std::cout << "flux = " << flux << std::endl;
    // std::cout << "fac_overall = " << fac_overall << std::endl;
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
    // std::cout << "ee = " << ee << " ee^2 = " << ee*ee << std::endl;
    // std::cout << "cLzl = " << cLzl << std::endl;
    // std::cout << "gLzQ = " << gLzQ[0] << " " << gLzQ[1] << " " << gLzQ[2] << " " << gLzQ[3] << " " << gLzQ[4] << std::endl;
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
    // double PropDen_y = 1./q2;
    // std::complex<double> PropDen_Z = 1./(q2 - mz*mz + II*mz*zw);
    // std::complex<double> PropDen_Zp = 1./(q2 - mzp*mzp + II*mzp*zpw);
    //TCRIDGE CORRECTED PROPAGATORS WITH S(3,4) IN PLACE OF q2:
    double PropDen_y = 1./(2*ss(3,4));
    std::complex<double> PropDen_Z = 1./(2*ss(3,4) - mz*mz + II*mz*zw);
    std::complex<double> PropDen_Zp = 1./(2*ss(3,4) - mzp*mzp + II*mzp*zpw);
    // std::cout << "q2 = " << q2 << " 2*s(3,4) = " << 2*ss(3,4) << std::endl;
    // std::cout << "PropDen_y = " << PropDen_y << std::endl;
    // std::cout << "PropDen_Z = " << PropDen_Z << std::endl;

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

      // std::vector<std::complex<double> > ampkinfac_qqb(8), ampkinfac_qbq(8), ampkinfac_qg(8), ampkinfac_gq(8), ampkinfac_gqb(8), ampkinfac_qbg(8);
// See member function definition for an explanation of the notation
      // this->ampkinfacs(1, 2, -3, -4, -5, ampkinfac_qqb);
      // this->ampkinfacs(2, 1, -3, -4, -5, ampkinfac_qbq);
      // this->ampkinfacs(1, 5, -3, -4, -2, ampkinfac_qg);
      // this->ampkinfacs(5, 1, -3, -4, -2, ampkinfac_gq);
      // this->ampkinfacs(5, 2, -3, -4, -1, ampkinfac_gqb);
      // this->ampkinfacs(2, 5, -3, -4, -1, ampkinfac_qbg);

      std::vector<std::complex<double>> ampkinfactom_qqbar(8), ampkinfactom_qbarq(8), ampkinfactom_gq(8), ampkinfactom_gqbar(8), ampkinfactom_qg(8), ampkinfactom_qbarg(8);  
      this->ampkinfacs_tom(ampkinfactom_qqbar,1);
      this->ampkinfacs_tom(ampkinfactom_qbarq,2);
      this->ampkinfacs_tom(ampkinfactom_gq,3);
      this->ampkinfacs_tom(ampkinfactom_qg,4);
      this->ampkinfacs_tom(ampkinfactom_gqbar,5);
      this->ampkinfacs_tom(ampkinfactom_qbarg,6);


      // for (int i=0; i<8;i++) {
      // 	// std::cout << "ampkinfactom_qqbar" << "[" << i << "] = " << ampkinfactom_qqbar[i] << std::endl;
      // 	// std::cout << "ampkinfactom_qqbar[" << i << "]/(2sqrt(2))=" << ampkinfactom_qqbar[i]/(2*pow(2,0.5)) << std::endl;
      // }
      // // for (int i=0; i<8;i++) {
      // // 	std::cout << "ampkinfac_qbarq" << "[" << i << "] = " << ampkinfac_qbarq[i] << std::endl;
      // // }
      // for (int i=0; i<8;i++) {
      // 	// std::cout << "ampkinfactom_qbarq" << "[" << i << "]/(2sqrt(2)) = " << ampkinfactom_qbarq[i]/(2*pow(2,0.5)) << std::endl;
      // }
      // // for (int i=0; i<8;i++) {
      // // 	std::cout << "ampkinfac_gq" << "[" << i << "] = " << ampkinfac_gq[i] << std::endl;
      // // }
      // for (int i=0; i<8;i++) {
      // 	// std::cout << "ampkinfactom_gq" << "[" << i << "]/(2sqrt(2)) = " << ampkinfactom_gq[i]/(2*pow(2,0.5)) << std::endl;
      // }
      // // for (int i=0; i<8;i++) {
      // // 	std::cout << "ampkinfac_qg" << "[" << i << "] = " << ampkinfac_qg[i] << std::endl;
      // // }
      // for (int i=0; i<8;i++) {
      // 	std::cout << "ampkinfactom_qg" << "[" << i << "]/(2sqrt(2)) = " << ampkinfactom_qg[i]/(2*pow(2,0.5)) << std::endl;
      // }
      // // for (int i=0; i<8;i++) {
      // // 	std::cout << "ampkinfac_gqbar" << "[" << i << "] = " << ampkinfac_gqb[i] << std::endl;
      // // }
      // for (int i=0; i<8;i++) {
      // 	std::cout << "ampkinfactom_gqbar" << "[" << i << "]/(2sqrt(2)) = " << ampkinfactom_gqbar[i]/(2*pow(2,0.5)) << std::endl;
      // }
      // // for (int i=0; i<8;i++) {
      // // 	std::cout << "ampkinfac_qbarg" << "[" << i << "] = " << ampkinfac_qbg[i] << std::endl;
      // // }
      // for (int i=0; i<8;i++) {
      // 	std::cout << "ampkinfactom_qbarg" << "[" << i << "]/(2sqrt(2)) = " << ampkinfactom_qbarg[i]/(2*pow(2,0.5)) << std::endl;
      // }

      std::vector<std::complex<double> > amp_dynLL(Nf);
      std::vector<std::complex<double> > amp_dynRR(Nf);
      std::vector<std::complex<double> > amp_dynLR(Nf);
      std::vector<std::complex<double> > amp_dynRL(Nf);
      double pre_amp;
      double pre_ampsq, col_avg;

      if (DYproc == 4) {
        for(unsigned int ii=0; ii<Nf; ii++){
          amp_dynLL[ii] = gLzQ[ii]*cLzl*PropDen_Z;
          amp_dynRR[ii] = gRzQ[ii]*cRzl*PropDen_Z;
          amp_dynLR[ii] = gLzQ[ii]*cRzl*PropDen_Z;
          amp_dynRL[ii] = gRzQ[ii]*cLzl*PropDen_Z;
        }
      }
// gamma* + Z -- subprocess 5
      else if (DYproc == 5) {
        for(unsigned int ii=0; ii<Nf; ii++){
          amp_dynLL[ii] = -QQ[ii]*ee*ee*PropDen_y + gLzQ[ii]*cLzl*PropDen_Z;
          amp_dynRR[ii] = -QQ[ii]*ee*ee*PropDen_y + gRzQ[ii]*cRzl*PropDen_Z;
          amp_dynLR[ii] = -QQ[ii]*ee*ee*PropDen_y + gLzQ[ii]*cRzl*PropDen_Z;
          amp_dynRL[ii] = -QQ[ii]*ee*ee*PropDen_y + gRzQ[ii]*cLzl*PropDen_Z;
	  // std::cout << "amp_dyn[" << ii << "]=" << amp_dyn[ii] << std::endl;
	  // std::cout << "amp_dyn[" << ii << "]/ee^2=" << amp_dyn[ii]/(ee*ee) << std::endl;
	  // std::cout << "-QQ*ee*ee*PropDen_y = " << -QQ[ii]*ee*ee*PropDen_y << std::endl;
	  // std::cout << "-QQ*PropDen_y = " << -QQ[ii]*PropDen_y << std::endl;
	  // std::cout << "QQ[" << ii << "]=" << QQ[ii] << " ee^2 = " << ee*ee << std::endl;
	  // std::cout << "gLzQ[ii]*cLzl*PropDen_Z = " << gLzQ[ii]*cLzl*PropDen_Z << std::endl;
	  // std::cout << "gLzQ[ii]*cLzl*PropDen_Z/ee^2 = " << gLzQ[ii]*cLzl*PropDen_Z/(ee*ee) << std::endl;
	  // std::cout << "gLzQ[" << ii << "] = " << gLzQ[ii] << "cLzl = " << cLzl << std::endl;
        }
      }
//
// gamma* + Z + Z' -- subprocess 6
      else if (DYproc == 6) {
        for(unsigned int ii=0; ii<Nf; ii++){
          amp_dynLL[ii] = -QQ[ii]*ee*ee*PropDen_y + gLzQ[ii]*cLzl*PropDen_Z + gLZpQ[ii]*cLZpl*PropDen_Zp;
          amp_dynRR[ii] = -QQ[ii]*ee*ee*PropDen_y + gRzQ[ii]*cRzl*PropDen_Z + gRZpQ[ii]*cRZpl*PropDen_Zp;
          amp_dynLR[ii] = -QQ[ii]*ee*ee*PropDen_y + gLzQ[ii]*cRzl*PropDen_Z + gLZpQ[ii]*cRZpl*PropDen_Zp;
          amp_dynRL[ii] = -QQ[ii]*ee*ee*PropDen_y + gRzQ[ii]*cLzl*PropDen_Z + gRZpQ[ii]*cLZpl*PropDen_Zp;

        }
      }

      //Must only combine the LL pieces with the LL coupling propagator etc, e.g. qLqRbar -> lL-lR+ g- and g+ only can go with the gLzQ*cLzl coupling dynamics piece etc
      std::vector<double> ampkinfacsq_tom_qqbar(4),ampkinfacsq_tom_qbarq(4),ampkinfacsq_tom_qg(4),ampkinfacsq_tom_gq(4),ampkinfacsq_tom_gqbar(4),ampkinfacsq_tom_qbarg(4);
      ampkinfacssq_tom(ampkinfacsq_tom_qqbar,ampkinfactom_qqbar);
      ampkinfacssq_tom(ampkinfacsq_tom_qbarq,ampkinfactom_qbarq);
      ampkinfacssq_tom(ampkinfacsq_tom_qg,ampkinfactom_qg);
      ampkinfacssq_tom(ampkinfacsq_tom_gq,ampkinfactom_gq);
      ampkinfacssq_tom(ampkinfacsq_tom_qbarg,ampkinfactom_qbarg);
      ampkinfacssq_tom(ampkinfacsq_tom_gqbar,ampkinfactom_gqbar);
      //Above 6 lines assign LL,LR,RL,RR components to elements 0,1,2,3 respectively

      // // std::cout << "ampkinfacsq_tom_qqbar[0]=" << ampkinfacsq_tom_qqbar[0] << std::endl;
      // // std::cout << "ampkinfacsq_tom_qqbar[1]=" << ampkinfacsq_tom_qqbar[1] << std::endl;
      // // std::cout << "ampkinfacsq_tom_qqbar[2]=" << ampkinfacsq_tom_qqbar[2] << std::endl;
      // // std::cout << "ampkinfacsq_tom_qqbar[3]=" << ampkinfacsq_tom_qqbar[3] << std::endl;
      // std::cout << "ampkinfacsq_tom_qqbar[0]/8=" << ampkinfacsq_tom_qqbar[0]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qqbar[1]/8=" << ampkinfacsq_tom_qqbar[1]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qqbar[2]/8=" << ampkinfacsq_tom_qqbar[2]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qqbar[3]/8=" << ampkinfacsq_tom_qqbar[3]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qg[0]/8=" << ampkinfacsq_tom_qg[0]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qg[1]/8=" << ampkinfacsq_tom_qg[1]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qg[2]/8=" << ampkinfacsq_tom_qg[2]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qg[3]/8=" << ampkinfacsq_tom_qg[3]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_gq[0]/8=" << ampkinfacsq_tom_gq[0]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_gq[1]/8=" << ampkinfacsq_tom_gq[1]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_gq[2]/8=" << ampkinfacsq_tom_gq[2]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_gq[3]/8=" << ampkinfacsq_tom_gq[3]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qbarq[0]/8=" << ampkinfacsq_tom_qbarq[0]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qbarq[1]/8=" << ampkinfacsq_tom_qbarq[1]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qbarq[2]/8=" << ampkinfacsq_tom_qbarq[2]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qbarq[3]/8=" << ampkinfacsq_tom_qbarq[3]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qbarg[0]/8=" << ampkinfacsq_tom_qbarg[0]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qbarg[1]/8=" << ampkinfacsq_tom_qbarg[1]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qbarg[2]/8=" << ampkinfacsq_tom_qbarg[2]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_qbarg[3]/8=" << ampkinfacsq_tom_qbarg[3]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_gqbar[0]/8=" << ampkinfacsq_tom_gqbar[0]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_gqbar[1]/8=" << ampkinfacsq_tom_gqbar[1]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_gqbar[2]/8=" << ampkinfacsq_tom_gqbar[2]/8 << std::endl;
      // std::cout << "ampkinfacsq_tom_gqbar[3]/8=" << ampkinfacsq_tom_gqbar[3]/8 << std::endl;

      

// Put pieces together
      // std::cout << "Nf = " << Nf << std::endl;
      for(unsigned int ii=1; ii<=Nf; ii++){
	// std::cout << "STARTING ii: " << ii << std::endl;
// qqb channel
        pre_ampsq = 0.;
        col_avg = 1./(Nc*Nc);
        // for(unsigned int jj=0; jj<8; jj++){
	  // std::cout << "amp_dynLL[" << ii << "]/ee^2=" << amp_dynLL[ii]/(ee*ee) << std::endl;
	  // std::cout << "amp_dynRR[" << ii << "]/ee^2=" << amp_dynRR[ii]/(ee*ee) << std::endl;
	  // std::cout << "amp_dynLR[" << ii << "]/ee^2=" << amp_dynLR[ii]/(ee*ee) << std::endl;
	  // std::cout << "amp_dynRL[" << ii << "]/ee^2=" << amp_dynRL[ii]/(ee*ee) << std::endl;
	  // std::cout << "ampkinfacsq_tom_qqbar[0]*norm(amp_dynLL[" << ii << "]) = " << ampkinfacsq_tom_qqbar[0]*norm(amp_dynLL[ii]) << std::endl;
	  // std::cout << "ampkinfacsq_tom_qqbar[1]*norm(amp_dynLR[" << ii << "]) = " << ampkinfacsq_tom_qqbar[1]*norm(amp_dynLR[ii]) << std::endl;
	  // std::cout << "ampkinfacsq_tom_qqbar[2]*norm(amp_dynRL[" << ii << "]) = " << ampkinfacsq_tom_qqbar[2]*norm(amp_dynRL[ii]) << std::endl;
	  // std::cout << "ampkinfacsq_tom_qqbar[3]*norm(amp_dynRR[" << ii << "]) = " << ampkinfacsq_tom_qqbar[3]*norm(amp_dynRR[ii]) << std::endl;
	  // std::cout << "ampkinfacsq_tom_qqbar[0]*norm(amp_dynLL[" << ii << "])/8ee^4 = " << ampkinfacsq_tom_qqbar[0]*norm(amp_dynLL[ii])/(8*ee*ee*ee*ee) << std::endl;
	  // std::cout << "ampkinfacsq_tom_qqbar[1]*norm(amp_dynLR[" << ii << "])/8ee^4 = " << ampkinfacsq_tom_qqbar[1]*norm(amp_dynLR[ii])/(8*ee*ee*ee*ee) << std::endl;
	  // std::cout << "ampkinfacsq_tom_qqbar[2]*norm(amp_dynRL[" << ii << "])/8ee^4 = " << ampkinfacsq_tom_qqbar[2]*norm(amp_dynRL[ii])/(8*ee*ee*ee*ee) << std::endl;
	  // std::cout << "ampkinfacsq_tom_qqbar[3]*norm(amp_dynRR[" << ii << "])/8ee^4 = " << ampkinfacsq_tom_qqbar[3]*norm(amp_dynRR[ii])/(8*ee*ee*ee*ee) << std::endl;
	// std::cout << "qqbar channel:" << std::endl;
	pre_amp = ampkinfacsq_tom_qqbar[0]*std::norm(amp_dynLL[ii-1])+ampkinfacsq_tom_qqbar[1]*std::norm(amp_dynLR[ii-1])+ampkinfacsq_tom_qqbar[2]*std::norm(amp_dynRL[ii-1])+ampkinfacsq_tom_qqbar[3]*std::norm(amp_dynRR[ii-1]);
	// std::cout << "ii: " << ii << " pre_amp/(8*ee^4) = (" << pre_amp/(8*ee*ee*ee*ee) << std::endl;
	pre_ampsq += 4.*pi*as*pre_amp;
	// std::cout << "as = " << as << std::endl;
	// std::cout << "as*pi = " << as*pi << std::endl;
	// std::cout << "pre_ampsq = " << pre_ampsq << std::endl;
        // }
	// std::cout << "fac_pverall = " << fac_overall << std::endl;
	// std::cout << "4*pi*alphas*gevpb/2q2 = " << 4*pi*as*gevpb/(2*q2) << std::endl;
	// std::cout << "flux = " << flux << std::endl;
	// std::cout << "4*pi*alphas = " << 4*pi*as << std::endl;
        // sigmaij[Nf+ii][Nf-ii] = fac_overall*pre_ampsq*col_avg;
	// std::cout << "sigmaij[" << Nf+ii << "][" << Nf-ii << "]=" << sigmaij[Nf+ii][Nf-ii] << std::endl;
	// std::cout << "sigmaij[" << Nf+ii << "][" << Nf-ii << "]/(flux*gevpb)*pi=" << sigmaij[Nf+ii][Nf-ii]/(flux*gevpb)*pi << std::endl;

// qbq channel
	// std::cout << "qbarq channel:" << std::endl;
        pre_ampsq = 0.;
        col_avg = 1./(Nc*Nc);
        // for(unsigned int jj=0; jj<8; jj++){
	pre_amp = ampkinfacsq_tom_qbarq[0]*std::norm(amp_dynLL[ii-1])+ampkinfacsq_tom_qbarq[1]*std::norm(amp_dynLR[ii-1])+ampkinfacsq_tom_qbarq[2]*std::norm(amp_dynRL[ii-1])+ampkinfacsq_tom_qbarq[3]*std::norm(amp_dynRR[ii-1]);
	  // std::cout << "ii: " << ii << " pre_amp/(8*ee^4) = (" << pre_amp/(8*ee*ee*ee*ee) << std::endl;
	pre_ampsq += 4.*pi*as*pre_amp;
        // }
        sigmaij[Nf-ii][Nf+ii] = fac_overall*pre_ampsq*col_avg;
	// std::cout << "sigmaij[" << Nf-ii << "][" << Nf+ii << "]=" << sigmaij[Nf-ii][Nf+ii] << std::endl;
	// std::cout << "sigmaij[" << Nf-ii << "][" << Nf+ii << "]/(flux*gevpb)*pi=" << sigmaij[Nf-ii][Nf+ii]/(flux*gevpb)*pi << std::endl;
// qg channel
	// std::cout << "qg channel:" << std::endl;
        pre_ampsq = 0.;
        col_avg = 1./(Nc*(Nc*Nc-1.));
	// for(unsigned int jj=0; jj<8; jj++){
	pre_amp = ampkinfacsq_tom_qg[0]*std::norm(amp_dynLL[ii-1])+ampkinfacsq_tom_qg[1]*std::norm(amp_dynLR[ii-1])+ampkinfacsq_tom_qg[2]*std::norm(amp_dynRL[ii-1])+ampkinfacsq_tom_qg[3]*std::norm(amp_dynRR[ii-1]);
	// std::cout << "ii: " << ii << " pre_amp/(8*ee^4) = (" << pre_amp/(8*ee*ee*ee*ee) << std::endl;
	pre_ampsq += 4.*pi*as*pre_amp;
        // }
        sigmaij[Nf+ii][Nf] = fac_overall*pre_ampsq*col_avg;
	// std::cout << "sigmaij[" << Nf+ii << "][" << Nf << "]=" << sigmaij[Nf+ii][Nf] << std::endl;
	// std::cout << "sigmaij[" << Nf+ii << "][" << Nf << "]/(flux*gevpb)*pi=" << sigmaij[Nf+ii][Nf]/(flux*gevpb)*pi << std::endl;
// gq channel
	// std::cout << "gq channel:" << std::endl;
        pre_ampsq = 0.;
        col_avg = 1./(Nc*(Nc*Nc-1.));
	// for(unsigned int jj=0; jj<8; jj++){
          pre_amp = ampkinfacsq_tom_gq[0]*std::norm(amp_dynLL[ii-1])+ampkinfacsq_tom_gq[1]*std::norm(amp_dynLR[ii-1])+ampkinfacsq_tom_gq[2]*std::norm(amp_dynRL[ii-1])+ampkinfacsq_tom_gq[3]*std::norm(amp_dynRR[ii-1]);
	  // std::cout << "ii: " << ii << " pre_amp/(8*ee^4) = (" << pre_amp/(8*ee*ee*ee*ee) << std::endl;
	pre_ampsq += 4.*pi*as*pre_amp;
        // }
        sigmaij[Nf][Nf+ii] = fac_overall*pre_ampsq*col_avg;
	// std::cout << "sigmaij[" << Nf << "][" << Nf+ii << "]=" << sigmaij[Nf][Nf+ii] << std::endl;
	// std::cout << "sigmaij[" << Nf << "][" << Nf+ii << "]/(flux*gevpb)*pi=" << sigmaij[Nf][Nf+ii]/(flux*gevpb)*pi << std::endl;
// gqb channel
	// std::cout << "gqb channel:" << std::endl;
        pre_ampsq = 0.;
        col_avg = 1./(Nc*(Nc*Nc-1.));
	// for(unsigned int jj=0; jj<8; jj++){
	pre_amp = ampkinfacsq_tom_gqbar[0]*std::norm(amp_dynLL[ii-1])+ampkinfacsq_tom_gqbar[1]*std::norm(amp_dynLR[ii-1])+ampkinfacsq_tom_gqbar[2]*std::norm(amp_dynRL[ii-1])+ampkinfacsq_tom_gqbar[3]*std::norm(amp_dynRR[ii]);
	// std::cout << "ii: " << ii << " pre_amp/(8*ee^4) = (" << pre_amp/(8*ee*ee*ee*ee) << std::endl;
	pre_ampsq += 4.*pi*as*pre_amp;
        // }
        sigmaij[Nf][Nf-ii] = fac_overall*pre_ampsq*col_avg;
	// std::cout << "sigmaij[" << Nf << "][" << Nf-ii << "]=" << sigmaij[Nf][Nf-ii] << std::endl;
	// std::cout << "sigmaij[" << Nf << "][" << Nf-ii << "]/(flux*gevpb)*pi=" << sigmaij[Nf][Nf-ii]/(flux*gevpb)*pi << std::endl;
// qbg channel
	// std::cout << "qbg channel:" << std::endl;
        pre_ampsq = 0.;
        col_avg = 1./(Nc*(Nc*Nc-1.));
	// for(unsigned int jj=0; jj<8; jj++){
	pre_amp = ampkinfacsq_tom_qbarg[0]*std::norm(amp_dynLL[ii-1])+ampkinfacsq_tom_qbarg[1]*std::norm(amp_dynLR[ii-1])+ampkinfacsq_tom_qbarg[2]*std::norm(amp_dynRL[ii-1])+ampkinfacsq_tom_qbarg[3]*std::norm(amp_dynRR[ii-1]);
	// std::cout << "ii: " << ii << " pre_amp/(8*ee^4) = (" << pre_amp/(8*ee*ee*ee*ee) << std::endl;
	pre_ampsq += 4.*pi*as*pre_amp;
        // }
        sigmaij[Nf-ii][Nf] = fac_overall*pre_ampsq*col_avg;
	// std::cout << "sigmaij[" << Nf-ii << "][" << Nf << "]=" << sigmaij[Nf-ii][Nf] << std::endl;
	// std::cout << "sigmaij[" << Nf-ii << "][" << Nf << "]/(flux*gevpb)*pi=" << sigmaij[Nf-ii][Nf]/(flux*gevpb)*pi << std::endl;
	// std::cout << "ENDING ii: " << ii << std::endl;
      }	
    }
      
    else{
      std::cout << "Wrong DYprocess flag -- STOPPING" << std::endl;
      exit(EXIT_FAILURE);
    }

    //     //For now set sigmaij to just include the PS fac * the spincolavg as want to test PS vol
    // for(int ii=0; ii<=2*Nf; ii++){
    //   for(int ij=0; ij<=2*Nf; ij++){
    // // 	// std::cout << "Nf = " << Nf << std::endl;
    // // 	// std::cout << ii << " " << ij << std::endl;
    // // 	sigmaij[ii][ij] = 1.0;
    // // 	//sigmaij[ii][ij] = 1/(32*pi*pi);
    // 	if (sigmaij[ii][ij]!=0) {
    // 	  // std::cout << "sigmaij[" << ii << "][" << ij << "]=" << sigmaij[ii][ij] << std::endl;
    // 	  std::cout << "sigmaij[" << ii << "][" << ij << "]/(flux*gevpb)*pi=" << sigmaij[ii][ij]/(flux*gevpb)*pi << std::endl;
    // 	}
    //   }
    // }


}

void drellyan_amplitude::ampkinfacssq_tom(std::vector<double> & ampfacsq_tom_vec, std::vector<std::complex<double>> ampfac_tom_vec) {
  ampfacsq_tom_vec[0] = std::norm(ampfac_tom_vec[0])+std::norm(ampfac_tom_vec[4]); //LL CONTRIBUTION
  ampfacsq_tom_vec[3] = std::norm(ampfac_tom_vec[3])+std::norm(ampfac_tom_vec[7]); //RR CONTRIBUTION
  ampfacsq_tom_vec[1] = std::norm(ampfac_tom_vec[2])+std::norm(ampfac_tom_vec[6]); //LR CONTRIBUTION 
  ampfacsq_tom_vec[2] = std::norm(ampfac_tom_vec[1])+std::norm(ampfac_tom_vec[5]); //RL CONTRIBUTION

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

// helicity - qqb: -1, 1, l-l+: -1, 1, g: -1; -WRONG!
// helicity - qqb: 1, -1 l-l+: -1, 1, g: -1; -RIGHT! SO RL CONTRIBUTION
      ampfac[1] = (sa(i3,i5)*sb(i2,i4))/sb(i1,i5) + (sa(i1,i3)*sb(i1,i2)*sb(i2,i4))/(sb(i1,i5)*sb(i2,i5));

// helicity - qqb: -1, 1, l-l+: -1, 1, g: 1; -WRONG!
// helicity - qqb: 1, -1, l-l+: -1, 1, g: 1; -RIGHT! SO RL CONTRIBUTION
      ampfac[2] = -((sa(i1,i2)*sa(i1,i3)*sb(i2,i4))/(sa(i1,i5)*sa(i2,i5))) + (sa(i1,i3)*sb(i4,i5))/sa(i2,i5);

// helicity - qqb: -1, 1, l-l+: 1, -1, g: -1; -WRONG!
// helicity - qqb: 1, -1, l-l+: 1, -1, g: -1; -RIGHT! SO RR CONTRIBUTION
      ampfac[3] = -((sa(i4,i5)*sb(i2,i3))/sb(i1,i5)) - (sa(i1,i4)*sb(i1,i2)*sb(i2,i3))/(sb(i1,i5)*sb(i2,i5));

// helicity - qqb: -1, 1, l-l+: 1, -1, g: 1; -WRONG!
// helicity - qqb: 1, -1, l-l+: 1, -1, g: 1; -RIGHT! SO RR CONTRIBUTION
      ampfac[0] = (sa(i1,i2)*sa(i1,i4)*sb(i2,i3))/(sa(i1,i5)*sa(i2,i5)) - (sa(i1,i4)*sb(i3,i5))/sa(i2,i5);

// helicity - qqb: 1, -1, l-l+: -1, 1, g: -1; -WRONG!
// helicity - qqb: -1, 1, l-l+: -1, 1, g: -1; -RIGHT! SO LL CONTRIBUTION
      ampfac[7] = (sa(i3,i5)*sb(i1,i4))/sb(i2,i5) - (sa(i2,i3)*sb(i1,i2)*sb(i1,i4))/(sb(i1,i5)*sb(i2,i5));

// helicity - qqb: 1, -1, l-l+: -1, 1, g: 1; -WRONG!
// helicity - qqb: -1, 1, l-l+: -1, 1, g: 1; -RIGHT! SO LL CONTRIBUTION
      ampfac[4] = (sa(i1,i2)*sa(i2,i3)*sb(i1,i4))/(sa(i1,i5)*sa(i2,i5)) + (sa(i2,i3)*sb(i4,i5))/sa(i1,i5);

// helicity - qqb: 1, -1, l-l+: 1, -1, g: -1; -WRONG!
// helicity - qqb: -1, 1, l-l+: 1, -1, g: -1; -RIGHT! SO LR CONTRIBUTION
      ampfac[5] = -((sa(i4,i5)*sb(i1,i3))/sb(i2,i5)) + (sa(i2,i4)*sb(i1,i2)*sb(i1,i3))/(sb(i1,i5)*sb(i2,i5));

// helicity - qqb: 1, -1, l-l+: 1, -1, g: 1; -WRONG!
// helicity - qqb: -1, 1, l-l+: 1, -1, g: 1; -RIGHT! SO LR CONTRIBUTION
      ampfac[6] = -((sa(i1,i2)*sa(i2,i4)*sb(i1,i3))/(sa(i1,i5)*sa(i2,i5))) - (sa(i2,i4)*sb(i3,i5))/sa(i1,i5);

      for(unsigned int ii = 0; ii<8; ii++) ampfac[ii] *= 2*std::sqrt(2.);

      
  // std::cout << "sa(1,2) = " << sa(1,2) << " sb(1,2) = " << sb(1,2) << std::endl;
      
}

void drellyan_amplitude::ampkinfacs_tom(std::vector<std::complex<double>> & ampfac_tom_vec, int proc){
  if (proc == 1) { //qqbar initiated 
    //M1 qLqRbar -> lL- lR+ g-: i1 = 1, i2 = 2, i3 = 3, i4 = 4, i5 = 5
    ampfac_tom_vec[0] = ampfac_tom(1,2,3,4,5,1,0); //M1 qLqRbar->lL- lR+ g- -LL
    ampfac_tom_vec[1] = ampfac_tom(2,1,3,4,5,1,0); //M2 qRqLbar->lL- lR+ g- -RL
    ampfac_tom_vec[2] = ampfac_tom(1,2,4,3,5,1,0); //M3 qLqRbar->lR- lL+ g- -LR
    ampfac_tom_vec[3] = ampfac_tom(2,1,4,3,5,1,0); //M4 qRqLbar->lR- lL+ g- -RR
    ampfac_tom_vec[4] = ampfac_tom(2,1,4,3,5,-1,0); //M5 qLqRbar->lL- lR+ g+ -LL
    ampfac_tom_vec[5] = ampfac_tom(1,2,4,3,5,-1,0); //M6 qRqLbar->lL- lR+ g+ -RL
    ampfac_tom_vec[6] = ampfac_tom(2,1,3,4,5,-1,0); //M7 qLqRbar->lR- lL+ g+ -LR
    ampfac_tom_vec[7] = ampfac_tom(1,2,3,4,5,-1,0); //M8 qRqLbar->lR- lL+ g+ -RR
  }
  else if (proc == 2) {//qbarq initiated
    ampfac_tom_vec[0] = ampfac_tom(2,1,3,4,5,1,0); //M1 qRbarqL->lL- lR+ g-
    ampfac_tom_vec[1] = ampfac_tom(1,2,3,4,5,1,0); //M2 qLbarqR->lL- lR+ g-
    ampfac_tom_vec[2] = ampfac_tom(2,1,4,3,5,1,0); //M3 qRbarqL->lR- lL+ g-
    ampfac_tom_vec[3] = ampfac_tom(1,2,4,3,5,1,0); //M4 qLbarqR->lR- lL+ g-
    ampfac_tom_vec[4] = ampfac_tom(1,2,4,3,5,-1,0); //M5 qRbarqL->lL- lR+ g+
    ampfac_tom_vec[5] = ampfac_tom(2,1,4,3,5,-1,0); //M6 qLbarqR->lL- lR+ g+
    ampfac_tom_vec[6] = ampfac_tom(1,2,3,4,5,-1,0); //M7 qRbarqL->lR- lL+ g+
    ampfac_tom_vec[7] = ampfac_tom(2,1,3,4,5,-1,0); //M8 qLbarqR->lR- lL+ g+
  }
  else if (proc == 3) {//gq initiated, just as qbarq but 1<->5
    ampfac_tom_vec[0] = ampfac_tom(2,5,3,4,1,1,0);
    ampfac_tom_vec[1] = ampfac_tom(5,2,3,4,1,1,0);
    ampfac_tom_vec[2] = ampfac_tom(2,5,4,3,1,1,0);
    ampfac_tom_vec[3] = ampfac_tom(5,2,4,3,1,1,0);
    ampfac_tom_vec[4] = ampfac_tom(5,2,4,3,1,-1,0);
    ampfac_tom_vec[5] = ampfac_tom(2,5,4,3,1,-1,0);
    ampfac_tom_vec[6] = ampfac_tom(5,2,3,4,1,-1,0);
    ampfac_tom_vec[7] = ampfac_tom(2,5,3,4,1,-1,0);

  }
  else if (proc == 4) {//qg initiated new, just as qqbar but 2<->5
    ampfac_tom_vec[0] = ampfac_tom(1,5,3,4,2,1,0); 
    ampfac_tom_vec[1] = ampfac_tom(5,1,3,4,2,1,0); 
    ampfac_tom_vec[2] = ampfac_tom(1,5,4,3,2,1,0); 
    ampfac_tom_vec[3] = ampfac_tom(5,1,4,3,2,1,0); 
    ampfac_tom_vec[4] = ampfac_tom(5,1,4,3,2,-1,0); 
    ampfac_tom_vec[5] = ampfac_tom(1,5,4,3,2,-1,0); 
    ampfac_tom_vec[6] = ampfac_tom(5,1,3,4,2,-1,0); 
    ampfac_tom_vec[7] = ampfac_tom(1,5,3,4,2,-1,0); 

  }
  else if (proc == 5) {//gqbar get from qqbar with 1<->5
    ampfac_tom_vec[0] = ampfac_tom(5,2,3,4,1,1,0); //M1 qLqRbar->lL- lR+ g- -LL
    ampfac_tom_vec[1] = ampfac_tom(2,5,3,4,1,1,0); //M2 qRqLbar->lL- lR+ g- -RL
    ampfac_tom_vec[2] = ampfac_tom(5,2,4,3,1,1,0); //M3 qLqRbar->lR- lL+ g- -LR
    ampfac_tom_vec[3] = ampfac_tom(2,5,4,3,1,1,0); //M4 qRqLbar->lR- lL+ g- -RR
    ampfac_tom_vec[4] = ampfac_tom(2,5,4,3,1,-1,0); //M5 qLqRbar->lL- lR+ g+ -LL
    ampfac_tom_vec[5] = ampfac_tom(5,2,4,3,1,-1,0); //M6 qRqLbar->lL- lR+ g+ -RL
    ampfac_tom_vec[6] = ampfac_tom(2,5,3,4,1,-1,0); //M7 qLqRbar->lR- lL+ g+ -LR
    ampfac_tom_vec[7] = ampfac_tom(5,2,3,4,1,-1,0); //M8 qRqLbar->lR- lL+ g+ -RR

  }  
  else if (proc == 6) {//qbarg get from qbarq with 2<->5
    ampfac_tom_vec[0] = ampfac_tom(5,1,3,4,2,1,0); //M1 qRbarqL->lL- lR+ g-
    ampfac_tom_vec[1] = ampfac_tom(1,5,3,4,2,1,0); //M2 qLbarqR->lL- lR+ g-
    ampfac_tom_vec[2] = ampfac_tom(5,1,4,3,2,1,0); //M3 qRbarqL->lR- lL+ g-
    ampfac_tom_vec[3] = ampfac_tom(1,5,4,3,2,1,0); //M4 qLbarqR->lR- lL+ g-
    ampfac_tom_vec[4] = ampfac_tom(1,5,4,3,2,-1,0); //M5 qRbarqL->lL- lR+ g+
    ampfac_tom_vec[5] = ampfac_tom(5,1,4,3,2,-1,0); //M6 qLbarqR->lL- lR+ g+
    ampfac_tom_vec[6] = ampfac_tom(1,5,3,4,2,-1,0); //M7 qRbarqL->lR- lL+ g+
    ampfac_tom_vec[7] = ampfac_tom(5,1,3,4,2,-1,0); //M8 qLbarqR->lR- lL+ g+
  }
  else {
    std::cout << "ISSUE! - proc must be 1 to 6 in ampfac_tom_vec in DYjet_hard.cc" << std::endl;
  }
  //OLD
  // else if (proc == 3 || proc == 5) {//gq initiated or gqbar initiated
  //   ampfac_tom_vec[0] = ampfac_tom(1,2,3,4,5,1,1); //M1' g+ qL->lL- lR+ qL
  //   ampfac_tom_vec[1] = ampfac_tom(2,1,3,4,5,1,1); //M2' g+ qR->lL- lR+ qR
  //   ampfac_tom_vec[2] = ampfac_tom(1,2,4,3,5,1,1); //M3' g+ qL->lR- lL+ qL
  //   ampfac_tom_vec[3] = ampfac_tom(2,1,4,3,5,1,1); //M4' g+ qR->lR- lL+ qR
  //   ampfac_tom_vec[4] = ampfac_tom(2,1,4,3,5,-1,1); //M5' g- qL->lL- lR+ qL
  //   ampfac_tom_vec[5] = ampfac_tom(1,2,4,3,5,-1,1); //M6' g- qR->lL- lR+ qR
  //   ampfac_tom_vec[6] = ampfac_tom(2,1,3,4,5,-1,1); //M7' g- qL->lR- lL+ qL
  //   ampfac_tom_vec[7] = ampfac_tom(1,2,3,4,5,-1,1); //M8' g- qR->lR- lL+ qR
  // }
  // else if (proc == 4 || proc == 6) {//qg initiated or qbarg initiated
  //   ampfac_tom_vec[0] = ampfac_tom(2,1,3,4,5,1,1); //M1' g+ qL->lL- lR+ qL
  //   ampfac_tom_vec[1] = ampfac_tom(1,2,3,4,5,1,1); //M2' g+ qR->lL- lR+ qR
  //   ampfac_tom_vec[2] = ampfac_tom(2,1,4,3,5,1,1); //M3' g+ qL->lR- lL+ qL
  //   ampfac_tom_vec[3] = ampfac_tom(1,2,4,3,5,1,1); //M4' g+ qR->lR- lL+ qR
  //   ampfac_tom_vec[4] = ampfac_tom(1,2,4,3,5,-1,1); //M5' g- qL->lL- lR+ qL
  //   ampfac_tom_vec[5] = ampfac_tom(2,1,4,3,5,-1,1); //M6' g- qR->lL- lR+ qR
  //   ampfac_tom_vec[6] = ampfac_tom(1,2,3,4,5,-1,1); //M7' g- qL->lR- lL+ qL
  //   ampfac_tom_vec[7] = ampfac_tom(2,1,3,4,5,-1,1); //M8' g- qR->lR- lL+ qR
  // } 

}

std::complex<double> drellyan_amplitude::ampfac_tom(int i1, int i2, int i3, int i4, int i5,int bra_ket, int gq){
  std::complex<double> ampfac_tomout = std::complex<double>(0.0,0.0);

  // std::cout << "sa(1,2) tom = " << sa(1,2) << " sb(1,2) tom = " << sb(1,2) << std::endl;
  double i1temp = 0, i2temp = 0, i5temp = 0;
  //double sign = 0; //drop sign as not needed as once conservation of momentum used to reduce to a single expression (1 term in numerator not 2) then it gives just an overall global minus which is not physical

  for (int ii=1;ii<=5;ii++) {
    for (int jj=1;jj<=5;jj++) {
      // std::cout << "sa(" << ii << ")(" << jj << ")=" << sa(ii,jj) << std::endl;
      // std::cout << "sb(" << ii << ")(" << jj << ")=" << sb(ii,jj) << std::endl;
    }
  }

  if (gq == 1) { //Need to do crossing 1->5, 5->2 and 2->1 and also introduce minus signs for one type of bracket (choose []=sb) for 1 and 5
    //sign = -1;
    i1temp = i1;
    i2temp = i2;
    i5temp = i5;
    i5 = i2temp; //5->2
    i2 = i1temp; //2->1
    i1 = i5temp; //1->5
  }
  else {
    //sign = 1;
  }
  if (bra_ket ==1) {
    //M1 qLqRbar -> lL- lR+ g-: i1 = 1, i2 = 2, i3 = 3, i4 = 4, i5 = 5
    ampfac_tomout = sb(i1,i4)*sb(i1,i4)/(sb(i1,i5)*sb(i2,i5)*sb(i4,i3)); //M1=prefacA*[14]^2/([15]*[25]*[43])
    // std::cout << "M1' = " << sb(4,5)*sb(4,5)/(sb(1,2)*sb(5,1)*sb(4,3))*sa(3,4)*sb(4,3)*2.0*std::sqrt(2.) << std::endl;
  }
  else if (bra_ket ==-1) {
    ampfac_tomout = sa(i1,i4)*sa(i1,i4)/(sa(i1,i5)*sa(i2,i5)*sa(i4,i3)); //M8=prefacA*<14>^2/(<15>*<25>*<43>)
  }
  else {
    std::cout << "ERROR bra_ket must be 1 or -1" << std::endl;
    exit(EXIT_FAILURE);
  }

  //REMOVE <34>[43] as that is in Z propagator elsewhere
  ampfac_tomout = ampfac_tomout*sa(i3,i4)*sb(i4,i3);

  // std::cout << "i1 = " << i1 << " i2 = " << i2 << " i3 = " << i3 << " i4 = " << i4 << " i5 = " << i5 << std::endl;
  // std::cout << "bra_ket = " << bra_ket << " gq = " << gq << std::endl;
  // std::cout << "ampfac_tomout (no 2sqrt(2) factor yet) = " << ampfac_tomout << std::endl;

  //*2sqrt(2) as factor of sqrt(2) in epsilon numerator(not denominator!) and factor of 2 from Fierz identity
  ampfac_tomout *= 2*std::sqrt(2.);
  return ampfac_tomout;
}
