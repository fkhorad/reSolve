
#include "drellyan_res_hard.h"


#include "constants.h"
#include "phase_space.h"
#include "drellyan_input.h"
#include "diphoton_hard.h" //for Asmitad etc



void sigmaijdrellyancalc(drellyan_input* drellyan_in, PSpoint& PS, std::vector<std::vector<double> >& sigmaij, double alphas, double& q2) {
// Born partonic cross-sections
// Normalisation is such that the COMPLETE dsigma/dOmega double differential cross-sections are returned
//
// Constants
    double gevpb = drellyan_in->res_1.gevm2topb;
    double pi = k_constants::pi;
// Overall factor for amplitudes
    double flux = 1./2./q2;
    double spincolavg = 1/12.;
    double PSfac = 1/32./pi/pi; // == 1/16./pi/pi * |p_CM|/|E^tot_CM| in the massles limit
    double totfac = flux*spincolavg*PSfac*gevpb;

    //TCRIDGE FOR TESTING OCT20 VS MCFM LO
    //q2 = 8204.6921208607964;

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
// Kinematics
    double s = 2.*PS.ss(1,2);
    double t = -2.*PS.ss(1,3);
    double costheta = 1. + 2.*t/s; // Note that theta here is the CM angle between parton 1 and the negatively charged lepton -- regardless of whether parton 1 is a quark or a antiquark

    std::cout << "s = " << s << " t = " << t << " costheta = " << costheta << std::endl;
    std::cout << "q2 = " << q2 << std::endl;
//
// Propagator factors for W, Z and Z'
    double facW = q2*q2/( (q2-mw*mw)*(q2-mw*mw)+mw*mw*ww*ww ); // is W on-shell only??
    double facZ = q2*q2/( (q2-mz*mz)*(q2-mz*mz)+mz*mz*zw*zw );
    double facZ1 = ( q2*(q2-mz*mz) ) / ( (q2-mz*mz)*(q2-mz*mz)+mz*mz*zw*zw );
    double facZ2 = ( q2*mz*zw ) / ( (q2-mz*mz)*(q2-mz*mz)+mz*mz*zw*zw );
//
    double facZp1 = ( q2*(q2-mzp*mzp) ) / ( (q2-mzp*mzp)*(q2-mzp*mzp)+mzp*mzp*zpw*zpw );
    double facZp2 = ( q2*mzp*zpw ) / ( (q2-mzp*mzp)*(q2-mzp*mzp)+mzp*mzp*zpw*zpw );

////////////////////////


// Charged current: subprocesses 1-3
//
    if(DYproc == 1 || DYproc == 2 || DYproc == 3) {
      if(DYproc == 1 || DYproc ==3){
// W+ boson contribution (for processes 1 and 3)
        sigmaij[Nf+1][Nf-2] = totfac*facW* Vud*Vud*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //u dbar
        sigmaij[Nf-2][Nf+1] = totfac*facW* Vud*Vud*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //dbar u
        if(Nf>=3) sigmaij[Nf+1][Nf-3] = totfac*facW* Vus*Vus*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //u sbar
        if(Nf>=3) sigmaij[Nf-3][Nf+1] = totfac*facW* Vus*Vus*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //sbar u
        if(Nf>=5) sigmaij[Nf+1][Nf-5] = totfac*facW* Vub*Vub*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //u bbar
        if(Nf>=5) sigmaij[Nf+1][Nf-5] = totfac*facW* Vub*Vub*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //bbar u
        if(Nf>=4) sigmaij[Nf+4][Nf-3] = totfac*facW* Vcs*Vcs*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //c sbar
        if(Nf>=4) sigmaij[Nf-3][Nf+4] = totfac*facW* Vcs*Vcs*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //sbar c
        if(Nf>=4) sigmaij[Nf+4][Nf-2] = totfac*facW* Vcd*Vcd*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //c dbar
        if(Nf>=4) sigmaij[Nf+4][Nf-2] = totfac*facW* Vcd*Vcd*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //dbar c
        if(Nf>=5) sigmaij[Nf+4][Nf-5] = totfac*facW* Vcb*Vcb*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //c bbar
        if(Nf>=5) sigmaij[Nf-5][Nf+4] = totfac*facW* Vcb*Vcb*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //bbar c
      }
      if(DYproc == 2 || DYproc == 3) {
// W- boson contribution (for subprocesses 2 and 3)
        sigmaij[Nf+2][Nf-1] = totfac*facW* Vud*Vud*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //d ubar
        sigmaij[Nf-1][Nf+2] = totfac*facW* Vud*Vud*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //ubar d
        if(Nf>=3) sigmaij[Nf+3][Nf-1] = totfac*facW* Vus*Vus*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //s ubar
        if(Nf>=3) sigmaij[Nf-1][Nf+3] = totfac*facW* Vus*Vus*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //ubar s
        if(Nf>=5) sigmaij[Nf+5][Nf-1] = totfac*facW* Vub*Vub*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //b ubar
        if(Nf>=5) sigmaij[Nf-1][Nf+5] = totfac*facW* Vub*Vub*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //ubar b
        if(Nf>=4) sigmaij[Nf+3][Nf-4] = totfac*facW* Vcs*Vcs*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //s cbar
        if(Nf>=4) sigmaij[Nf-4][Nf+3] = totfac*facW* Vcs*Vcs*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //cbar s
        if(Nf>=4) sigmaij[Nf+2][Nf-4] = totfac*facW* Vcd*Vcd*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //d cbar
        if(Nf>=4) sigmaij[Nf-4][Nf+2] = totfac*facW* Vcd*Vcd*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //cbar d
        if(Nf>=5) sigmaij[Nf+5][Nf-4] = totfac*facW* Vcb*Vcb*cLw*cLw*cLw*cLw*(1.0+costheta)*(1.0+costheta); //b cbar
        if(Nf>=5) sigmaij[Nf-4][Nf+5] = totfac*facW* Vcb*Vcb*cLw*cLw*cLw*cLw*(1.0-costheta)*(1.0-costheta); //cbar b
      }
    }
//
// Neutral current: subprocesses 4-6
//
// Z boson only (no photon): subprocess 4
    else if (DYproc == 4) {
      double FLL[Nf], FLR[Nf], FRL[Nf], FRR[Nf];
      for(unsigned int ii=0; ii<Nf; ii++){
        FLL[ii] = std::pow( gLzQ[ii]*cLzl, 2)*facZ;
        FLR[ii] = std::pow( gLzQ[ii]*cRzl, 2)*facZ;
        FRL[ii] = std::pow( gRzQ[ii]*cLzl, 2)*facZ;
        FRR[ii] = std::pow( gRzQ[ii]*cRzl, 2)*facZ;
//
        sigmaij[Nf+(ii+1)][Nf-(ii+1)] = totfac*( (FLL[ii]+FRR[ii])*(1+costheta)*(1+costheta) + (FLR[ii]+FRL[ii])*(1-costheta)*(1-costheta) );
//
        sigmaij[Nf-(ii+1)][Nf+(ii+1)] = totfac*( (FLL[ii]+FRR[ii])*(1-costheta)*(1-costheta) + (FLR[ii]+FRL[ii])*(1+costheta)*(1+costheta) );
      }
    }
// gamma* + Z -- subprocess 5
    else if (DYproc == 5) {
      double FLL[Nf], FLR[Nf], FRL[Nf], FRR[Nf];
      for(unsigned int ii=0; ii<Nf; ii++){
        FLL[ii] = std::pow( -QQ[ii]*ee*ee + gLzQ[ii]*cLzl*facZ1, 2) + std::pow( gLzQ[ii]*cLzl*facZ2, 2);
        FLR[ii] = std::pow( -QQ[ii]*ee*ee + gLzQ[ii]*cRzl*facZ1, 2) + std::pow( gLzQ[ii]*cRzl*facZ2, 2);
        FRL[ii] = std::pow( -QQ[ii]*ee*ee + gRzQ[ii]*cLzl*facZ1, 2) + std::pow( gRzQ[ii]*cLzl*facZ2, 2);
        FRR[ii] = std::pow( -QQ[ii]*ee*ee + gRzQ[ii]*cRzl*facZ1, 2) + std::pow( gRzQ[ii]*cRzl*facZ2, 2);
//
        sigmaij[Nf+(ii+1)][Nf-(ii+1)] = totfac*( (FLL[ii]+FRR[ii])*(1+costheta)*(1+costheta) + (FLR[ii]+FRL[ii])*(1-costheta)*(1-costheta) );
//
        sigmaij[Nf-(ii+1)][Nf+(ii+1)] = totfac*( (FLL[ii]+FRR[ii])*(1-costheta)*(1-costheta) + (FLR[ii]+FRL[ii])*(1+costheta)*(1+costheta) );
      }
    }
//
// gamma* + Z + Z' -- subprocess 6
    else if (DYproc == 6) {
      double FLL[Nf], FLR[Nf], FRL[Nf], FRR[Nf];
      for(unsigned int ii=0; ii<Nf; ii++){
        FLL[ii] = std::pow( -QQ[ii]*ee*ee + gLzQ[ii]*cLzl*facZ1 + gLZpQ[ii]*cLZpl*facZp1, 2) + std::pow( gLzQ[ii]*cLzl*facZ2 + gLZpQ[ii]*cLZpl*facZp2, 2);
        FLR[ii] = std::pow( -QQ[ii]*ee*ee + gLzQ[ii]*cRzl*facZ1 + gLZpQ[ii]*cRZpl*facZp1, 2) + std::pow( gLzQ[ii]*cRzl*facZ2 + gLZpQ[ii]*cRZpl*facZp2, 2);
        FRL[ii] = std::pow( -QQ[ii]*ee*ee + gRzQ[ii]*cLzl*facZ1 + gRZpQ[ii]*cLZpl*facZp1, 2) + std::pow( gRzQ[ii]*cLzl*facZ2 + gRZpQ[ii]*cLZpl*facZp2, 2);
        FRR[ii] = std::pow( -QQ[ii]*ee*ee + gRzQ[ii]*cRzl*facZ1 + gRZpQ[ii]*cRZpl*facZp1, 2) + std::pow( gRzQ[ii]*cRzl*facZ2 + gRZpQ[ii]*cRZpl*facZp2, 2);
//
        sigmaij[Nf+(ii+1)][Nf-(ii+1)] = totfac*( (FLL[ii]+FRR[ii])*(1+costheta)*(1+costheta) + (FLR[ii]+FRL[ii])*(1-costheta)*(1-costheta) );
//
        sigmaij[Nf-(ii+1)][Nf+(ii+1)] = totfac*( (FLL[ii]+FRR[ii])*(1-costheta)*(1-costheta) + (FLR[ii]+FRL[ii])*(1+costheta)*(1+costheta) );
      }
    }
    else{
      std::cout << "Wrong DYprocess flag -- STOPPING" << std::endl;
      exit(EXIT_FAILURE);
    }

    // //For now set sigmaij to just include the PS fac * the spincolavg as want to test PS vol
    for(int ii=0; ii<=2*Nf; ii++){
      for(int ij=0; ij<=2*Nf; ij++){
    // 	// std::cout << "Nf = " << Nf << std::endl;
    // 	// std::cout << ii << " " << ij << std::endl;
    // 	sigmaij[ii][ij] = PSfac;
    // 	// std::cout << "PSfac*spincolavg = " << PSfac*spincolavg << std::endl;
	if (sigmaij[ii][ij]!=0) {
	  std::cout << "sigmaij[" << ii << "][" << ij << "]=" << sigmaij[ii][ij] << std::endl;
	  std::cout << "sigmaij[" << ii << "][" << ij << "]/(flux*gevpb)*pi=" << sigmaij[ii][ij]/(flux*gevpb)*pi << std::endl;
	  
	}
      }
    }

}
