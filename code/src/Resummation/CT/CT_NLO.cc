#include "CT_NLO.h"

#include <vector>
#include <cmath>

#include "constants.h"
#include "AP_LO.h"
#include "resu_PS.h"
#include "resummation_input.h"
#include "Itilde.h"
#include "xsection.h"


double CT_NLO(const resu_PS& resuvars, ResummationInfo* info, const double alpha[]){

  //Set momentum fractions and introduce Monte Carlo extra variables alpha[0], alpha[1]
  double ax, ax1, ax2, x1, x2;
  ax = std::log(resuvars.x);
  ax1 = (ax+2*resuvars.eta)/2.0;
  ax2 = (ax-2*resuvars.eta)/2.0;
  x1 = std::exp(ax1);
  x2 = std::exp(ax2);
  double z1 = std::pow(x1, alpha[0]);
  double z2 = std::pow(x2, alpha[1]);

  int Nf = info->Nf;
  double mures2 = resuvars.mures2;
  double q2 = resuvars.q2;
  double LQ = std::log(q2/mures2);

  // PDFs -- MSTW08 for now
  int PDFlen = 2*Nf+1;
  double muf = std::sqrt(resuvars.muf2);
  double pdf_1[PDFlen], pdf_2[PDFlen], pdf0_1[PDFlen], pdf0_2[PDFlen];
  //proton or antiproton
  int ih1 = info->ih1;
  int ih2 = info->ih2;
  //momentum fractions pre z1, z2 splittings
  double x1p = x1/z1, x2p = x2/z2;
  k_getPDF(info->pdf_res, pdf0_1, Nf, ih1, x1, muf);
  k_getPDF(info->pdf_res, pdf0_2, Nf, ih2, x2, muf);
  k_getPDF(info->pdf_res, pdf_1, Nf, ih1, x1p, muf);
  k_getPDF(info->pdf_res, pdf_2, Nf, ih2, x2p, muf);
  //Note PDF labelling is 0 -> 10 with 5(=Nf) as the gluon; 0, 1, 2, 3, 4 are antiquarks and 6, 7, 8, 9, 10 are quarks

  // Born x-sec
  std::vector<std::vector<double> > sigmaBorn = resuvars.sigmaij;

  // std::cout << "sigmaBorn = " << std::endl;
  //   for(int j=0; j<=2*Nf; j++){
  //     for(int k=0; k<=2*Nf; k++){
  // 	std::cout << "sigmaBorn[" << j << "][" << k << "]=" << sigmaBorn[j][k] << std::endl;
  //     }
  //   }

  //Calculate the Sigmas up to NLO for now
  double Sigma11 = 0., Sigma12 = 0., CT_out=0.;

  if( info->order != 1 ){
    std::cout << "Only NLO CT implemented so far" << std::endl;
    exit(EXIT_FAILURE);
  }

  if(info->qq_order == 1){
// qqbar-initiated
    // std::cout << "Doing qq initiated bit...." << std::endl;
    double A1q = info->resuNindep.A1q;
    double B1q = info->resuNindep.B1q;

    for(int j=0; j<=2*Nf; j++){
      for(int k=0; k<=2*Nf; k++){
        if(j!=Nf && k!=Nf){ //i.e. summing over all quarks and antiquarks but not gluons

//Sigma12 first as no Mellin dependence (as no splitting functions)
          Sigma12 = Sigma12 + pdf0_1[j]*pdf0_2[k]*-0.5*A1q*sigmaBorn[j][k];

//Sigma11 is more difficult because as well as the constant term it also has Mellin space dependence via the splitting functions
          Sigma11 = Sigma11 + (-(B1q+A1q*LQ)*pdf0_1[j]*pdf0_2[k] + std::log(x1)*Pqg(z1)*pdf0_2[k]*pdf_1[Nf] + std::log(x2)*Pqg(z2)*pdf0_1[j]*pdf_2[Nf] + pdf0_1[j]*pdf0_2[k]*Pqqint(x1) + pdf0_1[j]*pdf0_2[k]*Pqqint(x2) + pdf0_1[j]*std::log(x2)*Pqqreg(z2)*(pdf_2[k]-z2*pdf0_2[k]) + pdf0_2[k]*std::log(x1)*Pqqreg(z1)*(pdf_1[j]-z1*pdf0_1[j]))*sigmaBorn[j][k];
	  // std::cout << "B1q = " << B1q << " A1q = " << A1q << " LQ = " << LQ << " pdf0_1[j] = " << pdf0_1[j] << " pdf0_2[k] = " << pdf0_2[k] << " x1 = " << x1 << " z1 = " << z1 << " Pqg(z1) = " << Pqg(z1) << std::endl;
	  // std::cout << "sigmaBorn[j][k] = " << sigmaBorn[j][k] << std::endl;
        }
      }
    }
  }
  if(info->gg_order == 1) {
// gg-initiated
    // std::cout << "Doing gg initiated bit...." << std::endl;

    double A1g = info->resuNindep.A1g;
    double B1g = info->resuNindep.B1g;
    double beta0 = info->resuNindep.beta0;

    //Note only summing needed now is when considering Pgq splittings, otherwise only gluon contributions so no sum

    //Sigma12 first as easier as no Mellin space dependence as now splittin functions
    Sigma12 = Sigma12 + pdf0_1[Nf]*pdf0_2[Nf]*-0.5*A1g*sigmaBorn[Nf][Nf];

    //Now Sigma11 -> harder as splitting functions and also as need to sum over quarks and antiquarks for Pgq case
    //First the constant contribution and the Pgg contributions so no sum required
    Sigma11 = Sigma11 + (pdf0_1[Nf]*pdf0_2[Nf]*-(B1g+A1g*LQ)-beta0/2*pdf0_1[Nf]*pdf0_2[Nf]+6*D0int(x2)*pdf0_1[Nf]*pdf0_2[Nf] + 6*D0int(x1)*pdf0_1[Nf]*pdf0_2[Nf] + pdf0_1[Nf]*std::log(x2)*Pggreg(z2)*pdf_2[Nf] + pdf0_2[Nf]*std::log(x1)*Pggreg(z1)*pdf_1[Nf] + pdf0_1[Nf]*std::log(x2)*6/(1-z2)*(pdf_2[Nf]-z2*pdf0_1[Nf]) + pdf0_2[Nf]*std::log(x1)*6/(1-z1)*(pdf_1[Nf]-z1*pdf0_1[Nf]))*sigmaBorn[Nf][Nf];

    //Now the Pgq, therefore we must sum of the quarks and antiquarks - but only one sum as only one index - the g is still fixed as the outcome as it initiates
    for(int j=0; j<=2*Nf; j++){
      if(j!=Nf){ //As looking at q,qbar -> g splitting so j must not be a gluon (that is Pgg and is above)
	Sigma11 = Sigma11 + (pdf0_1[Nf]*pdf_2[j]*std::log(x2)*Pgq(z2) + pdf0_2[Nf]*pdf_1[j]*std::log(x1)*Pgq(z1))*sigmaBorn[j][Nf];
      }
    }
  }
  //Now in order to obtain the counter term we also need the Itilde integrals which do the b space integrals involving the powers of logs
  // qT Logs
  double Itarg = std::sqrt(resuvars.qt2/q2);

  //double CT_1 = Sigma11*Itilde(1,Itarg)/(q2*q2) + Sigma12*Itilde(2,Itarg)/(q2*q2);
  //The Itilde (at least in the qT/Q > 1.5 region where they are done numerically, are identical to the Itilde in 0508068, therefore no 1/(q2*q2) are needed)
  double CT_1 = Sigma11*Itilde(1,Itarg) + Sigma12*Itilde(2,Itarg);
  std::cout << "CT_1 = " << CT_1 << std::endl;
  std::cout << "Sigma11*Itilde(1,Itarg) = " << Sigma11*Itilde(1,Itarg) << " Sigma12*Itilde(2,Itarg) = " << Sigma12*Itilde(2,Itarg) << std::endl;
  std::cout << "Sigma11 = " << Sigma11 << " Itilde(1,Itarg) = " << Itilde(1,Itarg) << " Sigma12 = " << Sigma12 << " Itilde(2,Itarg) = " << Itilde(2,Itarg) << std::endl;
  // std::cout << "CT_1*q^4 = " << CT_1*q2*q2 << std::endl;
  // std::cout << "Check the Itilde agree at 1.5- and 1.5+ given they have two different forms near there...." << std::endl;
  // std::cout << "Itilde(1,1.49999999) = " << Itilde(1,1.49999999) << std::endl;
  // std::cout << "Itilde(1,1.50000001) = " << Itilde(1,1.50000001) << std::endl;
  // std::cout << "Itilde(2,1.49999999) = " << Itilde(2,1.49999999) << std::endl;
  // std::cout << "Itilde(2,1.50000001) = " << Itilde(2,1.50000001) << std::endl;
  // std::cout << "Itilde(1,0.001) = " << Itilde(1,0.001) << std::endl;
  // std::cout << "Itilde(1,0.01) = " << Itilde(1,0.01) << std::endl;
  // std::cout << "Itilde(1,0.1) = " << Itilde(1,0.1) << std::endl;
  // std::cout << "Itilde(1,1) = " << Itilde(1,1) << std::endl;
  // std::cout << "Itilde(1,10) = " << Itilde(1,10) << std::endl;
  // std::cout << "Itilde(1,0.01) = " << Itilde(1,0.01) << std::endl;
  // std::cout << "Itilde(1,0.01) = " << Itilde(1,0.01) << std::endl;
  // std::cout << "Itilde(1,0.01) = " << Itilde(1,0.01) << std::endl;
  // std::cout << "Itilde(1,0.01) = " << Itilde(1,0.01) << std::endl;
  // std::cout << "Itilde(1,0.01) = " << Itilde(1,0.01) << std::endl;
  // std::cout << "Itilde(1,0.01) = " << Itilde(1,0.01) << std::endl;

      
    //TCRIDGE I THINK THAT WHEN WE DO THE CT WE EXCHANGE THE Q2 AND Y INTEGRALS FOR Z1 (->ALPHA1) and Z2 (->ALPHA2) INTEGRALS, THEREFORE PERHAPS RANDSJACOB FOR THE CT SHOULD NOT THEN HAVE THE PIECES FOR THE Q2 AND Y? 
  //TCRIDGE ALTERED - I DON'T THINK YOU ADD THE BORN PIECE UNTIL THE NLO COUNTER TERM IS NEEDED, WE CURRENTLY CALCULATE THE LO COUNTER TERM AS LO IN BOSON + JET
  // CT_out = resuvars.alphas*CT_1*resuvars.x; //Note only resuvars.alphas needed here as in resu_PS.cc alphas = alphas(mur)/pi and I used splitting functions with alphas/pi normalisation so ok, the x factor comes from the x/q2 flux factor but the q2 cancels with that from the jacobian of dq^2dy ->dz1dz2, don't think q2*q2 necessary here unless to cancel some factors in the Itilde, note got rid of *q2*q2 here as dropped 1/(q2*q2) in CT_1
  CT_out = resuvars.alphas*CT_1*resuvars.x/q2;
  std::cout << "qt = " << std::sqrt(resuvars.qt2) << " x = " << resuvars.x << " q2 = " << resuvars.q2 << " CT_out = " << CT_out << " z = " << z1*z2 << std::endl;
  // CT_out = resuvars.alphas*CT_1*resuvars.x*resuvars.x/q2; //try *x/q2 extra as seems to be about correct factor
  // std::cout << "qt = " << std::sqrt(resuvars.qt2) << " x = " << resuvars.x << " q2 = " << resuvars.q2 << " CT_out = " << CT_out << std::endl;

  
  //Add the Born cross-section (the LO piece) to obtain the NLO counter term
  // CT_out = xsection(resuvars, info) + 4*resuvars.alphas*CT_1*q2*q2;
  // CT_out = xsection(resuvars, info) + 4*resuvars.alphas*CT_1*q2*q2;
  //CT_out = xsection(resuvars, info)*resuvars.x/resuvars.q2 + 4*resuvars.alphas*CT_1; //factor 4 here accounts for fact splittin functions are alphas/pi normalisation but that we use alphas from resu_procindep.cc which has a 1/4pi normalisation
  //std::cout << "qt = " << std::sqrt(resuvars.qt2) << " xsec = " << xsection(resuvars,info) << " xsec*x/q2 = " << xsection(resuvars,info)*resuvars.x/resuvars.q2 << " CT = " << 4*resuvars.alphas*CT_1 << " x = " << resuvars.x << " q2 = " << resuvars.q2 << " CT_out = " << CT_out << std::endl;
  // std::cout << "qt = " << std::sqrt(resuvars.qt2) << " xsec = " << xsection(resuvars,info) << " CT = " << 4*resuvars.alphas*CT_1 << " CT*q^4 = " << 4*resuvars.alphas*CT_1*q2*q2 << std::endl;
  // std::cout << "qt = " << std::sqrt(resuvars.qt2) << " xsec = " << xsection(resuvars, info) << " CTtemp = CT_1*4as*q^4 = " << CT_1*4*resuvars.alphas*q2*q2 << " xsec/CT_temp = " << xsection(resuvars, info)/(CT_1*4*resuvars.alphas*q2*q2) << std::endl;  
  // std::cout << "qt = " << std::sqrt(resuvars.qt2) << " q = " << std::sqrt(q2) << " xsec = " << xsection(resuvars, info) << " CTtemp = CT_1*4as*q^4 = " << CT_1*4*resuvars.alphas*q2*q2 << " xsec/CT_temp = " << xsection(resuvars, info)/(CT_1*4*resuvars.alphas*q2*q2) <<  " xsec/(CT_temp/q) = " << xsection(resuvars, info)/(CT_1*4*resuvars.alphas*q2*std::sqrt(q2)) << std::endl;
  // CT_out = xsection(resuvars, info)*(1 + 4*resuvars.alphas*CT_1);
  // std::cout << "qt = " << std::sqrt(resuvars.qt2) << " xsec = " << xsection(resuvars,info) << " CT = " << xsection(resuvars,info)*4*resuvars.alphas*CT_1 << std::endl;
  
  
  return CT_out;
}
