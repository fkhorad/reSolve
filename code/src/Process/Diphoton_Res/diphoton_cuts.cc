
#include "diphoton_cuts.h"
//Function to do all diphoton cuts via phase space and kinematical limit calls
bool diph_cuts(diphoton_input* diph_in, double q2, double qt2, double eta, double mures2, PSpoint* PS_) {
    bool gencuts = false; //i.e. not cut
    double xx = 0, ax = 0;
    xx = q2/std::pow(diph_in->CM_energy,2);
    ax = std::log(xx);
// Check that y does not exceed kinematical limits
    if (std::abs(eta)>(-0.5*ax)) {
	gencuts = true; //i.e. cut
    }
// Check that qT does not exceed kinematical limits
    if (qt2 > 4*mures2) {
	gencuts = true; //i.e. cut
    }

//Now phase space cuts
    bool pscuts = false; //i.e. not cut
    pscuts = PScuts_1(diph_in, PS_);

// END CUTS BLOCK

    bool cuts = false;
    if (gencuts == true || pscuts == true) {
	cuts = true; //i.e. fail cuts if either cut on kinematical limits (gencuts) or on phase space (pscuts)
    }
    else {
	cuts = false;
    }

    return cuts;
}


// Standard set of phase space diphoton cuts
bool PScuts_1(diphoton_input* diph_1, PSpoint* PS_){

  four_momentum MOM3 = PS_->mom(2);
  four_momentum MOM4 = PS_->mom(3);


// pT1 and pT2 cuts
  double pT1 = std::sqrt(std::pow(MOM3[1],2) + std::pow(MOM3[2],2));
  double pT2 = std::sqrt(std::pow(MOM4[1],2) + std::pow(MOM4[2],2));


  double pT1cut = diph_1->pT1cut;
  double pT2cut = diph_1->pT2cut;
  if(std::max(pT1,pT2) < pT1cut) {
      return true;
  }
  if(std::min(pT1,pT2) < pT2cut) { 
      return true;
  }
  
// eta_1 and eta_2 cuts
  double max_eta = diph_1->etaCut;
  double crack1 = diph_1->crack1; double crack2 = diph_1->crack2;
  double eta_1 = 1./2. * std::log((MOM3[0] + MOM3[3])/(MOM3[0] - MOM3[3]));
  double aeta1 = std::abs(eta_1);
  double eta_2 = 1./2. * std::log((MOM4[0] + MOM4[3])/(MOM4[0] - MOM4[3]));
  double aeta2 = std::abs(eta_2);
  if (diph_1->verbosity >= 12) {
      std::cout << "aeta1 = " << aeta1 << std::endl;
      std::cout << "aeta2 = " << aeta2 << std::endl;
  }
  

  if(aeta1 > max_eta or aeta2 > max_eta) {
      return true;
  }
  if(crack1 < aeta1 and aeta1 < crack2 ) {
      return true;
  }
  if(crack1 < aeta2 and aeta2 < crack2 ) {
      return true;
  }
 
// Delta_R
  double phi_1 = std::atan2(MOM3[1], MOM3[2]);
  double phi_2 = std::atan2(MOM4[1], MOM4[2]);
  double Delta_R = std::sqrt(std::pow(eta_1-eta_2,2) + std::pow(phi_1-phi_2,2));
  double Rcut = diph_1->Rcut;
 
  if(Delta_R < Rcut) {
      return true;
  }
  return false;
};
