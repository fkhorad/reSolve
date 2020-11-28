
#include "drellyan_cuts.h"

#include "lorentz.h"
#include "phase_space.h"


//Function to do all drellyan cuts via phase space and kinematical limit calls
bool drellyan_cuts(drellyan_input* drellyan_in, double q2, double qt2, double eta, double mures2, PSpoint* PS_) {
  // std::cout << "q2 at start of cuts = " << q2 << std::endl;
    bool gencuts = false; //i.e. not cut
    double xx = 0, ax = 0;
    xx = q2/std::pow(drellyan_in->res_1.CM_energy,2);
    ax = std::log(xx);

// Check that y does not exceed kinematical limits
    if (std::abs(eta)>(-0.5*ax)) {
      gencuts = true; //i.e. cut
    }
// Check that qT does not exceed kinematical limits
//    if (qt2 > 4*mures2) {
//      gencuts = true; //i.e. cut
//    }

//Now phase space cuts
    bool pscuts = false; //i.e. not cut

    pscuts = PSDYcuts_1(drellyan_in, PS_);

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


// Standard set of phase space drellyan cuts
bool PSDYcuts_1(drellyan_input* drellyan_1, PSpoint* PS_){

  four_momentum MOM3 = PS_->mom(2);
  four_momentum MOM4 = PS_->mom(3);

  int DYprocess = drellyan_1->DYprocess;

// pT1 and pT2 cuts
  double pT1 = std::sqrt(std::pow(MOM3[1],2) + std::pow(MOM3[2],2));
  double pT2 = std::sqrt(std::pow(MOM4[1],2) + std::pow(MOM4[2],2));

//rapidities and pseudorapidities
  double y_1 = 1./2. * std::log((MOM3[0] + MOM3[3])/(MOM3[0] - MOM3[3])); //actual rapidity
  double ay1 = std::abs(y_1);
  double y_2 = 1./2. * std::log((MOM4[0] + MOM4[3])/(MOM4[0] - MOM4[3])); //actual rapidity
  double ay2 = std::abs(y_2);
//rapidities
  double eta_1 = 1./2. * std::log((MOM3[3] + std::sqrt(std::pow(MOM3[3],2)+std::pow(MOM3[2],2)+std::pow(MOM3[1],2)))/(std::sqrt(std::pow(MOM3[3],2)+std::pow(MOM3[2],2)+std::pow(MOM3[1],2))-MOM3[3])); //actual pseudorapidity
  double aeta1 = std::abs(eta_1);
  double eta_2 = 1./2. * std::log((MOM4[3] + std::sqrt(std::pow(MOM4[3],2)+std::pow(MOM4[2],2)+std::pow(MOM4[1],2)))/(std::sqrt(std::pow(MOM4[3],2)+std::pow(MOM4[2],2)+std::pow(MOM4[1],2))-MOM4[3]));  //actual pseudorapidity
  double aeta2 = std::abs(eta_2);

  double totpt = std::sqrt(std::pow(MOM3[1]+MOM4[1],2)+std::pow(MOM3[2]+MOM4[2],2)); //transverse mass in m->0 limit as is just (pxtot)^2 + (pytot)^2
  double totm = std::sqrt(std::pow(MOM3[0]+MOM4[0],2)-std::pow(MOM3[1]+MOM4[1],2)-std::pow(MOM3[2]+MOM4[2],2)-std::pow(MOM3[3]+MOM4[3],2));
  double toty = 0.5*std::log((MOM3[0]+MOM4[0]+MOM3[3]+MOM4[3])/(MOM3[0]+MOM4[0]-MOM3[3]-MOM4[3])); //actual total rapidity

  double pT1dotpT2 = MOM3[1]*MOM4[1]+MOM3[2]*MOM4[2]; //dot product of 2 transverse momenta

  double tmass2 = 2.0*(pT1*pT2-pT1dotpT2); //experimental definition of transverse mass as is just (ET1 + ET2)^2 - (pT1+pT2)^2, which in m=0 limit is 2(ET1*ET2 - pT1dotpT2) = 2(|pT1|*|pT2| - pT1dotpT2)
  double tmass = 0.0;

  if (tmass2 <0) {
    std::cout << "WARNING! tmass^2 = " << tmass2 << " is <0!" << std::endl;
  }
  else { tmass = std::sqrt(tmass2);}

  double pte = 0.0, ptmiss = 0.0, etae = 0.0;

  double crack1 = drellyan_1->crack1; double crack2 = drellyan_1->crack2;
  if(crack1 < aeta1 and aeta1 < crack2 ) {
    return true;
  }
  if(crack1 < aeta2 and aeta2 < crack2 ) {
    return true;
  }


// Z - Zp cuts
  if (DYprocess == 4 || DYprocess == 5 || DYprocess == 6) {
    if (pT1 < drellyan_1->pT1cut) {
      return true;
    }
    if (pT2 < drellyan_1->pT2cut) {
      return true;
    }
    if (aeta1 > drellyan_1->eta1cut) {
      return true;
    }
    if (aeta2 > drellyan_1->eta2cut) {
      return true;
    }
    if (std::max(aeta1,aeta2) > std::max(drellyan_1->etacutNEW1, drellyan_1->etacutNEW2) ) {
      return true;
    }
    if (std::min(aeta1,aeta2) > std::min(drellyan_1->etacutNEW1, drellyan_1->etacutNEW2) ) {
      return true;
    }

  }
// W cuts
  else if (DYprocess == 1 || DYprocess == 2 || DYprocess == 3) {
      pte = pT1;
      etae = eta_1;
      ptmiss = pT2;
      if (pte < drellyan_1->pTecut) {
        return true;
      }
      if (ptmiss < drellyan_1->pTmisscut) {
        return true;
      }
      if (std::abs(etae) > drellyan_1->etaecut) {
        return true;
      }
      if (tmass < drellyan_1->tmasscut) {
        return true;
      }
  }

// Forward, Backward
  double s = 2.*PS_->ss(1,2), t = -2.*PS_->ss(1,3);
  double costheta = 1. + 2.*t/s;
  double eta = 1./2.*(y_1 + y_2);
  if(drellyan_1->forward){
    if(costheta*eta < 0) return true;
  }
  if(drellyan_1->backward){
    if(costheta*eta > 0) return true;
  }

  return false;
}
