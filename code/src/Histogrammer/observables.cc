
#include "observables.h"

#include "phase_space.h"
#include "lorentz.h"


int obs_values(std::string obs, const PSpoint& PS_, double& value){
  int return_code = 0;

  if(obs=="qT") {
    value = qT_obs(PS_);
    // std::cout << "making qT histogram" << std::endl;
  }
  else if(obs=="qq") {
    value = qq_obs(PS_);
     // std::cout << "making qq histogram" << std::endl;
  }
  else if (obs=="eta") {
    value = eta_obs(PS_);
     // std::cout << "making eta histogram" << std::endl;
  }
  else if (obs=="mT") {
    value = mT_obs(PS_);
    // std::cout << "making mT histogram" << std::endl;
  }
  else if (obs=="pTmin") {
    value = pTmin_obs(PS_);
     // std::cout << "making pTmin histogram" << std::endl;
  }
  else if (obs=="pTmax") {
    value = pTmax_obs(PS_);
     // std::cout << "making pTmax histogram" << std::endl;
  }
  else{
    value = 0.;
    return_code = 1;
  }

  return return_code;
}



double qT_obs(const PSpoint& PS_){

  double qTx = PS_.mom(2)[1] + PS_.mom(3)[1];
  double qTy = PS_.mom(2)[2] + PS_.mom(3)[2];
  double qT2 = qTx*qTx + qTy*qTy;

  return std::sqrt(qT2);

};

double qq_obs(const PSpoint& PS_){

  four_momentum q_vec(4);
  for(int ii=0; ii<4; ii++) q_vec[ii] = PS_.mom(2)[ii] + PS_.mom(3)[ii];
  double qq = LorNorm(q_vec);

  return std::sqrt(qq);

};

double eta_obs(const PSpoint& PS_){

  double E = PS_.mom(2)[0] + PS_.mom(3)[0];
  double pZ = PS_.mom(2)[3] + PS_.mom(3)[3];
  double eta = 0.5*log((E+pZ)/(E-pZ));
  // std::cout << "eta sys = " << eta << std::endl;
  // std::cout << "E = " << E << " pZ = " << pZ << std::endl;

  return eta;
};

double mT_obs(const PSpoint& PS_){

  double qT1sq = PS_.mom(2)[1]*PS_.mom(2)[1] + PS_.mom(2)[2]*PS_.mom(2)[2];
  double qT2sq = PS_.mom(3)[1]*PS_.mom(3)[1] + PS_.mom(3)[2]*PS_.mom(3)[2];
  double qT1dotqT2 = PS_.mom(2)[1]*PS_.mom(3)[1] + PS_.mom(2)[2]*PS_.mom(3)[2];
  double qT1 = std::sqrt(qT1sq);
  double qT2 = std::sqrt(qT2sq);

  double mT = std::sqrt(2.0*(qT1*qT2-qT1dotqT2));


  return mT;
};

double pTmin_obs(const PSpoint& PS_){

  double qT1sq = PS_.mom(2)[1]*PS_.mom(2)[1] + PS_.mom(2)[2]*PS_.mom(2)[2];
  double qT2sq = PS_.mom(3)[1]*PS_.mom(3)[1] + PS_.mom(3)[2]*PS_.mom(3)[2];
  double pTminsq = 0;
  if (qT1sq < qT2sq)
    pTminsq = qT1sq;
  else
    pTminsq = qT2sq;

  return std::sqrt(pTminsq);
};

double pTmax_obs(const PSpoint& PS_){

  double qT1sq = PS_.mom(2)[1]*PS_.mom(2)[1] + PS_.mom(2)[2]*PS_.mom(2)[2];
  double qT2sq = PS_.mom(3)[1]*PS_.mom(3)[1] + PS_.mom(3)[2]*PS_.mom(3)[2];
  double pTmaxsq = 0;
  if (qT1sq < qT2sq)
    pTmaxsq = qT2sq;
  else
    pTmaxsq = qT1sq;

  return std::sqrt(pTmaxsq);
}
