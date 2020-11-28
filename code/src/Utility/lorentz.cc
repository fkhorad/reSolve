
#include "lorentz.h"

// LORENTZ SPACE RELATED FUNCTIONS

// While there several libraries implementing Lorentz vector and the like around,
// the idea here is to keep everything as simple as possible: I would prefer to avoid
// linking any "big" libraries unless it is strictly necessary.

// IN PARTICULAR: I use std::vector<double> for 4-momenta, aliased as "four_momentum" to allow easier replacing in the future. There is no "4-momentum class"

double LorDot(const four_momentum& mom1, const four_momentum& mom2){

    double res;

    res = mom1[0]*mom2[0]-mom1[1]*mom2[1]-mom1[2]*mom2[2]-mom1[3]*mom2[3];

    return res;

};
four_momentum LorDot(const std::vector<four_momentum>& MM, const four_momentum& vv){

  four_momentum res(4, 0.);
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      res[i] += MM[i][j]*vv[j];
    }
  }
  return res;
};
std::vector<four_momentum> LorDot(const std::vector<four_momentum>& M1, const std::vector<four_momentum>& M2){
  std::vector<four_momentum> res(4, four_momentum(4, 0.) );
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      res[i][j] = 0.;
      for(int k=0; k<4; k++){
        res[i][j] += M1[i][k]*M2[k][j];
      }
    }
  }
  return res;

};

double LorNorm(const four_momentum& mom){

  double res;

  res = mom[0]*mom[0]-mom[1]*mom[1]-mom[2]*mom[2]-mom[3]*mom[3];

  return res;

};


// Using Mangano-Parke conventions
std::complex<double> SpinProdA(const four_momentum& mom1, const four_momentum& mom2){

  std::complex<double> res;
  std::complex<double> II = std::complex<double>(0.,1.);

  double mom1P = mom1[0] + mom1[3];
  double mom2P = mom2[0] + mom2[3];
  // std::cout << "mom1 = " << mom1[0] << " " << mom1[1] << " " << mom1[2] << " " << mom1[3] << std::endl;
  // std::cout << "mom2 = " << mom2[0] << " " << mom2[1] << " " << mom2[2] << " " << mom2[3] << std::endl;
  // std::cout << "mom1P = " << mom1P << std::endl;
  // std::cout << "mom2P = " << mom2P << std::endl;

  double prod = mom1P*mom2P;
  // std::cout << "prod  = " << prod  << std::endl;

  if( prod == 0.){
    res = std::complex<double>(std::sqrt(2.*std::abs(LorDot(mom1,mom2))), 0.);
    if(mom1P<0.) res*=-1.;
    if(mom2P>0.) res*=-1.;
  }
  else{
    double fac = 1./std::sqrt(std::abs(prod));
    // std::cout << "fac   = " << fac   << std::endl;
    double res_R = fac*(mom1[1]*mom2P - mom2[1]*mom1P);
    double res_I = fac*(mom1[2]*mom2P - mom2[2]*mom1P);
    // std::cout << "res_R  = " << res_R  << std::endl;
    // std::cout << "res_I  = " << res_I  << std::endl;
    res = std::complex<double>(res_R, res_I);
  }

  if(mom1[0]<0 and mom2[0]>0) res*=-II;
  if(mom1[0]>0 and mom2[0]<0) res*=-II;
  if(mom1[0]<0 and mom2[0]<0) res*=-1.;

  return res;
}

std::vector<four_momentum> SetBoost(const four_momentum& pp, bool invert){

  std::vector<four_momentum> lor1(4, four_momentum(4, 0.) );
  double beta[3];

  if (invert){
    beta[0]=-pp[1]/pp[0]; beta[1]=-pp[2]/pp[0]; beta[2]=-pp[3]/pp[0];
  }
  else{
    beta[0]=pp[1]/pp[0]; beta[1]=pp[2]/pp[0]; beta[2]=pp[3]/pp[0];
  }

  double beta2, denom, gamma;
  beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
  gamma = 1/std::sqrt(1-beta2);
  denom = beta2;
  if(denom==0.) denom = 1.;

  lor1[0][0]=gamma;
  lor1[0][1]=beta[0]*gamma;
  lor1[0][2]=beta[1]*gamma;
  lor1[0][3]=beta[2]*gamma;
  lor1[1][0]=beta[0]*gamma;
  lor1[1][1]=1.+(gamma-1.)*beta[0]*beta[0]/denom;
  lor1[1][2]=(gamma-1.)*beta[0]*beta[1]/denom;
  lor1[1][3]=(gamma-1.)*beta[0]*beta[2]/denom;
  lor1[2][0]=beta[1]*gamma;
  lor1[2][1]=(gamma-1.)*beta[1]*beta[0]/denom;
  lor1[2][2]=1.+(gamma-1.)*beta[1]*beta[1]/denom;
  lor1[2][3]=(gamma-1.)*beta[1]*beta[2]/denom;
  lor1[3][0]=beta[2]*gamma;
  lor1[3][1]=(gamma-1.)*beta[2]*beta[0]/denom;
  lor1[3][2]=(gamma-1.)*beta[2]*beta[1]/denom;
  lor1[3][3]=1.+(gamma-1.)*beta[2]*beta[2]/denom;

  return lor1;

};

// Useful for debugging
void LorPrint(const four_momentum& vec){
  std::cout << vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3] << "\n";
};
void LorPrint(const four_momentum& vec, std::ofstream& file){
  file << vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3] << "\n";
};
