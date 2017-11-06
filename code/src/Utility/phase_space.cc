
#include "phase_space.h"

// PSpoint class methods

// Set:
int PSpoint::set_mom(int num, double* comps){

  int dimension = (int) momenta.size();

  if (num<dimension){
    for (int i=0;i<4;i++) momenta[num].at(i) = comps[i];
    return 0;
  }
  else
    return 1;
};
int PSpoint::set_mom(int num, four_momentum mom){

  int dimension = (int) momenta.size();

  if (num<dimension){
    momenta[num] = mom;
    return 0;
  }
  else
    return 1;
};

void PSpoint::set_products(){

  int dimension = momenta.size();

  std::vector<double> empty_r (dimension,0.);
  std::vector<std::complex<double> > empty_c (dimension,std::complex<double>(0.,0.));

  DotProducts.assign(dimension,empty_r);
  SpinProductsA.assign(dimension,empty_c);
  SpinProductsB.assign(dimension,empty_c);

  for(int i=0;i<dimension;i++){
    for(int j=0;j<dimension;j++){
      DotProducts[i][j] = LorDot(momenta[i],momenta[j]);
      SpinProductsA[i][j] = SpinProdA(momenta[i],momenta[j]);
      if(SpinProductsA[i][j] == std::complex<double>(0.,0.))
        SpinProductsB[i][j] = 0.;
      else
        SpinProductsB[i][j] = - DotProducts[i][j]/SpinProductsA[i][j];
    };
  };

};


void PSpoint::set_dim(int new_dim){

  four_momentum empty = four_momentum(4, 0.);
  momenta.resize(new_dim, empty);

};



// Get

int PSpoint::dim(){
   int res = momenta.size();
   return res;
};


four_momentum PSpoint::mom(int i){
  four_momentum res = momenta.at(i);
  return res;
};
double PSpoint::ss(int i, int j){
  double res = DotProducts.at(i).at(j);
  return res;
};
std::complex<double> PSpoint::sa(int i, int j){
  std::complex<double> res = SpinProductsA.at(i).at(j);
  return res;
};
std::complex<double> PSpoint::sb(int i, int j){
  std::complex<double> res = SpinProductsB.at(i).at(j);
  return res;
};


void PSpoint::PS_checker(int n_incoming){

  int DIM_ = momenta.size();

  std::cout << "PS CHECK" << "\n";
  for(int ii=0; ii<DIM_; ii++){
    std::cout << "Momentum " << ii+1 << " : ";
    LorPrint(momenta[ii]);
    std::cout << "Square of momentum " << ii+1 << " : " << LorNorm(momenta[ii]) << "\n";
  }

  four_momentum mom_temp = four_momentum(4, 0.);
  for(int ii=0; ii<n_incoming; ii++){
    for(int jj=0; jj<4; jj++){
      mom_temp[jj] += momenta[ii][jj];
    }
  }
  std::cout << "Sum of incoming momenta: ";
  LorPrint(mom_temp);
//
  mom_temp = four_momentum(4, 0.);
  for(int ii=n_incoming; ii<DIM_; ii++){
    for(int jj=0; jj<4; jj++){
      mom_temp[jj] += momenta[ii][jj];
    }
  }
  std::cout << "Sum of outgoing momenta: ";
  LorPrint(mom_temp);
  std::cout << "END OF PS CHECK" << "\n\n";

};


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

// Non-member functions: various basic PS-setting routines

double set_PS_twobody(const double* randoms, const double* masses, PSpoint& PS){

// Randoms to physical assignment
  double pi = k_constants::pi;
  double costheta = 2.*randoms[0]-1.;
  double phi = 2.*pi*randoms[2];
// Aux parameters
  double sintheta = std::sqrt(1.-pow(costheta,2));
//
  double lambda;
  lambda = std::pow(masses[0],4) + std::pow(masses[1],4) + std::pow(masses[2],4) -
    2.*std::pow(masses[0],2) * std::pow(masses[1],2) - 2.*std::pow(masses[0],2) * std::pow(masses[2],2) -
    2.*std::pow(masses[1],2) * std::pow(masses[2],2);
  double p12 = std::sqrt(lambda)/(2.*masses[0]);
  double vers[3] = {sintheta*std::cos(phi), std::sqrt(1.-std::pow(costheta,2))*std::sin(phi), costheta};
  double en1 = std::sqrt(std::pow(p12,2) + std::pow(masses[1],2));
  double en2 = std::sqrt(std::pow(p12,2) + std::pow(masses[2],2));
//
// The 2 4-momenta
  double mom_1[4] = {en1, p12*vers[0], p12*vers[1], p12*vers[2]};
  double mom_2[4] = {en2, -p12*vers[0], -p12*vers[1], -p12*vers[2]};
// Transfer them to the PS object
  PS.set_dim(2);
  PS.set_mom(0,mom_1);
  PS.set_mom(1,mom_2);
  double weight = 1./(16.*std::pow(pi,2)) * p12/masses[0] *4.*pi;
  return weight;
};


double set_PS_threebody(const double* randoms, const double* masses, PSpoint& PS){

  double pi = k_constants::pi;

  double mQ = std::sqrt( pow(masses[2]+masses[3],2) +
    ( std::pow(masses[0]-masses[1],2) - std::pow(masses[2]+masses[3],2) ) *
    randoms[0] );

  double massA[] = {masses[0], masses[1], mQ};
  double randA[] = {randoms[1], randoms[2]};
  double massB[] = {mQ, masses[2], masses[3]};
  double randB[] = {randoms[3], randoms[4]};

  PSpoint psp_aux;

  double weightA = set_PS_twobody (randA, massA, PS);

  double weightB = set_PS_twobody (randB, massB, psp_aux);

  std::vector<four_momentum> boost= SetBoost(PS.mom(1),false);

  PS.set_mom(1, LorDot(boost,psp_aux.mom(0)) );
  PS.set_dim(3);
  PS.set_mom(2, LorDot(boost,psp_aux.mom(1)) );

  double weight = weightA * weightB *
    (std::pow(masses[0]-masses[1],2) - std::pow(masses[2]+masses[3],2)) / (2.*pi);
  return weight;

};
double set_PS_threebody(const double* randoms, const double* masses,
  double mQ2aa, double mQ2bb, PSpoint& PS){

  double pi = k_constants::pi;

  double mQ2min = std::max(std::pow(masses[2]+masses[3],2), mQ2aa);
  double mQ2max = std::min(std::pow(masses[0]-masses[1],2), mQ2bb);
  if(mQ2max < mQ2min) mQ2max = mQ2min;

  double mQ = std::sqrt( mQ2min + (mQ2max - mQ2min) * randoms[0] );

  double massA[] = {masses[0], masses[1], mQ};
  double randA[] = {randoms[1], randoms[2]};
  double massB[] = {mQ, masses[2], masses[3]};
  double randB[] = {randoms[3], randoms[4]};

  PSpoint psp_aux;

  double weightA = set_PS_twobody (randA, massA, PS);

  double weightB = set_PS_twobody (randB, massB, psp_aux);

  std::vector<four_momentum> boost= SetBoost(PS.mom(1),false);

  PS.set_mom(1, LorDot(boost,psp_aux.mom(0)) );
  PS.set_dim(3);
  PS.set_mom(2, LorDot(boost,psp_aux.mom(1)) );

  double weight = weightA * weightB *(mQ2max - mQ2min) / (2.*pi);

  return weight;

};


double nonQCD2to2PS(const double* x, double en, PSpoint& PS){

  PSpoint PSfin;
  const double energy[3] = {en, 0., 0.};

  double weight = set_PS_twobody(x, energy, PSfin);

  four_momentum mom1(4, 0.);
  four_momentum mom2(4, 0.);
  mom1[0] = mom1[3] = mom2[0] = en/2.;
  mom2[3] = -en/2.;
  PS.set_dim(4);
  PS.set_mom(0, mom1);
  PS.set_mom(1, mom2);
  PS.set_mom(2, PSfin.mom(0));
  PS.set_mom(3, PSfin.mom(1));

  PS.set_products();

  return weight;
};

double nonQCD2to3PS(const double* x, double en, PSpoint& PS){

  PSpoint PSfin;
  const double energy[3] = {en, 0., 0.};

  double weight = set_PS_threebody(x, energy, PSfin);

  four_momentum mom1(4, 0.);
  four_momentum mom2(4, 0.);
  mom1[0] = mom1[3] = mom2[0] = en/2.;
  mom2[3] = -en/2.;
  PS.set_dim(5);
  PS.set_mom(0, mom1);
  PS.set_mom(1, mom2);
  for(int ii=2; ii<5; ii++)
    PS.set_mom(ii, PSfin.mom(ii));

  PS.set_products();

  return weight;

};
double nonQCD2to3PS(const double* x, double en, double mQ2aa, double mQ2bb, PSpoint& PS){

  PSpoint PSfin;
  const double energy[3] = {en, 0., 0.};

  double weight = set_PS_threebody(x, energy, mQ2aa, mQ2bb, PSfin);

  four_momentum mom1(4, 0.);
  four_momentum mom2(4, 0.);
  mom1[0] = mom1[3] = mom2[0] = en/2.;
  mom2[3] = -en/2.;

  PS.set_dim(5);
  PS.set_mom(0, mom1);
  PS.set_mom(1, mom2);
  for(int ii=2; ii<5; ii++)
    PS.set_mom(ii, PSfin.mom(ii));

  PS.set_products();

  return weight;

};

double QCD_0_2to3PS(const double* x, double en, PSpoint& PS){

  PSpoint PSfin;

// For the moment, I use the first two randoms as momentum fractions

  double x1 = x[0]; double x2 = x[1];
  double en_hat = std::sqrt(x1*x2)*en;


// First generate final state momenta in partonic CM frame

  const double energy[3] = {en_hat, 0., 0.};

  double weight = set_PS_threebody(x+2, energy, PSfin);


// Then boost back to hadronic CM frame and set the full PS
// Note that I'm using all-outgoing formalism, so initial state momenta have
// negative energies

  four_momentum mom1(4, 0.);
  four_momentum mom2(4, 0.);
  mom1[0] = -x1*en/2.; mom1[3] = x1*en/2.;
  mom2[0] = -x2*en/2.; mom2[3] = -x2*en/2.;

  PS.set_dim(5);
  PS.set_mom(0, mom1);
  PS.set_mom(1, mom2);

  four_momentum mom_temp;
  for(int ii=0; ii<4; ii++) mom_temp[ii] = -mom1[ii] - mom2[ii];
  std::vector<four_momentum> boost = SetBoost(mom_temp, false);

  for(int ii=0; ii<3; ii++){
    mom_temp = LorDot(boost, PSfin.mom(ii));
    PS.set_mom(ii+2, mom_temp);
  }

  PS.set_products();

//  std::cout << "\n" << "x1 and x2: " << x1 << " - " << x2 << "\n";
//  PS_checker(0);

  return weight;

};

double QCD_0_2to2PS(const double* x, double en, PSpoint& PS){

  PSpoint PSfin;

// For the moment, I use the first two randoms as momentum fractions

  double x1 = x[0]; double x2 = x[1];
  double en_hat = std::sqrt(x1*x2)*en;


// First generate final state momenta in partonic CM frame

  const double energy[3] = {en_hat, 0., 0.};

  double weight = set_PS_twobody(x+2, energy, PSfin);


// Then boost back to hadronic CM frame and set the full PS
// Note that I'm using all-outgoing formalism, so initial state momenta have
// negative energies

  four_momentum mom1(4, 0.);
  four_momentum mom2(4, 0.);
  mom1[0] = -x1*en/2.; mom1[3] = x1*en/2.;
  mom2[0] = -x2*en/2.; mom2[3] = -x2*en/2.;
  PS.set_dim(4);
  PS.set_mom(0, mom1);
  PS.set_mom(1, mom2);

  four_momentum mom_temp;
  for(int ii=0; ii<4; ii++) mom_temp[ii] = -mom1[ii] - mom2[ii];
  std::vector<four_momentum> boost = SetBoost(mom_temp, false);

  for(int ii=0; ii<2; ii++){
    mom_temp = LorDot(boost, PSfin.mom(ii));
    PS.set_mom(ii+2, mom_temp);
  }

  PS.set_products();

//  std::cout << "\n" << "x1 and x2: " << x1 << " - " << x2 << "\n";
//  PS_checker(0);

  return weight;

};


double QCD_1_2to3PS(const double* x, double en, double mQ2aa, double mQ2bb, PSpoint& PS){

  PSpoint PSfin;

// For the moment, I use the first two randoms as momentum fractions

  double x1 = x[0]; double x2 = x[1];
  double en_hat = std::sqrt(x1*x2)*en;


// First generate final state momenta in partonic CM frame

  const double energy[3] = {en_hat, 0., 0.};

  double weight = set_PS_threebody(x+2, energy, mQ2aa, mQ2bb, PSfin);


// Then boost back to hadronic CM frame and set the full PS
// Note that I'm using all-outgoing formalism, so initial state momenta have
// negative energies

  four_momentum mom1(4, 0.);
  four_momentum mom2(4, 0.);
  mom1[0] = -x1*en/2.; mom1[3] = x1*en/2.;
  mom2[0] = -x2*en/2.; mom2[3] = -x2*en/2.;
  PS.set_dim(5);
  PS.set_mom(0, mom1);
  PS.set_mom(1, mom2);

  four_momentum mom_temp;
  for(int ii=0; ii<4; ii++) mom_temp[ii] = -mom1[ii] - mom2[ii];
  std::vector<four_momentum> boost = SetBoost(mom_temp, false);

  for(int ii=0; ii<3; ii++){
    mom_temp = LorDot(boost, PSfin.mom(ii));
    PS.set_mom(ii+2, mom_temp);
  }

  PS.set_products();

//  std::cout << "\n" << "x1 and x2: " << x1 << " - " << x2 << "\n";
//  PS_checker(0);

  return weight;

};



void set_PS_fromfile(const char* filename, PSpoint& PS){

    std::fstream myfile(filename, std::ios_base::in);

    double aa1;
    double grid[50];
    int ii = 0;

    while (myfile >> aa1)
    {
        grid[ii]=aa1; ii++;
    }

    int dim = ii/4;

    if (dim > 0){

      PS.set_dim(dim);
      for(int jj=0; jj<dim; jj++){
        four_momentum mom;
        mom[0]=grid[4*jj+3]; mom[1]=grid[4*jj+0]; mom[2]=grid[4*jj+1]; mom[3]=grid[4*jj+2];
        PS.set_mom(ii, mom);
      }
    }

};
