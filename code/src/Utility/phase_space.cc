#include "phase_space.h"

#include <iostream>
#include <fstream>


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
void PSpoint::add_mom(four_momentum mom){
  momenta.push_back(mom);
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
        SpinProductsB[i][j] = - 2.0*DotProducts[i][j]/SpinProductsA[i][j];
    };
  };

};


void PSpoint::set_dim(int new_dim){

  four_momentum empty = four_momentum(4, 0.);
  momenta.resize(new_dim, empty);

};



// Get

int PSpoint::dim() const{
   int res = momenta.size();
   return res;
};


four_momentum PSpoint::mom(int i) const{
  four_momentum res = momenta.at(i);
  return res;
};
// Readers of momenta and spinor products; allow negative indices to account for crossings.
// TAKE NOTICE OF THE OFFSET!
double PSpoint::ss(int i, int j) const{
  double res = DotProducts.at(std::abs(i)-1).at(std::abs(j)-1);
  if(i<0) res*=-1.;
  if(j<0) res*=-1.;
  return res;
};
std::complex<double> PSpoint::sa(int i, int j) const{
  std::complex<double> res = SpinProductsA.at(std::abs(i)-1).at(std::abs(j)-1);
  // std::cout << "res = " << std::real(res) << "+i*" << std::imag(res) << std::endl;
// Choose spinors phases so that i> spinors are unchanged upon crossing, while i] change sign.
  return res;
};
std::complex<double> PSpoint::sb(int i, int j) const{
  std::complex<double> res = SpinProductsB.at(std::abs(i)-1).at(std::abs(j)-1);
// Choose spinors phases so that i> spinors are unchanged upon crossing, while i] change sign.
  if(i<0) res*=-1.;
  if(j<0) res*=-1.;
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

// void set_PS_fromfile_jets(const char* filename, PSpoint& PS){

//     std::fstream myfile(filename, std::ios_base::in);

//     double aa1;
//     double grid[50];
//     int ii = 0;

//     while (myfile >> aa1)
//     {
//         grid[ii]=aa1; ii++;
//     }

//     int dim = ii/5; //divide 5 now as 5 particles per event not 4 as also jet

//     if (dim > 0){

//       PS.set_dim(dim);
//       for(int jj=0; jj<dim; jj++){
//         four_momentum mom;
//         mom[0]=grid[4*jj+3]; mom[1]=grid[4*jj+0]; mom[2]=grid[4*jj+1]; mom[3]=grid[4*jj+2];
//         PS.set_mom(ii, mom);
//       }
//     }

// };


///////////////////////////////////////////////////////////////////////////////

void set_PS_twobody(double costheta, double phi, const four_momentum P_in, double m1, double m2, PSpoint& PS){

  double sintheta = std::sqrt(1.-std::pow(costheta,2));
  double M0_2 = LorDot(P_in,P_in);
  double m1_2 = m1*m1;
  double m2_2 = m2*m2;

//
  double lambda;
  lambda = M0_2*M0_2 + m1_2*m1_2 + m2_2*m2_2 - 2.*M0_2*m1_2 - 2.*M0_2*m2_2 - 2.*m1_2*m2_2;
  double p12_2 = lambda/(4.*M0_2);
  double p12 = std::sqrt(p12_2);
  double vers[3] = {sintheta*std::cos(phi), std::sqrt(1.-std::pow(costheta,2))*std::sin(phi), costheta};
  double en1 = std::sqrt(p12_2 + m1_2);
  double en2 = std::sqrt(p12_2 + m2_2);
//
// The 4-momenta in the CM
  four_momentum qq1 = {en1, p12*vers[0], p12*vers[1], p12*vers[2]};
  four_momentum qq2 = {en2, -p12*vers[0], -p12*vers[1], -p12*vers[2]};

  // std::cout << "qq1 = " << qq1[0] << " " << qq1[1] << " " << qq1[2] << " " << qq1[3] << std::endl;
  // std::cout << "qq2 = " << qq2[0] << " " << qq2[1] << " " << qq2[2] << " " << qq2[3] << std::endl;

// Boost them to original frame
  std::vector<four_momentum> lor1 = SetBoost(P_in, false);
  four_momentum mom_1 = LorDot(lor1, qq1);
  four_momentum mom_2 = LorDot(lor1, qq2);

  // std::cout << "lor1 = " << std::endl;
  // std::cout << lor1[0][0] << " " << lor1[0][1] << " " << lor1[0][2] << " " << lor1[0][3] << std::endl;
  // std::cout << lor1[1][0] << " " << lor1[1][1] << " " << lor1[1][2] << " " << lor1[1][3] << std::endl;
  // std::cout << lor1[2][0] << " " << lor1[2][1] << " " << lor1[2][2] << " " << lor1[2][3] << std::endl;
  // std::cout << lor1[3][0] << " " << lor1[3][1] << " " << lor1[3][2] << " " << lor1[3][3] << std::endl;


  // std::cout << "mom_1 = " << mom_1[0] << " " << mom_1[1] << " " << mom_1[2] << " " << mom_1[3] << std::endl;
  // std::cout << "mom_2 = " << mom_2[0] << " " << mom_2[1] << " " << mom_2[2] << " " << mom_2[3] << std::endl;
  // std::cout << "mom_1.mom_1 = " << LorDot(mom_1,mom_1) << std::endl;
  // std::cout << "mom_2.mom_2 = " << LorDot(mom_2,mom_2) << std::endl;


// Transfer them to the PS object
  PS.add_mom(mom_1);
  PS.add_mom(mom_2);
}
//
void set_PS_threebody(double cos1, double phi1, double mQ, double cos2, double phi2, const four_momentum P_in, double m1, double m2, double m3, PSpoint& PS){

  set_PS_twobody (cos1, phi1, P_in, m1, mQ, PS);

  PSpoint psp_aux;
  set_PS_twobody (cos2, phi2, PS.mom(1), m2, m3, psp_aux);

  PS.set_mom(PS.dim()-1, psp_aux.mom(0));
  PS.add_mom(psp_aux.mom(1));

}
