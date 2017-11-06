//Routine to read in pdf fit parameters from file pdffitparams.dat

#include "pdffit_in.h"

#include <fstream>
#include <sstream>
#include <complex>
#include <iostream>


void pdffitread_in (std::string infilename, int energysector,
                    pdffit& pdfbeam1fit, pdffit& pdfbeam2fit){

  std::ifstream myfile (infilename.c_str());
  if(myfile.fail()){
    std::cout << "Error opening input pdf file (maybe file missing or wrong dir?)"
      << " - Stopping.\n";
      exit(EXIT_FAILURE);
  }

  int nfitmax = k_constants::nfitmax;
  int nfitpars = k_constants::nfitpars;
  int kk=energysector;

  if(pdfbeam1fit.A_UV.size() < kk+1){
    pdfbeam1fit.A_UV.resize(kk+1);
    pdfbeam1fit.A_DV.resize(kk+1);
    pdfbeam1fit.A_US.resize(kk+1);
    pdfbeam1fit.A_DS.resize(kk+1);
    pdfbeam1fit.A_SS.resize(kk+1);
    pdfbeam1fit.A_GL.resize(kk+1);
    pdfbeam1fit.A_CH.resize(kk+1);
    pdfbeam1fit.A_BO.resize(kk+1);
//
    pdfbeam2fit.A_UV.resize(kk+1);
    pdfbeam2fit.A_DV.resize(kk+1);
    pdfbeam2fit.A_US.resize(kk+1);
    pdfbeam2fit.A_DS.resize(kk+1);
    pdfbeam2fit.A_SS.resize(kk+1);
    pdfbeam2fit.A_GL.resize(kk+1);
    pdfbeam2fit.A_CH.resize(kk+1);
    pdfbeam2fit.A_BO.resize(kk+1);
  }

  
  std::string line;
  if (myfile.is_open()){
// Ignore the first five lines
    std::getline( myfile, line);
    std::getline( myfile, line);
    std::getline( myfile, line);
    std::getline( myfile, line);
    std::getline( myfile, line);
//
    std::getline( myfile, line);
    std::stringstream ss(line);
    for(int jj=0; jj<nfitpars; jj++){
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam1fit.A_UV[kk][jj][ii];
      }
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
    }
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam1fit.A_DV[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam1fit.A_US[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam1fit.A_DS[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam1fit.A_SS[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam1fit.A_GL[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam1fit.A_CH[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam1fit.A_BO[kk][jj][ii];
      }
    }
// Beam 2
    std::getline( myfile, line);
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam2fit.A_UV[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam2fit.A_DV[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam2fit.A_US[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam2fit.A_DS[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam2fit.A_SS[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam2fit.A_GL[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam2fit.A_CH[kk][jj][ii];
      }
    }
    std::getline( myfile, line);
    for(int jj=0; jj<nfitpars; jj++){
      std::getline( myfile, line);
      ss.clear(); ss.str(line);
      for(int ii=0; ii<nfitmax; ii++){
        ss >> pdfbeam2fit.A_BO[kk][jj][ii];
      }
    }
  }

}



void fitparams::fitparamscalc(double mu_min, double en_sec_mul,
  double qmin, double qmax, double mu_F, int muF_flag){

//pdfs are fitted as 1D functions (of x) at fixed energy. If a range of energies with nonconstant muf is considered (muF_flag=1), what we actually do is do the fit on a logarithmically equispatiated grid of energies.
    if(muF_flag == 1){
      minmuf = std::max(mu_min,qmin);
      double q_temp = minmuf;
      double muf_temp;
      if(q_temp >= qmax) mufgrid.push_back(q_temp);
      while (q_temp < qmax) {
        muf_temp =
          std::min(q_temp*std::sqrt(en_sec_mul),std::sqrt(q_temp*qmax));
        q_temp = en_sec_mul*q_temp;
        mufgrid.push_back(mu_F*muf_temp);
      }
    }
    else mufgrid.push_back(mu_F);
    en_sec_multiplier = en_sec_mul;
    nenergysectors = mufgrid.size();
}
