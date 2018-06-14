//Header file for routine to read in pdf fit parameters from file pdffitparams.dat

#ifndef PDFFIT_IN_H
#define PDFFIT_IN_H

#include <vector>
#include <array>
#include <string>

#include "constants.h"


struct pdffit {
 public:
  std::vector<std::array<std::array<double, k_constants::nfitmax>, k_constants::nfitpars> >
    A_UV, A_DV, A_US, A_DS, A_SS, A_GL, A_CH, A_BO;

};


void pdffitread_in (std::string infilename, int energysector,
                    pdffit& pdfbeam1fit, pdffit& pdfbeam2fit);

class fitparams {
 public:
//  int ifit, ifit2;
  int nenergysectors;
  std::vector<double> mufgrid;
  double minmuf;
  double en_sec_multiplier;

  void fitparamscalc(double mu_min, double en_sec_mul,
  double qmin, double qmax, double mu_F, int muF_flag);
};


#endif
