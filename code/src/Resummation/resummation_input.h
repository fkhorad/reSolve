#ifndef _resummation_inputH_
#define _resummation_inputH_

#include "mellinspace_fns.h"
#include "pdfmellin.h"
#include "resu_procindep.h"
#include "mstwpdf.h"

struct InputPars;


extern "C" {
  double alphas_(double*);
}

extern "C" {
  void initalphas_(int*, double*, double*, double*, double*, double*, double*);
}


struct ResummationInfo {

// PDF fit coefficients for hadron 1 & 2
//  pdffit pdfbeam1fita;
//  pdffit pdfbeam2fita;
// PDF values in Mellin space along the contour
//  pdfmellin pdfbeam1pos;
//  pdfmellin pdfbeam2pos;
//  pdfmellin pdfbeam1min;
//  pdfmellin pdfbeam2min;
// Additional, constant PDF fit parameters (NOT the coefficients ifit and
// ifit2, which are PS-dependent and introduced deeper into the code)
//  fitparams fitparamsa;

// Weights & points for the inverse Mellin transform
  contour_info contoursa;
// N-independent resummation coefficients: Ac, Bc, diagonal Cc and beta functions
  resummationIndepfns resuNindep;
// N-dependent resummation coefficients: Cc and anomalous dimensions
  std::array<resuDep_N, k_constants::mellin_points> PosBranch, NegBranch;

// Workspace for inverse Fourier transform in intdeo
  double* aw;

// Input or PreProc parameters relevant for resummation
  double QQ_Min, QQ_Max, QT_Min, QT_Max, eta_Min, eta_Max;
  bool auto_qtlim, auto_etalim;
  double mu_R, mu_F, mu_S;
  int muR_flag, muF_flag, verbosity, ih1, ih2, Nf, Nc;
  double ggnp, gqnp;
  c_mstwpdf* pdf; // Dynamic memory allocation, at least for now
  int lenaw;
  double tiny, de_eps;
  double Ca, Cf, mu_min, en_sec_multiplier;
  double CM_energy, gevm2topb;
  double alpha_QED;
  int pdf_flag;
  int resum_flag, fitonly;
  std::string pdffit_file;
  double pdf_fussiness;
  double M_p, mTop;

  int tot_em_charge; // signed, single digit; +10 if sign not measured, i.e. 11 means +-1
  int qq_order, gg_order, order, pcF;
  bool qq_initiated;

};


void resummation_input(std::string filename, ResummationInfo&);


#endif
