#ifndef _resu_preprocH_
#define _resu_preprocH_

// For N-indep resummation parameters initialisation
#include "resu_procindep.h"
// For pdffit
class c_mstwpdf;
extern c_mstwpdf* pdf_tofit;
// For N-dep resummation parameters + Mellin space PDFs initialisation
#include "mellinspace_fns.h"
#include "pdffit_in.h"
#include "pdfmellin.h"


namespace global_fitpars{
  
// PDF fit coefficients for hadron 1 & 2
  extern pdffit pdfbeam1fita;
  extern pdffit pdfbeam2fita;
// PDF values in Mellin space along the contour
  extern pdfmellin pdfbeam1pos;
  extern pdfmellin pdfbeam2pos;
  extern pdfmellin pdfbeam1min;
  extern pdfmellin pdfbeam2min;
// Additional, constant PDF fit parameters (NOT the coefficients ifit and
// ifit2, which are PS-dependent and introduced deeper into the code)
  extern fitparams fitparamsa;
  
}


  
//class vegas_resummationinfo {
class ResummationInfo {
public:

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
  double mu_R, mu_F, mu_S;
  int muR_flag, muF_flag, verbosity, order, ih1, ih2, Nf;
  double ggnp, gqnp;
  c_mstwpdf* pdf;
  int lenaw;
  double tiny, de_eps;
  double Ca, Cf, mu_min, en_sec_multiplier;
  double CM_energy;
  int pdf_flag;
  int process; //only passed to allow process = -1 to mean Born only for testing
  /* std::string prefix */
  std::string pdffit_file;
  double pdf_fussiness;
//
  int multi_machine;
  std::string machine_tag;
  

};


void resu_preproc(ResummationInfo*, c_mstwpdf*, int process);

void pdffit_interface(ResummationInfo*);

#endif
