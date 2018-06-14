#ifndef _resu_preprocH_
#define _resu_preprocH_

struct ResummationInfo;
struct event_dumper_info;

// For pdffit
class c_mstwpdf;
extern c_mstwpdf* pdf_tofit;

// For N-indep resummation parameters initialisation
#include "resu_procindep.h"
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


void resu_preproc(const event_dumper_info&, ResummationInfo&);

void MSTW_init(ResummationInfo& resuminfo);

void pdffit_interface(const event_dumper_info& ev_info, ResummationInfo& resuminfo);

#endif
