#ifndef _PreProcH_
#define _PreProcH_

#include "mstwpdf.h"
class all_info;


// Wrapper around Fortran code for alpha_S initialisation.
extern "C" {
  void initalphas_(int *IORD, double *FR2, double *MUR, double *ASMUR,
                   double *MC, double *MB, double *MT);
}


void PreProc_basic(all_info&);

#endif
