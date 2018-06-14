#ifndef _partons_ccH_
#define _partons_ccH_

#include <cmath>
#include <iostream>

#include "mstwpdf.h"


extern "C" {

  void fiteador_(double& xtauf, double& muf2, int& energy_sector, int& pdf_label, int& aa_in);

  void writepdfout_(const char* filename, double& xx, double& muf, int& energy_sector,
                    int& pdf_label);

  void partons_cc_(double& sq2, double& sx,
                 double& fx0,
                 double& fx1, double& fx2, double& fx3, double& fx4, double& fx5,
                 double& fxm1, double& fxm2, double& fxm3, double& fxm4, double& fxm5,
                 int& nf, int& isetproton, int& ippbar);

}



extern c_mstwpdf* pdf_tofit;

#endif
