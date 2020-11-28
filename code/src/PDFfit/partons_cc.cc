
#include "partons_cc.h"


PDF_res_interface* pdf_tofit;

void partons_cc_(double& sq2, double& sx,
                 double& fx0,
                 double& fx1, double& fx2, double& fx3, double& fx4, double& fx5,
                 double& fxm1, double& fxm2, double& fxm3, double& fxm4, double& fxm5,
                 int& nf, int& isetproton, int& ippbar){

  double q = std::sqrt(sq2);

    fx0 = pdf_tofit->xfxQ(0,sx,q)/sx;
//
    fxm5 = pdf_tofit->xfxQ(-5,sx,q)/sx;
    fxm4 = pdf_tofit->xfxQ(-4,sx,q)/sx;
    fxm3 = pdf_tofit->xfxQ(-3,sx,q)/sx;
    fxm2 = pdf_tofit->xfxQ(-2,sx,q)/sx;
    fxm1 = pdf_tofit->xfxQ(-1,sx,q)/sx;
//
    fx1 = pdf_tofit->xfxQ(1,sx,q)/sx;
    fx2 = pdf_tofit->xfxQ(2,sx,q)/sx;
    fx3 = pdf_tofit->xfxQ(3,sx,q)/sx;
    fx4 = pdf_tofit->xfxQ(4,sx,q)/sx;
    fx5 = pdf_tofit->xfxQ(5,sx,q)/sx;

}
