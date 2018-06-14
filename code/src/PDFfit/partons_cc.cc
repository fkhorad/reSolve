
#include "partons_cc.h"

c_mstwpdf* pdf_tofit;

void partons_cc_(double& sq2, double& sx,
                 double& fx0,
                 double& fx1, double& fx2, double& fx3, double& fx4, double& fx5,
                 double& fxm1, double& fxm2, double& fxm3, double& fxm4, double& fxm5,
                 int& nf, int& isetproton, int& ippbar){

  double q = std::sqrt(sq2);

  if(isetproton==80 || isetproton==81 || isetproton == 82){
    if(ippbar==1){

// A call to the method
//   c_mstwpdf::update(x,q)
// updates the parton content to the values at (x,q^2).
// The parton contents can then be accessed in
//   c_mstwpdf::cont.upv etc.
//double upv,dnv,usea,dsea,str,sbar,chm,cbar,bot,bbar,glu,phot;

      pdf_tofit->update(sx,q);

      fx0 = pdf_tofit->cont.glu/sx;
//
      fxm5 = pdf_tofit->cont.bbar/sx;
      fxm4 = pdf_tofit->cont.cbar/sx;
      fxm3 = pdf_tofit->cont.sbar/sx;
      fxm2 = pdf_tofit->cont.dsea/sx;
      fxm1 = pdf_tofit->cont.usea/sx;
//
      fx1 = (pdf_tofit->cont.upv+pdf_tofit->cont.usea)/sx;
      fx2 = (pdf_tofit->cont.dnv+pdf_tofit->cont.dsea)/sx;
      fx3 = pdf_tofit->cont.str/sx;
      fx4 = pdf_tofit->cont.chm/sx;
      fx5 = pdf_tofit->cont.bot/sx;

    }
    else{
      std::cout << "pbar not currently supported\n";
    }
  }
  else{
    std::cout << "Only MSTW08nnlo central value for now\n";
  }
}
