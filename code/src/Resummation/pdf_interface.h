#ifndef _pdf_interfaceH_
#define _pdf_interfaceH_

#include "mstwpdf.h"

namespace LHAPDF{
  class PDF;
};
struct ResummationInfo;
//class c_mstwpdf;

extern "C" {
  double alphas_(double *MUR);
}
extern "C" {
  void initalphas_(int*, double*, double*, double*, double*, double*, double*);
}




class PDF_res_interface {

public:

  PDF_res_interface(){MSTW_set = 0; LHA_set = 0;}
  ~PDF_res_interface(){ if(MSTW_set==1) delete pdf_MSTW; }
// "Rule of 3" note: shallow copy is intended

  double xfxQ(int pid, double x, double Q) const;
  friend void LHA_init(ResummationInfo& resuminfo);
  friend void MSTW_init(ResummationInfo& resuminfo);
  double alphasQ(double) const;


private:

  c_mstwpdf* pdf_MSTW;
  int MSTW_set;
//
  LHAPDF::PDF* pdf_LHA;
  int LHA_set;

};

void MSTW_init(ResummationInfo& resuminfo);
void LHA_init(ResummationInfo& resuminfo);
void PDF_res_init(ResummationInfo& resuminfo);




#endif
