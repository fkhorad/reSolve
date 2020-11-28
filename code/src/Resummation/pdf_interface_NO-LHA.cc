
#include "pdf_interface.h"

#include <iostream>

#include "resummation_input.h"
#include "mstwpdf.h"


double PDF_res_interface::alphasQ(double Q) const{

  double res;
  if(LHA_set == 1){
    res = 0.;
  }
  else if(MSTW_set == 1){
    res = alphas_(&Q);
  }
  else res = 0.;

  return res;

}


void PDF_res_init(ResummationInfo& resuminfo){

  if(resuminfo.lha_flag == 1){
    std::cout << "Error: LHAPDFs requested, but the code was compiled without LHAPDF support" << std::endl;
    exit(EXIT_FAILURE);
  }
  else if(resuminfo.lha_flag == 0){
    MSTW_init(resuminfo);
  }
  else{
    std::cout << "Error -- LHA flag not set" << std::endl;
    exit(EXIT_FAILURE);
  }

}

void LHA_init(ResummationInfo& resuminfo){

// LHA PDF inizialisation: empty as LHA not used

}


void MSTW_init(ResummationInfo& resuminfo){

// QCD inizialisation: PDFs and alpha_s

  // MSTW08 central value -- ONLY AVAILABLE PDF OPTION FOR NOW, BUILT-IN
  // Open the central PDF grid corresponding to the prefix
    char filename[100];
    const char * prefix = resuminfo.pdf_setname.c_str();
    sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);

  // Declare a pdf object (with "new" as it'll need to survive out of here)
    std::cout << "PDF grid filename = " << filename << std::endl;
    resuminfo.pdf_res.pdf_MSTW = new c_mstwpdf(filename);

  // Call the Fortran initialisation routine with alpha_S(Q_0).
  // FR2 and mTop in principle should be taken from input; FR2 = (muF/muR)^2 Q0
  // is essentially a dumb parameter (as long as alphaSQ0 is
  // used along with it)
    double FR2 = 1.,  Q0 = 1., mTop = resuminfo.mTop;
    if( (resuminfo.muF_flag==1 && resuminfo.muR_flag==1) || (resuminfo.muF_flag==1 && resuminfo.muR_flag==1) ){
      FR2 = std::pow(resuminfo.mu_F/resuminfo.mu_R,2);
    }
    else if(resuminfo.muF_flag==0 && resuminfo.muR_flag==1){
      double muRtemp = std::sqrt(resuminfo.QQ_Max * resuminfo.QQ_Min)*resuminfo.mu_R;
      FR2 = std::pow(resuminfo.mu_F/muRtemp,2);
    }
    else if(resuminfo.muF_flag==1 && resuminfo.muR_flag==0){
      double muFtemp = std::sqrt(resuminfo.QQ_Max * resuminfo.QQ_Min)*resuminfo.mu_F;
      FR2 = std::pow(muFtemp/resuminfo.mu_R,2);
    }
    else if(resuminfo.muF_flag==0 && resuminfo.muR_flag==0){
      FR2 = std::pow(resuminfo.mu_F/resuminfo.mu_R,2);
    }

    initalphas_(&resuminfo.pdf_res.pdf_MSTW->alphaSorder,&FR2,&Q0,&resuminfo.pdf_res.pdf_MSTW->alphaSQ0,&resuminfo.pdf_res.pdf_MSTW->mCharm,&resuminfo.pdf_res.pdf_MSTW->mBottom,&mTop);

    resuminfo.pdf_res.MSTW_set = 1;

  // End MSTW08

}



double PDF_res_interface::xfxQ(int pid, double x, double Q) const{

  double res;
  if(LHA_set == 1){
    0.;
  }
  else if(MSTW_set == 1){
    int ii = pid;
    if(pid == 21) ii = 0;
    pdf_MSTW->parton(ii, x, Q);
  }
  else res = 0.;

}
