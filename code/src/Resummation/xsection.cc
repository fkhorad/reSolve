//xsection and xsection2 do the cross-section calculations to determine the normalisation required to account for pdf fits. Called in resummed in inv_fourier.cc

#include "xsection.h"

#include "pdffit_in.h"
#include "resu_preproc.h"
#include "hardfns.h"
#include "mstwpdf.h"


double xsection(PSdep_variables& resuvars, ResummationInfo* resuminfo){

    double ax, ax1, ax2, x1, x2, muf, xsec=0;
    int Nf = resuminfo->Nf;

    ax = std::log(resuvars.x);
    ax1 = (ax+2*resuvars.eta)/2.0;
    ax2 = (ax-2*resuvars.eta)/2.0;
    x1 = std::exp(ax1);
    x2 = std::exp(ax2);
    muf = std::sqrt(resuvars.muf2);

    if (resuminfo->verbosity >= 12) {
      std::cout << "In xsection: " << std::endl;
      std::cout << "muf = " << muf << std::endl;
      std::cout << "ax1 = " << ax1 << std::endl;
      std::cout << "ax2 = " << ax2 << std::endl;
      std::cout << "x1 = " << x1 << std::endl;
      std::cout << "x2 = " << x2 << std::endl;
    }

// Get PDFs
    int PDFlen = 2*Nf+1;
    double pdf1[PDFlen],pdf2[PDFlen];
    for(int ii=-1*Nf; ii<=Nf; ii++){
      pdf1[ii+Nf] = resuminfo->pdf->parton(ii,x1,muf)/x1;
      pdf2[ii+Nf] = resuminfo->pdf->parton(ii,x2,muf)/x2;
      if (resuminfo->verbosity >= 12) {
      std::cout << ii << " " << x1 << " " << muf << " " << pdf1[ii+Nf] 
        << " " <<  pdf2[ii+Nf] << std::endl;
      }
    }
// Note in parton then use PDG notation for parton flavour (apart from gluon=0, not 21):
//f = -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6 = tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t.

    double temp1ub = pdf1[4]; //the up and down contributions are wrong way round from the partons script in mstwpdf.cc so need to swap them
    double temp2u = pdf1[6];
    double temp3db = pdf1[3];
    double temp4d = pdf1[7];
    pdf1[4] = temp3db;
    pdf1[3] = temp1ub;
    pdf1[6] = temp4d;
    pdf1[7] = temp2u;
    double temp5ub = pdf2[4]; //the up and down contributions are wrong way round from the partons script in mstwpdf.cc so need to swap them
    double temp6u = pdf2[6];
    double temp7db = pdf2[3];
    double temp8d = pdf2[7];
    pdf2[4] = temp7db;
    pdf2[3] = temp5ub;
    pdf2[6] = temp8d;
    pdf2[7] = temp6u;


// Put contributions together
    double QQBN = 0.;
    for (int si = 0; si < PDFlen; si++) {
      for (int sj = 0; sj < PDFlen; sj++) {
        QQBN = QQBN + pdf1[si]*pdf2[sj]*resuvars.sigmaij[si][sj];
        if (resuminfo->verbosity >= 50) {
          std::cout << "pdf1[" << si << "]*pdf2[" << sj << "]"
            << "*sigmaij[" << si << "][" << sj << "] = "
            << pdf1[si]*pdf2[sj]*resuvars.sigmaij[si][sj] << std::endl;
          std::cout << "pieces: " << pdf1[si] << " " << pdf2[sj] << " "
            << resuvars.sigmaij[si][sj] << std::endl;
        }
      }
    }
    xsec = QQBN;

    return xsec;
}

double xsection2 (PSdep_variables& resuvars, ResummationInfo* resuminfo) {

    double ax, ax1, ax2, x1, x2, xsec2=0;
    int Nf, ifit, ifit2, ih1, ih2;

    ax = std::log(resuvars.x);
    ax1 = (ax+2*resuvars.eta)/2.0;
    ax2 = (ax-2*resuvars.eta)/2.0;
    x1 = std::exp(ax1);
    x2 = std::exp(ax2);
    Nf = resuminfo->Nf;
    ifit = resuvars.ifit;
    ifit2 = resuvars.ifit2;
    ih1 = resuminfo->ih1;
    ih2 = resuminfo->ih2;

    if (resuminfo->verbosity >= 12) {
      std::cout << "xsection2: " << std::endl;
      std::cout << "ax = " << ax << std::endl;
      std::cout << "eta = " << resuvars.eta << std::endl;
      std::cout << "sigmaij = " << std::endl;
      for (int i=0; i<11; i++) {
        for (int j=0; j<11; j++) {
          std::cout << "sigmaij[" << i << "][" << j << " ] = "
            << resuvars.sigmaij[i][j] << std::endl;
        }
      }
      std::cout << "ax1 = " << ax1 << std::endl;
      std::cout << "ax2 = " << ax2 << std::endl;
      std::cout << "x1 = " << x1 << std::endl;
      std::cout << "x2 = " << x2 << std::endl;
      std::cout << "ifit = " << ifit << std::endl;
      std::cout << "ifit2 = " << ifit2 << std::endl;
    }


// Get fitted PDFs
// Note that currently the fit procedure effectively limits Nf<=5.
    int PDFlen = 2*Nf+1;
    double FX1[11] = {0.}, FX2[11] ={0.};
    double aapower = k_constants::aa;
    if(PDFlen>11){
      std::cout << "Nf > 5 currently not allowed due to PDF fit" << std::endl;
      exit(EXIT_FAILURE);
    }

    distributions(x1, FX1[6], FX1[7], FX1[4], FX1[3], FX1[8], FX1[5], FX1[9], FX1[10], ih1, ifit, ifit2, global_fitpars::pdfbeam1fita, aapower,
                  resuminfo->verbosity);
    FX1[2]=FX1[8];
    FX1[1]=FX1[9];
    FX1[0]=FX1[10];

    distributions(x2, FX2[6], FX2[7], FX2[4], FX2[3], FX2[8], FX2[5], FX2[9], FX2[10], ih2, ifit, ifit2, global_fitpars::pdfbeam2fita, aapower,
      resuminfo->verbosity);
    FX2[2]=FX2[8];
    FX2[1]=FX2[9];
    FX2[0]=FX2[10];

    if (resuminfo->verbosity >= 12) {
      for(int ii=0; ii<11; ii++){
        std::cout << "FX1[" << ii << "] = " << FX1[ii] << "  ";
      }
      std::cout << std::endl;
      for(int ii=0; ii<11; ii++){
        std::cout << "FX2[" << ii << "] = " << FX2[ii] << "  ";
      }
      std::cout << std::endl;
    }

// Put contributions together
    double QQBN = 0.;
    for (int si = 0; si < PDFlen; si++) {
      for (int sj = 0; sj < PDFlen; sj++) {
        QQBN = QQBN + FX1[si]*FX2[sj]*resuvars.sigmaij[si][sj]/(x1*x2);
        if (resuminfo->verbosity >= 12) {
          std::cout << "QQBN bit: "
            << FX1[si]*FX2[sj]*resuvars.sigmaij[si][sj]/(x1*x2) << std::endl;
        }
      }
    }

    xsec2 = QQBN;

    return xsec2;
}


void distributions(double x, double& U, double& D, double& US, double& DS, double& SS, double& GL, double& CH, double& BO, int ih, int ifit1, int ifit2, pdffit& pdfbeamfit, double aapow, int verbosity){

    double UV = 0., DV = 0.;
    UV = pdf_func(pdfbeamfit.A_UV,x,ifit1,ifit2,aapow);
    DV = pdf_func(pdfbeamfit.A_DV,x,ifit1,ifit2,aapow);
    US = pdf_func(pdfbeamfit.A_US,x,ifit1,ifit2,aapow);
    DS = pdf_func(pdfbeamfit.A_DS,x,ifit1,ifit2,aapow);

    U = UV + US;
    D = DV + DS;

    SS = pdf_func(pdfbeamfit.A_SS,x,ifit1,ifit2,aapow);
    GL = pdf_func(pdfbeamfit.A_GL,x,ifit1,ifit2,aapow);
    CH = pdf_func(pdfbeamfit.A_CH,x,ifit1,ifit2,aapow);
    BO = pdf_func(pdfbeamfit.A_BO,x,ifit1,ifit2,aapow);

    if (ih == -1) {
      double Utemp = 0., Dtemp = 0.;
      Utemp = U;
      U = US;
      US = Utemp;
      Dtemp = D;
      D = DS;
      DS = Dtemp;
    }
    if (verbosity >= 12) {
      std::cout << "U = " << U << std::endl;
      std::cout << "D = " << D << std::endl;
      std::cout << "US = " << US << std::endl;
      std::cout << "DS = " << DS << std::endl;
      std::cout << "SS = " << SS << std::endl;
      std::cout << "GL = " << GL << std::endl;
      std::cout << "CH = " << CH << std::endl;
      std::cout << "BO = " << BO << std::endl;
    }

}

double pdf_func(std::vector<std::array<std::array<double, k_constants::nfitmax>, k_constants::nfitpars> >& A_FL,
   double x, int ifit1, int ifit2, double aapow){

  double res;

  res = A_FL[ifit2][0][ifit1]*std::pow(x,A_FL[ifit2][1][ifit1])*std::pow(1-x,A_FL[ifit2][2][ifit1])*(1+A_FL[ifit2][3][ifit1]*x + A_FL[ifit2][4][ifit1]*std::pow(x,0.5)+A_FL[ifit2][5][ifit1]*std::pow(x,1.5)+A_FL[ifit2][6][ifit1]*std::pow(x,2)+A_FL[ifit2][7][ifit1]*std::pow(x,aapow));

  return res;

}
