//xsection and xsection2 do the cross-section calculations to determine the normalisation required to account for pdf fits. Called in resummed in inv_fourier.cc

#include "xsection.h"

#include "pdffit_in.h"
#include "resummation_input.h"
#include "resu_preproc.h"
#include "resu_PS.h"
#include "pdf_interface.h"


double xsection(const resu_PS& resuvars, ResummationInfo* resuminfo){

    double ax, ax1, ax2, x1, x2, muf, xsec=0;
    int Nf = resuminfo->Nf;
    int ih1 = 0, ih2 = 0;
    if (resuminfo->order==1) //so finite piece for process+jet at NLO
      {
	std::cout << "problem, should call xsection_nlojet not xsection as otherwise will call PDFs with wrong x1, x2" << std::endl;
	// avoids finite piece issue - for nlo jet case the x and eta corresponding to the overall process (and hence those needed for the PDFs) are different to those for the final state being resummed (e.g. dilepton, diphoton) as the jet carries away energy, so whilst for all the resummation etc the x and eta are needed to be used, for the PDFs need essentially x_hat and eta_hat to extract PDFs at correct x1 and x2, see below

	exit(EXIT_FAILURE);
      }
    else{

	ax = std::log(resuvars.x);
	ax1 = (ax+2*resuvars.eta)/2.0;
	ax2 = (ax-2*resuvars.eta)/2.0;
	x1 = std::exp(ax1);
	x2 = std::exp(ax2);
	muf = std::sqrt(resuvars.muf2);
	ih1 = resuminfo->ih1;
	ih2 = resuminfo->ih2;
    }
	
    
// Get PDFs
    int PDFlen = 2*Nf+1;
    double pdf1[PDFlen], pdf2[PDFlen];
	
    k_getPDF(resuminfo->pdf_res, pdf1, Nf, ih1, x1, muf);
    k_getPDF(resuminfo->pdf_res, pdf2, Nf, ih2, x2, muf);

// Put contributions together
    double QQBN = 0.;
    for (int si = 0; si < PDFlen; si++) {
      for (int sj = 0; sj < PDFlen; sj++) {
        QQBN = QQBN + pdf1[si]*pdf2[sj]*resuvars.sigmaij[si][sj];
      }
    }
    xsec = QQBN;
    //Test PS vol
    // xsec = resuvars.sigmaij[0][0]/resuvars.x; //i.e. just PSfac*spincolavg

    return xsec;
}

double xsection_nlojet(const resu_PS& resuvars, ResummationInfo* resuminfo, double ss_hat, double etaa_hat){

    double ax, ax1, ax2, x1, x2, muf, xsec=0;
    int Nf = resuminfo->Nf;
    int ih1 = 0, ih2 = 0;
    if (resuminfo->order !=1) //i.e. either Born level finite piece or resummed or CT piece
      {
      	std::cout << "problem, should call xsection not xsection_nlojet as otherwise will call PDFs with wrong x1, x2" << std::endl;
	// avoids finite piece issue - for nlo jet case the x and eta corresponding to the overall process (and hence those needed for the PDFs) are different to those for the final state being resummed (e.g. dilepton, diphoton) as the jet carries away energy, so whilst for all the resummation etc the x and eta are needed to be used, for the PDFs need essentially x_hat and eta_hat to extract PDFs at correct x1 and x2, see below
	exit(EXIT_FAILURE);
      }
    else{
	// avoids finite piece issue - for nlo jet case the x and eta corresponding to the overall process (and hence those needed for the PDFs) are different to those for the final state being resummed (e.g. dilepton, diphoton) as the jet carries away energy, so whilst for all the resummation etc the x and eta are needed to be used, for the PDFs need essentially x_hat and eta_hat to extract PDFs at correct x1 and x2:
      ax = std::log(ss_hat/pow(resuminfo->CM_energy,2)); //x_hat is now s_hat (calculated in DYjet_ps.cc for the Jacobian) divided by s rather than simply q2/s which it is without the jet, x_hat > x as the jet carries away energy, similarly then need eta_hat not eta
	ax1 = (ax+2*etaa_hat)/2.0;
	ax2 = (ax-2*etaa_hat)/2.0;
	x1 = std::exp(ax1);
	x2 = std::exp(ax2);
	// std::cout << "x_hat = " << x1*x2 << std::endl;
	muf = std::sqrt(resuvars.muf2);
	ih1 = resuminfo->ih1;
	ih2 = resuminfo->ih2;
    }
    // Get PDFs
    int PDFlen = 2*Nf+1;
    double pdf1[PDFlen], pdf2[PDFlen];

    k_getPDF(resuminfo->pdf_res, pdf1, Nf, ih1, x1, muf);
    k_getPDF(resuminfo->pdf_res, pdf2, Nf, ih2, x2, muf);

// Put contributions together
    double QQBN = 0.;
    for (int si = 0; si < PDFlen; si++) {
      for (int sj = 0; sj < PDFlen; sj++) {
        QQBN = QQBN + pdf1[si]*pdf2[sj]*resuvars.sigmaij[si][sj];
      }
    }
    xsec = QQBN;
    //Test PS vol
    // xsec = resuvars.sigmaij[0][0]/resuvars.x; //i.e. just PSfac*spincolavg

    return xsec;
}


void k_getPDF(PDF_res_interface& PDF_, double* pdf, int Nf, int ih, double x, double muf){

// PDG notation:
//f = -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6 = tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t.
// By contrast, we currently use u = Nf+ 1 and d = Nf + 2
// and g = 0.

  for(int jj=-Nf; jj<=Nf; jj++){
    int pid = jj;
    if(jj==1) pid=2;  // Up and Down are swapped w.r.t. PDG
    if(jj==2) pid=1;  // Up and Down are swapped w.r.t. PDG
    if(jj==-1) pid=-2;  // Up and Down are swapped w.r.t. PDG
    if(jj==-2) pid=-1;  // Up and Down are swapped w.r.t. PDG
    if(jj==0) pid=21; // gluon
    if (ih == 1) { //p
      pdf[Nf+jj] = PDF_.xfxQ(pid,x,muf)/x;
    }
    else if (ih == -1) { //pbar - swap all antiquark and quark pdfs
      int npid = -pid;
      if(pid==21) npid=21;
      pdf[Nf+jj] = PDF_.xfxQ(npid,x,muf)/x;
    }
    else {
      std::cout << "WARNING! ih must be 1 or -1, i.e. p or pbar!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

}



double xsection2 (const resu_PS& resuvars, ResummationInfo* resuminfo) {

    double ax, ax1, ax2, x1, x2, xsec2=0;
    int Nf, ifit, ifit2, ih1, ih2;
    if (resuminfo->order==1) //so finite piece for process+jet at NLO
      {
	std::cout << "problem, should call xsection_nlojet not xsection as otherwise will call PDFs with wrong x1, x2" << std::endl;
	// avoids finite piece issue - for nlo jet case the x and eta corresponding to the overall process (and hence those needed for the PDFs) are different to those for the final state being resummed (e.g. dilepton, diphoton) as the jet carries away energy, so whilst for all the resummation etc the x and eta are needed to be used, for the PDFs need essentially x_hat and eta_hat to extract PDFs at correct x1 and x2, see below

	exit(EXIT_FAILURE);
      }
    else{

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
    }

// Get fitted PDFs
// Note that currently the fit procedure effectively limits Nf<=5 (which gives 12 PDFs when considering gluon and photon)
    int PDFlen = 2*Nf+1;
    double FX1[PDFlen], FX2[PDFlen];
    double aapower = k_constants::aa;
    if(Nf>5){
      std::cout << "Nf > 5 currently not allowed due to PDF fit" << std::endl;
      exit(EXIT_FAILURE);
    }

    distributions(x1, FX1, Nf, ih1, ifit, ifit2, global_fitpars::pdfbeam1fita, aapower);

    distributions(x2, FX2, Nf, ih2, ifit, ifit2, global_fitpars::pdfbeam2fita, aapower);


// Put contributions together
    double QQBN = 0.;
    for (int si = 0; si < PDFlen; si++) {
      for (int sj = 0; sj < PDFlen; sj++) {
        QQBN = QQBN + FX1[si]*FX2[sj]*resuvars.sigmaij[si][sj];
      }
    }

    xsec2 = QQBN;

    return xsec2;
}

double xsection2_nlojet(const resu_PS& resuvars, ResummationInfo* resuminfo, double ss_hat, double etaa_hat) {
  double ax, ax1, ax2, x1, x2, xsec2=0;
  int Nf, ifit, ifit2, ih1, ih2;
  
  if (resuminfo->order !=1) //i.e. either Born level finite piece or resummed or CT piece
      {
      	std::cout << "problem, should call xsection not xsection_nlojet as otherwise will call PDFs with wrong x1, x2" << std::endl;
	// avoids finite piece issue - for nlo jet case the x and eta corresponding to the overall process (and hence those needed for the PDFs) are different to those for the final state being resummed (e.g. dilepton, diphoton) as the jet carries away energy, so whilst for all the resummation etc the x and eta are needed to be used, for the PDFs need essentially x_hat and eta_hat to extract PDFs at correct x1 and x2, see below
	exit(EXIT_FAILURE);
      }
    else{
	// avoids finite piece issue - for nlo jet case the x and eta corresponding to the overall process (and hence those needed for the PDFs) are different to those for the final state being resummed (e.g. dilepton, diphoton) as the jet carries away energy, so whilst for all the resummation etc the x and eta are needed to be used, for the PDFs need essentially x_hat and eta_hat to extract PDFs at correct x1 and x2:
      ax = std::log(ss_hat/pow(resuminfo->CM_energy,2)); //x_hat is now s_hat (calculated in DYjet_ps.cc for the Jacobian) divided by s rather than simply q2/s which it is without the jet, x_hat > x as the jet carries away energy, similarly then need eta_hat not eta
      ax1 = (ax+2*etaa_hat)/2.0;
      ax2 = (ax-2*etaa_hat)/2.0;
      x1 = std::exp(ax1);
      x2 = std::exp(ax2);
      Nf = resuminfo->Nf;
      ifit = resuvars.ifit;
      ifit2 = resuvars.ifit2;
      ih1 = resuminfo->ih1;
      ih2 = resuminfo->ih2;

    }

// Get fitted PDFs
// Note that currently the fit procedure effectively limits Nf<=5 (which gives 12 PDFs when considering gluon and photon)
    int PDFlen = 2*Nf+1;
    double FX1[PDFlen], FX2[PDFlen];
    double aapower = k_constants::aa;
    if(Nf>5){
      std::cout << "Nf > 5 currently not allowed due to PDF fit" << std::endl;
      exit(EXIT_FAILURE);
    }

    distributions(x1, FX1, Nf, ih1, ifit, ifit2, global_fitpars::pdfbeam1fita, aapower);

    distributions(x2, FX2, Nf, ih2, ifit, ifit2, global_fitpars::pdfbeam2fita, aapower);


// Put contributions together
    double QQBN = 0.;
    for (int si = 0; si < PDFlen; si++) {
      for (int sj = 0; sj < PDFlen; sj++) {
        QQBN = QQBN + FX1[si]*FX2[sj]*resuvars.sigmaij[si][sj];
      }
    }

    xsec2 = QQBN;

    return xsec2;
}


void distributions(double x, double* FX, int Nf, int ih, int ifit1, int ifit2, pdffit& pdfbeamfit, double aapow){

    double UV, DV;
    UV = pdf_func(pdfbeamfit.A_UV,x,ifit1,ifit2,aapow)/x;
    DV = pdf_func(pdfbeamfit.A_DV,x,ifit1,ifit2,aapow)/x;
    FX[Nf-1] = pdf_func(pdfbeamfit.A_US,x,ifit1,ifit2,aapow)/x;
    FX[Nf-2] = pdf_func(pdfbeamfit.A_DS,x,ifit1,ifit2,aapow)/x;

    FX[Nf+1] = UV + FX[Nf-1];
    FX[Nf+2] = DV + FX[Nf-2];

    FX[Nf] = pdf_func(pdfbeamfit.A_GL,x,ifit1,ifit2,aapow)/x;
    if(Nf>=3){
      FX[Nf+3] = pdf_func(pdfbeamfit.A_SS,x,ifit1,ifit2,aapow)/x;
      FX[Nf-3] = FX[Nf+3];
    }
    if(Nf>=4){
      FX[Nf+4] = pdf_func(pdfbeamfit.A_CH,x,ifit1,ifit2,aapow)/x;
      FX[Nf-4] = FX[Nf+4];
    }
    if(Nf>=5){
      FX[Nf+5] = pdf_func(pdfbeamfit.A_BO,x,ifit1,ifit2,aapow)/x;
      FX[Nf-5] = FX[Nf+5];
    }

    if (ih == -1) {
      double Utemp, Dtemp;
      Utemp = FX[Nf+1];
      FX[Nf+1] = FX[Nf-1];
      FX[Nf-1] = Utemp;
      Dtemp = FX[Nf+2];
      FX[Nf+2] = FX[Nf-2];
      FX[Nf-2] = Dtemp;
    }

}

double pdf_func(std::vector<std::array<std::array<double, k_constants::nfitmax>, k_constants::nfitpars> >& A_FL,
   double x, int ifit1, int ifit2, double aapow){

  double res;

  res = A_FL[ifit2][0][ifit1]*std::pow(x,A_FL[ifit2][1][ifit1])*std::pow(1-x,A_FL[ifit2][2][ifit1])*(1+A_FL[ifit2][3][ifit1]*x + A_FL[ifit2][4][ifit1]*std::pow(x,0.5)+A_FL[ifit2][5][ifit1]*std::pow(x,1.5)+A_FL[ifit2][6][ifit1]*std::pow(x,2)+A_FL[ifit2][7][ifit1]*std::pow(x,aapow));

  return res;

}
