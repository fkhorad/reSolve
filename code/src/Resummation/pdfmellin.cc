//Routine to calculate the pdf mellin moments, this should be done once at start for all and then matrix referred to after or will slow down code if continually reevaluated.

#include "pdfmellin.h"

#include <iostream>

#include "pdffit_in.h"
#include "mellinspace_fns.h"


void pdfmomentsoverallcalc(contour_info& contours, pdfmellin &pdf1p, pdfmellin &pdf2p, pdfmellin &pdf1m, pdfmellin &pdf2m, pdffit& pdfbeam1fit, pdffit& pdfbeam2fit, int verbosity, int nenergysectors) {

    std::complex<double> Xn;
    Xn = std::complex<double>(0.0,0.0);
    int m_points = k_constants::mellin_points;

//Set up positive branch of contour
    for (int k=0; k<m_points; k++) {
	Xn = contours.Np[k];
//Beam1
	pdfmomentscalc(pdf1p, Xn, k, pdfbeam1fit, pdfbeam2fit, 1, verbosity,
                       nenergysectors);
//Beam2
	pdfmomentscalc(pdf2p, Xn, k, pdfbeam1fit, pdfbeam2fit, 2, verbosity,
                       nenergysectors
        );
    }
//Set up negative branch of contour
   for (int k=0; k<m_points; k++) {
	Xn = contours.Nm[k];
//Beam1
	pdfmomentscalc(pdf1m, Xn, k, pdfbeam1fit, pdfbeam2fit, 1, verbosity,
                       nenergysectors
        );
//Beam2
	pdfmomentscalc(pdf2m, Xn, k, pdfbeam1fit, pdfbeam2fit, 2, verbosity,
                       nenergysectors
        );
   }

}

void pdfmomentscalc(pdfmellin& a, std::complex<double> x, int i, pdffit& pdff1, pdffit& pdff2, int beam, int verbosity, int nenergysectors) {
    if (beam == 1) {}
    else if (beam == 2) {}
    else {
      std::cout << "beam must be 1 or 2 in pdfmomentscalc" << std::endl;
      exit(EXIT_FAILURE);
    }

    std::complex<double> N = x;

    int nfitmax = k_constants::nfitmax;
    int nfitpars = k_constants::nfitpars;
    double aa = k_constants::aa;

    std::complex<double> ca[8];

    a.uv.resize(nenergysectors);
    a.dv.resize(nenergysectors);
    a.us.resize(nenergysectors);
    a.ds.resize(nenergysectors);
    a.ss.resize(nenergysectors);
    a.gl.resize(nenergysectors);
    a.ch.resize(nenergysectors);
    a.bo.resize(nenergysectors);
    for (int ik = 0; ik<nfitmax; ik++) {
      for (int jk = 0; jk<nenergysectors; jk++) {
       if (beam == 1) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff1.A_UV[jk][kkk][ik],0.0);
         }
       }
       else if (beam == 2) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff2.A_UV[jk][kkk][ik],0.0);
         }
       }
       a.uv[jk][ik][i] = pdf_func_mellin(N, aa, verbosity, ca);

       if (beam == 1) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff1.A_DV[jk][kkk][ik],0.0);
         }
       }
       else if (beam == 2) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff2.A_DV[jk][kkk][ik],0.0);
         }
       }
       a.dv[jk][ik][i] = pdf_func_mellin(N, aa, verbosity, ca);

       if (beam == 1) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff1.A_US[jk][kkk][ik],0.0);
         }
       }
       else if (beam == 2) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff2.A_US[jk][kkk][ik],0.0);
         }
       }
       a.us[jk][ik][i] = pdf_func_mellin(N, aa, verbosity, ca);

       if (beam == 1) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff1.A_DS[jk][kkk][ik],0.0);
         }
       }
       else if (beam == 2) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff2.A_DS[jk][kkk][ik],0.0);
         }
       }
       a.ds[jk][ik][i] = pdf_func_mellin(N, aa, verbosity, ca);

       if (beam == 1) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff1.A_SS[jk][kkk][ik],0.0);
         }
       }
       else if (beam == 2) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff2.A_SS[jk][kkk][ik],0.0);
         }
       }
       a.ss[jk][ik][i] = pdf_func_mellin(N, aa, verbosity, ca);

       if (beam == 1) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff1.A_GL[jk][kkk][ik],0.0);
         }
       }
       else if (beam == 2) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff2.A_GL[jk][kkk][ik],0.0);
         }
       }
       a.gl[jk][ik][i] = pdf_func_mellin(N, aa, verbosity, ca);

       if (beam == 1) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff1.A_CH[jk][kkk][ik],0.0);
         }
       }
       else if (beam == 2) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff2.A_CH[jk][kkk][ik],0.0);
         }
       }
       a.ch[jk][ik][i] = pdf_func_mellin(N, aa, verbosity, ca);

       if (beam == 1) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff1.A_BO[jk][kkk][ik],0.0);
         }
       }
       else if (beam == 2) {
         for(int kkk=0; kkk<nfitpars; kkk++){
           ca[kkk] = std::complex<double>(pdff2.A_BO[jk][kkk][ik],0.0);
         }
       }
       a.bo[jk][ik][i] = pdf_func_mellin(N, aa, verbosity, ca);

     }  
    }
}


std::complex<double> pdf_func_mellin(std::complex<double> N, double aa,
                                     int verbosity, std::complex<double>* ca){

  std::complex<double> res;

  std::complex<double> aa_c = std::complex<double>(aa,0.);
  std::complex<double> OneComplex = std::complex<double>(1.0,0.0);
  res = ca[0]*cbeta(ca[1]-OneComplex+N,ca[2]+OneComplex,verbosity) +  ca[0]*ca[3]*cbeta(ca[1]+N,ca[2]+OneComplex,verbosity) + ca[0]*ca[4]*cbeta(ca[1]-0.5*OneComplex+N,ca[2]+OneComplex,verbosity) + ca[0]*ca[5]*cbeta(ca[1]+0.5*OneComplex+N,ca[2]+OneComplex,verbosity) + ca[0]*ca[6]*cbeta(ca[1]+OneComplex+N,ca[2]+OneComplex,verbosity) + ca[0]*ca[7]*cbeta(ca[1]-OneComplex+aa_c+N,ca[2]+OneComplex,verbosity);

  return res;
}
