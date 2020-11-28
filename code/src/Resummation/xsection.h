//header file for xsection.cc

#ifndef XSECTION_H
#define XSECTION_H

#include <array>
#include <vector>

#include "constants.h"

struct resu_PS;
struct ResummationInfo;
struct pdffit;
class PDF_res_interface;

double xsection(const resu_PS&, ResummationInfo*);

double xsection_nlojet(const resu_PS&, ResummationInfo*, double ss_hat, double etaa_hat);

void k_getPDF(PDF_res_interface& PDF_, double* pdf1, int Nf, int ih1, double x1, double muf);

double xsection2(const resu_PS&, ResummationInfo*);

double xsection2_nlojet(const resu_PS&, ResummationInfo*, double ss_hat, double etaa_hat);

void distributions(double x, double* FX, int Nf, int ih, int ifit1, int ifit2, pdffit& pdfbeamfit, double aapow);

double pdf_func(std::vector<std::array<std::array<double, k_constants::nfitmax>, k_constants::nfitpars> >& A_FL, double x, int ifit1, int ifit2, double aapow);

#endif
