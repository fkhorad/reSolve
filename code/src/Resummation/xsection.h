//header file for xsection.cc

#ifndef XSECTION_H
#define XSECTION_H

#include <array>
#include <vector>

#include "constants.h"

struct resu_PS;
struct ResummationInfo;
struct pdffit;
class c_mstwpdf;

double xsection(const resu_PS&, ResummationInfo*);

void k_getMSTWpdf(c_mstwpdf* PDF_, double* pdf1, int Nf, int ih1, double x1, double muf);

double xsection2(const resu_PS&, ResummationInfo*);

void distributions(double x, double* FX, int Nf, int ih, int ifit1, int ifit2, pdffit& pdfbeamfit, double aapow);

double pdf_func(std::vector<std::array<std::array<double, k_constants::nfitmax>, k_constants::nfitpars> >& A_FL, double x, int ifit1, int ifit2, double aapow);

#endif
