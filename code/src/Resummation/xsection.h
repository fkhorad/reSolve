//header file for xsection.cc

#ifndef XSECTION_H
#define XSECTION_H

#include <array>
#include <vector>

#include "constants.h"

class PSdep_variables;
class ResummationInfo;
struct pdffit;


double xsection(PSdep_variables&, ResummationInfo*);

double xsection2(PSdep_variables&, ResummationInfo*);

void distributions(double x, double& U, double& D, double& US, double& DS, double& SS, double& GL, double& CH, double& BO, int ih, int ifit1, int ifit2, pdffit& pdfbeam1fit, double aapow, int verbosity);

double pdf_func(std::vector<std::array<std::array<double, k_constants::nfitmax>, k_constants::nfitpars> >& A_FL,
   double x, int ifit1, int ifit2, double aapow);

#endif
