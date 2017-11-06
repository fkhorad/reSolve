//Routine to calculate the pdf mellin moments, this should be done once at start for all and then matrix referred to after or will slow down code if continually reevaluated.

#ifndef PDFMELLIN_H
#define PDFMELLIN_H

#include <complex>
#include <array>
#include <vector>


#include "constants.h"

struct pdffit;
class contour_info;


struct pdfmellin {
 public:
  std::vector<std::array<std::array< std::complex<double>, k_constants::mellin_points>, k_constants::nfitmax> >
    uv, dv, us, ds, ss, gl, ch, bo;
};

void pdfmomentsoverallcalc(contour_info& contours, pdfmellin &pdf1p, pdfmellin &pdf2p, pdfmellin &pdf1m, pdfmellin &pdf2m, pdffit& pdfbeam1fit, pdffit& pdfbeam2fit, int verbosity, int nenergysectors);

void pdfmomentscalc(pdfmellin& a, std::complex<double> x, int i, pdffit& pdff1, pdffit& pdff2, int beam, int verbosity, int nenergysectors);

std::complex<double> pdf_func_mellin(std::complex<double> N, double aa,
                                     int verbosity, std::complex<double>* ca);

#endif
