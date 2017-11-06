#ifndef _lorentzH_
#define _lorentzH_

#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>

// LORENTZ SPACE RELATED FUNCTIONS, KISS VERSION

// I use std::vector<double> for 4-momenta, aliased as "four_momentum"
// to allow easier replacing in the future. There is no "4-momentum class",
// just a (small) bunch of functions which operate on 4-dimensional arrays.

// The space-time dimension is HARDCODED for now


typedef std::vector<double> four_momentum;

double LorDot(const four_momentum&, const four_momentum&);
four_momentum LorDot(const std::vector<four_momentum>&, const four_momentum&);
std::vector<four_momentum> LorDot(const std::vector<four_momentum>&,
                                          const std::vector<four_momentum>&);

double LorNorm(const four_momentum&);

std::complex<double> SpinProdA(const four_momentum&, const four_momentum&);

std::vector<four_momentum > SetBoost(const four_momentum&, bool);

void LorPrint(const four_momentum&);
void LorPrint(const four_momentum&, std::ofstream&);

#endif
