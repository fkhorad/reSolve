#ifndef _drellyan_inputH_
#define _drellyan_inputH_

#include <string>
#include <map>
#include <fstream>

#include "read_data.h"
#include "events_out.h"
#include "resummation_input.h"


struct drellyan_input{

  ResummationInfo res_1;
  event_dumper_info event_info;

  int ndim;

// Process flags
  int DYprocess; //Selection of mediator boson if DY, 1 = W+, 2 = W-, 3 = W+ and W-, 4 = Z only, 5 = Z+gamma, 6 = Z prime+Z+gamma
  int DYnarrowwidthapprox; //Use narrow width approx for DY so q = mw/mz

// Generic drellyan cuts
  double pT1cut, pT2cut;
  double eta1cut, eta2cut, etacutNEW1, etacutNEW2;
  double crack1, crack2;
  double Rcut; // relative diphoton separation
  double pTecut, pTmisscut;
  double etaecut;
  double tmasscut;

// CKM
  double Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb;

// EW sector
  double gf;
  double mz,zw;
  double mw,ww;
  double gz;
  double gw;
  double sw2;
  double cw2;

// TEMP
  int AFB_flag; // AFB calculation
  int forward, backward;
//

// Z'
  int Zlike_flag, width_as_ratio;
//
  double mzp, zpw; // Z prime mass, width, width/mass ratio
  double gzp; // Z prime coupling
//
// Z'-fermions: no FCNC for now
  double cLu1, cLu2; // No t quark
  double cLd1, cLd2, cLd3;
  double cLl1; // We consider only one lepton species at a time
  double cRu1, cRu2; // No t quark
  double cRd1, cRd2, cRd3;
  double cRl1; // We consider only one lepton species at a time
};

void drellyan_setup(std::string filename, const event_dumper_info&, drellyan_input&);

void drellyan_ReadInput(std::string filename, drellyan_input&);

void drellyan_preproc(drellyan_input& drellyan_in);


#endif
