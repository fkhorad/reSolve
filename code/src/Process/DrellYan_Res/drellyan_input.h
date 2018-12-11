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

// Generic drellyan cuts
  double pT1cut, pT2cut;
  double eta1cut, eta2cut;
  double crack1, crack2;
  double Rcut; // relative diphoton separation
  int DYprocess; //Selection of mediator boson if DY, 1 = W+, 2 = W-, 3 = W+ and W-, 4 = Z only, 5 = Z+gamma, 6 = Z prime+Z+gamma
  int DYnarrowwidthapprox; //Use narrow width approx for DY so q = mw/mz
  double pTecut, pTmisscut;
  double etaecut;
  double tmasscut;

// CKM
  double Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb;
// EW
  double gf;
  double mz,zw;
  double mw,ww;
  double gz;
  double gw;
  double sw2;
  double cw2;
  double mzp, zpw; // Z prime mass and width
  double gzp; // Z prime coupling

};

void drellyan_setup(std::string filename, const event_dumper_info&, drellyan_input&);

void drellyan_ReadInput(std::string filename, drellyan_input&);

void drellyan_preproc(drellyan_input& drellyan_in);


#endif
