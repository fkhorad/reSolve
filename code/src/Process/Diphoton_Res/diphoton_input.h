#ifndef _diphoton_inputH_
#define _diphoton_inputH_

#include <string>
#include <map>
#include <fstream>

#include "read_data.h"
#include "events_out.h"
#include "resummation_input.h"

struct diphoton_input{

  ResummationInfo res_1;
  event_dumper_info event_info;

// Generic diphoton cuts
  double pT1cut, pT2cut;
  double etaCut;
  double crack1, crack2;
  double Rcut; // relative diphoton separation

// Box yes/no
  int boxflag; //include  gg -> yy box (0 = no, 1 = yes, 2 = box only

// General info
  int ndim;

};


void diphoton_setup(std::string filename, const event_dumper_info&, diphoton_input& diph_in);

void diphoton_ReadInput(std::string filename, diphoton_input& diph_in);

#endif
