#ifndef _diphoton_inputH_
#define _diphoton_inputH_

#include <string>
#include <map>
#include <fstream>

#include "read_data.h"


//class diphoton_input : public user_input{
class diphoton_input{

public:

// Generic diphoton cuts
  double pT1cut, pT2cut;
  double etaCut;
  double crack1, crack2;
  double Rcut; // relative diphoton separation

// Box yes/no
  int boxflag; //include  gg -> yy box (0 = no, 1 = yes, 2 = box only

// General info
  double gevm2topb;
  double M_p;
  double CM_energy;
  double alpha_QED;
  std::string filename;
  std::string workdir;
  std::string machine_tag;
  int verbosity;
  int Nf, Nc;
  
  int integrator_flag, save_events;

};


void ReadInput_diphoton(diphoton_input*,
                        std::map<std::string, int>& input_status,
                        std::map<std::string, const int>& default_ints,
                        std::map<std::string, const double>& default_reals);

#endif
