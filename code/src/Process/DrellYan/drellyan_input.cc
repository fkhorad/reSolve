#include "drellyan_input.h"

#include <cmath>
#include "constants.h"
#include "resu_preproc.h"
#include "resu_PS.h"

// Routine that reads the additional input parameters for this specific process
// Drell-Yan. Written in the same style of the generic input reader function in InputPars.cc.

void drellyan_setup(std::string filename, const event_dumper_info& event_info, drellyan_input& drellyan_in){

// Read common input data for all resummation processes
  ResummationInfo res_1;
  resummation_input(filename, res_1);
  drellyan_in.res_1 = res_1;

// Read extra input specific to the Drell-Yan process
  drellyan_in.event_info = event_info;
  drellyan_ReadInput(filename, drellyan_in);

// PreProcessing - 1 (resummation)
  resu_preproc(event_info, drellyan_in.res_1);

// PreProcessing - 2 (DY-specific)
  drellyan_preproc(drellyan_in);

}

void drellyan_ReadInput(std::string filename, drellyan_input& drellyan_in){

// Open the additional parameters input file, by default the same as the basic one.

  std::ifstream infile(filename.c_str());

// The additional parameters in this case are mainly for kinematic cuts

  std::map<std::string, int> input_status;
  std::map<std::string, const int> default_ints;
  std::map<std::string, const double> default_reals;
  std::map<std::string, const std::string> default_strings;


// Do the reading
  std::string line;
  do{

    std::getline(infile, line);

    if(line[0]!='#'){
//
// Process flags
      set_def<int>(default_ints,"DYprocess",5);
      read_data<int>(input_status,default_ints,line,"DYprocess",drellyan_in.DYprocess);
      set_def<int>(default_ints,"DYnarrowwidthapprox",0);
      read_data<int>(input_status,default_ints,line,"DYnarrowwidthapprox",drellyan_in.DYnarrowwidthapprox);
// TEMP
      set_def<int>(default_ints,"forward",0);
      read_data<int>(input_status,default_ints,line,"forward",drellyan_in.forward);
      set_def<int>(default_ints,"backward",0);
      read_data<int>(input_status,default_ints,line,"backward",drellyan_in.backward);
      set_def<int>(default_ints,"AFB_flag",0);
      read_data<int>(input_status,default_ints,line,"AFB_flag",drellyan_in.AFB_flag);
// Calculate sigma^F - sigma^B instead of sigma -- for AFB
//

//
// Basic inputs in EW scheme
      set_def<double>(default_reals,"gf",0.000011663787);
      read_data<double>(input_status,default_reals,line,"gf",drellyan_in.gf);
      set_def<double>(default_reals,"mz",91.1876);
      read_data<double>(input_status,default_reals,line,"mz",drellyan_in.mz);
      set_def<double>(default_reals,"mw",80.398);
      read_data<double>(input_status,default_reals,line,"mw",drellyan_in.mw);
      set_def<double>(default_reals,"zw",2.4952);
      read_data<double>(input_status,default_reals,line,"zw",drellyan_in.zw);
      set_def<double>(default_reals,"ww",2.141);
      read_data<double>(input_status,default_reals,line,"ww",drellyan_in.ww);
//
// Z' parameters -- 1
      set_def<int>(default_ints,"Zlike_flag",0);
      read_data<int>(input_status,default_ints,line,"Zlike_flag",drellyan_in.Zlike_flag);
      set_def<int>(default_ints,"width_as_ratio",1);
      read_data<int>(input_status,default_ints,line,"width_as_ratio",drellyan_in.width_as_ratio);
      set_def<double>(default_reals,"mzp",10000.); // Z' mass - set as 10 TeV as default
      read_data<double>(input_status,default_reals,line,"mzp",drellyan_in.mzp);
      set_def<double>(default_reals,"gzp",1.); // Z' coupling constant
      read_data<double>(input_status,default_reals,line,"gzp",drellyan_in.gzp);
//
// CUTS
      set_def<double>(default_reals,"crack1",1.0);
      read_data<double>(input_status,default_reals,line,"crack1",drellyan_in.crack1); //general cut
      set_def<double>(default_reals,"crack2",1.0);
      read_data<double>(input_status,default_reals,line,"crack2",drellyan_in.crack2); //general cut
      set_def<double>(default_reals,"pT1cut",0.0);
      read_data<double>(input_status,default_reals,line,"pT1cut",drellyan_in.pT1cut);  //general cut for DY Z
      set_def<double>(default_reals,"pT2cut",0.0);
      read_data<double>(input_status,default_reals,line,"pT2cut",drellyan_in.pT2cut);  //general cut for DY Z
      set_def<double>(default_reals,"eta1cut",10.);
      read_data<double>(input_status,default_reals,line,"eta1cut",drellyan_in.eta1cut);  //general cut for DY Z
      set_def<double>(default_reals,"eta2cut",10.);
      read_data<double>(input_status,default_reals,line,"eta2cut",drellyan_in.eta2cut);  //general cut for DY Z
      set_def<double>(default_reals,"etacutNEW1",10.);
      read_data<double>(input_status,default_reals,line,"etacutNEW1",drellyan_in.etacutNEW1);
      set_def<double>(default_reals,"etacutNEW2",10.);
      read_data<double>(input_status,default_reals,line,"etacutNEW2",drellyan_in.etacutNEW2);
      set_def<double>(default_reals,"pTecut",0.);
      read_data<double>(input_status,default_reals,line,"pTecut",drellyan_in.pTecut);  //general cut for DY W
      set_def<double>(default_reals,"pTmisscut",0.);
      read_data<double>(input_status,default_reals,line,"pTmisscut",drellyan_in.pTmisscut);  //general cut for DY W
      set_def<double>(default_reals,"etaecut",5.);
      read_data<double>(input_status,default_reals,line,"etaecut",drellyan_in.etaecut);  //general cut for DY W
      set_def<double>(default_reals,"tmasscut",0.);
      read_data<double>(input_status,default_reals,line,"tmasscut",drellyan_in.tmasscut);  //general cut for DY W
    }
  }
  while(!infile.eof());
//
// Do a second infile sweep because some defaults are conditional
  std::ifstream infile2(filename.c_str());
  double def_zpw, def_zpcoup;
  if(drellyan_in.Zlike_flag == 1) def_zpcoup = 1.;
  else def_zpcoup = 0.;
  if(drellyan_in.width_as_ratio == 0) def_zpw = 300.;
  else def_zpw = 0.03;
  do{

    std::getline(infile2, line);

    if(line[0]!='#'){
// Z' parameters -- 2
      set_def<double>(default_reals,"zpw",def_zpw); // Z' width - set to 3% of the mass
      read_data<double>(input_status,default_reals,line,"zpw",drellyan_in.zpw);
// Z'-fermion couplings
      set_def<double>(default_reals,"cLu1",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cLu1",drellyan_in.cLu1);
      set_def<double>(default_reals,"cLu2",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cLu2",drellyan_in.cLu2);
      set_def<double>(default_reals,"cLd1",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cLd1",drellyan_in.cLd1);
      set_def<double>(default_reals,"cLd2",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cLd2",drellyan_in.cLd2);
      set_def<double>(default_reals,"cLd3",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cLd3",drellyan_in.cLd3);
      set_def<double>(default_reals,"cLl1",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cLl1",drellyan_in.cLl1);
//
      set_def<double>(default_reals,"cRu1",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cRu1",drellyan_in.cRu1);
      set_def<double>(default_reals,"cRu2",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cRu2",drellyan_in.cRu2);
      set_def<double>(default_reals,"cRd1",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cRd1",drellyan_in.cRd1);
      set_def<double>(default_reals,"cRd2",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cRd2",drellyan_in.cRd2);
      set_def<double>(default_reals,"cRd3",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cRd3",drellyan_in.cRd3);
      set_def<double>(default_reals,"cRl1",def_zpcoup);
      read_data<double>(input_status,default_reals,line,"cRl1",drellyan_in.cRl1);
    }
  }
  while(!infile2.eof());


  dump_default_parameters(input_status, default_ints, default_reals, default_strings);

}

void drellyan_preproc(drellyan_input& drellyan_in){

// Dimension of the MC integral (# of randoms)
  drellyan_in.ndim = 4; // Basic # of random for 2->2 QCD phase space
  if(drellyan_in.res_1.resum_flag == 1 || drellyan_in.res_1.CT_flag == 1){
// CT needs 4 extra vars for qT and 2 extra momentum fractions, and CT_flag takes precedence over resum_flag
    if(drellyan_in.res_1.CT_flag) drellyan_in.ndim += 4;
// Resummation needs 2 extra vars for qT (as a 2d vector)
    else if(drellyan_in.res_1.resum_flag) drellyan_in.ndim +=2;
  }
  else if(drellyan_in.res_1.order == 1){
// If no resummation/CT, adjust ndim according to order; only NLO implemented, and 3 particle PS needs 3 extra vars
    drellyan_in.ndim += 3;
  }
//
  drellyan_in.res_1.qq_order = drellyan_in.res_1.order;
  drellyan_in.res_1.gg_order = -1; // gg channel disabled
  drellyan_in.res_1.pcF = 0;

  int npart = 4;
  drellyan_in.event_info.n_particles = npart;
  double mass[] = {drellyan_in.res_1.M_p, drellyan_in.res_1.M_p, 0., 0.};
  std::vector<double> mass2 (mass, mass + npart);
  drellyan_in.event_info.mass = mass2;
  int inout[] = {-1, -1, 1, 1};
  std::vector<int> inout2 (inout, inout + npart);
  drellyan_in.event_info.inout = inout2;
  int mother1[] = {0, 0, 1, 1};
  std::vector<int> mother1_2 (mother1, mother1 + npart);
  drellyan_in.event_info.mother1 = mother1_2;
  int mother2[] = {0, 0, 2, 2};
  std::vector<int> mother2_2 (mother2, mother2 + npart);
  drellyan_in.event_info.mother2 = mother2_2;
  int PDG_num[] = {drellyan_in.res_1.ih1*2212,drellyan_in.res_1.ih2*2212,0,0};
  if(drellyan_in.DYprocess==1){
    PDG_num[2]=11;
    PDG_num[3]=-12;
    drellyan_in.res_1.tot_em_charge=1;
  }
  if(drellyan_in.DYprocess==3){
    PDG_num[2]=11;
    PDG_num[3]=-12;
    drellyan_in.res_1.tot_em_charge=11;
  }
  if(drellyan_in.DYprocess==2){
    PDG_num[2]=-11;
    PDG_num[3]=12;
    drellyan_in.res_1.tot_em_charge=-1;
  }
  if(drellyan_in.DYprocess==4 || drellyan_in.DYprocess==5 || drellyan_in.DYprocess==6){
    PDG_num[2]=11;
    PDG_num[3]=-11;
  }
  std::vector<int> PDG_num2 (PDG_num, PDG_num + npart);
  drellyan_in.event_info.PDG_num = PDG_num2;
//
  double pi = k_constants::pi;
//
// "EW scheme": using mw, mz and GF as inputs for EW sectors
  std::cout << "DY process currently uses mw, mz and GF as inputs\n";
  std::cout << "alpha_QED and sW are calculated\n\n";
//
  double mw = drellyan_in.mw;
  double mz = drellyan_in.mz;
  double gf = drellyan_in.gf;
//
  drellyan_in.cw2 = mw*mw/mz/mz;
  double sw2 = 1. - drellyan_in.cw2;
  drellyan_in.sw2 = sw2;
  std::cout << "sW^2 = " << drellyan_in.sw2 << std::endl;
//
  std::cout << "(Re)calculating alpha_QED value" << std::endl;
  drellyan_in.res_1.alpha_QED = std::sqrt(2.)*gf*mw*mw*sw2/pi;
  std::cout << "alpha_QED = " << drellyan_in.res_1.alpha_QED << std::endl;
  double gz = std::sqrt(4.*std::sqrt(2.0)*gf*mz*mz);
  drellyan_in.gz = gz;
  drellyan_in.gw = std::sqrt(4.*std::sqrt(2.0)*gf*mw*mw);
//
//
//CKM matrix (for W boson in DY) -- currently hardcoded, values from 0903.2120
  drellyan_in.Vud = 0.97419;
  drellyan_in.Vus = 0.2257;
  drellyan_in.Vub = 0.00359;
  drellyan_in.Vcd = 0.2256;
  drellyan_in.Vcs = 0.97334;
  drellyan_in.Vcb = 0.0415;
  drellyan_in.Vtd = 0.;
  drellyan_in.Vts = 0.;
  drellyan_in.Vtb = 1.;


// Z' parameters
  double zpw = drellyan_in.zpw, mzp = drellyan_in.mzp, gzp = drellyan_in.gzp;
//
  if(drellyan_in.width_as_ratio != 0) drellyan_in.zpw = zpw*mzp;
//
  if(drellyan_in.Zlike_flag == 1){
    drellyan_in.gzp = gzp*gz;
//
    double Qu = 2./3.;
    double Qd = -1./3.;
    double cLzu=(0.5-Qu*sw2);
    double cLzd=(-0.5-Qd*sw2);
    double cRzu=-Qu*sw2;
    double cRzd=-Qd*sw2;
    double cLzl=(-0.5+sw2);
    double cRzl=sw2;
    drellyan_in.cLu1 = drellyan_in.cLu1*cLzu;
    drellyan_in.cLu2 = drellyan_in.cLu2*cLzu;
    drellyan_in.cLd1 = drellyan_in.cLd1*cLzd;
    drellyan_in.cLd2 = drellyan_in.cLd2*cLzd;
    drellyan_in.cLd3 = drellyan_in.cLd3*cLzd;
    drellyan_in.cLl1 = drellyan_in.cLl1*cLzl;
//
    drellyan_in.cRu1 = drellyan_in.cRu1*cRzu;
    drellyan_in.cRu2 = drellyan_in.cRu2*cRzu;
    drellyan_in.cRd1 = drellyan_in.cRd1*cRzd;
    drellyan_in.cRd2 = drellyan_in.cRd2*cRzd;
    drellyan_in.cRd3 = drellyan_in.cRd3*cRzd;
    drellyan_in.cRl1 = drellyan_in.cRl1*cRzl;
  }

}
