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

// PreProcessing - 1
  resu_preproc(event_info, drellyan_in.res_1);

// PreProcessing - 2
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

      set_def<int>(default_ints,"DYprocess",5);
      read_data<int>(input_status,default_ints,line,"DYprocess",drellyan_in.DYprocess);
      set_def<int>(default_ints,"DYnarrowwidthapprox",0);
      read_data<int>(input_status,default_ints,line,"DYnarrowwidthapprox",drellyan_in.DYnarrowwidthapprox);

      set_def<double>(default_reals,"gf",0.000011663787);
      read_data<double>(input_status,default_reals,line,"gf",drellyan_in.gf);
      set_def<double>(default_reals,"mz",91.1876);
      read_data<double>(input_status,default_reals,line,"mz",drellyan_in.mz);
      set_def<double>(default_reals,"mw",80.398);
      read_data<double>(input_status,default_reals,line,"mw",drellyan_in.mw);

      set_def<double>(default_reals,"mzp",10000); // Z prime mass - set as 10 TeV as default
      read_data<double>(input_status,default_reals,line,"mzp",drellyan_in.mzp);

      set_def<double>(default_reals,"crack1",1.0);
      read_data<double>(input_status,default_reals,line,"crack1",drellyan_in.crack1); //general cut
      set_def<double>(default_reals,"crack2",1.0);
      read_data<double>(input_status,default_reals,line,"crack2",drellyan_in.crack2); //general cut
      set_def<double>(default_reals,"pT1cut",0.0);
      read_data<double>(input_status,default_reals,line,"pT1cut",drellyan_in.pT1cut);  //general cut for DY Z
      set_def<double>(default_reals,"pT2cut",0.0);
      read_data<double>(input_status,default_reals,line,"pT2cut",drellyan_in.pT2cut);  //general cut for DY Z
      set_def<double>(default_reals,"eta1cut",5.);
      read_data<double>(input_status,default_reals,line,"eta1cut",drellyan_in.eta1cut);  //general cut for DY Z
      set_def<double>(default_reals,"eta2cut",5.);
      read_data<double>(input_status,default_reals,line,"eta2cut",drellyan_in.eta2cut);  //general cut for DY Z
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

  dump_default_parameters(input_status, default_ints, default_reals, default_strings);

}

void drellyan_preproc(drellyan_input& drellyan_in){

  drellyan_in.ndim = 6;
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

  std::cout << "DY process currently uses mw, mz and GF as inputs\n";
  std::cout << "alpha_QED and sW are calculated\n\n";

  double mw = drellyan_in.mw;
  double mz = drellyan_in.mz;
  double mzp = drellyan_in.mzp; // Z prime mass
  double gf = drellyan_in.gf;
  double pi = k_constants::pi;

  drellyan_in.mw = mw;
  drellyan_in.mz = mz;

  drellyan_in.cw2 = mw*mw/mz/mz;
  double sw2 = 1. - drellyan_in.cw2;
  drellyan_in.sw2 = sw2;
//
  std::cout << "(Re)calculating alpha_QED value" << std::endl;
  drellyan_in.res_1.alpha_QED = std::sqrt(2.)*gf*mw*mw*sw2/pi;
  std::cout << "alpha_QED = " << drellyan_in.res_1.alpha_QED << std::endl;
  drellyan_in.gz = std::sqrt(std::sqrt(2.0)*gf*mz*mz);
  drellyan_in.gw = std::sqrt(4.0*std::sqrt(2.0)*gf*mw*mw);
  drellyan_in.gzp = std::sqrt(std::sqrt(2.0)*gf*mzp*mzp); // Set Z prime coupling using same form as Z coupling 

// W and Z widths -- hardcoded for now
  drellyan_in.zw = 2.4952;
  drellyan_in.ww = 2.141;
  drellyan_in.zpw = 0.03*mzp; // Set default Z prime width to 3 percent of Z prime mass

//CKM matrix (for W boson in DY), from 0903.2120 -- hardcoded for now
  drellyan_in.Vud = 0.97419;
  drellyan_in.Vus = 0.2257;
  drellyan_in.Vub = 0.00359;
  drellyan_in.Vcd = 0.2256;
  drellyan_in.Vcs = 0.97334;
  drellyan_in.Vcb = 0.0415;
  drellyan_in.Vtd = 0.;
  drellyan_in.Vts = 0.;
  drellyan_in.Vtb = 1.;

}
