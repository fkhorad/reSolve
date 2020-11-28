#include "diphoton_input.h"

#include "resu_preproc.h"

void diphoton_setup(std::string filename, const event_dumper_info& event_info, diphoton_input& diph_in){

// Read common input data for all resummation processes
  ResummationInfo res_1;
  resummation_input(filename, res_1);
  diph_in.res_1 = res_1;

// Read extra input specific to the diphoton process
  diph_in.event_info = event_info;
  diphoton_ReadInput(filename, diph_in);

// PreProcessing - 1 (resummation)
  resu_preproc(event_info, diph_in.res_1);

// PreProcessing - 2 (diphoton-specific)
  diph_in.res_1.qq_order = diph_in.res_1.order;
  diph_in.res_1.gg_order = diph_in.res_1.order-2;
  diph_in.res_1.pcF = 0;
  diph_in.ndim = 6;
//
  int npart = 4;
  diph_in.event_info.n_particles = npart;
  double mass[] = {res_1.M_p, res_1.M_p, 0., 0.};
  std::vector<double> mass2 (mass, mass + npart);
  diph_in.event_info.mass = mass2;
  int PDG_num[] = {res_1.ih1*2212,res_1.ih2*2212,22,22};
  std::vector<int> PDG_num2 (PDG_num, PDG_num + npart);
  diph_in.event_info.PDG_num = PDG_num2;
  int inout[] = {-1, -1, 1, 1};
  std::vector<int> inout2 (inout, inout + npart);
  diph_in.event_info.inout = inout2;
  int mother1[] = {0, 0, 1, 1};
  std::vector<int> mother1_2 (mother1, mother1 + npart);
  diph_in.event_info.mother1 = mother1_2;
  int mother2[] = {0, 0, 2, 2};
  std::vector<int> mother2_2 (mother2, mother2 + npart);
  diph_in.event_info.mother2 = mother2_2;

}


// Routine that reads the additional input parameters for this specific process
// (pp->yy). Written in the same style of the generic input reader function in InputPars.cc.

void diphoton_ReadInput(std::string filename, diphoton_input& diph_in){

// Open the additional parameters input file, by default the same as the basic one.

  std::ifstream infile(filename.c_str());

  std::cout << "Reading diphoton input" << std::endl;

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

      set_def<int>(default_ints,"boxflag",1);
      read_data<int>(input_status,default_ints,line,"boxflag",diph_in.boxflag);

      set_def<double>(default_reals,"etaCut",5.);
      read_data<double>(input_status,default_reals,line,"etaCut",diph_in.etaCut);
      set_def<double>(default_reals,"crack1",1.);
      read_data<double>(input_status,default_reals,line,"crack1",diph_in.crack1);
      set_def<double>(default_reals,"crack2",1.);
      read_data<double>(input_status,default_reals,line,"crack2",diph_in.crack2);
      set_def<double>(default_reals,"pT1cut",0.);
      read_data<double>(input_status,default_reals,line,"pT1cut",diph_in.pT1cut);
      set_def<double>(default_reals,"pT2cut",0.);
      read_data<double>(input_status,default_reals,line,"pT2cut",diph_in.pT2cut);
      set_def<double>(default_reals,"Rcut",0.);
      read_data<double>(input_status,default_reals,line,"Rcut",diph_in.Rcut);

    }

  }
  while(!infile.eof());

  dump_default_parameters(input_status, default_ints, default_reals, default_strings);

}
