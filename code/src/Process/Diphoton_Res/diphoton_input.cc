#include "diphoton_input.h"


// Routine that reads the additional input parameters for this specific process
// (pp->yy). Written in the same style of the generic input reader function in InputPars.cc.

void ReadInput_diphoton(diphoton_input* diph_in,
                        std::map<std::string, int>& input_status,
                        std::map<std::string, const int>& default_ints,
                        std::map<std::string, const double>& default_reals){

// Open the additional parameters input file, by default the same as the basic one.

  std::ifstream infile(diph_in->filename.c_str());

// The additional parameters in this case are mainly for kinematic cuts


// Do the reading
  std::string line;
  do{

    std::getline(infile, line);

    set_def<int>(default_ints,"boxflag",1);
    read_data<int>(input_status,default_ints,line,"boxflag",diph_in->boxflag);

    set_def<double>(default_reals,"etaCut",5.);
    read_data<double>(input_status,default_reals,line,"etaCut",diph_in->etaCut);
    set_def<double>(default_reals,"crack1",1.);
    read_data<double>(input_status,default_reals,line,"crack1",diph_in->crack1);
    set_def<double>(default_reals,"crack1",1.);
    read_data<double>(input_status,default_reals,line,"crack2",diph_in->crack2);
    set_def<double>(default_reals,"pT1cut",0.);
    read_data<double>(input_status,default_reals,line,"pT1cut",diph_in->pT1cut);
    set_def<double>(default_reals,"pT2cut",0.);
    read_data<double>(input_status,default_reals,line,"pT2cut",diph_in->pT2cut);
    set_def<double>(default_reals,"Rcut",0.);
    read_data<double>(input_status,default_reals,line,"Rcut",diph_in->Rcut);

  }
  while(!infile.eof());

}
