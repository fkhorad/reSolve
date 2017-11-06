
#include "InputPars.h"

#include <fstream>

#include "read_data.h"


/*
  InputPars object initializator: reads inputs from file "infilename"
  (assigns default values, if they exist, to pars not found on file).
  Does the reading line by line. Looks for a specific string TAG for each
  parameter. The order of parameters in the input file is not significant.
*/


void InputPars::basic_input(int argc, char* argv[]){

// Process command-line arguments

  std::string infilename;
  if(argc==1){
    infilename.assign("input.dat");
  }
  else{
    infilename.assign(argv[1]);
  }
  if(argc>2){
    machine_tag.assign(argv[2]);
    if(argc>3) std::cout << "WARNING: command arguments beyond the 2nd currently ignored";
  }


// Open input file
  std::ifstream infile(infilename.c_str());
  if(infile.fail()){
    std::cout << "Error opening input file (maybe file missing or wrong dir?) - Stopping.\n";
    exit(EXIT_FAILURE);
  }

// Record the name of the input file
  filename_0 = infilename;


// NOW THE PARAMETERS
// TRY to keep the same ordering as in the header file
// The sequence is: assign a default value (if not already assigned),
// try to read the parameter. Any succesfully read value is dumped to screen (the program is rather verbose for now)
// All is inside a loop, so ordering in the input file does not matter

  std::string line;
  std::cout << "Main input file: " << filename_0 << "\n\nInput parameters:\n\n";
  do{

    std::getline(infile, line);

    if(line[0]!='#'){

// Basic parameters: process flag, mu_R, mu_F, total energy ...

// process, order, initial state hadrons
      set_def<int>(default_ints,"process",1);
      read_data<int>(input_status,default_ints,line,"process",process);
//
      set_def<int>(default_ints,"order",0);
      read_data<int>(input_status,default_ints,line,"order",order);
//
      set_def<int>(default_ints,"ih1",1);
      read_data<int>(input_status,default_ints,line,"ih1",ih1);
      set_def<int>(default_ints,"ih2",1);
      read_data<int>(input_status,default_ints,line,"ih2",ih2);
      set_def<int>(default_ints,"resum_flag",1);
      read_data<int>(input_status,default_ints,line,"resum_flag",resum_flag);

// ALL ENERGIES IN GEV!!

// Renormalization & factorization scales
      set_def<int>(default_ints,"muR_flag",1);
      read_data<int>(input_status,default_ints,line,"muR_flag",muR_flag);
      set_def<int>(default_ints,"muF_flag",1);
      read_data<int>(input_status,default_ints,line,"muF_flag",muF_flag);
//
      set_def<double>(default_reals,"mu_R",1.);
      read_data<double>(input_status,default_reals,line,"mu_R",mu_R);
      set_def<double>(default_reals,"mu_F",1.);
      read_data<double>(input_status,default_reals,line,"mu_F",mu_F);
      set_def<double>(default_reals,"mu_min",20.);
      read_data<double>(input_status,default_reals,line,"mu_min",mu_min);

// CM energy
      set_def<double>(default_reals,"CM_energy",13000.);
      read_data<double>(input_status,default_reals,line,"CM_energy",CM_energy);

// Verbosity
      set_def<int>(default_ints,"verbosity",1);
      read_data<int>(input_status,default_ints,line,"verbosity",verbosity);


// Resummation-specific stuff that probably should be moved elsewhere
      set_def<double>(default_reals,"mu_S",1.);
      read_data<double>(input_status,default_reals,line,"mu_S",mu_S);
      set_def<double>(default_reals,"en_sec_multiplier",2.);
      read_data<double>(input_status,default_reals,line,"en_sec_multiplier",en_sec_multiplier);
      set_def<double>(default_reals,"QQ_Min",0.);
      read_data<double>(input_status,default_reals,line,"QQ_Min",QQ_Min);
      set_def<double>(default_reals,"QQ_Max",500.);
      read_data<double>(input_status,default_reals,line,"QQ_Max",QQ_Max);
      set_def<double>(default_reals,"QT_Min",0.);
      read_data<double>(input_status,default_reals,line,"QT_Min",QT_Min);
      set_def<double>(default_reals,"QT_Max",100.);
      read_data<double>(input_status,default_reals,line,"QT_Max",QT_Max);
      set_def<double>(default_reals,"eta_Min",-2.5);
      read_data<double>(input_status,default_reals,line,"eta_Min",eta_Min);
      set_def<double>(default_reals,"eta_Max",2.5);
      read_data<double>(input_status,default_reals,line,"eta_Max",eta_Max);
      set_def<double>(default_reals,"ggnp",0.);
      read_data<double>(input_status,default_reals,line,"ggnp",ggnp);
      set_def<double>(default_reals,"gqnp",0.);
      read_data<double>(input_status,default_reals,line,"gqnp",gqnp);
// For De-Quadrature (intde2.cc)
      set_def<int>(default_ints,"lenaw",8000);
      read_data<int>(input_status,default_ints,line,"lenaw",lenaw);
      set_def<double>(default_reals,"tiny",1.e-307);
      read_data<double>(input_status,default_reals,line,"tiny",tiny);
      set_def<double>(default_reals,"de_eps",0.01);
      read_data<double>(input_status,default_reals,line,"de_eps",de_eps);
// Specific PDF fit file (default none)
      set_def<std::string>(default_strings,"pdffit_file","");
      read_data<std::string>(input_status,default_strings,line,"pdffit_file",pdffit_file);
      set_def<double>(default_reals,"pdf_fussiness",0.01);
      read_data<double>(input_status,default_reals,line,"pdf_fussiness",pdf_fussiness);
      set_def<int>(default_ints,"pdf_fitonly",0);
      read_data<int>(input_status,default_ints,line,"pdf_fitonly",pdf_fitonly);
      


// PARALLELIZATION: Multi-machine runs -- UNDER TESTING NOW
      set_def<int>(default_ints,"multi_machine",0.);
      read_data<int>(input_status,default_ints,line,"multi_machine",multi_machine);


// PDFs
      set_def<int>(default_ints,"pdf_flag",82);
      read_data<int>(input_status,default_ints,line,"pdf_flag",pdf_flag);
      // set_def<std::string>(default_strings,"prefix","Grids/mstw2008nnlo");
      // read_data<std::string>(input_status,default_strings,line,"prefix",prefix);


// Some basic constants in SM
      set_def<int>(default_ints,"Nc",3);
      read_data<int>(input_status,default_ints,line,"Nc",Nc);
      set_def<int>(default_ints,"Nf",5); // Flavours for fixed flavour scheme
      read_data<int>(input_status,default_ints,line,"Nf",Nf);
      set_def<double>(default_reals,"alpha_QED",1./137.0359997);
      read_data<double>(input_status,default_reals,line,"alpha_QED",alpha_QED);
      set_def<double>(default_reals,"gevm2topb",389379365.6);
      read_data<double>(input_status,default_reals,line,"gevm2topb",gevm2topb);
      set_def<double>(default_reals,"M_p",0.9382720813);
      read_data<double>(input_status,default_reals,line,"M_p",M_p);


// INTEGRATOR-RELATED STUFF
      set_def<int>(default_ints,"integrator_flag",1);
      read_data<int>(input_status,default_ints,line,"integrator_flag",integrator_flag);
      set_def<int>(default_ints,"save_events",1);
      read_data<int>(input_status,default_ints,line,"save_events",save_events);

// Needed for option: integrator_flag=0
      set_def<int>(default_ints,"no_eventsin",100);
      read_data<int>(input_status,default_ints,line,"no_eventsin",no_eventsin);
      set_def<std::string>(default_strings,"randoms_file","");
      read_data<std::string>(input_status,default_strings,line,"randoms_file",randoms_file);



// CUBA (Vegas)-specific inputs (see CUBA library documentation)
      set_def<int>(default_ints,"ndim",1);
      read_data<int>(input_status,default_ints,line,"ndim",ndim);
      set_def<int>(default_ints,"ncomp",1);
      read_data<int>(input_status,default_ints,line,"ncomp",ncomp);
      set_def<int>(default_ints,"nvec",1);
      read_data<int>(input_status,default_ints,line,"nvec",nvec);
      set_def<double>(default_reals,"epsrel",0.);
      read_data<double>(input_status,default_reals,line,"epsrel",epsrel);
      set_def<double>(default_reals,"epsabs",0.);
      read_data<double>(input_status,default_reals,line,"epsabs",epsabs);
      set_def<int>(default_ints,"VegFlags",3);
      read_data<int>(input_status,default_ints,line,"VegFlags",VegFlags);
      set_def<int>(default_ints,"cuba_seed",0);
      read_data<int>(input_status,default_ints,line,"cuba_seed",cuba_seed);
// Approx. min and max # of evaluations
      set_def<int>(default_ints,"mineval",0);
      read_data<int>(input_status,default_ints,line,"mineval",mineval);
      set_def<int>(default_ints,"maxeval",10000);
      read_data<int>(input_status,default_ints,line,"maxeval",maxeval);
// # of evaluation at first run, # of increases at each run and
      set_def<int>(default_ints,"nstart",1000);
      read_data<int>(input_status,default_ints,line,"nstart",nstart);
      set_def<int>(default_ints,"nincrease",1000);
      read_data<int>(input_status,default_ints,line,"nincrease",nincrease);
// # of "batch" evaluations - see CUBA documentation
      set_def<int>(default_ints,"nbatch",1000);
      read_data<int>(input_status,default_ints,line,"nbatch",nbatch);
// 0: grid is deleted at end of run; 1-10: grid is stored and the number servesas an identifier
      set_def<int>(default_ints,"gridno",0);
      read_data<int>(input_status,default_ints,line,"gridno",gridno);
// Text prefix for the state file
      set_def<std::string>(default_strings,"statefile","");
      read_data<std::string>(input_status,default_strings,line,"statefile",statefile);

// For k_vegas
      set_def<int>(default_ints,"seed",-1);
      read_data<int>(input_status,default_ints,line,"seed",seed);


// Directory where output events are stored
      set_def<std::string>(default_strings,"workdir","events/");
      read_data<std::string>(input_status,default_strings,line,"workdir",workdir);

    }
  }
  while(!infile.eof());

}
// End of initializator


// Auxiliary function:
// dumps the values of any parameters which still have their default values at the end
// of reading. If any par has status '-1' - which means that the code attempted to read it,
// but found no value, either default or user-defined - exit with an error.

void InputPars::dump_default_parameters(){

  for (std::map<std::string, int>::iterator it=input_status.begin();
       it!=input_status.end(); it++){

    std::string tag;
    std::map<std::string, const int>::iterator it1;
    std::map<std::string, const double>::iterator it2;
    std::map<std::string, const std::string>::iterator it3;

    tag = it->first;
    if(it->second == 1){
      std::cout << tag << ": ";
      it1 = default_ints.find(tag);
      it2 = default_reals.find(tag);
      it3 = default_strings.find(tag);
      if(it1 != default_ints.end()) std::cout << default_ints[tag];
      if(it2 != default_reals.end()) std::cout << default_reals[tag];
      if(it3 != default_strings.end()) std::cout << default_strings[tag];
      std::cout << " (default value)\n";
    }
    if(it->second == -1){
      std::cout << "Parameter "<< tag << " missing in input file and no default value found: STOPPING" << "\n";
      exit(EXIT_FAILURE);
    }
  }
  std::cout << "\n";

}
