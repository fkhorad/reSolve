
#include "InputPars.h"

#include <fstream>
#include <sstream>

#include "read_data.h"


/*
  InputPars object initializator: reads inputs from file "infilename"
  (assigns default values, if they exist, to pars not found on file).
  Does the reading line by line. Looks for a specific string TAG for each
  parameter. The order of parameters in the input file is not significant.
*/


void InputPars::ReadInput(int argc, char* argv[]){

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

// process
      set_def<int>(default_ints,"process",1);
      read_data<int>(input_status,default_ints,line,"process",process);

// General "verbosity" level
      set_def<int>(default_ints,"verbosity",1);
      read_data<int>(input_status,default_ints,line,"verbosity",verbosity);

// ALL ENERGIES IN GEV!!


// PARALLELISATION: Multi-machine runs -- UNDER TESTING NOW
      set_def<int>(default_ints,"multi_machine",0.);
      read_data<int>(input_status,default_ints,line,"multi_machine",multi_machine);


// PDFfit-only mode
      set_def<int>(default_ints,"pdf_fitonly",0);
      read_data<int>(input_status,default_ints,line,"pdf_fitonly",pdf_fitonly);


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
// Text prefix for the integration state file (for reusing)
      set_def<std::string>(default_strings,"statefile","");
      read_data<std::string>(input_status,default_strings,line,"statefile",statefile);
      set_def<int>(default_ints,"resume_integration",0);
      read_data<int>(input_status,default_ints,line,"resume_integration",resume_integration);
//
// CUBA (Vegas)-specific inputs (see CUBA library documentation)
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
// # of "batch" evaluations - see CUBA documentation
      set_def<int>(default_ints,"nbatch",1000);
      read_data<int>(input_status,default_ints,line,"nbatch",nbatch);
// 0: grid is deleted at end of run; 1-10: grid is stored and the number servesas an identifier
      set_def<int>(default_ints,"gridno",0);
      read_data<int>(input_status,default_ints,line,"gridno",gridno);
//
// For k_vegas
      set_def<int>(default_ints,"seed",-1);
      read_data<int>(input_status,default_ints,line,"seed",seed);


// Directory where output files are stored
      set_def<std::string>(default_strings,"workdir","events/");
      read_data<std::string>(input_status,default_strings,line,"workdir",workdir);

// Histograms
      set_def<int>(default_ints,"hist_only",0);
      read_data<int>(input_status,default_ints,line,"hist_only",hist_only);
      std::string tag_0("histo");
      std::size_t pos = line.find(tag_0);
      std::stringstream linestream;
      if(pos!=std::string::npos){
        linestream << line.substr(pos + tag_0.length() + 1);
        histo_data.push_back(linestream.str());
      }

    }
  }
  while(!infile.eof());

}
// End of initializator
