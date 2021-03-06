
#include "cuba_interface.h"

#include <fstream>
#include <cmath>
#include <sys/types.h> // For pid, used for CUBA library
#include <unistd.h>    // For pid, used for CUBA library
#include <sstream>

#include "InputPars.h"
#include "phase_space.h"
#include "events_out.h"


void cuba_vegas_call(const InputPars& input_1, int ndim, int ncomp, integrand_t integrand, void* userdata){


// Vegas input parameters
  const int nvec = input_1.nvec,
            flags = input_1.VegFlags,
            seed = input_1.seed, mineval = input_1.mineval,
            maxeval = input_1.maxeval, nstart = input_1.nstart,
            nincrease = input_1.nincrease, nbatch = input_1.nbatch,
            gridno = input_1.gridno;
  const double epsabs = input_1.epsabs, epsrel = input_1.epsrel;
// Define statefile
  std::string pre_statefile;
  if(input_1.statefile==""){
    pre_statefile = "";
  }
  else{
    std::stringstream statefile_stream;
    statefile_stream << input_1.workdir << input_1.statefile;
    pre_statefile = statefile_stream.str();
  }
  const char* statefile = pre_statefile.c_str();
  std::cout << "statefile " << statefile << std::endl;
  int* spin = NULL;   // What is this?

  std::cout << "nincrease = " << nincrease << std::endl;

// Vegas outputs
  int neval, fail;
  double integral, error, prob;

// call
  Vegas(ndim, ncomp,
    integrand, userdata, nvec,
    epsrel, epsabs,
    flags, seed,
    mineval, maxeval,
    nstart, nincrease, nbatch,
    gridno, statefile, spin,
    &neval, &fail,
    &integral, &error, &prob);

  std::cout << "\nIntegral: " << integral << "+-" << error << "\n\n";
  std::stringstream main_outfile;
  main_outfile << input_1.workdir << "reSolve_main_out.dat";
  std::ofstream main_file(main_outfile.str().c_str(), std::ios::out);
  main_file << "Final result:\n";
  main_file << "\nIntegral: " << integral << "+-" << error << "\n\n";
  main_file.close();

  collect_events(input_1.workdir,maxeval,nstart,nincrease,input_1.machine_tag,input_1.save_events);

}


void collect_events(std::string workdir, int maxeval, int nstart, int nincrease,
                    std::string machine_tag, int save_events){

  int maxiter = 1;
  if(maxeval > nstart){
    do{maxiter++;}
    while(maxeval > maxiter*nstart + ((maxiter-1)*maxiter)/2*nincrease);
  }

  int pid;
  std::string line;
  std::stringstream pids_filename;
  pids_filename << workdir << "pids.dat" ;

  if(save_events != 1 && save_events != 2) return;

  for(int ii=1; ii<=maxiter; ii++){
    std::ifstream pids_file_in(pids_filename.str().c_str());
    while( (pids_file_in.is_open()) and (!pids_file_in.eof()) ){
// Get names
      pids_file_in >> pid;
      std::string filename = get_event_filename(save_events, workdir, machine_tag, ii);
      std::string part_filename = get_pids_filename(pid, save_events, workdir, machine_tag, ii);
// Define streams
      std::ifstream part_event_file(part_filename.c_str());
      std::ofstream event_file(filename.c_str(), std::ios::out | std::ios::app);
// Copy files
      event_file << part_event_file.rdbuf();
// Close streams, remove redundant file
      part_event_file.close();
      event_file.close();
      std::remove(part_filename.c_str());
    }
    pids_file_in.close();
  }
  remove(pids_filename.str().c_str());

};
