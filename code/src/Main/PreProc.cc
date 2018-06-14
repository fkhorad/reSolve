
#include "PreProc.h"

#include <string>
#include <sstream>
#include <fstream>
#include <complex>

#include "InputPars.h"
#include "k_vegas_interface.h"
#include "events_out.h"


// Basic preprocessing -- COMPLETELY AGNOSTIC regarding the nature of the process

int PreProc_basic(InputPars& input_1){

  int outcode = 0;


// For histo-only mode

  if(input_1.hist_only==1) outcode=1;
  else if(input_1.pdf_fitonly){
    std::cout << "PDF fit only mode" << std::endl;
    input_1.event_info.workdir = input_1.workdir;
    input_1.event_info.multi_machine = 0;
  }
  else{

// Attempt to create a main output file, which also automatically checks if WORKDIR exists.
    std::stringstream main_outfile;

// Non k_vegas-parallelized version
    if(input_1.multi_machine==0){
// If the output file already exists, stops to avoid overwriting of data.
      main_outfile << input_1.workdir << "reSolve_main_out.dat";
      std::ifstream main_file_test(main_outfile.str().c_str(), std::ios::in);
      if(!main_file_test.fail()){
        bool continuation = false;
        if(input_1.resume_integration == 1){
          std::stringstream grid_file;
          grid_file << input_1.workdir << input_1.statefile;
          std::ifstream gridin(grid_file.str().c_str(), std::ios::in);
          if( !(gridin.fail()) ){
            continuation = true;
          }
        }
        if(!continuation) {
          std::cout << "Main output file already existed - stopping to avoid overwriting data. Set a different working directory." << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      main_file_test.close();
    }
// k_vegas-parallelized version
    else if(input_1.multi_machine==1){
      if(input_1.integrator_flag!=1){
        std::cout << "Multi-machine mode only allowed for k_vegas -- stopping" << std::endl;
        exit(EXIT_FAILURE);
      }
      if(input_1.machine_tag == ""){
        std::cout << "\"machine_tag\" variable must be non-empty for parallel running with k_vegas - STOPPING" << std::endl;
        exit(EXIT_FAILURE);
      }
//
      main_outfile << input_1.workdir << "reSolve_main_out_" << input_1.machine_tag << ".dat";
    }
    else{
      std::cout << "Unrecognized value for multi_machine flag " << input_1.multi_machine << " - STOPPING" << std::endl;
      exit(EXIT_FAILURE);
    }
// If basic tests are passed, create the file
    std::ofstream main_file(main_outfile.str().c_str(), std::ios::out);
    if(main_file.is_open()){
      main_file << "This is reSolve main out file" << std::endl;
    }
// Failure condition: stops.
    else{
      std::cout << "Problem in accessing working directory \"" << input_1.workdir << "\" (probably does not exist?): STOPPING" << std::endl;
      exit(EXIT_FAILURE);
    }

// If multi-machine and iteration end, combine previous k_vegas iterations and exit.
    if(input_1.multi_machine==1 && input_1.machine_tag == "END_ITER"){
      k_vegas_combiner(input_1);
      outcode = 1;
    }

// Initial definition of the structured variable for event saving
    input_1.event_info.save_events = input_1.save_events;
    input_1.event_info.integrator_flag = input_1.integrator_flag;
    input_1.event_info.workdir = input_1.workdir;
    input_1.event_info.multi_machine = input_1.multi_machine;
    input_1.event_info.machine_tag = input_1.machine_tag;
  }


  return outcode;

};
