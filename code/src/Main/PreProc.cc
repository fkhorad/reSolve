
#include "PreProc.h"

#include "InputPars.h"
#include <string>
#include <sstream>
#include <cstdio>
#include <complex>

#include "USER.h"


void PreProc_basic(all_info& info){

  InputPars input_1 = info.input_1;


// Parallelization check
  if(input_1.multi_machine==1 && input_1.machine_tag == ""){
    std::cout << "\"machine_tag\" variable must be non-empty for parallel running with k_vegas - STOPPING"
      << std::endl;
    exit(EXIT_FAILURE);
  }

  
// Try to create a main output file, at the same time checking if WORKDIR exists and if it's clean or already contains data.
  std::stringstream main_outfile;
  if(input_1.multi_machine==0){
    main_outfile << input_1.workdir << "reSolve_main_out.dat";
    std::ifstream main_file_test(main_outfile.str().c_str(), std::ios::in);
    if(!main_file_test.fail()){
      std::cout << "Main output file already existed - stopping to avoid overwriting data. Set a different working directory." << std::endl;
      exit(EXIT_FAILURE);
    }
    main_file_test.close();
  }
  else if(input_1.multi_machine==1){
    main_outfile << input_1.workdir << "reSolve_main_out_" << input_1.machine_tag << ".dat";
  }
  else{
    std::cout << "Unrecognized value for multi_machine flag " << input_1.multi_machine << " - stopping." << std::endl;
      exit(EXIT_FAILURE);    
  }
  std::ofstream main_file(main_outfile.str().c_str(), std::ios::out);
  if(main_file.is_open()){
    main_file << "This is reSolve main out file" << std::endl;
  }
  else{
    std::cout << "Problem in accessing working directory \"" << input_1.workdir << "\" (probably does not exist?): STOPPING" << std::endl;
    exit(EXIT_FAILURE);
  }



// FOR MSTW08 PDFs (default at least for now)

// Open the central PDF grid corresponding to the prefix
  char filename[100];
  if (input_1.pdf_flag == 80) {
      const char * prefix = "mstw2008lo";
      sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);
  }
  else if (input_1.pdf_flag == 81) {
      const char * prefix = "mstw2008nlo";
      sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);
  }
  else if (input_1.pdf_flag == 82) {
      const char * prefix = "mstw2008nnlo";
      sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);
  }
  // Declare a pdf object
  std::cout << "filename = " << filename << std::endl;
  info.pdf = new c_mstwpdf(filename);


// Call the Fortran initialisation routine with alpha_S(Q_0).
// FR2 and mTop in principle should be taken from input; Q0
// is essentially a dumb parameter (as long as alphaSQ0 is
// used along with it)
  double FR2 = 1., Q0 = 1., mTop = 1.e10;
  initalphas_(&info.pdf->alphaSorder,&FR2,&Q0,&info.pdf->alphaSQ0,&info.pdf->mCharm,&info.pdf->mBottom,&mTop);

// End MSTW08


};
