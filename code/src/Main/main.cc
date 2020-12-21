
#include <iostream>

#include "InputPars.h"
#include "PreProc.h"
#include "USER.h"
#include "PostProc.h"


int main(int argc, char* argv[]) {

/*
Read PROCESS-INDEPENDENT input (like basic functioning flags, process choice, parameters for the integrator, ...) from the input card.
*/
  InputPars input_basic;
  std::cout << "Reading input " << std::endl;
  int puppa = argc;
  
  input_basic.ReadInput(argc, argv);

// Process-independent setup & preprocessing
  int out1 = PreProc_basic(input_basic);

// Now call the main function, in USER.cc
  if(out1==0) MonteCarlo(input_basic);

// Process-independent setup & postprocessing
  PostProc_basic(input_basic);

  return 0;
}
