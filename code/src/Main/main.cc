
#include "main.h"

#include "USER.h"
#include "InputPars.h"
#include "PreProc.h"
#include "PostProc.h"


int main(int argc, char* argv[]) {

/*
  TO CHANGE (COMMENT):
  EVERY input parameter is recorded in a single structured variable, class
  "InputPars".
  The code very first instruction is defining a variable of InputPars class;
  the variable is then initialized by the function ReadInput, which
  reads data from an input file (or possibly more than one), and which includes
  a "generic" part plus a user-defined (process-dependent) one.
*/
  all_info info;
  ReadInput(argc, argv, info);
  /**/

/*
  PreProc: a routine that does any pre-processing, such as setting weights and points for deterministic integrations routines.
  In general, every calculation which does NOT depend on the phase space point
  should be done once and for all here.
  Again, this is split in a "generic" part plus a user-defined (process-dependent)
  one.
*/

  PreProc(info);


/**/

/*
  MonteCarlo_user: the main integration/histogramming routine, takes as an argument
  the InputPars structure.
  Optionally could return an additional output to be passed to PostProc - currently none.
  It is strictly a process-dependent function, declared and defined in USER.h and USER.cc
*/

  if (info.input_1.pdf_fitonly != 1) {
      MonteCarlo_user(info);
  }
/**/

/*
  PostProc: routine doing any necessary post-processing. Again, in principle there is a
  "generic" part plus a user-defined (process-dependent) one, though both are currently
  empty.
*/
  PostProc(info);

}





// MAIN FUNCTIONS


// INPUT

void ReadInput(int argc, char* argv[], all_info& info){

///////////////////////////////////////////////////////////////////////////////
// 1) Standard, process-independent part

// Now read "basic" (generic) input

  info.input_1.basic_input(argc, argv);


///////////////////////////////////////////////////////////////////////////////
// 2) User-defined (process-dependent) input (-> USER.cc)
  ReadInput_user(argc, argv, info);


///////////////////////////////////////////////////////////////////////////////
// 3) FINAL: a routine dumps all parameters that were not read, but assigned default values (the others have already been printed as they were read).
  info.input_1.dump_default_parameters();

}



// PREPROCESSING

void PreProc(all_info& info){

// Process-independent part
  PreProc_basic(info);
// ... and process-dependent part (->USER.cc)
  PreProc_user(info);

}


// POSTPROCESSING

void PostProc(all_info& info){

// Process-independent part
  PostProc_basic(info);
// ... and process-dependent part (->USER.cc)
  PostProc_user(info);

}
