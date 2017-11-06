#ifndef _USERH_
#define _USERH_

#include "InputPars.h"

class c_mstwpdf;
class ResummationInfo;
class diphoton_input;

// This is where basic user-defined, or process dependent, definitions are made.
// Functions defined here are actually meant to be WRAPPERS: the low-level definitions
// should be in additional files. The idea is that this (along with USER.cc) should be
// the only "core" file one needs to modify when adding new processes - all other process-
// dependent definitions should go in additional files. See how the diphoton example works

// NOTE that the basic structure for the process-dependent definitions follows the same
// pattern as the main program: read input, do any pre-processing, do the main Monte Carlo
// calculation, do any post-processing, exit.


// Parameters collection

class all_info{

public:

  InputPars input_1;

// pointer needed for the MSTW pdf set
  c_mstwpdf* pdf;

// Process-dependent additional objects
  ResummationInfo* res_1;
  diphoton_input* diph_1;

};


// INPUT

// Wrapper function for process-specific input reading. Ultimately, this will just
// redirected to another process-specific, user-defined function provided in an
// additional file.
//

void ReadInput_user(int argc, char* argv[], all_info& info);


// PRE-PROCESSING WRAPPER

void PreProc_user(all_info& info);


// MAIN CALCULATION WRAPPER

void MonteCarlo_user(all_info& info);


// POST-PROCESSING WRAPPER

void PostProc_user(all_info&);

//

//////////////////////////////////
// INTERFACE TO DIPHOTON PROCESS
void interface_diphoton(all_info& info);

///////////////////////////////////////////
// INTERFACE TO QT-RESUMMATION PROCESSES
void interface_resummation(all_info& info);


#endif
