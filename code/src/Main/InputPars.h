#ifndef _InputParsH_
#define _InputParsH_

#include <string>
#include <map>
#include <vector>

#include "events_out.h"


/* Here the InputPars structure is defined; it contains all input parameters
(which are INTENDED to be constant throughout the execution, even if not actively enforced to be const).
Input are either read from a file or assigned a default/standard value.
TRY TO KEEP THEM ORDERED!
*/

struct InputPars{

/* Define 'input_status': a MAP object containing input status during reading.
  Key:
    - tag absent: the program did not attempt to read the parameter yet
    - value -1: the program attempted to read at least once, did not match the tag,
      and found no default value (not a problem unless it still is =-1 at the end of the in
      file!)
    - value  1: the program attempted to read at least once, did not match the tag,
      but found and assigned a default value
    - value  0: a valid value was read for the parameter (can in principle be overwritten
      if the same tag is present twice in the in file by mistake)
*/
  std::map<std::string, int> input_status;

// Containers (MAP objects) for the default values (one per data type)
  std::map<std::string, const int> default_ints;
  std::map<std::string, const double> default_reals;
  std::map<std::string, const std::string> default_strings;


// Record the name of the input file

  std::string filename_0;

// Basic parameters: process flag, mu_R, mu_F, total energy ...

// process, order, initial state hadrons
  int process; // Process selector - each process should have a unique integer id, 1 = diphoton, 2 = drellyan
  int order; //flag for calculation order, flag should be -1 for LO with resummmation turned off, and 0, 1, 2 for LL, NLL and NNLL
  int resum_flag; // Resummation turned on or off

// CM energy
  double CM_energy; // kept in GeV throughout the code

// Verbosity
  int verbosity; //0=output nothing but answer, 1=output key test variables (i.e. resu, hard factors, etc), 2=output some more important couts, >=10 = important variables evaluated once per phase space point, >=12 less important variables evaluated once per phase space point, >=13 even less important variables evaluated once per phase space point, >=15 important variables evaluated once per b space point, >=16 less important variables evaluated/used per b space point, >=20 important variables evaluated once per mellin space point, >=25 less important variables used/evaluated once per mellin space point, >=50 is all couts I ever used for testing(LOADS!)


// PARALLELIZATION: Multi-machine runs -- UNDER TESTING NOW
  int multi_machine;
  std::string machine_tag;

// INTEGRATOR-RELATED STUFF
  int integrator_flag; // 1 = k_vegas, 2 = CUBA, 0 = read randoms in one by one
  int save_events; // 0 = No, 1 = "easy" version, 2 = lhe version (not available in cuba for now)
  std::string statefile; // if nonempty, state of integration is stored to "statefile" (binary file, allowing future resuming of integration)
  int resume_integration;

// Needed for option: integrator_flag=0
  int no_eventsin; //number of events read in from external file if events_flag is set to 1.
  std::string randoms_file;

// CUBA (Vegas)-specific inputs (see CUBA library documentation)
  int nvec; // Used for vectorization
  double epsrel, epsabs; // relative and absolute target uncertainties
  int VegFlags, cuba_seed; // various options + integer seed
  int mineval, maxeval; // min and approximate max # of evaluations
  int nstart, nincrease, nbatch; // # of evaluation at first run,
                                 // # of increases at each run and
                                 // # of "batch" evaluations
  int gridno; // 0: grid is deleted at end of run;
              // 1-10: grid is stored and the number serves as an identifier

// For k_vegas
  int seed;

// Directory where output events are stored
  std::string  workdir;

// histograms
  std::vector<std::string> histo_data;
  int hist_only;

// PDF fit-only mode
  int pdf_fitonly;

// events
  event_dumper_info event_info;

// METHODS

// Basic initialization function
  void ReadInput(int argc, char* argv[]);


};

#endif
