#ifndef _InputParsH_
#define _InputParsH_

#include <string>
#include <map>


/* Here the InputPars structure is defined; it contains all input parameters
(which are INTENDED to be constant throughout the execution, even if not actively enforced to be const).
Input are either read from a file or assigned a default/standard value.
TRY TO KEEP THEM ORDERED!
*/

class InputPars{

public:


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
  int process; // Process selector - each process should have a unique integer id
  int order; //flag for calculation order, flag should be -1 for LO with resummmation turned off, and 0, 1, 2 for LL, NLL and NNLL
  int ih1, ih2; //flag for changing beam1 and beam2 initial hadron from proton (1) to antiproton (-1)
  int resum_flag; // Resummation turned on or off

// Renormalization & factorization scales
  int muR_flag, muF_flag; // If = 1, they set mu_F, mu_R to the final state invariant mass
  double mu_R, mu_F, mu_min; // mu_R, mu_F right now have a double meaning:
// if they are read as a constant value for mu_X. If muF_flag=1, they are MULTIPLIERS, that is
// the actual scale is M*mu_R, with M the final state invariant mass.
// mu_min is the mininum value used for mu_F - needed to avoid troublesome zeroes.

// CM energy
  double CM_energy; // kept in GeV throughout the code

// Verbosity
  int verbosity; //0=output nothing but answer, 1=output key test variables (i.e. resu, hard factors, etc), 2=output some more important couts, >=10 = important variables evaluated once per phase space point, >=12 less important variables evaluated once per phase space point, >=13 even less important variables evaluated once per phase space point, >=15 important variables evaluated once per b space point, >=16 less important variables evaluated/used per b space point, >=20 important variables evaluated once per mellin space point, >=25 less important variables used/evaluated once per mellin space point, >=50 is all couts I ever used for testing(LOADS!)


// Resummation-specific stuff that probably should be moved elsewhere
  double mu_S; // Resummation scale
  double en_sec_multiplier; // Used for old-style PDF fit
  double QQ_Min, QQ_Max; //qmin and qmax for integration
  double QT_Min, QT_Max; //qtmin and qtmax for integration
  double eta_Min, eta_Max; //etamin and etamax for integration
  double ggnp, gqnp; //factors to account for Non-perturbative corrections to sudakovs, to factor in uncertainty from low qT of order lambdaQCD or less
// For De-Quadrature (intde2.cc)
  int lenaw;
  double tiny;
  double de_eps;
// If a specific PDF fit file is needed
  std::string pdffit_file;
  double pdf_fussiness;
  int pdf_fitonly;

// PARALLELIZATION: Multi-machine runs -- UNDER TESTING NOW
  int multi_machine;
  std::string machine_tag;

// PDFs
// This part is to be changed to make it LHA-compliant - we can discuss it later.
  int pdf_flag; // an identifier for the pdf set, but currently only one is present
  /* std::string prefix; // prefix for the pdf grid file */

// Some SM basic constants part 1
  int Nc, Nf;
  double alpha_QED;
  double gevm2topb;
  double M_p;

// INTEGRATOR-RELATED STUFF

  int integrator_flag; // 1 = k_vegas, 2 = CUBA, 0 = read randoms in one by one
  int save_events; // 0 = No, 1 = "easy" version, 2 = lhe version (not available in cuba for now)

// Needed for option: integrator_flag=0
  int no_eventsin; //number of events read in from external file if events_flag is set to 1.
  std::string randoms_file;


// CUBA (Vegas)-specific inputs (see CUBA library documentation)
  int ndim, ncomp; // Number of randoms, Number if integrand functions
  int nvec; // Used for vectorization
  double epsrel, epsabs; // relative and absolute target uncertainties
  int VegFlags, cuba_seed; // various options + integer seed
  int mineval, maxeval; // min and approximate max # of evaluations
  int nstart, nincrease, nbatch; // # of evaluation at first run,
                                 // # of increases at each run and
                                 // # of "batch" evaluations
  int gridno; // 0: grid is deleted at end of run;
              // 1-10: grid is stored and the number serves as an identifier
  std::string statefile; // if nonempty, state of integration is stored to "statefile" (binary file, allowing future resuming of integration)


// For k_vegas
  int seed;

// Directory where output events are stored
  std::string  workdir;


// LAST:
// Generic pointer to additional process-dependent input
//  user_input* user_1;


// METHODS

// Basic initialization function
  void basic_input(int argc, char* argv[]);

// Dump parameter which were assigned default values to std::out
  void dump_default_parameters();


};


#endif
