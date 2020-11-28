
#include "USER.h"

#include <map>
#include <vector>
#include <sstream>

#include "cuba_interface.h"
#include "k_vegas_interface.h"
#include "constants.h"
#include "InputPars.h"

// Process-specific includes

// Generic Resummation
#include "resummation_input.h"
#include "resu_preproc.h"
#include "pdf_interface.h"
// Diphoton
#include "diphoton_input.h"
#include "diphoton_integrand.h"
// Drell-Yan
#include "drellyan_input.h"
#include "drellyan_integrand.h"


//extern std::ifstream debbie;
//std::ifstream debbie("debbie.dat");


// THE MAIN BODY: SELECTION, INITIALISATION AND EXECUTION OF PROCESSES
// USER-CUSTOMIZABLE

void MonteCarlo(const InputPars& info) {

//////////////////////
//  THE PROCESSES   //
//////////////////////

  if(info.pdf_fitonly==1){
    ResummationInfo resuminfo;
    resuminfo.fitonly = 1;
    resummation_input(info.filename_0, resuminfo);
    resu_preproc(info.event_info, resuminfo);
  }

// p(b)-p(b) --> Diphoton + X + qT resummation: process=1
  else if(info.process==1){
    diphoton_input data;
    diphoton_setup(info.filename_0, info.event_info, data);
//
    if(info.pdf_fitonly!=1){
      if(info.integrator_flag == 0){
// For debugging: instead of generating randoms, read them in from an events.dat file and evaluate the integrand function on those
        debugger(info.randoms_file, info.no_eventsin, data.ndim, diphoton_integrand, &data);
      }
      else if(info.integrator_flag == 1) {
        k_vegas_call(info, data.ndim, diphoton_integrand, &data);
      }
      else if(info.integrator_flag == 2) {
        integrand_t int_cuba = (integrand_t) diphoton_integrand_cuba;
        cuba_vegas_call(info, data.ndim, 1, int_cuba, &data);
      }
    }
  }

// Drell-Yan @ hadron (p or pb) collider + qT resummation: process=2
  else if(info.process==2){
    drellyan_input data;
    drellyan_setup(info.filename_0, info.event_info, data);
//
    if(info.pdf_fitonly!=1){
      if(info.integrator_flag == 0){
// For debugging: instead of generating randoms, read them in from an events.dat file and evaluate the integrand function on those
        debugger(info.randoms_file, info.no_eventsin, data.ndim, drellyan_integrand, &data);
      }
      else if(info.integrator_flag == 1) {
        k_vegas_call(info, data.ndim, drellyan_integrand, &data);
      }
      else if(info.integrator_flag == 2) {
        integrand_t int_cuba = (integrand_t) drellyan_integrand_cuba;
        cuba_vegas_call(info, data.ndim, 1, int_cuba, &data);
      }
    }

  }

  else{}


//  debbie.close();

}
