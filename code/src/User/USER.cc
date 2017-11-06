
#include "USER.h"

#include <map>
#include <vector>
#include <sstream>

#include "cuba_interface.h"
#include "k_vegas_interface.h"
#include "eventsread_in.h"

// Process-specific includes

// Resummation
#include "resu_preproc.h"
// Diphoton
#include "diphoton_input.h"
#include "diphoton_integrand.h"


// USER INPUT WRAPPER FUNCTION

void ReadInput_user(int argc, char* argv[], all_info& info){

// The idea is to add an IF block similar to the default diphoton one for each process
// or generically "mode of use" one wants to define. Each process should have an unique
// integer tag


  if(info.input_1.resum_flag == 1){

    interface_resummation(info);

    if(info.input_1.process==1 || info.input_1.process==-1){
// Diphoton (with resummation): process=1
      interface_diphoton(info);
      ReadInput_diphoton(info.diph_1,info.input_1.input_status,
        info.input_1.default_ints,info.input_1.default_reals);
    }
  }

}


// USER PREPROCESSING WRAPPER FUNCTION

void PreProc_user(all_info& info) {

// qT resummation
  if(info.input_1.resum_flag == 1){

      resu_preproc(info.res_1, info.pdf, info.input_1.process);

  }
  else{
// EMPTY FOR NOW
  }

}

// USER MAIN MONTE CARLO WRAPPER FUNCTION

void MonteCarlo_user(all_info& info) {

/*

 COMMENT TO BE UPDATED

*/

// Combine previous k_vegas iterations
  if(info.input_1.multi_machine == 1 &&
     info.input_1.machine_tag == "END_ITER"){
    k_vegas_combiner(info.input_1);
  }

// Diphoton (with resummation): process=1
  else if(info.input_1.process==1||info.input_1.process==-1){
//process = -1 option is just evaluate Born - for testing
    DiphRes data;
    data.diph_in = info.diph_1;
    data.res_1 = info.res_1;

    if(info.input_1.integrator_flag == 0){
// For debugging: instead of generating randoms, read them in from an events.dat file
// and evaluate the integrand function on those
      std::string inputname = info.input_1.randoms_file;
      int eventsets = info.input_1.no_eventsin;
      int ndim = info.input_1.ndim;

      std::vector<std::vector<double> > xvec;
      eventsread_in (0, inputname, xvec, eventsets, ndim);
      double res;
      for(int k = 0; k<eventsets; k++) {
         double xx[ndim];
         for(int jj=0; jj<ndim; jj++) xx[jj] = xvec[k][jj];
         diphoton_integrand(ndim, xx, res, &data, 1, 1.);
      }
    }
    else if(info.input_1.integrator_flag == 1) {
      k_vegas_call(info.input_1, diphoton_integrand, &data);
    }
    else if(info.input_1.integrator_flag == 2) {
      integrand_t int_cuba = (integrand_t) diphoton_integrand_cuba;
      cuba_vegas_call(info.input_1, int_cuba, &data);
//      cuba_vegas_call(info.input_1, diphoton_integrand_cuba, &data);
    }
  }

  else{}

}


// USER POSTPROCESSING WRAPPER FUNCTION

void PostProc_user(all_info& info){

// qT resummation
  if(info.input_1.resum_flag == 1){
    delete info.res_1;
  }
// Diphoton (with resummation): process=1
  if(info.input_1.process==1){
    delete info.diph_1;
  }
  else{}
}


//////////////////////////////////
// INTERFACE TO DIPHOTON PROCESS

void interface_diphoton(all_info& info){

// Construct a diphoton_input object
// (with 'new' because it needs to still be allocated when ReadInput exits)
// and connect it to the general parameter structure
  diphoton_input* diph_in = new diphoton_input;
  info.diph_1 = diph_in;

// Set MonteCarlo dimensions
  info.input_1.ndim = 5;

// Copy relevant information from main input to diphoton input structure

    info.diph_1->gevm2topb = info.input_1.gevm2topb;
    info.diph_1->CM_energy = info.input_1.CM_energy;
    info.diph_1->M_p = info.input_1.M_p;
    info.diph_1->alpha_QED = info.input_1.alpha_QED;
    info.diph_1->filename = info.input_1.filename_0;
    info.diph_1->workdir = info.input_1.workdir;
    info.diph_1->machine_tag = info.input_1.machine_tag;
    info.diph_1->verbosity = info.input_1.verbosity;
    info.diph_1->Nf = info.input_1.Nf;
    info.diph_1->Nc = info.input_1.Nc;
    info.diph_1->save_events = info.input_1.save_events;
    info.diph_1->integrator_flag = info.input_1.integrator_flag;

}


///////////////////////////////////////////
// INTERFACE TO QT-RESUMMATION PROCESSES

void interface_resummation(all_info& info){

// Construct a ResummationInfo object
// (with 'new' because it needs to still be allocated when PreProc exits)
// and connect it to the general parameter structure
    ResummationInfo* res_1 = new ResummationInfo;
    info.res_1 = res_1;

// Copy relevant values from general input to resuminfo
    info.res_1->QQ_Min = info.input_1.QQ_Min;
    info.res_1->QQ_Max = info.input_1.QQ_Max;
    info.res_1->QT_Min = info.input_1.QT_Min;
    info.res_1->QT_Max = info.input_1.QT_Max;
    info.res_1->eta_Min = info.input_1.eta_Min;
    info.res_1->eta_Max = info.input_1.eta_Max;
    info.res_1->mu_R = info.input_1.mu_R;
    info.res_1->mu_F = info.input_1.mu_F;
    info.res_1->mu_S = info.input_1.mu_S;
    info.res_1->verbosity = info.input_1.verbosity;
    info.res_1->order = info.input_1.order;
    info.res_1->ih1 = info.input_1.ih1;
    info.res_1->ih2 = info.input_1.ih2;
    info.res_1->ggnp = info.input_1.ggnp;
    info.res_1->gqnp = info.input_1.gqnp;
    info.res_1->Nf = info.input_1.Nf;

    info.res_1->muF_flag = info.input_1.muF_flag;
    info.res_1->muR_flag = info.input_1.muR_flag;
    info.res_1->lenaw = info.input_1.lenaw;
    info.res_1->tiny = info.input_1.tiny;
    info.res_1->de_eps = info.input_1.de_eps;
    info.res_1->en_sec_multiplier = info.input_1.en_sec_multiplier;
    info.res_1->CM_energy = info.input_1.CM_energy;
    info.res_1->mu_min = info.input_1.mu_min;
    info.res_1->pdf_flag = info.input_1.pdf_flag;
    // info.res_1->prefix = info.input_1.prefix;
    info.res_1->pdffit_file = info.input_1.pdffit_file;
    info.res_1->pdf_fussiness = info.input_1.pdf_fussiness;
    
    info.res_1-> multi_machine = info.input_1.multi_machine;
    info.res_1-> machine_tag = info.input_1.machine_tag;

    int Nc = info.input_1.Nc;
    info.res_1->Cf = (Nc*Nc-1.0)/(2.0*Nc);
    info.res_1->Ca = Nc;

    if(info.res_1->pdf_fussiness <= 0.){
      std::cout << "Parameter \"pdf_fussiness\" is " << info.res_1->pdf_fussiness
        << ", while it is supposed to be > 0. Using default value "
        << info.input_1.default_reals["pdf_fussiness"] << " instead" << std::endl;
      info.res_1->pdf_fussiness = info.input_1.default_reals["pdf_fussiness"];
    }
    if(info.res_1->en_sec_multiplier <= 1.){
      std::cout << "Parameter \"en_sec_multiplier\" is " << info.res_1->en_sec_multiplier
        << ", while it is supposed to be > 1. Using default value "
        << info.input_1.default_reals["en_sec_multiplier"] << " instead" << std::endl;
      info.res_1->en_sec_multiplier = info.input_1.default_reals["en_sec_multiplier"];
    }

}
