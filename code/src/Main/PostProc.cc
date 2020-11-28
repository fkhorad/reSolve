
#include "PostProc.h"

#include "InputPars.h"
#include "histogrammer.h"


void PostProc_basic(const InputPars& input_1){

// Do histograms
  if(input_1.save_events==1 || input_1.save_events==2 || input_1.hist_only == 1 || input_1.pdf_fitonly!=1){
    std::string file_prefix ("events");
    if(input_1.save_events==2) file_prefix = file_prefix.append("_lhe");
    histogrammer(input_1.histo_data, input_1.save_events, 1, file_prefix, input_1.workdir);
  }

};
