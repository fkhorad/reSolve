
#include "events_to_histograms.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include "phase_space.h"
#include "histogram.h"
#include "observables.h"
#include "event_reader.h"


int events_to_histograms(std::string workdir, std::string file_name, int event_type, std::vector<histogram>& all_histos){


// Open events file
  std::stringstream infilename_ev;
  infilename_ev << workdir << file_name;
  std::ifstream infile_ev(infilename_ev.str().c_str());
  if(infile_ev.fail()){
    std::cout << "Event file not found" << std::endl;
    exit(2);
  }

// Main loop: read event by event from file, update histograms
  int nevents = 0, uncut_events = 0;
  double tot_weight = 0., tot_weight2 = 0.;
  PSpoint PS_;
  std::string line;
  do{
    int outcode = 0;
    double weight = 0.;

    std::getline(infile_ev, line);
    if(event_type==1){
      if(line != ""){
        nevents++;
//        std::cout << "nevent  " << nevents << "\n";
        easy_event event(0);
        outcode = read_easy_event(line, infile_ev, event);
// If the event size is not at least 3, there can be no PS so the event is
// not well-formed, so assign weight 0.
        if(outcode==0 && event.size()>=3){
          unsigned int nn = event.size()-1;
          weight = event[nn][0]*event[nn][1];
          easy_to_PS(event, PS_);
        }
        else weight = 0.;
      }
    }
    else if(event_type==2){
      if(line == "<event>"){
        nevents++;
        lhe_event event(0);
        outcode = read_lhe_event(infile_ev, event);
        if(outcode==0){
          lhe_to_PS(event, PS_);
          weight = event[0][2];
        }
        else weight = 0.;
      }
    }

    tot_weight += weight;
    if(std::abs(weight)>1e-300 && outcode==0){
      uncut_events++;
      tot_weight2 = tot_weight2 + weight;
      for(unsigned int jj=0; jj<all_histos.size(); jj++){
        update_histogram(all_histos[jj], PS_, weight);
      }
    }
    else if(outcode!=0){
      std::cout << "Something fishy happened, outcode: " << outcode << std::endl;
      std::cout << "ignoring this event and trying to continue reading" << std::endl << std::endl;
    }
  }while(!infile_ev.eof());



// Finalize histograms by calculating uncertainties
  for(unsigned int jj=0; jj<all_histos.size(); jj++){
    all_histos[jj].calc_sd();
// Save histograms to files
    std::stringstream outfilename;
    outfilename << workdir << "histo_" << jj << "_" << all_histos[jj].get_obs() << ".dat";
    std::ofstream outfile(outfilename.str().c_str());
    all_histos[jj].histo_print(outfile);
  }

// DEBUG Output

  std::cout << std::endl;
  std::cout << "Total events: " << nevents << std::endl;
  std::cout << "Uncut events: " << uncut_events << std::endl << std::endl;
  std::cout << "Total weight: " << tot_weight << std::endl;
  std::cout << "Total weight x-check: " << tot_weight2 << std::endl;
  std::cout << std::endl << std::endl;

  return 0;

}

int events_to_histos_multi_iter(int iter_min, int iter_max, std::string workdir, std::string file_prefix, int event_type, std::vector<histogram>& all_histos){

  int outcode = 0;
  std::string file_ext;
  if(event_type==1) file_ext = ".dat";
  else if(event_type==2) file_ext = ".lhe";
  else outcode = -1;

  std::vector<histogram> all_histos_old;

  if(outcode==0){
    for(int kk=iter_min; kk<=iter_max; kk++){
      std::stringstream filename_stream;
      filename_stream << file_prefix << "_"  << kk << file_ext;
      std::string filename = filename_stream.str();
      std::cout << filename << std::endl;

      all_histos_old = all_histos;
      for(unsigned int mm=0; mm<all_histos.size(); mm++){
        all_histos[mm].reset();
      }

      outcode = events_to_histograms(workdir, filename, event_type, all_histos);
      if(kk>0){
        for(unsigned int mm=0; mm<all_histos.size(); mm++)
          all_histos[mm].weighted_avg(all_histos_old[mm]);
      }

    }
  }

  return outcode;

}


void update_histogram(histogram& histo, const PSpoint& PS_, double weight){
  double value = 0., modifier = 1.;
  obs_values(histo.get_obs(), PS_, value, modifier);
  if(weight==weight){  // Naive NaN-guard
    histo.add_to_histogram(value, weight*modifier);
  }
}
