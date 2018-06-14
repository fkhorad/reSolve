
#include <fstream>
#include <vector>
#include <sstream>

#include "histogram.h"
#include "events_to_histograms.h"


int histogrammer(std::vector<std::string> histo_data, int event_type, int multi_iter, std::string file_prefix, std::string workdir){


  std::vector<histogram> all_histos(0);
  std::string obs;
  std::vector<double> bin_edges;
  for(unsigned int ii=0; ii<histo_data.size(); ii++){
    read_histo_data(histo_data[ii], obs, bin_edges);
    histogram histo(bin_edges, obs);
    all_histos.push_back(histo);
  }


// DEBUG output
  std::cout << "# of histos found: ";
  std::cout << all_histos.size() << std::endl << std::endl;

// Main routine

  if(all_histos.size()>0){

  if(multi_iter == 0){
    events_to_histograms(workdir, file_prefix, event_type, all_histos);
  }
  else if(multi_iter == 1){

// Check for presence of event files and # of iterations
    int iter_min=0, iter_max=0;
    const int max_iter_check = 1000;
    std::string file_ext;
    if(event_type==1) file_ext = "dat";
    if(event_type==2) file_ext = "lhe";
    for(int ii=1; ii<=max_iter_check; ii++){
      std::stringstream infile_stream;
      infile_stream << workdir << file_prefix << "_" << ii << "." << file_ext;
      std::ifstream in_file(infile_stream.str().c_str());
      if(in_file.is_open()){
        if(iter_min==0) iter_min=ii;
        iter_max=ii;
      }
      else if(iter_min>0) break;
    }
//
    if(iter_min>0){
      events_to_histos_multi_iter(iter_min, iter_max, workdir, file_prefix, event_type, all_histos);
    }
    else{
      std::cout << "NO EVENTS FILES FOUND in workdir " << workdir << "!" << std::endl;
    }
  }


// Output
  for(unsigned int jj=0; jj<all_histos.size(); jj++){
    std::cout << "Histo " << jj << " - " << all_histos[jj].get_obs() << std::endl;
    std::cout << "Integral: " << all_histos[jj].integral() << std::endl;
//    all_histos[jj].histo_print(std::cout);
    std::cout << std::endl;
  }

  std::cout << std::endl;

  }

  return 0;

}
