#ifndef _events_to_histogramsH_
#define _events_to_histogramsH_

#include <string>
#include <vector>
class histogram;
class PSpoint;


int events_to_histograms(std::string workdir, std::string file_name, int event_type, std::vector<histogram>& all_histos);

int events_to_histos_multi_iter(int iter_min, int iter_max, std::string workdir, std::string file_prefix, int event_type, std::vector<histogram>& all_histos);

void update_histogram(histogram& histo, const PSpoint& PS_, double weight);


#endif
