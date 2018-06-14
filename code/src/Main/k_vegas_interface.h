#ifndef _k_vegas_interfaceH_
#define _k_vegas_interfaceH_

#include "k_vegas.h"
struct InputPars;
#include <string>
#include <vector>
#include <fstream>


void k_vegas_call(const InputPars& input_1, int ndim, k_integrand_t integrand, void* userdata);

void k_vegas_combiner(InputPars& input_1);

void debugger(std::string inputname, int eventsets, int ndim, k_integrand_t integrand, void* userdata);

void combine_reweight(std::string part_filename, std::string filename, int event_type, int old_nevent, int new_nevent);
void reweight_easy_event(std::string line, std::ifstream& part_event_file, std::vector<std::string>& event, int old_nevent, int new_nevent);
void reweight_lhe_event(std::string line, std::ifstream& part_event_file, std::vector<std::string>& event, int old_nevent, int new_nevent);

#endif
