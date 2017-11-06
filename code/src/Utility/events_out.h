#ifndef _events_outH_
#define _events_outH_

#include <iostream>
#include <sstream>
#include <cmath>
#include <unistd.h>

#include "phase_space.h"


void dumper_easy(std::string filename, double f, const double weight, PSpoint* PS_, int dim, const double x[]);

void dumper_lhe(std::string filename, double f_weight, int n_particles, double* mass,
                int* PDG_num, int* inout, int* mother1, int* mother2,
                double mu_R, double alpha_QED, double alpha_s, PSpoint* PS_);

void printRowPDG(int x);
void printRowmom(double x);


std::string get_pids_filename(int pid_in, std::string workdir, std::string machine_tag,
                              int iter, int save_events);


#endif