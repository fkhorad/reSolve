#ifndef _events_outH_
#define _events_outH_

#include <sstream>
#include <vector>

class PSpoint;

struct event_dumper_info{
  int integrator_flag, save_events, multi_machine;
  std::string workdir, machine_tag;
//
  int n_particles;
  std::vector<double> mass;
  std::vector<int> PDG_num, inout, mother1, mother2;

};

void events_out(const event_dumper_info& info, int iter, double f_value, double weight, double mu_R, double alpha_QED, double alpha_s, const PSpoint& PS_, int dim, const double x[]);

void dumper_easy(std::string filename, const event_dumper_info& info, double f, double weight, const PSpoint& PS_, int ndim, const double x[]);

void dumper_lhe(std::string filename, const event_dumper_info&, double f_weight, double mu_R, double alpha_QED, double alpha_s, const PSpoint& PS_);

void printRowPDG(int x);
void printRowmom(double x);


std::string get_event_filename(int event_type, std::string workdir, std::string machine_tag, int iter);

std::string get_pids_filename(int pid_in, int event_type, std::string workdir, std::string machine_tag, int iter);



#endif
