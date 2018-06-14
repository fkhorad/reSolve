
#include "events_out.h"

#include <iostream>
#include <cmath>
#include <unistd.h> // nonstandard include for CUBA

#include "phase_space.h"


// OUTPUT: SAVE EVENTS TO FILE


void events_out(const event_dumper_info& info, int iter, double f_value, double weight, double mu_R, double alpha_QED, double alpha_s, const PSpoint& PS_, int dim, const double x[]){

  std::string filename;
  if(info.integrator_flag == 0 || info.integrator_flag == 1){
    filename = get_event_filename(info.save_events, info.workdir, info.machine_tag, iter);
  }
  else if(info.integrator_flag == 2){
    filename = get_pids_filename(0, info.save_events, info.workdir, info.machine_tag, iter);
  }
  else{}

  if(info.save_events==1){
    dumper_easy(filename, info, f_value, weight, PS_, dim, x);
  }
  else if(info.save_events==2){
    dumper_lhe(filename, info, f_value*weight, mu_R, alpha_QED, alpha_s, PS_);
  }

}


// V1: naive version

void dumper_easy(std::string filename, const event_dumper_info& info, double f_value, double weight, const PSpoint& PS_, int ndim, const double x[]){

  std::ofstream event_file(filename.c_str(), std::ios::out | std::ios::app);

// New version with momenta
  int PS_size = PS_.dim();
  for(int ii=0; ii<PS_size; ii++) LorPrint(PS_.mom(ii), event_file);

// Old version with randoms

  for(int ii=0; ii<ndim; ii++){
    event_file << x[ii] << "  ";
    if(ii==ndim-1) event_file << std::endl;
  }

// Function value, MC weight & close

  event_file << f_value << "  " << weight << "\n\n";
  event_file.close();

}


// V2: QUASI-LHE format

void dumper_lhe(std::string filename, const event_dumper_info& info, double f_weight, double mu_R, double alpha_QED, double alpha_s, const PSpoint& PS_){

  std::ofstream event_file(filename.c_str(), std::ios::out | std::ios::app);

  event_file << "<event>" << "\n";
  event_file << " " << info.n_particles << "   " << 0 << "  " << f_weight << "  " << mu_R << "  "
    << alpha_QED << "  " << alpha_s << "\n";

  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(event_file.rdbuf()); //redirect std::cout to out.txt
  std::streamsize ss = std::cout.precision();
  for(int ii=0; ii<info.n_particles; ii++){
    printRowPDG(info.PDG_num[ii]);
    std::string filler = "   ";
    if(info.inout[ii]<0) filler = "  ";
    std::cout << filler << info.inout[ii] << "    " << info.mother1[ii] << "    " << info.mother2[ii] << "    " << 0 << "    " << 0;
    printRowmom(PS_.mom(ii)[1]); printRowmom(PS_.mom(ii)[2]); printRowmom(PS_.mom(ii)[3]); printRowmom(PS_.mom(ii)[0]);
    std::cout << "  " << info.mass[ii]; std::cout << std::fixed; std::cout.precision(1);
    std::cout << "  " << 0. << "  " << 0. << "\n"; std::cout.precision(ss);
  }
  std::cout.rdbuf(coutbuf);

  event_file << "</event>" << "\n";
  event_file.close();

}

void printRowPDG(int x) {

  /// make it return a character when you've worked out the equivalent of printf
  // inserts spaces so as to align columns for positive and negative numbers and numbers with different numbers of digits

  int xa;
  if (x >= 0){
    xa = x;
    std::cout << " ";
  }
  else xa = -x;
  if (xa<10) std::cout << "       " << x;
  else if (xa<100) std::cout << "      " << x;
  else if (xa<1000) std::cout << "     " << x;
  else if (xa<10000) std::cout << "    " << x;
  else {}
}

void printRowmom(double x) {

  /// make it return a character when you've worked out the equivalent of printf
  // inserts spaces so as to align columns for positive and negative numbers and numbers with different numbers of digits
    std::streamsize ss = std::cout.precision();
    std::cout.precision(11);
    std::cout << std::scientific;
    if (x >= 0) std::cout << "  "  << x;
    else std::cout << " " << x;
    std::cout << std::fixed;
    std::cout.precision(ss);
}


std::string get_event_filename(int event_type, std::string workdir, std::string machine_tag, int iter){

  std::stringstream filename_stream;
  std::string filler = "";
  if(machine_tag != "") filler = "_";
  if(event_type == 1)
    filename_stream << workdir << "events_" << iter << filler << machine_tag << ".dat";
  if(event_type == 2)
    filename_stream << workdir << "events_lhe_" << iter << filler << machine_tag << ".lhe";

  return filename_stream.str();
}


// THIS IS UNIX/LINUX ONLY (BUT SO IS CUBA, MOSTLY)
std::string get_pids_filename(int pid_in, int event_type, std::string workdir, std::string machine_tag, int iter){

  int pid;
  if(pid_in == 0) pid = getpid(); // ONLY OS-DEPENDENT INSTRUCTION
  else pid = pid_in;

  int oldpid = 0;
  std::stringstream pids_filename;
  pids_filename << workdir << "pids.dat" ;

  std::ifstream pids_file_in(pids_filename.str().c_str());
  while( (oldpid != pid) and (pids_file_in.is_open()) and (!pids_file_in.eof()) ){
    pids_file_in >> oldpid;
  }
  pids_file_in.close();

  if(oldpid != pid){
    std::ofstream pids_file_out(pids_filename.str().c_str(), std::ios::out | std::ios::app);
    pids_file_out << pid << "\n";
    pids_file_out.close();
  }

  std::stringstream filename;
  std::string filler = "";
  if(machine_tag != "") filler = "_";
  if(event_type==1)
    filename << workdir << "events_" << iter << filler << machine_tag << "_" << pid << ".dat";
  else if(event_type==2)
    filename << workdir << "events_lhe_" << iter << filler << machine_tag << "_" << pid << ".lhe";
  else
    filename << "";

  return filename.str();

}
