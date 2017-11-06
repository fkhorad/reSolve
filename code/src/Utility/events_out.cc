
#include "events_out.h"


// OUTPUT: SAVE EVENTS TO FILE

/*
  std::stringstream filename;
  filename << workdir << "events_" << iter << "_" << machine_tag << ".dat";
  std::ofstream event_file(filename.str().c_str(), ios::out | ios::app);
*/

// V1: naive version

void dumper_easy(std::string filename, double f, const double weight, PSpoint* PS_, int ndim, const double x[]){

  std::ofstream event_file(filename.c_str(), std::ios::out | std::ios::app);
  
// New version with momenta
  int PS_size = PS_->dim();
  for(int ii=0; ii<PS_size; ii++) LorPrint(PS_->mom(ii), event_file);

// Old version with randoms

  for(int ii=0; ii<ndim; ii++){
    event_file << x[ii] << "  ";
    if(ii==ndim-1) event_file << std::endl;
  }

// Function value, MC weight & close

  event_file << f << "  " << weight << "\n\n";
  event_file.close();

}


// V2: PSEUDO-LHE format

void dumper_lhe(std::string filename, double f_weight, int n_particles, double* mass,
                int* PDG_num, int* inout, int* mother1, int* mother2,
                double mu_R, double alpha_QED, double alpha_s, PSpoint* PS_){

  std::ofstream event_file(filename.c_str(), std::ios::out | std::ios::app);


  event_file << "<event>" << "\n";
  event_file << " " << n_particles << "   " << 0 << "  " << f_weight << "  " << mu_R << "  "
    << alpha_QED << "  " << alpha_s << "\n";

  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(event_file.rdbuf()); //redirect std::cout to out.txt
  std::streamsize ss = std::cout.precision();
  for(int ii=0; ii<n_particles; ii++){
    printRowPDG(PDG_num[ii]);
    std::string filler = "   ";
    if(inout[ii]<0) filler = "  ";
    std::cout << filler << inout[ii] << "    " << mother1[ii] << "    " << mother2[ii] << "    " << 0 << "    " << 0;
    printRowmom(PS_->mom(ii)[1]); printRowmom(PS_->mom(ii)[2]); printRowmom(PS_->mom(ii)[3]); printRowmom(PS_->mom(ii)[0]);
    std::cout << "  " << mass[ii]; std::cout << std::fixed; std::cout.precision(1);
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



// THIS IS UNIX/LINUX ONLY (BUT SO IS CUBA, MOSTLY)
std::string get_pids_filename(int pid_in, std::string workdir, std::string machine_tag, int iter, int save_events){

  int pid;
  if(pid_in == 0) pid = getpid();
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
  if(save_events==1)
    filename << workdir << "events_" << iter << filler << machine_tag << "_" << pid << ".dat";
  else if(save_events==2)
    filename << workdir << "events_lhe_" << iter << filler << machine_tag << "_" << pid << ".lhe";
  else
    filename << "";
    
  return filename.str();

}


