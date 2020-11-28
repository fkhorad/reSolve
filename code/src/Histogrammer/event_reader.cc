
#include "event_reader.h"

#include <sstream>
#include <fstream>
#include "phase_space.h"


// LHE and pseudo-LHE

int read_lhe_event(std::ifstream& infile, lhe_event& event){

// Start reading event
// Read and check 1st line
      int outcode = 0;
      std::string line;
      std::getline(infile,line);
      if(infile.fail()) outcode = 1;
      int n_particles;
      event.resize(0);
      if(outcode == 0){
        event.push_back(std::vector<double>(0));
        line_to_vec(line, event[0]);
        if(event[0].size() != 6) outcode = 10;
        n_particles = (int) event[0][0];
        if(n_particles<1) outcode = 10;
      }

// Read and check the other lines
      if(outcode==0){
        for(int ii=0; ii<n_particles; ii++){
          std::getline(infile, line);
          if(infile.fail()) outcode = 2;
          if(outcode == 0){
            event.push_back(std::vector<double>(0));
            line_to_vec(line, event[ii+1]);
          }
          if(event[ii+1].size() != 13) outcode = 20;
          if(outcode > 0) break;
        }
      }

// Read and check last line
      if(outcode == 0){
        std::getline(infile, line);
        if(line != "</event>") outcode = 3;
      }
// End reading event

      return outcode;

}

void lhe_to_PS(lhe_event event, PSpoint& PS_){

  int dim = event.size();
  PS_.set_dim(dim-1);
  for(int ii=1; ii<dim; ii++){
    four_momentum mom(4);
    mom[0] = event[ii][9]; mom[1]=event[ii][6]; mom[2]=event[ii][7]; mom[3]=event[ii][8];
    PS_.set_mom(ii-1,mom);
  }
  PS_.set_products();

}


// The "easy" format

int read_easy_event(std::string line, std::ifstream& infile, easy_event& event){

// Start reading event
      int outcode = 0;
      event.resize(0);
// Read the lines
      int ii=0;
      do{
        if(outcode == 0){
          event.push_back(std::vector<double>(0));
          line_to_vec(line, event[ii]);
	  // std::cout << line << std::endl;
          ii++;
        }
        std::getline(infile, line);
        if(infile.fail() && !infile.eof()) outcode = 1;
        if(outcode > 0) break;
      } while(line!="");

      if(event.size()>=2){
	// std::cout << "event.size() = " << event.size() << std::endl;
// The next line should have -2 --> -1 if the randoms line is not present
        int n_particles = event.size()-2;
	// std::cout << "n_particles = " << n_particles << std::endl;
// Check the lines
        for(int ii=0; ii<n_particles; ii++){
          if(event[ii].size()!=4) outcode = 2;
        }
        if(event[n_particles+1].size()!=2){
	  // std::cout << "event[n_particles+1].size() = " << event[n_particles+1].size() << std::endl;
          outcode = 3;
        }
      }
// End reading event

      return outcode;

}

void easy_to_PS(easy_event event, PSpoint& PS_){

// The next line should have -2 --> -1 if the randoms line is not present
  int dim = event.size()-2;
  PS_.set_dim(dim);
  for(int ii=0; ii<dim; ii++){
    four_momentum mom = event[ii];
    PS_.set_mom(ii,mom);
  }
  PS_.set_products();

}


// Generic

void line_to_vec(std::string line, std::vector<double>& vec){

  std::stringstream linestr;
  linestr << line;
  do{
    double val;
    linestr >> val;
    vec.push_back(val);
  }while(!linestr.fail());
  vec.pop_back();

}
