
#include "histogram.h"

#include <algorithm>
#include <sstream>
#include <cmath>


// First and last bins reserved for underflow / overflow values
//Constructors: set empty histograms
histogram::histogram(){
  double bin_edges_arr[]= {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};
  std::vector<double> bin_edges_(bin_edges_arr, bin_edges_arr+9);
  this->set(bin_edges_, "qT");
}
histogram::histogram(std::vector<double> bin_edges_, std::string obs){
  this->set(bin_edges_, obs);
}
// Set: used by constructors too, empties the histo if pre-esisting
void histogram::set(std::vector<double> bin_edges_, std::string obs){

  n_events = 0;
  obs_name = obs;
  bin_edges = bin_edges_;
  n_bins = bin_edges.size()-1;
  bins_are_sorted=std::is_sorted(bin_edges.begin(),bin_edges.end());
  bin_values.resize(n_bins+2);
  bin_values.assign(n_bins+2, 0.);
  bin_squared_vals.resize(n_bins+2);
  bin_squared_vals.assign(n_bins+2, 0.);
  bin_sd.resize(n_bins+2);
  bin_sd.assign(n_bins+2, 0.);

}
// Reset
void histogram::reset(){
  for(int ii=0; ii<n_bins+1; ii++){
     bin_values.at(ii) = 0.;
     bin_squared_vals.at(ii) = 0.;
     bin_sd.at(ii) = 0.;
  }
  bin_values.at(n_bins+1) = 0.;
}


// Main setter function (public)
void histogram::add_to_histogram(double obs_value, double weight){
  
  n_events++;
  int ii = 0;
  ii = which_bin(obs_value);

  if(ii==0 or ii>n_bins){
    bin_values[ii] += weight;
    bin_squared_vals[ii] += weight*weight;
  }
  else{
    bin_values[ii] += weight/(bin_edges[ii]-bin_edges[ii-1]);
    bin_squared_vals[ii] += weight*weight/(bin_edges[ii]-bin_edges[ii-1])/(bin_edges[ii]-bin_edges[ii-1]);
  }
  
}

void histogram::calc_sd(){
  if(n_events>=2){
    for(int ii=0; ii<n_bins+2; ii++){
      bin_sd[ii] = std::sqrt( (n_events*bin_squared_vals[ii] - bin_values[ii]*bin_values[ii])/(n_events-1.) );
    }
  }
}

// Public const getter functions
int histogram::get_num_of_bins() const{
  int res = n_bins;
  return res;
}
double histogram::get_bin_value(int the_bin) const{
  double res = bin_values[the_bin];
  return res;
}
double histogram::get_bin_value(double obs_value) const{
  double res = bin_values[which_bin(obs_value)];
  return res;
}
double histogram::get_bin_sd(int the_bin) const{
  double res = bin_sd[the_bin];
  return res;
}
double histogram::get_bin_sd(double obs_value) const{
  double res = bin_sd[which_bin(obs_value)];
  return res;
}
std::string histogram::get_obs() const{
  return obs_name;
}
int histogram::which_bin(double obs_value) const{
  int ii = 0;
  while(obs_value > bin_edges[ii]){
    ii++;
    if(ii==n_bins+1) break;
  }
  return ii;
}

void histogram::histo_print(std::ostream& out1){
  // 0s here are for x error bars for gnuplot plotting
  out1 << "<" << bin_edges[0] << " " << bin_values[0] << " 0 " << bin_sd[0] << std::endl;
  for(int ii=0; ii<n_bins; ii++){
    out1 << (bin_edges[ii]+bin_edges[ii+1])/2. << " " << bin_values[ii+1] << " 0 " << bin_sd[ii+1] << std::endl;
  }
  out1 << ">" << bin_edges[n_bins] << " " << bin_values[n_bins+1] << " 0 " << bin_sd[n_bins+1] << std::endl;

}

double histogram::max() const{
  return bin_edges[n_bins+1];
}
double histogram::min() const{
  return bin_edges[0];
}
double histogram::integral() const{
  double res;
  res = bin_values[0];
  for(int ii=1; ii<n_bins+1; ii++)
    res+=bin_values[ii]*(bin_edges[ii]-bin_edges[ii-1]);
  res+=bin_values[n_bins+1];
  return res;
}



int read_histo_data(const std::string line, std::string& obs, std::vector<double>& bin_edges){

  int outcode = 0;

  int bin_int;
  double min, max;
  bin_edges.resize(0);
  std::stringstream linestr;
  linestr << line;
  linestr >> obs;
  linestr >> bin_int;
  if(linestr.fail()) outcode = 1;

  if(outcode==0){
    if(bin_int==0){
      do{
        double val;
        linestr >> val;
        bin_edges.push_back(val);
      }while(!linestr.fail());
      bin_edges.pop_back();
    }
    else{
      linestr >> min;
      linestr >> max;
      if(linestr.fail()) outcode = 2;
      if(outcode==0){
        for(int ii=0; ii<bin_int+1; ii++){
          bin_edges.push_back(min+ii*(max-min)/bin_int);
        }
      }
    }
  }

  return outcode;
}


void histogram::weighted_avg(histogram old_histo){

  double wgt1, wgt2, wgt_tot;

  if(old_histo.get_num_of_bins() == n_bins){
    for(int ii=0; ii<n_bins+2; ii++){
      if(bin_sd[ii]>0.)
        wgt1 = 1/(bin_sd[ii]*bin_sd[ii]);
      else
        wgt1 = 0.;
      if(old_histo.get_bin_sd(ii))
        wgt2 = 1/(old_histo.get_bin_sd(ii)*old_histo.get_bin_sd(ii));
      else
        wgt2 = 0.;
      wgt_tot = wgt1+wgt2;
      if(wgt_tot>0.){
        bin_values[ii] = (wgt1*bin_values[ii] + wgt2*old_histo.get_bin_value(ii))/wgt_tot;
        bin_sd[ii] = 1/std::sqrt(wgt_tot);
      }
      else{
        bin_values[ii] = 0.;
        bin_sd[ii] = 0.;
      }
    }
  }

}
