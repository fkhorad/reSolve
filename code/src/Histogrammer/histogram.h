#ifndef _histogramH_
#define _histogramH_

#include <vector>
#include <string>
#include <iostream>


class histogram{

public:

// The default constructor builds an histo with 10 evenly spaced bins of width 1.
  histogram();
// Serious constructor : needs a std::vector of bin edges plus an observable name as an input
  histogram(std::vector<double> bin_edges, std::string obs);
// Reset
  void reset();
  void set(std::vector<double> bin_edges, std::string obs);
// Public setter functions
  void add_to_histogram(double obs_value, double weight);
  void calc_sd();
  void weighted_avg(histogram old_histo);
// Getter functions
  int get_num_of_bins() const;
  double get_bin_value(int the_bin) const;
  double get_bin_value(double obs_value) const;
  double get_bin_sd(int the_bin) const;
  double get_bin_sd(double obs_value) const;
  std::string get_obs() const;
  int which_bin(double value) const;
  double max() const;
  double min() const;
  double integral() const;

  void histo_print(std::ostream& out1=std::cout);


private:

// Main data members
  int n_events;
  int n_bins;
  std::vector<double> bin_edges;
  std::vector<double> bin_values;
  std::vector<double> bin_squared_vals;
  std::vector<double> bin_sd;
  bool bins_are_sorted;
  std::string obs_name;

};


int read_histo_data(const std::string line, std::string& obs, std::vector<double>& bin_edges);


#endif
