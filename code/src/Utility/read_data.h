#ifndef _read_dataH_
#define _read_dataH_

#include <iostream>
#include <sstream>
#include <map>
#include <cstdlib>


/*
Header for Template functions for data reading/setting. Details:
- read_data: parses one line of a file stream looking for data of type 'type'
tagged by a std::string 'tag'. Also has the additional functionality of checking
if a default value exists and assigning that unless it is overridden. It also prints
to std::out every value as it is succesfully read.
- set_def: set a new default value for a data of type 'type' tagged by a
std::string 'tag'.
*/


template <typename type>
void read_data(std::map<std::string, int>& input_status,
               std::map<std::string, const type>& defaults,
               std::string line, std::string tag, type& value){

// Check if it is the first time that the program looks for this parameter
  std::map<std::string, int>::iterator it1;

  it1 = input_status.find(tag);
  if(it1 == input_status.end()){
// If so, record the first attempt to read the parameter
    input_status.insert( std::pair<std::string, int>(tag, -1) );
// and assign the default value, if any, to the parameter
    if(defaults.find(tag) != defaults.end()){
      value = defaults[tag];
      input_status[tag] = 1;
    }
  }

// Check if parameter was already read
  if(input_status[tag] == 0) return;
  else{
// If not, try to read the parameter from the input line which is being read
    std::stringstream linestream;
    std::size_t pos;

    std::stringstream tagP;
    tagP << tag << ":";
    pos = line.find(tagP.str());
    if(pos!=std::string::npos){
      linestream << line.substr(pos + tag.length() + 1);
      linestream >> value;
// Check on possible reading errors
      if(linestream.fail()){
        std::cout << "Error while reading parameter" << tag << "\n";
        exit(EXIT_FAILURE);
      }
// Record that the parameter was actually read and dump the value
      input_status[tag] = 0;
      std::cout << tag << ": " << value << "\n";
    }
  }
}


// Similar to read_data, but sets default values
template <typename type>
void set_def( std::map<std::string, const type>& defaults,
              std::string tag, type value){

  if(defaults.find(tag) == defaults.end())
    defaults.insert(std::pair<std::string, const type>(tag, value));

}


// Auxiliary function:
// dumps the values of any parameters which still have their default values at the end
// of reading. If any par has status '-1' - which means that the code attempted to read it,
// but found no value, either default or user-defined - exit with an error.
void dump_default_parameters(std::map<std::string, int> input_status, std::map<std::string, const int> default_ints, std::map<std::string, const double> default_reals, std::map<std::string, const std::string> default_strings);

#endif
