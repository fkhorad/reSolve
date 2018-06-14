
#include "read_data.h"

#include <string>

// Auxiliary function:
// dumps the values of any parameters which still have their default values at the end
// of reading. If any par has status '-1' - which means that the code attempted to read it,
// but found no value, either default or user-defined - exit with an error.

void dump_default_parameters(std::map<std::string, int> input_status, std::map<std::string, const int> default_ints, std::map<std::string, const double> default_reals, std::map<std::string, const std::string> default_strings){

  for (std::map<std::string, int>::iterator it=input_status.begin();
       it!=input_status.end(); it++){

    std::string tag;
    std::map<std::string, const int>::iterator it1;
    std::map<std::string, const double>::iterator it2;
    std::map<std::string, const std::string>::iterator it3;

    tag = it->first;
    if(it->second == 1){
      std::cout << tag << ": ";
      it1 = default_ints.find(tag);
      it2 = default_reals.find(tag);
      it3 = default_strings.find(tag);
      if(it1 != default_ints.end()) std::cout << default_ints[tag];
      if(it2 != default_reals.end()) std::cout << default_reals[tag];
      if(it3 != default_strings.end()) std::cout << default_strings[tag];
      std::cout << " (default value)\n";
    }
    if(it->second == -1){
      std::cout << "Parameter "<< tag << " missing in input file and no default value found: STOPPING" << "\n";
      exit(EXIT_FAILURE);
    }
  }
  std::cout << "\n";

}
