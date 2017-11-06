//Routine to read in event file from old code

#include "eventsread_in.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>


void eventsread_in (int argc, std::string argv, std::vector<std::vector<double> >& xvec, int eventsets, int ndim)
{
    double drop = 0, drop1 = 0, drop2 = 0;

    std::string infilename;
    if(argc==1){
      infilename.assign("events.dat");
    }
    else{
      infilename.assign(argv);
    }

    xvec.resize(eventsets);
    for(int kk=0; kk<eventsets; kk++) xvec[kk].resize(ndim);


    std::ifstream myfile (infilename.c_str());
    if(myfile.fail()){
      std::cout << "Error opening input events file (maybe file missing or wrong dir?) - Stopping.\n";
      exit(EXIT_FAILURE);
    }
    std::string line, line1;
    if (myfile.is_open()){
      const int blanks =2; // # of useless lines before 1st set of randoms
      for(int i=0; i < blanks; i++)
        std::getline( myfile, line, '\n'); //blank line

      const int blanks2 =4; // # of useless lines between 2 sets of randoms (hopefully constant)     
      for (int j = 0; j < eventsets; j++) {
        std::getline( myfile, line, '\n'); //randoms line
        std::istringstream ss(line);
        ss >> xvec[j][0] >> xvec[j][1] >> xvec[j][2] >> xvec[j][3] >> xvec[j][4];
        std::cout << xvec[j][0] << " " << xvec[j][1] << " " << xvec[j][2] << " " << xvec[j][3] << " " << xvec[j][4] << std::endl;
        for(int i=0; i < blanks2; i++)
          std::getline( myfile, line, '\n'); //blank line
      }
    }

}
