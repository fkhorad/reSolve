#ifndef _lhe_readerH_
#define _lhe_readerH_

#include <vector>
#include <string>

class PSpoint;


// LHE and pseudo-LHE

typedef std::vector<std::vector<double> > lhe_event;

int read_lhe_event(std::ifstream& infile, lhe_event& event);
void lhe_to_PS(lhe_event event, PSpoint& PS_);


// The "easy" format

typedef std::vector<std::vector<double> > easy_event;

int read_easy_event(std::string line, std::ifstream& infile, easy_event& event);
void easy_to_PS(easy_event event, PSpoint& PS_);


// Generic

void line_to_vec(std::string line, std::vector<double>& vec);



#endif
