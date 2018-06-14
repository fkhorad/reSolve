#ifndef _observablesH_
#define _observablesH_

#include <string>
class PSpoint;


int obs_values(std::string obs, const PSpoint& PS_, double& value);

double qT_obs(const PSpoint&);
double qq_obs(const PSpoint&);
double eta_obs(const PSpoint&);
double mT_obs(const PSpoint&);
double pTmin_obs(const PSpoint&);
double pTmax_obs(const PSpoint&);

#endif
