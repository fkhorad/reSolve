
#include "AP_LO.h"

#include <cmath>

double Pqg(double y) {// LO qg splitting function (with as/pi normalisation as in Peskin and Schroeder)

  //  return 0.5*(y*y+(1-y)*(1-y));
  return 0.25*(y*y+(1-y)*(1-y));

}

double Pqqint(double y) {//Integral from 0 to y of the regular part of Pqq (again this is the LO qq splitting function with as/pi normalisation as in Peskin and Schroder), i.e. it is 4/3*\int_0^y ((1+y^2)/(1-y))

  //  return -4.0/3.0*(2*std::log(1-y)+y*y/2.0+y);
  return -2.0/3.0*(2*std::log(1-y)+y*y/2.0+y);

}

double Pqqreg(double y) {// LO qq splitting function regular part (i.e. with plus distribution already dealt with, with as/pi normalisation as in Peskin and Schroder)

  //return 4.0/3.0*(1+y*y)/(1-y);
  return 2.0/3.0*(1+y*y)/(1-y);

}

double D0int(double y) {//Integral of 1/(1-y) from 0 to y as occurs after rearraging and algebra of Pgg pieces, could move this into the CT_NLO code easily, just here for consistency

  // return std::log(1-y);
  return 0.5*std::log(1-y);

}

double Pggreg(double y) {//Regular part of the LO Pgg splitting function with an extra -1 piece (this is generated as the usual Pgg regular part of (1-y)/y + y*(1-y) has a -1+1 added with the -1 giving the Pggreg and the +1 altering the form of the 1/(1-zi)*pdfs term in the Sigma12

  // return 6*((1-2*y)/y + y*(1-y));
  return 0.5*6*((1-2*y)/y + y*(1-y));
 
}

double Pgq(double y) {// LO gq splitting function (with as/pi normalisation as in Peskin and Schroeder)

  // return 4.0/3.0*((1+(1-y)*(1-y))/y);
  return 2.0/3.0*((1+(1-y)*(1-y))/y);
}
