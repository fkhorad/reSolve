
#include "ADPINT.h"

#include <iostream>
#include <vector>
#include <cmath>

const int MAXINT = 500;
double FV[MAXINT],FU[MAXINT],U[MAXINT],V[MAXINT],RESULT[MAXINT],ERR[MAXINT];


void ADPINT(double A, double B, double z, int n, double& adp_res, double& adp_err, double (* F)(double, double, int, void*), void* userdata){

//c.....Integral of F(X) from A to B, with error
//c.....less than ABS(AERR) + ABS(RERR*INTEGRAL)
//c.....Best estimate of error returned in ERREST.


  std::vector< std::vector<double> > Cadpcal(MAXINT,std::vector<double>(2));
//     Work space:

  int NUMINT = 499;
  double DX = (B-A)/ NUMINT;
  for(int i = 1; i < NUMINT+1 ; i++){
    if (i == 1) {
      U[i] = A;
      FU[i] = F(U[i],z,n,userdata);
    }
    else {
      U[i] = V[i-1];
      FU[i] = FV[i-1];
    }
    if (i == NUMINT) {
      V[i] = B ;}
    else {
      V[i] = A + DX * i;}

    FV[i] = F(V[i],z,n,userdata);

    double WRE,WER;
    ADPCAL(i, U[i], V[i], FU[i], FV[i], z, n, WRE, WER, F, userdata);

    Cadpcal[i][0]= WRE;
    Cadpcal[i][1]= WER;
  }
//    FWe[i] = tmp.WFW();                                }

  double  ADPINT1 = 0.0;
  double   ERREST1 = 0.0;
  for(int i= 1; i< NUMINT + 1; i++){
    ADPINT1 += Cadpcal[i][0];
    ERREST1 += Cadpcal[i][1];
  }
  if((ERREST1/ADPINT1) > 0.01){
    std::cout << "Error greater than 1 per cent\n";
  }

  adp_res = ADPINT1;
  adp_err = ERREST1;

}



void ADPCAL(int ii,double Uii, double Vii, double FUii, double FVii,double z,int n,  double adp_res1, double adp_err1, double (* F)(double, double, int, void*), void* userdata){


//c.....Fill in details of interval I given endpoints

     double FW1 = F( (Uii + Vii) /2.0,z,n, userdata);
     double DX = Vii - Uii;
     double RESULT1 = DX * (FUii + 4.0 * FW1 + FVii) / 6.0;
     double ERR1 = fabs(DX * (FUii - 2.0 * FW1 + FVii )/ 12.0);

     adp_res1 = RESULT1;
     adp_err1 = ERR1;

}
