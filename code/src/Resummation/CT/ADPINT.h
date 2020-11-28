#ifndef _ADPINTH_
#define _ADPINTH_

void ADPINT(double A, double B, double z, int n, double& adp_res, double& adp_err, double (* F)(double, double, int, void*), void* userdata);
void ADPCAL(int ii,double Uii, double Vii, double FUii, double FVii,double z,int n,  double adp_res1, double adp_err1, double (* F)(double, double, int, void*), void* userdata);


#endif
