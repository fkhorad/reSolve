
#include "Itilde.h"

#include "constants.h"


//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//c.....Function BK(n,z)
//c.....BK(n,z) is the n-derivative of BesselK[nu,z]
//c.....with respect to nu in nu=1

//c.....Itilde defined as in the paper

double Itilde(int m,double xmio){

  double argum,Itildee;
  double b0,logx;

  double z2, z3, Egamma;
  z2 = k_constants::zeta2;
  z3 = k_constants::zeta3;
  Egamma = k_constants::Eulerconst;


  b0=2.0*std::exp(-Egamma);

  argum=b0*xmio;
  logx=std::log(xmio);

  if(m == 1){
    Itildee=-zBK(0,argum)/(xmio*xmio);
    // std::cout << "argum = " << argum << std::endl;
    // std::cout << "zBK(0,0.2) = " << zBK(0,0.2) << std::endl;
  }
  else if (m == 2){
    Itildee=2.0/(xmio*xmio)*(zBK(0,argum)*logx-zBK(1,argum));
  }
  else if (m == 3){
    Itildee=-3.0/(xmio*xmio)*(zBK(0,argum)*(logx*logx-z2)
            -2.0*zBK(1,argum)*logx+zBK(2,argum));
  }
  else if (m == 4){
    Itildee=4.0/(xmio*xmio)*(zBK(0,argum)*(std::pow(logx,3)-3.0*z2*logx+2.0*z3)
              -3.0*zBK(1,argum)*(logx*logx-z2)+3.0*zBK(2,argum)*logx
              -zBK(3,argum));
  }

  return Itildee;

}





//C     n-derivative of the function BesselK[nu,z]
//C     with respect to nu for nu=1
//C     NOTE: IT IS MULTIPLIED by z
//C     NOTE: needs ADPINT


double zBK(int n,double z){

  double zbk;
  double max=10.0;

// Use power-series approximated form only for z<1.5

  if(z < 1.5){

    if(n == 0){
      zbk=zbesselk0(z);}
    else if(n == 1){
      zbk=zbesselk1(z);}
    else if(n == 2){
      zbk=zbesselk2(z);}
    else if(n == 3){
      zbk=zbesselk3(z);}
  }

// Use integral representation (doing integral numerically) for z>1.5

  else{
    double integral, err;
    void* data;
    ADPINT(0.0,max,z,n, integral,err, F0, data);
    zbk=z*integral;
  }

  return zbk;

}

//C     Approximated forms of BesselK[nu,z] for nu=1 and
//C     derivatives with respect to nu of BesselK[nu,z] at nu=1
//C     All functions multiplied by z

double zbesselk0(double z){

  double zm,loz,zbesselk00;
  double Egamma = k_constants::Eulerconst;

  zm=z/2.0;
  loz=std::log(zm);

  zbesselk00=1.0+z*zm*(loz-0.5*(1.0-2.0*Egamma))
             +std::pow(zm,4)*(loz-0.5*(2.5-2.0*Egamma))
             +std::pow(zm,6)/6.0*(loz-0.5*(10.0/3.0-2.0*Egamma))
             +std::pow(zm,8)/72.0*(loz-0.5*(47.0/12.0-2.0*Egamma))
             +std::pow(zm,10)/1440.0*(loz-0.5*(131.0/30.0-2.0*Egamma));

      return zbesselk00;
  }


double zbesselk1(double z){

  double zbesselk11,zm,loz;

  zm=z/2.0;
  loz=std::log(zm);
  double Egamma = k_constants::Eulerconst;

  zbesselk11=-(loz+Egamma)-zm*zm*(loz-1.0+Egamma)
             -0.25*std::pow(zm,4)*(loz-1.5+Egamma)
             -std::pow(zm,6)/36.0*(loz-11.0/6.0+Egamma)
             -std::pow(zm,8)/576.0*(loz-25.0/12.0+Egamma);


      return zbesselk11;
  }


double zbesselk2(double z){

  double a[14];
  double zbesselk22,loz,zm;
  a[0] = 1.15443132980306572;
  a[1] = 1.97811199065594511;
  a[2] = 0.154431329803065721;
  a[3] = 4.801792651508824500;
  a[4] = 0.806235643470665767;
  a[5] =-0.672784335098467139;
  a[6] = 3.285072828402112960;
  a[7] =-1.945338757678943440;
  a[8] =-0.181575166960855634;
  a[9] = 0.694195147571435559;
  a[10]=-0.607655744858515573;
  a[11]=-0.019182189839330562;
  a[12]= 0.068894530444636532;
  a[13]=-0.070514317816328185;

  zm=z/2;
  loz=std::log(zm);

  zbesselk22=loz*loz+a[0]*loz+a[1]
       +zm*zm*(2.0*std::pow(loz,3)/3.0+a[2]*loz*loz+a[3]*loz+a[4])
       +std::pow(zm,4)*(std::pow(loz,3)/3.0+a[5]*loz*loz+a[6]*loz+a[7])
       +std::pow(zm,6)*(std::pow(loz,3)/18.0+a[8]*loz*loz+a[9]*loz+a[10])
       +std::pow(zm,8)*(std::pow(loz,3)/216.0+a[11]*loz*loz+a[12]*loz+a[13]);

  return zbesselk22;
}


double zbesselk3(double z){

  double b[15];
  double zbesselk33,loz,zm;

  b[0] = 1.731646994704598580;
  b[1] = 5.934335971967835330;
  b[2] = 5.444874456485317730;
  b[3] =-1.268353005295401420 ;
  b[4] = 8.471041982558638170;
  b[5] =-3.026167526073320430;
  b[6] =-0.692088251323850355 ;
  b[7] = 2.809848746963509900;
  b[8] =-2.161466255000085060;
  b[9] =-0.104676472369316706;
  b[10]= 0.381989731242156681;
  b[11]=-0.367492827636283900;
  b[12]=-0.007844362856415627;
  b[13]= 0.027796539630842606;
  b[14]=-0.029917436634978395;

  zm=z/2.;
  loz=std::log(zm);

  zbesselk33=std::pow(loz,3)+b[0]*loz*loz+b[1]*loz+b[2]
        +zm*zm*(std::pow(loz,3)+b[3]*loz*loz+b[4]*loz+b[5])
        +std::pow(zm,4)*(std::pow(loz,3)/4.0+b[6]*loz*loz+b[7]*loz+b[8])
        +std::pow(zm,6)*(std::pow(loz,3)/36.0+b[9]*loz*loz+b[10]*loz+b[11])
        +std::pow(zm,8)*(std::pow(loz,3)/576.0+b[12]*loz*loz+b[13]*loz+b[14]);

  zbesselk33 *=-1.0;

  return zbesselk33;

}



double F0(double variable, double z, int n, void* data){

  double t = variable;
  int nu,nn;
  double fbb;
  double zz=z;
  nu=1;
  nn=n;
  if(nn == 0){
    fbb=std::exp(-zz*std::cosh(t))*std::cosh(nu*t);}
  else if(nn == 1){
    fbb=std::exp(-zz*std::cosh(t))*t*std::sinh(nu*t);}
  else if(nn == 2){
    fbb=std::exp(-zz*std::cosh(t))*t*t*std::cosh(nu*t);}
  else if(nn == 3){
    fbb=std::exp(-zz*std::cosh(t))*t*t*t*std::sinh(nu*t);}

  return fbb;

}