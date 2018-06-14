#include "polylogs.h"

#include "constants.h"


double Li4(double x){
    double Li4 = 0.0;
    double Z3 = k_constants::zeta3;
    double Z4 = k_constants::zeta4;
    double pi = k_constants::pi;
    double PolyLog2m2 = k_constants::PolyLog2m2;
    double PolyLog3m2 = k_constants::PolyLog3m2;
    double PolyLog4m2 = k_constants::PolyLog4m2;

    if(x<=-2.0){
      Li4 = 1.0/360.0*(-60.0/x*(6.0+ 6.0/625.0/(x*x*x*x) + 3.0/128.0/(x*x*x) + 2.0/27.0/(x*x) + 3.0/8.0/x )-7.0*pi*pi*pi*pi
        - 30.0*pi*pi*std::log(-x)*std::log(-x) - 15.0*std::log(-x)*std::log(-x)*std::log(-x)*std::log(-x));
    }
    else if(x>-2.0 && x<-1.0){
      Li4 = (std::pow((2.0 + x),10.0)* (-25437292.0 + 71241525.0* std::log(3.0) + 62364492.0* PolyLog2m2 - 22044960.0 *PolyLog3m2))/225740390400.0
        + (std::pow((2.0 + x),9.0)* (-2404252.0 + 7176033.0* std::log(3.0) + 6657228.0* PolyLog2m2 - 2449440.0* PolyLog3m2))/11287019520.0
        + (std::pow((2.0 + x),8.0)* (-82123.0 + 265923.0* std::log(3.0) + 264627.0* PolyLog2m2 - 102060.0* PolyLog3m2))/209018880.0
        + (std::pow((2.0 + x),6.0)* (-1438.0 + 6075.0* std::log(3.0) + 7398.0* PolyLog2m2 - 3240.0* PolyLog3m2))/1244160.0
        + (std::pow((2.0 + x),5.0)* (-58.0 + 315.0* std::log(3.0) + 450.0* PolyLog2m2 - 216.0* PolyLog3m2))/34560.0
        + (std::pow((2.0 + x),4.0)* (-2.0 + 18.0* std::log(3) + 33.0* PolyLog2m2 - 18.0* PolyLog3m2))/1152.0
        + 1.0/48.0 *std::pow((2.0 + x),3.0)* (std::log(3.0) + 3.0* PolyLog2m2 - 2.0* PolyLog3m2) + 1.0/8.0* std::pow((2.0 + x),2.0)* (PolyLog2m2 - PolyLog3m2)
        - 1.0/2.0* (2.0 + x)* PolyLog3m2 + (std::pow((2.0 + x),7.0)* (10962.0* std::log(3.0) + 11907.0* PolyLog2m2 - 5.0* (607.0 + 972.0* PolyLog3m2)))/4354560.0 + PolyLog4m2;
    }
    else if(x<=0.7 && x>=-1.0){
      Li4 = x+x*x/16+std::pow(x,3.0)/81+std::pow(x,4.0)/256+std::pow(x,5.0)/625
      +std::pow(x,6.0)/1296+std::pow(x,7.0)/2401+std::pow(x,8.0)/4096+std::pow(x,9.0)/6561
      +std::pow(x,10.0)/10000+std::pow(x,11.0)/14641+std::pow(x,12.0)/20736+std::pow(x,13.0)/28561
      +std::pow(x,14.0)/38416+std::pow(x,15.0)/50625+std::pow(x,16.0)/65536+std::pow(x,17.0)/83521
      +std::pow(x,18.0)/104976+std::pow(x,19.0)/130321+std::pow(x,20.0)/160000+std::pow(x,21.0)/194481
      +std::pow(x,22.0)/234256+std::pow(x,23.0)/279841+std::pow(x,24.0)/331776+std::pow(x,25.0)/390625
      +std::pow(x,26.0)/456976+std::pow(x,27.0)/531441+std::pow(x,28.0)/614656+std::pow(x,29.0)/707281+std::pow(x,30.0)/810000
      +std::pow(x,31.0)/923521+std::pow(x,32.0)/1048576+std::pow(x,33.0)/1185921+std::pow(x,34.0)/1336336
      +std::pow(x,35.0)/1500625+std::pow(x,36.)/1679616+std::pow(x,37.0)/1874161+std::pow(x,38.0)/2085136
      +std::pow(x,39.0)/2313441+std::pow(x,40.0)/2560000+std::pow(x,41.0)/2825761+std::pow(x,42.0)/3111696
      +std::pow(x,43.0)/3418801+std::pow(x,44.0)/3748096+std::pow(x,45.0)/4100625+std::pow(x,46.0)/4477456
      +std::pow(x,47.0)/4879681+std::pow(x,48.0)/5308416+std::pow(x,49.0)/5764801+std::pow(x,50.0)/6250000;
    }
    else if(x<1.0 && x>0.7){
      Li4 = 1.0/457228800.0*(5080320.0*pi*pi*pi*pi-std::pow((-1.0+x),3.0)
      *(-1381393255.0+x* (4840853127.0+x* (-9435621561.0+x* (11568105449.0
        +x* (-9128211801.0+x* (4518682089.0+x* (-1281356743.0+159233895* x)))))))
        +1512.0*pi*pi* std::pow((-1.0+x),2.0)* (177133.0+x* (-617934.0
          +x* (1341449.0+x* (-1931968.0+x* (1883165.0+x* (-1233718.0+x* (522099.0+2.0* x* (-64642.0+7129.0* x))))))))
          +2520.0* std::pow((-1.0+x),3.0)* (-420475.0+x* (1615443.0+x* (-3282009.0+x* (4114961.0+x *(-3292089.0
            +x* (1644801.0+x *(-469507.0+58635.0* x)))))))* std::log(1.0-x)-181440.0* (-1.0+x)* (-7381.0
              +x* (17819.0+x* (-38881.0+x* (61919.0+x* (-70381.0+x* (56627.0
                +x* (-31573.0+7.0* x* (1661.0+4.0* x* (-91.0+9.0* x)))))))))* Z3);
    }
    else if(x==1){
      Li4 = Z4;
    }
    else {
      std::cout<< "Bad value of x in Li4 - " << x << std::endl;
      Li4 = 0.0;
    }

// Only Li4 receives values of Z, Li2 and Li3 are only evaluated in x and y,
// their arguments are bounded between 0 and 1 so they do not throw errors

    return Li4;
}


double Li3(double x) {
    double Li3 = 0.0;
    double xx = 0.0;
    double Z2 = k_constants::zeta2;

    if (x>1.0 && x < 1.0+1e-8)	{
      xx = 1.0;
    }
    else {
      xx = x;
    }

    if (xx > 0 && xx <1.0) {
      Li3 = Li3fn(xx);
    }
    else if (xx>-1.0 && xx<0.0) {
      Li3 = -Li3fn(-xx) + Li3fn(xx*xx)/4.0;
    }
    else if (xx < -1.0) {
      Li3 = -Li3fn(1.0/xx) + Li3fn(1.0/(xx*xx))/4.0 + Z2*std::log(-1.0/xx)
        + std::log(-1.0/xx)*std::log(-1.0/xx)*std::log(-1.0/xx)/6.0;
    }
    else {
      std::cout << "Wrong argument in Li3" << std::endl;
      exit(EXIT_FAILURE);
    }
    return Li3;
}

double Li3fn(double xx) {
    double Li3_0 = 0.0, yy = 0.0;
    double Z2 = k_constants::zeta2;
    double Z3 = k_constants::zeta3;

    if (xx <0.0 || xx>1.0) {
      std::cout << "wrong argument in Li3fn" << std::endl;
      exit(EXIT_FAILURE);
    }
    else if (xx < 0.35) {
      yy = std::log(1.0-xx);
      Li3_0 = -yy - (3*yy*yy)/8.0 - (17*yy*yy*yy)/216.0 - (5*yy*yy*yy*yy)/576.0
        - (7*yy*yy*yy*yy*yy)/54000.0 + (7*yy*yy*yy*yy*yy*yy)/86400.0 + 19*yy*yy*yy*yy*yy*yy*yy/5556600.0
        - yy*yy*yy*yy*yy*yy*yy*yy/752640.0 - 11*yy*yy*yy*yy*yy*yy*yy*yy*yy/127008000.0
        + 11*yy*yy*yy*yy*yy*yy*yy*yy*yy*yy/435456000.0;
    }
    else if (xx>0.35 && xx<1.0) {
      yy = std::log(xx);
      Li3_0 = Z3 + Z2*yy - (yy*yy*yy)/12.0 - (yy*yy*yy*yy)/288.0 + (yy*yy*yy*yy*yy*yy)/86400.0
        - (yy*yy*yy*yy*yy*yy*yy*yy)/10160640.0 + (yy*yy*yy*yy*yy*yy*yy*yy*yy*yy)/870912000.0
        + (yy*yy)/4.0*(3.0-2.0*std::log(-yy));
    }
    else if (xx == 1.0) {
      Li3_0 = Z3;
    }
    return Li3_0;
}


double Li2 (double x) {
    double Li2 = 0.0, xx =0.0;
    double Z2 = k_constants::zeta2;

    if (x > 1.0 && x < 1.0+1e-8) {
      xx = 1.0;
    }
    else {
      xx = x;
    }

    if (xx>0 && x<1.0) {
      Li2 = Li2fn(xx);
    }
    else if (xx>-1.0 && xx<0.0) {
      Li2 = -Li2fn(-xx) + Li2fn(xx*xx)/2.0;
    }
    else if (xx < -1.0) {
      Li2 = Li2fn(-1.0/xx) - Li2fn(1.0/(xx*xx))/2.0 - Z2 - std::log(-1.0/xx)*std::log(-1.0/xx)/2.0;
    }
    else {
      std::cout << "Wrong argument in Li2" << std::endl;
      exit(EXIT_FAILURE);
    }
    return Li2;
}

double Li2fn (double xx) {
    double Li2_0 = 0.0, yy = 0.0;
    double Z2 = k_constants::zeta2;

    if (xx < 0.0 || xx > 1.0) {
      std::cout << "Wrong argument in Z2 in Li2fn" << std::endl;
      exit(EXIT_FAILURE);
    }
    else if (xx < 0.35) {
      yy = std::log(1.0-xx);
      Li2_0 = -yy - (yy*yy)/4.0 - (yy*yy*yy)/36.0 + (yy*yy*yy*yy*yy)/3600.0
        - (yy*yy*yy*yy*yy*yy*yy)/211680.0 + (yy*yy*yy*yy*yy*yy*yy*yy*yy)/10886400.0;
    }
    else if (xx>0.35 && xx<1.0) {
      yy = std::log(xx);
      Li2_0 = Z2 - (yy*yy)/4.0 - (yy*yy*yy)/72.0 + (yy*yy*yy*yy*yy)/14400.0 - (yy*yy*yy*yy*yy*yy*yy)/1270080.0
        + (yy*yy*yy*yy*yy*yy*yy*yy*yy)/87091200.0 + yy*(1.0-std::log(-yy));
    }
    else if (xx == 1.0) {
      Li2_0 = Z2;
    }
    return Li2_0;
}

std::complex<double> Li2(double x, double y){
    double r = x/y;
    double DEB_PI6 = k_constants::pi2_6;
    if (r >= 0.0) return Li2(1.0-r);
    else return DEB_PI6 - Li2(r) - std::log(1.0-r)*(Log(x) - Log(y));
}


std::complex<double> Log(double x){
    double pi = k_constants::pi;
    return std::complex<double> ( std::log(std::abs(x)), (x < 0. ? pi : 0.) );
}
std::complex<double> L0(double x, double y){
    return (Log(x) - Log(y))/(1.0 - x/y);
}
std::complex<double> L1(double x, double y){
    return (L0(x, y) + 1.0)/(1.0 - x/y);
}
std::complex<double> L2(double x, double y){
    double r = x/y;
    double omr = 1-r;
    double omr3 = omr*omr*omr;

    return (Log(x)-Log(y)-(r-1.0/r)*0.5)/omr3;
}
std::complex<double> Ls1(double x, double y, double z, double w) {
    return (Ls0(x,y,z,w) + L0(x,y) +  L0(z,w))/(1.0-x/y-z/w);
}
std::complex<double> Ls0(double x, double y, double z, double w) {
    return Ls_1(x,y,z,w)/(1.0-x/y-z/w);
}
std::complex<double> Ls_1(double x, double y, double z, double w) {

    double DEB_PI6 = k_constants::pi2_6;

    return Li2(x, y) + Li2(z, w) + (Log(x)-Log(y))*(Log(z)-Log(w)) - DEB_PI6;
}
