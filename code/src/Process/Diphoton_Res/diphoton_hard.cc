// Calculates the process dependent parts for the hard function calculations for the diphoton process
// from H2f.cpp Copyright 2010 Leandro <leandro@ubuntu>

#include "diphoton_hard.h"

#include "polylogs.h"


void sigmaijdiphotoncalc(diphoton_input* diph_in, PSpoint& PS, std::vector<std::vector<double> >& sigmaij, double alphas) {
// Born partonic cross-sections
// Normalisation is such that the COMPLETE dsigma/dOmega double differential cross-sections are returned
// There is also the option to include the gg channel Box contribution

    int Nf = diph_in->res_1.Nf;
    double gevpb = diph_in->res_1.gevm2topb;
    double alphaem = diph_in->res_1.alpha_QED;
    double Qu = 2./3.;
    double Qd = -1./3.;
    double pi = k_constants::pi;

    double s = 2.*PS.ss(0,1);
    double t = -2.*PS.ss(0,2);
    double u = -2.*PS.ss(1,2);
    double costheta = 1. + 2.*t/s;

    double BOX = 0.;
    double sigma0 = 0.;
//Box switch
    int boxflag = diph_in->boxflag;
    double flux = 1./(2.*s);
    if (boxflag == 0 || boxflag == 1){
      double KinFac = 8*(t/u + u/t); // |M|^2 summed over spins stripped of charges and coupling constants
      double spincolavg = 1/12.;
      sigma0 = gevpb*spincolavg*flux*(alphaem*alphaem)/2.*KinFac; // This includes the 1/(32pi^2) PS factor for d/domega, with 1/(16pi^2) absorbed into the definition of alphaem
      // std::cout << "sigma0 = " << sigma0 << std::endl;
      // std::cout << "gevpb = " << gevpb << std::endl;
      // std::cout << "spincolavg = " << spincolavg << std::endl;
      // std::cout << "flux = " << flux << std::endl;
      // std::cout << "alphaem^2 = " << alphaem*alphaem << std::endl;
      // std::cout << "KinFac = " << KinFac << std::endl;
    }
    if (boxflag == 1 || boxflag == 2){
      double spincolavgGG = 1/32.;
      BOX=gevpb*spincolavgGG*flux*(alphas*alphas*16.*pi*pi)*(alphaem*alphaem)/2.*ggBoxdiphotoncalc(costheta, diph_in->res_1.verbosity); // This has an additional 16pi^2 to compensate for the introduction of alpha_s
    }
    if (boxflag<0 || boxflag>2){
      std::cout << "boxflag must be 0,1,2" << std::endl;
      exit(EXIT_FAILURE);
    }

// Up quarks
    sigmaij[Nf+1][Nf-1] = sigma0*(std::pow(Qu,4));
    if (Nf>=4) sigmaij[Nf+4][Nf-4] = sigmaij[Nf+1][Nf-1];
    sigmaij[Nf-1][Nf+1] = sigmaij[Nf+1][Nf-1];
    if (Nf>=4) sigmaij[Nf-4][Nf+4] = sigmaij[Nf+1][Nf-1];
// Down quarks
    sigmaij[Nf+2][Nf-2] = sigma0*std::pow(Qd,4);
    if (Nf>=3) sigmaij[Nf+3][Nf-3] = sigmaij[Nf+2][Nf-2];
    if (Nf>=5) sigmaij[Nf+5][Nf-5] = sigmaij[Nf+2][Nf-2];
    sigmaij[Nf-2][Nf+2] = sigmaij[Nf+2][Nf-2];
    if (Nf>=3) sigmaij[Nf-3][Nf+3] = sigmaij[Nf+2][Nf-2];
    if (Nf>=5) sigmaij[Nf-5][Nf+5] = sigmaij[Nf+2][Nf-2];
// Gluons
    sigmaij[Nf][Nf] = BOX;

}


double ggBoxdiphotoncalc(double costheta, int verbosity) {
// Defined in terms of costheta, the cosine of the azimuthal angle of photon 1 in the partonic CM frame,
// which is just a convenient way of dimensionlessly parametrizing the Mandelstams invariants, since t/s = -1 + costheta.
    double s = 0.0, t = 0.0, u = 0.0;
    s = 1.;
    t = -1./2.*(1.+costheta);
    u = -1./2.*(1.-costheta);

    double sumQq2 = 11./9.;
    double pi = k_constants::pi;

    double sumQqto4 = sumQq2*sumQq2;
    double logmsot, logmsou, logtou, ratiostou, ratiosuot, ratiotuos;
    logmsot = std::log(-s/t);
    logmsou = std::log(-s/u);
    logtou = std::log(t/u);
    ratiostou = (s*s+t*t)/(u*u);
    ratiosuot = (s*s+u*u)/(t*t);
    ratiotuos = (t*t+u*u)/(s*s);
    double factor = 0.0, Box = 0.0;
    factor = 1./(4.*pi*pi)*(1.0/8.0*((ratiostou*logmsot*logmsot+2.0*(s-t)/u*logmsot)*
    (ratiostou*logmsot*logmsot+2.0*(s-t)/u*logmsot) + (ratiosuot*logmsou*logmsou+2.0*(s-u)/t*logmsou)*
    (ratiosuot*logmsou*logmsou+2.0*(s-u)/t*logmsou) + (ratiotuos*(logtou*logtou+pi*pi)+2.0*(t-u)/s*logtou)*
    (ratiotuos*(logtou*logtou+pi*pi)+2.0*(t-u)/s*logtou)) +
    1.0/2.0*(ratiostou*logmsot*logmsot+2.0*(s-t)/u*logmsot + ratiosuot*logmsou*logmsou+2.0*(s-u)/t*logmsou + ratiotuos*(logtou*logtou+pi*pi)+2.0*(t-u)/s*logtou)
    + pi*pi/2*((ratiostou*logmsot+(s-t)/u)*(ratiostou*logmsot+(s-t)/u) + (ratiosuot*logmsou+(s-u)/t)*(ratiosuot*logmsou+(s-u)/t))+4.0);

    Box = sumQqto4*factor;

    return Box;
}

double H1qdiphoton_DY(double costheta, int verbosity) {
// H1q^gammagamma in the scheme H1q^DY=0, as a function of yy == c_theta_CM,
// alpha_s normalization is alpha_s/Pi. H1qGaGa_DY in old cdoe.
// Done separately w. respect to H1diphotocalc below as that's in a different scheme

    double H1qdiphoton_DY = 0.0;
    H1qdiphoton_DY = (2.0/3.0 + (1.0/3.0)*(1.0/(1.0+costheta*costheta))*
      (std::log(1.0/2.0*(1.0+costheta))*((5.0-costheta)*(1.0-costheta)) +
      std::log(1.0/2.0*(1.0-costheta))*((5.0+costheta)*(1.0+costheta)) +
      (std::log(1.0/2.0*(1+costheta))*std::log(1.0/2.0*(1.0+costheta)))*
      (5.0+costheta*(2.0+costheta)) + (std::log(1.0/2.0*(1-costheta))*
      std::log(1.0/2.0*(1-costheta)))*(5.0+costheta*-(2.0-costheta))));
    if (verbosity >= 13) {
      std::cout << "H1qdiphoton_DY(costheta) = "<< H1qdiphoton_DY << std::endl;
    }

    return H1qdiphoton_DY;
}

double H1diphotoncalc(PSpoint* PS, double Cf, int verbosity) {
// H1q^gammagamma in "hard" scheme, as a function of ve=-u/s.
// alpha_s normalization is alpha_S/2Pi (!!). See 1311.1654. H1Gian in old code.
    double H1 = 0.;

    double s = 2.*PS->ss(0,1);
    double t = -2.*PS->ss(0,2);
    double u = -2.*PS->ss(1,2);
    double pi = k_constants::pi;
    double ve = 0.0;
    ve = -u/s;
    if (verbosity >= 13) {
      std::cout << "s = " << s << std::endl;
      std::cout << "t = " << t << std::endl;
      std::cout << "u = " << u << std::endl;
      std::cout << "In H1diphotoncalc, ve = " << ve << std::endl;
    }

    H1 = Cf*((pi*pi-7.0)+1.0/((1-ve)*(1-ve)+ve*ve)*(((1.0-ve)*(1.0-ve)+1.0)*
      std::log(1.0-ve)*std::log(1.0-ve)+ve*(ve+2.0)*std::log(1.0-ve) +
      (ve*ve+1.0)*std::log(ve)*std::log(ve) + (1.0-ve)*(3.0-ve)*std::log(ve)));
    if (verbosity >= 13) {
      std::cout << "H1 = " << H1 << " for ve = " << ve << std::endl;
    }
    return H1;
}

double H2qdiphotoncalc_DY (PSpoint* PS, double costheta, int quark,
                           double Cf, int Nf, int Nc, int verbosity) {
// As in H2qGaGa_DY of old code. Hard scheme -> DY scheme conversion for H2q coefficient,
// for retro-compatibility reasons I keep for now alpha_s/Pi normalization,
// take care with the different (alpha_s/2Pi) normalization of Giancarlo's
// (H1diphotoncalc) and Leandro's coefficients (H1qdiphoton_DY) etc.

    double ve = 0.0, H2qDY_hard = 0.0, C1qqDY = 0.0, H2qdiphoton_DY = 0.0;
    ve = 1.0/2.0*(1-costheta);
    double pi = k_constants::pi;
    double Z3 = k_constants::zeta3;

    C1qqDY = Cf/4.0*(pi*pi-8.0);
    H2qDY_hard = -2561.0/144.0 + 127.0*Nf/72.0 + 3.0*pi*pi/2.0
      - 19.0*Nf*pi*pi/81.0 + 49.0*pi*pi*pi*pi/324.0 + 58.0*Z3/9.0 + 8.0*Nf*Z3/27.0;
    if (verbosity >= 13) {
      std::cout << "costheta = " << costheta << std::endl;
      std::cout << "H2qDY_hard = " << H2qDY_hard << std::endl;
      std::cout << "C1qqDY = " << C1qqDY << std::endl;
      std::cout << "H1diphotoncalc(PS) = " << H1diphotoncalc(PS, Cf, verbosity) << std::endl;
      std::cout << "H2stqqdiphotoncalc(PS,quark) = " << H2stqqdiphotoncalc(PS, quark, Cf, Nf, Nc, verbosity) << std::endl;
    }
    H2qdiphoton_DY = H2stqqdiphotoncalc(PS, quark, Cf, Nf, Nc, verbosity)/4.0 - H2qDY_hard
      - 2*C1qqDY*(H1diphotoncalc(PS, Cf, verbosity)/2.0-2.0*C1qqDY);
    if (verbosity >= 13) {
	std::cout << "H2qdiphoton_DY = " << H2qdiphoton_DY << std::endl;
    }
    return H2qdiphoton_DY;

}


double H2stqqdiphotoncalc(PSpoint* PS, int quark, double Cf, int Nf, int Nc, int verbosity) {
    double muR2 = 0.0;
    double H2 = 0.0;
    double QQ[6], Qsqd[6];
    double s = 2.*PS->ss(0,1);
    double t = -2.*PS->ss(0,2);
    double u = -2.*PS->ss(1,2);
    double ve = 0.0, ALO = 0.0;
    double Qu = 2./3.;
    double Qd = -1./3.;
    double Z2 = k_constants::zeta2;
    double Z3 = k_constants::zeta3;
    double Z4 = k_constants::zeta4;

    ALO = 24.0*(u/t+t/u);
    ve = -u/s;
    muR2 = s; //Set as will always be the case in our calculations
    if (verbosity >= 13) {
      std::cout << "s = " << s << std::endl;
      std::cout << "t = " << t << std::endl;
      std::cout << "u = " << u << std::endl;
      std::cout << "ve = " << ve << std::endl;
      std::cout << "quark = " << quark << std::endl;
      std::cout << "muR2 = " << muR2 << std::endl;
      std::cout << "Qsqd = " << Qsqd << std::endl;
    }

    for (int i =0; i<5; i++) {
      QQ[i] = 0.0;
      Qsqd[i] = 0.0;
    }
    QQ[0] = 1.0;
    QQ[1] = Qu;
    QQ[2] = Qd;
    QQ[3] = Qd;
    QQ[4] = Qu;
    QQ[5] = Qd;
    Qsqd[0] = QQ[0]*QQ[0];
    Qsqd[1] = QQ[1]*QQ[1];
    Qsqd[2] = QQ[2]*QQ[2];
    Qsqd[3] = QQ[3]*QQ[3];
    Qsqd[4] = QQ[4]*QQ[4];
    Qsqd[5] = QQ[5]*QQ[5];

    if (verbosity >= 13) {
      std::cout << "quark = " << quark << " Qsqd[quark] = " << Qsqd[quark] << std::endl;
      std::cout << "Finite0x2 = " << Finite0x2(t,u,s,muR2,Qsqd[quark], Cf,Nf,Nc,verbosity) << std::endl;
      std::cout << "Finite1x1 = " << Finite1x1(t,u,s,muR2, Cf,Nc,verbosity) << std::endl;
      std::cout << "H1diphotoncalc(PS) = " << H1diphotoncalc(PS, Cf, verbosity) << std::endl;
    }

//No log(muR2/s) pieces as s = muR2 for us here.
    H2 = 1.0/ALO*(Finite0x2(t,u,s,muR2,Qsqd[quark], Cf,Nf,Nc,verbosity)
      + Finite1x1(t,u,s,muR2, Cf,Nc,verbosity)) + 6.0*Z2*Cf*H1diphotoncalc(PS, Cf,verbosity)
      +Cf*Cf*(-45.0*Z4+3.0/4.0*(1.0-8.0*Z2+16.0*Z3)*std::log(muR2/s))
      +Cf*Nc/648.0*(4856.0+21258.0*Z2-3366.0*Z3-8505.0*Z4) -Cf*Nf/162.0*(164.0+873*Z2-153.0*Z3);


    if (verbosity >= 13) {
	std::cout << "H2stqqdiphotoncalc = " << H2 << std::endl;
    }

    return H2;
}

double Finite1x1(double t, double u, double s,double muR2,
                 double Cf, int Nc, int verbosity) {
    double finite1x1 = 0.0;
    finite1x1 = Nc*Cf*Cf*(G1smitad(u,t,s,muR2, verbosity)+G1smitad(t,u,s,muR2, verbosity));
    if (verbosity >= 13) {
      std::cout << "finite1x1 = " << finite1x1 << std::endl;
    }
    return finite1x1;
}


double Finite0x2(double t, double u, double s, double muR2, double Qsqd,
                 double Cf, int Nf, int Nc, int verbosity) {
    double finite0x2 = 0.0;
    double Tr = 1./2.;
    double sumQq2 = 11./9.;

    if (verbosity >= 13) {
      std::cout << "t = " << t << " u = " << u << " s = " << s << " muR2 = " << muR2 << " Qsqd = " << Qsqd << std::endl;
    }
    finite0x2= 2.0*Nc*((sumQq2/Qsqd)*Tr*Cf*(Asmitad(u,t,s, verbosity)
      +Asmitad(t,u,s, verbosity))+Cf*Cf*(Bsmitad(u,t,s,muR2, verbosity)
      +Bsmitad(t,u,s,muR2, verbosity)) +Cf*Nc*(D2smitad(u,t,s,muR2, verbosity)
      +D2smitad(t,u,s,muR2, verbosity))+Nf*Cf*(E3smitad(u,t,s,muR2, verbosity)
      +E3smitad(t,u,s,muR2, verbosity)));
   if (verbosity >= 13) {
     std::cout << "finite0x2 = " << finite0x2 << std::endl;
     std::cout << 2.0*Nc*((sumQq2/Qsqd)*Tr*Cf) << " " << (Asmitad(u,t,s, verbosity)+Asmitad(t,u,s, verbosity))
       << " " << Cf*Cf*(Bsmitad(u,t,s,muR2, verbosity)+Bsmitad(t,u,s,muR2, verbosity))
       << " " << Cf*Nc*(D2smitad(u,t,s,muR2, verbosity)+D2smitad(t,u,s,muR2, verbosity)) << " "
       << Nf*Cf*(E3smitad(u,t,s,muR2, verbosity)+E3smitad(t,u,s,muR2, verbosity)) << std::endl;
     std::cout << "Tr = " << Tr << " sumQq2 = " << sumQq2 << " Qsqd = " << Qsqd << std::endl;
    }

    return finite0x2;
}

double Asmitad(double u, double t, double s, int verbosity) {
    double x = 0.0, y = 0.0, z = 0.0, X = 0.0, Y = 0.0, Asm = 0.0;
    x = -t/s;
    y = -u/s;
    z = -u/t;
    double Z3 = k_constants::zeta3;
    double pi = k_constants::pi;

    X = std::log(x);
    Y = std::log(y);
    Asm = (128.0*Li4(z)-128.0*Li4(x)+128.0*Li4(y)+(-64.0/3.0+128*Y)*Li3(x)
      +(64.0/3.0*X-64.0/3.0*pi*pi)*Li2(x)+16.0/3.0*X*X*X*X-64.0/3.0*X*X*X*Y
      +(-16.0+32.0/3.0*pi*pi+32.0*Y*Y)*X*X+(-64/3.0*pi*pi*Y+48.0+160.0/9.0*pi*pi)
      *X+64.0/3.0*Z3+224.0/45.0*pi*pi*pi*pi-128.0*Y*Z3)*(t/u)+(32.0/3.0*Li3(x)-32.0/3.0*Li3(y)
      +(-32.0/3.0*X-32.0/3.0*Y)*Li2(x)+(-32.0/3.0*Y*Y-80.0/9.0*pi*pi-64.0/3.0)*X+(64.0/3.0+32.0/3.0*pi*pi)*Y)
      *(t*t/s/s)+24.0*X*X*t*t/u/u+(416.0/3.0*Li3(x)+64.0*Li3(y)*X-416.0/3.0*Li2(x)*X+(8.0*Y*Y+16)*X*X
      +(-8.0/3.0*Y+80.0/3.0+112.0/9.0*pi*pi-64.0*Z3-64.0*Y*Y)*X
      -416.0/3.0*Z3-148.0/9.0*pi*pi+44.0/45.0*pi*pi*pi*pi);
    return Asm;
}

double Bsmitad(double u,double t,double s,double muR2, int verbosity) {
    double x = 0.0, y = 0.0, z = 0.0, X = 0.0, Y = 0.0, S = 0.0, Bsm = 0.0;
    x = -t/s;
    y = -u/s;
    z = -u/t;
    double Z3 = k_constants::zeta3;
    double pi = k_constants::pi;

    X = std::log(x);
    Y = std::log(y);
    S = std::log(s/muR2);
    Bsm = (-112.0*Li4(z)-88.0*Li4(y)+(-128.0*Y+48.0*X-64.0)*Li3(x)
      + (-16.0*Y-16.0*X+12.0)*Li3(y) + (12.0*Y-4.0*Y*Y + 8.0*X*X-8.0*pi*pi+64.0*X)*Li2(x)
      + 2.0/3.0*X*X*X*X+56.0/3.0*X*X*X*Y + (44.0*Y-4.0*pi*pi+2.0-32.0*Y*Y)*X*X
      + (-4.0*Y*Y*Y-8.0-32.0*Z3-80.0/3.0*pi*pi+6.0*Y*Y + 56.0/3.0*pi*pi*Y)*X + Y*Y*Y*Y
      + 6.0*Y*Y*Y + (-10.0/3.0*pi*pi-5.0)*Y*Y + (-39.0-18.0*pi*pi+144.0*Z3)*Y + 3*S
      + 187.0/4.0-4.0*pi*pi*S + 4.0/45.0*pi*pi*pi*pi-5.0*pi*pi-20.0*Z3+48.0*Z3*S)*t/u
      + (-12.0*X*X+(24.0*Y+24.0)*X-12.0*Y*Y-24.0*Y-12.0*pi*pi)*t*t/(s*s)+8.0*X*X*t*t/(u*u)
      + (-80.0*Li4(y)+32.0*X*Li3(x) + (-128.0*X-152.0)*Li3(y)+152.0*Li2(x)*X + 8.0*Y*Y*Li2(y)
      + (-16.0*Y*Y-24.0)*X*X + (60.0*Y*Y + (28.0+32.0/3.0*pi*pi)*Y-58.0)*X + 14.0/3.0*Y*Y*Y*Y
      + 44.0/3.0*Y*Y*Y + 8.0/3.0*Y*Y*pi*pi + (96.0*Z3-32.0/3.0*pi*pi)*Y
      +32.0/45.0*pi*pi*pi*pi+16.0*Z3-86.0/3.0*pi*pi-2.0);
    return Bsm;
}

double D2smitad(double u,double t,double s,double muR2, int verbosity) {
    double x = 0.0, y = 0.0, z = 0.0, X = 0.0, Y = 0.0, S = 0.0, D2sm =0.0;
    x = -t/s;
    y = -u/s;
    z = -u/t;
    double Z3 = k_constants::zeta3;
    double pi = k_constants::pi;

    X = std::log(x);
    Y = std::log(y);
    S = std::log(s/muR2);
    D2sm = (48.0*Li4(z)-16.0*Li4(x)+24.0*Li4(y)+(56.0*Y-8.0*X+20.0)*Li3(x)
      +(8.0*X-12.0+16.0*Y)*Li3(y)+(16.0/3.0*pi*pi-20.0*X-12.0*Y-8.0*X*X+4.0*Y*Y)*Li2(x)
      +1.0/3.0*X*X*X*X+(-8.0*Y-70.0/9.0)*X*X*X+(-4.0*pi*pi+286.0/9.0-16.0*Y+14.0*Y*Y-44.0/3.0*S)*X*X
      +(-22.0/9.0*pi*pi+4.0*Y*Y*Y-8.0*pi*pi*Y-6.0*Y*Y)*X-44.0/9.0*Y*Y*Y
      +(-4.0/3.0*pi*pi+35.0/9.0-22.0/3.0*S)*Y*Y +(57-26.0/9.0*pi*pi-72.0*Z3-22.0*S)*Y
      +479.0/9.0*Z3+19.0/60.0*pi*pi*pi*pi-52.0*Z3*S+1141.0/27.0*S -215.0/18.0*pi*pi-43417.0/324.0+23.0/6.0*pi*pi*S)*t/u
      +(6.0*X*X+(-12.0-12.0*Y)*X+6.0*Y*Y+12*Y+6.0*pi*pi)*t*t/(s*s)-6.0*X*X*t*t/(u*u) +(16.0*Li4(y)+48.0*Li3(x)*Y+64.0*Li3(y)
      -8.0*Y*Y*Li2(y)-64.0*Li2(x)*X -4.0/3.0*X*X*X*X+(-20.0/3.0*pi*pi+6.0*Y*Y)*X*X+(-24.0*Y*Y+(-16.0/3.0*pi*pi-14.0)*Y-148.0/9.0*pi*pi)*X
      -112.0/9.0*Y*Y*Y+(-44.0/3.0*S+298.0/9.0)*Y*Y+(538.0/9.0-48.0*Z3-44.0/3.0*S)*Y-8.0*Z3-1.0/3.0*pi*pi*pi*pi+61.0/9.0*pi*pi);
    return D2sm;
}

double E3smitad(double u,double t,double s,double muR2, int verbosity) {
    double X = 0.0, Y = 0.0, S = 0.0, E3sm = 0.0;
    double Z3 = k_constants::zeta3;
    double pi = k_constants::pi;

    X = std::log(-t/s);
    Y = std::log(-u/s);
    S = std::log(s/muR2);
    E3sm = (16.0/9.0*X*X*X+(-76.0/9.0+8.0/3.0*S)*X*X+16.0/9.0*pi*pi*X
      +8.0/9.0*Y*Y*Y+(4.0/3.0*S-2.0/9.0)*Y*Y +(8.0/9.0*pi*pi+4.0*S-10.0)*Y
      -1.0/3.0*pi*pi*S-202.0/27.0*S+19.0/9.0*pi*pi-2.0/9.0*Z3+3401.0/162.0)*t/u
      +(16.0/9.0*pi*pi*X+16.0/9.0*Y*Y*Y+(8.0/3.0*S-52.0/9.0)*Y*Y+(-76.0/9.0+8.0/3.0*S)*Y
      +8.0/9.0*pi*pi);
    return E3sm;
}

double G1smitad(double u,double t,double s,double muR2, int verbosity) {
    double X = 0.0, Y = 0.0, S = 0.0, G1sm = 0.0;
    double pi = k_constants::pi;

    X = std::log(-t/s);
    Y = std::log(-u/s);
    S = std::log(s/muR2);
    G1sm= (14.0*X*X*X*X+28.0*X*X*X+8.0*X*X*Y*Y+56.0*X*X*pi*pi
      -48.0*X*X+12.0*X*X*Y+32.0*X*Y*pi*pi+80.0*pi*pi*X+2.0*Y*Y*Y*Y+12.0*Y*Y*Y-10.0*Y*Y
      +8.0*Y*Y*pi*pi+26.0*pi*pi+24.0*pi*pi*Y-84.0*Y+102.0)*t/u+8.0*X*(X*X*X+X*X+4.0*pi*pi*X+2.0*pi*pi)*t*t/u/u
      +2.0*X*X*(X*X+4.0*pi*pi)*t*t*t/u/u/u
      +(32.0*X*X*X+8.0*X*X*Y*Y+80.0*pi*pi*X+32.0*X*Y*pi*pi+8.0*Y*Y*X+8.0*Y*Y*Y*Y+32.0*Y*Y*pi*pi
        -32.0*Y*Y-4.0-56.0*Y+24.0*pi*pi);
    return G1sm;
}
