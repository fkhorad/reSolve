// Calculates the process dependent parts for the hard function calculations,
// here only diphoton is included for now, much of code is repackaged and rewritten
// from H2f.cpp Copyright 2010 Leandro <leandro@ubuntu>

#include "diphoton_hard.h"


void sigmaijdiphotoncalc(diphoton_input* diph_in, PSpoint& PS, double jacob,
                         std::vector<std::vector<double> >& sigmaij, double alphas) {
// Calculates the Born cross section and the complete modified cross-section which appears in the "universal"
// resummation formalism, in the DY resummation scheme. It returns the Born cross-section along with the H1q and H2q coefficients.
// It includes also the gg Box contribution. N.B. Without the 1/s "flux factor".

//  Note in old code sigmaij was indexed -5 to 5 in each direction but here it's 0 to 10 in each for simplicity.
    sigmaij.resize(11);
    for (int a = 0; a<11; a++) {
        sigmaij[a].resize(11,0.);
    }
    if (diph_in->verbosity >= 12) {
      std::cout << "sigmaijdiphotoncalc: " << std::endl;
    }
    double gevpb = diph_in->gevm2topb;
    double alphaem = diph_in->alpha_QED;
    double Qu = 2./3.;
    double Qd = -1./3.;
    double pi = k_constants::pi;

    double s = 2.*PS.ss(0,1);
    double t = -2.*PS.ss(0,2);
    double u = -2.*PS.ss(1,2);
    double costheta = 1. + 2.*t/s;

    if (diph_in->verbosity >= 12) {
      std::cout << "procdephard_diphoton.cc" << std::endl;
      std::cout << "costheta = " << costheta << std::endl;
      std::cout << "jacob = " << jacob << std::endl;
    }

    double KinFac = jacob*(t/u + u/t);
    if (diph_in->verbosity >= 12) {
      std::cout << "KinFac = " << KinFac << std::endl;
      std::cout << "for kinfac: " << std::endl;
    }



    double BOX = 0.;
    double sigma0 = 0.;
//Box switch
    int boxflag = diph_in->boxflag;
    if (boxflag == 0 || boxflag == 1){
      sigma0 = gevpb*pi*(alphaem*alphaem)/3.0;
    }
    if (boxflag == 1 || boxflag == 2){
      BOX=jacob*ggBoxdiphotoncalc(costheta, diph_in->verbosity);
    }
    if (boxflag<0 || boxflag>2){
      std::cout << "boxflag must be 0,1,2" << std::endl;
      exit(EXIT_FAILURE);
    }


    sigmaij[6][4] = sigma0*(std::pow(Qu,4))*KinFac;
    sigmaij[9][1] = sigmaij[6][4];
    sigmaij[4][6] = sigmaij[6][4];
    sigmaij[1][9] = sigmaij[6][4];
    sigmaij[7][3] = sigma0*std::pow(Qd,4)*KinFac;
    sigmaij[8][2] = sigmaij[7][3];
    sigmaij[10][0] = sigmaij[7][3];
    sigmaij[3][7] = sigmaij[7][3];
    sigmaij[2][8] = sigmaij[7][3];
    sigmaij[0][10] = sigmaij[7][3];
    sigmaij[5][5] = gevpb*((alphas*alphas)*(alphaem*alphaem)/(64*pi))*BOX;


      if (diph_in->verbosity >= 12) {
        std::cout << "BOX = " << BOX << std::endl;
      }

      if (diph_in->verbosity >= 13) {
        std::cout << "sigma0 = " << sigma0 << std::endl;
        std::cout << "gevpb = " << gevpb << std::endl;
        std::cout << "alphaem = " << alphaem << std::endl;
        std::cout << "KinFac = " << KinFac << std::endl;
      }
      if (diph_in->verbosity >= 13) {
        std::cout << "sigmaij = "<< std::endl;
        for (int i = 0; i <11; i++) {
          for (int j = 0; j <11; j++) {
          std::cout << "sigmaij[" << i<< "][" << j << "] = " << sigmaij[i][j] << std:: endl;
          }
        }
      }

    if (diph_in->verbosity >= 12) {
      std::cout << "KinFac = " << KinFac << std::endl;
      std::cout << "BOX = " << BOX << std::endl;
    }

}


double ggBoxdiphotoncalc(double costheta, int verbosity) {
// Defined in terms of costheta, the cosine of the azimuthal angle of photon 1 in the partonic CM frame,
// which is just a convenient way of dimensionlessly parametrizing the Mandelstams invariants, since t/s = -1 + costheta.
    double s = 0.0, t = 0.0, u = 0.0;
    s = 1.0;
    t = -1.0/2.0*(1+costheta);
    u = -1.0/2.0*(1-costheta);
    if (verbosity >= 13) {
	std::cout << "ggboxdiphotoncalc: " << std::endl;
	std::cout << "s = " << s << std::endl;
	std::cout << "t = " << t << std::endl;
	std::cout << "u = " << u << std::endl;
    }

    double sumQq2 = 11./9.;
    double pi = k_constants::pi;

    double sumQqto4 = 0.0;
    sumQqto4 = sumQq2*sumQq2;
    double av = 0.0, logmsot = 0.0, logmsou = 0.0, logtou = 0.0,
    ratiostou = 0.0, ratiosuot = 0.0, ratiotuos = 0.0;
    av = 1.0/256.0;
    logmsot = std::log(-s/t);
    logmsou = std::log(-s/u);
    logtou = std::log(t/u);
    ratiostou = (s*s+t*t)/(u*u);
    ratiosuot = (s*s+u*u)/(t*t);
    ratiotuos = (t*t+u*u)/(s*s);
    double factor = 0.0, Box = 0.0;
    factor = av*16*2.0/(pi*pi)*(1.0/8.0*((ratiostou*logmsot*logmsot+2.0*(s-t)/u*logmsot)*
    (ratiostou*logmsot*logmsot+2.0*(s-t)/u*logmsot) + (ratiosuot*logmsou*logmsou+2.0*(s-u)/t*logmsou)*
    (ratiosuot*logmsou*logmsou+2.0*(s-u)/t*logmsou) + (ratiotuos*(logtou*logtou+pi*pi)+2.0*(t-u)/s*logtou)*
    (ratiotuos*(logtou*logtou+pi*pi)+2.0*(t-u)/s*logtou)) +
    1.0/2.0*(ratiostou*logmsot*logmsot+2.0*(s-t)/u*logmsot + ratiosuot*logmsou*logmsou+2.0*(s-u)/t*logmsou + ratiotuos*(logtou*logtou+pi*pi)+2.0*(t-u)/s*logtou)
    + pi*pi/2*((ratiostou*logmsot+(s-t)/u)*(ratiostou*logmsot+(s-t)/u) + (ratiosuot*logmsou+(s-u)/t)*(ratiosuot*logmsou+(s-u)/t))+4.0);

    Box = 1.0/2.0*4*sumQqto4*(2.0*pi)*(2.0*pi)*(2.0*pi)*(2.0*pi)*factor; //Factor 1/2 for identical particles

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
