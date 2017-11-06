
#include "evolution.h"




// Much of the stuff that's inside here is redundant and sometimes slightly inconsistent


//
// PROVIDES MOMENTS OF DENSITIES AT A GIVEN SCALE (INCLUDES EVOLUTION)
// AND ANOMALOUS DIMENSIONS
// EVERYTHING IN MELLIN SPACE
//

void RENO2 (std::complex<double> XN, int I, int ISIGN, int IBEAM, int flag1,
             int nf,
             std::complex<double> QQI, std::complex<double> QGF, std::complex<double> GQI,
             std::complex<double> GGI, std::complex<double> GGF, std::complex<double> NS1MI,
             std::complex<double> NS1PI, std::complex<double> NS1F, std::complex<double> QQ1F,
             std::complex<double> QG1F, std::complex<double> GQ1I, std::complex<double> GQ1F,
             std::complex<double> GG1I,std::complex<double>GG1F, std::complex<double> UVI,
             std::complex<double> DVI, std::complex<double> USI, std::complex<double> DSI,
             std::complex<double> SSI, std::complex<double>GLI, std::complex<double>CHI,
             std::complex<double> BOI,
	    double ALPS, double q2, double b0p, std::complex<double> ALPQ,
	    double a_param, double as, double mur2, int verbosity,
            std::complex<double>* FN){
//
//
      double XL, XL1, S, ALP;
      int F, iord, naord;
      std::complex<double> UVN,DVN,NS3N,NS8N,GLN,SIN,NS15N,NS24N,NS35N,
       ALPr,BON,CHN;
      std::complex<double> ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,
                 RGQ, RGG, C2Q, C2G, CDYQ, CDYG,
                 C2QI, C2GF, CDYQI, CDYGI;
      std::complex<double> ENS,EM,EP,AC,EMP,EPM,NMP,NPM,DMQQ,DMQG,DMGQ,DMGG,
       DPQQ,DPQG,DPGQ,DPGG,RMMQQ,RMMQG,RMMGQ,RMMGG,RMPQQ,RMPQG,RMPGQ,
       RMPGG,RPMQQ,RPMQG,RPMGQ,RPMGG,RPPQQ,RPPQG,RPPGQ,RPPGG,GL,SG,
       SSN,USN,DSN;
//
//
      if (verbosity >= 50) {
	  std::cout << "In RENO2: " << std::endl;
	  std::cout << "qqi = " << QQI << std::endl;
	  std::cout << "qgf = " << QGF << std::endl;
	  std::cout << "gqi = " << GQI << std::endl;
	  std::cout << "ggi = " << GGI << std::endl;
	  std::cout << "ggf = " << GGF << std::endl;
	  std::cout << "ns1mi = " << NS1MI << std::endl;
	  std::cout << "ns1pi = " << NS1PI << std::endl;
	  std::cout << "ns1f = " << NS1F << std::endl;
	  std::cout << "qq1f = " << QQ1F << std::endl;
	  std::cout << "qg1f = " << QG1F << std::endl;
	  std::cout << "gq1i = " << GQ1I << std::endl;
	  std::cout << "gq1f = " << GQ1F << std::endl;
	  std::cout << "gg1i = " << GG1I << std::endl;
	  std::cout << "gg1f = " << GG1F << std::endl;
      }

      naord = flag1;
      if (flag1 < 1) naord=1;

      UVN = UVI;
      DVN = DVI;
      NS3N = UVI + 2.*USI - DVI - 2.*DSI;
      NS8N = UVI + 2.*USI + DVI + 2.*DSI - 4.*SSI;
      GLN = GLI;
      SIN = UVI + DVI + 2.*USI + 2.*DSI
          + 2.*SSI + 2.*CHI + 2.*BOI;
      NS15N = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI - 6.*CHI;
      NS24N = UVI + DVI + 2.*USI + 2.*DSI +
          2.*SSI + 2.*CHI - 8.*BOI;
      NS35N = SIN;

      if (nf==3){
        SIN = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI;
        NS15N = SIN;
        NS24N = SIN;
        NS35N = SIN;
      }
      F = 5;

      XL = ALPS / ALPQ.real(); //CHECK IMPLICIT TAKING OF REAL PART IN OLD CODE 
      ALP=ALPQ.real(); //CHECK IMPLICIT TAKING OF REAL PART IN OLD CODE

      if (verbosity >= 50) {
	  std::cout << "I = " << I << std::endl;
	  std::cout << "ISIGN = " << ISIGN << std::endl;
	  std::cout << "IBEAM = " << IBEAM << std::endl;
	  std::cout << "XL = " << XL << std::endl;
	  std::cout << "UVN = " << UVN << std::endl;
	  std::cout << "DVN = " << DVN << std::endl;
	  std::cout << "NS3N = " << NS3N << std::endl;
	  std::cout << "NS8N = " << NS8N << std::endl;  
	  std::cout << "GLN = " << GLN << std::endl; 
	  std::cout << "SIN = " << SIN << std::endl;  
	  std::cout << "NS15N = " << NS15N << std::endl;
	  std::cout << "NS24N = " << NS24N << std::endl;
	  std::cout << "NS35N = " << NS35N << std::endl;
	  std::cout << "ALPS = " << ALPS << std::endl;
	  std::cout << "ALPQ = " << ALPQ << std::endl;
	  std::cout << "XL = " << XL << std::endl;
      }


      S   = std::log(XL);
      XL1 = 1.- XL;
// SELECT ORDER FOR EVOLUTION LO/NLO
      double naordm1 = naord-1.;
      ALPr= ALP * std::complex<double>(1.,0.)*naordm1;

      ANOM (ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,
               RGQ, RGG, C2Q, C2G, CDYQ, CDYG, XN, F,
               QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
               QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF,
	    CDYQI, CDYGI, verbosity, nf);


       ENS = std::exp (-ANS*S);

       if (verbosity >= 50) {
	   if (ENS.real() > 1e50) {
	       std::cout << "ANS = " << ANS << " S = " << S << std::endl;
	   }
       }
       
//
       EM  = std::exp (-AM*S);
       EP  = std::exp (-AP*S);

       if (verbosity >= 50) {
	   std::cout << "AP = " << AP << " S = " << S << std::endl;
       }

       AC  = 1.- AL;
       EMP = std::exp (S * (-AM+AP));
       EPM = std::exp (S * (+AM-AP));
       NMP = 1.- AM + AP;
       NPM = 1.- AP + AM;
       DMQQ =  AL * RQQ + BE * RGQ;
       DMQG =  AL * RQG + BE * RGG;
       DMGQ =  AB * RQQ + AC * RGQ;
       DMGG =  AB * RQG + AC * RGG;
       DPQQ =  AC * RQQ - BE * RGQ;
       DPQG =  AC * RQG - BE * RGG;
       DPGQ = -AB * RQQ + AL * RGQ;
       DPGG = -AB * RQG + AL * RGG;
       RMMQQ =   AL * DMQQ + AB * DMQG;
       RMMQG =   BE * DMQQ + AC * DMQG;
       RMMGQ =   AL * DMGQ + AB * DMGG;
       RMMGG =   BE * DMGQ + AC * DMGG;
       RMPQQ =  (AC * DMQQ - AB * DMQG) / NMP;
       RMPQG = (-BE * DMQQ + AL * DMQG) / NMP;
       RMPGQ =  (AC * DMGQ - AB * DMGG) / NMP;
       RMPGG = (-BE * DMGQ + AL * DMGG) / NMP;
       RPMQQ =  (AL * DPQQ + AB * DPQG) / NPM;
       RPMQG =  (BE * DPQQ + AC * DPQG) / NPM;
       RPMGQ =  (AL * DPGQ + AB * DPGG) / NPM;
       RPMGG =  (BE * DPGQ + AC * DPGG) / NPM;
       RPPQQ =   AC * DPQQ - AB * DPQG;
       RPPQG =  -BE * DPQQ + AL * DPQG;
       RPPGQ =   AC * DPGQ - AB * DPGG;
       RPPGG =  -BE * DPGQ + AL * DPGG;

//...EVOLUTION OF LIGHT PARTON DENSITIES
         UVN  = UVN  * ENS * (1.+  ALPr * XL1 * RMIN);
         DVN  = DVN  * ENS * (1.+  ALPr * XL1 * RMIN);
         NS3N = NS3N * ENS * (1.+  ALPr * XL1 * RPLUS);
         NS8N = NS8N * ENS * (1.+  ALPr * XL1 * RPLUS);

	 if (verbosity >= 50) {
	     if (UVN.real() > 1e100) {
		 std::cout << "UVN = " << UVN << std::endl;
		 std::cout << "ENS = " << ENS << std::endl;
		 std::cout << "ALPr = " << ALPr << " XL1 = " << XL1 << " RMIN = " << RMIN << std::endl;
	     }
	 }

//
       SG = SIN;
       GL = GLN;
       SIN = EM * ((AL + ALPr * (RMMQQ * XL1 + RMPQQ * (EPM-XL)))* SG
               + (BE + ALPr * (RMMQG * XL1 + RMPQG * (EPM-XL))) * GL)
         + EP * ((AC + ALPr * (RPPQQ * XL1 + RPMQQ * (EMP-XL))) * SG
               +(-BE + ALPr * (RPPQG * XL1 + RPMQG * (EMP-XL))) * GL);
       GLN = EM * ((AB + ALPr * (RMMGQ * XL1 + RMPGQ * (EPM-XL))) * SG
               + (AC + ALPr * (RMMGG * XL1 + RMPGG * (EPM-XL))) * GL)
         + EP *((-AB + ALPr * (RPPGQ * XL1 + RPMGQ * (EMP-XL))) * SG
               + (AL + ALPr * (RPPGG * XL1 + RPMGG * (EMP-XL))) * GL);
//
      NS15N = NS15N * ENS * (1.+  ALPr * XL1 * RPLUS);
      NS24N = NS24N * ENS * (1.+  ALPr * XL1 * RPLUS);
      NS35N = SIN;

      if (verbosity >= 50) {
	  std::cout << "SIN testing: "<< std::endl;;
	  std::cout << "EM = " << EM << std::endl;
	  std::cout << "AL = " << AL << std::endl;
	  std::cout << "ALPr = " << ALPr << std::endl;
	  std::cout << "RMMQQ = " << RMMQQ << std::endl;
	  std::cout << "XL1 = " << XL1 << std::endl;
	  std::cout << "RMPQQ = " << RMPQQ << std::endl;
	  std::cout << "EPM = " << EPM << std::endl;
	  std::cout << "XL = " << XL << std::endl;
	  std::cout << "SG = " << SG << std::endl;
	  std::cout << "BE = " << BE << std::endl;
	  std::cout << "RMMQG = " << RMMQG << std::endl;
	  std::cout << "RMPQG = " << RMPQG << std::endl;
	  std::cout << "GL = " << GL << std::endl;
	  std::cout << "EP = " << EP << std::endl;
	  std::cout << "AC = " << AC << std::endl;
	  std::cout << "RPPQQ = " << RPPQQ << std::endl;
	  std::cout << "RPMQQ = " << RPMQQ << std::endl;
	  std::cout << "EMP = " << EMP << std::endl;
	  std::cout << "EMP = " << EMP << std::endl;
	  std::cout << "RPPQG = " << RPPQG << std::endl;
	  std::cout << "RPMQG = " << RPMQG << std::endl;
      }


//...  FLAVOUR DECOMPOSITION OF THE QUARK SEA :
      SSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N - 20.* NS8N)/ 120.;
      DSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N
       - 30.* NS3N - 60.* DVN) / 120.;
      USN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N
       + 30.* NS3N - 60.* UVN) / 120.;
      CHN = (UVN + DVN + 2.*USN + 2.*DSN + 2.*SSN - NS15N)/6.;
      BON = (UVN + DVN + 2.*USN + 2.*DSN + 2.*SSN + 2.*CHN - NS24N)/8.;

      if (nf==3){                                     //GRV
         SSN= (20.* SIN - 20.* NS8N)/120.;
         DSN = (20.* SIN + 10.* NS8N - 30.* NS3N - 60.* DVN) / 120.;
         USN = (20.* SIN + 10.* NS8N + 30.* NS3N - 60.* UVN) / 120.;
         CHN=std::complex<double>(0.,0.);
         BON=std::complex<double>(0.,0.);
      }

      if (verbosity >= 50) {
	  std::cout << "ENS = " << ENS << " ALPr = " << ALPr << " XL1 = " << XL1 << " RMIN = " << RMIN << std::endl;
	  std::cout << "ANS = " << ANS << " S = " << S << std::endl;
      }
      if (verbosity >= 16) {	  
	  std::cout << "UVN = " << UVN << std::endl;
	  std::cout << "USN = " << USN << std::endl;
	  std::cout << "DVN = " << DVN << std::endl;
	  std::cout << "DSN = " << DSN << std::endl;
	  std::cout << "SSN = " << SSN << std::endl;
	  std::cout << "CHN = " << CHN << std::endl;
	  std::cout << "BON = " << BON << std::endl;
	  std::cout << "GLN = " << GLN << std::endl;
      }
      if (verbosity >= 50) {
	  std::cout << "UV evolution factor: " << ENS * (1.+  ALPr * XL1 * RMIN) << std::endl;
	  std::cout << "SIN = " << SIN << std::endl;
	  std::cout << "NS35N = " << NS35N << std::endl;
	  std::cout << "NS24N = " << NS24N << std::endl;
	  std::cout << "NS15N = " << NS15N << std::endl;
	  std::cout << "NS8N = " << NS8N << std::endl;
	  std::cout << "NS3N = " << NS3N << std::endl;
	  std::cout << "UVN = " << UVN << std::endl;
	  std::cout << "USN pieces: " << 10.* SIN << " " << 2.* NS35N << " " << 3.* NS24N << " " << 5.* NS15N << " " << 10.* NS8N << " " << 30.* NS3N << " " << -60.* UVN << std::endl;
      }


//...  OUTPUT
//  U   UB   D   DB   S   C    B    G
//  0    1   2   3    4   5    6    7
      FN[0] = UVN + USN;
      FN[1] = USN;
      FN[2] = DVN + DSN;
      FN[3] = DSN;
      FN[4] = SSN;
      FN[5] = CHN;
      FN[6] = BON;
      FN[7] = GLN;
      FN[8] = 0.;

}


/*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/

//
//...ANOMALOUS DIMENSIONS FOR LEADING AND NEXT TO LEADING ORDER
//...EVOLUTION OF PARTON DENSITIES AND WILSON COEFFICIENTS FOR
//...NLO STRUCTURE FUNCTIONS :
void ANOM (std::complex<double>& ANS, std::complex<double>& AM, std::complex<double>& AP, std::complex<double>& AL,
           std::complex<double>& BE, std::complex<double>& AB, std::complex<double>& RMIN, std::complex<double>& RPLUS,
           std::complex<double>& RQQ, std::complex<double>& RQG, std::complex<double>& RGQ, std::complex<double>& RGG,
           std::complex<double>& C2Q, std::complex<double>& C2G, std::complex<double>& CDYQ,
           std::complex<double>& CDYG, std::complex<double>& XN, int& FR, std::complex<double>& QQI,
           std::complex<double>& QGF, std::complex<double>& GQI, std::complex<double>& GGI, std::complex<double>& GGF,
           std::complex<double>& NS1MI, std::complex<double>& NS1PI, std::complex<double>& NS1F,
           std::complex<double>& QQ1F, std::complex<double>& QG1F, std::complex<double>& GQ1I,
           std::complex<double>& GQ1F, std::complex<double>& GG1I, std::complex<double>& GG1F,
           std::complex<double>& C2QI, std::complex<double>& C2GF, std::complex<double>& CDYQI,
           std::complex<double>& CDYGI, int verbosity, const int nf){

       double B0, B1, B02, B10;
       int N;
       double F;
       std::complex<double> QQ,QG,GQ,GG,SQ,GP,GM,NS1M,NS1P,QQ1,QG1,GQ1,GG1,DEL,XN1;
//...ANOMALOUS DIMENSIONS AND RELATED QUANTITIES IN LEADING ORDER :

// --->  COMMON/NF/nF

       if (verbosity >= 50) {
	   std::cout << "ANOM inputs: " << std::endl;
	   std::cout << "FR = " << FR << std::endl;
	   std::cout << "QQI = " << QQI << std::endl;
	   std::cout << "QGF = " << QGF  << std::endl;
	   std::cout << "GQI = " << GQI << std::endl;
	   std::cout << "GGI = " << GGI << std::endl;
	   std::cout << "GGF = " << GGF << std::endl;
	   std::cout << "NS1MI = " << NS1MI << std::endl;
	   std::cout << "NS1PI = " << NS1PI << std::endl;
	   std::cout << "NS1F = " << NS1F << std::endl;
	   std::cout << "QQ1F = " << QQ1F << std::endl;
	   std::cout << "QG1F = " << QG1F << std::endl;
	   std::cout << "GQ1I = " << GQ1I << std::endl;
	   std::cout << "GQ1F = " << GQ1F << std::endl;
	   std::cout << "GG1I = " << GG1I << std::endl;
	   std::cout << "GG1F = " << GG1F << std::endl;
       }

       F = FR;
       if (nf==3) F = 3;
       N = 1;
       B0 = 11.- 2./3.* FR;
       B02 = 2.* B0;
       QQ = QQI;
       QG = F * QGF;
       GQ = GQI;
       GG = GGI + F * GGF;
       SQ = std::sqrt ((GG - QQ) * (GG - QQ) + 4.* QG * GQ);
       GP = 0.5 * (QQ + GG + SQ);
       GM = 0.5 * (QQ + GG - SQ);
       ANS = QQ / B02;
       AM = GM / B02;
       AP = GP / B02;

       if (verbosity >= 50) {
	   std::cout << "GP = " << GP << " B02 = " << B02 << std::endl;
	   std::cout << "QQ = " << QQ << " GG = " << GG << " SQ = " << SQ << std::endl;
       }

       AL = (QQ - GP) / (GM - GP);
       BE = QG / (GM - GP);
       AB = GQ / (GM - GP);

       if (verbosity >= 50) {
	   std::cout << "B0 = " << B0 << std::endl;
	   std::cout << "QG = " << QG << std::endl;
	   std::cout << "GG = " << GG << std::endl;
	   std::cout << "SQ = " << SQ << std::endl;
	   std::cout << "GP = " << GP << std::endl;
	   std::cout << "GM = " << GM << std::endl;
	   std::cout << "ANS = " << ANS << std::endl;
	   std::cout << "AM = " << AM << std::endl;
	   std::cout << "AP = " << AP << std::endl;
	   std::cout << "AL = " << AL << std::endl;
	   std::cout << "BE = " << BE << std::endl;
	   std::cout << "AB = " << AB << std::endl;
       }


//...NEXT TO LEADING ORDER : ANOMALOUS DIMENSIONS AND WILSON COEFFICIENTS
//...IN THE MS-BAR FACTORIZATION SCHEME OF BARDEEN ET AL. (1981) :
       NS1M = NS1MI + F * NS1F;
       NS1P = NS1PI + F * NS1F;
       QQ1 = NS1P + F * QQ1F;
       QG1 = F * QG1F;
       GQ1 = GQ1I + F * GQ1F;
       GG1 = GG1I + F * GG1F;
       C2Q = C2QI;
       C2G = F * C2GF;
       CDYQ = CDYQI;
       CDYG = CDYGI;
//...CHANGE TO THE SCHEME OF ALTARELLI ET AL. (1979) FOR N = 2 OR TO THE
//...DIS SCHEME FOR N = 3 :
       if (N == 2){
          DEL = -0.5 * QG;
          C2G = C2G + DEL;
          QQ1 = QQ1 - DEL * GQ;
          QG1 = QG1 + DEL * (QQ - GG - B02 - QG);
          GQ1 = GQ1 + DEL * GQ;
          GG1 = GG1 + DEL * (GQ + B02);
       }
       else if (N == 3){
          NS1P = NS1P + B02 * C2Q;
          NS1M = NS1M + B02 * C2Q;
          QQ1 = QQ1 + C2Q * (QG + B02) + C2G * GQ;
          QG1 = QG1 + C2Q * QG + C2G * (GG - QQ + QG + B02);
          GQ1 = GQ1 - C2Q * (QQ - GG + GQ + B02) - C2G * GQ;
          GG1 = GG1 - C2Q * QG - C2G * (GQ + B02);
          C2Q = 0.;
          C2G = 0.;
       }
       XN1 = XN + 1.;
//       C3Q = C2Q - 8./3.* (1./ XN + 1./ XN1)
//       CLQ = 16./ (3.* XN1)
//       CLG = 8.* F / (XN1 * (XN + 2.))

       if (verbosity >= 50) {
	   std::cout << "NS1M = " << NS1M << std::endl;
	   std::cout << "NS1P = " << NS1P << std::endl;
	   std::cout << "QQ1 = " << QQ1 << std::endl;
	   std::cout << "QG1 = " << QG1 << std::endl;
	   std::cout << "GQ1 = " << GQ1 << std::endl;
	   std::cout << "GG1 = " << GG1 << std::endl;
	   std::cout << "C2Q = " << C2Q << std::endl;
	   std::cout << "C2G = " << C2G << std::endl;
	   std::cout << "CDYQ = " << CDYQ << std::endl;
	   std::cout << "CDYG = " << CDYG << std::endl;
       }

//...COMBINATIONS OF ANOMALOUS DIMENSIONS FOR NLO - SINGLET EVOLUTION :
       B1 = 102. - 38./3.* FR;
       B10 = B1 / B0;
       RMIN = (NS1M - QQ * B10) / B02;
       RPLUS = (NS1P - QQ * B10) / B02;
       RQQ = (QQ1 - QQ * B10) / B02;
       RQG = (QG1 - QG * B10) / B02;
       RGQ = (GQ1 - GQ * B10) / B02;
       RGG = (GG1 - GG * B10) / B02;

       if (verbosity >= 50) {
	   std::cout << "B1 = " << B1 << std::endl;
	   std::cout << "B10 = " << B10 << std::endl;
	   std::cout << "RMIN = " << RMIN << std::endl;
	   std::cout << "RPLUS = " << RPLUS << std::endl;
	   std::cout << "RQQ = " << RQQ << std::endl;
	   std::cout << "RQG = " << RQG << std::endl;
	   std::cout << "RGQ = " << RGQ << std::endl;
	   std::cout << "RGG = " << RGG << std::endl;
       }

}


/*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/

//... alphasl gives the LL/NLL evolution of alpha from
//    Qres^2 = (Q/a_param)^2
//    to
//    q2=bo^2/b^2
// alps=alpha(Qres)=alpQF

void alphaslcalc(double q2, double b0p, double a_param, double as, double mur2,
                 std::complex<double> nq2, resummationIndepfns* resufns,
                 std::complex<double>& alphasl, std::complex<double>& aexp, std::complex<double>& aexpB, int verbosity, 
                 int flag1){
//
      std::complex<double> xlambda,aa1,all,qq,t,xlt,bstar,b,blog,aa2;
      double blim,xlp,q;
      int iord;

      iord=flag1-1;
      if (iord<0) iord=0;

//.....Here computes NLL expression for alphas
//     here nq2=b0^2/b^2 and the result is now  alpha(nq2)/alpha(Qres)

//     HERE CHANGE: order of alphas related to order of evolution
      if(iord==1) xlp=1.;
      else if(iord==0) xlp=0.;

      q = std::sqrt(q2);

      int flagrealcomplex = 0, imod = 1;
      double beta0 = resufns->beta0, beta1 = resufns->beta1;

      b=std::pow(b0p*b0p/nq2, 0.5);
      blim=b0p*(1./q)*std::exp(1./(2.*as*beta0));                            // avoid Landau pole
      bstar=b;

//.....choose bstar (b) for real axis (complex plane) integration

      if (flagrealcomplex==0) bstar=b/sqrt(1. + (b*b)/(blim*blim));
      if (imod==1) blog=log( std::pow(q*bstar/b0p,2) + 1.);                  // modified sudakov
//      if (imod==0) blog= log( std::pow(q*bstar/b0p,2) );                     // normal sudakov

      // std::cout << "blog = " << blog << std::endl;
      // std::cout << "blim = " << blim << std::endl;
      // std::cout << "bstar = " << bstar << std::endl;

      xlambda=beta0*as*blog;

      if (xlambda.real() > 1) {
	  std::cout << "Warning: xlambda >1 means complex alphas as go near Landau pole - don't trust output" << std::endl;
      }

//     HERE now a dependence (without constant term)!
      aa1=log(1.-xlambda)+as*xlp*
        (beta1/beta0*log(1.-xlambda)/(1.-xlambda)
        + beta0*xlambda/(1.-xlambda)*log(q2/mur2)
//        +   beta0*log(q2/muf2)
        -2.*beta0*xlambda*log(a_param)/(1.-xlambda)   );
      alphasl=std::exp(-aa1);

      if (verbosity >= 50) {
	  std::cout << "aa1 = " << aa1 << std::endl;
	  std::cout << "xlambda = " << xlambda << std::endl;
	  std::cout << "log(1.-xlamda) = " << log(1.-xlambda) << std::endl;
	  std::cout << "as = " << as << std::endl;
	  std::cout << "xlp = " << xlp << std::endl;
	  std::cout << "beta1 = " << beta1 << std::endl;
	  std::cout << "beta0 = " << beta0 << std::endl;
	  std::cout << "q2 = " << q2 << std::endl;
	  std::cout << "mur2 = " << mur2 << std::endl;
	  std::cout << "a_param = " << a_param << std::endl;
	  std::cout << "blog = " << blog << std::endl;
	  std::cout << "q = " << q << std::endl;
	  std::cout << "bstar = " << bstar << std::endl;
	  std::cout << "b0p = " << b0p << std::endl;
	  std::cout << "b = " << b << std::endl;
	  std::cout << "nq2 = " << nq2 << std::endl;
	  std::cout << "blim = " << blim << std::endl;
      }

//Now compute the factors
//needed to resum the logs which multiply the N-dependent part
//of the C coefficients
//the limit below implies xlambda<1/2 and then aa2<= 1
//blim=b0p*(1/q)*exp(1/(2*as*beta0))
//blim=b0p*(1/q)*exp(1/(4*as*beta0))
//Set a limit to avoid very large values of b (= very small scales ~1/b)
      blim=0.5;

      if (flagrealcomplex==0) bstar=b/sqrt(1. + (b*b)/(blim*blim));
      if (imod==1) blog=log( std::pow(q*bstar/b0p,2) + 1.);                  // modified sudakov
//      if (imod==0) blog= log( std::pow(q*bstar/b0p,2) );                     // normal sudakov

      xlambda=beta0*as*blog;

      if (verbosity > 15) {
	  std::cout << "blog = " << blog << std::endl;
	  std::cout << "bstar = " << bstar << std::endl;
	  std::cout << "blim = " << blim << std::endl;
      }

      aa1=log(1.-xlambda)+as*xlp*(beta1/beta0*log(1.-xlambda)/(1.-xlambda) + beta0*xlambda/(1.-xlambda)*log(q2/mur2)-2.0*beta0*xlambda*log(a_param)/(1.-xlambda));
      if (verbosity > 50) {
	  std::cout << "aa1 = " << aa1 << std::endl;
	  std::cout << "xlambda = " << xlambda << std::endl;
	  std::cout << "beta0 = " << beta0 << " as = " << as << " blog = " << blog << std::endl;
	  std::cout << "xlp = " << xlp << std::endl;
      }

      aexp = std::exp(-aa1);
	  
      aa2 = xlambda/(1. - xlambda);
      aexpB = std::exp(aa2);

      if (verbosity > 15) {
	  std::cout << "aa1 = " << aa1 << std::endl;
	  std::cout << "aexp = " << aexp << std::endl;
	  std::cout << "aa2 = " << aa2 << std::endl;
	  std::cout << "aexpb = " << aexpB << std::endl;
      }

}
