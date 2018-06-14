//Calculates Hard functions for inverse mellin transform

#include "hardfns.h"

#include "constants.h"
#include "resummation_input.h"
#include "resu_preproc.h"
#include "resu_PS.h"
#include "evolution.h"
#include "resu_procindep.h"


//DEBUG
#include <chrono>
extern int k_count, k_count0, k_count2;
extern std::chrono::high_resolution_clock::time_point start, end_;


//Calculates Hard functions for inverse mellin transform argument, sudakq and sudakg are passed in on end for cases where sudakovs begin very small (typically less than 1e-300 so stored as 0) and hard factors very large (typically more than 1e300 so stored as inf and then get nans).

void Hardfns_calc (int qq_order, int gg_order, int pcF, int em_charge, resu_PS* resu, resummationIndepfns* indep, AllDep_N PR1, AllDep_N PR2, resuDep_N N1, resuDep_N N2, std::complex<double> alpq, std::complex<double> aexp, std::complex<double> aexpb, std::complex<double> sudakq, std::complex<double> sudakg, std::complex<double>& HCRNqq, std::complex<double>& HCRNgg){

//NOTE - all if statements commented out here as they add significant time to the running of this program. They can be uncommented for debugging. Hardfns_calc is called up to 136 (number of points on mellin contour) times per b point, and there are typically 20 (although can be lot more) b points per invbtoqt inverse fourier transform and invbtoqt is called for each of the Monte Carlo iterations - therefore for 100000 iterations then Hardfns_calc may be called 100 million times, therefore this is the time critical part of the code, it currently takes approx 10^-5 seconds, so 100000 iterations should take 1000s-10000s of seconds of cpu_time (i.e. per core), of course all this is very dependent on the computer used.

/*
// DEBUG
if((k_count0==407||k_count0==408) && k_count==2 && k_count2==10){
  end_ = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//  std::cout<<"hardfns start chkpoint: "<< time_span.count() <<'\n';
}
*/

  // std::cout << "em_charge in hardfns.cc = " << em_charge << std::endl;
  
    double pi = k_constants::pi;
    int Nf = (PR1.FP.size()-1)/2; // Get Nf from PDF vector size (2*Nf+1)

    std::complex<double> OneComplex = std::complex<double>(1.,0.);
    std::complex<double> zero_C = std::complex<double>(0.,0.);
    std::complex<double> cqg1, cqg2;
    cqg1 = alpq*2.0*OneComplex*N1.C1qg;
    cqg2 = alpq*2.0*OneComplex*N2.C1qg;

    double as2 = resu->alphas/2.0;
    double loga = std::log(resu->a);
    double logmuf = std::log(resu->muf2/resu->q2);
    double logmur = std::log(resu->mur2/resu->q2);

    std::complex<double> H1stgg, H1stqg_1, H1stqg_2, H1stqqb;
    std::complex<double> H2stgg, H2stqg_1, H2stqg_2, H2stqq_1, H2stqq_2;
    std::complex<double> H2stqqp_1, H2stqqp_2;
    std::complex<double> H2stqqb[5], Hqqb[5];
    std::complex<double> Hgg, Hqg_1, Hqg_2, Hgq_1, Hgq_2, Hqq, Hqq_1, Hqq_2, Hqqp_1, Hqqp_2;
    H1stgg = zero_C;
    H1stqg_1 = zero_C;
    H1stqg_2 = zero_C;
    H1stqqb = zero_C;
    H2stgg = zero_C;
    H2stqg_1 = zero_C;
    H2stqg_2 = zero_C;
    H2stqq_1 = zero_C;
    H2stqq_2 = zero_C;
    H2stqqp_1 = zero_C;
    H2stqqp_2 = zero_C;
    Hgg = zero_C;
    Hqg_1 = zero_C;
    Hqg_2 = zero_C;
    Hgq_1 = zero_C;
    Hgq_2 = zero_C;
    Hqq = zero_C;
    Hqq_1 = zero_C;
    Hqq_2 = zero_C;
    Hqqp_1 = zero_C;
    Hqqp_2 = zero_C;
    for(int i = 0; i<5; i++) {
      H2stqqb[i] = zero_C;
      Hqqb[i] = zero_C;
    }

    HCRNqq = zero_C, HCRNgg = zero_C;


    //Quark contribution - final state in Hst is always Q Qb or Qb Q here
    //Leading log - trivial but useful for testing
    if (qq_order == 0) {
      for (int i = 0; i<5; i++) {Hqqb[i] = OneComplex;}
// All non-diagonal channels are 0 at LL
//End LL
    }

//NLL
    else if (qq_order == 1) {

      for (int i = 0; i<5; i++) {
        Hqqb[i] = OneComplex + resu->alphas/2.0*(resu->H1q + N1.C1qq + N2.C1qq) - resu->alphas/2.0*(N1.gamma1qq + N2.gamma1qq)*(logmuf+2*loga) + resu->alphas/2.0*(-4*loga)*(indep->B1q + indep->A1q*loga) + as2*2.*pcF*indep->beta0*logmur;
      }
      Hqg_1 = cqg1 + as2*(-N1.gamma1qg)*(logmuf + 2.0*loga);
      Hqg_2 = cqg2 + as2*(-N2.gamma1qg)*(logmuf + 2.0*loga);

// all other channels are still = 0 at NLL
      Hgg = 0.0;
      Hqq_1 = 0.0;
      Hqq_2 = 0.0;
      Hqq= 0.0;
      Hqqp_1 = 0.0;
      Hqqp_2 = 0.0;

    }



//NNLL
    else if (qq_order == 2) {

// DEBUG
if((k_count0==407||k_count0==408) && k_count==2 && k_count2==10){
  end_ = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//  std::cout<<"hardfns 0- chkpoint: "<< time_span.count() <<'\n';
}

      H1stqqb = (resu->H1q + N1.C1qq + N2.C1qq)-(N1.gamma1qq + N2.gamma1qq)*(logmuf+2*loga)+(-4*loga)*(indep->B1q+indep->A1q*loga)+ 2.*pcF*indep->beta0*logmur;
      H1stqg_1 = N1.C1qg + (-N1.gamma1qg)*(logmuf+2*loga);
      H1stqg_2 = N2.C1qg + (-N2.gamma1qg)*(logmuf+2*loga);
      H1stgg = 0.0;

// All *4 because of normalization (as/2pi)^2 in stead of as/pi, normalization of as/2pi  implies  gamma^1->gamma^1/2, C^1, H^1 -> C^1/2, H^1/2, gamma^2-> gamma^2/4, beta,A,B -> beta,A,B are in as/pi already

      H2stqq_1 = N1.C2NSqqb + N1.C2Sqqb;
      H2stqq_2 = N2.C2NSqqb + N2.C2Sqqb;

      H2stqqp_1 = N1.C2Sqqb;
      H2stqqp_2 = N2.C2Sqqb;

      for (int l = 0; l<5; l++) {
        H2stqqb[l] = resu->H2q[l] + N1.C2NSqq + N2.C2NSqq + N1.C2Sqqb + N2.C2Sqqb + N1.C1qq*N2.C1qq + resu->H1q*(N1.C1qq+N2.C1qq);
        H2stqqb[l] = H2stqqb[l] + 4.0*(indep->A1q*indep->beta0*8.0/6.0*loga*loga*loga + 2.0*loga*loga*(indep->A2q-indep->beta0*(indep->B1q+2.0*indep->A1q*loga + N1.gamma1qq/2.0 + N2.gamma1qq/2.0)) -2.0*loga*(indep->B2q + 2.0*indep->A2q*loga - indep->beta0*(N1.C1qq + N2.C1qq)/2.0 + N1.gamma2qqV/4.0 + N1.gamma2qqS/4.0 + N2.gamma2qqV/4.0 + N2.gamma2qqS/4.0) + indep->beta0/2.0*(N1.gamma1qq + N2.gamma1qq)/2.0*logmuf*logmuf - (N1.gamma2qqV/4.0 + N1.gamma2qqS/4.0 +
        N2.gamma2qqV/4.0 + N2.gamma2qqS/4.0)*logmuf + H1stqqb/2.0*indep->beta0*logmur + 1.0/2.0*(H1stqqb + resu->H1q + N1.C1qq + N2.C1qq)/2.0*((N1.gamma1qq + N2.gamma1qq)/2.0*(-logmuf - 2.0*loga) - ((indep->B1q+indep->A1q*loga)*2.0*loga)) + 1.0/4.0*(H1stqg_1+N1.C1qg)*N1.gamma1gq/2.0*(-logmuf-2.0*loga) + 1.0/4.0*(H1stqg_2+N2.C1qg)*N2.gamma1gq/2.0*(-logmuf-2.0*loga))
        - 4.*pcF*(1./2.*std::pow(indep->beta0*logmur,2) - indep->beta1*logmur) ;
      }

      H2stqg_1 = N1.C2qg + N1.C1qg*N2.C1qq + 4.0*(indep->beta0*2.0*loga*loga*(-N1.gamma1qg/2.0) - 2.0*loga*(-indep->beta0*N1.C1qg/2.0 + N1.gamma2qg/4.0) + indep->beta0/2.0*logmuf*logmuf*(N1.gamma1qg/2.0) - N1.gamma2qg/4.0*logmuf + indep->beta0*logmur*H1stqg_1/2.0 + 0.5*(H1stqqb + resu->H1q + N1.C1qq + N2.C1qq)/2.0*(-logmuf-2.0*loga)*N1.gamma1qg/2.0 + 0.5*(H1stqg_1 + N1.C1qg)/2.0*((-logmuf-2.0*loga)*(N2.gamma1qq + N1.gamma1gg)/2.0 - ((indep->B1q+indep->A1q*loga)*2.0*loga)));

      H2stqg_2 = N2.C2qg + N2.C1qg*N1.C1qq + 4.0*(indep->beta0*2.0*loga*loga*(-N2.gamma1qg/2.0) - 2.0*loga*(-indep->beta0*N2.C1qg/2.0 + N2.gamma2qg/4.0) + indep->beta0/2.0*logmuf*logmuf*(N2.gamma1qg/2.0) - N2.gamma2qg/4.0*logmuf + indep->beta0*logmur*H1stqg_2/2.0 + 0.5*(H1stqqb + resu->H1q + N2.C1qq + N1.C1qq)/2.0*(-logmuf-2.0*loga)*N2.gamma1qg/2.0 + 0.5*(H1stqg_2 + N2.C1qg)/2.0*((-logmuf-2.0*loga)*(N1.gamma1qq + N2.gamma1gg)/2.0 - ((indep->B1q+indep->A1q*loga)*2.0*loga)));

//GG
      H2stgg = N2.C1qg*N1.C1qg-4.0*(0.5*(logmuf + 2*loga))*((H1stqg_1 + N1.C1qg)/2.0*N2.gamma1qg/2.0 + (H1stqg_2 + N2.C1qg)/2.0*N1.gamma1qg/2.0);

//qq_1 i.e. QQ -> Qb Q = Qb Qb -> Q Qb
      H2stqq_1 = H2stqq_1 + 4.0*(-2.0*loga*((0.25*(N1.gamma2qqbV + N1.gamma2qqbS))) - (0.25*(N1.gamma2qqbV + N1.gamma2qqbS))*logmuf + 0.5*(H1stqg_1 + N1.C1qg)/2.0*(0.5*(N1.gamma1gq)*(-logmuf - 2.0*loga)));

//qq_2 i.e. Qb Qb -> Qb Q = Q Q -> Q Qb
      H2stqq_2 = H2stqq_2 + 4.0*(-2.0*loga*((0.25*(N2.gamma2qqbV+N2.gamma2qqbS))) - (0.25*(N2.gamma2qqbV + N2.gamma2qqbS))*logmuf + 0.5*(H1stqg_2 + N2.C1qg)/2.0*(0.5*(N2.gamma1gq)*(-logmuf-2.0*loga)));

//qqp_1 i.e. Q' Q -> Qb Q flavour in sigmaQQb determined by second parton
      H2stqqp_1 = H2stqqp_1 + 4.0*(-2.0*loga*(0.25*N1.gamma2qqbS) - (0.25*N1.gamma2qqbS)*logmuf + 0.5*(H1stqg_1 + N1.C1qg)/2.0*(0.5*N1.gamma1gq*(-logmuf -2.0*loga)));

//qqp_2 i.e. Q Q; -> Q Qb flavour in sigmaQQb determined by first parton
      H2stqqp_2 = H2stqqp_2 + 4.0*(-2.0*loga*(0.25*N2.gamma2qqbS) - 0.25*N2.gamma2qqbS*logmuf + 0.5*(H1stqg_2 + N2.C1qg)/2.0*(0.5*N2.gamma1gq*(-logmuf-2.0*loga)));

// DEBUG
if((k_count0==407||k_count0==408) && k_count==2 && k_count2==10){
  end_ = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//  std::cout<<"hardfns 1 chkpoint: "<< time_span.count() <<'\n';
}

// std::cout << "H1stqqb = " << H1stqqb << std::endl;
// std::cout << "H1stqg_1 = " << H1stqg_1 << std::endl;
// std::cout << "H1stqg_2 = " << H1stqg_2 << std::endl;
// std::cout << "H1stgg = " << H1stgg << std::endl;
// std::cout << "H2stqqb = " << H2stqqb[0] << " " << H2stqqb[1] << " " <<H2stqqb[2] << " " <<H2stqqb[3] << " " <<H2stqqb[4] << " " << std::endl;
// std::cout << "H2stqg_1 = " << H2stqg_1 << std::endl;
// std::cout << "H2stqg_2 = " << H2stqg_2 << std::endl;
// std::cout << "H2stgg = " << H2stgg << std::endl;
// std::cout << "H2stqq_1 = " << H2stqq_1 << std::endl;
// std::cout << "H2stqq_2 = " << H2stqq_2 << std::endl;
// std::cout << "H2stqqp_1 = " << H2stqqp_1 << std::endl;
// std::cout << "H2stqqp_2 = " << H2stqqp_2 << std::endl;
 
// Factor in N-dependent part of resummation


     std::complex<double> leg1qq_exp = std::pow(aexpb,resu->alphas/2.0*(N1.C1qq - 2.0*pi*pi/3.0 + 16.0/3.0));
     std::complex<double> leg2qq_exp = std::pow(aexpb,resu->alphas/2.0*(N2.C1qq - 2.0*pi*pi/3.0 + 16.0/3.0));

     std::complex<double> leg1qg_exp = std::pow(aexpb,resu->alphas/2.0*(N1.C2qg/N1.C1qg-2.0*pi*pi/3.0+16.0/3.0));
     std::complex<double> leg2qg_exp = std::pow(aexpb,resu->alphas/2.0*(N2.C2qg/N2.C1qg-2.0*pi*pi/3.0+16.0/3.0));

//QQb
//first pow is q leg, second pow is qb leg
      for (int l = 0; l<5; l++) {
        Hqqb[l] = (1.0 + as2*H1stqqb + as2*as2*H2stqqb[l])*leg1qq_exp*leg2qq_exp;
      }

//qg_1 i.e. gq initial state
      Hqg_1 = (as2*H1stqg_1 + as2*as2*H2stqg_1)*aexp*leg1qg_exp*leg1qq_exp;

//qg_2 i.e. qg initial state
      Hqg_2 = (as2*H1stqg_2 + as2*as2*H2stqg_2)*aexp*leg2qg_exp*leg2qq_exp;

//GG
      Hgg = as2*as2*H2stgg*aexp*leg1qg_exp*aexp*leg2qg_exp;
      // std::cout << "as2 = " << as2 << std::endl;
      // std::cout << "H2stgg = " << H2stgg << std::endl;
      // std::cout << "aexp = " << aexp << std::endl;
      // std::cout << "leg1qg_exp = " << leg1qg_exp << std::endl;
      // std::cout << "leg2qg_exp = " << leg2qg_exp << std::endl;
      // std::cout << "Hgg = " << Hgg << std::endl;
      
//QQ -> QbQ = QbQb -> Q Qb
      Hqq_1 = (as2*as2*H2stqq_1)*leg1qq_exp*std::pow(aexp,2);

//QbQb -> QbQ = QQ -> Q Qb
      Hqq_2 = (as2*as2*H2stqq_2)*leg2qq_exp*std::pow(aexp,2);

//Average QQ -> QQb and QbQb -> QQb
      Hqq = (Hqq_1 + Hqq_2)/2.0;

//qqp_1 i.e. Q' Q -> Qb Q flavour in sigmaQQb determined by second parton
      Hqqp_1 = (as2*as2*H2stqqp_1)*leg1qq_exp*aexp*aexp; // aexpb power is q leg, aexp power is qb leg

//qqp_2 i.e. Q Q' -> Q Qb flavour in sigmaQQb determined by first parton
      Hqqp_2 = (as2*as2*H2stqqp_2)*leg2qq_exp*aexp*aexp; // aexpb power is q leg, aexp power is qb leg

    }
/*
// DEBUG
if((k_count0==407||k_count0==408) && k_count==2 && k_count2==10){
  end_ = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//  std::cout<<"hardfns 2 chkpoint: "<< time_span.count() <<'\n';
}
*/

//Final: combine pieces

    std::complex<double> GGN, QGN_1, QGN_2, QQBN;
    GGN = zero_C;
    QGN_1 = zero_C;
    QGN_2 = zero_C;
    QQBN = zero_C;

    int m = 0;

    //Sums for diphoton or Z drellyan case as here only symmetric quarks are linked by the boson

    if (em_charge == 0 && qq_order >= 0) {

      for (int n = 0; n<2.0*Nf+1.0; n++) {
        if (n != 5) {
        if (n < 5) {m = 4 - n;}
        else {m = n - 6;}
//GG
        GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[n][10-n]*Hgg;
	// std::cout << "PR1.FP[5] = " << PR1.FP[5] << std::endl;
	// std::cout << "PR2.FP[5] = " << PR2.FP[5] << std::endl;
	// std::cout << "resu->sigmaij[" << n << "][" << 10-n << "] = " << resu->sigmaij[n][10-n] << std::endl;
	// std::cout << "Hgg = " << Hgg << std::endl;
	// std::cout << "sudakq = " << sudakq << std::endl;
//GQ + GQb
		QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[10-n]*resu->sigmaij[n][10-n]*Hqg_1;
//QG + QbG
		QGN_2 = QGN_2 + PR1.FP[n]*sudakq*PR2.FP[5]*resu->sigmaij[n][10-n]*Hqg_2;
//QQb + QbQ
		QQBN = QQBN + PR1.FP[n]*sudakq*PR2.FP[10-n]*resu->sigmaij[n][10-n]*Hqqb[m];
//QQ + QbQb
		QQBN = QQBN + PR1.FP[n]*sudakq*PR2.FP[n]*resu->sigmaij[n][10-n]*Hqq;
		for (int p = 0; p<2.0*Nf+1.0; p++) {
		    int q = 0;
		    q = 10 - n;
		    if (p != n && p != 5 && p != q) {
			QQBN = QQBN + PR1.FP[p]*sudakq*PR2.FP[10-n]*resu->sigmaij[n][10-n]*Hqqp_1;
			QQBN = QQBN + PR1.FP[n]*sudakq*PR2.FP[p]*resu->sigmaij[n][10-n]*Hqqp_2;
			// std::cout << "resu->sigmaij[n][10-n] = " << resu->sigmaij[n][10-n] << std::endl;
		    }
		}
	    }
	}
    }

    else if (em_charge == 1) {
      // W+ so can only go to udbar, usbar, ubbar, cdbar, csbar, cbbar, dbaru, sbaru, bbaru, dbarc, sbarc, bbarc
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[6][3]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[6][2]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[6][0]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[9][3]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[9][2]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[9][0]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[3][6]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[2][6]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[0][6]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[3][9]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[2][9]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[0][9]*Hgg;

      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[3]*resu->sigmaij[6][3]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[2]*resu->sigmaij[6][2]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[0]*resu->sigmaij[6][0]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[3]*resu->sigmaij[9][3]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[2]*resu->sigmaij[9][3]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[0]*resu->sigmaij[9][3]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[6]*resu->sigmaij[3][6]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[6]*resu->sigmaij[2][6]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[6]*resu->sigmaij[0][6]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[9]*resu->sigmaij[3][9]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[9]*resu->sigmaij[2][9]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[9]*resu->sigmaij[0][9]*Hqg_1;

      QGN_2 = QGN_2 + PR1.FP[6]*sudakq*PR2.FP[5]*resu->sigmaij[6][3]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[6]*sudakq*PR2.FP[5]*resu->sigmaij[6][2]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[6]*sudakq*PR2.FP[5]*resu->sigmaij[6][0]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[9]*sudakq*PR2.FP[5]*resu->sigmaij[9][3]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[9]*sudakq*PR2.FP[5]*resu->sigmaij[9][2]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[9]*sudakq*PR2.FP[5]*resu->sigmaij[9][0]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[3]*sudakq*PR2.FP[5]*resu->sigmaij[3][6]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[2]*sudakq*PR2.FP[5]*resu->sigmaij[2][6]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[0]*sudakq*PR2.FP[5]*resu->sigmaij[0][6]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[3]*sudakq*PR2.FP[5]*resu->sigmaij[3][9]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[2]*sudakq*PR2.FP[5]*resu->sigmaij[2][9]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[0]*sudakq*PR2.FP[5]*resu->sigmaij[0][9]*Hqg_2;

      QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[3]*resu->sigmaij[6][3]*Hqqb[0];
      QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[2]*resu->sigmaij[6][2]*Hqqb[0];
      QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[0]*resu->sigmaij[6][0]*Hqqb[0];
      QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[3]*resu->sigmaij[9][3]*Hqqb[3];
      QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[2]*resu->sigmaij[9][2]*Hqqb[3];
      QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[0]*resu->sigmaij[9][0]*Hqqb[3];
      QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[6]*resu->sigmaij[3][6]*Hqqb[1];
      QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[6]*resu->sigmaij[2][6]*Hqqb[2];
      QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[6]*resu->sigmaij[0][6]*Hqqb[4];
      QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[9]*resu->sigmaij[3][9]*Hqqb[1];
      QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[9]*resu->sigmaij[2][9]*Hqqb[2];
      QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[9]*resu->sigmaij[0][9]*Hqqb[4];

      QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[7]*resu->sigmaij[6][3]*Hqqb[0];
      QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[8]*resu->sigmaij[6][2]*Hqqb[0];
      QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[10]*resu->sigmaij[6][0]*Hqqb[0];
      QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[7]*resu->sigmaij[9][3]*Hqqb[3];
      QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[8]*resu->sigmaij[9][2]*Hqqb[3];
      QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[10]*resu->sigmaij[9][0]*Hqqb[3];
      QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[4]*resu->sigmaij[3][6]*Hqqb[1];
      QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[4]*resu->sigmaij[2][6]*Hqqb[2];
      QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[4]*resu->sigmaij[0][6]*Hqqb[4];
      QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[1]*resu->sigmaij[3][9]*Hqqb[1];
      QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[1]*resu->sigmaij[2][9]*Hqqb[2];
      QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[1]*resu->sigmaij[0][9]*Hqqb[4];

      for (int qk = 0; qk < 2.0*Nf+1.0;qk++) {
	if (qk != 5) {
	  if (qk != 6 && qk != 4) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[3]*resu->sigmaij[6][3]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[2]*resu->sigmaij[6][2]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[0]*resu->sigmaij[6][0]*Hqqp_1;
	  }
	  if (qk != 9 && qk !=1) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[3]*resu->sigmaij[9][3]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[2]*resu->sigmaij[9][3]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[0]*resu->sigmaij[9][3]*Hqqp_1;
	  }
	  if (qk != 3 && qk != 7) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[6]*resu->sigmaij[3][6]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[9]*resu->sigmaij[3][9]*Hqqp_1;
	  }
	  if (qk != 2 && qk!= 8) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[6]*resu->sigmaij[2][6]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[9]*resu->sigmaij[2][9]*Hqqp_1;
	  }
	  if (qk != 0 && qk!= 10) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[6]*resu->sigmaij[0][6]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[9]*resu->sigmaij[0][9]*Hqqp_1;
	  }
	  if (qk != 3 && qk !=7) {
	    QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[qk]*resu->sigmaij[6][3]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[qk]*resu->sigmaij[9][3]*Hqqp_2;
	  }
	  if( qk != 2 && qk != 8) {
	    QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[qk]*resu->sigmaij[6][2]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[qk]*resu->sigmaij[9][2]*Hqqp_2;
	  }
	  if (qk != 0 && qk!= 10) {
	    QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[qk]*resu->sigmaij[6][0]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[qk]*resu->sigmaij[9][0]*Hqqp_2;
	  }
	  if (qk != 6 && qk != 4) {
	    QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[qk]*resu->sigmaij[3][6]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[qk]*resu->sigmaij[2][6]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[qk]*resu->sigmaij[0][6]*Hqqp_2;
	  }
	  if (qk != 9 && qk != 1) {
	    QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[qk]*resu->sigmaij[3][9]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[qk]*resu->sigmaij[2][9]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[qk]*resu->sigmaij[0][9]*Hqqp_2;
	  }
	}
      }
    }
    else if (em_charge == -1) {
      //W- so can only go to ubard, ubars, ubarb, cbard, cbars, cbarb, dubar, subar, bubar, dcbar, scbar, bcbar
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[4][7]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[4][8]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[4][10]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[1][7]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[1][8]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[1][10]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[7][4]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[8][4]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[10][4]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[7][1]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[8][1]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[10][1]*Hgg;

      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[7]*resu->sigmaij[4][7]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[8]*resu->sigmaij[4][8]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[10]*resu->sigmaij[4][10]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[7]*resu->sigmaij[1][7]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[8]*resu->sigmaij[1][8]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[10]*resu->sigmaij[1][10]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[4]*resu->sigmaij[7][4]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[4]*resu->sigmaij[8][4]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[4]*resu->sigmaij[10][4]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[1]*resu->sigmaij[7][1]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[1]*resu->sigmaij[8][1]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[1]*resu->sigmaij[10][1]*Hqg_1;

      QGN_2 = QGN_2 + PR1.FP[4]*sudakq*PR2.FP[5]*resu->sigmaij[4][7]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[4]*sudakq*PR2.FP[5]*resu->sigmaij[4][8]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[4]*sudakq*PR2.FP[5]*resu->sigmaij[4][10]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[1]*sudakq*PR2.FP[5]*resu->sigmaij[1][7]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[1]*sudakq*PR2.FP[5]*resu->sigmaij[1][8]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[1]*sudakq*PR2.FP[5]*resu->sigmaij[1][10]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[7]*sudakq*PR2.FP[5]*resu->sigmaij[7][4]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[8]*sudakq*PR2.FP[5]*resu->sigmaij[8][4]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[10]*sudakq*PR2.FP[5]*resu->sigmaij[10][4]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[7]*sudakq*PR2.FP[5]*resu->sigmaij[7][1]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[8]*sudakq*PR2.FP[5]*resu->sigmaij[8][1]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[10]*sudakq*PR2.FP[5]*resu->sigmaij[10][1]*Hqg_2;

      QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[7]*resu->sigmaij[4][7]*Hqqb[0];
      QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[8]*resu->sigmaij[4][8]*Hqqb[0];
      QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[10]*resu->sigmaij[4][10]*Hqqb[0];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[7]*resu->sigmaij[1][7]*Hqqb[3];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[8]*resu->sigmaij[1][8]*Hqqb[3];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[10]*resu->sigmaij[1][10]*Hqqb[3];
      QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[4]*resu->sigmaij[7][4]*Hqqb[1];
      QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[4]*resu->sigmaij[8][4]*Hqqb[2];
      QQBN = QQBN + PR1.FP[10]*sudakq*PR2.FP[4]*resu->sigmaij[10][4]*Hqqb[4];
      QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[1]*resu->sigmaij[7][1]*Hqqb[1];
      QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[1]*resu->sigmaij[8][1]*Hqqb[2];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[1]*resu->sigmaij[10][1]*Hqqb[4];

      QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[3]*resu->sigmaij[4][7]*Hqqb[0];
      QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[2]*resu->sigmaij[4][8]*Hqqb[0];
      QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[0]*resu->sigmaij[4][10]*Hqqb[0];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[3]*resu->sigmaij[1][7]*Hqqb[3];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[2]*resu->sigmaij[1][8]*Hqqb[3];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[0]*resu->sigmaij[1][10]*Hqqb[3];
      QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[6]*resu->sigmaij[7][4]*Hqqb[1];
      QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[6]*resu->sigmaij[8][4]*Hqqb[2];
      QQBN = QQBN + PR1.FP[10]*sudakq*PR2.FP[6]*resu->sigmaij[10][4]*Hqqb[4];
      QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[9]*resu->sigmaij[7][1]*Hqqb[1];
      QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[9]*resu->sigmaij[8][1]*Hqqb[2];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[9]*resu->sigmaij[10][1]*Hqqb[4];


      for (int qk = 0; qk < 2.0*Nf+1.0;qk++) {
	if (qk != 5) {
	  if (qk != 6 && qk != 4) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[7]*resu->sigmaij[4][7]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[8]*resu->sigmaij[4][8]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[10]*resu->sigmaij[4][10]*Hqqp_1;
	  }
	  if (qk != 9 && qk !=1) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[7]*resu->sigmaij[1][7]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[8]*resu->sigmaij[1][8]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[10]*resu->sigmaij[1][10]*Hqqp_1;
	  }
	  if (qk != 3 && qk != 7) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[4]*resu->sigmaij[7][4]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[1]*resu->sigmaij[7][1]*Hqqp_1;
	  }
	  if (qk != 2 && qk!= 8) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[4]*resu->sigmaij[8][4]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[4]*resu->sigmaij[8][1]*Hqqp_1;
	  }
	  if (qk != 0 && qk!= 10) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[4]*resu->sigmaij[10][4]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[1]*resu->sigmaij[10][1]*Hqqp_1;
	  }
	  if (qk != 3 && qk !=7) {
	    QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[qk]*resu->sigmaij[4][7]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[qk]*resu->sigmaij[1][7]*Hqqp_2;
	  }
	  if( qk != 2 && qk != 8) {
	    QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[qk]*resu->sigmaij[4][8]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[qk]*resu->sigmaij[1][8]*Hqqp_2;
	  }
	  if (qk != 0 && qk!= 10) {
	    QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[qk]*resu->sigmaij[4][10]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[qk]*resu->sigmaij[1][10]*Hqqp_2;
	  }
	  if (qk != 6 && qk != 4) {
	    QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[qk]*resu->sigmaij[7][4]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[qk]*resu->sigmaij[8][4]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[10]*sudakq*PR2.FP[qk]*resu->sigmaij[10][4]*Hqqp_2;
	  }
	  if (qk != 9 && qk != 1) {
	    QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[qk]*resu->sigmaij[7][1]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[qk]*resu->sigmaij[8][1]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[10]*sudakq*PR2.FP[qk]*resu->sigmaij[10][1]*Hqqp_2;
	  }
	}
      }
    }
    else if (em_charge == 11) {
      //W+ and W-
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[6][3]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[6][2]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[6][0]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[9][3]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[9][2]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[9][0]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[3][6]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[2][6]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[0][6]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[3][9]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[2][9]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[0][9]*Hgg;

      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[3]*resu->sigmaij[6][3]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[2]*resu->sigmaij[6][2]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[0]*resu->sigmaij[6][0]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[3]*resu->sigmaij[9][3]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[2]*resu->sigmaij[9][3]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[0]*resu->sigmaij[9][3]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[6]*resu->sigmaij[3][6]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[6]*resu->sigmaij[2][6]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[6]*resu->sigmaij[0][6]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[9]*resu->sigmaij[3][9]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[9]*resu->sigmaij[2][9]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[9]*resu->sigmaij[0][9]*Hqg_1;

      QGN_2 = QGN_2 + PR1.FP[6]*sudakq*PR2.FP[5]*resu->sigmaij[6][3]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[6]*sudakq*PR2.FP[5]*resu->sigmaij[6][2]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[6]*sudakq*PR2.FP[5]*resu->sigmaij[6][0]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[9]*sudakq*PR2.FP[5]*resu->sigmaij[9][3]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[9]*sudakq*PR2.FP[5]*resu->sigmaij[9][2]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[9]*sudakq*PR2.FP[5]*resu->sigmaij[9][0]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[3]*sudakq*PR2.FP[5]*resu->sigmaij[3][6]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[2]*sudakq*PR2.FP[5]*resu->sigmaij[2][6]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[0]*sudakq*PR2.FP[5]*resu->sigmaij[0][6]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[3]*sudakq*PR2.FP[5]*resu->sigmaij[3][9]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[2]*sudakq*PR2.FP[5]*resu->sigmaij[2][9]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[0]*sudakq*PR2.FP[5]*resu->sigmaij[0][9]*Hqg_2;

      QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[3]*resu->sigmaij[6][3]*Hqqb[0];
      QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[2]*resu->sigmaij[6][2]*Hqqb[0];
      QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[0]*resu->sigmaij[6][0]*Hqqb[0];
      QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[3]*resu->sigmaij[9][3]*Hqqb[3];
      QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[2]*resu->sigmaij[9][2]*Hqqb[3];
      QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[0]*resu->sigmaij[9][0]*Hqqb[3];
      QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[6]*resu->sigmaij[3][6]*Hqqb[1];
      QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[6]*resu->sigmaij[2][6]*Hqqb[2];
      QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[6]*resu->sigmaij[0][6]*Hqqb[4];
      QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[9]*resu->sigmaij[3][9]*Hqqb[1];
      QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[9]*resu->sigmaij[2][9]*Hqqb[2];
      QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[9]*resu->sigmaij[0][9]*Hqqb[4];

      for (int qk = 0; qk < 2.0*Nf+1.0;qk++) {
        if (qk != 5) {
	  if (qk != 6 && qk != 4) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[3]*resu->sigmaij[6][3]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[2]*resu->sigmaij[6][2]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[0]*resu->sigmaij[6][0]*Hqqp_1;
	  }
	  if (qk != 9 && qk !=1) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[3]*resu->sigmaij[9][3]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[2]*resu->sigmaij[9][3]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[0]*resu->sigmaij[9][3]*Hqqp_1;
	  }
	  if (qk != 3 && qk != 7) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[6]*resu->sigmaij[3][6]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[9]*resu->sigmaij[3][9]*Hqqp_1;
	  }
	  if (qk != 2 && qk!= 8) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[6]*resu->sigmaij[2][6]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[9]*resu->sigmaij[2][9]*Hqqp_1;
	  }
	  if (qk != 0 && qk!= 10) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[6]*resu->sigmaij[0][6]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[9]*resu->sigmaij[0][9]*Hqqp_1;
	  }
	  if (qk != 3 && qk !=7) {
	    QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[qk]*resu->sigmaij[6][3]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[qk]*resu->sigmaij[9][3]*Hqqp_2;
	  }
	  if( qk != 2 && qk != 8) {
	    QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[qk]*resu->sigmaij[6][2]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[qk]*resu->sigmaij[9][2]*Hqqp_2;
	  }
	  if (qk != 0 && qk!= 10) {
	    QQBN = QQBN + PR1.FP[6]*sudakq*PR2.FP[qk]*resu->sigmaij[6][0]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[9]*sudakq*PR2.FP[qk]*resu->sigmaij[9][0]*Hqqp_2;
	  }
	  if (qk != 6 && qk != 4) {
	    QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[qk]*resu->sigmaij[3][6]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[qk]*resu->sigmaij[2][6]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[qk]*resu->sigmaij[0][6]*Hqqp_2;
	  }
	  if (qk != 9 && qk != 1) {
	    QQBN = QQBN + PR1.FP[3]*sudakq*PR2.FP[qk]*resu->sigmaij[3][9]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[2]*sudakq*PR2.FP[qk]*resu->sigmaij[2][9]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[0]*sudakq*PR2.FP[qk]*resu->sigmaij[0][9]*Hqqp_2;
	  }
	}
      }
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[4][7]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[4][8]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[4][10]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[1][7]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[1][8]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[1][10]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[7][4]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[8][4]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[10][4]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[7][1]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[8][1]*Hgg;
      GGN = GGN + PR1.FP[5]*sudakq*PR2.FP[5]*resu->sigmaij[10][1]*Hgg;

      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[7]*resu->sigmaij[4][7]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[8]*resu->sigmaij[4][8]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[10]*resu->sigmaij[4][10]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[7]*resu->sigmaij[1][7]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[8]*resu->sigmaij[1][8]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[10]*resu->sigmaij[1][10]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[4]*resu->sigmaij[7][4]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[4]*resu->sigmaij[8][4]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[4]*resu->sigmaij[10][4]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[1]*resu->sigmaij[7][1]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[1]*resu->sigmaij[8][1]*Hqg_1;
      QGN_1 = QGN_1 + PR1.FP[5]*sudakq*PR2.FP[1]*resu->sigmaij[10][1]*Hqg_1;

      QGN_2 = QGN_2 + PR1.FP[4]*sudakq*PR2.FP[5]*resu->sigmaij[4][7]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[4]*sudakq*PR2.FP[5]*resu->sigmaij[4][8]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[4]*sudakq*PR2.FP[5]*resu->sigmaij[4][10]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[1]*sudakq*PR2.FP[5]*resu->sigmaij[1][7]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[1]*sudakq*PR2.FP[5]*resu->sigmaij[1][8]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[1]*sudakq*PR2.FP[5]*resu->sigmaij[1][10]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[7]*sudakq*PR2.FP[5]*resu->sigmaij[7][4]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[8]*sudakq*PR2.FP[5]*resu->sigmaij[8][4]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[10]*sudakq*PR2.FP[5]*resu->sigmaij[10][4]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[7]*sudakq*PR2.FP[5]*resu->sigmaij[7][1]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[8]*sudakq*PR2.FP[5]*resu->sigmaij[8][1]*Hqg_2;
      QGN_2 = QGN_2 + PR1.FP[10]*sudakq*PR2.FP[5]*resu->sigmaij[10][1]*Hqg_2;

      QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[7]*resu->sigmaij[4][7]*Hqqb[0];
      QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[8]*resu->sigmaij[4][8]*Hqqb[0];
      QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[10]*resu->sigmaij[4][10]*Hqqb[0];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[7]*resu->sigmaij[1][7]*Hqqb[3];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[8]*resu->sigmaij[1][8]*Hqqb[3];
      QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[10]*resu->sigmaij[1][10]*Hqqb[3];
      QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[4]*resu->sigmaij[7][4]*Hqqb[1];
      QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[4]*resu->sigmaij[8][4]*Hqqb[2];
      QQBN = QQBN + PR1.FP[10]*sudakq*PR2.FP[4]*resu->sigmaij[10][4]*Hqqb[4];
      QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[1]*resu->sigmaij[7][1]*Hqqb[1];
      QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[1]*resu->sigmaij[8][1]*Hqqb[2];
      QQBN = QQBN + PR1.FP[10]*sudakq*PR2.FP[1]*resu->sigmaij[10][1]*Hqqb[4];


      for (int qk = 0; qk < 2.0*Nf+1.0;qk++) {
	if (qk != 5) {
	  if (qk != 6 && qk != 4) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[7]*resu->sigmaij[4][7]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[8]*resu->sigmaij[4][8]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[10]*resu->sigmaij[4][10]*Hqqp_1;
	  }
	  if (qk != 9 && qk !=1) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[7]*resu->sigmaij[1][7]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[8]*resu->sigmaij[1][8]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[10]*resu->sigmaij[1][10]*Hqqp_1;
	  }
	  if (qk != 3 && qk != 7) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[4]*resu->sigmaij[7][4]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[1]*resu->sigmaij[7][1]*Hqqp_1;
	  }
	  if (qk != 2 && qk!= 8) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[4]*resu->sigmaij[8][4]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[4]*resu->sigmaij[8][1]*Hqqp_1;
	  }
	  if (qk != 0 && qk!= 10) {
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[4]*resu->sigmaij[10][4]*Hqqp_1;
	    QQBN = QQBN + PR1.FP[qk]*sudakq*PR2.FP[1]*resu->sigmaij[10][1]*Hqqp_1;
	  }
	  if (qk != 3 && qk !=7) {
	    QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[qk]*resu->sigmaij[4][7]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[qk]*resu->sigmaij[1][7]*Hqqp_2;
	  }
	  if( qk != 2 && qk != 8) {
	    QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[qk]*resu->sigmaij[4][8]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[qk]*resu->sigmaij[1][8]*Hqqp_2;
	  }
	  if (qk != 0 && qk!= 10) {
	    QQBN = QQBN + PR1.FP[4]*sudakq*PR2.FP[qk]*resu->sigmaij[4][10]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[1]*sudakq*PR2.FP[qk]*resu->sigmaij[1][10]*Hqqp_2;
	  }
	  if (qk != 6 && qk != 4) {
	    QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[qk]*resu->sigmaij[7][4]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[qk]*resu->sigmaij[8][4]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[10]*sudakq*PR2.FP[qk]*resu->sigmaij[10][4]*Hqqp_2;
	  }
	  if (qk != 9 && qk != 1) {
	    QQBN = QQBN + PR1.FP[7]*sudakq*PR2.FP[qk]*resu->sigmaij[7][1]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[8]*sudakq*PR2.FP[qk]*resu->sigmaij[8][1]*Hqqp_2;
	    QQBN = QQBN + PR1.FP[10]*sudakq*PR2.FP[qk]*resu->sigmaij[10][1]*Hqqp_2;
	  }
	}
      }
    }

    HCRNqq = GGN + QGN_1 + QGN_2 + QQBN;
//   std::cout << "HCRNqq = " << HCRNqq << std::endl;
//   std::cout << "GGN = " << GGN << " QGN_1 = " << QGN_1 << std::endl;
//   std::cout << "QGN_2 = " << QGN_2 << " QQBN = " << QQBN << std::endl;

//Gluon contribution - final state in Hst is always g g here, in development

//LL

    if (gg_order == 0){
      Hgg = 1.;
      Hqg_1 = zero_C;
      Hqg_2 = zero_C;
      Hqq = zero_C;
      Hqq_1 = zero_C;
      Hqq_2 = zero_C;
      Hqqp_1 = zero_C;
      Hqqp_2 = zero_C;
      for(int i = 0; i<5; i++) Hqqb[i] = zero_C;
    }

//NLL

    else if (gg_order == 1) {
      for(int i; i<Nf; i++) Hqqb[i] = zero_C;
//
      Hgq_1 = as2*N1.C1gq + as2*(-N1.gamma1gq)*(logmuf+2.*loga);
      Hgq_2 = as2*N2.C1gq + as2*(-N2.gamma1gq)*(logmuf+2.*loga);
//
      Hgg=1.+as2*(resu->H1g + N1.C1gg + N2.C1gg) - as2*(N1.gamma1gg + N2.gamma1gg)*(logmuf+2.*loga) + as2*(-4.*loga)*(indep->B1g+indep->A1g*loga) + as2*2.*pcF*indep->beta0*logmur;
    }

// Combine GG
    if(gg_order >= 0){
      std::complex<double> GGN = zero_C, GQN_1 = zero_C, GQN_2 = zero_C;

      GGN = GGN + Hgg*PR1.FP[Nf]*PR2.FP[Nf]*resu->sigmaij[Nf][Nf];

      for(int ii=-Nf; ii<=Nf; ii++){
        if (ii!=0){
          GQN_1 = GQN_1 + PR1.FP[Nf+ii]*PR2.FP[Nf]*resu->sigmaij[Nf][Nf]*Hgq_1;
          GQN_2 = GQN_2 + PR1.FP[Nf]*PR1.FP[Nf+ii]*resu->sigmaij[Nf][Nf]*Hgq_2;
        }
      }

      HCRNgg = sudakg*(GGN + GQN_1 + GQN_2);
    }

/*
// DEBUG
if((k_count0==407||k_count0==408) && k_count==2 && k_count2==10){
  end_ = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
}
*/

};



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//Sets the C1, C2, gamma (renormalised anomalous dimensions) and evolved (up to mures)
// pdf mellin moment values for use in calculating the Hst coefficients
void GetResuPars (int i, int isign, int ibeam, ResummationInfo* resuminfo, resu_PS* resu, std::complex<double> alpq, AllDep_N& ParamsforResu){

    std::complex<double> xn, qqi, qgf, gqi, ggi, ggf, ns1mi, ns1pi, ns1f,
      qq1f, qg1f, gq1i, gq1f, gg1i, gg1f;
    std::complex<double> uvi,dvi,usi,dsi,ssi,gli,chi,boi;

    int Nf = resuminfo->Nf;
    int ifit = resu->ifit;
    int ifit2 = resu->ifit2;


// RECOVER PARAMETERS from "resuminfo" structure

    if (isign == 1) {
      xn = resuminfo->contoursa.Np[i];
      qqi = resuminfo->PosBranch[i].qqi;
      qgf = resuminfo->PosBranch[i].qgf;
      gqi = resuminfo->PosBranch[i].gqi;
      ggi = resuminfo->PosBranch[i].ggi;
      ggf = resuminfo->PosBranch[i].ggf;
      ns1mi = resuminfo->PosBranch[i].ns1mi;
      ns1pi = resuminfo->PosBranch[i].ns1pi;
      ns1f = resuminfo->PosBranch[i].ns1f;
      qq1f = resuminfo->PosBranch[i].qq1f;
      qg1f = resuminfo->PosBranch[i].qg1f;
      gq1i = resuminfo->PosBranch[i].gq1i;
      gq1f = resuminfo->PosBranch[i].gq1f;
      gg1i = resuminfo->PosBranch[i].gg1i;
      gg1f = resuminfo->PosBranch[i].gg1f;

      if (ibeam == 1) {
        uvi = global_fitpars::pdfbeam1pos.uv[ifit2][ifit][i];
        dvi = global_fitpars::pdfbeam1pos.dv[ifit2][ifit][i];
        usi = global_fitpars::pdfbeam1pos.us[ifit2][ifit][i];
        dsi = global_fitpars::pdfbeam1pos.ds[ifit2][ifit][i];
        ssi = global_fitpars::pdfbeam1pos.ss[ifit2][ifit][i];
        gli = global_fitpars::pdfbeam1pos.gl[ifit2][ifit][i];
        chi = global_fitpars::pdfbeam1pos.ch[ifit2][ifit][i];
        boi = global_fitpars::pdfbeam1pos.bo[ifit2][ifit][i];
      }
      else if (ibeam == 2) {
        uvi = global_fitpars::pdfbeam2pos.uv[ifit2][ifit][i];
        dvi = global_fitpars::pdfbeam2pos.dv[ifit2][ifit][i];
        usi = global_fitpars::pdfbeam2pos.us[ifit2][ifit][i];
        dsi = global_fitpars::pdfbeam2pos.ds[ifit2][ifit][i];
        ssi = global_fitpars::pdfbeam2pos.ss[ifit2][ifit][i];
        gli = global_fitpars::pdfbeam2pos.gl[ifit2][ifit][i];
        chi = global_fitpars::pdfbeam2pos.ch[ifit2][ifit][i];
        boi = global_fitpars::pdfbeam2pos.bo[ifit2][ifit][i];
      }

    }
    else if (isign == -1) {
      xn = resuminfo->contoursa.Nm[i];
      qqi = resuminfo->NegBranch[i].qqi;
      qgf = resuminfo->NegBranch[i].qgf;
      gqi = resuminfo->NegBranch[i].gqi;
      ggi = resuminfo->NegBranch[i].ggi;
      ggf = resuminfo->NegBranch[i].ggf;
      ns1mi = resuminfo->NegBranch[i].ns1mi;
      ns1pi = resuminfo->NegBranch[i].ns1pi;
      ns1f = resuminfo->NegBranch[i].ns1f;
      qq1f = resuminfo->NegBranch[i].qq1f;
      qg1f = resuminfo->NegBranch[i].qg1f;
      gq1i = resuminfo->NegBranch[i].gq1i;
      gq1f = resuminfo->NegBranch[i].gq1f;
      gg1i = resuminfo->NegBranch[i].gg1i;
      gg1f = resuminfo->NegBranch[i].gg1f;
/*
      ParamsforResu.c1qg = resuminfo->NegBranch[i].C1qg;
      ParamsforResu.c1gq = resuminfo->NegBranch[i].C1gq;
      ParamsforResu.c1qq = resuminfo->NegBranch[i].C1qq;
      ParamsforResu.c1gg = resuminfo->NegBranch[i].C1gg;
      ParamsforResu.c2qgM = resuminfo->NegBranch[i].C2qg;
      ParamsforResu.c2NSqqM = resuminfo->NegBranch[i].C2NSqq;
      ParamsforResu.c2SqqbM = resuminfo->NegBranch[i].C2Sqqb;
      ParamsforResu.c2NSqqbM = resuminfo->NegBranch[i].C2NSqqb;
      ParamsforResu.gamma1qq = resuminfo->NegBranch[i].gamma1qq;
      ParamsforResu.gamma1qg = resuminfo->NegBranch[i].gamma1qg;
      ParamsforResu.gamma1gq = resuminfo->NegBranch[i].gamma1gq;
      ParamsforResu.gamma1gg = resuminfo->NegBranch[i].gamma1gg;
      ParamsforResu.gamma2qq = resuminfo->NegBranch[i].gamma2qq;
      ParamsforResu.gamma2qqV = resuminfo->NegBranch[i].gamma2qqV;
      ParamsforResu.gamma2qqbV = resuminfo->NegBranch[i].gamma2qqbV;
      ParamsforResu.gamma2qqS = resuminfo->NegBranch[i].gamma2qqS;
      ParamsforResu.gamma2qqbS = resuminfo->NegBranch[i].gamma2qqbS;
      ParamsforResu.gamma2qg = resuminfo->NegBranch[i].gamma2qg;
      ParamsforResu.gamma2gq = resuminfo->NegBranch[i].gamma2gq;
      ParamsforResu.gamma2gg = resuminfo->NegBranch[i].gamma2gg;
*/
      if (ibeam == 1) {
        uvi = global_fitpars::pdfbeam1min.uv[ifit2][ifit][i];
        dvi = global_fitpars::pdfbeam1min.dv[ifit2][ifit][i];
        usi = global_fitpars::pdfbeam1min.us[ifit2][ifit][i];
        dsi = global_fitpars::pdfbeam1min.ds[ifit2][ifit][i];
        ssi = global_fitpars::pdfbeam1min.ss[ifit2][ifit][i];
        gli = global_fitpars::pdfbeam1min.gl[ifit2][ifit][i];
        chi = global_fitpars::pdfbeam1min.ch[ifit2][ifit][i];
        boi = global_fitpars::pdfbeam1min.bo[ifit2][ifit][i];
      }
      else if (ibeam == 2) {
        uvi = global_fitpars::pdfbeam2min.uv[ifit2][ifit][i];
        dvi = global_fitpars::pdfbeam2min.dv[ifit2][ifit][i];
        usi = global_fitpars::pdfbeam2min.us[ifit2][ifit][i];
        dsi = global_fitpars::pdfbeam2min.ds[ifit2][ifit][i];
        ssi = global_fitpars::pdfbeam2min.ss[ifit2][ifit][i];
        gli = global_fitpars::pdfbeam2min.gl[ifit2][ifit][i];
        chi = global_fitpars::pdfbeam2min.ch[ifit2][ifit][i];
        boi = global_fitpars::pdfbeam2min.bo[ifit2][ifit][i];
      }

    }

/*
// Define renormalised anomalous dimensions from basic coefficients
    std::complex<double> Nfcomplex = std::complex<double>(Nf,0.0);
    std::complex<double> OneComplex = std::complex<double>(1.0,0.0);
    ParamsforResu.gamma1qq = -(qqi/4.0);
    ParamsforResu.gamma1qg = -(qgf/8.0);
    ParamsforResu.gamma1gq = -(gqi/4.0);
    ParamsforResu.gamma1gg = -((ggi + Nfcomplex*ggf)/4.0);
    ParamsforResu.gamma2qq = -(((ns1pi + Nfcomplex*ns1f) + Nfcomplex*qq1f)/8.0);
    ParamsforResu.gamma2qqV = -(ns1pi + 2.0*Nfcomplex*ns1f + ns1mi)/16.0;
    ParamsforResu.gamma2qqbV = -(ns1pi-ns1mi)/16.0;
    ParamsforResu.gamma2qqS = -(qq1f/16.0);
    ParamsforResu.gamma2qqbS = ParamsforResu.gamma2qqS;
    ParamsforResu.gamma2qg = -(qg1f/16.0);
    ParamsforResu.gamma2gq = -((gq1i + Nfcomplex*gq1f)/8.0);
    ParamsforResu.gamma2gg = -((gg1i + Nfcomplex*gg1f)/8.0);
*/


// Call PDF evolution routine RENO2
// As of current version, RENO2 is hard-coded to work in conjunction with the legacy PDF
// fitting routine. The output is the FP vector. Flavour code:
//  0   1   2   3   4   5   6   7   8
//  u   ub  d   db  s   ch  bo  gl  empty
// There is no top and it is assumed strange = antistrange, charm = anticharm, bottom = antibottom
//    std::cout << "resu->alphaqf = " << resu->alphaqf  << std::endl;
//    std::cout << "resu->a = " << resu->a  << std::endl;
//    std::cout << "resu->alphas = " << resu->alphas << std::endl;
//    std::cout << "resu->mur2 = " << resu->mur2 << std::endl;
//    std::cout << "resu->q2 = " << resu->q2 << std::endl;
//    std::cout << "resu->b0p = " << resu->b0p << std::endl;

    
    std::complex<double> FP[9];
    RENO2(xn, i, isign, ibeam, resuminfo->order, Nf, qqi, qgf, gqi, ggi, ggf, ns1mi, ns1pi, ns1f, qq1f, qg1f, gq1i, gq1f, gg1i, gg1f, uvi, dvi, usi, dsi, ssi, gli, chi, boi, resu->alphaqf, resu->q2, resu->b0p, alpq, resu->a, resu->alphas, resu->mur2, resuminfo->verbosity, FP);

     if (resuminfo->verbosity >= 50) {
      std::cout << "Post RENO2: " << std::endl; //RENO2 calculates the FP from all these inputs
      std::cout << "FP = " << FP[0] << " " << FP[1] << " " << FP[2] << " " << FP[3] << " " << FP[4] << " " << FP[5] << " " << FP[6] << " " << FP[7] << " " << FP[8] << std::endl;
     }


// Swap u <-> ubar and d <-> dbar for antiprotons in beam
    int ih;
    if(ibeam==1) ih = resuminfo->ih1;
    if(ibeam==2) ih = resuminfo->ih2;
    if (ih == -1) {
      std::complex<double> utemp, dtemp;
      utemp = FP[0];
      dtemp = FP[2];
      FP[0] = FP[1];
      FP[1] = utemp;
      FP[2] = FP[3];
      FP[3] = dtemp;
    }

// Convert to a more flexible flavour notation (Nf active flavours):
// FX[0] to FX[Nf-1]      : antiquarks
// FX[Nf]                 : gluons
// FX[Nf+1] to FX[2*Nf+1] : quarks

      ParamsforResu.FP.resize(2*Nf+1);

      // std::cout << "Nf = " << Nf << std::endl;

// It is implicitly assumed that Nf >= 2
      ParamsforResu.FP[Nf] = FP[7];   // gluon
      ParamsforResu.FP[Nf+1] = FP[0]; // up
      ParamsforResu.FP[Nf+2] = FP[2]; // down
      if (Nf >= 3) ParamsforResu.FP[Nf+3] = FP[4]; // strange
      else ParamsforResu.FP[Nf+3] = std::complex<double>(0., 0.);
      if (Nf >= 4) ParamsforResu.FP[Nf+4] = FP[5]; // charm
      else ParamsforResu.FP[Nf+4] = std::complex<double>(0., 0.);
      if (Nf >= 5) ParamsforResu.FP[Nf+5] = FP[6]; // bottom
      else ParamsforResu.FP[Nf+5] = std::complex<double>(0., 0.);

      ParamsforResu.FP[Nf-1] = FP[1]; // antiup
      ParamsforResu.FP[Nf-2] = FP[3]; // antidown
      if (Nf >= 3) ParamsforResu.FP[Nf-3] = FP[4]; // antistrange
      else ParamsforResu.FP[Nf-3] = std::complex<double>(0., 0.);
      if (Nf >= 4) ParamsforResu.FP[Nf-4] = FP[5]; // anticharm
      else ParamsforResu.FP[Nf-4] = std::complex<double>(0., 0.);
      if (Nf >= 5) ParamsforResu.FP[Nf-5] = FP[6]; // antibottom
      else ParamsforResu.FP[Nf-5] = std::complex<double>(0., 0.);

// Do to current restrictions on RENO2, this can't work for Nf>5.

}






/* DEBUG */

void AllDep_N::preent(){

/*
  std::cout << "Local Mellin C1 pars\n";
  std::cout << c1qg << std::endl;
  std::cout << c1gq << std::endl;
  std::cout << c1qq << std::endl;
  std::cout << c1gg << std::endl;
  std::cout << "Local Mellin C2 pars\n";
  std::cout << c2qgM << std::endl;
  std::cout << c2NSqqM << std::endl;
  std::cout << c2SqqbM << std::endl;
  std::cout << c2NSqqbM << std::endl;

  std::cout << "Local Mellin Anomalous Dimensions 1\n";
  std::cout << gamma1qq << std::endl;
  std::cout << gamma1qg << std::endl;
  std::cout << gamma1gq << std::endl;
  std::cout << gamma1gg << std::endl;
  std::cout << "Local Mellin Anomalous Dimensions 2\n";
  std::cout << gamma2qq << std::endl;
  std::cout << gamma2qqV << std::endl;
  std::cout << gamma2qqbV << std::endl;
  std::cout << gamma2qqS << std::endl;
  std::cout << gamma2qqbS << std::endl;
  std::cout << gamma2qg << std::endl;
  std::cout << gamma2gq << std::endl;
  std::cout << gamma2gg << std::endl;
*/

  std::cout << "Local Mellin PDFs, evolved 1\n";
  for(int ii=0; ii<FP.size(); ii++) std::cout << FP[ii] << std::endl;

  std::cout << std::endl;

}
