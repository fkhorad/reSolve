//Takes the resummed form factor, multiplies by the sudakov factor and does the inverse mellin transform. The double mellin integral and all its calls (e.g. Hardfns_calc) are the time critical parts of the code.

#include "inv_mellin.h"

#include <ctime>
#include <limits>
#include <array>

#include "constants.h"

#include "resummation_input.h"
#include "sudakov.h"
#include "hardfns.h"
#include "resu_PS.h"
#include "evolution.h"


//DEBUG
#include <chrono>
extern int k_count, k_count2, k_count0, k_count00;
int k_count2;

extern std::chrono::high_resolution_clock::time_point start, end_;
std::chrono::high_resolution_clock::time_point start, end_;


double inversemellin_resummed(std::complex<double> b, resu_PS* resu, ResummationInfo* resuminfo) {


//DEBUG
    k_count = k_count+1;
    int ktest = 1, i1test = 8, i2test = 3;

    double invres;

// N1, N2-INDEPENDENT PART

// Basic variables

    double pi = k_constants::pi;
    const int m_points = k_constants::mellin_points;
    std::complex<double> OneComplex = std::complex<double>(1.0,0.0);


//Debug
    start = std::chrono::high_resolution_clock::now();


    double q2 = resu->q2;
    double qt2 = resu->qt2;
    double eta = resu->eta;
    double mur2 = resu->mur2;
    double a = resu->a;
    double b0p = resu->b0p;
    double x = resu->x;
    double alphas = resu->alphas;
    double alphaqf = resu->alphaqf;

    int ih1 = 0, ih2 = 0;
    ih1 = resuminfo->ih1;
    ih2 = resuminfo->ih2;
    int order = resuminfo->order;
    int gg_order = resuminfo->gg_order;
    int qq_order = resuminfo->qq_order;
    double ggnp = resuminfo->ggnp;
    double gqnp = resuminfo->gqnp;
    int ifit = resu->ifit;
    int ifit2 = resu->ifit2;
    double Ax = std::log(x);
    resummationIndepfns* resuNindep = &resuminfo->resuNindep;

    if (resuminfo->verbosity >= 16) {
      std::cout << "inv_mellin.cc" << std::endl;
    }

// Regulated b
// To resolve singularities - I want to pass bstar to Hardfns_calc to resolve singularities seen at large b
    std::complex<double> bstar, blog, blim;
    blim = b0p*1.0/std::pow(q2,0.5)*std::exp(1.0/(2.0*alphas*resuNindep->beta0));
    bstar=b/std::pow(1.0+(b*b)/(blim*blim),0.5); //resolves singularities at large b
    blog = std::log(q2*bstar*bstar/(b0p*b0p)+1.0); //modified sudakov resolves singularities at small b

//Use bstar instead to avoid issues of evolution over large scale differences in alphasl -> FX1/2 matrices -> Hardfns_calc, bstar is same as b for small and medium b but cutsoff at blim for large b (essentially a small qT cutoff)
    std::complex<double> scale2 = std::complex<double>(b0p*b0p/(bstar.real()*bstar.real()),0.0);

// SUDAKOV
    std::complex<double> sudakq, sudakg;
    sudakq = Sc(b, q2, b0p, alphas, mur2, a, resuNindep, order, gqnp, ggnp, 'q',
                resuminfo->verbosity);
    sudakg = Sc(b, q2, b0p, alphas, mur2, a, resuNindep, order, gqnp, ggnp, 'g',
                resuminfo->verbosity);

// Use Bessel function as integrand is only dependent on b magnitude, not direction
// Must use actual b here as is part of measure for inverse fourier transform
    double xj0 = 0.5*j0(std::pow(qt2,0.5)*b.real()); //bessel function used, must use actual b here as is part of measure for inverse fourier transform

    double factorfin = b.real()*xj0;
// Get alphasl, aexp, aexpb
    std::complex<double> alphasl = std::complex<double>(1.0,0.0);
    std::complex<double> aexp, aexpb;
    aexp = std::complex<double>(0.0,0.0);
    aexpb = std::complex<double>(0.0,0.0);
// (Can take real part of scale2 as it is real anyway, just defined as complex)
    alphaslcalc(resu->q2, resu->b0p, resu->a, resu->alphas, resu->mur2,
                scale2.real(), resuNindep, alphasl, aexp, aexpb, resuminfo->verbosity, resuminfo->order);
    std::complex<double> alpq = resu->alphaqf*alphasl;


// BEGIN SETUP of inverse N1, N2, Mellin integrals

    double Ax1 = (Ax+2*eta)/2.0;
    double Ax2 = (Ax-2*eta)/2.0;
    double xx1 = std::exp(Ax1);
    double xx2 = std::exp(Ax2);
    double xnormal = 1.0;


// Setup # of points in the contours and extra normalization (xnormal):
// to avoid issues at x>0.87d0 (very large rapidities) suppress the result
    int nmax1, nmax2;
    if (xx1 > 0.87) {xnormal=(1-xx1)*(1-xx1)*(1-xx1)/((1-0.87)*(1-0.87)*(1-0.87));}
    if (xx2 > 0.87) {xnormal=(1-xx2)*(1-xx2)*(1-xx2)/((1-0.87)*(1-0.87)*(1-0.87));}
    if (xx1 < 0.001) { nmax1 = 40;}
    else if (xx1 < 0.05) { nmax1 = 40;}
    else if (xx1 < 0.2) { nmax1 = 56;}
    else if (xx1 < 0.4) { nmax1 = 72;}
    else if (xx1 < 0.7) { nmax1 = 88;}
    else { nmax1 = 88;}

    if (xx2 < 0.001) { nmax2 = 40;}
    else if (xx2 < 0.05) { nmax2 = 40;}
    else if (xx2 < 0.2) { nmax2 = 56;}
    else if (xx2 < 0.4) { nmax2 = 72;}
    else if (xx2 < 0.7) { nmax2 = 88;}
    else { nmax2 = 88;}

// END SETUP, start calculating quantities



// Split up double mellin integration into N1 dependent parts, N2 dependent aparts and N1 and N2 dependent parts to save time

    contour_info* contours = &resuminfo->contoursa;

    AllDep_N ParamsforResu1, ParamsforResu2p, ParamsforResu2n;
    std::array<AllDep_N, m_points> P1, P2p, P2n;
    std::complex<double> hcrnqq, hcrngg;

//    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    end_ = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//    std::cout<<"invres setup chkpoint: "<< time_span.count() <<'\n';


//N1 dependent only
    std::complex<double> xn1pos, xnm1pos, xn1neg, xnm1neg;
    std::complex<double> cex1pos[m_points];
    std::complex<double> cex1neg[m_points];
    for (int i1 = 0; i1 < nmax1; i1++) {
// Contour: Positive branch
      xn1pos = contours->Np[i1];
      xnm1pos = -xn1pos*Ax1;
      cex1pos[i1] = std::exp(xnm1pos)/pi*contours->phicomplexp;
      xn1pos = xn1pos + OneComplex;
//Negative branch
      xn1neg = contours->Np[i1]; //Np not Nm as do over both parts here
      xnm1neg = -xn1neg*Ax1;
      cex1neg[i1] = std::exp(xnm1neg)/pi*contours->phicomplexp; //phicomplexp not phicomplexm as do over both parts here
      xn1neg = xn1neg + OneComplex;

      GetResuPars(i1, 1, 1, resuminfo, resu, alpq, ParamsforResu1);
      // std::cout << "alpq = " << alpq << std::endl;
      P1[i1] = ParamsforResu1;
    }

    end_ = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//    std::cout<<"invres N1 chkpoint: "<< time_span.count() <<'\n';


//N2 dependent only
    std::complex<double> cex2pos[m_points];
    std::complex<double> cex2neg[m_points];
    std::complex<double> xn2pos, xnm2pos, xn2neg, xnm2neg;
    for (int i2 = 0; i2 < nmax2; i2++) {
//Positive branch
      xn2pos = contours->Np[i2];
      xnm2pos = -xn2pos*Ax2;
      cex2pos[i2] = std::exp(xnm2pos)/pi*contours->phicomplexp;
//Negative branch
      xn2neg = contours->Nm[i2];
      xnm2neg = -xn2neg*Ax2;
      cex2neg[i2] = std::exp(xnm2neg)/pi*contours->phicomplexm;

//Need to call GetResuPars 2x for beam 2, once for isign + 1 and once for isign -1.
      GetResuPars(i2, 1, 2, resuminfo, resu, alpq, ParamsforResu2p);
      GetResuPars(i2, -1, 2, resuminfo, resu, alpq, ParamsforResu2n);

      P2p[i2] = ParamsforResu2p;
      P2n[i2] = ParamsforResu2n;
    }

    end_ = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//    std::cout<<"invres N2 chkpoint: "<< time_span.count() <<'\n';


// Now the main integration loop

    double fz = 0., fun = 0.;
    k_count2=0;
    for (int i1 = 0; i1 < nmax1; i1++) {
      for (int i2 = 0; i2 < nmax2; i2++) {
        k_count2++;

// Positive branch contribution

// Note: sudakovs passed so they combine with hardfactors at this stage, the hcrnqq and hcrngg out then are actually hcrnqq*sudakq and hcrngg*sudakg to avoid issues when hard factors become very large at same time as sudakovs become very small or vice versa

if((k_count0==407||k_count0==408) && k_count==2 && k_count2==10){
    end_ = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//    std::cout<<"invres pre-hardfns chkpoint: "<< time_span.count() <<'\n';
}

 // std::cout << "resuminfo->pcF = " << resuminfo->pcF << std::endl;
 // std::cout << "resuminfo->tot_em_charge = " << resuminfo->tot_em_charge << std::endl;

Hardfns_calc (qq_order, gg_order, resuminfo->pcF, resuminfo->tot_em_charge, resu, resuNindep, P1[i1], P2p[i2], resuminfo->PosBranch[i1], resuminfo->PosBranch[i2], alpq, aexp, aexpb, sudakq, sudakg, hcrnqq, hcrngg);

// std::cout << "hcrnqq/sudakq = " << hcrnqq/sudakq << std::endl;

if((k_count0==407||k_count0==408) && k_count==2 && k_count2==10){
    end_ = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//    std::cout<<"invres mid-hardfns chkpoint: "<< time_span.count() <<'\n';
}


// The hcrnqq and hcrngg out then are actually hcrnqq*sudakq and hcrngg*sudakg to avoid issues when hard factors become very large at same time as sudakovs become very small or vice versa
        std::complex<double> int1 = (hcrnqq + hcrngg)
          *cex1pos[i1]*cex2pos[i2]*contours->weights[i1]*contours->weights[i2]*factorfin;
	// std::cout << "hcrnqq = " << hcrnqq << " hcrngg = " << hcrngg << std::endl;
	// std::cout << "cex1pos[i1] = " << cex1pos[i1] << " cex2pos[i2] = " << cex2pos[i2]  << std::endl;
	// std::cout << "contours->weights[i1] = " << contours->weights[i1] << " contours->weights[i2] = " << contours->weights[i2]  << std::endl;
	


// Negative branch contribution

// Reset main functions to zero before negative branch
        hcrnqq = std::complex<double>(0.0,0.0);
        hcrngg = std::complex<double>(0.0,0.0);


// Note: sudakovs passed so they combine with hardfactors at this stage, the hcrnqq and hcrngg out then are actually hcrnqq*sudakq and hcrngg*sudakg to avoid issues when hard factors become very large at same time as sudakovs become very small or vice versa
       Hardfns_calc (qq_order, gg_order, resuminfo->pcF, resuminfo->tot_em_charge, resu, resuNindep, P1[i1], P2n[i2], resuminfo->PosBranch[i1], resuminfo->NegBranch[i2], alpq, aexp, aexpb, sudakq, sudakg, hcrnqq, hcrngg);


//The hcrnqq and hcrngg out then are actually hcrnqq*sudakq and hcrngg*sudakg to avoid issues when hard factors become very large at same time as sudakovs become very small or vice versa
       std::complex<double> int2 = (hcrnqq + hcrngg)
         *cex1neg[i1]*cex2neg[i2]*contours->weights[i1]*contours->weights[i2]*factorfin;

//old code used dble in fortran to take real part here I think.....
       fz = -0.5*(int1.real()-int2.real());
       // std::cout << "int1.real() = " << int1.real() << std::endl;
       // std::cout << "int2.real() = " << int2.real() << std::endl;


       fun = fun + fz;
       // std::cout << "fz = " << fz << std::endl;
       // std::cout << "fun = " << fun << std::endl;


       if(i1==nmax1/4 && i2==nmax2/4){
         end_ = std::chrono::high_resolution_clock::now();
         time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//         std::cout<<"quarter int chkpoint: "<< time_span.count() <<'\n';
       }
       if(i1==nmax1/2 && i2==nmax2/2){
         end_ = std::chrono::high_resolution_clock::now();
         time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//         std::cout<<"half int chkpoint: "<< time_span.count() <<'\n';
       }
       if(i1==3*nmax1/4 && i2==3*nmax2/4){
         end_ = std::chrono::high_resolution_clock::now();
         time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);
//         std::cout<<"3/4 int chkpoint: "<< time_span.count() <<'\n';
       }


      }
    }


    invres = fun*xnormal;

    // std::cout << "invres = " << invres << std::endl;
    // std::cout << "fun = " << fun << std::endl;
    // std::cout << "xnormal = " << xnormal << std::endl;


    end_ = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end_ - start);

//     if(time_span.count()>0.2){
// //      std::cout << "WARNING!  " << time_span.count() << std::endl;
//       invres = 0.;
//     }


    return invres;

}
