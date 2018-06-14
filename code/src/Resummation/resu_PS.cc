#include "resu_PS.h"

#include "resummation_input.h"
#include "constants.h"


void in_state_w_recoil(four_momentum qvec, double sqs, double alpha, four_momentum& k1, four_momentum& k2){
// To enforce momentum conservation with qT != 0 and unresolved radiation, we make the two initial state partons recoil. There's an intrinsic amiguity in the way the qT is distributed between parton 1 and parton 2.

    double M2, MT2, qT2;

    k1.resize(4);
    k2.resize(4);

    qT2 = qvec[1]*qvec[1] + qvec[2]*qvec[2];
    M2 = LorDot(qvec,qvec);
    MT2 = M2 + qT2;

// kT1 PRESCRIPTION (alpha = 0 --> CS frame)
    double kT1 = 0., kT2 = 0.;
    kT1 = qvec[1]/2.*(1+alpha);
    kT2 = qvec[2]/2.*(1+alpha);
    double qTkT, kTsq;
    qTkT = qvec[1]*kT1+qvec[2]*kT2;
    kTsq = std::pow(kT1,2)+std::pow(kT2,2);

    double zeta1, qP1, qP2, P1P2;
    zeta1 = 1./(2.*M2)*(M2+2.*(qTkT) + std::sqrt(std::pow(M2+2.0*qTkT,2)-4.*MT2*kTsq));
    qP1 = (qvec[0]-qvec[3])*sqs/2.;
    qP2 = (qvec[0]+qvec[3])*sqs/2.;

    // std::cout << "qP1 = " << qP1 << std::endl;
    
    P1P2 = sqs*sqs/2.;

    k1[0] = sqs/2.*(zeta1*M2/(2.*qP1) + kTsq/zeta1*qP1/(M2*P1P2));
    k1[1] = kT1;
    k1[2] = kT2;
    k1[3] = sqs/2.*(zeta1*M2/(2.*qP1) - kTsq/zeta1*qP1/(M2*P1P2));

    // std::cout << "k1 = " << k1[0] << " " << k1[1] << " " << k1[2] << " " << k1[3] << std::endl;

    for(int i=0; i<4; i++) k2[i] = qvec[i]-k1[i];

}


void resu_PS::set(const ResummationInfo& res_1, double q2_in, double eta_in, double qT2_in){

    double pi = k_constants::pi;

    double mur = res_1.mu_R;
    double muf = res_1.mu_F;
    double mures = res_1.mu_S;
    double qq = std::sqrt(q2_in);

    q2 = q2_in;
    qt2 = qT2_in;
    eta = eta_in;

//If muR_flag == 1 set renormalisation and resummation scales to q2 multiplied by fraction set in usual direct input (times 1/2 for resummation)
    if (res_1.muR_flag == 1) {
      mur = qq*res_1.mu_R;
      mures = qq*res_1.mu_S/2.;
    }
//If muF_flag == 1 set factorisation scale to q2 multiplied by fraction set in usual direct input
    if (res_1.muF_flag == 1) {
      muf = qq*res_1.mu_F;
    }

//Set to values now flags accounted for to pass on through resuminfo
    mur2 = std::pow(mur,2);
    mures2 = std::pow(mures,2);
    muf2 = std::pow(muf,2);

    // alphas = alphas_(&mures)/pi;
    // alphaqf = alphas_(&mur)/(4*pi);
    //TEMP switch mur and mures as I believe we have it opposite here to how it was previously and could be causing discrepancy
    alphas = alphas_(&mur)/pi;
    alphaqf = alphas_(&mures)/(4*pi);
    double sqs = res_1.CM_energy;
    x = q2_in/std::pow(sqs,2);
    a = qq/mures;
    double Euler = 0.57721566;
    b0p = 2*std::exp(-Euler)*a;

    int Nf = res_1.Nf;
    H2q.resize(Nf);
    sigmaij.resize(2*Nf+1);
    for(int ii=0; ii<2*Nf+1; ii++){
      sigmaij[ii].resize(2*Nf+1);
    }

}
