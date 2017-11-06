
#include "diphoton_ps.h"

#include "phase_space.h"
#include "constants.h"

#include "diphoton_input.h"

#include "resu_preproc.h"
#include "hardfns.h"


double diphoton_ps(const double x[], double CM_energy,
                 ResummationInfo* resuminfo, double& randsjacob,
		   PSpoint& pp1, PSdep_variables& resuvars){
    // NOTE PROCESS ONLY PASSED HERE TO ALLOW PROCESS = -1 TO MEAN BORN ONLY FOR TESTING
    double pi = k_constants::pi;

//Assign q2, qt2, eta, thetaCM, phiCM with randoms from vegas
    double q2 = std::pow(resuminfo->QQ_Min,2) + x[0]*(std::pow(resuminfo->QQ_Max,2) - std::pow(resuminfo->QQ_Min,2));
    double qt2 = std::pow(resuminfo->QT_Min,2) + x[1]*(std::pow(resuminfo->QT_Max,2) - std::pow(resuminfo->QT_Min,2));
    double eta = resuminfo->eta_Min + x[2]*(resuminfo->eta_Max - resuminfo->eta_Min);
    double thetaCM = pi/2*x[3]; //thetaCM between 0 and pi/2
                                //(takes into account the 1/2 factor for
                                // identical photons)
    double phiCM = 2*pi*x[4];   //phiCM between 0 and 2pi
    randsjacob = (std::pow(resuminfo->QQ_Max,2)-std::pow(resuminfo->QQ_Min,2))*(std::pow(resuminfo->QT_Max,2)-std::pow(resuminfo->QT_Min,2))*(resuminfo->eta_Max-resuminfo->eta_Min)*pi/2.0*2.0*pi;
    if (resuminfo->process == -1) {
	qt2 = 0.0; //for Born only option for testing
	randsjacob = (std::pow(resuminfo->QQ_Max,2)-std::pow(resuminfo->QQ_Min,2))*(resuminfo->eta_Max-resuminfo->eta_Min)*pi/2.0*2.0*pi; //no qT integration in Born only case as qT = 0 at LO
    }

    if (resuminfo->verbosity >= 10) {
        std::cout << "x = " << x[0] << " " << x[1] << " " << x[2]
        << " " << x[3] << " " << x[4] << std::endl;
        std::cout << "q2 = " << q2 << std::endl;
        std::cout << "qt2 = " << qt2 << std::endl;
        std::cout << "eta = " << eta << std::endl;
        std::cout << "thetaCM = " << thetaCM << std::endl;
        std::cout << "phiCM = " << phiCM << std::endl;
    }

    four_momentum k1(4), k2(4), qq1(4), qq2(4);
    double costheta = 0., jacob = 0.;

    costheta = kinematics(q2, qt2, eta, thetaCM, phiCM, CM_energy,
                          k1, k2, qq1, qq2, resuminfo->verbosity); //calculates costheta and qq1 and qq2 in lab frame (i.e. boosted back to lab frame from CM fram
    jacob = std::pow(1-std::pow(costheta,2),0.5)/(2*pi); //phase space jacobian

    if (resuminfo->verbosity >= 12) {
        LorPrint(k1);
        LorPrint(k2);
        LorPrint(qq1);
        LorPrint(qq2);
    }

    PSpoint pp1fin;
    double x1 = x[0]; double x2 = x[1];
    double en_hat = std::sqrt(x1*x2)*CM_energy;
    const double energy[3] = {en_hat, 0., 0.};
    set_PS_twobody(x+2, energy, pp1fin); //x+2 here is memory address of randoms[2] onwards, i.e randoms[2],[3] and [4] used not [0] and [1]

    pp1.set_dim(4);
    //k1, k2, qq1, qq2 are all lab frame momenta
    pp1.set_mom(0, k1); //Quark1 momentum
    pp1.set_mom(1, k2); //Quark2 momentum
    pp1.set_mom(2, qq1); //Diphoton1 momentum, Save as 2 as 0 and 1 are quark momenta
    pp1.set_mom(3, qq2); //Diphoton2 momentum, Save as 3 as 0 and 1 are quark momenta
    pp1.set_products();

    if (resuminfo->verbosity >= 12) {
        LorPrint(pp1.mom(0));
        LorPrint(pp1.mom(1));
        LorPrint(pp1.mom(2));
        LorPrint(pp1.mom(3));
    }

// PSdep_Variables calculation -- not really ideal, but ok

    double mur = resuminfo->mu_R;
    double muf = resuminfo->mu_F;
    double mures = resuminfo->mu_S;

    // std::cout << "mur = " << mur << std::endl;
    // std::cout << "muf = " << muf << std::endl;
    // std::cout << "mures = " << mures << std::endl;
    // std::cout << "muR_flag = " << resuminfo->muR_flag << std::endl;
    // std::cout << "muF_flag = " << resuminfo->muF_flag << std::endl;

//If muR_flag == 1 set renormalisation and resummation scales to q2 multiplied by fraction set in usual direct input (times 1/2 for resummation)
    if (resuminfo->muR_flag == 1) {
      mur = std::pow(q2,0.5)*resuminfo->mu_R;
      mures = std::pow(q2,0.5)*resuminfo->mu_S/2;
    }
//If muF_flag == 1 set factorisation scale to q2 multiplied by fraction set in usual direct input
    if (resuminfo->muF_flag == 1) {
      muf = std::pow(q2,0.5)*resuminfo->mu_F;
    }

    // std::cout << "mur = " << mur << std::endl;
    // std::cout << "muf = " << muf << std::endl;
    // std::cout << "mures = " << mures << std::endl;

    //Set to values now flags accounted for to pass on through resuminfo
    resuvars.mur2 = pow(mur,2);
    resuvars.mures2 = pow(mures,2);
    resuvars.muf2 = pow(muf,2);


    double sqs = CM_energy;
    double alpqf = alphas_(&mures)/(4*pi);
    double as = alphas_(&mur)/pi;
    double xx = q2/std::pow(sqs,2);
    double mur2 = std::pow(mur,2);
    double muf2 = std::pow(muf,2);
    double a = std::pow(q2,0.5)/mures;
    double mures2 = std::pow(mures,2);
    double Euler = 0.57721566;
    double b0p = 2*std::exp(-Euler)*a;

    resuvars.x = xx;
    resuvars.q2 = q2;
    resuvars.qt2 = qt2;
    resuvars.eta = eta;
    resuvars.mur2 = mur2;
    resuvars.muf2 = muf2;
    resuvars.mures2 = mures2;
    resuvars.a = a;
    resuvars.b0p = b0p;
    resuvars.alphas = as;
    resuvars.alphaqf = alpqf;

    return jacob;

}


double kinematics(double M2, double qT2, double eta, double thetaCM,
                  double phiCM, double sqs, four_momentum& k1,
                  four_momentum& k2, four_momentum& qq1, four_momentum& qq2,
                  int verbosity) {
//Full exclusive final state kinematics. Due to an intrinsic kinematic ambiguity of the resummation formula, we need to specify an additional (abitrary, unphysical) qT1 in order to completely define the initial state kinematics with momentum conservation. We work in the hadron CM frame.
    double MT2 = 0.0, MT = 0.0, MM = 0.0, qT = 0.0, qP1 = 0.0, k1TqT = 0.0, x1 = 0.0, y1 = 0.0, theta1CM = 0.0, phi1CM = 0.0, k1T_2 = 0.0, t =0.0, costheta = 0.0 , kin = 0.0;
    four_momentum qVec(4), P1(4), P2(4), k1T(4), q1CM(4), q2CM(4);
    for (int i = 0; i<4; i++) {
        qVec[i] = 0.0;
        P1[i] = 0.0;
        P2[i] = 0.0;
        k1T[i] = 0.0;
        k1[i] = 0.0;
        k2[i] = 0.0;
        q1CM[i] = 0.0;
        q2CM[i] = 0.0;
    }

    MT2 = M2 + qT2;
    MT = std::pow(MT2,0.5);
    MM = std::pow(M2,0.5);
    qT = std::pow(qT2, 0.5);
    qVec[0] = MT*cosh(eta);
    qVec[1] = qT;
    qVec[2] = 0.0;
    qVec[3] = MT*sinh(eta);

    if (verbosity >= 13) {
        std::cout << "MT2 = " << MT2 << std::endl;
        std::cout << "MT = " << MT << std::endl;
        std::cout << "MM = " << MM << std::endl;
        std::cout << "qT = " << qT << std::endl;
        std::cout << "qVec = " << qVec[0] << " " << qVec[1] << " " << qVec[2] << " " << qVec[3] << std::endl;
    }
    P1[0] = sqs/2.0;
    P1[1] = 0.0;
    P1[2] = 0.0;
    P1[3] = sqs/2.0;
    P2[0] = sqs/2.0;
    P2[1] = 0.0;
    P2[2] = 0.0;
    P2[3] = -sqs/2.0;
    qP1 = LorDot(qVec, P1);
    if (verbosity >= 13) {
        std::cout << "qVec = " << qVec[0] << " " <<qVec[1] << " " <<qVec[2] << " " <<qVec[3] << std::endl;
        std::cout << "P1 = " << P1[0] << " " <<P1[1] << " " <<P1[2] << " " <<P1[3] << std::endl;
        std::cout << "P2 = " << P2[0] << " " <<P2[1] << " " <<P2[2] << " " <<P2[3] << std::endl;
        std::cout << "qP1 = " << qP1 << std::endl;
    }
    //Simplest prescription for qT1

    k1T[0] = 0.0 ;
    k1T[1] = qT/2.0;
    k1T[2] = 0.0;
    k1T[3] = 0.0;// By construction, the x and y components of k1 are k1T
    k1T_2 = -LorDot(k1T,k1T);
    k1TqT = LorDot(k1T,qVec);

    if (verbosity >= 13) {
        std::cout << "k1T = " << k1T[0] << " " << k1T[1] << " " << k1T[2] << " " << k1T[3] << std::endl;
        std::cout << "k1T_2 =  " << k1T_2 << std::endl;
        std::cout << "k1TqT = " << k1TqT << std::endl;
    }

    //Initial state:
    x1 = (M2 - 2.0*k1TqT + pow((M2-2.0*k1TqT)*(M2-2.0*k1TqT)-4.0*MT2*k1T_2,0.5))/(4.0*qP1); //P1 coefficient of k1
    y1 = k1T_2/(x1*sqs*sqs); //P2 coefficient of k1

    //As these are four_momentum objects you have to tell it how to do the calculations component wise.
    for(int ii=0; ii<4; ii++) {
        k1[ii] = x1*P1[ii] + y1*P2[ii] + k1T[ii];
        k2[ii] = qVec[ii] - k1[ii];
    }

    if (verbosity >= 13) {
        std::cout << "x1 = " << x1 << std::endl;
        std::cout << "y1 = " << y1 << std::endl;
        std::cout << "k1 = " << k1[0] << " " << k1[1] << " " << k1[2] << " " << k1[3] << std::endl;
        std::cout << "k2 = " << k2[0] << " " << k2[1] << " " << k2[2] << " " << k2[3] << std::endl;
    }

    //Final state
    q1CM[0] = MM/2.0;
    q1CM[1] = MM/2.0*sin(thetaCM)*cos(phiCM);
    q1CM[2] = MM/2.0*sin(thetaCM)*sin(phiCM);
    q1CM[3] = MM/2.0*cos(thetaCM);
    q2CM[0] = MM/2.0;
    q2CM[1] = -MM/2.0*sin(thetaCM)*cos(phiCM);
    q2CM[2] = -MM/2.0*sin(thetaCM)*sin(phiCM);
    q2CM[3] = -MM/2.0*cos(thetaCM);

    if (verbosity >= 13) {
        std::cout << "thetaCM = " << thetaCM << std::endl;
        std::cout << "phiCM = " << phiCM << std::endl;
        std::cout << "q1CM = " << q1CM[0] << " " << q1CM[1] << " " << q1CM[2] << " " << q1CM[3] << std::endl;
        std::cout << "q2CM = " << q2CM[0] << " " << q2CM[1] << " " << q2CM[2] << " " << q2CM[3] << std::endl;
    }

    std::vector<four_momentum> lor1;
    
    lor1 = SetBoost(qVec, false);

    if (verbosity >= 13) {
        std::cout << "lor1 = " << std::endl;
        std::cout << lor1[0][0] << " " << lor1[0][1] << " " << lor1[0][2] << " " << lor1[0][3] << std::endl;
        std::cout << lor1[1][0] << " " << lor1[1][1] << " " << lor1[1][2] << " " << lor1[1][3] << std::endl;
        std::cout << lor1[2][0] << " " << lor1[2][1] << " " << lor1[2][2] << " " << lor1[2][3] << std::endl;
        std::cout << lor1[3][0] << " " << lor1[3][1] << " " << lor1[3][2] << " " << lor1[3][3] << std::endl;
    }
    qq1 = LorDot(lor1,q1CM);
    qq2 = LorDot(lor1,q2CM);
//    MatrixMult4xvec(q1CM,lor1,qq1);
//    MatrixMult4xvec(q2CM,lor1,qq2);
    t = -2*LorDot(k1,qq1);
    costheta = 1 + 2*t/M2;

    return costheta;
}
