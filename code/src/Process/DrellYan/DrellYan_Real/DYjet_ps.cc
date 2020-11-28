#include "DYjet_ps.h"

#include <fstream> // TEMP
#include <sstream> // TEMP
#include <iomanip>

#include "phase_space.h"
#include "constants.h"

#include "drellyan_input.h"
#include "DYjet_hard.h"

#include "resu_PS.h"


double DYjet_ps(const double x[], drellyan_input* drellyan_in, drellyan_amplitude& amp1, resu_PS& resuPS_0, double &ss_hat, double &etaa_hat){

    double pi = k_constants::pi;

// Dimension check: this routine needs a total of 7 randoms to generate the phase space -- stop if dimension is incorrect
    if(!(drellyan_in->ndim == 7) ){
      std::cout << "Something's wrong: ndim does not match phase-space generation routine" << std::endl;
      exit(EXIT_FAILURE);
    }

// Assign q2, qt2, eta, thetaCM, phiCM with randoms from vegas
    double q2, qt2, eta, thetaCM, c_thetaCM, phiCM;
    double randsjacob = 1.;
    double CM_energy = drellyan_in->res_1.CM_energy;

    //Trying a different order of integration variables for the randoms
    //First set s_hat since it is involved in many of the limits, sets x = s_hat/s
    double QQ_Min = std::max(0., drellyan_in->res_1.QQ_Min);
    double s_hat = QQ_Min*QQ_Min + x[0]*(CM_energy*CM_energy-QQ_Min*QQ_Min);
    double xx = s_hat/(CM_energy*CM_energy);
    // std::cout << "s_hat = " << s_hat << " i.e. between" << QQ_Min << " and " << CM_energy*CM_energy << std::endl;
    // std::cout << "xx = " << xx << std::endl;
    // std::cout << "randsjacob = " << randsjacob << std::endl;
    randsjacob = randsjacob*CM_energy*CM_energy;
    // std::cout << "randsjacob after s_hat = " << randsjacob << std::endl;

    //Next set q^2 in range 0 to s_hat
    double QQ_Max = std::min(drellyan_in->res_1.QQ_Max,std::sqrt(s_hat));
    if (drellyan_in->DYnarrowwidthapprox!=1) {
      q2 = std::pow(drellyan_in->res_1.QQ_Min,2) + x[1]*(std::pow(QQ_Max,2) - std::pow(QQ_Min,2));
      randsjacob = randsjacob*(std::pow(QQ_Max,2)-std::pow(QQ_Min,2));
    }
    // std::cout << "q2 = " << q2 << " i.e. between " << QQ_Min*QQ_Min << " and " << QQ_Max*QQ_Max << std::endl;
    // std::cout << "randsjacob after q2 = " << randsjacob << std::endl;
    // else if (drellyan_in->DYprocess == 4 || drellyan_in->DYprocess == 5) {
    //   q2 = drellyan_in->mz*drellyan_in->mz;
    //   randsjacob = randsjacob*pi*drellyan_in->mz*drellyan_in->zw;
    // }
    // else if (drellyan_in->DYprocess == 1 || drellyan_in->DYprocess == 2 || drellyan_in->DYprocess == 3) {
    //   q2 = drellyan_in->mw*drellyan_in->mw;
    //   randsjacob = randsjacob*pi*drellyan_in->mw*drellyan_in->ww;
    // }

    //Next set qt^2, can be set anywhere in range 0 to (s_hat-q^2)/(2*sqrt(s_hat))
    double QT_Min = std::max(0., drellyan_in->res_1.QT_Min);
    // std::cout << "drellyan_in->res_1.QT_Min = " << drellyan_in->res_1.QT_Min << std::endl;
    // std::cout << "QT_Min = " << QT_Min << std::endl;
    double QT_kinmax = (s_hat-q2)/(2*std::sqrt(s_hat));
    // std::cout << "QT_kinmax = " << QT_kinmax << std::endl;
    // std::cout << "QT_kinmax^2 = " << QT_kinmax*QT_kinmax << std::endl;
    double QT_Max = std::min(drellyan_in->res_1.QT_Max, QT_kinmax);
    // std::cout << "QT_Max = " << QT_Max << std::endl;
    if (QT_Min > QT_Max) {
      randsjacob = 0;
      //      qt2 = 100; //large to avoid any nan issues in events output, contribution will be zero anyway as randsjacob = 0 set and multiplied from here on
      // std::cout << "QT_Min = " << QT_Min << std::endl;
      // std::cout << "QT_Max = " << QT_Max << std::endl;
    }
    else {
      qt2 = std::pow(QT_Min,2) + x[2]*(std::pow(QT_Max,2) - std::pow(QT_Min,2));
      // std::cout << "qt2 = " << qt2 << " i.e. between " << QT_Min*QT_Min << " and " << QT_Max*QT_Max << std::endl;
      // std::cout << "qt2 = " << qt2 << std::endl;
      randsjacob = randsjacob*(std::pow(QT_Max,2)-std::pow(QT_Min,2));
    }
    // std::cout << "randsjacob after qt2 = " << randsjacob << std::endl;
    //Choosing qT set sinthetaqt=qT/qTmax,
    double costhetaqt = std::sqrt(std::pow(s_hat-q2,2)-4*s_hat*qt2)/(s_hat-q2); // at this stage I have chosen the positive costhetaqt, could also choice minus this value as costhetaqt=+-sqrt(1-sin^2thetaqt) and only sinthetaqt = qT/qTmax is set, choosing positive sign means that expect forward event rather than backward event in parton CM frame and ultimately therefore generally we will have x1 > x2 (although not guaranteed as boost to lab could overwhelm that effect for individual events) (will I therefore need a factor of 2?)
    // std::cout << "costhetaqt = " << costhetaqt << std::endl;


    //Now set y, y-eta_hat is the dilepton rapidity in the parton CM frame and is fixed by s_hat, q2 and qt2, however eta_hat is yet to be fixed as whilst x = x1*x2 is fixed x1/x2 is not which is what sets eta_hat
    //N.B given sign choice of costhetaqt expect x1>x2 and so eta_hat positive
    double eta_Min = drellyan_in->res_1.eta_Min, eta_Max = drellyan_in->res_1.eta_Max;
    double arg = (s_hat-q2)/(s_hat+q2)*costhetaqt;
    double yminuseta_hat = 0.5*std::log((1+arg)/(1-arg));
    // std::cout << "y-eta_hat = " << yminuseta_hat << std::endl;
    //Given y-eta_hat is fixed we know y must be between eta_hat_min + yminuseta_hat and eta_hat_max + yminuseta_hat
    double eta_hat_min = 0.5*std::log(xx);
    double eta_hat_max = -0.5*std::log(xx);
    // std::cout << "eta_hat_min = " << eta_hat_min << " eta_hat_max = " << eta_hat_max << std::endl;
    double y_Min = std::max(eta_Min,eta_hat_min+yminuseta_hat);
    double y_Max = std::min(eta_Max,eta_hat_max+yminuseta_hat);
    eta = y_Min + x[3]*(y_Max-y_Min);
    // std::cout << "eta = " << eta << " i.e. between " << y_Min << " and " << y_Max << std::endl;
    randsjacob = randsjacob*(y_Max-y_Min);
    // std::cout << "randsjacob after y = " << randsjacob << std::endl;
    thetaCM = pi*x[4]; //thetaCM between 0 and pi
    c_thetaCM = std::cos(thetaCM);
    // std::cout << "costhetaCM = " << c_thetaCM << std::endl;
    randsjacob = randsjacob*pi*std::sin(thetaCM); // this Jacobian is actually for d(costheta)
    // std::cout << "randsjacob after costhetaCM = " << randsjacob << std::endl;
//
    phiCM = 2*pi*x[5];   //phiCM between 0 and 2pi
    // std::cout << "phiCM = " << phiCM << std::endl;
    randsjacob = randsjacob*2*pi;
    // std::cout << "randsjacob after phiCM = " << randsjacob << std::endl;

// azimuthal angle of the qT vector
    double phisys = 2*pi*x[6];
    // std::cout << "phisys = " << phisys << std::endl;
    randsjacob = randsjacob*2*pi;
    // std::cout << "randsjacob after phisys = " << randsjacob << std::endl;

    //Now we have set all 7 variables we have fully defined the phase space, therefore any other variables we may wish to make the defining of momenta etc easier are already fixed in terms of our variables
    double x1 = std::sqrt(s_hat/(CM_energy*CM_energy))*std::exp(eta)*std::sqrt((1-arg)/(1+arg));
    double x2 = std::sqrt(s_hat/(CM_energy*CM_energy))*std::exp(-eta)*std::sqrt((1+arg)/(1-arg));
    // std::cout << "x1 = " << x1 << " x2 = " << x2 << std::endl;
    double eta_hat = 0.5*std::log(x1/x2);
    // std::cout << "eta_hat = " << eta_hat << std::endl;
    double MT = std::sqrt(q2+qt2);
    // std::cout << "MT = " << MT << std::endl;
    //Now kjet_plus is actually fixed by s_hat, q2, qt2 and y(called eta in code)
    double kjet_plus = std::exp(eta)*(std::sqrt(s_hat)*std::sqrt((1-arg)/(1+arg))-MT);
    //    if (kjet_plus == 0 || kjet_plus!=kjet_plus) {
    if (kjet_plus!=kjet_plus) {
      std::cout << "kjet_plus NAN ISSUE!" << std::endl;
      std::cout << "kjet_plus = " << kjet_plus << std::endl;
      std::cout << "arg = " << arg << std::endl;
      std::cout << "sqrt(s_hat) = " << std::sqrt(s_hat) << std::endl;
      std::cout << "std::sqrt(s_hat)*std::sqrt((1-arg)/(1+arg)) = " << std::sqrt(s_hat)*std::sqrt((1-arg)/(1+arg)) << " MT = " << MT << std::endl;
      std::cout << "Now kjet_plus should be in range " << qt2/(std::sqrt(s_hat)-MT*std::exp(-eta)) << " to " << std::sqrt(s_hat)-MT*std::exp(eta) << std::endl;
      std::cout << "for this point we have: " << std::endl;
      std::cout << "q2 = " << q2 << " qt2 = " << qt2 << " eta = " << eta << " thetaCM = " << thetaCM << " phiCM = " << phiCM << " phisys = " << phisys << std::endl;
    std::cout << "alternative kjet_plus form = " << qt2*std::exp(eta)/(std::sqrt(s_hat)*std::sqrt((1+arg)/(1-arg))-MT) << std::endl;
    std::cout << "alternative kjet_plus form 2 = " << x1*CM_energy-MT*std::exp(eta) << std::endl;
    std::cout << "alternative kjet_plus form 3 = " << qt2/(x2*CM_energy-MT*std::exp(-eta)) << std::endl;
    }

    //std::cout << "kjet_plus = " << kjet_plus << std::endl;
    if (kjet_plus < 0.00000000001 && randsjacob > 0) {
      std::cout << "kjet_plus = " << kjet_plus << std::endl;
      std::cout << "eta = " << eta << std::endl;
      std::cout << "qt2 = " << qt2 << std::endl;
      std::cout << "sqrt(s_hat) = " << std::sqrt(s_hat) << std::endl;
      std::cout << "arg = " << arg << std::endl;
      std::cout << "MT = " << MT << std::endl;
      std::cout << "exp(eta) = " << std::exp(eta) << std::endl;
      std::cout << "std::sqrt(s_hat)*std::sqrt((1-arg)/(1+arg))-MT = " << std::sqrt(s_hat)*std::sqrt((1-arg)/(1+arg))-MT << std::endl;
      std::cout << "std::sqrt((1-arg)/(1+arg))-MT) = " << std::sqrt((1-arg)/(1+arg))-MT << std::endl;
      std::cout << "kjet_plus = " << kjet_plus << std::endl;
      std::cout << "randsjacob here= " << randsjacob << std::endl;
    }
    // std::cout << "kjet_plus = " << kjet_plus << std::endl;
    // std::cout << "alternative kjet_plus form = " << qt2*std::exp(eta)/(std::sqrt(s_hat)*std::sqrt((1+arg)/(1-arg))-MT) << std::endl;
    // std::cout << "alternative kjet_plus form 2 = " << x1*CM_energy-MT*std::exp(eta) << std::endl;
    // std::cout << "alternative kjet_plus form 3 = " << qt2/(x2*CM_energy-MT*std::exp(-eta)) << std::endl;
    //double kjet_plus;
    // kjet_plus_Min = qt2 / ( CM_energy - std::sqrt(q2 + qt2)*std::exp(-eta) );
    // kjet_plus_Max = CM_energy - std::sqrt(q2 + qt2)*std::exp(eta);
    // kjet_plus = kjet_plus_Min + x[6]*(kjet_plus_Max - kjet_plus_Min);
    // randsjacob = randsjacob*(kjet_plus_Max - kjet_plus_Min);

// Total momentum of dilepton system
    four_momentum qvec(4);
    // double MT = std::sqrt(q2 + qt2);
    qvec[0] = MT*std::cosh(eta); qvec[1] = std::sqrt(qt2)*std::cos(phisys);
    qvec[2] = std::sqrt(qt2)*std::sin(phisys); qvec[3] = MT*std::sinh(eta);

    // std::cout << "qvec = " << qvec[0] << " " << qvec[1] << " " << qvec[2] << " " << qvec[3] << std::endl;
    // double kjet_plus1 = -std::exp(y)*(2*qt2+q2-s_hat
// PS setting: initial state
    //TCRIDGE NEED TO COMMENT OUT NEXT TWO LINES FOR INPUTTING MOMENTA BY HAND SEP 20
       four_momentum k1(4), k2(4), kjet(4);
       double kjet_plus_1, kjet_plus_2, kjet_plus_3;
       kjet_plus_1 = qt2*std::exp(eta)/(std::sqrt(s_hat)*std::sqrt((1+arg)/(1-arg))-MT);
       kjet_plus_2 = x1*CM_energy-MT*std::exp(eta);
       kjet_plus_3 = qt2/(x2*CM_energy-MT*std::exp(-eta));
       // std::cout << "kjet_plus_1 = " << kjet_plus_1 << std::endl;
       // std::cout << "kjet_plus_2 = " << kjet_plus_2 << std::endl;
       // std::cout << "kjet_plus_3 = " << kjet_plus_3 << std::endl;
       if (kjet_plus == 0 || kjet_plus!=kjet_plus) {
       //       if (kjet_plus!=kjet_plus) {
	 double kjet_plus_temp = 0;
	 kjet_plus_temp = std::max(kjet_plus_1,kjet_plus_2);
	 kjet_plus_temp = std::max(kjet_plus_temp, kjet_plus_3);
	 kjet_plus = kjet_plus_temp;
	 // if (kjet_plus_1!=0) {
	 //   kjet_plus = kjet_plus_1;
	 // }
	 // else if (kjet_plus_2!=0) {
	 //   kjet_plus = kjet_plus_2;
	 // }
	 // else if (kjet_plus_3!=0) {
	 //   kjet_plus = kjet_plus_3;
	 // }
	 // else {
	 //   std::cout << "kjet_plus = 0 in all forms!" << std::endl;
	 // }
	 // std::cout << "kjet_plus_1 = " << kjet_plus_1 << std::endl;
	 // std::cout << "kjet_plus_2 = " << kjet_plus_2 << std::endl;
	 // std::cout << "kjet_plus_3 = " << kjet_plus_3 << std::endl;
	 // std::cout << "kjet_plus now = " << kjet_plus << std::endl;
       }


       in_state_plus_jet(qvec, CM_energy, kjet_plus, k1, k2, kjet);
       if (k1 !=k1 || k2!=k2 || kjet!=kjet) {
	 // std::cout << "qvec = " << qvec[0] << " " << qvec[1] << " " << qvec[2] << " " << qvec[3] << std::endl;
	 // std::cout << "CM_energy = " << CM_energy << std::endl;
	 // std::cout << "kjet_plus = " << kjet_plus << std::endl;
	 // std::cout << "k1 = " << k1[0] << " " << k1[1] << " " << k1[2] << " " << k1[3] << std::endl;
	 // std::cout << "k2 = " << k2[0] << " " << k2[1] << " " << k2[2] << " " << k2[3] << std::endl;
	 // std::cout << "kjet = " << kjet[0] << " " << kjet[1] << " " << kjet[2] << " " << kjet[3] << std::endl;
       }
//     //TCRIDGE SEP 20 FOR READING IN MOMENTA FOR TESTING
//    PSpoint pp1_MCFM;
//    four_momentum k1_MCFM(4),k2_MCFM(4),kjet_MCFM(4), qvec_MCFM(4), pp1_MCFM_mom1(4), pp1_MCFM_mom2(4);
// //     //FIRST POINT IN 8TeV finite piece run for reSolve
// //     // k1_MCFM[1] = 0;
// //     // k1_MCFM[2] = 0;
// //     // k1_MCFM[3] = 3380.97;
// //     // k1_MCFM[0] = 3380.97;
// //     // k2_MCFM[1] = 0;
// //     // k2_MCFM[2] = 0;
// //     // k2_MCFM[3] = -2466.61;
// //     // k2_MCFM[0] = 2466.61;
// //     // kjet_MCFM[1] = 11.7355;
// //     // kjet_MCFM[2] = -52.7137;
// //     // kjet_MCFM[3] = -2465.98;
// //     // kjet_MCFM[0] = 2466.57;
// //     // // std::cout << "here 1... " << std::endl;
// //     // pp1_MCFM_mom1[1] = -10.4336;
// //     // pp1_MCFM_mom1[2] = 34.0296;
// //     // pp1_MCFM_mom1[3] = 2998.27;
// //     // pp1_MCFM_mom1[0] = 2998.48;
// //     // pp1_MCFM_mom2[1] = -1.30194;
// //     // pp1_MCFM_mom2[2] = 18.6842;
// //     // pp1_MCFM_mom2[3] = 382.071;
// //     // pp1_MCFM_mom2[0] = 382.53;

// //     //FIRST POINT OF MCFM Zgam_NLO_for_reSolve_real_vshort run
// //     //NOTE TO GO FROM MCFM P TO THESE YOU NEED MCFM (x,y,z,E) -> reSolve (3,-2,1,0)
//    k1_MCFM[3] = 0;
//    k1_MCFM[2] = 0;
//    k1_MCFM[1] = -502.27662862141653;
//    k1_MCFM[0] = -502.27662862141653;
//    k2_MCFM[3] = 0;
//    k2_MCFM[2] = 0;
//    k2_MCFM[1] = 188.56891264578434;
//    k2_MCFM[0] = -188.56891264578434;
//    pp1_MCFM_mom1[3] = -91.701893323247063;
//    pp1_MCFM_mom1[2] = -47.009765785378235;
//    pp1_MCFM_mom1[1] = 19.158571234597737;
//    pp1_MCFM_mom1[0] = 104.81510468446550;
//    // std::cout << "here 1... " << std::endl;
//    pp1_MCFM_mom2[3] = -148.63221580738792;
//    pp1_MCFM_mom2[2] = 0.72611346519357767;
//    pp1_MCFM_mom2[1] = -55.710311799636919;
//    pp1_MCFM_mom2[0] = 158.73153957985488;
//    kjet_MCFM[3] = 240.33410913063500;
//    kjet_MCFM[2] = 46.283652320184657;
//    kjet_MCFM[1] = 350.25945654067141;
//    kjet_MCFM[0] = 427.29889700288044;
//    qvec_MCFM[0] = pp1_MCFM_mom1[0] + pp1_MCFM_mom2[0];
//    qvec_MCFM[1] = pp1_MCFM_mom1[1] + pp1_MCFM_mom2[1];
//    qvec_MCFM[2] = pp1_MCFM_mom1[2] + pp1_MCFM_mom2[2];
//    qvec_MCFM[3] = pp1_MCFM_mom1[3] + pp1_MCFM_mom2[3];
//    // std::cout << "here 2... " << std::endl;
//    // pp1_MCFM.mom(0)[1] = 0;
//    // pp1_MCFM.mom(0)[3] = 0;
//    // pp1_MCFM.mom(0)[4] = -1.9194959036398892E-002;
//    // pp1_MCFM.mom(0)[0] = -1.9194959036398892E-002;
//    // pp1_MCFM.mom(1)[1] = 0;
//    // pp1_MCFM.mom(1)[2] = 0;
//    // pp1_MCFM.mom(1)[3] = 242.64459090955395;
//    // pp1_MCFM.mom(1)[0] = -242.64459090955395;
//    // pp1_MCFM.mom(2)[1] = 0.15821094057113416;
//    // pp1_MCFM.mom(2)[2] = -1.5126752823951315;
//    // pp1_MCFM.mom(2)[3] = -52.380408917549317;
//    // pp1_MCFM.mom(2)[0] = 52.402485204247377;
//    // pp1_MCFM.mom(3)[1] = 0.58473788783420766;
//    // pp1_MCFM.mom(3)[2] = 1.3969767846953922;
//    // pp1_MCFM.mom(3)[3] = -159.22643989316549;
//    // pp1_MCFM.mom(3)[0] = 159.23364162006183;
//    // pp1_MCFM.mom(4)[1] = -0.74294882840534182;
//    // pp1_MCFM.mom(4)[2] = 0.11569849769973921;
//    // pp1_MCFM.mom(4)[3] = -31.018547139802735;
//    // pp1_MCFM.mom(4)[0] = 31.027659044281165;
//    amp1.add_mom(k1_MCFM);
//    amp1.add_mom(k2_MCFM);
//    // std::cout << "here 3... " << std::endl;
//    amp1.add_mom(pp1_MCFM_mom1);
//    amp1.add_mom(pp1_MCFM_mom2);
//    //Add the jet
//    amp1.add_mom(kjet_MCFM);
//    // std::cout << "here 4... " << std::endl;
//    //Finalise PS
//    amp1.set_products();

//    resuPS_0.set(drellyan_in->res_1, q2, eta, qt2);

//    qvec_MCFM[0] = pp1_MCFM_mom1[0] + pp1_MCFM_mom2[0];
//    qvec_MCFM[1] = pp1_MCFM_mom1[1] + pp1_MCFM_mom2[1];
//    qvec_MCFM[2] = pp1_MCFM_mom1[2] + pp1_MCFM_mom2[2];
//    qvec_MCFM[3] = pp1_MCFM_mom1[3] + pp1_MCFM_mom2[3];
//    // std::cout << "here 2... " << std::endl;
//    // pp1_MCFM.mom(0)[1] = 0;
//    // pp1_MCFM.mom(0)[3] = 0;
//    // pp1_MCFM.mom(0)[4] = -1.9194959036398892E-002;
//    // pp1_MCFM.mom(0)[0] = -1.9194959036398892E-002;
//    // pp1_MCFM.mom(1)[1] = 0;
//    // pp1_MCFM.mom(1)[2] = 0;
//    // pp1_MCFM.mom(1)[3] = 242.64459090955395;
//    // pp1_MCFM.mom(1)[0] = -242.64459090955395;
//    // pp1_MCFM.mom(2)[1] = 0.15821094057113416;
//    // pp1_MCFM.mom(2)[2] = -1.5126752823951315;
//    // pp1_MCFM.mom(2)[3] = -52.380408917549317;
//    // pp1_MCFM.mom(2)[0] = 52.402485204247377;
//    // pp1_MCFM.mom(3)[1] = 0.58473788783420766;
//    // pp1_MCFM.mom(3)[2] = 1.3969767846953922;
//    // pp1_MCFM.mom(3)[3] = -159.22643989316549;
//    // pp1_MCFM.mom(3)[0] = 159.23364162006183;
//    // pp1_MCFM.mom(4)[1] = -0.74294882840534182;
//    // pp1_MCFM.mom(4)[2] = 0.11569849769973921;
//    // pp1_MCFM.mom(4)[3] = -31.018547139802735;
//    // pp1_MCFM.mom(4)[0] = 31.027659044281165;
//    amp1.add_mom(k1_MCFM);
//    amp1.add_mom(k2_MCFM);
//    // std::cout << "here 3... " << std::endl;
//    amp1.add_mom(pp1_MCFM_mom1);
//    amp1.add_mom(pp1_MCFM_mom2);
//    //Add the jet
//    amp1.add_mom(kjet_MCFM);
//    // std::cout << "here 4... " << std::endl;
// //Finalise PS
//    amp1.set_products();

//    resuPS_0.set(drellyan_in->res_1, q2, eta, qt2);
//    // std::cout << "q2 in ps = " << q2 << std::endl;
//   // std::cout << "here 5... " << std::endl;
   //UP TO HERE IS THE SETTING OF MOMENTA BY HAND SEP 20 - NEED TO UNCOMMENT BELOW 10 LINES TO RESET TO NORMAL


    //TCRIDGE COMMENT OUT THESE NEXT APPROX 10 LINES AS INPUT MOMENTA BY HAND - SEP20
     amp1.add_mom(k1);
     amp1.add_mom(k2);
     //PS setting: final state (not counting the jet)
     PSpoint pp1;
     set_PS_twobody(c_thetaCM, phiCM, qvec, 0., 0., pp1);


     amp1.add_mom(pp1.mom(0));
     amp1.add_mom(pp1.mom(1));
 // Add the jet
     amp1.add_mom(kjet);
 // Finalise PS
     amp1.set_products();
 //
     resuPS_0.set(drellyan_in->res_1, q2, eta, qt2);



// Final: insert the Jacobian factor for the full PS measure as a function of the randoms-related physical variables, that is:
//
//  dx1 dx2 dP_3 = Jac * dq2 deta dqt2 dc_thetaCM dphiCM dphisys dkjet_plus
//
// where dP_3 is the PS of 3 massless particles with total 4-momentum
// x1 P1 + x2 P2, with P1 and P2 the hadron momenta.

    // double s_hat = q2 + std::exp(-eta)*MT*kjet_plus + qt2*(2. + std::exp(-eta)*MT/kjet_plus);
    //Commented out as s_hat now generated as one of MC variables
    // double s_hat = q2 + std::exp(-eta)*MT*kjet_plus + qt2*(2. + std::exp(eta)*MT/kjet_plus);
    ss_hat = s_hat; //load into ss_hat so can pass into xsection later to read correct x1 and x2 values for PDFs
    // std::cout << "set ss_hat to " << ss_hat << std::endl;
    // std::cout << "s_hat = " << s_hat << std::endl;
    // double eta_hat = 1./2.*std::log( (kjet_plus*kjet_plus + std::exp(eta)*MT*kjet_plus)/(q2 + std::exp(-eta)*MT*kjet_plus) );
    // double eta_hat = 1./2.*std::log( (kjet_plus*kjet_plus + std::exp(eta)*MT*kjet_plus)/(qt2 + std::exp(-eta)*MT*kjet_plus) );
    etaa_hat = eta_hat; //needed later to call the PDFs at the right x1 and x2 in xsection_nlojet
    // std::cout << "set etaa_hat to " << etaa_hat << std::endl;
    // std::cout << "eta_hat = " << eta_hat << std::endl;
    // double c_theta_qt = MT*std::sinh(eta - eta_hat)/std::sqrt(qt2 + MT*MT*std::pow(std::sinh(eta - eta_hat),2));
    // std::cout << "eta = " << eta << " eta_hat = " << eta_hat << std::endl;
    // std::cout << "x = " << q2/pow(CM_energy,2) << " x_hat = " << s_hat/pow(CM_energy,2) << std::endl;
    // std::cout << "c_theta_qt = " << c_theta_qt << std::endl;
    // std::cout << "MT = "<< MT << std::endl;
    // std::cout << "sinh(eta - eta_hat) = "<< std::sinh(eta - eta_hat) << std::endl;
    // std::cout << "eta = " << eta << " eta_hat = " << eta_hat << std::endl;
    // std::cout << "sqrt(qt2 + MT*MT*std::sinh(eta - eta_hat)) = " << sqrt(qt2 + MT*MT*std::sinh(eta - eta_hat)) << std::endl;
    // std::cout << "qt2 = " << qt2 << " MT*MT*std::sinh(eta - eta_hat) = " << MT*MT*std::sinh(eta - eta_hat) << std::endl;
    // double Jac = std::exp(-eta_hat)*s_hat/(512.*std::pow(pi,5))/(1.-c_theta_qt)/std::sqrt(s_hat)/(s_hat - q2)/CM_energy/CM_energy;
    //New Jacobian
    // double Atilde = -qt2*MT*(std::cosh(eta-eta_hat))/std::pow(qt2+MT*MT*std::pow(std::sinh(eta-eta_hat),2),1.5)*0.5*(1/kjet_plus+1/(kjet_plus+std::exp(eta)*MT) - MT/(qt2*std::exp(eta)+kjet_plus*MT)); //jacobian for dcosthetaqt -> dkjet_plus change
    // double Atilde = qt2*MT*(std::cosh(eta-eta_hat))/std::pow(qt2+MT*MT*std::pow(std::sinh(eta-eta_hat),2),1.5)*0.5*(1/kjet_plus+1/(kjet_plus+std::exp(eta)*MT) - MT/(qt2*std::exp(eta)+kjet_plus*MT)); //jacobian for dcosthetaqt -> dkjet_plus change, sign unimportant?
    // double Btilde = 1/CM_energy*MT*(std::exp(eta)/kjet_plus+1/MT); //jacobian for dx1dx2 -> dqT2 dy change
    // double Jac = 1/(512*std::pow(pi,5)*std::sqrt(q2)*std::sqrt(s_hat))*Btilde*Atilde;
    // std::cout << "randsjacob Jac piece = " << Jac << std::endl;
    // double dx1dx2dcosthetaqt_dqT2dydkjetplus = (std::exp(-eta)*(qt2*(kjet_plus*kjet_plus+std::exp(2*eta)*qt2+2*std::exp(eta)*kjet_plus*MT)*std::cosh(eta-eta_hat)-q2*(kjet_plus*kjet_plus-std::exp(2*eta)*qt2)*std::sinh(eta-eta_hat)))/(2*kjet_plus*kjet_plus*CM_energy*CM_energy*std::pow(qt2+std::pow(MT,2)*std::sinh(eta-eta_hat)*std::sinh(eta-eta_hat),1.5));
    // std::cout << "dx1dx2dcosthetaqt_dqT2dydkjetplus = " << dx1dx2dcosthetaqt_dqT2dydkjetplus << std::endl;
    // std::cout << "num = " << (std::exp(-eta)*(qt2*(kjet_plus*kjet_plus+std::exp(2*eta)*qt2+2*std::exp(eta)*kjet_plus*MT)*std::cosh(eta-eta_hat)-q2*(kjet_plus*kjet_plus-std::exp(2*eta)*qt2)*std::sinh(eta-eta_hat))) << " denom = " << (2*kjet_plus*kjet_plus*CM_energy*CM_energy*std::pow(qt2+std::pow(MT,2)*std::sinh(eta-eta_hat)*std::sinh(eta-eta_hat),1.5)) << std::endl;
    // double Jac = 1/(512*std::pow(pi,5)*std::sqrt(q2)*std::sqrt(s_hat)*dx1dx2dcosthetaqt_dqT2dydkjetplus);
    double dx1dx2dcosthetaqt_dshatdqT2dy = std::abs(2*s_hat/(CM_energy*CM_energy*std::pow(s_hat-q2,2)*costhetaqt));
    // std::cout << "extra Jac piece from x1, x2, costhetaqt ->s_hat, qt2, y = " << dx1dx2dcosthetaqt_dshatdqT2dy << std::endl;
    double Jac = 1/(512*std::pow(pi,5)*std::sqrt(q2)*std::sqrt(s_hat))*dx1dx2dcosthetaqt_dshatdqT2dy;
    // std::cout << "extra total Jac = " << Jac << std::endl;

    // double Jac = (1+c_theta_qt)*(CM_energy*CM_energy)*(q2-s_hat)*(q2-s_hat)*(-q2+s_hat)*(q2+s_hat)/(16384*pow(pi,5)*pow(s_hat,4));
    // std::cout << "Atilde = " << Atilde << " Btilde = " << Btilde << " rest = " << 1/(512*std::pow(pi,5)*std::sqrt(q2)*std::sqrt(s_hat)) <<  std::endl;

    // std::cout << "randsjacob = " << randsjacob << std::endl;
    // std::cout << "randsjacob Jac piece = " << Jac << std::endl;
    // std::cout << "Jac = " << Jac << std::endl;
    randsjacob *= Jac;
    //finally multiply by |k| and |q1| which are the norms of the 3-momenta of the jet in the parton CM frame and that of either lepton in the Z frame, which are given by Kallen functions, and add a factor of 2 for fact I have chosen sign of costhetaqt
    // std::cout << "randsjacob after Jac piece = " << randsjacob << std::endl;
    // std::cout << "mz = " << drellyan_in->mz << std::endl;
    randsjacob *= 2*(s_hat-q2)/(2*std::sqrt(s_hat))*std::sqrt(q2)/2;
    // std::cout << "randsjacob |k||q1| and factor of 2 piece = " << 2*(s_hat-q2)/(2*std::sqrt(s_hat))*std::sqrt(q2)/2 << std::endl;
    // std::cout << "randsjacob |k||q1| piece and factor of 2piece /q2 = " << 2*(s_hat-q2)/(2*std::sqrt(s_hat))*std::sqrt(q2)/2*1/q2 << std::endl;
    // std::cout << "randsjacob final = "<< randsjacob << std::endl;
    //randsjacob = randsjacob/(Jac*(kjet_plus_Max - kjet_plus_Min));
    if (randsjacob!=randsjacob) {
      std::cout << "randsjacob nan issue= " << randsjacob << std::endl;
      std::cout << "set randsjacob 0 for now to avoid issues if it comes as nan...." << std::endl;
      randsjacob = 0; //nan catcher!
      std::cout << "Point had s_hat, q2, qt2, y, thetaCM, phiCM, phisys = " << s_hat << " " << q2 << " " << qt2 << " " << eta << " " << thetaCM << " " << phiCM << " " << phisys << std::endl;
      //Mainly catches cases where it had QT_Min > QT_Max which have randsjacob explicitly set to 0 anyway
    }

    //TCRIDGE OCT20 TEMP
    //randsjacob = randsjacob*std::sqrt(qt2);
    //std::cout << "s_hat = " << s_hat << " qt = " << std::sqrt(qt2) << std::endl;

    return randsjacob;

}


void in_state_plus_jet(four_momentum qvec, double CM_energy, double kjet_plus, four_momentum& k1, four_momentum& k2, four_momentum& kjet){

// Transverse parts
  kjet[1] = -qvec[1]; kjet[2] = -qvec[2];
  double kjetT2 = kjet[1]*kjet[1] + kjet[2]*kjet[2];
  k1[1] = 0.; k1[2] = 0.;
  k2[1] = 0.; k2[2] = 0.;

// Longitudinal parts and energies
  double kjet_minus = kjetT2/kjet_plus;
  kjet[0] = (kjet_plus + kjet_minus)/2.;
  kjet[3] = (kjet_plus - kjet_minus)/2.;
  // k1[0] = (kjet_plus)/2.; k1[3] = (kjet_plus)/2.;
  // k2[0] = (kjet_minus)/2.; k2[3] = -(kjet_minus)/2.;
  k1[0] = (kjet_plus+qvec[0]+qvec[3])/2.; k1[3] = (kjet_plus+qvec[0]+qvec[3])/2.;
  k2[0] = (kjet_minus+qvec[0]-qvec[3])/2.; k2[3] = -(kjet_minus+qvec[0]-qvec[3])/2.;

  if(k2[0]!=k2[0] || k2[0] > 1000000000000) {
    std::cout << "kjet_plus = " << kjet_plus << " kjet_minus = " <<kjet_minus << " kjetT2 = " << kjetT2 << std::endl;
  }

}
