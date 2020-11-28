
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

    double QQ_Max = std::max(0., drellyan_in->res_1.QQ_Max);
    QQ_Max = std::max(QQ_Max,CM_energy*CM_energy);
//
    if (drellyan_in->DYnarrowwidthapprox!=1) {
      q2 = std::pow(drellyan_in->res_1.QQ_Min,2) + x[0]*(std::pow(drellyan_in->res_1.QQ_Max,2) - std::pow(drellyan_in->res_1.QQ_Min,2));
      randsjacob = randsjacob*(std::pow(drellyan_in->res_1.QQ_Max,2)-std::pow(drellyan_in->res_1.QQ_Min,2));
    }
    else if (drellyan_in->DYprocess == 4 || drellyan_in->DYprocess == 5) {
      q2 = drellyan_in->mz*drellyan_in->mz;
      randsjacob = randsjacob*pi*drellyan_in->mz*drellyan_in->zw;
    }
    else if (drellyan_in->DYprocess == 1 || drellyan_in->DYprocess == 2 || drellyan_in->DYprocess == 3) {
      q2 = drellyan_in->mw*drellyan_in->mw;
      randsjacob = randsjacob*pi*drellyan_in->mw*drellyan_in->ww;
    }
    // std::cout << "randoms = " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << std::endl;
    // std::cout << "randsjacob qq = " << randsjacob << std::endl;
    // std::cout << "CM_energy = " << CM_energy << std::endl;
    // std::cout << "QQ_Max^2 = " << drellyan_in->res_1.QQ_Max*drellyan_in->res_1.QQ_Max << std::endl;
    //
    double eta_Min = drellyan_in->res_1.eta_Min, eta_Max = drellyan_in->res_1.eta_Max;
    //TCRIDGE need etalim all the time now as otherwise kjet_plus_Max can go negative and cause problems of nans etc in eta_hat and beyond
    // if(drellyan_in->res_1.auto_etalim){
      eta_Max = 1./2. * std::log(CM_energy*CM_energy/q2);
      // std::cout << "eta_Max = " << eta_Max << std::endl;
      eta_Min = -eta_Max;
    // }
    eta = eta_Min + x[2]*(eta_Max - eta_Min);
    randsjacob = randsjacob*(eta_Max - eta_Min);
    // std::cout << "eta_lim = " << 1./2. * std::log(CM_energy*CM_energy/q2) << std::endl;

    // std::cout << "eta_Max = " << eta_Max << std::endl;
	  
    // std::cout << "randsjacob eta = " << randsjacob << std::endl;
//
    double QT_Min = std::max(0., drellyan_in->res_1.QT_Min);
    double QT_Max = std::max(0., drellyan_in->res_1.QT_Max);
    double QT_kinmax = std::sqrt( std::pow( CM_energy*CM_energy + q2, 2) / (4.*CM_energy*CM_energy*std::cosh(eta)*std::cosh(eta)) - q2);
    // std::cout << "QT_Max = " << QT_Max << " QT_kinmax = " << QT_kinmax << std::endl;
    // std::cout << "std::pow( CM_energy*CM_energy + q2, 2)= " << std::pow( CM_energy*CM_energy + q2, 2) << std::endl;
    // std::cout << "((4.*CM_energy*CM_energy*std::cosh(eta)*std::cosh(eta)) - q2 ) = " << ((4.*CM_energy*CM_energy*std::cosh(eta)*std::cosh(eta)) - q2 ) << std::endl;
    // std::cout << "4.*CM_energy*CM_energy*std::cosh(eta)*std::cosh(eta) = " << 4.*CM_energy*CM_energy*std::cosh(eta)*std::cosh(eta) << " q2 = " << q2 << std::endl;
    QT_Max = std::min(QT_Max, QT_kinmax);
    // std::cout << "QT_Max = " << QT_Max << std::endl;
    qt2 = std::pow(QT_Min,2) + x[1]*(std::pow(QT_Max,2) - std::pow(QT_Min,2));
    randsjacob = randsjacob*(std::pow(QT_Max,2)-std::pow(QT_Min,2));
    // std::cout << "randsjacob qt = " << randsjacob << std::endl;
    // std::cout << "randsjacob qt piece = " << (std::pow(QT_Max,2)-std::pow(QT_Min,2)) << std::endl;
//
    thetaCM = pi*x[3]; //thetaCM between 0 and pi
    c_thetaCM = std::cos(thetaCM);
    randsjacob = randsjacob*pi*std::sin(thetaCM); // this Jacobian is actually for d(costheta)
    // std::cout << "randsjacob theta = " << randsjacob << std::endl;
//
    phiCM = 2*pi*x[4];   //phiCM between 0 and 2pi
    randsjacob = randsjacob*2*pi;
//
// azimuthal angle of the qT vector
    double phisys = 2*pi*x[5];
    randsjacob = randsjacob*2*pi;
    // std::cout << "randsjacob phi, az = " << randsjacob << std::endl;
//
    // std::cout << "q2 = " << q2 << " qt2 = " << qt2 << " eta = " << eta << " CM_energy = " << CM_energy << std::endl;
    // std::cout << "std::sqrt(q2 + qt2)*std::exp(eta) = " << std::sqrt(q2 + qt2)*std::exp(eta) << std::endl;
// Extra variable for the jet
    double kjet_plus, kjet_plus_Min, kjet_plus_Max;
    kjet_plus_Min = qt2 / ( CM_energy - std::sqrt(q2 + qt2)*std::exp(-eta) );
    kjet_plus_Max = CM_energy - std::sqrt(q2 + qt2)*std::exp(eta);
    kjet_plus = kjet_plus_Min + x[6]*(kjet_plus_Max - kjet_plus_Min);
    randsjacob = randsjacob*(kjet_plus_Max - kjet_plus_Min);
    //randsjacob = randsjacob;
    //*(kjet_plus_Max - kjet_plus_Min);
    // std::cout << "kjet_plus_Min = " << kjet_plus_Min << std::endl;
    // std::cout << "kjet_plus_Max = " << kjet_plus_Max << std::endl;
    // std::cout << "randsjacob kjet_plus piece = " << (kjet_plus_Max - kjet_plus_Min) << std::endl;
    // std::cout << "randsjacob kjet = " << randsjacob << std::endl;

// Total momentum
    four_momentum qvec(4);
    double MT = std::sqrt(q2 + qt2);
    qvec[0] = MT*std::cosh(eta); qvec[1] = std::sqrt(qt2)*std::cos(phisys);
    qvec[2] = std::sqrt(qt2)*std::sin(phisys); qvec[3] = MT*std::sinh(eta);

    // std::cout << "qvec = " << qvec[0] << " " << qvec[1] << " " << qvec[2] << " " << qvec[3] << std::endl;
// PS setting: initial state
    four_momentum k1(4), k2(4), kjet(4);
    in_state_plus_jet(qvec, CM_energy, kjet_plus, k1, k2, kjet);
    // std::cout << "k1 = " << k1[0] << " " << k1[1] << " " << k1[2] << " " << k1[3] << std::endl;
    // std::cout << "k2 = " << k2[0] << " " << k2[1] << " " << k2[2] << " " << k2[3] << std::endl;
    // std::cout << "kjet = " << kjet[0] << " " << kjet[1] << " " << kjet[2] << " " << kjet[3] << std::endl;
    amp1.add_mom(k1);
    amp1.add_mom(k2);
// PS setting: final state (not counting the jet)
    PSpoint pp1;
    set_PS_twobody(c_thetaCM, phiCM, qvec, 0., 0., pp1);
    // std::cout << "pp1.mom(0) = " << pp1.mom(0)[0] << " " << pp1.mom(0)[1] << " " << pp1.mom(0)[2] << " " << pp1.mom(0)[3] << std::endl;
    // std::cout << "pp1.mom(1) = " << pp1.mom(1)[0] << " " << pp1.mom(1)[1] << " " << pp1.mom(1)[2] << " " << pp1.mom(1)[3] << std::endl;
    // std::cout << "k1.k2 = " << LorDot(k1,k2) << std::endl;
    // std::cout << "k1.pp1.mom(0) = " << LorDot(k1,pp1.mom(0)) << std::endl;
    // std::cout << "k1.pp1.mom(1) = " << LorDot(k1,pp1.mom(1)) << std::endl;
    // std::cout << "k1.kjet = " << LorDot(k1,kjet) << std::endl;
    // std::cout << "k2.pp1.mom(0) = " << LorDot(k2,pp1.mom(0)) << std::endl;
    // std::cout << "k2.pp1.mom(1) = " << LorDot(k2,pp1.mom(1)) << std::endl;
    // std::cout << "k2.kjet = " << LorDot(k2,kjet) << std::endl;
    // std::cout << "pp1.mom(0).pp1.mom(1) = " << LorDot(pp1.mom(0),pp1.mom(1)) << std::endl;
    // std::cout << "pp1.mom(0).kjet = " << LorDot(pp1.mom(0),kjet) << std::endl;
    // std::cout << "pp1.mom(1).kjet = " << LorDot(pp1.mom(1),kjet) << std::endl;  
    
    
    // std::cout << "Check momentum conservation: k1+k2 = qvec+kjet: " << std::endl;
    // std::cout << "k1+k2 = " << k1[0]+k2[0] << " " << k1[1]+k2[1] << " " << k1[2]+k2[2] << " " << k1[3]+k2[3] << std::endl;
    // std::cout << "qvec+kjet = " << qvec[0]+kjet[0] << " " << qvec[1]+kjet[1] << " " << qvec[2]+kjet[2] << " " << qvec[3]+kjet[3] << std::endl;
    //Test for overall momentum conservation of parton1 + parton2 -> jet + dilepton
    if (k1[0]+k2[0]-qvec[0]-kjet[0]>=0.00000000001 || k1[1]+k2[1]-qvec[1]-kjet[1]>=0.00000000001 || k1[2]+k2[2]-qvec[2]-kjet[2]>=0.00000000001 || k1[3]+k2[3]-qvec[3]-kjet[3]>=0.00000000001) {
      std::cout << "Momentum conservation issue!!!! " << std::endl;
      std::cout << k1[0]+k2[0]-qvec[0]-kjet[0] << " " << k1[1]+k2[1]-qvec[1]-kjet[1] << " " << k1[2]+k2[2]-qvec[2]-kjet[2] << " " << k1[3]+k2[3]-qvec[3]-kjet[3] << std::endl;
    }
    //Test for on-shellness of each of parton1, parton2, jet, dilepton, individual leptons
    // std::cout << "Testing on-shellness..." << std::endl;
    // std::cout << "parton1: " << k1[0]*k1[0]-k1[1]*k1[1]-k1[2]*k1[2]-k1[3]*k1[3] << std::endl;
    // std::cout << "parton2: " << k2[0]*k2[0]-k2[1]*k2[1]-k2[2]*k2[2]-k2[3]*k2[3] << std::endl;
    // std::cout << "jet: " << kjet[0]*kjet[0]-kjet[1]*kjet[1]-kjet[2]*kjet[2]-kjet[3]*kjet[3] << std::endl;
    // std::cout << "dilepton: " << qvec[0]*qvec[0]-qvec[1]*qvec[1]-qvec[2]*qvec[2]-qvec[3]*qvec[3] << std::endl;
    // std::cout << "lepton1: " << pp1.mom(0)[0]*pp1.mom(0)[0] - pp1.mom(0)[1]*pp1.mom(0)[1]-pp1.mom(0)[2]*pp1.mom(0)[2]-pp1.mom(0)[3]*pp1.mom(0)[3] << std::endl;
    // std::cout << "lepton2: " << pp1.mom(1)[0]*pp1.mom(1)[0] - pp1.mom(1)[1]*pp1.mom(1)[1]-pp1.mom(1)[2]*pp1.mom(1)[2]-pp1.mom(1)[3]*pp1.mom(1)[3] << std::endl;
    if (k1[0]*k1[0]-k1[1]*k1[1]-k1[2]*k1[2]-k1[3]*k1[3]>=0.0000001 || k2[0]*k2[0]-k2[1]*k2[1]-k2[2]*k2[2]-k2[3]*k2[3]>=0.0000001 || kjet[0]*kjet[0]-kjet[1]*kjet[1]-kjet[2]*kjet[2]-kjet[3]*kjet[3]>=0.0000001 || pp1.mom(0)[0]*pp1.mom(0)[0] - pp1.mom(0)[1]*pp1.mom(0)[1]-pp1.mom(0)[2]*pp1.mom(0)[2]-pp1.mom(0)[3]*pp1.mom(0)[3]>=0.0000001 || pp1.mom(1)[0]*pp1.mom(1)[0] - pp1.mom(1)[1]*pp1.mom(1)[1]-pp1.mom(1)[2]*pp1.mom(1)[2]-pp1.mom(1)[3]*pp1.mom(1)[3]>=0.0000001) {
      std::cout << "On-shellness issue!!!!! " << std::endl;
      std::cout << k1[0]*k1[0]-k1[1]*k1[1]-k1[2]*k1[2]-k1[3]*k1[3] << " " << k2[0]*k2[0]-k2[1]*k2[1]-k2[2]*k2[2]-k2[3]*k2[3] << " " << kjet[0]*kjet[0]-kjet[1]*kjet[1]-kjet[2]*kjet[2]-kjet[3]*kjet[3] << " " << pp1.mom(0)[0]*pp1.mom(0)[0] - pp1.mom(0)[1]*pp1.mom(0)[1]-pp1.mom(0)[2]*pp1.mom(0)[2]-pp1.mom(0)[3]*pp1.mom(0)[3] << " " << pp1.mom(1)[0]*pp1.mom(1)[0] - pp1.mom(1)[1]*pp1.mom(1)[1]-pp1.mom(1)[2]*pp1.mom(1)[2]-pp1.mom(1)[3]*pp1.mom(1)[3] << std::endl;
    }
    //Test to ensure total energy of each event is appropriate fraction of s_hadronic, i.e. x1*x2*s_hadronic
    double x1 = 1/CM_energy*(MT*std::exp(eta)+kjet_plus);
    double x2 = 1/CM_energy*(MT*std::exp(-eta)+qt2/kjet_plus);
    // std::cout << "CM_energy = " << CM_energy << std::endl;
    // std::cout << "x1 = " << x1 << " x2 = " << x2 << std::endl;
    // std::cout << "x1*x2 = " << x1*x2 << std::endl;
    // std::cout << "0.5(x1+x2)*CM_energy = " << 0.5*(x1+x2)*CM_energy << std::endl;
    // std::cout << "Incoming energy = " << k1[0]+k2[0] << std::endl;
    // std::cout << "Outgoing energy = " << kjet[0]+qvec[0] << std::endl;
    if ((k1[0]+k2[0] - 0.5*(x1+x2)*CM_energy)/(0.5*(x1+x2)*CM_energy) >=0.00000001 || (kjet[0]+qvec[0] - 0.5*(x1+x2)*CM_energy)/(0.5*(x1+x2)*CM_energy)>=0.00000001) {
      std::cout << "Incoming or outgoing energy not correct issue!!!!!!" << std::endl;
      std::cout << "0.5(x1+x2)*CM_energy = " << 0.5*(x1+x2)*CM_energy << std::endl;
      std::cout << "Incoming energy = " << k1[0]+k2[0] << std::endl;
      std::cout << "Outgoing energy = " << kjet[0]+qvec[0] << std::endl;
    }
    //Test to ensure total momentum squared of each event is same before and after, i.e. x1*x2*s
    // std::cout << "Total momentum squared should be x1*x2*CM_energy*CM_energy = " << x1*x2*CM_energy*CM_energy << std::endl;
    // std::cout << "Ingoing momentum squared = " << std::pow(k1[0]+k2[0],2) - std::pow(k1[1]+k2[1],2) - std::pow(k1[2]+k2[2],2) - std::pow(k1[3]+k2[3],2) << std::endl;
    // std::cout << "Outgoing momentum squared = " << std::pow(kjet[0]+qvec[0],2) - std::pow(kjet[1]+qvec[1],2) - std::pow(kjet[2]+qvec[2],2) - std::pow(kjet[3]+qvec[3],2) << std::endl;
    if ((std::pow(k1[0]+k2[0],2) - std::pow(k1[1]+k2[1],2) - std::pow(k1[2]+k2[2],2) - std::pow(k1[3]+k2[3],2) - x1*x2*CM_energy*CM_energy)/(x1*x2*CM_energy*CM_energy) >=0.00000001 ||(std::pow(kjet[0]+qvec[0],2) - std::pow(kjet[1]+qvec[1],2) - std::pow(kjet[2]+qvec[2],2) - std::pow(kjet[3]+qvec[3],2) - x1*x2*CM_energy*CM_energy)/(x1*x2*CM_energy*CM_energy)>=0.00000001) {
      std::cout << "Total momentum squared incoming or outgoing not correct issue!!!!!!!" << std::endl;
      std::cout << "Total momentum squared should be x1*x2*CM_energy*CM_energy = " << x1*x2*CM_energy*CM_energy << std::endl;
      std::cout << "Ingoing momentum squared = " << std::pow(k1[0]+k2[0],2) - std::pow(k1[1]+k2[1],2) - std::pow(k1[2]+k2[2],2) - std::pow(k1[3]+k2[3],2) << std::endl;
      std::cout << "Outgoing momentum squared = " << std::pow(kjet[0]+qvec[0],2) - std::pow(kjet[1]+qvec[1],2) - std::pow(kjet[2]+qvec[2],2) - std::pow(kjet[3]+qvec[3],2) << std::endl;
    }

    
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
    double s_hat = q2 + std::exp(-eta)*MT*kjet_plus + qt2*(2. + std::exp(eta)*MT/kjet_plus);
    ss_hat = s_hat; //load into ss_hat so can pass into xsection later to read correct x1 and x2 values for PDFs
    // std::cout << "set ss_hat to " << ss_hat << std::endl;
    // std::cout << "s_hat = " << s_hat << std::endl;
    // double eta_hat = 1./2.*std::log( (kjet_plus*kjet_plus + std::exp(eta)*MT*kjet_plus)/(q2 + std::exp(-eta)*MT*kjet_plus) );
    double eta_hat = 1./2.*std::log( (kjet_plus*kjet_plus + std::exp(eta)*MT*kjet_plus)/(qt2 + std::exp(-eta)*MT*kjet_plus) );
    etaa_hat = eta_hat; //needed later to call the PDFs at the right x1 and x2 in xsection_nlojet
    // std::cout << "set etaa_hat to " << etaa_hat << std::endl;
    // std::cout << "eta_hat = " << eta_hat << std::endl;
    double c_theta_qt = MT*std::sinh(eta - eta_hat)/std::sqrt(qt2 + MT*MT*std::pow(std::sinh(eta - eta_hat),2));
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
    double dx1dx2dcosthetaqt_dqT2dydkjetplus = (std::exp(-eta)*(qt2*(kjet_plus*kjet_plus+std::exp(2*eta)*qt2+2*std::exp(eta)*kjet_plus*MT)*std::cosh(eta-eta_hat)-q2*(kjet_plus*kjet_plus-std::exp(2*eta)*qt2)*std::sinh(eta-eta_hat)))/(2*kjet_plus*kjet_plus*CM_energy*CM_energy*std::pow(qt2+std::pow(MT,2)*std::sinh(eta-eta_hat)*std::sinh(eta-eta_hat),1.5));
    std::cout << "dx1dx2dcosthetaqt_dqT2dydkjetplus = " << dx1dx2dcosthetaqt_dqT2dydkjetplus << std::endl;
    std::cout << "num = " << (std::exp(-eta)*(qt2*(kjet_plus*kjet_plus+std::exp(2*eta)*qt2+2*std::exp(eta)*kjet_plus*MT)*std::cosh(eta-eta_hat)-q2*(kjet_plus*kjet_plus-std::exp(2*eta)*qt2)*std::sinh(eta-eta_hat))) << " denom = " << (2*kjet_plus*kjet_plus*CM_energy*CM_energy*std::pow(qt2+std::pow(MT,2)*std::sinh(eta-eta_hat)*std::sinh(eta-eta_hat),1.5)) << std::endl;
    double Jac = 1/(512*std::pow(pi,5)*std::sqrt(q2)*std::sqrt(s_hat)*dx1dx2dcosthetaqt_dqT2dydkjetplus);
      
    // double Jac = (1+c_theta_qt)*(CM_energy*CM_energy)*(q2-s_hat)*(q2-s_hat)*(-q2+s_hat)*(q2+s_hat)/(16384*pow(pi,5)*pow(s_hat,4));
    // std::cout << "Atilde = " << Atilde << " Btilde = " << Btilde << " rest = " << 1/(512*std::pow(pi,5)*std::sqrt(q2)*std::sqrt(s_hat)) <<  std::endl;
    
    std::cout << "randsjacob Jac piece = " << Jac << std::endl;
    std::cout << "randsjacob = " << randsjacob << std::endl;
    // std::cout << "Jac = " << Jac << std::endl;
    randsjacob *= Jac;
    // std::cout << "randsjacob final = "<< randsjacob << std::endl;
    //randsjacob = randsjacob/(Jac*(kjet_plus_Max - kjet_plus_Min));

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

}
