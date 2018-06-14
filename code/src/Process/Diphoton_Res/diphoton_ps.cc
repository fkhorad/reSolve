
#include "diphoton_ps.h"

#include "phase_space.h"
#include "constants.h"

#include "diphoton_input.h"

#include "resu_preproc.h"
#include "resu_PS.h"


double diphoton_ps(const double x[], ResummationInfo* resuminfo, PSpoint& pp1, resu_PS& resuPS_0){

    double pi = k_constants::pi;

//Assign q2, qt2, eta, thetaCM, phiCM with randoms from vegas
    double q2, qt2, eta, thetaCM, phiCM;
    double randsjacob = 1.;
    double CM_energy = resuminfo->CM_energy;
//
    q2 = std::pow(resuminfo->QQ_Min,2) + x[0]*(std::pow(resuminfo->QQ_Max,2) - std::pow(resuminfo->QQ_Min,2));
    randsjacob = randsjacob*(std::pow(resuminfo->QQ_Max,2)-std::pow(resuminfo->QQ_Min,2));
    // std::cout << "QQ^2 randsjacob piece: " << (std::pow(resuminfo->QQ_Max,2)-std::pow(resuminfo->QQ_Min,2)) << std::endl;
//
    double eta_Min = resuminfo->eta_Min, eta_Max = resuminfo->eta_Max;
    if(resuminfo->auto_etalim){
      eta_Max = 1./2. * std::log(CM_energy*CM_energy/q2);
      eta_Min = -eta_Max;
    }
    eta = eta_Min + x[2]*(eta_Max - eta_Min);
    // std::cout << "eta = " << eta << std::endl;
    randsjacob = randsjacob*(eta_Max - eta_Min);
    // std::cout << "eta randsjacob piece: " << (eta_Max - eta_Min) << std::endl;
//
    if (resuminfo->resum_flag != 0) {
      double QT_Min = resuminfo->QT_Min, QT_Max = resuminfo->QT_Max;
      if(resuminfo->auto_qtlim){
        QT_Min = 0.;
        QT_Max = std::sqrt( std::pow( CM_energy*CM_energy + q2, 2) / (4.*CM_energy*CM_energy*std::cosh(eta)*std::cosh(eta)) - q2 );
      }
      qt2 = std::pow(QT_Min,2) + x[1]*(std::pow(QT_Max,2) - std::pow(QT_Min,2));
      randsjacob = randsjacob*(std::pow(QT_Max,2)-std::pow(QT_Min,2));
      // std::cout << "qt2 randsjacob piece: " << (std::pow(QT_Max,2)-std::pow(QT_Min,2)) << std::endl;
    }
    else {
      qt2 = 0;
    }
//
    thetaCM = pi*x[3]; //thetaCM between 0 and pi
    // randsjacob = randsjacob*pi*std::sin(thetaCM); // the Jacobian is actually for d(costheta)
    // std::cout << "thetaCM randsjacob piece: " << pi*std::sin(thetaCM) << std::endl;
    randsjacob = randsjacob*pi; // the Jacobian is actually for d(costheta), for now no sinthetaCM as it is accounted for as sintheta elsewhere (diphoton_integrand.cc) for diphoton case
    // std::cout << "sin(thetaCM) = " << sin(thetaCM) << std::endl;
//
    phiCM = 2*pi*x[4];   //phiCM between 0 and 2pi
    randsjacob = randsjacob*2*pi;
    // std::cout << "phiCM randsjacob piece: " << 2*pi << std::endl;
//
    randsjacob = randsjacob/2.; // Factor of 2 for identical final state particles
    // std::cout << "Divide randsjacob by 2 as identical photons : " << 0.5 << std::endl;
// azimuthal angle of the  total final state momentum: not part of the integration
    double phisys = 2*pi*x[5];

// Total momentum
    four_momentum qvec(4);
    double MT = std::sqrt(q2 + qt2);
    // qvec[0] = MT*std::cosh(eta); qvec[1] = std::sqrt(qt2)*std::cos(phisys);
    // qvec[2] = std::sqrt(qt2)*std::sin(phisys); qvec[3] = MT*std::sinh(eta);
    //TEMP bug search
    qvec[0] = MT*std::cosh(eta); qvec[1] = std::sqrt(qt2);
    qvec[2] = 0.0; qvec[3] = MT*std::sinh(eta);

    
    // std::cout << "q2 = " << q2 << " qt2 = " << qt2 << std::endl;
    // std::cout << "MT = " << MT << " eta = " << eta << std::endl;
    // std::cout << "qvec = " << qvec[0] << " " << qvec[1] << " " << qvec[2] << " " << qvec[3] << std::endl;

// PS setting: initial state
    four_momentum k1(4), k2(4);
    in_state_w_recoil(qvec, CM_energy, 0., k1, k2);
    pp1.add_mom(k1);
    pp1.add_mom(k2);
// PS setting: final state
    set_PS_twobody(std::cos(thetaCM), phiCM, qvec, 0., 0., pp1);
    pp1.set_products();
// PS setting: reduced version for resummation
    resuPS_0.set(*resuminfo, q2, eta, qt2);

    return randsjacob;

}
