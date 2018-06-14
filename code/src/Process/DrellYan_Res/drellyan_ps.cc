
#include "drellyan_ps.h"

#include "phase_space.h"
#include "constants.h"

#include "drellyan_input.h"

#include "resu_preproc.h"
#include "resu_PS.h"

#include <fstream> // TEMP
#include <sstream> // TEMP
#include <iomanip>

//extern std::ifstream debbie;

double drellyan_ps(const double x[], drellyan_input* drellyan_in, PSpoint& pp1, resu_PS& resuPS_0){
    // NOTE PROCESS ONLY PASSED HERE TO ALLOW PROCESS = -1 TO MEAN BORN ONLY FOR TESTING
    double pi = k_constants::pi;


//Assign q2, qt2, eta, thetaCM, phiCM with randoms from vegas
    double q2, qt2, eta, thetaCM, c_thetaCM, phiCM;
    double randsjacob = 1.;
    double CM_energy = drellyan_in->res_1.CM_energy;
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
//
    double eta_Min = drellyan_in->res_1.eta_Min, eta_Max = drellyan_in->res_1.eta_Max;
    if(drellyan_in->res_1.auto_etalim){
      eta_Max = 1./2. * std::log(CM_energy*CM_energy/q2);
      eta_Min = -eta_Max;
    }
    eta = eta_Min + x[2]*(eta_Max - eta_Min);
    randsjacob = randsjacob*(eta_Max - eta_Min);
//
    if (drellyan_in->res_1.resum_flag != 0) {
      double QT_Min = drellyan_in->res_1.QT_Min, QT_Max = drellyan_in->res_1.QT_Max;
      if(drellyan_in->res_1.auto_qtlim){
        QT_Min = 0.;
        QT_Max = std::sqrt( std::pow( CM_energy*CM_energy + q2, 2) / (4.*CM_energy*CM_energy*std::cosh(eta)*std::cosh(eta)) - q2 );
      }
      qt2 = std::pow(QT_Min,2) + x[1]*(std::pow(QT_Max,2) - std::pow(QT_Min,2));
      randsjacob = randsjacob*(std::pow(QT_Max,2)-std::pow(QT_Min,2));
    }
    else {
      qt2 = 0;
    }
//
    thetaCM = pi*x[3]; //thetaCM between 0 and pi
    c_thetaCM = std::cos(thetaCM);
    randsjacob = randsjacob*pi*std::sin(thetaCM); // the Jacobian is actually for d(costheta)
//
    phiCM = 2*pi*x[4];   //phiCM between 0 and 2pi
    randsjacob = randsjacob*2*pi;
//
// azimuthal angle of the  total final state momentum: not part of the integration
    double phisys = 2*pi*x[5];

/*
// TEMP
    std::string line;
    std::stringstream ss;
    std::getline(debbie, line);
    ss << line;
    std::getline(debbie, line);
    ss << line;
    std::getline(debbie, line);
    ss << line;
    std::getline(debbie, line);
    ss << line;
    std::getline(debbie, line);
    ss << line;
    ss >> qt2;
    ss >> eta;
    ss >> phisys;
    ss >> c_thetaCM;
    ss >> phiCM;
    std::ofstream debborah; 
    debborah.open("debborah.dat", std::ios_base::app);
    debborah << std::setprecision(17) << qt2 << "\n";
    debborah << std::setprecision(17) << eta << "\n";
    debborah << std::setprecision(17) << phisys << "\n";
    debborah << std::setprecision(17) << c_thetaCM << "\n";
    debborah << std::setprecision(17) << phiCM << "\n";
    debborah << "----------\n";
    debborah.close();

    phiCM = pi/2. - phiCM;
    std::getline(debbie, line);
    std::getline(debbie, line);
    std::getline(debbie, line);
    std::getline(debbie, line);
*/

// Total momentum
    four_momentum qvec(4);
    double MT = std::sqrt(q2 + qt2);
    qvec[0] = MT*std::cosh(eta); qvec[1] = std::sqrt(qt2)*std::cos(phisys);
    qvec[2] = std::sqrt(qt2)*std::sin(phisys); qvec[3] = MT*std::sinh(eta);


// PS setting: initial state
    four_momentum k1(4), k2(4);
    in_state_w_recoil(qvec, CM_energy, 0., k1, k2);
    pp1.add_mom(k1);
    pp1.add_mom(k2);
// PS setting: final state
    set_PS_twobody(c_thetaCM, phiCM, qvec, 0., 0., pp1);
    pp1.set_products();
// PS setting: reduced version for resummation
    resuPS_0.set(drellyan_in->res_1, q2, eta, qt2);

/*
// Temp to match DYres usage
    double mur = std::sqrt(resuPS_0.mur2);
    double mures = std::sqrt(resuPS_0.mures2);
    double cmass = 2.25, bmass = 21.3444;
    setqmass_(&cmass, &bmass);
    double as, alpqf, amz;
    int nlooprun;
    if(drellyan_in->res_1.order == 0){
      amz=0.13939;
      nlooprun=1;
    }
    else if(drellyan_in->res_1.order == 1){
      amz=0.12018;
      nlooprun=2;
    }
    else if(drellyan_in->res_1.order == 2){
      amz=0.11706999689340591;
      nlooprun=3;
    }
    as = alphas_ellis_(&mur, &amz, &nlooprun)/pi;
    alpqf = alphas_ellis_(&mures, &amz, &nlooprun)/(4*pi);
    resuPS_0.alphas = as;
    resuPS_0.alphaqf = alpqf;
//  SWITCHING FUNCTION
    double switch0 = 1.;
    double q = std::sqrt(q2), qt = std::sqrt(qt2);
    if(qt >= q*3./4. && q > 0) switch0=std::exp( -std::pow((q*3./4.-qt)/(q/2.),2) ); // GAUSS SWITCH
    if(switch0<0.01) switch0 = 0.;
    randsjacob = randsjacob*switch0;
*/

    return randsjacob;

}
