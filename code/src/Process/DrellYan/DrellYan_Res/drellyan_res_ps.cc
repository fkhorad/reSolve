
#include "drellyan_res_ps.h"

#include "phase_space.h"
#include "constants.h"

#include "drellyan_input.h"

#include "resu_preproc.h"
#include "resu_PS.h"

#include <fstream> // TEMP
#include <sstream> // TEMP
#include <iomanip>

//extern std::ifstream debbie;

double drellyan_res_ps(const double x[], drellyan_input* drellyan_in, PSpoint& pp1, resu_PS& resuPS_0){

    double pi = k_constants::pi;

// Dimension check: this routine needs a total of 4, 6 or 8 -- depending on resum and CT flags -- randoms to generate the phase space -- stop if dimension is incorrect
    if(drellyan_in->res_1.CT_flag==1){
      if( !(drellyan_in->ndim == 8) ){
        std::cout << "Something's wrong: ndim does not match phase-space generation routine" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    else if(drellyan_in->res_1.resum_flag==1){
      if( !(drellyan_in->ndim == 6) ){
        std::cout << "Something's wrong: ndim does not match phase-space generation routine" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    else if(drellyan_in->res_1.order == 0){
      if( !(drellyan_in->ndim == 4) ){
        std::cout << "Something's wrong: ndim does not match phase-space generation routine" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    else{
      std::cout << "Something's wrong: order, resummation or CT flags uncorrectly set" << std::endl;
      exit(EXIT_FAILURE);
    }


//Assign q2, qt2, eta, thetaCM, phiCM with randoms from vegas
    double q2, qt2, eta, thetaCM, c_thetaCM, phiCM;
    double randsjacob = 1.;
    double CM_energy = drellyan_in->res_1.CM_energy;
//    //TCRIDGE I THINK THAT WHEN WE DO THE CT WE EXCHANGE THE Q2 AND Y INTEGRALS FOR Z1 (->ALPHA1) and Z2 (->ALPHA2) INTEGRALS, THEREFORE PERHAPS RANDSJACOB FOR THE CT SHOULD NOT THEN HAVE THE PIECES FOR THE Q2 AND Y? SHOULD I TRY REMOVING THESE OR IS IT THAT INSTEAD I SHOULD INCLUDE THE Q^2 FACTOR (AND SOMETHING FOR Y PRESUMABLY) STILL IN THE CT CALCULATION?
    if (drellyan_in->DYnarrowwidthapprox!=1) {
      q2 = std::pow(drellyan_in->res_1.QQ_Min,2) + x[0]*(std::pow(drellyan_in->res_1.QQ_Max,2) - std::pow(drellyan_in->res_1.QQ_Min,2));
      randsjacob = randsjacob*(std::pow(drellyan_in->res_1.QQ_Max,2)-std::pow(drellyan_in->res_1.QQ_Min,2));
    }
    else if ((drellyan_in->DYprocess == 4 || drellyan_in->DYprocess == 5)) {
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
    eta = eta_Min + x[1]*(eta_Max - eta_Min);
    randsjacob = randsjacob*(eta_Max - eta_Min);
//
    if (drellyan_in->res_1.resum_flag > 0 || drellyan_in->res_1.CT_flag > 0) {
      double QT_Min = drellyan_in->res_1.QT_Min, QT_Max = drellyan_in->res_1.QT_Max;
      if(drellyan_in->res_1.auto_qtlim){
        QT_Min = 0.;
        QT_Max = std::sqrt( std::pow( CM_energy*CM_energy + q2, 2) / (4.*CM_energy*CM_energy*std::cosh(eta)*std::cosh(eta)) - q2 );
      }
      qt2 = std::pow(QT_Min,2) + x[4]*(std::pow(QT_Max,2) - std::pow(QT_Min,2));
      randsjacob = randsjacob*(std::pow(QT_Max,2)-std::pow(QT_Min,2));
    }
    else {
      qt2 = 0;
    }
//
    thetaCM = pi*x[2]; //thetaCM between 0 and pi
    c_thetaCM = std::cos(thetaCM);
    randsjacob = randsjacob*pi*std::sin(thetaCM); // the Jacobian is actually for d(costheta)
//
    phiCM = 2*pi*x[3];   //phiCM between 0 and 2pi
    randsjacob = randsjacob*2*pi;
//
// azimuthal angle of the  total final state momentum: not part of the integration
    double phisys = 0.;
    if (drellyan_in->res_1.resum_flag > 0 || drellyan_in->res_1.CT_flag > 0) double phisys = 2*pi*x[5];


// Total momentum
    four_momentum qvec(4);
    double MT = std::sqrt(q2 + qt2);
    qvec[0] = MT*std::cosh(eta); qvec[1] = std::sqrt(qt2)*std::cos(phisys);
    qvec[2] = std::sqrt(qt2)*std::sin(phisys); qvec[3] = MT*std::sinh(eta);


// PS setting: initial state
    //TCRIDGE COMMENT OUT THESE NEXT 10 LINES WHEN SETTING SAME AS MCFM FOR TESTING - OCT20
    four_momentum k1(4), k2(4);
    in_state_w_recoil(qvec, CM_energy, 0., k1, k2);
    pp1.add_mom(k1);
    pp1.add_mom(k2);
// PS setting: final state
    set_PS_twobody(c_thetaCM, phiCM, qvec, 0., 0., pp1);
    pp1.set_products();
// PS setting: reduced version for resummation
    resuPS_0.set(drellyan_in->res_1, q2, eta, qt2);
    
    //TCRIDGE OCT 20 FOR READING IN MOMENTA FOR TESTING                      
    // PSpoint pp1_MCFM;
    // four_momentum k1_MCFM(4),k2_MCFM(4), qvec_MCFM(4), pp1_MCFM_mom1(4), pp1_MCFM_mom2(4);
    // k1_MCFM[2] = 0;
    // k1_MCFM[1] = -3.3314207577702164;
    // k1_MCFM[0] = -3.3314207577702164;
    // k2_MCFM[3] = 0;
    // k2_MCFM[2] = 0;
    // k2_MCFM[1] = 615.70518387118511;
    // k2_MCFM[0] = -615.70518387118511;
    // pp1_MCFM_mom1[3] = 34.756629925026012;
    // pp1_MCFM_mom1[2] = -14.093164897416626;
    // pp1_MCFM_mom1[1] = -132.68365083097132;
    // pp1_MCFM_mom1[0] = 137.88252905428632;
    // // std::cout << "here 1... " << std::endl;  
    // pp1_MCFM_mom2[3] = -34.756629925026012;
    // pp1_MCFM_mom2[2] = 14.093164897416626;
    // pp1_MCFM_mom2[1] = -479.69011228244364;
    // pp1_MCFM_mom2[0] = 481.15407557466904;
    // qvec_MCFM[0] = pp1_MCFM_mom1[0] + pp1_MCFM_mom2[0];                        
    // qvec_MCFM[1] = pp1_MCFM_mom1[1] + pp1_MCFM_mom2[1];                        
    // qvec_MCFM[2] = pp1_MCFM_mom1[2] + pp1_MCFM_mom2[2];                        
    // qvec_MCFM[3] = pp1_MCFM_mom1[3] + pp1_MCFM_mom2[3];   

    // pp1.add_mom(k1_MCFM);                                                     
    // pp1.add_mom(k2_MCFM);                                                     
    // pp1.add_mom(pp1_MCFM_mom1);                                              
    // pp1.add_mom(pp1_MCFM_mom2);                                              
    // pp1.set_products();                                                       
    // resuPS_0.set(drellyan_in->res_1, q2, eta, qt2);  


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
    //std::cout << "randsjacob = " << randsjacob << std::endl;

    return randsjacob;

}
