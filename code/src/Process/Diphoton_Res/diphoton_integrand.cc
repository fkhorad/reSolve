//Integrands for diphoton case

#include "diphoton_integrand.h"

extern int k_count00, k_count0;
int k_count00 = 0;


int diphoton_integrand(int ndim, const double x[], double& f, void *userdata,
                       double weight, int iter) {
    // NOTE process only passed here to allow process = -1 to do born only for testing
    DiphRes* data = (DiphRes*) userdata;
    ResummationInfo* resuminfo = data->res_1;
    diphoton_input* diph_in = data->diph_in;

    double pi = k_constants::pi;

// DEBUG
    k_count00++;


// Phase space
    PSpoint pp1;
    double randsjacob;
    PSdep_variables resuvars;

    double jacob = diphoton_ps(x, diph_in->CM_energy, resuminfo,
			       randsjacob, pp1, resuvars);

// DEBUG
   // std::cout << "Diphoton Momenta" << std::endl;
   // LorPrint(pp1.mom(2));
   // LorPrint(pp1.mom(3));
   // std::cout << std::endl;


// CUTS:
    bool cutcheck = false; //i.e not cut
    cutcheck = diph_cuts(diph_in, resuvars.q2, resuvars.qt2, resuvars.eta, resuvars.mures2, &pp1);


    double tempf;
    if(cutcheck==true){
      tempf = 0.;

// DEBUG
    // std::cout << "Main counters  " << k_count0 << "  " << k_count00 << std::endl;
    // std::cout << "xborn and co. == 0. here" << std::endl << std::endl << std::endl;

    }else{

// Calculate the process-dependent stuff for resummation
      std::vector<std::vector<double> > sigmaij;
      sigmaijdiphotoncalc(diph_in, pp1, jacob, sigmaij, resuvars.alphas);
      for (int j = 0; j<11; j++) {
        for (int k = 0; k < 11; k++) {
          resuvars.sigmaij[j][k] = sigmaij[j][k];
        }
      }

      double s = 2.*pp1.ss(0,1);
      double t = -2.*pp1.ss(0,2);
      double costheta = 1. + 2.*t/s;

      double H1q = 2*H1qdiphoton_DY(costheta, diph_in->verbosity);
      double H1g = 0.0;
      double H2q[5];
      for (int j = 0; j < 5; j++) {
        H2q[j] = 4*H2qdiphotoncalc_DY(&pp1, costheta, j+1,
                 resuminfo->Cf, diph_in->Nf, diph_in->Nc, diph_in->verbosity);
      }
      double H2g = 0.0;
      resuvars.H1q = std::complex<double>(H1q,0.0);
      resuvars.H1g = std::complex<double>(H1g,0.0);
      resuvars.H2g = std::complex<double>(H2g,0.0);
      for (int i = 0; i<5; i++) {
        resuvars.H2q[i] = H2q[i];
      }

// Call main resummation function

      tempf = resummed(resuvars, resuminfo);
    }

    f = tempf*randsjacob;
//
    if (diph_in->verbosity >= 10) {
      std::cout << "resummed = " << tempf << std::endl;
      std::cout << "f = " << f << std::endl;
    }

// Print event to file

    if(diph_in->save_events==1 || diph_in->save_events==2){
      std::string filename;
      if(diph_in->integrator_flag == 0 || diph_in->integrator_flag == 1){
        std::stringstream filename_stream;
        std::string filler = "";
        if(diph_in->machine_tag != "") filler = "_";
        if(diph_in->save_events == 1)
          filename_stream << diph_in->workdir << "events_" << iter << filler << diph_in->machine_tag << ".dat";
        if(diph_in->save_events == 2)
          filename_stream << diph_in->workdir << "events_lhe_" << iter << filler << diph_in->machine_tag << ".lhe";
        filename = filename_stream.str();
      }
      if(diph_in->integrator_flag == 2){
        filename = get_pids_filename(0, diph_in->workdir, diph_in->machine_tag, iter, diph_in->save_events);
      }
      if(diph_in->save_events==1){
        dumper_easy(filename, f, weight, &pp1, ndim, x);
      }
      else if(diph_in->save_events==2){
        int n_part = 4;
        double mass[] = {diph_in->M_p, diph_in->M_p, 0., 0.};
        int PDG_num[] = {resuminfo->ih1*2212,resuminfo->ih2*2212,22,22};
        int inout[] = {-1, -1, 1, 1};
        int mother1[] = {0, 0, 1, 1};
        int mother2[] = {0, 0, 2, 2};
        dumper_lhe(filename, f*weight, n_part, mass, PDG_num, inout, mother1, mother2,
                   std::sqrt(resuvars.mur2), diph_in->alpha_QED, pi*resuvars.alphas, &pp1);
      }
    }

    return 0;
}



int diphoton_integrand_cuba(const int* ndim, const double x[], const int *ncomp,
                            double f[], void *userdata,
                            const int* unused1, const int* unused2, const double* weight, const int* iter){
  
  double fres;
  int outcode;
  
  outcode = diphoton_integrand(*ndim, x, fres, userdata, *weight, *iter);

  if(*ncomp > 1) outcode = 1;
  
  f[0] = fres;
  
  return outcode;
  
}
