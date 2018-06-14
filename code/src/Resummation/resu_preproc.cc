
#include "resu_preproc.h"

#include "resummation_input.h"
// For PDFs and fit
#include "mstwpdf.h"
#include "partons_cc.h"
// For inverse Fourier transform initialisation
#include "intde2.h"
//
#include "events_out.h"
#include "constants.h"


namespace global_fitpars{

// PDF fit coefficients for hadron 1 & 2
  pdffit pdfbeam1fita;
  pdffit pdfbeam2fita;
// PDF values in Mellin space along the contour
  pdfmellin pdfbeam1pos;
  pdfmellin pdfbeam2pos;
  pdfmellin pdfbeam1min;
  pdfmellin pdfbeam2min;
// Additional, constant PDF fit parameters (NOT the coefficients ifit and
// ifit2, which are PS-dependent and introduced deeper inside the code)
  fitparams fitparamsa;

}


//void resu_preproc(ResummationInfo* resuminfo, c_mstwpdf* pdf, int process, int DYprocess, int DYnarrowwidthapprox, int resum_flag, int hist_only, double mz, double mw){
void resu_preproc(const event_dumper_info& ev_info, ResummationInfo& resuminfo){

// A small extra
  int Nc = resuminfo.Nc;
  resuminfo.Cf = (Nc*Nc-1.0)/(2.0*Nc);
  resuminfo.Ca = Nc;

  resuminfo.tot_em_charge = 0; // Default value

  MSTW_init(resuminfo);

  if(resuminfo.resum_flag==1){

    std::cout << "Initialisation of resummation coefficients" << std::endl;

// Initialisation for intdeo --> inverse fourier transform
// Dynamic allocation of aw for intdeo
    int lenaw = resuminfo.lenaw;
    double tiny = resuminfo.tiny;
    double eps = resuminfo.de_eps;
    double* aw = new double[lenaw];
    intdeoini(lenaw, tiny, eps, aw);
    resuminfo.aw = aw; // aw just a pointer!
//


// Initialisation for N-independent resummation parameters
    resu_init_0(resuminfo.resuNindep, resuminfo.Nf, resuminfo.Ca, resuminfo.Cf);


// Initialisation for N-dependent resummation parameter + Mellin contour
    inv_mel_init(resuminfo.contoursa, resuminfo.PosBranch, resuminfo.NegBranch, resuminfo.verbosity, resuminfo.Nf, resuminfo.Ca, resuminfo.Cf);

    pdffit_interface(ev_info, resuminfo);

    std::cout << "Initialisation ended" << std::endl << std::endl;

  }

}


void MSTW_init(ResummationInfo& resuminfo){

  // QCD inizialisation: PDFs and alpha_s

  // MSTW08 central value -- ONLY AVAILABLE PDF OPTION FOR NOW, BUILT-IN
  // Open the central PDF grid corresponding to the prefix
    char filename[100];
    if (resuminfo.pdf_flag == 80) {
        const char * prefix = "mstw2008lo";
        sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);
    }
    else if (resuminfo.pdf_flag == 81) {
        const char * prefix = "mstw2008nlo";
        sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);
    }
    else if (resuminfo.pdf_flag == 82) {
        const char * prefix = "mstw2008nnlo";
        sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);
    }
  // Declare a pdf object (with "new" as it'll need to survive out of here)
    std::cout << "PDF grid filename = " << filename << std::endl;
    resuminfo.pdf = new c_mstwpdf(filename);

  // Call the Fortran initialisation routine with alpha_S(Q_0).
  // FR2 and mTop in principle should be taken from input; FR2 = (muF/muR)^2 Q0
  // is essentially a dumb parameter (as long as alphaSQ0 is
  // used along with it)
    double FR2 = 1.,  Q0 = 1., mTop = resuminfo.mTop;
    if( (resuminfo.muF_flag==1 && resuminfo.muR_flag==1) || (resuminfo.muF_flag==1 && resuminfo.muR_flag==1) ){
      FR2 = std::pow(resuminfo.mu_F/resuminfo.mu_R,2);
    }
    else if(resuminfo.muF_flag==0 && resuminfo.muR_flag==1){
      double muRtemp = std::sqrt(resuminfo.QQ_Max * resuminfo.QQ_Min)*resuminfo.mu_R;
      FR2 = std::pow(resuminfo.mu_F/muRtemp,2);
    }
    else if(resuminfo.muF_flag==1 && resuminfo.muR_flag==0){
      double muFtemp = std::sqrt(resuminfo.QQ_Max * resuminfo.QQ_Min)*resuminfo.mu_F;
      FR2 = std::pow(muFtemp/resuminfo.mu_R,2);
    }
    else if(resuminfo.muF_flag==0 && resuminfo.muR_flag==0){
      FR2 = std::pow(resuminfo.mu_F/resuminfo.mu_R,2);
    }

    initalphas_(&resuminfo.pdf->alphaSorder,&FR2,&Q0,&resuminfo.pdf->alphaSQ0,&resuminfo.pdf->mCharm,&resuminfo.pdf->mBottom,&mTop);

  // End MSTW08
}


void pdffit_interface(const event_dumper_info& ev_info, ResummationInfo& resuminfo){


// First check if the automatic fit procedure is called for, or if the user desires for a
// specific pdf_fit file to be used


// Constant for now. Change? Include in constants.h?
  std::string pdf_fit_dir = "pdf_fits/"; // Slash included!

  if(resuminfo.pdffit_file != ""){

// Read in the user-requested ad-hoc PDF fit file
    std::stringstream inputname;
    inputname << pdf_fit_dir << resuminfo.pdffit_file;
    std::cout << "Overriding automatic fit parameters and reading from " << inputname.str() << std::endl;
    resuminfo.muF_flag = 0;
    global_fitpars::fitparamsa.nenergysectors = 1;
    pdffitread_in(inputname.str(), 0, global_fitpars::pdfbeam1fita, global_fitpars::pdfbeam2fita);
  }
  else{


// Get parameters for PDF fit
    global_fitpars::fitparamsa.fitparamscalc(resuminfo.mu_min, resuminfo.en_sec_multiplier, resuminfo.QQ_Min, resuminfo.QQ_Max, resuminfo.mu_F, resuminfo.muF_flag);

    double xtau = (resuminfo.QQ_Max*resuminfo.QQ_Max)/(resuminfo.CM_energy*resuminfo.CM_energy);
    double etam = std::abs(resuminfo.eta_Max);
    int nenergysectors = global_fitpars::fitparamsa.nenergysectors;
    int pdf_flag = resuminfo.pdf_flag;
    int nf = resuminfo.Nf;

    std::cout << "We need " << nenergysectors << " pdf_fit file(s) at mu_F = ";
    for(int ii=0; ii<nenergysectors; ii++)
      std::cout << global_fitpars::fitparamsa.mufgrid[ii] << "  ";
    std::cout << "." << std::endl;


// Before launching the (time-consuming) fit, check if the needed .dat files
// are already around.
    std::vector<int> pdffit_found (nenergysectors, 0);
    std::vector<std::string> pdffit_names(nenergysectors, "");
    int last_file = 0;
    double pdf_fussiness = resuminfo.pdf_fussiness;
    int max_fit_files = k_constants::max_fit_files;
    for(int kk=0; kk < max_fit_files; kk++){
      double xtau_temp, muf_temp;
      int pdf_flag_temp;
      std::string word;
      std::stringstream in_name;
      in_name << pdf_fit_dir << "pdf_fit_" << kk << ".dat";
      std::ifstream in_file(in_name.str().c_str());
      if(in_file.is_open()){
        last_file = kk;
        in_file >> word >> xtau_temp;
        in_file >> word >> muf_temp;
        in_file >> word >> pdf_flag_temp;
        in_file.close();
        if( (std::abs(xtau_temp/xtau-1.) < pdf_fussiness) && (pdf_flag==pdf_flag_temp) ){
          for(int jj=0;jj<nenergysectors;jj++){
            if( std::abs(muf_temp/global_fitpars::fitparamsa.mufgrid[jj]-1.) < pdf_fussiness ){
              pdffit_found[jj] = 1;
              pdffit_names[jj] = in_name.str();
              std::cout << "Acceptable PDF fit found at mu_F = " << muf_temp << std::endl;
              break;
            }
          }
        }
      }
    }

// PDF grid -- for now assumes use of MSTW08
    pdf_tofit = resuminfo.pdf; // Needs to be a global object
// I'm trying to avoid global vars, expecially big ones, so this one will
// only be used by the pdffit routine and deleted immediately after the fit.
    for(int jj=0; jj<nenergysectors; jj++){
      if(pdffit_found[jj]==0){
        double muf2 = global_fitpars::fitparamsa.mufgrid[jj]*global_fitpars::fitparamsa.mufgrid[jj];
        int aa_in = k_constants::aa;
        int jjp = jj+1;
        fiteador_(xtau,muf2,jjp,pdf_flag,aa_in);
        std::stringstream outname;
        if(ev_info.multi_machine == 0)
          outname << pdf_fit_dir << "pdf_fit_" << last_file+1+jj << ".dat";
        else if(ev_info.multi_machine == 1)
          outname << pdf_fit_dir << "pdf_fit_" << ev_info.machine_tag << "_" << last_file+1+jj << ".dat";
        else{
          std::cout << "Wrong multi_machine tag in pdffit_interface" << std::endl;
          exit(EXIT_FAILURE);
        }
        writepdfout_(outname.str().c_str(), xtau, global_fitpars::fitparamsa.mufgrid[jj], jjp, pdf_flag);
        pdffit_names[jj] = outname.str();
      }
    }
// End of PDF-FIT block


// Now read the PDF fit file(s) in
    if(resuminfo.pdffit_file == ""){
      for(int jj=0; jj<nenergysectors; jj++){
        pdffitread_in(pdffit_names[jj], jj, global_fitpars::pdfbeam1fita, global_fitpars::pdfbeam2fita);
        std::cout << "PDF fit file n. " << jj << ": " << pdffit_names[jj] << std::endl;
      }
    }
   }

// Define PDFs in Mellin space along the contour
   pdfmomentsoverallcalc(resuminfo.contoursa, global_fitpars::pdfbeam1pos,global_fitpars::pdfbeam2pos, global_fitpars::pdfbeam1min, global_fitpars::pdfbeam2min, global_fitpars::pdfbeam1fita, global_fitpars::pdfbeam2fita, resuminfo.verbosity, global_fitpars::fitparamsa.nenergysectors);

// That's it

   std::cout << "PDF fit and reading ended" << std::endl;


}
