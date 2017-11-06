
#include "resu_preproc.h"

// For pdffit
#include "mstwpdf.h"
#include "partons_cc.h"
// For inverse Fourier transform initialisation
#include "intde2.h"


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
// ifit2, which are PS-dependent and introduced deeper into the code)
  fitparams fitparamsa;
  
}


void resu_preproc(ResummationInfo* resuminfo, c_mstwpdf* pdf, int process){

    std::cout << "Initialisation of resummation coefficients" << std::endl;
  
// Copy PDF object into resuminfo
    resuminfo->pdf = pdf;

// Dynamic allocation of aw for intdeo
    int lenaw = resuminfo->lenaw;
    double* aw = new double[lenaw];
    resuminfo->aw = aw;


// Initialisation for intdeo --> inverse fourier transform
    double tiny = resuminfo->tiny;
    double eps = resuminfo->de_eps;
    intdeoini(lenaw, tiny, eps, aw);


// Initialisation for N-independent resummation parameters
    resu_init_0(resuminfo->resuNindep, resuminfo->Nf,
                resuminfo->Ca, resuminfo->Cf);


// Initialisation for N-dependent resummation parameter + Mellin contour
    inv_mel_init(resuminfo->contoursa, resuminfo->PosBranch, resuminfo->NegBranch,
                 resuminfo->verbosity, resuminfo->Nf, resuminfo->Ca, resuminfo->Cf);


    resuminfo->process = process; //NOTE ONLY SET TO ALLOW BORN ONLY OPTION FOR PROCESS = -1 FOR TESTING

// PDF FIT!!
    if(resuminfo->process!=-1) {
	pdffit_interface(resuminfo);
    }

    std::cout << "Initialisation ended" << std::endl << std::endl;

}



void pdffit_interface(ResummationInfo* resuminfo){
  

// First check if the automatic fit procedure is called for, or if the user desires for a
// specific pdf_fit file to be used

// Constant for now. Change? Include in constants.h?
   std::string pdf_fit_dir = "pdf_fits/"; // Slash included!

   if(resuminfo->pdffit_file != ""){

// Read in the user-requested ad-hoc PDF fit file
     std::stringstream inputname;
     inputname << pdf_fit_dir << resuminfo->pdffit_file;
     std::cout << "Overriding automatic fit parameters and reading from " <<
       inputname.str() << std::endl;
     resuminfo->muF_flag = 0;
     global_fitpars::fitparamsa.nenergysectors = 1;
     pdffitread_in(inputname.str(), 0,
                   global_fitpars::pdfbeam1fita, global_fitpars::pdfbeam2fita); 
   }
   else{

// Get parameters for PDF fit
      global_fitpars::fitparamsa.fitparamscalc(resuminfo->mu_min,
               resuminfo->en_sec_multiplier, resuminfo->QQ_Min, resuminfo->QQ_Max,
               resuminfo->mu_F, resuminfo->muF_flag);
      double xtau = (resuminfo->QQ_Max*resuminfo->QQ_Max)/
        (resuminfo->CM_energy*resuminfo->CM_energy);
      int nenergysectors = global_fitpars::fitparamsa.nenergysectors;
      int pdf_flag = resuminfo->pdf_flag;

      std::cout << "We need " << nenergysectors << " pdf_fit file(s) at mu_F = ";
      for(int ii=0; ii<nenergysectors; ii++)
        std::cout << global_fitpars::fitparamsa.mufgrid[ii] << "  ";
      std::cout << "." << std::endl;


// Before doing the fit (time-consuming), check if the needed .dat files
// are already around.
	  std::vector<int> pdffit_found (nenergysectors, 0);
	  std::vector<std::string> pdffit_names(nenergysectors, "");
	  int last_file = 0;
	  double pdf_fussiness = resuminfo->pdf_fussiness;
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
		  if( (std::abs(xtau_temp/xtau-1.) < pdf_fussiness) && (pdf_flag==pdf_flag) )
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
      
// PDF grid -- for now assumes use of MSTW08
// Declare a pdf object
      char filename[100];
      if (resuminfo->pdf_flag == 80) {
         const char * prefix = "mstw2008lo";
         sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);
       }
       else if (resuminfo->pdf_flag == 81) {
         const char * prefix = "mstw2008nlo";
         sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);
       }
      else if (resuminfo->pdf_flag == 82) {
         const char * prefix = "mstw2008nnlo";
         sprintf(filename,"Grids/%s.%2.2d.dat",prefix,0);
      }
      // const char * prefix = resuminfo->prefix.c_str();
      // sprintf(filename,"%s.%2.2d.dat",prefix,0);
	  pdf_tofit = new c_mstwpdf(filename); // Needs to be a global object
// I'm trying to avoid global vars, expecially big ones, so this one will
// only be used by the pdffit routine and deleted immediately after the fit.
	  for(int jj=0; jj<nenergysectors; jj++){
	      if(pdffit_found[jj]==0){ 
		  double muf2 = global_fitpars::fitparamsa.mufgrid[jj]*global_fitpars::fitparamsa.mufgrid[jj];
		  fiteador_(xtau,muf2,jj,pdf_flag);
		  std::stringstream outname;
		  if(resuminfo->multi_machine == 0)
		      outname << pdf_fit_dir << "pdf_fit_" << last_file+1+jj << ".dat";
		  else if(resuminfo->multi_machine == 1)
		      outname << pdf_fit_dir << "pdf_fit_" << resuminfo->machine_tag << "_" << last_file+1+jj << ".dat";
		  else{
		      std::cout << "Wrong multi_machine tag in pdffit_interface" << std::endl;
		      exit(EXIT_FAILURE);
		  }
// TEMP
		  std::cout << "PDF FIT NAME: " << outname.str() << std::endl;
		  writepdfout_(outname.str().c_str(), xtau, global_fitpars::fitparamsa.mufgrid[jj],
			       jj, pdf_flag);
		  pdffit_names[jj] = outname.str();
	      }
	  }
	  delete pdf_tofit;
// End of PDF-FIT block


// Now read the PDF fit file(s) in
      if(resuminfo->pdffit_file == ""){
        for(int jj=0; jj<nenergysectors; jj++){
          pdffitread_in(pdffit_names[jj], jj,
                      global_fitpars::pdfbeam1fita, global_fitpars::pdfbeam2fita);
        }
      }
   }

// Define PDFs in Mellin space along the contour
    pdfmomentsoverallcalc(resuminfo->contoursa, global_fitpars::pdfbeam1pos,
      global_fitpars::pdfbeam2pos, global_fitpars::pdfbeam1min, global_fitpars::pdfbeam2min,
      global_fitpars::pdfbeam1fita, global_fitpars::pdfbeam2fita, resuminfo->verbosity,
      global_fitpars::fitparamsa.nenergysectors);

// That's it

    std::cout << "PDF fit and reading ended" << std::endl;

  
}
