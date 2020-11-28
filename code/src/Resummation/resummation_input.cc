
#include "resummation_input.h"

#include "read_data.h"
#include "constants.h"


// Reads resummation-specific input and returns a ResummationInfo object
void resummation_input(std::string filename, ResummationInfo& res_1){

  std::ifstream infile(filename.c_str());

  std::cout << "Reading general resummation input" << std::endl;

// The additional parameters in this case are mainly for kinematic cuts

  std::map<std::string, int> input_status;
  std::map<std::string, const int> default_ints;
  std::map<std::string, const double> default_reals;
  std::map<std::string, const std::string> default_strings;

// Do the reading
  std::string line;
  do{

    std::getline(infile, line);

    if(line[0]!='#'){

// Define Initial state (p and pbar currently supported) and PDF set
      set_def<double>(default_reals,"CM_energy",13000.);
      read_data<double>(input_status,default_reals,line,"CM_energy",res_1.CM_energy);
      set_def<int>(default_ints,"lha_flag",1);
      read_data<int>(input_status,default_ints,line,"lha_flag",res_1.lha_flag);
      set_def<std::string>(default_strings,"pdf_setname","MSTW2008nnlo68cl");
      read_data<std::string>(input_status,default_strings,line,"pdf_setname",res_1.pdf_setname);
      set_def<int>(default_ints,"lha_imem",0);
      read_data<int>(input_status,default_ints,line,"lha_imem",res_1.lha_imem);
      set_def<int>(default_ints,"ih1",1);
      read_data<int>(input_status,default_ints,line,"ih1",res_1.ih1);
      set_def<int>(default_ints,"ih2",1);
      read_data<int>(input_status,default_ints,line,"ih2",res_1.ih2);

// Order of calculation
      set_def<int>(default_ints,"order",0);
      read_data<int>(input_status,default_ints,line,"order",res_1.order);
// resummation & CT on or off
      set_def<int>(default_ints,"resum_flag",1);
      read_data<int>(input_status,default_ints,line,"resum_flag",res_1.resum_flag);
      set_def<int>(default_ints,"CT_flag",0);
      read_data<int>(input_status,default_ints,line,"CT_flag",res_1.CT_flag);

// Scales
      set_def<int>(default_ints,"muR_flag",1);
      read_data<int>(input_status,default_ints,line,"muR_flag",res_1.muR_flag);
      set_def<int>(default_ints,"muF_flag",1);
      read_data<int>(input_status,default_ints,line,"muF_flag",res_1.muF_flag);
      set_def<double>(default_reals,"mu_R",1.);
      read_data<double>(input_status,default_reals,line,"mu_R",res_1.mu_R);
      set_def<double>(default_reals,"mu_F",1.);
      read_data<double>(input_status,default_reals,line,"mu_F",res_1.mu_F);
      set_def<double>(default_reals,"mu_S",1.);
      read_data<double>(input_status,default_reals,line,"mu_S",res_1.mu_S);

// For PDF fit
      set_def<double>(default_reals,"mu_min",20.);
      read_data<double>(input_status,default_reals,line,"mu_min",res_1.mu_min);
      set_def<double>(default_reals,"pdf_fussiness",0.01);
      read_data<double>(input_status,default_reals,line,"pdf_fussiness",res_1.pdf_fussiness);
      set_def<double>(default_reals,"en_sec_multiplier",2.);
      read_data<double>(input_status,default_reals,line,"en_sec_multiplier",res_1.en_sec_multiplier);
// Specific PDF fit file (default none)
      set_def<std::string>(default_strings,"pdffit_file","");
      read_data<std::string>(input_status,default_strings,line,"pdffit_file",res_1.pdffit_file);


// Basic parameters for event phase space definition
      set_def<double>(default_reals,"QQ_Min",0.);
      read_data<double>(input_status,default_reals,line,"QQ_Min",res_1.QQ_Min);
      set_def<double>(default_reals,"QQ_Max",500.);
      read_data<double>(input_status,default_reals,line,"QQ_Max",res_1.QQ_Max);
      set_def<double>(default_reals,"QT_Min",0.);
      read_data<double>(input_status,default_reals,line,"QT_Min",res_1.QT_Min);
      set_def<double>(default_reals,"QT_Max",0.);
      read_data<double>(input_status,default_reals,line,"QT_Max",res_1.QT_Max);
      set_def<double>(default_reals,"eta_Min",0.);
      read_data<double>(input_status,default_reals,line,"eta_Min",res_1.eta_Min);
      set_def<double>(default_reals,"eta_Max",0.);
      read_data<double>(input_status,default_reals,line,"eta_Max",res_1.eta_Max);

// NP parameters
      set_def<double>(default_reals,"ggnp",0.);
      read_data<double>(input_status,default_reals,line,"ggnp",res_1.ggnp);
      set_def<double>(default_reals,"gqnp",0.);
      read_data<double>(input_status,default_reals,line,"gqnp",res_1.gqnp);

// For De-Quadrature (intde2.cc)
      set_def<int>(default_ints,"lenaw",8000);
      read_data<int>(input_status,default_ints,line,"lenaw",res_1.lenaw);
      set_def<double>(default_reals,"tiny",1.e-307);
      read_data<double>(input_status,default_reals,line,"tiny",res_1.tiny);
      set_def<double>(default_reals,"de_eps",0.01);
      read_data<double>(input_status,default_reals,line,"de_eps",res_1.de_eps);

// Some basic SM constants
      set_def<int>(default_ints,"Nc",3);
      read_data<int>(input_status,default_ints,line,"Nc",res_1.Nc);
      set_def<int>(default_ints,"Nf",5); // Flavours for fixed flavour scheme
      read_data<int>(input_status,default_ints,line,"Nf",res_1.Nf);
      set_def<double>(default_reals,"alpha_QED",1./137.0359997);
      read_data<double>(input_status,default_reals,line,"alpha_QED",res_1.alpha_QED);
      set_def<double>(default_reals,"gevm2topb",389379660.);
      read_data<double>(input_status,default_reals,line,"gevm2topb",res_1.gevm2topb);
      set_def<double>(default_reals,"M_p",0.9382720813);
      read_data<double>(input_status,default_reals,line,"M_p",res_1.M_p);
      set_def<double>(default_reals,"mTop",1.e10);
      read_data<double>(input_status,default_reals,line,"mTop",res_1.mTop);

    }
  }while(!infile.eof());

  res_1.verbosity = 0; // To be fixed

  if(res_1.fitonly==1){
    res_1.resum_flag=1;
    input_status["resum_flag"]=0;
  }
  if(input_status["pdf_setname"] == 1){
    if(res_1.lha_flag==1){
      if(res_1.order==0){
        set_def<std::string>(default_strings,"pdf_setname","MSTW2008lo68cl");
        res_1.pdf_setname="MSTW2008lo68cl";
      }
      if(res_1.order==1){
        set_def<std::string>(default_strings,"pdf_setname","MSTW2008nlo68cl");
        res_1.pdf_setname="MSTW2008nlo68cl";
      }
      if(res_1.order==2){
        set_def<std::string>(default_strings,"pdf_setname","MSTW2008nnlo68cl");
        res_1.pdf_setname="MSTW2008nnlo68cl";
      }
    }
    else{
      if(res_1.order==0){
        set_def<std::string>(default_strings,"pdf_setname","mstw2008lo");
        res_1.pdf_setname="mstw2008lo";
      }
      if(res_1.order==1){
        set_def<std::string>(default_strings,"pdf_setname","mstw2008nlo");
        res_1.pdf_setname="mstw2008nlo";
      }
      if(res_1.order==2){
        set_def<std::string>(default_strings,"pdf_setname","mstw2008nnlo");
        res_1.pdf_setname="mstw2008nnlo";
      }
    }
  }

  dump_default_parameters(input_status, default_ints, default_reals, default_strings);

  res_1.auto_etalim = false;
  if(input_status["eta_Max"] == 1) res_1.auto_etalim = true;
  res_1.auto_qtlim = false;
  if(input_status["QT_Max"] == 1) res_1.auto_qtlim = true;


// A sanity check for PDF fit inputs
  if(res_1.pdf_fussiness <= 0.){
    std::cout << "Parameter \"pdf_fussiness\" is " << res_1.pdf_fussiness << ", while it is supposed to be > 0. Using default value " << default_reals["pdf_fussiness"] << " instead" << std::endl;
    res_1.pdf_fussiness = default_reals["pdf_fussiness"];
  }
  if(res_1.en_sec_multiplier <= 1.){
    std::cout << "Parameter \"en_sec_multiplier\" is " << res_1.en_sec_multiplier << ", while it is supposed to be > 1. Using default value " << default_reals["en_sec_multiplier"] << " instead" << std::endl;
      res_1.en_sec_multiplier = default_reals["en_sec_multiplier"];
  }

}
