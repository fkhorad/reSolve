
#include "k_vegas_interface.h"

#include <cstdio>

#include "InputPars.h"
#include "eventsread_in.h"
#include "event_reader.h"


void k_vegas_call(const InputPars& input_1, int ndim, k_integrand_t integrand, void* userdata) {

// nstart and nincrease are PER CORE in parallelisation mode
    int nstart = input_1.nstart; //Number of calls in first iteration
    int nincrease = input_1.nincrease; //Number of call increases in first iteration
// maxeval ignored in parallelisation mode
    int maxeval = input_1.maxeval; //Approx. max total number of calls across all iterations

// Initialise integral object
    vegas_integral prova(ndim);
// Try to recover vegas_integral object from statefile, if present;
    std::string statefile;
    if(input_1.statefile != "")
      statefile = input_1.statefile;
    else
      statefile = "grid.dat";
    std::stringstream filename_full;
    filename_full << input_1.workdir << statefile;
    std::ifstream gridin(filename_full.str().c_str(), std::ios::in);
    if( !(gridin.fail()) ){
      prova.copy_from_file(gridin);
    }
    gridin.close();


// SET THE RANDOM SEED
    if(input_1.seed == -1){
      time_t seed_0;
      time(&seed_0);
      vegas_random_init(seed_0);
    }
    else if(input_1.seed > 0){
      vegas_random_init(input_1.seed);
    }
    else if (input_1.seed == -2 && input_1.multi_machine == 1) //seed generated from tag and iteration number for reproducible different seeds for parallelisation
    {
      int seedfordebugging = 0;
      std::cout << "input_1.machine_tag = " << input_1.machine_tag << std::endl;

      std::string::size_type sz; // alias of size_t
      int y = 1;
      for (unsigned int i=0; i < input_1.machine_tag.length(); i++) {
        int x = input_1.machine_tag[i];
        std::cout << input_1.machine_tag[i] << " " << x << std::endl;
        y = y*x;
      }
//
      seedfordebugging = std::abs(y*(prova.r_iter()+1));
//
      vegas_random_init(seedfordebugging);
    }
    else if (input_1.seed == -2 && input_1.multi_machine != 1) {
      std::cout << "WARNING: multi_machine must be 1 for seed input -2 so it gets seed from machine tag!" << std::endl;
      std::cout << "->seed taken as 0" << std::endl;
    }
    if(input_1.seed >=0 && input_1.multi_machine == 1)
      std::cout << "WARNING: this seed choice is probably unsafe for parallel running" << std::endl;

//


// Integral
    std::cout << std::endl << " Vegas integration -- maxeval   = " << maxeval << std::endl
                             << "                   -- nstart    = " << nstart << std::endl
                             << "                   -- nincrease = " << nincrease << std::endl;
    int niter = prova.r_iter(); //Iteration number as read from grid file
    std::cout << "iter read: " << niter << std::endl;
    int ntotcalls = nstart*niter + nincrease*niter*(niter-1.)/2.; //Total number of calls so far carried out
    std::cout << "ntotcalls sofar: " << ntotcalls << std::endl;
    int ncalls = nstart + nincrease*niter;
//
// START ITERATION
    do{
      std::cout << "ncalls: " << ncalls << std::endl;
      std::cout << "nstart = " << nstart << std::endl;
      std::cout << "nincrease = " << nincrease << std::endl;
      std::cout << "niter = " << niter << std::endl;
// Main
      prova.integrate(integrand, ncalls, input_1.verbosity, userdata);
//      std::cout << "here\n";
// Save integrand object to statefile, if necessary
      if( input_1.multi_machine == 1 ){
        std::stringstream filename_part;
//        std::cout << "there\n";
        filename_part << input_1.workdir << "part_" << input_1.machine_tag << "_" << statefile;
        std::ofstream outfile(filename_part.str().c_str(), std::ios::out);
        prova.info_dump(outfile);
//        std::cout << "tthere\n";
        std::stringstream listname;
        listname << input_1.workdir << "tag_list.dat";
        std::ofstream listfile(listname.str().c_str(), std::ios::out | std::ios::app);
        listfile << input_1.machine_tag << "\n";
        outfile.close();
        listfile.close();
//        std::cout << "ttthere\n";
      }
      else if(input_1.statefile != ""){
        std::ofstream outfile(filename_full.str().c_str(), std::ios::out);
        prova.info_dump(outfile);
      }
// Do necessary updates if single-core mode
      if( input_1.multi_machine == 0 ){
        prova.refine_grid();
        ntotcalls += ncalls;
        niter = prova.r_iter();
        ncalls = nstart + nincrease*niter;
      }
//
    } while(ntotcalls < maxeval && input_1.multi_machine == 0);
// End iteration


    std::stringstream main_outfile;
    if(input_1.multi_machine==0){
      main_outfile << input_1.workdir << "reSolve_main_out.dat";
    }
    else if(input_1.multi_machine==1){
      main_outfile << input_1.workdir << "reSolve_main_out_" << input_1.machine_tag << ".dat";
    }
    std::ofstream main_file(main_outfile.str().c_str(), std::ios::out);
    if(main_file.is_open()){
      prova.iteration_output(main_file);
    }
    main_file.close();

}


// Combine several k_vegas integration grids from files (for our KISS-version parallelization)
void k_vegas_combiner(InputPars& input_1){

// Read list of partial runs in this iter
      std::stringstream listname;
      listname << input_1.workdir << "tag_list.dat";
      std::ifstream listfile(listname.str().c_str(), std::ios::in);
      int iter = 0;
      int ndim = 1;
      int new_nevent;
      std::vector<int> old_nevent(0);

// Read in the state files from each partial run
      std::string curr_tag;
      std::vector<vegas_integral> veg_parts;
      while(listfile >> curr_tag){
        std::stringstream gridname;
        if(input_1.statefile != "")
          gridname << input_1.workdir << "part_" << curr_tag << "_" << input_1.statefile;
        else
          gridname << input_1.workdir << "part_" << curr_tag << "_grid.dat";
        std::ifstream gridfile(gridname.str().c_str(), std::ios::in);
        vegas_integral curr_grid;
        int grid_code = curr_grid.copy_from_file(gridfile);

// Print some info
        std::cout << "Pieces\n";
        std::cout << gridname.str() << std::endl;
//        curr_grid.info_dump();
        std::cout << "code " << grid_code << std::endl;
//

        if(grid_code == 0) veg_parts.push_back(curr_grid);
        old_nevent.push_back(curr_grid.r_neval());
        gridfile.close();
        std::remove(gridname.str().c_str());
      }
      listfile.close();

//
// Initialise combined integral object
      vegas_integral provona(1);
// Combine the partial runs
      if(veg_parts.size() > 0){
        provona.combine_partial_iters(veg_parts.size(), &veg_parts[0]);
        provona.refine_grid();
        iter = provona.r_iter();
        new_nevent = provona.r_neval();

// Output
        std::stringstream gridout;
        if(input_1.statefile != "")
          gridout << input_1.workdir << input_1.statefile;
        else
          gridout << input_1.workdir << "grid.dat";
        std::remove(gridout.str().c_str());
        std::ofstream grid_out(gridout.str().c_str());
        provona.info_dump(grid_out);
        provona.iteration_output();
//
        std::stringstream main_outfile;
        main_outfile << input_1.workdir << "reSolve_main_out_" << input_1.machine_tag << ".dat";
        std::cout << "MAIN OUT FILE:  " << main_outfile.str() << std::endl;
        std::ofstream main_file(main_outfile.str().c_str(), std::ios_base::app);
        if(main_file.is_open()){
          provona.iteration_output(main_file);
        }
        main_file.close();

      }
      else{
        std::cout << "No well-defined partial grids found and no final grid produced" << std::endl;
      }

// Finally, combine event files

      int event_type = input_1.save_events;
      if(iter>0 && (event_type==1 || event_type==2) ){
        std::ifstream listfile_2(listname.str().c_str(), std::ios::in);
        std::string filename = get_event_filename(event_type, input_1.workdir, "", iter);
        int list_counter = 0;
        while(listfile_2 >> curr_tag){
          std::string part_filename = get_event_filename(event_type, input_1.workdir, curr_tag, iter);
          combine_reweight(part_filename, filename, event_type, old_nevent[list_counter], new_nevent);
          list_counter++;
        }
        listfile_2.close();
      }
      std::remove(listname.str().c_str());


}



void combine_reweight(std::string part_filename, std::string filename, int event_type, int old_nevent, int new_nevent){

// Define streams
   std::ifstream part_event_file(part_filename.c_str());
   std::ofstream event_file(filename.c_str(), std::ios::out | std::ios::app);
//
  std::string line;
  std::vector<std::string> event;
//
// Parse in-file
  do{
    std::getline(part_event_file, line);
//
// READ, REWEIGHT
    if(event_type==1){
      if(line != ""){
        reweight_easy_event(line, part_event_file, event, old_nevent, new_nevent);
      }
    }
    else if(event_type==2){
      if(line == "<event>"){
        reweight_lhe_event(line, part_event_file, event, old_nevent, new_nevent);
      }
    }
// WRITE
    for(unsigned int ii=0; ii<event.size(); ii++){
      event_file << event[ii] << std::endl;
    }

  }while(!part_event_file.eof());

// Close streams, remove redundant file
   part_event_file.close();
   event_file.close();
   std::remove(part_filename.c_str());
}
//
void reweight_easy_event(std::string line, std::ifstream& part_event_file, std::vector<std::string>& event, int old_nevent, int new_nevent){

  event.resize(0);
  event.push_back(line);
  std::string new_line;
  do{
    std::getline(part_event_file, new_line);
    event.push_back(new_line);
  } while(new_line != "");
//
  std::vector<double> weight_line;
  line_to_vec(event[event.size()-2], weight_line);
  double f_value = weight_line[0];
  double weight = weight_line[1]*old_nevent/new_nevent;
//
  std::stringstream replacement;
  replacement << f_value << "  " << weight;
  event[event.size()-2] = replacement.str();

}
void reweight_lhe_event(std::string line, std::ifstream& part_event_file, std::vector<std::string>& event, int old_nevent, int new_nevent){

  event.resize(0);
  event.push_back(line);
  std::string new_line;
  do{
    std::getline(part_event_file, new_line);
    event.push_back(new_line);
  } while(new_line != "</event>");
//
  std::vector<double> weight_line;
  line_to_vec(event[1], weight_line);
  double f_weight = (old_nevent*weight_line[2])/new_nevent;
  int n_particles = (int) weight_line[0];
  double mu_R = weight_line[3];
  double alpha_QED = weight_line[4];
  double alpha_s = weight_line[5];
//
  std::stringstream replacement;
  replacement << " " << n_particles << "   " << 0 << "  " << f_weight << "  " << mu_R << "  " << alpha_QED << "  " << alpha_s;
  event[1] = replacement.str();

}



void debugger(std::string inputname, int eventsets, int ndim, k_integrand_t integrand, void* userdata){

      std::vector<std::vector<double> > xvec;
      eventsread_in (0, inputname, xvec, eventsets, ndim);
      double res;
      for(int k = 0; k<eventsets; k++) {
         double xx[ndim];
         for(int jj=0; jj<ndim; jj++) xx[jj] = xvec[k][jj];
         integrand(ndim, xx, res, userdata, 1., 1);
      }

}
