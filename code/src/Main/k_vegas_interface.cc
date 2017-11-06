
#include "k_vegas_interface.h"

#include <fstream>
#include <cstdio>

#include "InputPars.h"


void k_vegas_call(InputPars& input_1, k_integrand_t integrand, void* userdata) {

// All numbers are per core in parallelisation mode
    int nstart = input_1.nstart; //Number of calls in first iteration
    int nincrease = input_1.nincrease; //Number of call increases in first iteration
    int maxeval = input_1.maxeval; //Approx. max total number of calls across all iterations

// Initialise integral object
    vegas_integral prova(input_1.ndim);

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
      for (int i=0; i < input_1.machine_tag.length(); i++) {
        int x = input_1.machine_tag[i];
        std::cout << input_1.machine_tag[i] << " " << x << std::endl;
        y = y*x;
      }

      seedfordebugging = std::abs(y*(prova.r_iter()+1));

      vegas_random_init(seedfordebugging);
    }
    else if (input_1.seed == -2 && input_1.multi_machine != 1) {
      std::cout << "WARNING: multi_machine must be 1 for seed input -2 so it gets seed from machine tag!" << std::endl;
      std::cout << "->seed taken as 0" << std::endl;
    }
    if(input_1.seed >=0 && input_1.multi_machine == 1)
      std::cout << "WARNING: this seed choice is probably unsafe for parallel running" << std::endl;

//

// Try to recover vegas_integral object from file;
// stick to standard if not found

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


    if(input_1.verbosity>=1){
      std::cout << std::endl << " Vegas integration -- maxeval   = " << maxeval << std::endl
                             << "                   -- nstart    = " << nstart << std::endl
                             << "                   -- nincrease = " << nincrease << std::endl;
    }
    int niter = prova.r_iter(); //Iteration number as read from grid file
    std::cout << "iter read: " << niter << std::endl;
    int ntotcalls = nstart*niter + nincrease*niter*(niter-1.)/2.; //Total number of calls so far carried out
    std::cout << "ntotcalls sofar: " << ntotcalls << std::endl;
    int ncalls = nstart + nincrease*niter;
    while (ntotcalls < maxeval) {

//      double avgi,sd,chi2a;

      std::cout << "ncalls: " << ncalls << std::endl;
      prova.integrate(integrand, ncalls, input_1.verbosity, userdata); //first argument will be integrand function

// Save integrand object to statefile
      if( input_1.multi_machine == 1 ){
        std::stringstream filename_part;
        filename_part << input_1.workdir << "part_" << input_1.machine_tag << "_" << statefile;
        std::ofstream outfile(filename_part.str().c_str(), std::ios::out);
        prova.info_dump(outfile);
        std::stringstream listname;
        listname << input_1.workdir << "tag_list.dat";
        std::ofstream listfile(listname.str().c_str(), std::ios::out | std::ios::app);
        listfile << input_1.machine_tag << "\n";
        outfile.close();
        listfile.close();
      }
      else if(input_1.statefile != ""){
        std::ofstream outfile(filename_full.str().c_str(), std::ios::out);
        prova.info_dump(outfile);
      }

// End iteration
      if( input_1.multi_machine == 1 ) break; // Only 1 iteration at a time and no grid reninement if multi-machine tag is set
      else{
        prova.refine_grid();
// Update # of total calls so far
        ntotcalls += ncalls;
      }


    }

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

// DEBUG
        std::cout << "Pieces\n";
        std::cout << gridname.str() << std::endl;
        curr_grid.info_dump();

        std::cout << "code " << grid_code << std::endl;

        if(grid_code == 0)
          veg_parts.push_back(curr_grid);
        gridfile.close();
        std::remove(gridname.str().c_str());
      }
      listfile.close();
      std::remove(listname.str().c_str());

// Combine the partial runs
      if(veg_parts.size() > 0){
        std::cout << veg_parts[0].r_dim() << std::endl;
        vegas_integral provona(veg_parts[0].r_dim());
        provona.combine_partial_iters(veg_parts.size(), &veg_parts[0]);
        provona.refine_grid();
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
      }
      else{
        std::cout << "No well-defined partial grids found and no final grid produced" << std::endl;
      }

 }
