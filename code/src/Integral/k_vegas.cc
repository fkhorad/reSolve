#include "k_vegas.h"

// Integrate on an unit hypercube!

int vegas_integral::integrate(k_integrand_t func, int neval, int verbose, void *userdata){

// Initialization

  int return_code;
  // std::cout << "niter = " << niter << std::endl;
  int niter = last_iter+1;


  if(verbose>1 || (verbose==1 && niter ==1))
    std::cout << "\n" << "****************************************************" << "\n"
      <<" Input parameters for Vegas:  ndim  = "<< int_dim <<"     g_dim = " << grid_dim << "\n"
      <<"                              alpha = " << alpha << "\n\n";

  double ti,ti2,tsi,sdi;
  ti = 0.; ti2 = 0.; // Single iteration total integral estimate and its square
  tsi = 0.;          // Single iteration total squared integral
  sdi = 0.;          // Single iteration variance estimate

  double weight; // Weight to associate to a function evaluation during a single iteration

  for(int ii=0; ii<int_dim; ii++){
    for(int jj=0; jj<grid_dim; jj++){
      dd[ii][jj] = 0.;
      di[ii][jj] = 0.;
    }
  }

  for(int nn=0; nn<neval; nn++){

    double rnd[int_dim], xx[int_dim], ff, ff2;
    int ia[int_dim];
    double x0,xn;
    weight = 1./neval;

    std::clock_t start;
    double duration;
    start = std::clock();

// Generate the point
    vegas_random(int_dim,rnd);                  // RANDOMS CALL
//    randa_(&int_dim,rnd);
    for(int ii=0; ii<int_dim; ii++){
      xn = (grid_dim)*(1.-rnd[ii]);
      ia[ii] = (int) xn;
      if(ia[ii]==0){
        x0 = grid[ii][ia[ii]];
        xx[ii] = (xn-ia[ii])*x0;
      }
      else{
        x0 = grid[ii][ia[ii]]-grid[ii][ia[ii]-1];
        xx[ii] = grid[ii][ia[ii]-1] + (xn-ia[ii])*x0;
      }
      weight = weight*grid_dim*x0;
    }

// Function evaluation + jacobian, standard dev., chi2...
//
    func(int_dim,xx,ff,userdata,weight,niter);          // FUNCTION CALL


    if (ff != ff) {
      std::cout << "WARNING - nan attained for integrand evaluation in k_vegas.cc for point:" << std::endl;
      std::cout << "nn(evaluation number) = " << nn << std::endl;
      std::cout << "ff(integrand evaluation) = " << ff << std::endl;
      std::cout << "randoms for this point: " << xx[0] << " " << xx[1] << " " << xx[2] << " " << xx[3] << " " << xx[4] << std::endl;
      std::cout << "weight = " << weight << std::endl;
      std::cout << "iteration number: " << niter << std::endl;
      ff = 0.0; //set ff = 0.0 so avoid nan infecting further points and overall answer
      std::cout << "ff taken as: " << ff << std::endl;
      exit(EXIT_FAILURE);
    }

    ff = ff*weight;
    ff2 = ff*ff;         // Functions & squares, always FULLY WEIGHTED
//
    ti = ti + ff;
    tsi = tsi + ff2;


// Data that will be used to refine the grid
    for(int ii=0; ii<int_dim; ii++){
      di[ii][ia[ii]] = di[ii][ia[ii]] + ff;
      dd[ii][ia[ii]] = dd[ii][ia[ii]] + ff2;
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

  }

  ti2 = ti*ti;
  sdi = tsi*neval - ti2;
  sdi = sdi/(neval-1.);

// Final results for this iteration

  double weight2;

  if(sdi==0.)
    weight2=1.e15;
  else
    weight2=ti2/sdi;  // (Relative uncertainty estimate on this iter)^(-2) --
                      // mainly used for averaging over multiple iterations

// TEMP
  std::cout << "PRE\n";
  std::cout << "si: " << si << "\n";
  std::cout << "si2: " << si2 << "\n";
  std::cout << "swgt: " << swgt << "\n";
  std::cout << "schi: " << schi << "\n";


// Store previous values for cross-checks
  pre_si = si;
  pre_si2 = si2;
  pre_swgt = swgt;
  pre_schi = schi;

  if(weight2<1.e14){
    si = si + ti*weight2;
    si2 = si2 + ti2;
    swgt = swgt + weight2;
    schi = schi + ti2*weight2;
//    avgi = si/swgt;
//    sd = swgt*niter/si2;
//    chi2a = sd*(schi/swgt-avgi*avgi)/(niter-.999);
//    sd = std::sqrt(1./sd);

    return_code = 0;
  }
  else{

    if(verbose >= 1) std::cout << " WARNING: zero error estimate!" << "\n\n";

    si = ti;
    si2 = ti*ti;
//    avgi = ti;
//    sd = 0.;
//    chi2a = -1.;

    return_code = -1;
  }

// TEMP
 std::cout << "POST\n";
 std::cout << "si: " << si << "\n";
 std::cout << "si2: " << si2 << "\n";
 std::cout << "swgt: " << swgt << "\n";
 std::cout << "schi: " << schi << "\n";


  last_neval = neval;
  last_iter = niter;
  last_ti = ti;
  last_tsi = tsi;
  last_sdi = sdi;

// Info dumping
  if(verbose  > 1) this->info_dump();

  this->iteration_output();

  return return_code;

}


void vegas_integral::iteration_output(std::ostream & out1){

  double ti = last_ti, sdi = last_sdi, avgi, sd, chi2a;
  double ti2 = ti*ti;
  int niter = last_iter;
  int neval = last_neval;
  double weight2 = ti2/sdi;

  if(weight2<1.e14){
    avgi = si/swgt;
    sd = swgt*niter/si2;
    chi2a = sd*(schi/swgt-avgi*avgi)/(niter-.999);
    sd = std::sqrt(1./sd);
  }
  else{
    out1 << " WARNING: zero error estimate in iteration " << niter << "\n\n";
    avgi = ti;
    sd = 0.;
    chi2a = -1.;
  }

  std::string spacing1 = "", spacing2 = "";
  if(niter<100) spacing1 = " ";
  if(niter<10) spacing1 = "  ";
  if(niter<1000000) spacing2 = " ";
  if(niter<100000) spacing2 = "  ";
  if(niter<10000) spacing2 = "   ";
  if(niter<1000) spacing2 = "    ";
  if(niter<100) spacing2 = "     ";
  if(niter<10) spacing2 = "      ";
  out1 << " Integration by Vegas" << "\n ----\n"
              << " Iteration no. " << niter << spacing1 << "     integral = " << ti << "\n"
              << "  neval      = " << neval << spacing2 <<      "std dev  = " << std::sqrt(sdi) << "\n"
              << " Accumulated results:" <<                   "  integral = " << avgi << "\n"
              << "               " <<                   "        std dev  = " << sd << "\n"
              << "               " <<                          " chi**2 per it'n = " << chi2a << "\n"
              << "\n";


}

void vegas_integral::refine_grid(){


  for(int jj=0; jj<int_dim; jj++){

    double x0, xn, rc, dr, dt[int_dim], rr[grid_dim], xin[grid_dim];

    // std::cout << "int_dim = " << int_dim << std::endl;
    // std::cout << "grid_dim = " << grid_dim << std::endl;

    x0 = dd[jj][0];
    xn = dd[jj][1];
    dd[jj][0] = (x0 + xn)/2.;
    dt[jj] = dd[jj][0];
    for(int ii=1; ii<grid_dim-1; ii++){
      dd[jj][ii] = x0 + xn;
      x0 = xn;
      xn = dd[jj][ii+1];
      dd[jj][ii] = (dd[jj][ii]+xn)/3.;
      dt[jj] = dt[jj] + dd[jj][ii];
    }
    dd[jj][grid_dim-1] = (x0+xn)/2.;
    dt[jj] = dt[jj] + dd[jj][grid_dim-1];

    rc = 0.;
    for(int ii=0; ii<grid_dim; ii++){
      rr[ii] = 0.;
      if(dd[jj][ii]>0.){
        x0 = dt[jj]/dd[jj][ii];
        rr[ii] = std::pow((x0-1.)/x0/std::log(x0), alpha);
      }
      rc = rc + rr[ii];
    }
    rc = rc/grid_dim;

/*
            k=0
            xn=0.
            dr=xn
            i=k
 25         k=k+1
            dr=dr+r(k)
            xo=xn
            xn=xi(k,j)
 26         if(rc.gt.dr)go to 25
            i=i+1
            dr=dr-rc
            xin(i)=xn-(xn-xo)*dr/r(k)
            if(i.lt.ndm)go to 26
*/
    xn = 0.; dr = 0.;
    int kk = -1, ii = -1;
   tag1: kk++;
    dr = dr + rr[kk];
    x0 = xn;
    xn = grid[jj][kk];
   tag2: if(rc>dr) goto tag1;
    ii++;
    dr=dr-rc;
    xin[ii]= xn - (xn-x0)*dr/rr[kk];


    if(ii < grid_dim-2) goto tag2;
/*
      do{
      kk++;
      dr = dr + rr[kk];
      x0 = xn;
      xn = grid[jj][kk];
      do{
        if(rc > dr) break;
        ii++;
        dr=dr-rc;
        xin[ii]= xn - (xn-x0)*dr/rr[kk];
      } while(ii < grid_dim-1);
    } while(rc > dr);
*/


    for(int ii=0; ii<grid_dim-1; ii++){
      grid[jj][ii] = xin[ii];
      // std::cout << "grid[" << jj << "][" << ii << "]=" << grid[jj][ii] << std::endl;
    }
    grid[jj][grid_dim-1] = 1.;

  }

}

void vegas_integral::combine_partial_iters(int nparts, vegas_integral* iter_parts){

 this->reset_resize(iter_parts[0].int_dim);
 this->set_alpha(iter_parts[0].alpha);

// Copy common info from first object in grid: grid+cumulatives
 for(int jj=0; jj<int_dim; jj++){
   for(int kk=0; kk<grid_dim; kk++){
     di[jj][kk] = 0.;
     dd[jj][kk] = 0.;
     grid[jj][kk] = iter_parts[0].grid[jj][kk];
   }
 }
 last_iter = iter_parts[0].last_iter;
//
 si = iter_parts[0].pre_si;
 si2 = iter_parts[0].pre_si2;
 swgt = iter_parts[0].pre_swgt;
 schi = iter_parts[0].pre_schi;

// TEMP
 std::cout << "PRE\n";
 std::cout << "si: " << si << "\n";
 std::cout << "si2: " << si2 << "\n";
 std::cout << "swgt: " << swgt << "\n";
 std::cout << "schi: " << schi << "\n";

 double nti = 0., ntsi = 0.;

 for(int ii=0; ii<nparts; ii++){

   int curr_neval = iter_parts[ii].last_neval;
   last_neval = last_neval + curr_neval;
   nti = nti + curr_neval*iter_parts[ii].last_ti;
   ntsi = ntsi + curr_neval*iter_parts[ii].last_tsi;
 }
 last_ti = nti/last_neval;
 last_tsi = ntsi/last_neval;
//
 for(int ii=0; ii<nparts; ii++){
   for(int jj=0; jj<int_dim; jj++){
     for(int kk=0; kk<grid_dim; kk++){
       double dcurr_neval = (double) iter_parts[ii].last_neval;
       di[jj][kk] = di[jj][kk] + (dcurr_neval/last_neval)*iter_parts[ii].di[jj][kk];
       dd[jj][kk] = dd[jj][kk] + (dcurr_neval/last_neval)*(dcurr_neval/last_neval)*iter_parts[ii].dd[jj][kk];
     }
   }

 }
  double ti2 = last_ti*last_ti;
  last_sdi = last_tsi*last_neval - ti2;
  last_sdi = last_sdi/(last_neval-1.);

  double weight2;
  if(last_sdi==0.)
    weight2=1.e15;
  else
    weight2=ti2/last_sdi;  // (Relative uncertainty estimate on this iter)^(-2) -- mainly used for averaging over multiple iterations

  if(weight2<1.e14){
    si = si + last_ti*weight2;
    si2 = si2 + ti2;
    swgt = swgt + weight2;
    schi = schi + ti2*weight2;
  }
  else{
    si = last_ti;
    si2 = last_ti*last_ti;
    swgt = 0.;
    schi = 0.;
  }

// TEMP
 std::cout << "POST\n";
 std::cout << "si: " << si << "\n";
 std::cout << "si2: " << si2 << "\n";
 std::cout << "swgt: " << swgt << "\n";
 std::cout << "schi: " << schi << "\n";

}



// Constructor

vegas_integral::vegas_integral(){

 int_dim = 1;
 alpha = def_alpha;
 grid_dim = def_grid_dim;

 this->reset();

}
vegas_integral::vegas_integral(int dim){

 int_dim = dim;
 alpha = def_alpha;
 grid_dim = def_grid_dim;

 this->reset();

}

void vegas_integral::reset(){

 this->reset_cumulatives();
 alpha = def_alpha;

// Initialize an evenly-spaced grid

 double step = 1./grid_dim;
 grid.resize(int_dim);
 dd.resize(int_dim);
 di.resize(int_dim);
 for(int i=0; i<int_dim; i++){
   grid.at(i).resize(grid_dim);
   dd.at(i).resize(grid_dim);
   di.at(i).resize(grid_dim);
  for(int j=0; j<grid_dim-1; j++){
     grid.at(i).at(j) = (j+1)*step;
   }
   grid.at(i).at(grid_dim-1) = 1.;
 }

}
void vegas_integral::reset_resize(int int_dim_new, int grid_dim_new){

 int_dim = int_dim_new;
 grid_dim = grid_dim_new;

 this->reset();

}
void vegas_integral::reset_cumulatives(){

 si = 0.; si2 = 0.; swgt = 0.; schi = 0.;
 pre_si = 0.; pre_si2 = 0.; pre_swgt = 0.; pre_schi = 0.;
 last_neval = 0; last_iter = 0;
 last_ti = 0.; last_tsi = 0.; last_sdi = 0.;

}
void vegas_integral::set_alpha(double new_alpha){
  alpha =  new_alpha;
}
void vegas_integral::copy(vegas_integral& vi1){
  si = vi1.si; si2 = vi1.si2; swgt = vi1.swgt; schi = vi1.schi;
  int_dim = vi1.int_dim; grid_dim = vi1.grid_dim; alpha = vi1.alpha;
  grid = vi1.grid;
}

int vegas_integral::r_dim(){
  return grid.size();
}
int vegas_integral::r_gdim(){
  return grid[0].size();
}
double vegas_integral::r_grid(int i,int j){
  return grid[i][j];
}
double vegas_integral::r_si(){
  return si;
}
double vegas_integral::r_si2(){
  return si2;
}
double vegas_integral::r_swgt(){
  return swgt;
}
double vegas_integral::r_schi(){
  return schi;
}
int vegas_integral::r_neval(){
  return last_neval;
}
int vegas_integral::r_iter(){
  return last_iter;
}
//
void vegas_integral::info_dump(std::ostream & out1){

  out1 << "Basics:\n";
  out1 << "  alpha: " << alpha << "  int_dim: " << int_dim
    << "  grid_dim: " << grid_dim << "\n\n";

  out1 << "Cumulatives:\n";
  out1 << "  si: " << si << "  si2: " << si2 << "  swgt: " << swgt << "  schi: " << schi << "  last_iter: " << last_iter << "\n";
  out1 << "  pre_si: " << pre_si << "  pre_si2: " << pre_si2 << "  pre_swgt: " << pre_swgt << "  pre_schi: " << pre_schi << "\n";
  out1 << "  last_neval: " << last_neval  << "  last_ti: " << last_ti << "  last_tsi: " << last_tsi << "  last_sdi: " << last_sdi << "\n\n";

  for(int ll=0; ll<int_dim; ll++){
      out1 << "\n data for axis  " << ll+1 << "\n"
          << "    x         di            d" << "\n";
    for(int hh=0; hh<grid_dim; hh++){
        out1 << "   " << std::fixed<<std::setprecision(6) << grid[ll][hh]
            << "  " << std::scientific << di[ll][hh] << "  " << dd[ll][hh] << "\n";
    }
  }
  out1.unsetf(std::ios_base::floatfield);
  out1 << "\n\n";
}

int vegas_integral::copy_from_file(std::istream & infile){

  int return_code = 0;

  std::string line, word;
//
  std::getline(infile, line); // BLANK
  std::getline(infile, line);
  std::istringstream iss1(line);
  iss1 >> word >> alpha >> word >> int_dim >> word >> grid_dim;
//
  std::getline(infile, line); // BLANK
  std::getline(infile, line); // BLANK
  std::getline(infile, line);
  std::istringstream iss2(line);
  iss2 >> word >> si >> word >> si2 >> word >> swgt >> word >> schi
     >> word >> last_iter;
  std::getline(infile, line);
  std::istringstream iss3(line);
  iss3 >> word >> pre_si >> word >> pre_si2 >> word >> pre_swgt >> word >> pre_schi;
  std::getline(infile, line);
  std::istringstream iss4(line);
  iss4 >> word >> last_neval >> word >> last_ti >> word >> last_tsi >> word >> last_sdi;
//
  grid.resize(int_dim);
  dd.resize(int_dim);
  di.resize(int_dim);
  for(int i=0; i<int_dim; i++){
    grid.at(i).resize(grid_dim);
    dd.at(i).resize(grid_dim);
    di.at(i).resize(grid_dim);
  }
  std::getline(infile, line); // BLANK
  for(int ll=0; ll<int_dim; ll++){
    std::getline(infile, line); // BLANK
    std::getline(infile, line); // BLANK
    std::getline(infile, line);
    for(int hh=0; hh<grid_dim; hh++){
      std::getline(infile, line);
      std::istringstream issL(line);
      issL >> grid[ll][hh] >> di[ll][hh] >> dd[ll][hh];
    }
  }

  if(infile.fail()) return_code = 1;
  return return_code;

}


//////////////////////////////////


//  std::default_random_engine generator;
//  std::uniform_real_distribution<double> distribution(0.0,1.0);

void vegas_random_init(int seed_1){
//  generator.seed(seed_1);
  randin_(&seed_1);
}

void vegas_random(int dim, double* xx){
//  for (int i=0; i<dim; i++) xx[i] = distribution(generator);
  randa_(&dim,xx);
  // std::cout << "dim = " << dim << std::endl;
  // std::cout << "xx = " << xx[0] << " " << xx[1] << " " << xx[2] << " " << xx[3] << " " << xx[4] << std::endl;
}
