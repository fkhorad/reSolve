#ifndef _constantsH_
#define _constantsH_


// An effort has been made to avoid the use global variables/constants in the code except where absolutely necessary. Collected here are those few numerical values that should have no business being changed at input file level.

namespace k_constants{

// Only reason to change these is if a different precision is desired;
// recompiling seems reasonable in this case
const double zeta2= 1.6449340668482262;
const double zeta3= 1.2020569031595942;
const double zeta4= 1.082323233711138;
const double Eulerconst = 0.57721566490153;
const double pi = 3.141592653589793238462643383279503;

// Parameters for PDF fit -- no point in changing them with current fit code
const int nfitpars = 8;
const int nfitmax = 14;
const double aa = 2.5;
const int max_fit_files = 100;

// Parameter for inverse Mellin transform -- no point in changing it with current Mellin integration code
const int mellin_points = 136;


//Constants for diphoton hard factor -- to be reviewed
//N[PolyLog[2, -25]]
const double PolyLog2m25 = -6.785907899712118;
// N[PolyLog[3, -25]]
const double PolyLog3m25 = -10.893189584831829;
//N[PolyLog[4, -25]]
const double PolyLog4m25 = -14.848948252983412;
// (*N[PolyLog[2,-10]]*)
const double PolyLog2m10 = 4.198277886858104;
// (*N[PolyLog[3,-10]]*)
const double PolyLog3m10 = -5.921064803756975;
// (*N[PolyLog[3,-5]]*)
const double PolyLog3m5 = -3.537511437618607;
// (*N[PolyLog[2,-5]]*)
const double PolyLog2m5 = -2.749279126060808;
// (*N[PolyLog[4,-5]]*)
const double PolyLog4m5 = -4.10646797909497;
// (*N[PolyLog[4,-10]]*)
const double PolyLog4m10 = -7.326570248027082;
// (*N[PolyLog[4,-2]]*)
const double PolyLog4m2 = -1.813126015328491;
// (* N[PolyLog[3,-2]]*)
const double PolyLog3m2 = -1.6682833639665713;
// (*N[PolyLog[2,-2]]*)
const double PolyLog2m2 = -1.4367463668836808;
/* To Hardcode (just fractions, man!)
const double Qu = 0.66666666666666666666;
const double Qd = -0.33333333333333333333;
const double Qu4 = 0.19753086419753086419;
const double Qd4 = 0.01234567901234567901;
const double Tr = 0.5;
const double sumQq2 = 1.2222222222222222222;
*/

const double gf=1.16639e-5; // To input

}

#endif
