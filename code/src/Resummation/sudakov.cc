//Calculates Quark and Gluon Sudakov form factors Sq, Sg, these are the N-independent, i.e. soft, parts of the exponential

#include "sudakov.h"

std::complex<double> Sc(std::complex<double> b, double qsquared, double b0p, double alphas, double muresumsquared, double a, resummationIndepfns* resufns, int order, double gqnp, double ggnp, char qorg, int verbosity) {
    double q = 0;
    double A1q=0, A1g=0, A2q=0, A2g=0, B1q=0, B1g=0, A3q=0, A3g=0, C1qqn=0, C1ggn=0, beta0=0, beta1=0, beta2=0, B2q=0, B2g=0, fac=0;
    std::complex<double> S, bstar, blog, blim, y, f0c, f1c, f2c;

    S = std::complex<double>(0.0,0.0);
    bstar = std::complex<double>(0.0,0.0);
    blog = std::complex<double>(0.0,0.0);
    blim = std::complex<double>(0.0,0.0);
    y = std::complex<double>(0.0,0.0);
    f0c = std::complex<double>(0.0,0.0);
    f1c = std::complex<double>(0.0,0.0);
    f2c = std::complex<double>(0.0,0.0);
    beta0 = resufns->beta0;
    q = std::pow(qsquared,0.5);

    blim = b0p*1.0/q*std::exp(1.0/(2.0*alphas*beta0));
    bstar=b/std::pow(1.0+(b*b)/(blim*blim),0.5); //resolves singularities at large b
    blog = std::log(q*q*bstar*bstar/(b0p*b0p)+1.0); //modified sudakov resolves singularities at small b

    y = beta0*alphas*blog;

    if (verbosity >= 15) {
      std::cout << "sudakov: " << std::endl;
      std::cout << "b = " << b << std::endl;
      std::cout << "qsquared = " << qsquared << std::endl;
      std::cout << "b0p = " << b0p << std::endl;
      std::cout << "alphas = " << alphas << std::endl;
      std::cout << "muresumsquared = " << muresumsquared << std::endl;
      std::cout << "a = " << a << std::endl;
      std::cout << "beta0 = " << beta0 << std::endl;
      std::cout << "alphas = " << alphas << std::endl;
      std::cout << "blog = " << blog << std::endl;
      std::cout << "b = " << b << std::endl;
      std::cout << "y = " << y << std::endl;
      std::cout << "ggnp = " << ggnp << std::endl;
      std::cout << "gqnp = " << gqnp << std::endl;
    }

    A1q = resufns->A1q;
    A1g = resufns->A1g;
    A2q = resufns->A2q;
    A2g = resufns->A2g;
    B1q = resufns->B1q;
    B1g = resufns->B1g;
    A3q = resufns->A3q;
    A3g = resufns->A3g;
    C1qqn = resufns->C1qqn;
    C1ggn = resufns->C1ggn;
    beta0 = resufns->beta0;
    beta1 = resufns->beta1;
    beta2 = resufns->beta2;
    B2q = resufns->B2q;
    B2g = resufns->B2g;
    f0c = f0(y, beta0, A1q, A1g, qorg, verbosity);
    f1c = f1(y, qsquared, a, beta0, beta1, A1q, A1g, A2q, A2g, B1q, B1g, muresumsquared, qorg, verbosity);
    f2c = f2(y, qsquared, a, beta0, beta1, beta2, A1q, A1g, A2q, A2g, B1q, B1g, A3q, A3g, B2q, B2g, C1qqn, C1ggn, muresumsquared, qorg, verbosity);

    if (verbosity >= 50) {
	std::cout << "f0c = " << f0c << std::endl;
	std::cout << "f1c = " << f1c << std::endl;
	std::cout << "f2c = " << f2c << std::endl;
    }

    if (order == 0) {
	S = std::exp(blog*f0c);
    }
    else if (order == 1) {
	S = std::exp(blog*f0c+f1c);
    }
    else if (order == 2) {
	S = std::exp(blog*f0c + f1c + alphas*f2c);
    }
    if (qorg == 'q') {fac = gqnp;}
    else if (qorg == 'g') {fac = ggnp;}
    S=S*std::exp(-fac*b*b);
    if (verbosity >= 50) {
	std::cout << "Sq exponent pieces: " << blog*f0c << " " << f1c << " " << alphas*f2c << " " << -fac*b*b << std::endl;
	std::cout << "S should = " << std::exp(blog*f0c + f1c + alphas*f2c-fac*b*b) << std::endl;
	std::cout << "fac = " << fac << " b*b = " << b*b << std::endl;
	std::cout << "S = " << S << std::endl;
    }

    return S;
}

std::complex<double> f0(std::complex<double> y, double beta0, double A1q, double A1g, char qorg, int verbosity) 
//soft resummation LL coefficient
{
    std::complex<double> f0, OneComplex;
    double A1 = 0;
    f0 = std::complex<double>(0.0,0.0);
    OneComplex = std::complex<double>(1.0,0.0);
    if (qorg == 'q') {A1 = A1q;}
    else if (qorg == 'g') {A1 = A1g;}
    f0 = A1/beta0*(y+std::log(OneComplex-y))/y;
    return f0;
}

std::complex<double> f1(std::complex<double> y, double q2, double a, double beta0, double beta1, double A1q, double A1g, double A2q, double A2g, double B1q, double B1g, double mur2, char qorg, int verbosity) 
//soft resummation NLL coefficient, now mur dependence
{
    std::complex<double> f1, OneComplex;
    double A1=0, A2=0, B1=0;
    f1 = std::complex<double>(0.0,0.0);
    OneComplex = std::complex<double>(1.0,0.0);
    if (qorg == 'q') {
	A1 = A1q;
	A2 = A2q;
	B1 = B1q;
    }
    else if (qorg == 'g') {
	A1 = A1g;
	A2 = A2g;
	B1 = B1g;
    }
    f1 = A1*beta1/(beta0*beta0*beta0)*(0.5*std::log(OneComplex-y)*std::log(OneComplex-y) + y/(OneComplex-y) + std::log(OneComplex-y)/(OneComplex-y)) - A2/(beta0*beta0)*(std::log(OneComplex-y)+y/(OneComplex-y)) + B1/beta0*std::log(OneComplex-y) + A1/beta0*(y/(OneComplex-y)+std::log(OneComplex-y))*std::log(q2/mur2);
    f1 = f1 - 2.0*OneComplex*std::log(a)*A1/beta0*y/(OneComplex-y); //a dependence
    return f1;
}

std::complex<double> f2(std::complex<double> y, double q2, double a, double beta0, double beta1, double beta2, double A1q, double A1g, double A2q, double A2g, double B1q, double B1g, double A3q, double A3g, double B2q, double B2g, double C1qqn, double C1ggn, double mur2, char qorg, int verbosity) 
//soft resummation NNLL coefficient, now mur dependence
{
    std::complex<double> f2, OneComplex;
    double A1=0, A2=0, B1=0, B2=0, A3=0, C1n=0;
    f2 = std::complex<double>(0.0,0.0);
    OneComplex = std::complex<double>(1.0,0.0);
    if (qorg == 'q') {
	A1 = A1q;
	A2 = A2q;
	B1 = B1q;
	B2 = B2q;
	A3 = A3q;
	C1n = C1qqn;
    }
    else if (qorg == 'g') {
	A1 = A1g;
	A2 = A2g;
	B1 = B1g;
	B2 = B2g;
	A3 = A3g;
	C1n = C1ggn;
    }
    f2 = ((A2*beta1)/(beta0*beta0*beta0))*((y/2.0)*((3.0*y-2.0)/((OneComplex-y)*(OneComplex-y))) - ((OneComplex-2.0*y)*log(OneComplex-y)/((OneComplex-y)*(OneComplex-y)))) - B2/beta0*(y/(OneComplex-y)) + B1*beta1/(beta0*beta0)*(y/(OneComplex-y) + log(OneComplex-y)/(OneComplex-y)) - A3/(2.0*beta0*beta0)*y*y/((OneComplex-y)*(OneComplex - y)) + A1*((beta1*beta1/(2.0*beta0*beta0*beta0*beta0)*(OneComplex-2.0*y)/((OneComplex-y)*(OneComplex-y))*log(OneComplex-y)*log(OneComplex-y) + log(OneComplex-y)*((beta0*beta2-beta1*beta1)/(beta0*beta0*beta0*beta0) + beta1*beta1/(beta0*beta0*beta0*beta0*(OneComplex-y)))) + y/(2.0*beta0*beta0*beta0*beta0*(OneComplex-y)*(OneComplex-y))*(beta0*beta2*(2.0-3.0*y)+beta1*beta1*y)) - A1/2.0*y*y/((OneComplex-y)*(OneComplex-y))*log(q2/mur2)*log(q2/mur2) + log(q2/mur2)*(B1*y/(OneComplex-y)+A2/beta0*y*y/((OneComplex-y)*(OneComplex-y)) + A1*beta1/(beta0*beta0)*(y/(OneComplex-y) + (OneComplex-2.0*y)/((OneComplex-y)*(OneComplex-y))*log(OneComplex-y))) + 2.0*C1n*y/(OneComplex-y);
    f2 = f2+2.0*A1*y*(y-2.0*OneComplex)/((OneComplex-y)*(OneComplex-y))*log(a)*log(a) - log(a)*(2.0*B1*y/(OneComplex-y)+2.0*y/beta0*A2/((OneComplex-y)*(OneComplex-y)) - 2.0*A1*beta1/(beta0*beta0)*y*log(OneComplex-y)/((OneComplex-y)*(OneComplex-y))) + A1*log(a)*log(q2/mur2)*y*2.0/((OneComplex-y)*(OneComplex-y));
    return f2;

}


