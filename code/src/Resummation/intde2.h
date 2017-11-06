//header file for intde2.c

#ifndef INTDE2_H
#define INTDE2_H

#include <cmath>


void intdeoini(int lenaw, double tiny, double eps, double *aw);

void intdeo(double (*f)(double, void* userdata),
            double a, double omega, double *aw, double *i, double *err,
            void* userdata);


#endif
