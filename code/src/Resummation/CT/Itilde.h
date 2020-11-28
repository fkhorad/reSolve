// MODIFIED BY KORA -- HAVE TO DECIDE WHAT TO DO WITH THE LICENCE

//      Itilde.cpp
//
//      Copyright 2010 Leandro <leandro@ubuntu>
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include <iostream>
//#include <cstdlib>
//#include "wadp.h"
#include <vector>
//#include <stdio.h>
#include <cmath>

#include "ADPINT.h"

double Itilde(int m,double xmio);

double zBK(int ,double );
double zbesselk0(double);
double zbesselk1(double);
double zbesselk2(double);
double zbesselk3(double);

double F0(double variable, double z, int n, void* data);
