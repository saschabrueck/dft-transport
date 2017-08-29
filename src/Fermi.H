/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __FERMI
#define __FERMI

#include <mpi.h>
#include "Utilities.H"
#include "Types.H"

class Fermi
{
public:
    
    Fermi();
    ~Fermi();

    double find_fermi(double,double,double*,double,int,double*,int,double*,int,int,double,\
		      double,double);
    void density(double*,double,double*,double,int,double*,int,double*,int,int,double,double,\
		 int,double,double,double);
    void derivate(double*,double,double*,double,int,double*,int,double*,int,int,double,double,\
		  int,double,double,double);

private:
    double Eps0,e,kB;

    double Newton(double,double,double,double*,double,int,double*,int,double*,int,int,double,\
		  double,double);
    double get_edge(double*,double,int,int,int);
};

#endif
