/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __MATERIAL
#define __MATERIAL

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "CSR.H"
#include "Utilities.H"

class Material{

public:

    int **neighbor_table,*no_orb,no_cation,no_anion;
    double *bond_length,*atomic_mass;
    double *alpha_ph,*beta_ph,*kappa_ph,*tau_ph,*gamma_ph;
    double *alphap_ph,*betap_ph,*taup_ph;
    double *delta1p_ph,*delta3p_ph,*delta4p_ph;
    double **phiR,**phiT,**phiZ;
    double **mid_gap_energy,**band_gap_table;
    
    Material(const char*,const char*,int,double*,int,double);
    ~Material();
    
    void initialize(int);
    void get_band_edge(double*,double*);
    void get_d_matrix(CSR*,int,int,int);
    void get_nd_matrix(CSR*,double,double,double,double,int,int,int);
    void get_n_and_nd_matrix_derivative(double*,CSR*,CSR*,double,double,double,double,\
					int,int,int,int);
    void get_strain_shift(double*,double*,double*,double*,double*,int,\
			  int,int);
    
    // For strain calculation
    double alpha,beta;
    double ideal_a0;
    int TB;
    int IS,IPX,IPY,IPZ,ISS,ID1,ID2,ID3,ID4,ID5;
    
protected:
    
    char *mat_name;
    char *table_file;
    double binary_x[3];
    double *Esa, *Epa, *Estara, *Eda12, *Eda15, *lambdaa, *Vsssa, *Vstarstarsa, *Vsstarsa;
    double *Epa1, *Epa2, *Epa3, *Eda1, *Eda2, *Eda3, *Eda4, *Eda5;
    double *Vspsa, *Vstarpsa, *Vsdsa, *Vstardsa, *Vppsa, *Vpppa, *Vpdsa, *Vpdpa, *Vddsa;
    double *Vddpa, *Vddda;
    double *Csasc, *Cstarastarc, *Csastarc, *Csapc, *Cstarapc, *Csadc, *Cstaradc, *Cpapc;
    double *Cpadc, *Cdadc;
    double *Esc, *Epc, *Estarc, *Edc12, *Edc15, *lambdac, *Vsssc, *Vstarstarsc, *Vsstarsc;
    double *Epc1, *Epc2, *Epc3, *Edc1, *Edc2, *Edc3, *Edc4, *Edc5;
    double *Vspsc, *Vstarpsc, *Vsdsc, *Vstardsc, *Vppsc, *Vpppc, *Vpdsc, *Vpdpc, *Vddsc;
    double *Vddpc, *Vdddc;
    double *Cscsa, *Cstarcstara, *Cscstara, *Cscpa, *Cstarcpa, *Cscda, *Cstarcda, *Cpcpa;
    double *Cpcda, *Cdcda;
    double *eta_sss, *eta_sstars, *eta_starstars, *eta_sps, *eta_starps, *eta_sds;
    double *eta_stards, *eta_pps, *eta_ppp, *eta_pds, *eta_pdp, *eta_dds, *eta_ddp;
    double *eta_ddd;
    double Eg, ECmin, EVmax, *Eshift;
    double *r2av_pp,*r2av_dd,*r4av_dd,*r_pd,*r3_pd,*Zeff;
    double amu;
    double mu_scal;
    double Temp;
    int N3D;
    int strain_model;
    int read_hamiltonian;
    int sc_dist_dep;
    int init_gap_tables;    
    int mod_scaling_fcn;    
    int table_dim, no_material, mat_code;
    int sp3d5ss,sp3d5,sp3ss,sp3,sorb,pzorb,pzdxzdyz;

    double si_temp_dep(double);
    double graphene_scaling_function(double,double);
    void calc_matrix_element(double*,double,double,double,double,int,int);
    void prepare_strain_constants(double*,int,int);
    void read_material_from_file(FILE*);
    void get_Hdd_shift(double*,double*,double,double,double,double,double);
    void get_Hpp_shift(double*,double*,double,double,double,double);
    void get_Hpd_shift(double*,double*,double,double,double,double,double);
    void Y_2_4_combos(double*,double*,double*);
    
};

#endif
