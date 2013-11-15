#include "Material.H"
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

Material::Material(const char* pmat_name,double *pmat_binary_x,int pstrain_model,double pTemp)
{

    amu          = 1.0/(6.022e26*scale_norm);
    N3D          = 3;

    strain_model = pstrain_model;

    Temp         = pTemp;

    mat_name     = new char[255];

    strcpy(mat_name,pmat_name);

    c_dcopy(3,pmat_binary_x,1,binary_x,1);
    
    mat_code = -1;
    
    if(strcmp(mat_name,"Si")==0){
        mat_code = 0;
    }
    if(strcmp(mat_name,"GaAs")==0){
        mat_code = 1;
    }
    if(strcmp(mat_name,"GaAs_AlGaAs")==0){
        mat_code = 2;
    }
    if(strcmp(mat_name,"Ge")==0){
        mat_code = 3;
    }
    if(strcmp(mat_name,"SiGe")==0){
        mat_code = 4;
    }
    if((strcmp(mat_name,"InAs_GaAs_AlAs")==0)||(strcmp(mat_name,"InAs_InGaAs_InAlAs")==0)||(strcmp(mat_name,"InAs")==0)){
        mat_code = 5;
    }
    if(strcmp(mat_name,"Si_P")==0){
        mat_code = 6;
    }
    if(strcmp(mat_name,"Si_SiO2")==0){
        mat_code = 7;
    }
    if(strcmp(mat_name,"Si_sp3ss_VB")==0){
        mat_code = 8;
    }
    if((strcmp(mat_name,"InAs_sp3ss")==0)||(strcmp(mat_name,"InGaAs_sp3ss")==0)){
        mat_code = 9;
    }
    if(strcmp(mat_name,"Si_sp3ss_CB")==0){
        mat_code = 10;
    }
    if(strcmp(mat_name,"Ge_Vogl_sp3ss")==0){
        mat_code = 11;
    }
    if(strcmp(mat_name,"GaAs_Vogl_sp3ss")==0){
        mat_code = 12;
    }
    if(strcmp(mat_name,"InAs_Vogl_sp3ss")==0){
        mat_code = 13;
    }
    if(strcmp(mat_name,"InSb")==0){
        mat_code = 14;
    }
    if(strcmp(mat_name,"InSb_sp3ss")==0){
        mat_code = 15;
    }
    if(strcmp(mat_name,"InAs_GaSb_sp3ss")==0){
        mat_code = 16;
    }
    if(strcmp(mat_name,"GaSb")==0){
        mat_code = 17;
    }
    if(strcmp(mat_name,"AlN")==0){
        mat_code = 20;
    }
    if((strcmp(mat_name,"GaN")==0)||(strcmp(mat_name,"GaN_AlGaN")==0)){
        mat_code = 21;
    }
    if(strcmp(mat_name,"InN")==0){
        mat_code = 22;
    }
    if((strcmp(mat_name,"graphene_pz")==0)||(strcmp(mat_name,"carbon_pz")==0)){
        mat_code = 50;
    }
    if((strcmp(mat_name,"graphene_sp3d5ss")==0)||(strcmp(mat_name,"cnt_sp3d5ss")==0)||\
       (strcmp(mat_name,"carbon_sp3d5ss")==0)){
        mat_code = 51;
    }
    if((strcmp(mat_name,"graphene_sp3")==0)||(strcmp(mat_name,"cnt_sp3")==0)||\
       (strcmp(mat_name,"carbon_sp3")==0)){
        mat_code = 52;
    }
    if((strcmp(mat_name,"graphene_pzdxzdyz")==0)||(strcmp(mat_name,"cnt_pzdxzdyz")==0)||\
       (strcmp(mat_name,"carbon_pzdxzdyz")==0)){
        mat_code = 53;
    }
    if(strcmp(mat_name,"cnt_pz")==0){
        mat_code = 54;
    }
    if((strcmp(mat_name,"graphene_pzdxzdyz_hydro_pz")==0)||\
       (strcmp(mat_name,"cnt_pzdxzdyz_hydro_pz")==0)||\
       (strcmp(mat_name,"carbon_pzdxzdyz_hydro_pz")==0)){
        mat_code = 55;
    }
    if((strcmp(mat_name,"graphene_py")==0)||(strcmp(mat_name,"carbon_py")==0)){
        mat_code = 56;
    }
    if((strcmp(mat_name,"graphene_sigma_pi")==0)||\
       (strcmp(mat_name,"bilayer_sigma_pi")==0)||\
       (strcmp(mat_name,"multilayer_sigma_pi")==0)||\
       (strcmp(mat_name,"cnt_sigma_pi")==0)||\
       (strcmp(mat_name,"carbon_sigma_pi")==0)){
        mat_code = 57;
    }
    if((strcmp(mat_name,"bilayer_pz")==0)||\
       (strcmp(mat_name,"multilayer_pz")==0)){
        mat_code = 58;
    }
    if(strcmp(mat_name,"Bi2Te3")==0){
        mat_code = 100;
    }
    if(strcmp(mat_name,"Si_new")==0){
        mat_code = 200;
    }
        
    TB       = 10;
    sp3d5ss  = 10;
    sp3d5    = 9;
    sp3ss    = 5;
    sp3      = 4;
    sorb     = 1;
    pzorb    = 1;
    pzdxzdyz = 3;

    IS       = 0;
    IPX      = 1;
    IPY      = 2;
    IPZ      = 3;
    ISS      = 4;
    ID1      = 5;
    ID2      = 6;
    ID3      = 7;
    ID4      = 8;
    ID5      = 9;

}

/************************************************************************************************/

Material::~Material()
{
    for(int i=0;i<table_dim;i++){
        delete[] neighbor_table[i];
	delete[] mid_gap_energy[i];
	delete[] band_gap_table[i];
    }
    for(int i=0;i<2*max(no_anion,no_cation);i++){
        delete[] phiR[i];
	delete[] phiT[i];
	delete[] phiZ[i];
    }
    delete [] mat_name;
    delete [] neighbor_table;
    delete [] mid_gap_energy;
    delete [] band_gap_table;
    delete [] phiR;
    delete [] phiT;
    delete [] phiZ;
    delete [] bond_length;
    delete [] no_orb;
    delete [] atomic_mass;
    delete [] alpha_ph;
    delete [] beta_ph;
    delete [] kappa_ph;
    delete [] tau_ph;
    delete [] gamma_ph;
    delete [] alphap_ph;
    delete [] betap_ph;
    delete [] taup_ph;
    delete [] delta1p_ph;
    delete [] delta3p_ph;
    delete [] delta4p_ph;
    delete [] Esa;
    delete [] Epa;
    delete [] Epa1;
    delete [] Epa2;
    delete [] Epa3;
    delete [] Estara;
    delete [] Eda12;
    delete [] Eda15;
    delete [] Eda1;
    delete [] Eda2;
    delete [] Eda3;
    delete [] Eda4;
    delete [] Eda5;
    delete [] lambdaa;
    delete [] Vsssa;
    delete [] Vstarstarsa;
    delete [] Vsstarsa;
    delete [] Vspsa;
    delete [] Vstarpsa;
    delete [] Vsdsa;
    delete [] Vstardsa;
    delete [] Vppsa;
    delete [] Vpppa;
    delete [] Vpdsa;
    delete [] Vpdpa;
    delete [] Vddsa;
    delete [] Vddpa;
    delete [] Vddda;
    delete [] Csasc;
    delete [] Cstarastarc;
    delete [] Csastarc;
    delete [] Csapc;
    delete [] Cstarapc;
    delete [] Csadc;
    delete [] Cstaradc;
    delete [] Cpapc;
    delete [] Cpadc;
    delete [] Cdadc;
    delete [] Esc;
    delete [] Epc;
    delete [] Epc1;
    delete [] Epc2;
    delete [] Epc3;
    delete [] Estarc;
    delete [] Edc12;
    delete [] Edc15;
    delete [] Edc1;
    delete [] Edc2;
    delete [] Edc3;
    delete [] Edc4;
    delete [] Edc5;
    delete [] lambdac;
    delete [] Vsssc;
    delete [] Vstarstarsc;
    delete [] Vsstarsc;
    delete [] Vspsc;
    delete [] Vstarpsc;
    delete [] Vsdsc;
    delete [] Vstardsc;
    delete [] Vppsc;
    delete [] Vpppc;
    delete [] Vpdsc;
    delete [] Vpdpc;
    delete [] Vddsc;
    delete [] Vddpc;
    delete [] Vdddc;
    delete [] Cscsa;
    delete [] Cstarcstara;
    delete [] Cscstara;
    delete [] Cscpa;
    delete [] Cstarcpa;
    delete [] Cscda;
    delete [] Cstarcda;
    delete [] Cpcpa;
    delete [] Cpcda;
    delete [] Cdcda;
    delete [] eta_sss;
    delete [] eta_starps;
    delete [] eta_sstars;
    delete [] eta_starstars;
    delete [] eta_sps;
    delete [] eta_sds;
    delete [] eta_stards;
    delete [] eta_pps;
    delete [] eta_ppp;
    delete [] eta_pds;
    delete [] eta_pdp;
    delete [] eta_dds;
    delete [] eta_ddp;
    delete [] eta_ddd;
    delete [] Eshift;
    delete [] r2av_pp;
    delete [] r2av_dd;
    delete [] r4av_dd;
    delete [] r_pd;
    delete [] r3_pd;
    delete [] Zeff;

}

/************************************************************************************************/

void Material::initialize(int psc_dist_dep)
{

    FILE *F;

    switch(mat_code){
    case 0:
        no_anion  = 2;
        no_cation = 2;
        break;
    case 1:
        no_anion  = 1;
        no_cation = 1;
        break;    
    case 2:
        no_anion  = 1;
        no_cation = 2;
        break;
    case 3:
        no_anion  = 1;
        no_cation = 1;
        break;
    case 4:
        no_anion  = 2;
        no_cation = 2;
        break;
    case 5:
        no_anion  = 1;
        no_cation = 3;
        break;
    case 6:
        no_anion  = 2;
        no_cation = 2;
        break;
    case 7:
        no_anion  = 2;
	no_cation = 2;
	break;
    case 8:
        no_anion  = 1;
        no_cation = 1;
	break;
    case 9:
        no_anion  = 1;
	no_cation = 1;
	break;
    case 10: 
        no_anion  = 1;
        no_cation = 1;
        break;   
    case 11: 
        no_anion  = 1;
        no_cation = 1;
        break;   
    case 12: 
        no_anion  = 1;
        no_cation = 1;
        break;
    case 13: 
        no_anion  = 1;
        no_cation = 1;
        break;
    case 14: 
        no_anion  = 1;
        no_cation = 1;
        break;
    case 15:
        no_anion  = 1;
	no_cation = 1;
	break;
    case 16:
        no_anion  = 2;
	no_cation = 2;
	break;
    case 17:
        no_anion  = 1;
	no_cation = 1;
	break;
    case 20: 
        no_anion  = 2;
        no_cation = 2;
        break;    
    case 21: 
        no_anion  = 2;
        no_cation = 4;
        break;    
    case 22: 
        no_anion  = 2;
        no_cation = 2;
        break;
    case 50: 
        no_anion  = 1;
        no_cation = 1;
        break;
    case 51: 
        no_anion  = 1;
        no_cation = 1;
        break;
    case 52: 
        no_anion  = 1;
        no_cation = 1;
        break;
    case 53: 
        no_anion  = 1;
        no_cation = 1;
        break;
    case 54: 
        no_anion  = 1;
        no_cation = 1;
        break;
    case 55:
        no_anion  = 2;
	no_cation = 2;
	break;
    case 56: 
        no_anion  = 1;
        no_cation = 1;
        break;
    case 57: 
        no_anion  = 2;
        no_cation = 2;
        break;
    case 58: 
        no_anion  = 2;
        no_cation = 2;
        break;	    
    case 100: 
        no_anion  = 2;
        no_cation = 3;
        break;
    case 200: 
        no_anion  = 1;
        no_cation = 1;
        break;
    default:
        F = fopen(mat_name,"r");
	if(F!=0){
	    fscanf(F,"%i",&no_anion);
	    fscanf(F,"%i",&no_cation);
        }else{
	    printf("This material is not defined and you have not provided a parameter file\n");
	    exit(0);
	}
    }	
 
    no_material    = no_anion*no_cation;
    table_dim      = Round(max(2*no_cation,2*(no_anion-1)+1));
    
    bond_length    = new double[no_material];
    neighbor_table = new int*[table_dim];
    mid_gap_energy = new double*[table_dim];
    band_gap_table = new double*[table_dim];
    no_orb         = new int[2*max(no_anion,no_cation)];
    atomic_mass    = new double[2*max(no_anion,no_cation)];
    alpha_ph       = new double[2*max(no_anion,no_cation)];
    beta_ph        = new double[2*max(no_anion,no_cation)];
    kappa_ph       = new double[2*max(no_anion,no_cation)];
    tau_ph         = new double[2*max(no_anion,no_cation)];
    gamma_ph       = new double[2*max(no_anion,no_cation)];
    alphap_ph      = new double[2*max(no_anion,no_cation)];
    betap_ph       = new double[2*max(no_anion,no_cation)];
    taup_ph        = new double[2*max(no_anion,no_cation)];
    delta1p_ph     = new double[2*max(no_anion,no_cation)];
    delta3p_ph     = new double[2*max(no_anion,no_cation)];
    delta4p_ph     = new double[2*max(no_anion,no_cation)];
    phiR           = new double*[2*max(no_anion,no_cation)];
    phiT           = new double*[2*max(no_anion,no_cation)];
    phiZ           = new double*[2*max(no_anion,no_cation)];

    for(int i=0;i<table_dim;i++){
        neighbor_table[i] = new int[table_dim];
	mid_gap_energy[i] = new double[table_dim];
	band_gap_table[i] = new double[table_dim];
    }

    for(int i=0;i<2*max(no_anion,no_cation);i++){
	phiR[i]           = new double[4];
	phiT[i]           = new double[4];
	phiZ[i]           = new double[4];
    }

    sc_dist_dep   = psc_dist_dep;
    
    Esa           = new double[no_material];
    Epa           = new double[no_material];
    Epa1          = new double[no_material];
    Epa2          = new double[no_material];
    Epa3          = new double[no_material];
    Estara        = new double[no_material];
    Eda12         = new double[no_material];
    Eda15         = new double[no_material];
    Eda1          = new double[no_material];
    Eda2          = new double[no_material];
    Eda3          = new double[no_material];
    Eda4          = new double[no_material];
    Eda5          = new double[no_material];
    lambdaa       = new double[no_material];
    Vsssa         = new double[no_material];
    Vstarstarsa   = new double[no_material];
    Vsstarsa      = new double[no_material];
    Vspsa         = new double[no_material];
    Vstarpsa      = new double[no_material];
    Vsdsa         = new double[no_material];
    Vstardsa      = new double[no_material];
    Vppsa         = new double[no_material];
    Vpppa         = new double[no_material];
    Vpdsa         = new double[no_material];
    Vddsa         = new double[no_material];
    Vpdpa         = new double[no_material];
    Vddpa         = new double[no_material];
    Vddda         = new double[no_material];
    Csasc         = new double[no_material];
    Cstarastarc   = new double[no_material];
    Csastarc      = new double[no_material];
    Csapc         = new double[no_material];
    Cstarapc      = new double[no_material];
    Csadc         = new double[no_material];
    Cstaradc      = new double[no_material];
    Cpapc         = new double[no_material];
    Cpadc         = new double[no_material];
    Cdadc         = new double[no_material];
    Esc           = new double[no_material];
    Epc           = new double[no_material];
    Epc1          = new double[no_material];
    Epc2          = new double[no_material];
    Epc3          = new double[no_material];
    Estarc        = new double[no_material];
    Edc12         = new double[no_material];
    Edc15         = new double[no_material];
    Edc1          = new double[no_material];
    Edc2          = new double[no_material];
    Edc3          = new double[no_material];
    Edc4          = new double[no_material];
    Edc5          = new double[no_material];
    lambdac       = new double[no_material];
    Vsssc         = new double[no_material];
    Vstarstarsc   = new double[no_material];
    Vsstarsc      = new double[no_material];
    Vspsc         = new double[no_material];
    Vstarpsc      = new double[no_material];
    Vsdsc         = new double[no_material];
    Vstardsc      = new double[no_material];
    Vppsc         = new double[no_material];
    Vpppc         = new double[no_material];
    Vpdsc         = new double[no_material];
    Vpdpc         = new double[no_material];
    Vddsc         = new double[no_material];
    Vddpc         = new double[no_material];
    Vdddc         = new double[no_material];
    Cscsa         = new double[no_material];
    Cstarcstara   = new double[no_material];
    Cscstara      = new double[no_material];
    Cscpa         = new double[no_material];
    Cstarcpa      = new double[no_material];
    Cscda         = new double[no_material];
    Cstarcda      = new double[no_material];
    Cpcpa         = new double[no_material];
    Cpcda         = new double[no_material];
    Cdcda         = new double[no_material];
    eta_sss       = new double[no_material];
    eta_sstars    = new double[no_material];
    eta_starstars = new double[no_material];
    eta_sps       = new double[no_material];
    eta_starps    = new double[no_material];
    eta_sds       = new double[no_material];
    eta_stards    = new double[no_material];
    eta_pps       = new double[no_material];
    eta_ppp       = new double[no_material];
    eta_pds       = new double[no_material];
    eta_pdp       = new double[no_material];
    eta_dds       = new double[no_material];
    eta_ddp       = new double[no_material];
    eta_ddd       = new double[no_material];
    Eshift        = new double[no_material];
    r2av_pp       = new double[no_material];
    r2av_dd       = new double[no_material];
    r4av_dd       = new double[no_material];
    r_pd          = new double[no_material];
    r3_pd         = new double[no_material];
    Zeff          = new double[no_material];

    init_var(r2av_pp,no_material);
    init_var(r2av_dd,no_material);
    init_var(r4av_dd,no_material);
    init_var(r_pd,no_material);
    init_var(r3_pd,no_material);
    init_var(Zeff,no_material);

    init_var(alpha_ph,2*max(no_anion,no_cation));
    init_var(beta_ph,2*max(no_anion,no_cation));
    init_var(kappa_ph,2*max(no_anion,no_cation));
    init_var(tau_ph,2*max(no_anion,no_cation));
    init_var(gamma_ph,2*max(no_anion,no_cation));

    init_var(alphap_ph,2*max(no_anion,no_cation));
    init_var(betap_ph,2*max(no_anion,no_cation));
    init_var(taup_ph,2*max(no_anion,no_cation));
    init_var(delta1p_ph,2*max(no_anion,no_cation));
    init_var(delta3p_ph,2*max(no_anion,no_cation));
    init_var(delta4p_ph,2*max(no_anion,no_cation));

    //Sui and Herman, Phys. Rev. B 48, 17938 (1993)

    //alphap_ph  = 2/(sqrt(3)*bond_length)*n_alpha*alpha_ph
    //betap_ph   = -4*sqrt(3)/bond_length*l_beta*beta_ph
    //taup_ph    = sqrt(3)/(2*bond_length)*(n_tau-2*l_tau)*tau_ph
    //delta1p_ph = sqrt(3)/bond_length*(m_beta-2*l_beta)*beta_ph
    //delta3p_ph = -12*sqrt(3)/bond_length*l_tau*tau_ph
    //delta4p_ph = 0 (much smaller than the others)
    
    init_gap_tables = 1;
    mod_scaling_fcn = 0;
    mu_scal         = 0.0;
       
    switch(mat_code){
        
    case 0:  // Si (with H passivation)

        Eg                   = 1.13;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 48.5;
        beta                 = 13.8;
        ideal_a0             = 0.543095;

//Anion and from anion to
        Esa[0]               = -2.15168;
        Epa[0]               = 4.22925;
        Estara[0]            = 19.11650;
        Eda12[0]             = 13.78950;
	Eda15[0]             = 13.78950;
        lambdaa[0]           = 0.01989;     //(Delta/3)
        Vsssa[0]             = -1.95933;
        Vstarstarsa[0]       = -4.24135;
        Vsstarsa[0]          = -1.52230;
        Vspsa[0]             = 3.02562;
        Vstarpsa[0]          = 3.15565;
        Vsdsa[0]             = -2.28485;
        Vstardsa[0]          = -0.80993;
        Vppsa[0]             = 4.10364;
        Vpppa[0]             = -1.51801;
        Vpdsa[0]             = -1.35554;
        Vpdpa[0]             = 2.38479;
        Vddsa[0]             = -1.68136;
        Vddpa[0]             = 2.58880;
        Vddda[0]             = -1.81400;
//Strain
        Csasc[0]             = 1.6788;
        Cstarastarc[0]       = 0.7777;
        Csastarc[0]          = 1.7843;
        Csapc[0]             = 0.4801;
        Cstarapc[0]          = 3.5888;
        Csadc[0]             = 0;
        Cstaradc[0]          = 0.3421;
        Cpapc[0]             = 4.0664;
        Cpadc[0]             = 0;
        Cdadc[0]             = 5.2973;

//Cation and from cation to anion
        Esc[0]               = -2.15168;
        Epc[0]               = 4.22925;
        Estarc[0]            = 19.11650;
        Edc12[0]             = 13.78950;
	Edc15[0]             = 13.78950;
        lambdac[0]           = 0.01989;     //(Delta/3)
        Vsssc[0]             = -1.95933;
        Vstarstarsc[0]       = -4.24135;
        Vsstarsc[0]          = -1.52230;
        Vspsc[0]             = 3.02562;
        Vstarpsc[0]          = 3.15565;
        Vsdsc[0]             = -2.28485;
        Vstardsc[0]          = -0.80993;
        Vppsc[0]             = 4.10364;
        Vpppc[0]             = -1.51801;
        Vpdsc[0]             = -1.35554;
        Vpdpc[0]             = 2.38479;
        Vddsc[0]             = -1.68136;
        Vddpc[0]             = 2.58880;
        Vdddc[0]             = -1.81400;
//Strain
        Cscsa[0]             = 1.6788;
        Cstarcstara[0]       = 0.7777;
        Cscstara[0]          = 1.7843;
        Cscpa[0]             = 0.4801;
        Cstarcpa[0]          = 3.5888;
        Cscda[0]             = 0;
        Cstarcda[0]          = 0.3421;
        Cpcpa[0]             = 4.0664;
        Cpcda[0]             = 0;
        Cdcda[0]             = 5.2973;

        eta_sss[0]           = 0.562469;
        eta_sstars[0]        = 0.132030;
        eta_starstars[0]     = 0.192369;
        eta_sps[0]           = 2.365484;
        eta_starps[0]        = 0.344918;
        eta_sds[0]           = 2.567199;
        eta_stards[0]        = 1.086006;
        eta_pps[0]           = 0.494354;
        eta_ppp[0]           = 1.843851;
        eta_pds[0]           = 2.236358;
        eta_pdp[0]           = 4.512500;
        eta_dds[0]           = 4.668357;
        eta_ddp[0]           = 2.302375;
        eta_ddd[0]           = 0.923911;

        Eshift[0]            = 27;

        bond_length[0]       = 0.543*sqrt(3.0)/4.0;

//Si and from Si to H
        Esa[1]               = -2.15168;
        Epa[1]               = 4.22925;
        Estara[1]            = 19.11650;
        Eda12[1]             = 13.78950;
	Eda15[1]             = 13.78950;
        lambdaa[1]           = 0.01989;     //(Delta/3)
        Vsssa[1]             = -3.99972;
        Vstarstarsa[1]       = 0.0;
        Vsstarsa[1]          = 0.0;
        Vspsa[1]             = 0.0;
        Vstarpsa[1]          = 0.0;
        Vsdsa[1]             = 0.0;
        Vstardsa[1]          = 0.0;
        Vppsa[1]             = 0.0;
        Vpppa[1]             = 0.0;
        Vpdsa[1]             = 0.0;
        Vpdpa[1]             = 0.0;
        Vddsa[1]             = 0.0;
        Vddpa[1]             = 0.0;
        Vddda[1]             = 0.0;
//Strain
        Csasc[1]             = 0.0;
        Cstarastarc[1]       = 0.0;
        Csastarc[1]          = 0.0;
        Csapc[1]             = 0.0;
        Cstarapc[1]          = 0.0;
        Csadc[1]             = 0.0;
        Cstaradc[1]          = 0.0;
        Cpapc[1]             = 0.0;
        Cpadc[1]             = 0.0;
        Cdadc[1]             = 0.0;

//H and from H to Si
        Esc[1]               = 0.99984;
        Epc[1]               = 0.0;
        Estarc[1]            = 0.0;
        Edc12[1]             = 0.0;
	Edc15[1]             = 0.0;
        lambdac[1]           = 0.0;     //(Delta/3)
        Vsssc[1]             = -3.99972;
        Vstarstarsc[1]       = 0.0;
        Vsstarsc[1]          = -1.6977;
        Vspsc[1]             = 4.25175;
        Vstarpsc[1]          = 0.0;
        Vsdsc[1]             = -2.10552;
        Vstardsc[1]          = 0.0;
        Vppsc[1]             = 0.0;
        Vpppc[1]             = 0.0;
        Vpdsc[1]             = 0.0;
        Vpdpc[1]             = 0.0;
        Vddsc[1]             = 0.0;
        Vddpc[1]             = 0.0;
        Vdddc[1]             = 0.0;
//Strain
        Cscsa[1]             = 0.0;
        Cstarcstara[1]       = 0.0;
        Cscstara[1]          = 0.0;
        Cscpa[1]             = 0.0;
        Cstarcpa[1]          = 0.0;
        Cscda[1]             = 0.0;
        Cstarcda[1]          = 0.0;
        Cpcpa[1]             = 0.0;
        Cpcda[1]             = 0.0;
        Cdcda[1]             = 0.0;

        eta_sss[1]           = 0.0;
        eta_sstars[1]        = 0.0;
        eta_starstars[1]     = 0.0;
        eta_sps[1]           = 0.0;
        eta_starps[1]        = 0.0;
        eta_sds[1]           = 0.0;
        eta_stards[1]        = 0.0;
        eta_pps[1]           = 0.0;
        eta_ppp[1]           = 0.0;
        eta_pds[1]           = 0.0;
        eta_pdp[1]           = 0.0;
        eta_dds[1]           = 0.0;
        eta_ddp[1]           = 0.0;
        eta_ddd[1]           = 0.0;

        Eshift[1]            = 27;

        bond_length[1]       = 0.543*sqrt(3.0)/4.0/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 11;
	neighbor_table[0][3] = 21;
        neighbor_table[3][0] = 22;
	neighbor_table[1][2] = 21;
        neighbor_table[2][1] = 22;

	no_orb[0]            = sp3d5ss; //Si anion
	no_orb[1]            = sp3d5ss; //Si cation
	no_orb[2]            = sorb;    //H anion
	no_orb[3]            = sorb;    //H cation

	atomic_mass[0]       = 28.0855*amu;  //Si
	atomic_mass[1]       = 28.0855*amu;  //Si
	atomic_mass[2]       = 1.00794*amu;  //H
	atomic_mass[3]       = 1.00794*amu;  //H

	alpha_ph[0]          = 49.4;
	beta_ph[0]           = 4.79;
	kappa_ph[0]          = 6.99;
	tau_ph[0]            = 5.2;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 49.4;
	beta_ph[1]           = 4.79;
	kappa_ph[1]          = 6.99;
	tau_ph[1]            = 5.2;
	gamma_ph[1]          = 0.0;

	alpha_ph[2]          = 49.4;
	beta_ph[2]           = 4.79;
	kappa_ph[2]          = 6.99;
	tau_ph[2]            = 5.2;
	gamma_ph[2]          = 0.0;

	alpha_ph[3]          = 49.4;
	beta_ph[3]           = 4.79;
	kappa_ph[3]          = 6.99;
	tau_ph[3]            = 5.2;
	gamma_ph[3]          = 0.0;

	alphap_ph[0]         = -2.3241e12/64.0/1.5*pow(Temp/300.0,0.9);
	betap_ph[0]          = -1.0388e12/64.0/1.5*pow(Temp/300.0,0.9);
	taup_ph[0]           = -2.2447e11/64.0/1.5*pow(Temp/300.0,0.9);
	delta1p_ph[0]        = -2.2300e11/64.0/1.5*pow(Temp/300.0,0.9);
	delta3p_ph[0]        = 6.3894e11/64.0/1.5*pow(Temp/300.0,0.9);

	alphap_ph[1]         = -2.3241e12/64.0/1.5*pow(Temp/300.0,0.9);
	betap_ph[1]          = -1.0388e12/64.0/1.5*pow(Temp/300.0,0.9);
	taup_ph[1]           = -2.2447e11/64.0/1.5*pow(Temp/300.0,0.9);
	delta1p_ph[1]        = -2.2300e11/64.0/1.5*pow(Temp/300.0,0.9);
	delta3p_ph[1]        = 6.3894e11/64.0/1.5*pow(Temp/300.0,0.9);
	
	alphap_ph[2]         = -2.3241e12/64.0/1.5*pow(Temp/300.0,0.9);
	betap_ph[2]          = -1.0388e12/64.0/1.5*pow(Temp/300.0,0.9);
	taup_ph[2]           = -2.2447e11/64.0/1.5*pow(Temp/300.0,0.9);
	delta1p_ph[2]        = -2.2300e11/64.0/1.5*pow(Temp/300.0,0.9);
	delta3p_ph[2]        = 6.3894e11/64.0/1.5*pow(Temp/300.0,0.9);

	alphap_ph[3]         = -2.3241e12/64.0/1.5*pow(Temp/300.0,0.9);
	betap_ph[3]          = -1.0388e12/64.0/1.5*pow(Temp/300.0,0.9);
	taup_ph[3]           = -2.2447e11/64.0/1.5*pow(Temp/300.0,0.9);
	delta1p_ph[3]        = -2.2300e11/64.0/1.5*pow(Temp/300.0,0.9);
	delta3p_ph[3]        = 6.3894e11/64.0/1.5*pow(Temp/300.0,0.9);

        break;

    case 1:  // GaAs

        Eg                   = 1.424;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 41.49;
        beta                 = 8.94;
        ideal_a0             = 0.56532;

//Anion and from anion to
        Esa[0]               = -5.50042;
        Epa[0]               = 4.15107;
        Estara[0]            = 19.71059;
        Eda12[0]             = 13.03169;
	Eda15[0]             = 13.03169;
        lambdaa[0]           = 0.17234;     //(Delta/3)
        Vsssa[0]             = -1.64508;
        Vstarstarsa[0]       = -3.67720;
        Vsstarsa[0]          = -1.31491;
        Vspsa[0]             = 2.66493;
        Vstarpsa[0]          = 1.97650;
        Vsdsa[0]             = -2.58357;
        Vstardsa[0]          = -0.62820;
        Vppsa[0]             = 4.15080;
        Vpppa[0]             = -1.42744;
        Vpdsa[0]             = -1.87428;
        Vpdpa[0]             = 2.52926;
        Vddsa[0]             = -1.26996;
        Vddpa[0]             = 2.50536;
        Vddda[0]             = -0.85174;
//Strain
        Csasc[0]             = 0.586960;
        Cstarastarc[0]       = 0.486090;
        Csastarc[0]          = 0.770950;
        Csapc[0]             = 0.759790;
        Cstarapc[0]          = 0.810790;
        Csadc[0]             = 1.070150;
        Cstaradc[0]          = 1.032560;
        Cpapc[0]             = 2;
        Cpadc[0]             = 1.613500;
        Cdadc[0]             = 1.262620;

//Cation and from cation to anion
        Esc[0]               = -0.24119;
        Epc[0]               = 6.70776;
        Estarc[0]            = 22.66352;
        Edc12[0]             = 12.74846;
	Edc15[0]             = 12.74846;
        lambdac[0]           = 0.02179;     //(Delta/3)
        Vsssc[0]             = -1.64508;
        Vstarstarsc[0]       = -3.67720;
        Vsstarsc[0]          = -2.20777;
        Vspsc[0]             = 2.96032;
        Vstarpsc[0]          = 1.02755;
        Vsdsc[0]             = -2.32059;
        Vstardsc[0]          = -0.13324;
        Vppsc[0]             = 4.15080;
        Vpppc[0]             = -1.42744;
        Vpdsc[0]             = -1.88964;
        Vpdpc[0]             = 2.54913;
        Vddsc[0]             = -1.26996;
        Vddpc[0]             = 2.50536;
        Vdddc[0]             = -0.85174;
//Strain
        Cscsa[0]             = 0.586960;
        Cstarcstara[0]       = 0.486090;
        Cscstara[0]          = 0.889210;
        Cscpa[0]             = 1.458910;
        Cstarcpa[0]          = 1.212020;
        Cscda[0]             = 0.580530;
        Cstarcda[0]          = 1.323850;
        Cpcpa[0]             = 2;
        Cpcda[0]             = 1.5;
        Cdcda[0]             = 1.26262;

        eta_sss[0]           = 2.060010;
        eta_sstars[0]        = 0;
        eta_starstars[0]     = 0.212660;
        eta_sps[0]           = 1.384980;
        eta_starps[0]        = 1.399300;
        eta_sds[0]           = 1.898890;
        eta_stards[0]        = 1.785400;
        eta_pps[0]           = 2.684970;
        eta_ppp[0]           = 1.314050;
        eta_pds[0]           = 1.812350;
        eta_pdp[0]           = 2.379640;
        eta_dds[0]           = 1.724430;
        eta_ddp[0]           = 1.972530;
        eta_ddd[0]           = 1.896720;

        Eshift[0]            = 27;

        bond_length[0]       = 0.56533*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;

	no_orb[0]            = sp3d5ss; //As
	no_orb[1]            = sp3d5ss; //Ga

	atomic_mass[0]       = 74.9216*amu;  //As
	atomic_mass[1]       = 69.7230*amu;  //Ga
        
        break;
    case 2:  // GaAs-AlxGa1-xAs

        Eg                   = 1.424;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 41.49;
        beta                 = 8.94;
        ideal_a0             = 0.56532;

//Anion 1 (As) and from anion 1 to cation 1 (Ga)
        Esa[0]               = -5.50042;
        Epa[0]               = 4.15107;
        Estara[0]            = 19.71059;
        Eda12[0]             = 13.03169;
	Eda15[0]             = 13.03169;
        lambdaa[0]           = 0.17234;     //(Delta/3)
        Vsssa[0]             = -1.64508;
        Vstarstarsa[0]       = -3.67720;
        Vsstarsa[0]          = -1.31491;
        Vspsa[0]             = 2.66493;
        Vstarpsa[0]          = 1.97650;
        Vsdsa[0]             = -2.58357;
        Vstardsa[0]          = -0.62820;
        Vppsa[0]             = 4.15080;
        Vpppa[0]             = -1.42744;
        Vpdsa[0]             = -1.87428;
        Vpdpa[0]             = 2.52926;
        Vddsa[0]             = -1.26996;
        Vddpa[0]             = 2.50536;
        Vddda[0]             = -0.85174;
//Strain
        Csasc[0]             = 0.586960;
        Cstarastarc[0]       = 0.486090;
        Csastarc[0]          = 0.770950;
        Csapc[0]             = 0.759790;
        Cstarapc[0]          = 0.810790;
        Csadc[0]             = 1.070150;
        Cstaradc[0]          = 1.032560;
        Cpapc[0]             = 2;
        Cpadc[0]             = 1.613500;
        Cdadc[0]             = 1.262620;

//Cation 1 (Ga) and from cation 1 to anion 1 (As)
        Esc[0]               = -0.24119;
        Epc[0]               = 6.70776;
        Estarc[0]            = 22.66352;
        Edc12[0]             = 12.74846;
	Edc15[0]             = 12.74846;
        lambdac[0]           = 0.02179;     //(Delta/3)
        Vsssc[0]             = -1.64508;
        Vstarstarsc[0]       = -3.67720;
        Vsstarsc[0]          = -2.20777;
        Vspsc[0]             = 2.96032;
        Vstarpsc[0]          = 1.02755;
        Vsdsc[0]             = -2.32059;
        Vstardsc[0]          = -0.13324;
        Vppsc[0]             = 4.15080;
        Vpppc[0]             = -1.42744;
        Vpdsc[0]             = -1.88964;
        Vpdpc[0]             = 2.54913;
        Vddsc[0]             = -1.26996;
        Vddpc[0]             = 2.50536;
        Vdddc[0]             = -0.85174;
//Strain
        Cscsa[0]             = 0.586960;
        Cstarcstara[0]       = 0.486090;
        Cscstara[0]          = 0.889210;
        Cscpa[0]             = 1.458910;
        Cstarcpa[0]          = 1.212020;
        Cscda[0]             = 0.580530;
        Cstarcda[0]          = 1.323850;
        Cpcpa[0]             = 2;
        Cpcda[0]             = 1.5;
        Cdcda[0]             = 1.26262;

        eta_sss[0]           = 2.060010;
        eta_sstars[0]        = 0;
        eta_starstars[0]     = 0.212660;
        eta_sps[0]           = 1.384980;
        eta_starps[0]        = 1.399300;
        eta_sds[0]           = 1.898890;
        eta_stards[0]        = 1.785400;
        eta_pps[0]           = 2.684970;
        eta_ppp[0]           = 1.314050;
        eta_pds[0]           = 1.812350;
        eta_pdp[0]           = 2.379640;
        eta_dds[0]           = 1.724430;
        eta_ddp[0]           = 1.972530;
        eta_ddd[0]           = 1.896720;

        Eshift[0]            = 27;

        bond_length[0]       = 0.56533*sqrt(3.0)/4.0;

//Anion 1 (As) and from anion 1 to cation 2 (AlxGa1-x)
        Esa[1]               = -5.17012*binary_x[1]-5.50042*(1-binary_x[1]);
        Epa[1]               = 4.39708*binary_x[1]+4.15107*(1-binary_x[1]);
        Estara[1]            = 19.80474*binary_x[1]+19.71059*(1-binary_x[1]);
        Eda12[1]             = 13.13880*binary_x[1]+13.03169*(1-binary_x[1]);
	Eda15[1]             = 13.13880*binary_x[1]+13.03169*(1-binary_x[1]);
        lambdaa[1]           = 0.17386*binary_x[1]+0.17234*(1-binary_x[1]);     //(Delta/3)
        Vsssa[1]             = -1.64584*binary_x[1]-1.64508*(1-binary_x[1]);
        Vstarstarsa[1]       = -2.84245*binary_x[1]-3.67720*(1-binary_x[1]);
        Vsstarsa[1]          = -2.78690*binary_x[1]-1.31491*(1-binary_x[1]);
        Vspsa[1]             = 3.02223*binary_x[1]+2.66493*(1-binary_x[1]);
        Vstarpsa[1]          = 1.92174*binary_x[1]+1.97650*(1-binary_x[1]);
        Vsdsa[1]             = -3.03196*binary_x[1]-2.58357*(1-binary_x[1]);
        Vstardsa[1]          = -1.84300*binary_x[1]-0.62820*(1-binary_x[1]);
        Vppsa[1]             = 4.53156*binary_x[1]+4.15080*(1-binary_x[1]);
        Vpppa[1]             = -1.86816*binary_x[1]-1.42744*(1-binary_x[1]);
        Vpdsa[1]             = -2.47345*binary_x[1]-1.87428*(1-binary_x[1]);
        Vpdpa[1]             = 2.52741*binary_x[1]+2.52926*(1-binary_x[1]);
        Vddsa[1]             = -1.97058*binary_x[1]-1.26996*(1-binary_x[1]);
        Vddpa[1]             = 1.67733 *binary_x[1]+2.50536*(1-binary_x[1]);
        Vddda[1]             = -1.58868 *binary_x[1]-0.85174*(1-binary_x[1]);
//Strain
	Csasc[1]             = 0.586960;
        Cstarastarc[1]       = 0.486090;
        Csastarc[1]          = 0.770950;
        Csapc[1]             = 0.759790;
        Cstarapc[1]          = 0.810790;
        Csadc[1]             = 1.070150;
        Cstaradc[1]          = 1.032560;
        Cpapc[1]             = 2;
        Cpadc[1]             = 1.613500;
        Cdadc[1]             = 1.262620;

//Cation 2 (AxGa1-x) and from cation 2 to anion 1 (As)
        Esc[1]               = 0.79695*binary_x[1]-0.24119*(1-binary_x[1]);
        Epc[1]               = 6.63291*binary_x[1]+6.70776*(1-binary_x[1]);
        Estarc[1]            = 24.16587*binary_x[1]+22.66352*(1-binary_x[1]);
        Edc12[1]             = 12.92122*binary_x[1]+12.74846*(1-binary_x[1]);
	Edc15[1]             = 12.92122*binary_x[1]+12.74846*(1-binary_x[1]);
        lambdac[1]           = 0.01589*binary_x[1]+0.02179*(1-binary_x[1]);     //(Delta/3)
        Vsssc[1]             = -1.64584*binary_x[1]-1.64508*(1-binary_x[1]);
        Vstarstarsc[1]       = -2.84245*binary_x[1]-3.67720*(1-binary_x[1]);
        Vsstarsc[1]          = -1.88341*binary_x[1]-2.20777*(1-binary_x[1]);
        Vspsc[1]             = 2.95309*binary_x[1]+2.96032*(1-binary_x[1]);
        Vstarpsc[1]          = 1.30469*binary_x[1]+1.02755*(1-binary_x[1]);
        Vsdsc[1]             = -2.64111*binary_x[1]-2.32059*(1-binary_x[1]);
        Vstardsc[1]          = -1.73510*binary_x[1]-0.13324*(1-binary_x[1]);
        Vppsc[1]             = 4.53156*binary_x[1]+4.15080*(1-binary_x[1]);
        Vpppc[1]             = -1.86816*binary_x[1]-1.42744*(1-binary_x[1]);
        Vpdsc[1]             = -1.02836*binary_x[1]-1.88964*(1-binary_x[1]);
        Vpdpc[1]             = 2.86419*binary_x[1]+2.54913*(1-binary_x[1]);
        Vddsc[1]             = -1.97058*binary_x[1]-1.26996*(1-binary_x[1]);
        Vddpc[1]             = 1.67733*binary_x[1]+2.50536*(1-binary_x[1]);
        Vdddc[1]             = -1.58868*binary_x[1]-0.85174*(1-binary_x[1]);
//Strain
	Cscsa[1]             = 0.586960;
        Cstarcstara[1]       = 0.486090;
        Cscstara[1]          = 0.889210;
        Cscpa[1]             = 1.458910;
        Cstarcpa[1]          = 1.212020;
        Cscda[1]             = 0.580530;
        Cstarcda[1]          = 1.323850;
        Cpcpa[1]             = 2;
        Cpcda[1]             = 1.5;
        Cdcda[1]             = 1.26262;

        eta_sss[1]           = 2.060010;
        eta_sstars[1]        = 0;
        eta_starstars[1]     = 0.212660;
        eta_sps[1]           = 1.384980;
        eta_starps[1]        = 1.399300;
        eta_sds[1]           = 1.898890;
        eta_stards[1]        = 1.785400;
        eta_pps[1]           = 2.684970;
        eta_ppp[1]           = 1.314050;
        eta_pds[1]           = 1.812350;
        eta_pdp[1]           = 2.379640;
        eta_dds[1]           = 1.724430;
        eta_ddp[1]           = 1.972530;
        eta_ddd[1]           = 1.896720;

        Eshift[1]            = 27;

        bond_length[1]       = 0.56533*sqrt(3.0)/4.0;

	neighbor_table[0][1] = 11;
        neighbor_table[0][3] = 21;
	neighbor_table[2][1] = 11;
        neighbor_table[2][3] = 21;
        neighbor_table[1][0] = 12;
        neighbor_table[3][0] = 22;
	neighbor_table[1][2] = 12;
        neighbor_table[3][2] = 22;

	no_orb[0]            = sp3d5ss; //As
	no_orb[1]            = sp3d5ss; //Ga
	no_orb[2]            = sp3d5ss; //As
	no_orb[3]            = sp3d5ss; //AlxGa1-x

	//atomic_mass[0]       = 74.9216*amu;  //As
	//atomic_mass[1]       = 69.7230*amu;  //Ga
	//atomic_mass[2]       = 74.9216*amu;  //As
	//atomic_mass[3]       = (26.981539*binary_x[1]+69.7230*(1-binary_x[1]))*amu; //AlxGa1-x

	//InAs Parameters
	atomic_mass[0]       = 74.9216*amu;
	atomic_mass[1]       = 114.8180*amu;
	atomic_mass[2]       = 74.9216*amu;
	atomic_mass[3]       = 114.8180*amu;

	alpha_ph[0]          = 34.6779;
	beta_ph[0]           = 2.6992;
	kappa_ph[0]          = 0.8682;
	tau_ph[0]            = 1.9896;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 34.6779;
	beta_ph[1]           = 2.6992;
	kappa_ph[1]          = 0.8682;
	tau_ph[1]            = 1.9896;
	gamma_ph[1]          = 0.0;

	alpha_ph[2]          = 34.6779;
	beta_ph[2]           = 2.6992;
	kappa_ph[2]          = 0.8682;
	tau_ph[2]            = 1.9896;
	gamma_ph[2]          = 0.0;

	alpha_ph[3]          = 34.6779;
	beta_ph[3]           = 2.6992;
	kappa_ph[3]          = 0.8682;
	tau_ph[3]            = 1.9896;
	gamma_ph[3]          = 0.0;

        break;

    case 3:  // Ge

        Eg                   = 0.664;
        ECmin                = 1.448;
        EVmax                = 0.77;

	// Strain constants
        alpha                = 39.0;
        beta                 = 12.0;
        ideal_a0             = 0.5657906;

//Anion and from anion to cation
        Esa[0]               = -1.95617;
        Epa[0]               = 5.30970;
        Estara[0]            = 19.29600;
        Eda12[0]             = 13.58060;
	Eda15[0]             = 13.58060;
        lambdaa[0]           = 0.10132;     //(Delta/3)
        Vsssa[0]             = -1.39456;
        Vstarstarsa[0]       = -3.56680;
        Vsstarsa[0]          = -2.01830;
        Vspsa[0]             = 2.73135;
        Vstarpsa[0]          = 2.68638;
        Vsdsa[0]             = -2.64779;
        Vstardsa[0]          = -1.12312;
        Vppsa[0]             = 4.28921;
        Vpppa[0]             = -1.73707;
        Vpdsa[0]             = -2.00115;
        Vpdpa[0]             = 2.10953;
        Vddsa[0]             = -1.32941;
        Vddpa[0]             = 2.56261;
        Vddda[0]             = -1.95120;
//Strain
        Csasc[0]             = 0;
        Cstarastarc[0]       = 6.28624;
        Csastarc[0]          = 1.86887;
        Csapc[0]             = 2.03278;
        Cstarapc[0]          = 6.28624;
        Csadc[0]             = 0.16396;
        Cstaradc[0]          = 1.98112;
        Cpapc[0]             = 0.42830;
        Cpadc[0]             = 0.12084;
        Cdadc[0]             = 3.85908;

//Cation and from cation to anion
        Esc[0]               = -1.95617;
        Epc[0]               = 5.30970;
        Estarc[0]            = 19.29600;
        Edc12[0]             = 13.58060;
	Edc15[0]             = 13.58060;
        lambdac[0]           = 0.10132;     //(Delta/3)
        Vsssc[0]             = -1.39456;
        Vstarstarsc[0]       = -3.56680;
        Vsstarsc[0]          = -2.01830;
        Vspsc[0]             = 2.73135;
        Vstarpsc[0]          = 2.68638;
        Vsdsc[0]             = -2.64779;
        Vstardsc[0]          = -1.12312;
        Vppsc[0]             = 4.28921;
        Vpppc[0]             = -1.73707;
        Vpdsc[0]             = -2.00115;
        Vpdpc[0]             = 2.10953;
        Vddsc[0]             = -1.32941;
        Vddpc[0]             = 2.56261;
        Vdddc[0]             = -1.95120;
//Strain
        Cscsa[0]             = 0;
        Cstarcstara[0]       = 6.28624;
        Cscstara[0]          = 1.86887;
        Cscpa[0]             = 2.03278;
        Cstarcpa[0]          = 6.28624;
        Cscda[0]             = 0.16396;
        Cstarcda[0]          = 1.98112;
        Cpcpa[0]             = 0.42830;
        Cpcda[0]             = 0.12084;
        Cdcda[0]             = 3.85908;

        eta_sss[0]           = 1.995511;
        eta_sstars[0]        = 0;
        eta_starstars[0]     = 2.388227;
        eta_sps[0]           = 1.293029;
        eta_starps[0]        = 5;
        eta_sds[0]           = 2.792438;
        eta_stards[0]        = 0.751342;
        eta_pps[0]           = 1.136409;
        eta_ppp[0]           = 1.748033;
        eta_pds[0]           = 2.687841;
        eta_pdp[0]           = 4.369211;
        eta_dds[0]           = 5;
        eta_ddp[0]           = 0.697690;
        eta_ddd[0]           = 3.062529;

        Eshift[0]            = 27.77;

        bond_length[0]       = 0.56579*sqrt(3.0)/4.0;

	neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;

	no_orb[0]            = sp3d5ss; //Ge anion
	no_orb[1]            = sp3d5ss; //Ge cation

	atomic_mass[0]       = 72.6100*amu;
	atomic_mass[1]       = 72.6100*amu;

	alpha_ph[0]          = 44.32;
	beta_ph[0]           = 3.68;
	kappa_ph[0]          = 6.13;
	tau_ph[0]            = 4.95;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 44.32;
	beta_ph[1]           = 3.68;
	kappa_ph[1]          = 6.13;
	tau_ph[1]            = 4.95;
	gamma_ph[1]          = 0.0;

	alphap_ph[0]         = -2.0283e12/64.0;
	betap_ph[0]          = -5.1097e11/64.0;
	taup_ph[0]           = -2.6841e11/64.0;
	delta1p_ph[0]        = -1.6651e10/64.0;
	delta3p_ph[0]        = 1.1758e11/64.0;

	alphap_ph[1]         = -2.0283e12/64.0;
	betap_ph[1]          = -5.1097e11/64.0;
	taup_ph[1]           = -2.6841e11/64.0;
	delta1p_ph[1]        = -1.6651e10/64.0;
	delta3p_ph[1]        = 1.1758e11/64.0;

        break;

    case 4:  // SiGe

        Eg                   = 1.13;
        ECmin                = Eg;
        EVmax                = 0.3;

	// Strain constants
        alpha                = 43.25;
        beta                 = 12.9;
        ideal_a0             = 0.55443;

//Anion 1 (Si) and from anion 1 to cation 1 (Si)
        Esa[0]               = -2.15168;
        Epa[0]               = 4.22925;
        Estara[0]            = 19.11650;
        Eda12[0]             = 13.78950;
	Eda15[0]             = 13.78950;
        lambdaa[0]           = 0.01989;     //(Delta/3)
        Vsssa[0]             = -1.95933;
        Vstarstarsa[0]       = -4.24135;
        Vsstarsa[0]          = -1.52230;
        Vspsa[0]             = 3.02562;
        Vstarpsa[0]          = 3.15565;
        Vsdsa[0]             = -2.28485;
        Vstardsa[0]          = -0.80993;
        Vppsa[0]             = 4.10364;
        Vpppa[0]             = -1.51801;
        Vpdsa[0]             = -1.35554;
        Vpdpa[0]             = 2.38479;
        Vddsa[0]             = -1.68136;
        Vddpa[0]             = 2.58880;
        Vddda[0]             = -1.81400;
//Strain
        Csasc[0]             = 1.68054;
        Cstarastarc[0]       = 0.77849;
        Csastarc[0]          = 1.78613;
        Csapc[0]             = 0.48057;
        Cstarapc[0]          = 3.59244;
        Csadc[0]             = 0;
        Cstaradc[0]          = 0.34243;
        Cpapc[0]             = 4.07053;
        Cpadc[0]             = 0;
        Cdadc[0]             = 5.30270;

//Cation 1 (Si) and from cation 1 to anion 1 (Si)
        Esc[0]               = -2.15168;
        Epc[0]               = 4.22925;
        Estarc[0]            = 19.11650;
        Edc12[0]             = 13.78950;
	Edc15[0]             = 13.78950;
        lambdac[0]           = 0.01989;     //(Delta/3)
        Vsssc[0]             = -1.95933;
        Vstarstarsc[0]       = -4.24135;
        Vsstarsc[0]          = -1.52230;
        Vspsc[0]             = 3.02562;
        Vstarpsc[0]          = 3.15565;
        Vsdsc[0]             = -2.28485;
        Vstardsc[0]          = -0.80993;
        Vppsc[0]             = 4.10364;
        Vpppc[0]             = -1.51801;
        Vpdsc[0]             = -1.35554;
        Vpdpc[0]             = 2.38479;
        Vddsc[0]             = -1.68136;
        Vddpc[0]             = 2.58880;
        Vdddc[0]             = -1.81400;
//Strain
        Cscsa[0]             = 1.68054;
        Cstarcstara[0]       = 0.77849;
        Cscstara[0]          = 1.78613;
        Cscpa[0]             = 0.48057;
        Cstarcpa[0]          = 3.59244;
        Cscda[0]             = 0;
        Cstarcda[0]          = 0.34243;
        Cpcpa[0]             = 4.07053;
        Cpcda[0]             = 0;
        Cdcda[0]             = 5.30270;

        eta_sss[0]           = 0.562469;
        eta_sstars[0]        = 0.132030;
        eta_starstars[0]     = 0.192369;
        eta_sps[0]           = 2.365484;
        eta_starps[0]        = 0.344918;
        eta_sds[0]           = 2.567199;
        eta_stards[0]        = 1.086006;
        eta_pps[0]           = 0.494354;
        eta_ppp[0]           = 1.843851;
        eta_pds[0]           = 2.236358;
        eta_pdp[0]           = 4.512500;
        eta_dds[0]           = 4.668357;
        eta_ddp[0]           = 2.302375;
        eta_ddd[0]           = 0.923911;

        Eshift[0]            = 27;

        bond_length[0]       = 0.543095*sqrt(3.0)/4.0;
	
//Anion 1 (Si) and from anion 1 to cation 2 (Ge)
        Esa[1]               = -2.15168;
        Epa[1]               = 4.22925;
        Estara[1]            = 19.11650;
        Eda12[1]             = 13.78950;
	Eda15[1]             = 13.78950;
        lambdaa[1]           = 0.01989;     //(Delta/3)
        Vsssa[1]             = -1.574309;
        Vstarstarsa[1]       = -3.999753;
        Vsstarsa[1]          = -1.440045;
        Vspsa[1]             = 2.878285;
        Vstarpsa[1]          = 2.747303;
        Vsdsa[1]             = -2.471630;
        Vstardsa[1]          = -1.045346;
        Vppsa[1]             = 4.201374;
        Vpppa[1]             = -1.626914;
        Vpdsa[1]             = -1.415419;
        Vpdpa[1]             = 2.023372;
        Vddsa[1]             = -1.317709;
        Vddpa[1]             = 2.501057;
        Vddda[1]             = -1.777856;
//Strain
        Csasc[1]             = 3.09934;
        Cstarastarc[1]       = 2.09695;
        Csastarc[1]          = 1.69301;
        Csapc[1]             = 2.21094;
        Cstarapc[1]          = 2.58548;
        Csadc[1]             = 1.88788;
        Cstaradc[1]          = 1.82653;
        Cpapc[1]             = 2.89710;
        Cpadc[1]             = 3.51673;
        Cdadc[1]             = 3.55811;

//Cation 2 (Ge) and from cation 2 to anion 1 (Si)
        Esc[1]               = -1.95617;
        Epc[1]               = 5.3097;
        Estarc[1]            = 19.295998;
        Edc12[1]             = 13.580600;
	Edc15[1]             = 13.580600;
        lambdac[1]           = 0.10132;     //(Delta/3)
        Vsssc[1]             = -1.574309;
        Vstarstarsc[1]       = -3.999753;
        Vsstarsc[1]          = -1.499091;
        Vspsc[1]             = 2.760826;
        Vstarpsc[1]          = 2.814474;
        Vsdsc[1]             = -2.59531;
        Vstardsc[1]          = -1.029594;
        Vppsc[1]             = 4.201374;
        Vpppc[1]             = -1.626914;
        Vpdsc[1]             = -2.022367;
        Vpdpc[1]             = 2.504041;
        Vddsc[1]             = -1.317709;
        Vddpc[1]             = 2.501057;
        Vdddc[1]             = -1.777856;
//Strain
        Cscsa[1]             = 3.09934;
        Cstarcstara[1]       = 2.09695;
        Cscstara[1]          = 3.16404;
        Cscpa[1]             = 2.25252;
        Cstarcpa[1]          = 2.39048;
        Cscda[1]             = 1.26252;
        Cstarcda[1]          = 1.77635;
        Cpcpa[1]             = 2.89710;
        Cpcda[1]             = 1.26053;
        Cdcda[1]             = 3.55811;

        eta_sss[1]           = 2.155355;
        eta_sstars[1]        = 0.154363;
        eta_starstars[1]     = 0;
        eta_sps[1]           = 2.157637;
        eta_starps[1]        = 1.413615;
        eta_sds[1]           = 1.831938;
        eta_stards[1]        = 2.538574;
        eta_pps[1]           = 2.621803;
        eta_ppp[1]           = 1.491652;
        eta_pds[1]           = 1.787961;
        eta_pdp[1]           = 2.257439;
        eta_dds[1]           = 1.814857;
        eta_ddp[1]           = 1.973166;
        eta_ddd[1]           = 1.805441;

        Eshift[1]            = 27;

        bond_length[1]       = (0.543095+0.56579)*sqrt(3.0)/8.0;

//Anion 2 (Ge) and from anion 2 to cation 1 (Si)
        Esa[2]               = -1.95617;
        Epa[2]               = 5.30970;
        Estara[2]            = 19.29600;
        Eda12[2]             = 13.58060;
	Eda15[2]             = 13.58060;
        lambdaa[2]           = 0.10132;     //(Delta/3)
        Vsssa[2]             = -1.574309;
        Vstarstarsa[2]       = -3.999753;
        Vsstarsa[2]          = -1.499091;
        Vspsa[2]             = 2.760826;
        Vstarpsa[2]          = 2.814474;
        Vsdsa[2]             = -2.59531;
        Vstardsa[2]          = -1.029594;
        Vppsa[2]             = 4.201374;
        Vpppa[2]             = -1.626914;
        Vpdsa[2]             = -2.022367;
        Vpdpa[2]             = 2.504041;
        Vddsa[2]             = -1.317709;
        Vddpa[2]             = 2.501057;
        Vddda[2]             = -1.777856;
//Strain
        Csasc[2]             = 3.09934;
        Cstarastarc[2]       = 2.09695;
        Csastarc[2]          = 3.16404;
        Csapc[2]             = 2.25252;
        Cstarapc[2]          = 2.39048;
        Csadc[2]             = 1.26252;
        Cstaradc[2]          = 1.77635;
        Cpapc[2]             = 2.89710;
        Cpadc[2]             = 1.26053;
        Cdadc[2]             = 3.55811;

//Cation 1 (Si) and from cation 1 to anion 2 (Ge)
        Esc[2]               = -2.151680;
        Epc[2]               = 4.22925;
        Estarc[2]            = 19.116498;
        Edc12[2]             = 13.7895;
	Edc15[2]             = 13.7895;
        lambdac[2]           = 0.01989;     //(Delta/3)
        Vsssc[2]             = -1.574309;
        Vstarstarsc[2]       = -3.999753;
        Vsstarsc[2]          = -1.440045;
        Vspsc[2]             = 2.878285;
        Vstarpsc[2]          = 2.747303;
        Vsdsc[2]             = -2.471630;
        Vstardsc[2]          = -1.045346;
        Vppsc[2]             = 4.201374;
        Vpppc[2]             = -1.626914;
        Vpdsc[2]             = -1.415419;
        Vpdpc[2]             = 2.023372;
        Vddsc[2]             = -1.317709;
        Vddpc[2]             = 2.501057;
        Vdddc[2]             = -1.777856;
//Strain
        Cscsa[2]             = 3.09934;
        Cstarcstara[2]       = 2.09695;
        Cscstara[2]          = 1.69301;
        Cscpa[2]             = 2.21094;
        Cstarcpa[2]          = 2.58548;
        Cscda[2]             = 1.88788;
        Cstarcda[2]          = 1.82653;
        Cpcpa[2]             = 2.89710;
        Cpcda[2]             = 3.51673;
        Cdcda[2]             = 3.55811;

        eta_sss[2]           = 2.155355;
        eta_sstars[2]        = 0.154363;
        eta_starstars[2]     = 0;
        eta_sps[2]           = 2.157637;
        eta_starps[2]        = 1.413615;
        eta_sds[2]           = 1.831938;
        eta_stards[2]        = 2.538574;
        eta_pps[2]           = 2.621803;
        eta_ppp[2]           = 1.491652;
        eta_pds[2]           = 1.787961;
        eta_pdp[2]           = 2.257439;
        eta_dds[2]           = 1.814857;
        eta_ddp[2]           = 1.973166;
        eta_ddd[2]           = 1.805441;

        Eshift[2]            = 27;

        bond_length[2]       = (0.56579+0.543095)*sqrt(3.0)/8.0;

        //Anion 2 (Ge) and from anion 2 to cation 2 (Ge)
        Esa[3]               = -1.95617;
        Epa[3]               = 5.30970;
        Estara[3]            = 19.29600;
        Eda12[3]             = 13.58060;
	Eda15[3]             = 13.58060;
        lambdaa[3]           = 0.10132;     //(Delta/3)
        Vsssa[3]             = -1.39456;
        Vstarstarsa[3]       = -3.56680;
        Vsstarsa[3]          = -2.01830;
        Vspsa[3]             = 2.73135;
        Vstarpsa[3]          = 2.68638;
        Vsdsa[3]             = -2.64779;
        Vstardsa[3]          = -1.12312;
        Vppsa[3]             = 4.28921;
        Vpppa[3]             = -1.73707;
        Vpdsa[3]             = -2.00115;
        Vpdpa[3]             = 2.10953;
        Vddsa[3]             = -1.32941;
        Vddpa[3]             = 2.56261;
        Vddda[3]             = -1.95120;
//Strain
        Csasc[3]             = 0;
        Cstarastarc[3]       = 6.28624;
        Csastarc[3]          = 1.86887;
        Csapc[3]             = 2.03278;
        Cstarapc[3]          = 6.28624;
        Csadc[3]             = 0.16396;
        Cstaradc[3]          = 1.98112;
        Cpapc[3]             = 0.42830;
        Cpadc[3]             = 0.12084;
        Cdadc[3]             = 3.85908;

//Cation 2 (Ge) and from cation 2 to anion 2 (Ge)
        Esc[3]               = -1.95617;
        Epc[3]               = 5.30970;
        Estarc[3]            = 19.29600;
        Edc12[3]             = 13.58060;
	Edc15[3]             = 13.58060;
        lambdac[3]           = 0.10132;     //(Delta/3)
        Vsssc[3]             = -1.39456;
        Vstarstarsc[3]       = -3.56680;
        Vsstarsc[3]          = -2.01830;
        Vspsc[3]             = 2.73135;
        Vstarpsc[3]          = 2.68638;
        Vsdsc[3]             = -2.64779;
        Vstardsc[3]          = -1.12312;
        Vppsc[3]             = 4.28921;
        Vpppc[3]             = -1.73707;
        Vpdsc[3]             = -2.00115;
        Vpdpc[3]             = 2.10953;
        Vddsc[3]             = -1.32941;
        Vddpc[3]             = 2.56261;
        Vdddc[3]             = -1.95120;
//Strain
        Cscsa[3]             = 0;
        Cstarcstara[3]       = 6.28624;
        Cscstara[3]          = 1.86887;
        Cscpa[3]             = 2.03278;
        Cstarcpa[3]          = 6.28624;
        Cscda[3]             = 0.16396;
        Cstarcda[3]          = 1.98112;
        Cpcpa[3]             = 0.42830;
        Cpcda[3]             = 0.12084;
        Cdcda[3]             = 3.85908;

        eta_sss[3]           = 1.995511;
        eta_sstars[3]        = 0;
        eta_starstars[3]     = 2.388227;
        eta_sps[3]           = 1.293029;
        eta_starps[3]        = 5;
        eta_sds[3]           = 2.792438;
        eta_stards[3]        = 0.751342;
        eta_pps[3]           = 1.136409;
        eta_ppp[3]           = 1.748033;
        eta_pds[3]           = 2.687841;
        eta_pdp[3]           = 4.369211;
        eta_dds[3]           = 5;
        eta_ddp[3]           = 0.697690;
        eta_ddd[3]           = 3.062529;

        Eshift[3]            = 27.77;

        bond_length[3]       = 0.56579*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[0][3] = 21;
        neighbor_table[1][0] = 12;
        neighbor_table[1][2] = 32;
        neighbor_table[2][1] = 31;
        neighbor_table[2][3] = 41;
        neighbor_table[3][0] = 22;
        neighbor_table[3][2] = 42;

	no_orb[0]            = sp3d5ss; //Si anion
	no_orb[1]            = sp3d5ss; //Si cation
	no_orb[2]            = sp3d5ss; //Ge anion
	no_orb[3]            = sp3d5ss; //Ge cation

	atomic_mass[0]       = 28.0855*amu;  //Si anion
	atomic_mass[1]       = 28.0855*amu;  //Si cation
	atomic_mass[2]       = 72.6100*amu;  //Ge anion
	atomic_mass[3]       = 72.6100*amu;  //Ge cation

	alpha_ph[0]          = 49.4;
	beta_ph[0]           = 4.79;
	kappa_ph[0]          = 6.99;
	tau_ph[0]            = 5.2;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 49.4;
	beta_ph[1]           = 4.79;
	kappa_ph[1]          = 6.99;
	tau_ph[1]            = 5.2;
	gamma_ph[1]          = 0.0;

	alpha_ph[2]          = 44.32;
	beta_ph[2]           = 3.68;
	kappa_ph[2]          = 6.13;
	tau_ph[2]            = 4.95;
	gamma_ph[2]          = 0.0;

	alpha_ph[3]          = 44.32;
	beta_ph[3]           = 3.68;
	kappa_ph[3]          = 6.13;
	tau_ph[3]            = 4.95;
	gamma_ph[3]          = 0.0;
	
	alphap_ph[0]         = -2.3241e12/64.0;
	betap_ph[0]          = -1.0388e12/64.0;
	taup_ph[0]           = -2.2447e11/64.0;
	delta1p_ph[0]        = -2.2300e11/64.0;
	delta3p_ph[0]        = 6.3894e11/64.0;

	alphap_ph[1]         = -2.3241e12/64.0;
	betap_ph[1]          = -1.0388e12/64.0;
	taup_ph[1]           = -2.2447e11/64.0;
	delta1p_ph[1]        = -2.2300e11/64.0;
	delta3p_ph[1]        = 6.3894e11/64.0;

	alphap_ph[2]         = -2.0283e12/64.0;
	betap_ph[2]          = -5.1097e11/64.0;
	taup_ph[2]           = -2.6841e11/64.0;
	delta1p_ph[2]        = -1.6651e10/64.0;
	delta3p_ph[2]        = 1.1758e11/64.0;

	alphap_ph[3]         = -2.0283e12/64.0;
	betap_ph[3]          = -5.1097e11/64.0;
	taup_ph[3]           = -2.6841e11/64.0;
	delta1p_ph[3]        = -1.6651e10/64.0;
	delta3p_ph[3]        = 1.1758e11/64.0;

        break;

    case 5:  // InAs_GaAs_AlAs

        Eg                   = 0.37;
        ECmin                = 0.59570;
        EVmax                = 0.22430;

        alpha                = 41.49;    // AlAs, GaAs
        beta                 = 8.94;     // AlAs, GaAs
        ideal_a0             = 0.56532;  // AlAs, GaAs
        //alpha              = 35.18;    // InAs
        //beta               = 5.49;     // InAs
        //ideal_a0           = 0.60583;  // InAs

//Anion (As) and from anion to cation 1 (In)
        Esa[0]               = -5.50042;
        Epa[0]               = 4.15107;
        Estara[0]            = 19.71059;
        Eda12[0]             = 13.03169;
	Eda15[0]             = 13.03169;
        lambdaa[0]           = 0.17234;     //(Delta/3)
        Vsssa[0]             = -1.69435;
        Vstarstarsa[0]       = -4.210450;
        Vsstarsa[0]          = -1.15987;
        Vspsa[0]             = 2.598230;
        Vstarpsa[0]          = 2.067660;
        Vsdsa[0]             = -2.268370;
        Vstardsa[0]          = -0.899370;
        Vppsa[0]             = 4.31064;
        Vpppa[0]             = -1.288950;
        Vpdsa[0]             = -1.73141;
        Vpdpa[0]             = 2.188860;
        Vddsa[0]             = -1.584610;
        Vddpa[0]             = 2.717930;
        Vddda[0]             = -0.505090;
//Strain
        Csasc[0]             = 1.258286;
        Cstarastarc[0]       = 2.481447;
        Csastarc[0]          = 4.557774;
        Csapc[0]             = 4.367575;
        Cstarapc[0]          = 3.298598;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 1.195042;
        Cpapc[0]             = 4.624438;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.246999;

//Cation 1 (In) and from cation 1 to anion (As)
        Esc[0]               = -0.581930;
        Epc[0]               = 6.971630;
        Estarc[0]            = 19.941380;
        Edc12[0]             = 13.307090;
	Edc15[0]             = 13.307090;
        lambdac[0]           = 0.131200;     //(Delta/3)
        Vsssc[0]             = -1.694350;
        Vstarstarsc[0]       = -4.210450;
        Vsstarsc[0]          = -2.426740;
        Vspsc[0]             = 2.809360;
        Vstarpsc[0]          = 0.937340;
        Vsdsc[0]             = -2.293090;
        Vstardsc[0]          = -0.488990;
        Vppsc[0]             = 4.310640;
        Vpppc[0]             = -1.288950;
        Vpdsc[0]             = -1.978420;
        Vpdpc[0]             = 2.456020;
        Vddsc[0]             = -1.584610;
        Vddpc[0]             = 2.717930;
        Vdddc[0]             = -0.505090;
//Strain
        Cscsa[0]             = 1.258286;
        Cstarcstara[0]       = 2.481447;
        Cscstara[0]          = 1.086223;
        Cscpa[0]             = 7.029660;
        Cstarcpa[0]          = 7.029496;
        Cscda[0]             = 0.187036;
        Cstarcda[0]          = 1.769483;
        Cpcpa[0]             = 4.624438;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.246999;

        eta_sss[0]           = 1.924940;
        eta_sstars[0]        = 0.060800;
        eta_starstars[0]     = 0.000810;
        eta_sps[0]           = 1.570030;
        eta_starps[0]        = 1.949370;
        eta_sds[0]           = 1.765660;
        eta_stards[0]        = 2.023870;
        eta_pps[0]           = 2.061510;
        eta_ppp[0]           = 1.602470;
        eta_pds[0]           = 2.383820;
        eta_pdp[0]           = 2.455600;
        eta_dds[0]           = 2.322910;
        eta_ddp[0]           = 1.615890;
        eta_ddd[0]           = 2.329600;

        Eshift[0]            = 27;

        bond_length[0]       = 0.60583*sqrt(3.0)/4.0;

	//InxGa1-xAs
//Anion (As) and from anion to cation 2 (Ga)
	Esa[1]               = -5.50042*binary_x[1]-5.50042*(1-binary_x[1]);
        Epa[1]               = 4.15107*binary_x[1]+4.15107*(1-binary_x[1]);
        Estara[1]            = 19.71059*binary_x[1]+19.71059*(1-binary_x[1]);
        Eda12[1]             = 13.03169*binary_x[1]+13.03169*(1-binary_x[1]);
	Eda15[1]             = 13.03169*binary_x[1]+13.03169*(1-binary_x[1]);
        lambdaa[1]           = 0.17234*binary_x[1]+0.17234*(1-binary_x[1]);     //(Delta/3)
        Vsssa[1]             = -1.69435*binary_x[1]-1.64508*(1-binary_x[1])-0.026152007*binary_x[1]*(1-binary_x[1]);
        Vstarstarsa[1]       = -4.21045*binary_x[1]-3.67720*(1-binary_x[1])-0.28501808*binary_x[1]*(1-binary_x[1]);
        Vsstarsa[1]          = -1.15987*binary_x[1]-1.31491*(1-binary_x[1])-0.10224072*binary_x[1]*(1-binary_x[1]);
        Vspsa[1]             = 2.59823*binary_x[1]+2.66493*(1-binary_x[1])-0.066696639*binary_x[1]*(1-binary_x[1]);
        Vstarpsa[1]          = 2.06766*binary_x[1]+1.97650*(1-binary_x[1])+0.09113023*binary_x[1]*(1-binary_x[1]);
        Vsdsa[1]             = -2.26837*binary_x[1]-2.58357*(1-binary_x[1])+0.31519991*binary_x[1]*(1-binary_x[1]);
        Vstardsa[1]          = -0.89937*binary_x[1]-0.62820*(1-binary_x[1])-0.27115*binary_x[1]*(1-binary_x[1]);
        Vppsa[1]             = 4.31064*binary_x[1]+4.15080*(1-binary_x[1])-0.13549584*binary_x[1]*(1-binary_x[1]);
        Vpppa[1]             = -1.28895*binary_x[1]-1.42744*(1-binary_x[1])+0.11847117*binary_x[1]*(1-binary_x[1]);
        Vpdsa[1]             = -1.73141*binary_x[1]-1.87428*(1-binary_x[1])+0.12099301*binary_x[1]*(1-binary_x[1]);
        Vpdpa[1]             = 2.18886*binary_x[1]+2.52926*(1-binary_x[1])-0.097900682*binary_x[1]*(1-binary_x[1]);
        Vddsa[1]             = -1.58461*binary_x[1]-1.26996*(1-binary_x[1])+0.032701635*binary_x[1]*(1-binary_x[1]);
        Vddpa[1]             = 2.71793*binary_x[1]+2.50536*(1-binary_x[1])+0.21167118*binary_x[1]*(1-binary_x[1]);
        Vddda[1]             = -0.50509*binary_x[1]-0.85174*(1-binary_x[1])+0.346428*binary_x[1]*(1-binary_x[1]);
//Strain
        Csasc[1]             = 1.258286*binary_x[1]+0.586960*(1-binary_x[1]);
        Cstarastarc[1]       = 2.481447*binary_x[1]+0.486090*(1-binary_x[1]);
        Csastarc[1]          = 4.557774*binary_x[1]+0.770950*(1-binary_x[1]);
        Csapc[1]             = 4.367575*binary_x[1]+0.759790*(1-binary_x[1]);
        Cstarapc[1]          = 3.298598*binary_x[1]+0.810790*(1-binary_x[1]);
        Csadc[1]             = 0.0*binary_x[1]+1.070150*(1-binary_x[1]);
        Cstaradc[1]          = 1.195042*binary_x[1]+1.032560*(1-binary_x[1]);
        Cpapc[1]             = 4.624438*binary_x[1]+2.0*(1-binary_x[1]);
        Cpadc[1]             = 0.0*binary_x[1]+1.613500*(1-binary_x[1]);
        Cdadc[1]             = 0.246999*binary_x[1]+1.262620*(1-binary_x[1]);

//Cation 2 (Ga) and from cation 2 to anion (As)
        Esc[1]               = -0.58193*binary_x[1]-0.24119*(1-binary_x[1])-0.20396862*binary_x[1]*(1-binary_x[1]);
        Epc[1]               = 6.97163*binary_x[1]+6.70776*(1-binary_x[1])+0.26205851*binary_x[1]*(1-binary_x[1]);
        Estarc[1]            = 19.94138*binary_x[1]+22.66352*(1-binary_x[1])-1.4771794*binary_x[1]*(1-binary_x[1]);
        Edc12[1]             = 13.30709*binary_x[1]+12.74846*(1-binary_x[1])+0.34268775*binary_x[1]*(1-binary_x[1]);
	Edc15[1]             = 13.30709*binary_x[1]+12.74846*(1-binary_x[1])+0.34268775*binary_x[1]*(1-binary_x[1]);
        lambdac[1]           = 0.1312*binary_x[1]+0.02179*(1-binary_x[1]);     //(Delta/3)
        Vsssc[1]             = -1.69435*binary_x[1]-1.64508*(1-binary_x[1])-0.026152007*binary_x[1]*(1-binary_x[1]);
        Vstarstarsc[1]       = -4.21045*binary_x[1]-3.67720*(1-binary_x[1])-0.28501808*binary_x[1]*(1-binary_x[1]);
        Vsstarsc[1]          = -2.42674*binary_x[1]-2.20777*(1-binary_x[1])-0.054524137*binary_x[1]*(1-binary_x[1]);
        Vspsc[1]             = 2.80936*binary_x[1]+2.96032*(1-binary_x[1])-0.14879334*binary_x[1]*(1-binary_x[1]);
        Vstarpsc[1]          = 0.93734*binary_x[1]+1.02755*(1-binary_x[1])-0.090199823*binary_x[1]*(1-binary_x[1]);
        Vsdsc[1]             = -2.29309*binary_x[1]-2.32059*(1-binary_x[1])-0.013103513*binary_x[1]*(1-binary_x[1]);
        Vstardsc[1]          = -0.48899*binary_x[1]-0.13324*(1-binary_x[1])+0.35564052*binary_x[1]*(1-binary_x[1]);
        Vppsc[1]             = 4.31064*binary_x[1]+4.15080*(1-binary_x[1])-0.13549584*binary_x[1]*(1-binary_x[1]);
        Vpppc[1]             = -1.28895*binary_x[1]-1.42744*(1-binary_x[1])+0.11847117*binary_x[1]*(1-binary_x[1]);
        Vpdsc[1]             = -1.97842*binary_x[1]-1.88964*(1-binary_x[1])+0.087554021*binary_x[1]*(1-binary_x[1]);
        Vpdpc[1]             = 2.45602*binary_x[1]+2.54913*(1-binary_x[1])-0.0931*binary_x[1]*(1-binary_x[1]);
        Vddsc[1]             = -1.58461*binary_x[1]-1.26996*(1-binary_x[1])+0.032701635*binary_x[1]*(1-binary_x[1]);
        Vddpc[1]             = 2.71793*binary_x[1]+2.50536*(1-binary_x[1])+0.21167118*binary_x[1]*(1-binary_x[1]);
        Vdddc[1]             = -0.50509*binary_x[1]-0.85174*(1-binary_x[1])+0.346428*binary_x[1]*(1-binary_x[1]);
//Strain
        Cscsa[1]             = 1.258286*binary_x[1]+0.586960*(1-binary_x[1]);
        Cstarcstara[1]       = 2.481447*binary_x[1]+0.486090*(1-binary_x[1]);
        Cscstara[1]          = 1.086223*binary_x[1]+0.889210*(1-binary_x[1]);
        Cscpa[1]             = 7.02966*binary_x[1]+1.458910*(1-binary_x[1]);
        Cstarcpa[1]          = 7.029496*binary_x[1]+1.212020*(1-binary_x[1]);
        Cscda[1]             = 0.187036*binary_x[1]+0.580530*(1-binary_x[1]);
        Cstarcda[1]          = 1.769483*binary_x[1]+1.323850*(1-binary_x[1]);
        Cpcpa[1]             = 4.624438*binary_x[1]+2.0*(1-binary_x[1]);
        Cpcda[1]             = 0.0*binary_x[1]+1.5*(1-binary_x[1]);
        Cdcda[1]             = 0.246999*binary_x[1]+1.26262*(1-binary_x[1]);

        eta_sss[1]           = 1.92494*binary_x[1]+2.060010*(1-binary_x[1]);
        eta_sstars[1]        = 0.0608*binary_x[1]+0.0*(1-binary_x[1]);
        eta_starstars[1]     = 0.00081*binary_x[1]+0.212660*(1-binary_x[1]);
        eta_sps[1]           = 1.57003*binary_x[1]+1.384980*(1-binary_x[1]);
        eta_starps[1]        = 1.94937*binary_x[1]+1.399300*(1-binary_x[1]);
        eta_sds[1]           = 1.76566*binary_x[1]+1.898890*(1-binary_x[1]);
        eta_stards[1]        = 2.02387*binary_x[1]+1.785400*(1-binary_x[1]);
        eta_pps[1]           = 2.06151*binary_x[1]+2.684970*(1-binary_x[1]);
        eta_ppp[1]           = 1.60247*binary_x[1]+1.314050*(1-binary_x[1]);
        eta_pds[1]           = 2.38382*binary_x[1]+1.812350*(1-binary_x[1]);
        eta_pdp[1]           = 2.4556*binary_x[1]+2.379640*(1-binary_x[1]);
        eta_dds[1]           = 2.32291*binary_x[1]+1.724430*(1-binary_x[1]);
        eta_ddp[1]           = 1.61589*binary_x[1]+1.972530*(1-binary_x[1]);
        eta_ddd[1]           = 2.3296*binary_x[1]+1.896720*(1-binary_x[1]);

        Eshift[1]            = 27;

        bond_length[1]       = (0.60583*binary_x[1]+0.56533*(1-binary_x[1]))*sqrt(3.0)/4.0;

	//InxAl1-xAs
//Anion (As) and from anion to cation 3 (Al)
        Esa[2]               = -5.50042*binary_x[2]-5.17012*(1-binary_x[2]);
        Epa[2]               = 4.15107*binary_x[2]+4.39708*(1-binary_x[2]);
        Estara[2]            = 19.71059*binary_x[2]+19.80474*(1-binary_x[2]);
        Eda12[2]             = 13.03169*binary_x[2]+13.13880*(1-binary_x[2]);
	Eda15[2]             = 13.03169*binary_x[2]+13.13880*(1-binary_x[2]);
        lambdaa[2]           = 0.17234*binary_x[2]+0.17386*(1-binary_x[2]);     //(Delta/3)
        Vsssa[2]             = -1.69435*binary_x[2]-1.64584*(1-binary_x[2]);
        Vstarstarsa[2]       = -4.21045*binary_x[2]-2.84245*(1-binary_x[2]);
        Vsstarsa[2]          = -1.15987*binary_x[2]-2.78690*(1-binary_x[2]);
        Vspsa[2]             = 2.59823*binary_x[2]+3.02223*(1-binary_x[2]);
        Vstarpsa[2]          = 2.06766*binary_x[2]+1.92174*(1-binary_x[2]);
        Vsdsa[2]             = -2.26837*binary_x[2]-3.03196*(1-binary_x[2]);
        Vstardsa[2]          = -0.89937*binary_x[2]-1.84300*(1-binary_x[2]);
        Vppsa[2]             = 4.31064*binary_x[2]+4.53156*(1-binary_x[2]);
        Vpppa[2]             = -1.28895*binary_x[2]-1.86816*(1-binary_x[2]);
        Vpdsa[2]             = -1.73141*binary_x[2]-2.47345*(1-binary_x[2]);
        Vpdpa[2]             = 2.18886*binary_x[2]+2.52741*(1-binary_x[2]);
        Vddsa[2]             = -1.58461*binary_x[2]-1.97058*(1-binary_x[2]);
        Vddpa[2]             = 2.71793*binary_x[2]+1.67733*(1-binary_x[2]);
        Vddda[2]             = -0.50509*binary_x[2]-1.58868*(1-binary_x[2]);
//Strain
	Csasc[2]             = 1.258286*binary_x[2]+0.586960*(1-binary_x[2]);
        Cstarastarc[2]       = 2.481447*binary_x[2]+0.486090*(1-binary_x[2]);
        Csastarc[2]          = 4.557774*binary_x[2]+0.770950*(1-binary_x[2]);
        Csapc[2]             = 4.367575*binary_x[2]+0.759790*(1-binary_x[2]);
        Cstarapc[2]          = 3.298598*binary_x[2]+0.810790*(1-binary_x[2]);
        Csadc[2]             = 0.0*binary_x[2]+1.070150*(1-binary_x[2]);
        Cstaradc[2]          = 1.195042*binary_x[2]+1.032560*(1-binary_x[2]);
        Cpapc[2]             = 4.624438*binary_x[2]+2.0*(1-binary_x[2]);
        Cpadc[2]             = 0.0*binary_x[2]+1.613500*(1-binary_x[2]);
        Cdadc[2]             = 0.246999*binary_x[2]+1.262620*(1-binary_x[2]);

//Cation 3 (Al) and from cation 3 to anion (As)
        Esc[2]               = -0.58193*binary_x[2]+0.79695*(1-binary_x[2]);
        Epc[2]               = 6.97163*binary_x[2]+6.63291*(1-binary_x[2]);
        Estarc[2]            = 19.94138*binary_x[2]+24.16587*(1-binary_x[2]);
        Edc12[2]             = 13.30709*binary_x[2]+12.92122*(1-binary_x[2]);
	Edc15[2]             = 13.30709*binary_x[2]+12.92122*(1-binary_x[2]);
        lambdac[2]           = 0.1312*binary_x[2]+0.01589*(1-binary_x[2]);     //(Delta/3)
        Vsssc[2]             = -1.69435*binary_x[2]-1.64584*(1-binary_x[2]);
        Vstarstarsc[2]       = -4.21045*binary_x[2]-2.84245*(1-binary_x[2]);
        Vsstarsc[2]          = -2.42674*binary_x[2]-1.88341*(1-binary_x[2]);
        Vspsc[2]             = 2.80936*binary_x[2]+2.95309*(1-binary_x[2]);
        Vstarpsc[2]          = 0.93734*binary_x[2]+1.30469*(1-binary_x[2]);
        Vsdsc[2]             = -2.29309*binary_x[2]-2.64111*(1-binary_x[2]);
        Vstardsc[2]          = -0.48899*binary_x[2]-1.73510*(1-binary_x[2]);
        Vppsc[2]             = 4.31064*binary_x[2]+4.53156*(1-binary_x[2]);
        Vpppc[2]             = -1.28895*binary_x[2]-1.86816*(1-binary_x[2]);
        Vpdsc[2]             = -1.97842*binary_x[2]-1.02836*(1-binary_x[2]);
        Vpdpc[2]             = 2.45602*binary_x[2]+2.86419*(1-binary_x[2]);
        Vddsc[2]             = -1.58461*binary_x[2]-1.97058*(1-binary_x[2]);
        Vddpc[2]             = 2.71793*binary_x[2]+1.67733*(1-binary_x[2]);
        Vdddc[2]             = -0.50509*binary_x[2]-1.58868*(1-binary_x[2]);
//Strain
        Cscsa[2]             = 1.258286*binary_x[2]+0.586960*(1-binary_x[2]);
        Cstarcstara[2]       = 2.481447*binary_x[2]+0.486090*(1-binary_x[2]);
        Cscstara[2]          = 1.086223*binary_x[2]+0.889210*(1-binary_x[2]);
        Cscpa[2]             = 7.02966*binary_x[2]+1.458910*(1-binary_x[2]);
        Cstarcpa[2]          = 7.029496*binary_x[2]+1.212020*(1-binary_x[2]);
        Cscda[2]             = 0.187036*binary_x[2]+0.580530*(1-binary_x[2]);
        Cstarcda[2]          = 1.769483*binary_x[2]+1.323850*(1-binary_x[2]);
        Cpcpa[2]             = 4.624438*binary_x[2]+2.0*(1-binary_x[2]);
        Cpcda[2]             = 0.0*binary_x[2]+1.5*(1-binary_x[2]);
        Cdcda[2]             = 0.246999*binary_x[2]+1.26262*(1-binary_x[2]);

        eta_sss[2]           = 1.92494*binary_x[2]+2.060010*(1-binary_x[2]);
        eta_sstars[2]        = 0.0608*binary_x[2]+0.0*(1-binary_x[2]);
        eta_starstars[2]     = 0.00081*binary_x[2]+0.212660*(1-binary_x[2]);
        eta_sps[2]           = 1.57003*binary_x[2]+1.384980*(1-binary_x[2]);
        eta_starps[2]        = 1.94937*binary_x[2]+1.399300*(1-binary_x[2]);
        eta_sds[2]           = 1.76566*binary_x[2]+1.898890*(1-binary_x[2]);
        eta_stards[2]        = 2.02387*binary_x[2]+1.785400*(1-binary_x[2]);
        eta_pps[2]           = 2.06151*binary_x[2]+2.684970*(1-binary_x[2]);
        eta_ppp[2]           = 1.60247*binary_x[2]+1.314050*(1-binary_x[2]);
        eta_pds[2]           = 2.38382*binary_x[2]+1.812350*(1-binary_x[2]);
        eta_pdp[2]           = 2.4556*binary_x[2]+2.379640*(1-binary_x[2]);
        eta_dds[2]           = 2.32291*binary_x[2]+1.724430*(1-binary_x[2]);
        eta_ddp[2]           = 1.61589*binary_x[2]+1.972530*(1-binary_x[2]);
        eta_ddd[2]           = 2.3296*binary_x[2]+1.896720*(1-binary_x[2]);

        Eshift[2]            = 27;

        bond_length[2]       = (0.60583*binary_x[2]+0.56611*(1-binary_x[2]))*sqrt(3.0)/4.0;

	neighbor_table[0][1] = 11;
        neighbor_table[0][3] = 21;
	neighbor_table[0][5] = 31;
	neighbor_table[2][1] = 11;
        neighbor_table[2][3] = 21;
	neighbor_table[2][5] = 31;
	neighbor_table[4][1] = 11;
        neighbor_table[4][3] = 21;
	neighbor_table[4][5] = 31;
        neighbor_table[1][0] = 12;
        neighbor_table[3][0] = 22;
        neighbor_table[5][0] = 32;
	neighbor_table[1][2] = 12;
        neighbor_table[3][2] = 22;
        neighbor_table[5][2] = 32;
	neighbor_table[1][4] = 12;
        neighbor_table[3][4] = 22;
        neighbor_table[5][4] = 32;

	no_orb[0]            = sp3d5ss; //As
	no_orb[1]            = sp3d5ss; //In
	no_orb[2]            = sp3d5ss; //As
	no_orb[3]            = sp3d5ss; //Ga
	no_orb[4]            = sp3d5ss; //As
	no_orb[5]            = sp3d5ss; //Al

	//atomic_mass[0]       = 74.9216*amu;
	//atomic_mass[1]       = 114.8180*amu;
	//atomic_mass[2]       = 74.9216*amu;
	//atomic_mass[3]       = (114.8180*binary_x[1]+69.7230*(1-binary_x[1]))*amu;
	//atomic_mass[4]       = 74.9216*amu;
	//atomic_mass[5]       = (114.8180*binary_x[2]+26.981539*(1-binary_x[2]))*amu;

	//InAs Parameters
	atomic_mass[0]       = 74.9216*amu;
	atomic_mass[1]       = 114.8180*amu;
	atomic_mass[2]       = 74.9216*amu;
	atomic_mass[3]       = 114.8180*amu;
	atomic_mass[4]       = 74.9216*amu;
	atomic_mass[5]       = 114.8180*amu;

	alpha_ph[0]          = 34.6779;
	beta_ph[0]           = 2.6992;
	kappa_ph[0]          = 0.8682;
	tau_ph[0]            = 1.9896;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 34.6779;
	beta_ph[1]           = 2.6992;
	kappa_ph[1]          = 0.8682;
	tau_ph[1]            = 1.9896;
	gamma_ph[1]          = 0.0;

	//InAs Parameters 
	alpha_ph[2]          = 34.6779;
	beta_ph[2]           = 2.6992;
	kappa_ph[2]          = 0.8682;
	tau_ph[2]            = 1.9896;
	gamma_ph[2]          = 0.0;

	alpha_ph[3]          = 34.6779;
	beta_ph[3]           = 2.6992;
	kappa_ph[3]          = 0.8682;
	tau_ph[3]            = 1.9896;
	gamma_ph[3]          = 0.0;

	alpha_ph[4]          = 34.6779;
	beta_ph[4]           = 2.6992;
	kappa_ph[4]          = 0.8682;
	tau_ph[4]            = 1.9896;
	gamma_ph[4]          = 0.0;

	alpha_ph[5]          = 34.6779;
	beta_ph[5]           = 2.6992;
	kappa_ph[5]          = 0.8682;
	tau_ph[5]            = 1.9896;
	gamma_ph[5]          = 0.0;

        break;

    case 6:  // SiP

        Eg                   = 1.13;
        ECmin                = Eg;
        EVmax                = 0;

        alpha                = 48.5;
        beta                 = 13.8;
        ideal_a0             = 0.543095;

//Si, Si-Si, and Si-P
        Esa[0]               = -2.15168;
        Epa[0]               = 4.22925;
        Estara[0]            = 19.11650;
        Eda12[0]             = 13.78950;
	Eda15[0]             = 13.78950;
        lambdaa[0]           = 0.01989;     //(Delta/3)
        Vsssa[0]             = -1.95933;
        Vstarstarsa[0]       = -4.24135;
        Vsstarsa[0]          = -1.52230;
        Vspsa[0]             = 3.02562;
        Vstarpsa[0]          = 3.15565;
        Vsdsa[0]             = -2.28485;
        Vstardsa[0]          = -0.80993;
        Vppsa[0]             = 4.10364;
        Vpppa[0]             = -1.51801;
        Vpdsa[0]             = -1.35554;
        Vpdpa[0]             = 2.38479;
        Vddsa[0]             = -1.68136;
        Vddpa[0]             = 2.58880;
        Vddda[0]             = -1.81400;
//Strain
        Csasc[0]             = 1.68054;
        Cstarastarc[0]       = 0.77849;
        Csastarc[0]          = 1.78613;
        Csapc[0]             = 0.48057;
        Cstarapc[0]          = 3.59244;
        Csadc[0]             = 0;
        Cstaradc[0]          = 0.34243;
        Cpapc[0]             = 4.07053;
        Cpadc[0]             = 0;
        Cdadc[0]             = 5.30270;

//P, P-Si
        Esc[0]               = -2.15168;
        Epc[0]               = 5.32925;
        Estarc[0]            = 19.11650;
        Edc12[0]             = 14.18950;
	Edc15[0]             = 14.18950;
        lambdac[0]           = 0.01989;     //(Delta/3)
        Vsssc[0]             = -1.95933;
        Vstarstarsc[0]       = -4.24135;
        Vsstarsc[0]          = -1.52230;
        Vspsc[0]             = 3.02562;
        Vstarpsc[0]          = 3.15565;
        Vsdsc[0]             = -2.28485;
        Vstardsc[0]          = -0.80993;
        Vppsc[0]             = 4.10364;
        Vpppc[0]             = -1.51801;
        Vpdsc[0]             = -1.35554;
        Vpdpc[0]             = 2.38479;
        Vddsc[0]             = -1.68136;
        Vddpc[0]             = 2.58880;
        Vdddc[0]             = -1.81400;
//Strain
        Cscsa[0]             = 1.68054;
        Cstarcstara[0]       = 0.77849;
        Cscstara[0]          = 1.78613;
        Cscpa[0]             = 0.48057;
        Cstarcpa[0]          = 3.59244;
        Cscda[0]             = 0;
        Cstarcda[0]          = 0.34243;
        Cpcpa[0]             = 4.07053;
        Cpcda[0]             = 0;
        Cdcda[0]             = 5.30270;

        eta_sss[0]           = 0.562469;
        eta_sstars[0]        = 0.132030;
        eta_starstars[0]     = 0.192369;
        eta_sps[0]           = 2.365484;
        eta_starps[0]        = 0.344918;
        eta_sds[0]           = 2.567199;
        eta_stards[0]        = 1.086006;
        eta_pps[0]           = 0.494354;
        eta_ppp[0]           = 1.843851;
        eta_pds[0]           = 2.236358;
        eta_pdp[0]           = 4.512500;
        eta_dds[0]           = 4.668357;
        eta_ddp[0]           = 2.302375;
        eta_ddd[0]           = 0.923911;

        Eshift[0]            = 27;

        bond_length[0]       = 0.543*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[0][3] = 11;
	neighbor_table[1][0] = 11;
	neighbor_table[1][2] = 11;
	neighbor_table[2][1] = 12;
        neighbor_table[3][0] = 12;

	no_orb[0]            = sp3d5ss; //Si anion
	no_orb[1]            = sp3d5ss; //Si cation
	no_orb[2]            = sp3d5ss; //P anion
	no_orb[3]            = sp3d5ss; //P cation

	atomic_mass[0]       = 28.0855*amu;  //Si 
	atomic_mass[1]       = 28.0855*amu;  //Si
	atomic_mass[2]       = 30.97376*amu; //P
	atomic_mass[3]       = 30.97376*amu; //P

        break;

    case 7: //Si_SiO2
      
        Eg                   = 1.13;
	ECmin                = Eg;
	EVmax                = 0;

	/*
	double hbar          = 1.0546e-34;
	double mo            = 9.1095e-31;
	double e             = 1.6022e-19;
	double a0            = 0.543e-9;

	double Es            = -4.8;
	double me            = 0.44;
	double mhh           = -1.5;
	double dECB          = 4.0;

	double EG            = 8.9;
	double ECB           = 4.28;

	double Vxx           = (EG+dECB)/2.0;
	double Ep            = ECB+dECB-Vxx;
	double Vss           = Es-ECB;
	double Vsp           = sqrt((hbar*hbar/(2.0*mo*me*e)*(4.0/a0)*(4.0/a0)-(Es-ECB)/2.0)*EG);
	double Vxy           = sqrt(Vxx*(Vxx-hbar*hbar/(mo*mhh*e)*(4.0/a0)*(4.0/a0)));

	double Vsss          = Vss/4.0;
	double Vsps          = sqrt(3.0)/4.0*Vsp;
	double Vpps          = (Vxx+2.0*Vxy)/4.0;
	double Vppp          = (Vxx-Vxy)/4.0;
	
	//Es                   = Es-3.15;
	//Ep                   = Ep-3.15;
	*/

 // Strain constants
	alpha                = 48.5;
	beta                 = 13.8;
	ideal_a0             = 0.543095;

//Anion 1 (Si) and from anion 1 to cation 1 (Si)
        Esa[0]               = -2.15168;
        Epa[0]               = 4.22925;
        Estara[0]            = 19.11650;
        Eda12[0]             = 13.78950;
	Eda15[0]             = 13.78950;
        lambdaa[0]           = 0.01989;     //(Delta/3)
        Vsssa[0]             = -1.95933;
        Vstarstarsa[0]       = -4.24135;
        Vsstarsa[0]          = -1.52230;
        Vspsa[0]             = 3.02562;
        Vstarpsa[0]          = 3.15565;
        Vsdsa[0]             = -2.28485;
        Vstardsa[0]          = -0.80993;
        Vppsa[0]             = 4.10364;
        Vpppa[0]             = -1.51801;
        Vpdsa[0]             = -1.35554;
        Vpdpa[0]             = 2.38479;
        Vddsa[0]             = -1.68136;
        Vddpa[0]             = 2.58880;
        Vddda[0]             = -1.81400;
//Strain
        Csasc[0]             = 1.68054;
        Cstarastarc[0]       = 0.77849;
        Csastarc[0]          = 1.78613;
        Csapc[0]             = 0.48057;
        Cstarapc[0]          = 3.59244;
        Csadc[0]             = 0;
        Cstaradc[0]          = 0.34243;
        Cpapc[0]             = 4.07053;
        Cpadc[0]             = 0;
        Cdadc[0]             = 5.30270;

//Cation 1 (Si )and from cation 1 (Si) to anion 1 (Si)
        Esc[0]               = -2.15168;
        Epc[0]               = 4.22925;
        Estarc[0]            = 19.11650;
        Edc12[0]             = 13.78950;
	Edc15[0]             = 13.78950;
        lambdac[0]           = 0.01989;     //(Delta/3)
        Vsssc[0]             = -1.95933;
        Vstarstarsc[0]       = -4.24135;
        Vsstarsc[0]          = -1.52230;
        Vspsc[0]             = 3.02562;
        Vstarpsc[0]          = 3.15565;
        Vsdsc[0]             = -2.28485;
        Vstardsc[0]          = -0.80993;
        Vppsc[0]             = 4.10364;
        Vpppc[0]             = -1.51801;
        Vpdsc[0]             = -1.35554;
        Vpdpc[0]             = 2.38479;
        Vddsc[0]             = -1.68136;
        Vddpc[0]             = 2.58880;
        Vdddc[0]             = -1.81400;
 //Strain
        Cscsa[0]             = 1.68054;
        Cstarcstara[0]       = 0.77849;
        Cscstara[0]          = 1.78613;
        Cscpa[0]             = 0.48057;
        Cstarcpa[0]          = 3.59244;
        Cscda[0]             = 0;
        Cstarcda[0]          = 0.34243;
        Cpcpa[0]             = 4.07053;
        Cpcda[0]             = 0;
        Cdcda[0]             = 5.30270;

        eta_sss[0]           = 0.562469;
        eta_sstars[0]        = 0.132030;
        eta_starstars[0]     = 0.192369;
        eta_sps[0]           = 2.365484;
        eta_starps[0]        = 0.344918;
        eta_sds[0]           = 2.567199;
        eta_stards[0]        = 1.086006;
        eta_pps[0]           = 0.494354;
        eta_ppp[0]           = 1.843851;
        eta_pds[0]           = 2.236358;
        eta_pdp[0]           = 4.512500;
        eta_dds[0]           = 4.668357;
        eta_ddp[0]           = 2.302375;
        eta_ddd[0]           = 0.923911;

        Eshift[0]            = 27;

        bond_length[0]       = 0.543*sqrt(3.0)/4.0;

//Anion 1 (Si) and from anion 1 (Si) to cation 2 (SiO2)
        Esa[1]               = -2.15168;
        Epa[1]               = 4.22925;
        Estara[1]            = 19.11650;
        Eda12[1]             = 13.78950;
	Eda15[1]             = 13.78950;
        lambdaa[1]           = 0.01989;     //(Delta/3)
        Vsssa[1]             = -1.95933;
        Vstarstarsa[1]       = 0.0;
        Vsstarsa[1]          = 0.0;
        Vspsa[1]             = 3.02562;
        Vstarpsa[1]          = 3.15565;
        Vsdsa[1]             = 0.0;
        Vstardsa[1]          = 0.0;
        Vppsa[1]             = 4.10364;
        Vpppa[1]             = -1.51801;
        Vpdsa[1]             = 0.0;
        Vpdpa[1]             = 0.0;
        Vddsa[1]             = 0.0;
        Vddpa[1]             = 0.0;
        Vddda[1]             = 0.0;
 //Strain
        Csasc[1]             = 0.0;
        Cstarastarc[1]       = 0.0;
        Csastarc[1]          = 0.0;
        Csapc[1]             = 0.0;
        Cstarapc[1]          = 0.0;
        Csadc[1]             = 0.0;
        Cstaradc[1]          = 0.0;
        Cpapc[1]             = 0.0;
        Cpadc[1]             = 0.0;
        Cdadc[1]             = 0.0;

 //Cation 2 (SiO2) and from cation 2 (SiO2) to anion 1 (Si)
        Esc[1]               = -4.8;
        Epc[1]               = 1.83;
        Estarc[1]            = 0.0;
        Edc12[1]             = 0.0;
	Edc15[1]             = 0.0;
        lambdac[1]           = 0.0;     //(Delta/3)
        Vsssc[1]             = -1.95933;
        Vstarstarsc[1]       = 0.0;
        Vsstarsc[1]          = -1.52230;
        Vspsc[1]             = 3.02562;
        Vstarpsc[1]          = 0.0;
        Vsdsc[1]             = -2.28485;
        Vstardsc[1]          = 0.0;
        Vppsc[1]             = 4.10364;
        Vpppc[1]             = -1.51801;
        Vpdsc[1]             = -1.35554;
        Vpdpc[1]             = 2.38479;
        Vddsc[1]             = 0.0;
        Vddpc[1]             = 0.0;
        Vdddc[1]             = 0.0;
	//Strain
        Cscsa[1]             = 0.0;
        Cstarcstara[1]       = 0.0;
        Cscstara[1]          = 0.0;
        Cscpa[1]             = 0.0;
        Cstarcpa[1]          = 0.0;
        Cscda[1]             = 0.0;
        Cstarcda[1]          = 0.0;
        Cpcpa[1]             = 0.0;
        Cpcda[1]             = 0.0;
        Cdcda[1]             = 0.0;

        eta_sss[1]           = 0.0;
        eta_sstars[1]        = 0.0;
        eta_starstars[1]     = 0.0;
        eta_sps[1]           = 0.0;
        eta_starps[1]        = 0.0;
        eta_sds[1]           = 0.0;
        eta_stards[1]        = 0.0;
        eta_pps[1]           = 0.0;
        eta_ppp[1]           = 0.0;
        eta_pds[1]           = 0.0;
        eta_pdp[1]           = 0.0;
        eta_dds[1]           = 0.0;
        eta_ddp[1]           = 0.0;
        eta_ddd[1]           = 0.0;

        Eshift[1]            = 27;

        bond_length[1]       = 0.543*sqrt(3.0)/4.0;

 //Anion 2 (SiO2) and from anion 2 (SiO2) to cation 1 (Si)
        Esa[2]               = -4.8;
        Epa[2]               = 1.83;
        Estara[2]            = 0.0;
        Eda12[2]             = 0.0;
	Eda15[2]             = 0.0;
        lambdaa[2]           = 0.0;     //(Delta/3)
        Vsssa[2]             = -1.95933;
        Vstarstarsa[2]       = 0.0;
        Vsstarsa[2]          = -1.52230;
        Vspsa[2]             = 3.02562;
        Vstarpsa[2]          = 0.0;
        Vsdsa[2]             = -2.28485;
        Vstardsa[2]          = 0.0;
        Vppsa[2]             = 4.10364;
        Vpppa[2]             = -1.51801;
        Vpdsa[2]             = -1.35554;
        Vpdpa[2]             = 2.38479;
        Vddsa[2]             = 0.0;
        Vddpa[2]             = 0.0;
        Vddda[2]             = 0.0;
 //Strain
        Csasc[2]             = 0.0;
        Cstarastarc[2]       = 0.0;
        Csastarc[2]          = 0.0;
        Csapc[2]             = 0.0;
        Cstarapc[2]          = 0.0;
        Csadc[2]             = 0.0;
        Cstaradc[2]          = 0.0;
        Cpapc[2]             = 0.0;
        Cpadc[2]             = 0.0;
        Cdadc[2]             = 0.0;

 //Cation 1 (Si) and from cation 1 (Si) to anion 2 (SiO2)
        Esc[2]               = -2.15168;
        Epc[2]               = 4.22925;
        Estarc[2]            = 19.11650;
        Edc12[2]             = 13.78950;
	Edc15[2]             = 13.78950;
        lambdac[2]           = 0.01989;     //(Delta/3)
        Vsssc[2]             = -1.95933;
        Vstarstarsc[2]       = 0.0;
        Vsstarsc[2]          = 0.0;
        Vspsc[2]             = 3.02562;
        Vstarpsc[2]          = 3.15565;
        Vsdsc[2]             = 0.0;
        Vstardsc[2]          = 0.0;
        Vppsc[2]             = 4.10364;
        Vpppc[2]             = -1.51801;
        Vpdsc[2]             = 0.0;
        Vpdpc[2]             = 0.0;
        Vddsc[2]             = 0.0;
        Vddpc[2]             = 0.0;
        Vdddc[2]             = 0.0;
	//Strain
        Cscsa[2]             = 0.0;
        Cstarcstara[2]       = 0.0;
        Cscstara[2]          = 0.0;
        Cscpa[2]             = 0.0;
        Cstarcpa[2]          = 0.0;
        Cscda[2]             = 0.0;
        Cstarcda[2]          = 0.0;
        Cpcpa[2]             = 0.0;
        Cpcda[2]             = 0.0;
        Cdcda[2]             = 0.0;

        eta_sss[2]           = 0.0;
        eta_sstars[2]        = 0.0;
        eta_starstars[2]     = 0.0;
        eta_sps[2]           = 0.0;
        eta_starps[2]        = 0.0;
        eta_sds[2]           = 0.0;
        eta_stards[2]        = 0.0;
        eta_pps[2]           = 0.0;
        eta_ppp[2]           = 0.0;
        eta_pds[2]           = 0.0;
        eta_pdp[2]           = 0.0;
        eta_dds[2]           = 0.0;
        eta_ddp[2]           = 0.0;
        eta_ddd[2]           = 0.0;

        Eshift[2]            = 27;

        bond_length[2]       = 0.543*sqrt(3.0)/4.0;

 //Anion 2 (SiO2) and from anion 2 (SiO2) to cation 2 (SiO2)
        Esa[3]               = -4.8;
        Epa[3]               = 1.83;
        Estara[3]            = 0.0;
        Eda12[3]             = 0.0;
	Eda15[3]             = 0.0;
        lambdaa[3]           = 0.0;     //(Delta/3)
        Vsssa[3]             = -2.27;
        Vstarstarsa[3]       = 0.0;
        Vsstarsa[3]          = 0.0;
        Vspsa[3]             = 3.92651;
        Vstarpsa[3]          = 0.0;
        Vsdsa[3]             = 0.0;
        Vstardsa[3]          = 0.0;
        Vppsa[3]             = 5.46553;
        Vpppa[3]             = -0.314015;
        Vpdsa[3]             = 0.0;
        Vpdpa[3]             = 0.0;
        Vddsa[3]             = 0.0;
        Vddpa[3]             = 0.0;
        Vddda[3]             = 0.0;
 //Strain
        Csasc[3]             = 0.0;
        Cstarastarc[3]       = 0.0;
        Csastarc[3]          = 0.0;
        Csapc[3]             = 0.0;
        Cstarapc[3]          = 0.0;
        Csadc[3]             = 0.0;
        Cstaradc[3]          = 0.0;
        Cpapc[3]             = 0.0;
        Cpadc[3]             = 0.0;
        Cdadc[3]             = 0.0;

 //Cation 2 (SiO2) and from cation 2 (SiO2) to anion 2 (SiO2)
        Esc[3]               = -4.8;
        Epc[3]               = 1.83;
        Estarc[3]            = 0.0;
        Edc12[3]             = 0.0;
	Edc15[3]             = 0.0;
        lambdac[3]           = 0.0;     //(Delta/3)
        Vsssc[3]             = -2.27;
        Vstarstarsc[3]       = 0.0;
        Vsstarsc[3]          = 0.0;
        Vspsc[3]             = 3.92651;
        Vstarpsc[3]          = 0.0;
        Vsdsc[3]             = 0.0;
        Vstardsc[3]          = 0.0;
        Vppsc[3]             = 5.46553;
        Vpppc[3]             = -0.314015;
        Vpdsc[3]             = 0.0;
        Vpdpc[3]             = 0.0;
        Vddsc[3]             = 0.0;
        Vddpc[3]             = 0.0;
        Vdddc[3]             = 0.0;
	//Strain
        Cscsa[3]             = 0.0;
        Cstarcstara[3]       = 0.0;
        Cscstara[3]          = 0.0;
        Cscpa[3]             = 0.0;
        Cstarcpa[3]          = 0.0;
        Cscda[3]             = 0.0;
        Cstarcda[3]          = 0.0;
        Cpcpa[3]             = 0.0;
        Cpcda[3]             = 0.0;
        Cdcda[3]             = 0.0;

        eta_sss[3]           = 0.0;
        eta_sstars[3]        = 0.0;
        eta_starstars[3]     = 0.0;
        eta_sps[3]           = 0.0;
        eta_starps[3]        = 0.0;
        eta_sds[3]           = 0.0;
        eta_stards[3]        = 0.0;
        eta_pps[3]           = 0.0;
        eta_ppp[3]           = 0.0;
        eta_pds[3]           = 0.0;
        eta_pdp[3]           = 0.0;
        eta_dds[3]           = 0.0;
        eta_ddp[3]           = 0.0;
        eta_ddd[3]           = 0.0;

        Eshift[3]            = 27;

        bond_length[3]       = 0.543*sqrt(3.0)/4.0;

	neighbor_table[0][1] = 11;
        neighbor_table[0][3] = 21;
        neighbor_table[1][0] = 12;
        neighbor_table[1][2] = 32;
        neighbor_table[2][1] = 31;
        neighbor_table[2][3] = 41;
        neighbor_table[3][0] = 22;
        neighbor_table[3][2] = 42;

	no_orb[0]            = sp3d5ss; //Si anion
	no_orb[1]            = sp3d5ss; //Si cation
	no_orb[2]            = sp3;     //SiO2 anion
	no_orb[3]            = sp3;     //SiO2 cation

	atomic_mass[0]       = 28.0855*amu;  //Si
	atomic_mass[1]       = 28.0855*amu;  //Si
	atomic_mass[2]       = 20.0278*amu;  //SiO2=(Si+2*O)/3 O:15.9994
	atomic_mass[3]       = 20.0278*amu;  //SiO2=(Si+2*O)/3

	alpha_ph[0]          = 49.4;
	beta_ph[0]           = 4.79;
	kappa_ph[0]          = 6.99;
	tau_ph[0]            = 5.2;
	gamma_ph[0]          = 0.0;	

	alpha_ph[1]          = 49.4;
	beta_ph[1]           = 4.79;
	kappa_ph[1]          = 6.99;
	tau_ph[1]            = 5.2;
	gamma_ph[1]          = 0.0;

	alpha_ph[2]          = 49.4;
	beta_ph[2]           = 4.79;
	kappa_ph[2]          = 6.99;
	tau_ph[2]            = 5.2;
	gamma_ph[2]          = 0.0;

	alpha_ph[3]          = 49.4;
	beta_ph[3]           = 4.79;
	kappa_ph[3]          = 6.99;
	tau_ph[3]            = 5.2;
	gamma_ph[3]          = 0.0;

	alphap_ph[0]         = -2.3241e12/64.0;
	betap_ph[0]          = -1.0388e12/64.0;
	taup_ph[0]           = -2.2447e11/64.0;
	delta1p_ph[0]        = -2.2300e11/64.0;
	delta3p_ph[0]        = 6.3894e11/64.0;

	alphap_ph[1]         = -2.3241e12/64.0;
	betap_ph[1]          = -1.0388e12/64.0;
	taup_ph[1]           = -2.2447e11/64.0;
	delta1p_ph[1]        = -2.2300e11/64.0;
	delta3p_ph[1]        = 6.3894e11/64.0;

	alphap_ph[2]         = -2.3241e12/64.0;
	betap_ph[2]          = -1.0388e12/64.0;
	taup_ph[2]           = -2.2447e11/64.0;
	delta1p_ph[2]        = -2.2300e11/64.0;
	delta3p_ph[2]        = 6.3894e11/64.0;

	alphap_ph[3]         = -2.3241e12/64.0;
	betap_ph[3]          = -1.0388e12/64.0;
	taup_ph[3]           = -2.2447e11/64.0;
	delta1p_ph[3]        = -2.2300e11/64.0;
	delta3p_ph[3]        = 6.3894e11/64.0;

        break;

    case 8:  // Si_sp3ss_VB
      
        Eg                   = 1.13;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 48.5;
        beta                 = 13.8;
        ideal_a0             = 0.543095;

//Anion and from anion to
        Esa[0]               = -3.31789;
        Epa[0]               = 1.67862;
        Estara[0]            = 8.23164;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           = 0.015;     //(Delta/3)
        Vsssa[0]             = -2.3997;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 3.0927;
        Vstarpsa[0]          = 3.1396;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 2.8117;
        Vpppa[0]             = -0.7701;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = -3.31789;
        Epc[0]               = 1.67862;
        Estarc[0]            = 8.23164;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.015;     //(Delta/3)
        Vsssc[0]             = -2.3997;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 3.0927;
        Vstarpsc[0]          = 3.1396;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 2.8117;
        Vpppc[0]             = -0.7701;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.543*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 11;

	no_orb[0]            = sp3ss;   //Si anion
	no_orb[1]            = sp3ss;   //Si cation
	
	atomic_mass[0]       = 28.0855*amu;
	atomic_mass[1]       = 28.0855*amu;

	alpha_ph[0]          = 49.4;
	beta_ph[0]           = 4.79;
	kappa_ph[0]          = 6.99;
	tau_ph[0]            = 5.2;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 49.4;
	beta_ph[1]           = 4.79;
	kappa_ph[1]          = 6.99;
	tau_ph[1]            = 5.2;
	gamma_ph[1]          = 0.0;

	alphap_ph[0]         = -2.3241e12/64.0;
	betap_ph[0]          = -1.0388e12/64.0;
	taup_ph[0]           = -2.2447e11/64.0;
	delta1p_ph[0]        = -2.2300e11/64.0;
	delta3p_ph[0]        = 6.3894e11/64.0;

	alphap_ph[1]         = -2.3241e12/64.0;
	betap_ph[1]          = -1.0388e12/64.0;
	taup_ph[1]           = -2.2447e11/64.0;
	delta1p_ph[1]        = -2.2300e11/64.0;
	delta3p_ph[1]        = 6.3894e11/64.0;

	break;

    case 9:  //InAs_sp3ss or InxGa1-xAs_sp3ss
    
        Eg                   = 0.37;
        ECmin                = 0.59570;
        EVmax                = 0.22430;

        alpha                = 35.18;    // InAs
        beta                 = 5.49;     // InAs
        ideal_a0             = 0.60583;  // InAs

//Anion (AS) and from anion to cation 1 (In)
        Esa[0]               = (-9.57566+EVmax)*binary_x[0]+(1-binary_x[0])*(-3.53284)+binary_x[0]*(1-binary_x[0])*(1.6691228);
        Epa[0]               = (0.02402+EVmax)*binary_x[0]+(1-binary_x[0])*0.27772+binary_x[0]*(1-binary_x[0])*(2.9399971e-02);
        Estara[0]            = (7.44461+EVmax)*binary_x[0]+(1-binary_x[0])*12.33930+binary_x[0]*(1-binary_x[0])*(4.6704);
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           = 0.1272*binary_x[0]+(1-binary_x[0])*0.10901+binary_x[0]*(1-binary_x[0])*(-8.1573390e-04);     //(Delta/3)
        Vsssa[0]             = -1.2671*binary_x[0]+(1-binary_x[0])*(-1.7191)+binary_x[0]*(1-binary_x[0])*(-1.8166900e-01);
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 1.0903*binary_x[0]+(1-binary_x[0])*1.2381+binary_x[0]*(1-binary_x[0])*(-1.4779451e-01);
        Vstarpsa[0]          = 1.6440*binary_x[0]+(1-binary_x[0])*2.7350+binary_x[0]*(1-binary_x[0])*(-1.091);
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 2.5550*binary_x[0]+(1-binary_x[0])*2.8719+binary_x[0]*(1-binary_x[0])*(-3.1679898e-01);
        Vpppa[0]             = -0.9591*binary_x[0]+(1-binary_x[0])*(-0.9351)+binary_x[0]*(1-binary_x[0])*(-2.3999787e-02);
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation 1 (In) and from cation 1 to anion (As)
        Esc[0]               = (-2.21525+EVmax)*binary_x[0]+(1-binary_x[0])*(-8.11499)+binary_x[0]*(1-binary_x[0])*(2.5562752);
        Epc[0]               = (4.62421+EVmax)*binary_x[0]+(1-binary_x[0])*4.57341+binary_x[0]*(1-binary_x[0])*(2.7519993e-01);
        Estarc[0]            = (4.12648+EVmax)*binary_x[0]+(1-binary_x[0])*4.31241+binary_x[0]*(1-binary_x[0])*(3.8399941e-02);
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.1251*binary_x[0]+(1-binary_x[0])*0.04+binary_x[0]*(1-binary_x[0])*(8.5000000e-02);     //(Delta/3)
        Vsssc[0]             = -1.2671*binary_x[0]+(1-binary_x[0])*(-1.7191)+binary_x[0]*(1-binary_x[0])*(-1.8166900e-01);
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 2.6762*binary_x[0]+(1-binary_x[0])*4.8055+binary_x[0]*(1-binary_x[0])*(-1.1767748);
        Vstarpsc[0]          = 1.0632*binary_x[0]+(1-binary_x[0])*2.1752+binary_x[0]*(1-binary_x[0])*(-1.1119995);
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 2.5550*binary_x[0]+(1-binary_x[0])*2.8719+binary_x[0]*(1-binary_x[0])*(-3.1679898e-01);
        Vpppc[0]             = -0.9591*binary_x[0]+(1-binary_x[0])*(-0.9351)+binary_x[0]*(1-binary_x[0])*(-2.3999787e-02);
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27.0;

        bond_length[0]       = (0.60583*binary_x[0]+0.56532*(1-binary_x[0]))*sqrt(3.0)/4.0;
	
        neighbor_table[0][0] = 11;
        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;
        neighbor_table[1][1] = 12;

	no_orb[0]            = sp3ss;   //As anion
	no_orb[1]            = sp3ss;   //InxGa1-x cation
	
	//InAs parameters
	atomic_mass[0]       = 74.9216*amu;  //As
	atomic_mass[1]       = 114.818*amu;  //In

	alpha_ph[0]          = 34.6779;
	beta_ph[0]           = 2.6992;
	kappa_ph[0]          = 0.8682;
	tau_ph[0]            = 1.9896;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 34.6779;
	beta_ph[1]           = 2.6992;
	kappa_ph[1]          = 0.8682;
	tau_ph[1]            = 1.9896;
	gamma_ph[1]          = 0.0;

	break;

    case 10:  // Si_sp3ss_CB
      
        Eg                   = 1.13;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 48.5;
        beta                 = 13.8;
        ideal_a0             = 0.543095;

//Anion and from anion to
        Esa[0]               = -3.6587;
        Epa[0]               = 1.67889;
        Estara[0]            = 3.87567;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           =  0.015;     //(Delta/3)
        Vsssa[0]             = -1.9929;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 3.8428;
        Vstarpsa[0]          = 2.3434;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 12.0859;
        Vpppa[0]             = -5.4071;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = -3.6587;
        Epc[0]               = 1.67889;
        Estarc[0]            = 3.87567;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.015;     //(Delta/3)
        Vsssc[0]             = -1.9929;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 3.8428;
        Vstarpsc[0]          = 2.3434 ;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 12.0859;
        Vpppc[0]             = -5.4071;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.543*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 11;

	no_orb[0]            = sp3ss;   //Si anion
	no_orb[1]            = sp3ss;   //Si cation

	atomic_mass[0]       = 28.0855*amu;
	atomic_mass[1]       = 28.0855*amu;

	alpha_ph[0]          = 49.4;
	beta_ph[0]           = 4.79;
	kappa_ph[0]          = 6.99;
	tau_ph[0]            = 5.2;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 49.4;
	beta_ph[1]           = 4.79;
	kappa_ph[1]          = 6.99;
	tau_ph[1]            = 5.2;
	gamma_ph[1]          = 0.0;

	alphap_ph[0]         = -2.3241e12/64.0;
	betap_ph[0]          = -1.0388e12/64.0;
	taup_ph[0]           = -2.2447e11/64.0;
	delta1p_ph[0]        = -2.2300e11/64.0;
	delta3p_ph[0]        = 6.3894e11/64.0;

	alphap_ph[1]         = -2.3241e12/64.0;
	betap_ph[1]          = -1.0388e12/64.0;
	taup_ph[1]           = -2.2447e11/64.0;
	delta1p_ph[1]        = -2.2300e11/64.0;
	delta3p_ph[1]        = 6.3894e11/64.0;
	
	break;

    case 11:  // Ge_Vogl_sp3ss
      
        Eg                   = 0.664;
        ECmin                = 1.448;
        EVmax                = 0.77;

	// Strain constants
        alpha                = 39.0;
        beta                 = 12.0;
        ideal_a0             = 0.56571;

//Anion and from anion to
        Esa[0]               = -5.8800;
        Epa[0]               = 1.610;
        Estara[0]            = 6.390;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           =  0.10132;     //(Delta/3)
        Vsssa[0]             = -1.6950;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 2.3664;
        Vstarpsa[0]          = 2.2599;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 2.8525;
        Vpppa[0]             = -0.8225;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = -5.8800;
        Epc[0]               = 1.610;
        Estarc[0]            = 6.390;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.10132;     //(Delta/3)
        Vsssc[0]             = -1.6950;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 2.3664;
        Vstarpsc[0]          = 2.2599;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 2.8525;
        Vpppc[0]             = -0.8225;
        
	Vsssc[0]             = -2.0750;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 2.4808;
        Vstarpsc[0]          = 2.3174;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 2.7162;
        Vpppc[0]             = -0.7150;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.565791*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 11;
	
	no_orb[0]            = sp3ss;   //Ge anion
	no_orb[1]            = sp3ss;   //Ge cation

	atomic_mass[0]       = 72.61*amu;
	atomic_mass[1]       = 72.61*amu;

	break ;

    case 12:  // GaAs_Vogl_sp3ss
      
        Eg                   = 1.424;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 41.49;
        beta                 = 8.94;
        ideal_a0             = 0.56532;

//Anion and from anion to
        Esa[0]               = -8.3431;
        Epa[0]               = 1.0414;
        Estara[0]            = 8.5914;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           =  0.17234;     //(Delta/3)
        Vsssa[0]             = -1.6128;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 1.9399;
        Vstarpsa[0]          = 2.0967;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 3.0276;
        Vpppa[0]             = -0.7808;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = -2.6569;
        Epc[0]               = 3.6686;
        Estarc[0]            = 6.3786;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.02179;     //(Delta/3)
        Vsssc[0]             = -1.6128;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 2.5045;
        Vstarpsc[0]          = 2.0818;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 3.0276;
        Vpppc[0]             = -0.7808;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.56532*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;
	
	no_orb[0]            = sp3ss;   //As anion
	no_orb[1]            = sp3ss;   //Ga cation

	atomic_mass[0]       = 74.9216*amu;  //As
	atomic_mass[1]       = 69.7230*amu;  //Ga

	break;

    case 13:  // InAs_Vogl_sp3ss
      
        Eg                   = 0.37;
        ECmin                = 0.59570;
        EVmax                = 0.22430;

	// Strain constants
        alpha                = 35.18;
        beta                 = 5.49;
        ideal_a0             = 0.60583;

//Anion and from anion to
        Esa[0]               = -9.5381;
        Epa[0]               = 0.9099;
        Estara[0]            = 7.4099;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           =  0.17234;     //(Delta/3)
        Vsssa[0]             = -1.4013;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 1.3144;
        Vstarpsa[0]          = 1.4612;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 2.6946;
        Vpppa[0]             = -0.6574;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = -2.7219;
        Epc[0]               = 3.7201;
        Estarc[0]            = 6.7401;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.131200;     //(Delta/3)        
        Vsssc[0]             = -1.4013;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 2.3551;
        Vstarpsc[0]          = 1.6929;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 2.6946;
        Vpppc[0]             = -0.6574;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.60583*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;

	no_orb[0]            = sp3ss;   //As anion
	no_orb[1]            = sp3ss;   //In cation

	atomic_mass[0]       = 74.9216*amu;  //As
	atomic_mass[1]       = 114.818*amu;  //In

	alpha_ph[0]          = 34.6779;
	beta_ph[0]           = 2.6992;
	kappa_ph[0]          = 0.8682;
	tau_ph[0]            = 1.9896;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 34.6779;
	beta_ph[1]           = 2.6992;
	kappa_ph[1]          = 0.8682;
	tau_ph[1]            = 1.9896;
	gamma_ph[1]          = 0.0;

	break;

    case 14:  // InSb

        Eg                   = 0.235;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.64794;

//Anion and from anion to
        Esa[0]               = -4.9527;
        Epa[0]               = 4.0797;
        Estara[0]            = 16.1664;
        Eda12[0]             = 11.2647;
	Eda15[0]             = 11.2647;
        lambdaa[0]           = 0.4495;     //(Delta/3)
        Vsssa[0]             = -1.1290;
        Vstarstarsa[0]       = -3.2248;
        Vsstarsa[0]          = -1.8819;
        Vspsa[0]             = 2.5362;
        Vstarpsa[0]          = 2.7380;
        Vsdsa[0]             = -2.5635;
        Vstardsa[0]          = -0.7371;
        Vppsa[0]             = 4.1830;
        Vpppa[0]             = -1.4688;
        Vpdsa[0]             = -2.1487;
        Vpdpa[0]             = 1.8462;
        Vddsa[0]             = -1.3052;
        Vddpa[0]             = 2.0784;
        Vddda[0]             = -1.4118;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = 0.3389;
        Epc[0]               = 6.4919;
        Estarc[0]            = 16.1664;
        Edc12[0]             = 11.2647;
	Edc15[0]             = 11.2647;
        lambdac[0]           = 0.1230;     //(Delta/3)
        Vsssc[0]             = -1.1290;
        Vstarstarsc[0]       = -3.2248;
        Vsstarsc[0]          = -2.0042;
        Vspsc[0]             = 2.6980;
        Vstarpsc[0]          = 2.3471;
        Vsdsc[0]             = -2.3085;
        Vstardsc[0]          = -0.8144;
        Vppsc[0]             = 4.1830;
        Vpppc[0]             = -1.4688;
        Vpdsc[0]             = -2.1652;
        Vpdpc[0]             = 1.8491;
        Vddsc[0]             = -1.3052;
        Vddpc[0]             = 2.0784;
        Vdddc[0]             = -1.4118;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 2.0;
        eta_stards[0]        = 2.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 2.0;
        eta_pdp[0]           = 2.0;
        eta_dds[0]           = 2.0;
        eta_ddp[0]           = 2.0;
        eta_ddd[0]           = 2.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.64794*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;

	no_orb[0]            = sp3d5ss; //Sb
	no_orb[1]            = sp3d5ss; //In

	atomic_mass[0]       = 121.760*amu;  //Sb
	atomic_mass[1]       = 114.818*amu;  //In
        
        break;

    case 15:  //InSb_sp3ss
    
        Eg                   = 0.169;
        ECmin                = Eg;
        EVmax                = 0.0;

        alpha                = 0.0;      // InSb
        beta                 = 0.0;      // InSb
        ideal_a0             = 0.64794;  // InSb

//Anion (Sb) and from anion to cation (In)
        Esa[0]               = -7.80905;
        Epa[0]               = -0.14734;
        Estara[0]            = 7.43195;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           = 0.85794/3.0;     //(Delta/3)
        Vsssa[0]             = -4.89637/4.0;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 3.33714*sqrt(3.0)/4.0;
        Vstarpsa[0]          = 4.59953*sqrt(3.0)/4.0;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = (0.75260+2*4.48030)/4.0;
        Vpppa[0]             = (0.75260-4.48030)/4.0;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation (In) and from cation to anion (Sb)
        Esc[0]               = -2.83599;
        Epc[0]               = 3.91522;
        Estarc[0]            = 3.54540;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.51000/3.0;     //(Delta/3)
        Vsssc[0]             = -4.89637/4.0;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 5.60426*sqrt(3.0)/4.0;
        Vstarpsc[0]          = 2.53756*sqrt(3.0)/4.0;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = (0.75260+2*4.48030)/4.0;
        Vpppc[0]             = (0.75260-4.48030)/4.0;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27.0;

        bond_length[0]       = 0.64794*sqrt(3.0)/4.0;
	
        neighbor_table[0][0] = 11;
        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;
        neighbor_table[1][1] = 12;

	no_orb[0]            = sp3ss;   //Sb anion
	no_orb[1]            = sp3ss;   //In cation

	/*
	atomic_mass[0]       = 121.760*amu;  //Sb
	atomic_mass[1]       = 114.818*amu;  //In
	*/

	//InAs Phonon Spectrum Parameters

	atomic_mass[0]       = 74.9216*amu;
	atomic_mass[1]       = 114.8180*amu;

	alpha_ph[0]          = 34.6779;
	beta_ph[0]           = 2.6992;
	kappa_ph[0]          = 0.8682;
	tau_ph[0]            = 1.9896;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 34.6779;
	beta_ph[1]           = 2.6992;
	kappa_ph[1]          = 0.8682;
	tau_ph[1]            = 1.9896;
	gamma_ph[1]          = 0.0;
	
	break;

    case 16:  //InAs_GaSb_sp3ss
    
        Eg                   = 0.37;
        ECmin                = 0.59570;
        EVmax                = 0.22430;

        alpha                = 35.18;    // InAs
        beta                 = 5.49;     // InAs
        ideal_a0             = 0.60583;  // InAs

	//Affinity (InAs)    = 4.90 eV
	//Affinity (GaSb)    = 4.06 eV (modified to 4.0 eV so that InAs-GaSb band gap offset = 150 meV)
	//Affinity (GaAs)    = 4.07 eV
	//Affinity (InSb)    = 4.59 eV

	//deltaV (InAs-GaSb) = 0.459 eV
	//EV(GaAs)           = 0.0 eV

//Anion 1 (As) and from anion 1 to cation 1 (In)
        Esa[0]               = -9.57566+0.2243;
        Epa[0]               = 0.02402+0.2243;
        Estara[0]            = 7.44461+0.2243;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           = 0.1272;     //(Delta/3)
        Vsssa[0]             = -1.2671;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 1.0903;
        Vstarpsa[0]          = 1.6440;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 2.5550;
        Vpppa[0]             = -0.9591;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation 1 (In) and from cation 1 to anion 1 (As)
        Esc[0]               = -2.21525+0.2243;
        Epc[0]               = 4.62421+0.2243;
        Estarc[0]            = 4.12648+0.2243;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.1251;     //(Delta/3)
        Vsssc[0]             = -1.2671;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 2.6762;
        Vstarpsc[0]          = 1.0632;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 2.5550;
        Vpppc[0]             = -0.9591;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 0.0;
        eta_sstars[0]        = 0.0;
        eta_starstars[0]     = 0.0;
        eta_sps[0]           = 0.0;
        eta_starps[0]        = 0.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 0.0;
        eta_ppp[0]           = 0.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27.0;

        bond_length[0]       = 0.60583*sqrt(3.0)/4.0;

//Anion 1 (As) and from anion 1 to cation 2 (Ga)
        Esa[1]               = -3.53284;
        Epa[1]               = 0.27772;
        Estara[1]            = 12.33930;
        Eda12[1]             = 0.0;
	Eda15[1]             = 0.0;
        lambdaa[1]           = 0.32703/3.0;     //(Delta/3)
        Vsssa[1]             = -6.87653/4.0;
        Vstarstarsa[1]       = 0.0;
        Vsstarsa[1]          = 0.0;
        Vspsa[1]             = 2.85929*sqrt(3.0)/4.0;
        Vstarpsa[1]          = 6.31619*sqrt(3.0)/4.0;
        Vsdsa[1]             = 0.0;
        Vstardsa[1]          = 0.0;
        Vppsa[1]             = (1.33572+2*5.07596)/4.0;
        Vpppa[1]             = (1.33572-5.07596)/4.0;
        Vpdsa[1]             = 0.0;
        Vpdpa[1]             = 0.0;
        Vddsa[1]             = 0.0;
        Vddpa[1]             = 0.0;
        Vddda[1]             = 0.0;
//Strain
        Csasc[1]             = 0.0;
        Cstarastarc[1]       = 0.0;
        Csastarc[1]          = 0.0;
        Csapc[1]             = 0.0;
        Cstarapc[1]          = 0.0;
        Csadc[1]             = 0.0;
        Cstaradc[1]          = 0.0;
        Cpapc[1]             = 0.0;
        Cpadc[1]             = 0.0;
        Cdadc[1]             = 0.0;

//Cation 2 (Ga) and from cation 2 to anion 1 (As)
        Esc[1]               = -8.11499;
        Epc[1]               = 4.57341;
        Estarc[1]            = 4.31241;
        Edc12[1]             = 0.0;
	Edc15[1]             = 0.0;
        lambdac[1]           = 0.120/3.0;     //(Delta/3)
        Vsssc[1]             = -6.87653/4.0;
        Vstarstarsc[1]       = 0.0;
        Vsstarsc[1]          = 0.0;
        Vspsc[1]             = 11.09774*sqrt(3.0)/4.0;
        Vstarpsc[1]          = 5.02335*sqrt(3.0)/4.0;
        Vsdsc[1]             = 0.0;
        Vstardsc[1]          = 0.0;
        Vppsc[1]             = (1.33572+2*5.07596)/4.0;
        Vpppc[1]             = (1.33572-5.07596)/4.0;
        Vpdsc[1]             = 0.0;
        Vpdpc[1]             = 0.0;
        Vddsc[1]             = 0.0;
        Vddpc[1]             = 0.0;
        Vdddc[1]             = 0.0;
//Strain
        Cscsa[1]             = 0.0;
        Cstarcstara[1]       = 0.0;
        Cscstara[1]          = 0.0;
        Cscpa[1]             = 0.0;
        Cstarcpa[1]          = 0.0;
        Cscda[1]             = 0.0;
        Cstarcda[1]          = 0.0;
        Cpcpa[1]             = 0.0;
        Cpcda[1]             = 0.0;
        Cdcda[1]             = 0.0;

        eta_sss[1]           = 0.0;
        eta_sstars[1]        = 0.0;
        eta_starstars[1]     = 0.0;
        eta_sps[1]           = 0.0;
        eta_starps[1]        = 0.0;
        eta_sds[1]           = 0.0;
        eta_stards[1]        = 0.0;
        eta_pps[1]           = 0.0;
        eta_ppp[1]           = 0.0;
        eta_pds[1]           = 0.0;
        eta_pdp[1]           = 0.0;
        eta_dds[1]           = 0.0;
        eta_ddp[1]           = 0.0;
        eta_ddd[1]           = 0.0;

        Eshift[1]            = 27.0;

        bond_length[1]       = 0.5666*sqrt(3.0)/4.0;

//Anion 2 (Sb) and from anion 2 to cation 1 (In)
        Esa[2]               = -7.80905+0.7353;
        Epa[2]               = -0.14734+0.7353;
        Estara[2]            = 7.43195+0.7353;
        Eda12[2]             = 0.0;
	Eda15[2]             = 0.0;
        lambdaa[2]           = 0.85794/3.0;     //(Delta/3)
        Vsssa[2]             = -4.89637/4.0;
        Vstarstarsa[2]       = 0.0;
        Vsstarsa[2]          = 0.0;
        Vspsa[2]             = 3.33714*sqrt(3.0)/4.0;
        Vstarpsa[2]          = 4.59953*sqrt(3.0)/4.0;
        Vsdsa[2]             = 0.0;
        Vstardsa[2]          = 0.0;
        Vppsa[2]             = (0.75260+2*4.48030)/4.0;
        Vpppa[2]             = (0.75260-4.48030)/4.0;
        Vpdsa[2]             = 0.0;
        Vpdpa[2]             = 0.0;
        Vddsa[2]             = 0.0;
        Vddpa[2]             = 0.0;
        Vddda[2]             = 0.0;
//Strain
        Csasc[2]             = 0.0;
        Cstarastarc[2]       = 0.0;
        Csastarc[2]          = 0.0;
        Csapc[2]             = 0.0;
        Cstarapc[2]          = 0.0;
        Csadc[2]             = 0.0;
        Cstaradc[2]          = 0.0;
        Cpapc[2]             = 0.0;
        Cpadc[2]             = 0.0;
        Cdadc[2]             = 0.0;

//Cation 1 (In) and from cation 1 to anion 2 (Sb)
        Esc[2]               = -2.83599+0.7353;
        Epc[2]               = 3.91522+0.7353;
        Estarc[2]            = 3.54540+0.7353;
        Edc12[2]             = 0.0;
	Edc15[2]             = 0.0;
        lambdac[2]           = 0.51000/3.0;     //(Delta/3)
        Vsssc[2]             = -4.89637/4.0;
        Vstarstarsc[2]       = 0.0;
        Vsstarsc[2]          = 0.0;
        Vspsc[2]             = 5.60426*sqrt(3.0)/4.0;
        Vstarpsc[2]          = 2.53756*sqrt(3.0)/4.0;
        Vsdsc[2]             = 0.0;
        Vstardsc[2]          = 0.0;
        Vppsc[2]             = (0.75260+2*4.48030)/4.0;
        Vpppc[2]             = (0.75260-4.48030)/4.0;
        Vpdsc[2]             = 0.0;
        Vpdpc[2]             = 0.0;
        Vddsc[2]             = 0.0;
        Vddpc[2]             = 0.0;
        Vdddc[2]             = 0.0;
//Strain
        Cscsa[2]             = 0.0;
        Cstarcstara[2]       = 0.0;
        Cscstara[2]          = 0.0;
        Cscpa[2]             = 0.0;
        Cstarcpa[2]          = 0.0;
        Cscda[2]             = 0.0;
        Cstarcda[2]          = 0.0;
        Cpcpa[2]             = 0.0;
        Cpcda[2]             = 0.0;
        Cdcda[2]             = 0.0;

        eta_sss[2]           = 0.0;
        eta_sstars[2]        = 0.0;
        eta_starstars[2]     = 0.0;
        eta_sps[2]           = 0.0;
        eta_starps[2]        = 0.0;
        eta_sds[2]           = 0.0;
        eta_stards[2]        = 0.0;
        eta_pps[2]           = 0.0;
        eta_ppp[2]           = 0.0;
        eta_pds[2]           = 0.0;
        eta_pdp[2]           = 0.0;
        eta_dds[2]           = 0.0;
        eta_ddp[2]           = 0.0;
        eta_ddd[2]           = 0.0;

        Eshift[2]            = 27.0;

        bond_length[2]       = 0.64794*sqrt(3.0)/4.0;
	
//Anion 2 (Sb) and from anion 2 to cation 2 (Ga)
        Esa[3]               = -7.16208+0.7443;
        Epa[3]               = -0.17071+0.7443;
        Estara[3]            = 7.32190+0.7443;
        Eda12[3]             = 0.0;
	Eda15[3]             = 0.0;
        lambdaa[3]           = 0.75773/3.0;     //(Delta/3)
        Vsssa[3]             = -6.60955/4.0;
        Vstarstarsa[3]       = 0.0;
        Vsstarsa[3]          = 0.0;
        Vspsa[3]             = 3.00325*sqrt(3.0)/4.0;
        Vstarpsa[3]          = 4.69778*sqrt(3.0)/4.0;
        Vsdsa[3]             = 0.0;
        Vstardsa[3]          = 0.0;
        Vppsa[3]             = (0.58073+2*4.76520)/4.0;
        Vpppa[3]             = (0.58073-4.76520)/4.0;
        Vpdsa[3]             = 0.0;
        Vpdpa[3]             = 0.0;
        Vddsa[3]             = 0.0;
        Vddpa[3]             = 0.0;
        Vddda[3]             = 0.0;
//Strain
        Csasc[3]             = 0.0;
        Cstarastarc[3]       = 0.0;
        Csastarc[3]          = 0.0;
        Csapc[3]             = 0.0;
        Cstarapc[3]          = 0.0;
        Csadc[3]             = 0.0;
        Cstaradc[3]          = 0.0;
        Cpapc[3]             = 0.0;
        Cpadc[3]             = 0.0;
        Cdadc[3]             = 0.0;

//Cation 2 (Ga) and from cation 2 to anion 2 (Sb)
        Esc[3]               = -4.77036+0.7443;
        Epc[3]               = 4.06643+0.7443;
        Estarc[3]            = 3.12330+0.7443;
        Edc12[3]             = 0.0;
	Edc15[3]             = 0.0;
        lambdac[3]           = 0.15778/3.0;     //(Delta/3)
        Vsssc[3]             = -6.60955/4.0;
        Vstarstarsc[3]       = 0.0;
        Vsstarsc[3]          = 0.0;
        Vspsc[3]             = 7.78033*sqrt(3.0)/4.0;
        Vstarpsc[3]          = 4.09285*sqrt(3.0)/4.0;
        Vsdsc[3]             = 0.0;
        Vstardsc[3]          = 0.0;
        Vppsc[3]             = (0.58073+2*4.76520)/4.0;
        Vpppc[3]             = (0.58073-4.76520)/4.0;
        Vpdsc[3]             = 0.0;
        Vpdpc[3]             = 0.0;
        Vddsc[3]             = 0.0;
        Vddpc[3]             = 0.0;
        Vdddc[3]             = 0.0;
//Strain
        Cscsa[3]             = 0.0;
        Cstarcstara[3]       = 0.0;
        Cscstara[3]          = 0.0;
        Cscpa[3]             = 0.0;
        Cstarcpa[3]          = 0.0;
        Cscda[3]             = 0.0;
        Cstarcda[3]          = 0.0;
        Cpcpa[3]             = 0.0;
        Cpcda[3]             = 0.0;
        Cdcda[3]             = 0.0;

        eta_sss[3]           = 0.0;
        eta_sstars[3]        = 0.0;
        eta_starstars[3]     = 0.0;
        eta_sps[3]           = 0.0;
        eta_starps[3]        = 0.0;
        eta_sds[3]           = 0.0;
        eta_stards[3]        = 0.0;
        eta_pps[3]           = 0.0;
        eta_ppp[3]           = 0.0;
        eta_pds[3]           = 0.0;
        eta_pdp[3]           = 0.0;
        eta_dds[3]           = 0.0;
        eta_ddp[3]           = 0.0;
        eta_ddd[3]           = 0.0;

        Eshift[3]            = 27.0;

        bond_length[3]       = 0.60959*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[0][3] = 21;
        neighbor_table[1][0] = 12;
        neighbor_table[1][2] = 32;
	neighbor_table[2][1] = 31;
        neighbor_table[2][3] = 41;
        neighbor_table[3][0] = 22;
        neighbor_table[3][2] = 42;

	mid_gap_energy[0][1] = 0.4093;
        mid_gap_energy[0][3] = 0.7120;
        mid_gap_energy[1][0] = 0.4093;
        mid_gap_energy[1][2] = 0.8198;
        mid_gap_energy[2][1] = 0.8198;
        mid_gap_energy[2][3] = 1.1198;
        mid_gap_energy[3][0] = 0.7120;
        mid_gap_energy[3][2] = 1.1198;

	band_gap_table[0][1] = 0.370;
        band_gap_table[0][3] = 1.424;
        band_gap_table[1][0] = 0.370;
        band_gap_table[1][2] = 0.169;
        band_gap_table[2][1] = 0.169;
        band_gap_table[2][3] = 0.751;
        band_gap_table[3][0] = 1.424;
        band_gap_table[3][2] = 0.751;

	no_orb[0]            = sp3ss;   //As anion
	no_orb[1]            = sp3ss;   //In cation
	no_orb[2]            = sp3ss;   //Sb anion
	no_orb[3]            = sp3ss;   //Ga cation
	
	/*
	atomic_mass[0]       = 74.9216*amu;  //As
	atomic_mass[1]       = 114.818*amu;  //In
	atomic_mass[2]       = 121.760*amu;  //Sb
	atomic_mass[3]       = 69.7230*amu;  //Ga
	*/

	//InAs Phonon Spectrum Parameters

	atomic_mass[0]       = 74.9216*amu;
	atomic_mass[1]       = 114.8180*amu;
	atomic_mass[2]       = 74.9216*amu;
	atomic_mass[3]       = 114.8180*amu;

	alpha_ph[0]          = 34.6779;
	beta_ph[0]           = 2.6992;
	kappa_ph[0]          = 0.8682;
	tau_ph[0]            = 1.9896;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 34.6779;
	beta_ph[1]           = 2.6992;
	kappa_ph[1]          = 0.8682;
	tau_ph[1]            = 1.9896;
	gamma_ph[1]          = 0.0;

	alpha_ph[2]          = 34.6779;
	beta_ph[2]           = 2.6992;
	kappa_ph[2]          = 0.8682;
	tau_ph[2]            = 1.9896;
	gamma_ph[2]          = 0.0;

	alpha_ph[3]          = 34.6779;
	beta_ph[3]           = 2.6992;
	kappa_ph[3]          = 0.8682;
	tau_ph[3]            = 1.9896;
	gamma_ph[3]          = 0.0;

	init_gap_tables      = 0;

	break;

    case 17:  // GaSb

        Eg                   = 0.7207;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.609593;

//Anion and from anion to
        Esa[0]               = -5.98375-0.7706;
        Epa[0]               = 4.676303-0.7706;
        Estara[0]            = 15.604254-0.7706;
        Eda12[0]             = 13.690765-0.7706;
	Eda15[0]             = 13.690765-0.7706;
        lambdaa[0]           = 0.419280;     //(Delta/3)
        Vsssa[0]             = -1.292943;
        Vstarstarsa[0]       = -5.190687;
        Vsstarsa[0]          = -2.646522;
        Vspsa[0]             = 3.001255;
        Vstarpsa[0]          = 3.814514;
        Vsdsa[0]             = -1.604796;
        Vstardsa[0]          = -1.992665;
        Vppsa[0]             = 4.595779;
        Vpppa[0]             = -1.424393;
        Vpdsa[0]             = -1.312547;
        Vpdpa[0]             = 1.653695;
        Vddsa[0]             = -0.444050;
        Vddpa[0]             = 2.170088;
        Vddda[0]             = -1.235843;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = -0.499224-0.7706;
        Epc[0]               = 7.147250-0.7706;
        Estarc[0]            = 20.849355-0.7706;
        Edc12[0]             = 9.297686-0.7706;
	Edc15[0]             = 9.297686-0.7706;
        lambdac[0]           = 0.170873;     //(Delta/3)
        Vsssc[0]             = -1.292943;
        Vstarstarsc[0]       = -5.190687;
        Vsstarsc[0]          = -1.056315;
        Vspsc[0]             = 2.517191;
        Vstarpsc[0]          = 2.326126;
        Vsdsc[0]             = -2.778422;
        Vstardsc[0]          = -1.075569;
        Vppsc[0]             = 4.595779;
        Vpppc[0]             = -1.424393;
        Vpdsc[0]             = -2.244886;
        Vpdpc[0]             = 2.411281;
        Vddsc[0]             = -0.444050;
        Vddpc[0]             = 2.170088;
        Vdddc[0]             = -1.235843;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 2.0;
        eta_stards[0]        = 2.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 2.0;
        eta_pdp[0]           = 2.0;
        eta_dds[0]           = 2.0;
        eta_ddp[0]           = 2.0;
        eta_ddd[0]           = 2.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.609593*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;

	no_orb[0]            = sp3d5ss; //Sb
	no_orb[1]            = sp3d5ss; //Ga

	/*
	atomic_mass[0]       = 121.760*amu;  //Sb
	atomic_mass[1]       = 69.7230*amu;  //Ga
	*/

	//InAs Phonon Spectrum Parameters

	atomic_mass[0]       = 74.9216*amu;
	atomic_mass[1]       = 114.8180*amu;

	alpha_ph[0]          = 34.6779;
	beta_ph[0]           = 2.6992;
	kappa_ph[0]          = 0.8682;
	tau_ph[0]            = 1.9896;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 34.6779;
	beta_ph[1]           = 2.6992;
	kappa_ph[1]          = 0.8682;
	tau_ph[1]            = 1.9896;
	gamma_ph[1]          = 0.0;
        
        break;

    case 21:  // GaN_AlGaN

        Eg                   = 3.53;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.0;

//Anion N and from anion N to cation Ga
        Esa[0]               = -8.9379-0.0554;
        Epa[0]               = 2.0626-0.0554;
        Estara[0]            = 28.1824-0.0554;
        Eda12[0]             = 29.4098-0.0554;
	Eda15[0]             = 27.9433-0.0554;
        lambdaa[0]           = 0.0035;     //(Delta/3)
        Vsssa[0]             = -2.5495;
        Vstarstarsa[0]       = -3.9997;
        Vsstarsa[0]          = -2.0860;
        Vspsa[0]             = 3.9210;
        Vstarpsa[0]          = 4.2911;
        Vsdsa[0]             = -3.9072;
        Vstardsa[0]          = -2.0963;
        Vppsa[0]             = 4.7429;
        Vpppa[0]             = -1.4302;
        Vpdsa[0]             = -1.9007;
        Vpdpa[0]             = 2.2761;
        Vddsa[0]             = -1.2016;
        Vddpa[0]             = 6.0757;
        Vddda[0]             = -4.4436;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation Ga and from cation Ga to anion N
        Esc[0]               = 4.7501-0.0554;
        Epc[0]               = 11.4501-0.0554;
        Estarc[0]            = 35.0507-0.0554;
        Edc12[0]             = 26.9898-0.0554;
	Edc15[0]             = 28.4088-0.0554;
        lambdac[0]           = 0.0410;     //(Delta/3)
        Vsssc[0]             = -2.5495;
        Vstarstarsc[0]       = -3.9997;
        Vsstarsc[0]          = -3.7569;
        Vspsc[0]             = 4.0489;
        Vstarpsc[0]          = 2.0861;
        Vsdsc[0]             = -1.2252;
        Vstardsc[0]          = -1.7553;
        Vppsc[0]             = 4.7429;
        Vpppc[0]             = -1.4302;
        Vpdsc[0]             = -1.3286;
        Vpdpc[0]             = 3.2195;
        Vddsc[0]             = -1.2016;
        Vddpc[0]             = 6.0757;
        Vdddc[0]             = -4.4436;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 3.65;
        eta_sstars[0]        = 1.98;
        eta_starstars[0]     = 0.0;
        eta_sps[0]           = 4.91;
        eta_starps[0]        = 1.05;
        eta_sds[0]           = 1.02;
        eta_stards[0]        = 4.32;
        eta_pps[0]           = 4.09;
        eta_ppp[0]           = 3.35;
        eta_pds[0]           = 1.24;
        eta_pdp[0]           = 1.02;
        eta_dds[0]           = 2.0;
        eta_ddp[0]           = 2.0;
        eta_ddd[0]           = 2.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.1954;

 //Anion N and from anion N to cation AlxGa1-x
        Esa[1]               = (-8.9717-0.1512+0.7)*binary_x[1]-(8.9379+0.0554)*(1-binary_x[1]);
        Epa[1]               = (2.1679-0.1512+0.7)*binary_x[1]+(2.0626-0.0554)*(1-binary_x[1]);
        Estara[1]            = (28.8253-0.1512+0.7)*binary_x[1]+(28.1824-0.0554)*(1-binary_x[1]);
        Eda12[1]             = (30.4332-0.1512+0.7)*binary_x[1]+(29.4098-0.0554)*(1-binary_x[1]);
	Eda15[1]             = (29.5599-0.1512+0.7)*binary_x[1]+(27.9433-0.0554)*(1-binary_x[1]);
        lambdaa[1]           = 0.0035*binary_x[1]+0.0035*(1-binary_x[1]);     //(Delta/3)
        Vsssa[1]             = -2.6261*binary_x[1]-2.5495*(1-binary_x[1]);
        Vstarstarsa[1]       = -4.4940*binary_x[1]-3.9997*(1-binary_x[1]);
        Vsstarsa[1]          = -1.8773*binary_x[1]-2.0860*(1-binary_x[1]);
        Vspsa[1]             = 3.7546*binary_x[1]+3.9210*(1-binary_x[1]);
        Vstarpsa[1]          = 4.0680*binary_x[1]+4.2911*(1-binary_x[1]);
        Vsdsa[1]             = -3.8483*binary_x[1]-3.9072*(1-binary_x[1]);
        Vstardsa[1]          = -2.1897*binary_x[1]-2.0963*(1-binary_x[1]);
        Vppsa[1]             = 4.3117*binary_x[1]+4.7429*(1-binary_x[1]);
        Vpppa[1]             = -1.2842*binary_x[1]-1.4302*(1-binary_x[1]);
        Vpdsa[1]             = -2.1306*binary_x[1]-1.9007*(1-binary_x[1]);
        Vpdpa[1]             = 2.5342*binary_x[1]+2.2761*(1-binary_x[1]);
        Vddsa[1]             = -1.0791*binary_x[1]-1.2016*(1-binary_x[1]);
        Vddpa[1]             = 6.4201*binary_x[1]+6.0757*(1-binary_x[1]);
        Vddda[1]             = -4.5012*binary_x[1]-4.4436*(1-binary_x[1]);
//Strain
        Csasc[1]             = 0.0;
        Cstarastarc[1]       = 0.0;
        Csastarc[1]          = 0.0;
        Csapc[1]             = 0.0;
        Cstarapc[1]          = 0.0;
        Csadc[1]             = 0.0;
        Cstaradc[1]          = 0.0;
        Cpapc[1]             = 0.0;
        Cpadc[1]             = 0.0;
        Cdadc[1]             = 0.0;

//Cation AlxGa1-x and from cation Ga to anion N
        Esc[1]               = (6.0728-0.1512+0.7)*binary_x[1]+(4.7501-0.0554)*(1-binary_x[1]);
        Epc[1]               = (11.5368-0.1512+0.7)*binary_x[1]+(11.4501-0.0554)*(1-binary_x[1]);
        Estarc[1]            = (35.0041-0.1512+0.7)*binary_x[1]+(35.0507-0.0554)*(1-binary_x[1]);
        Edc12[1]             = (28.6862-0.1512+0.7)*binary_x[1]+(26.9898-0.0554)*(1-binary_x[1]);
	Edc15[1]             = (29.9039-0.1512+0.7)*binary_x[1]+(28.4088-0.0554)*(1-binary_x[1]);
        lambdac[1]           = 0.0070*binary_x[1]+0.0410*(1-binary_x[1]);     //(Delta/3)
        Vsssc[1]             = -2.6261*binary_x[1]-2.5495*(1-binary_x[1]);
        Vstarstarsc[1]       = -4.4940*binary_x[1] -3.9997*(1-binary_x[1]);
        Vsstarsc[1]          = -2.9127*binary_x[1]-3.7569*(1-binary_x[1]);
        Vspsc[1]             = 3.9269*binary_x[1]+4.0489*(1-binary_x[1]);
        Vstarpsc[1]          = 1.8328*binary_x[1]+2.0861*(1-binary_x[1]);
        Vsdsc[1]             = -1.0331*binary_x[1]-1.2252*(1-binary_x[1]);
        Vstardsc[1]          = -1.5559*binary_x[1]-1.7553*(1-binary_x[1]);
        Vppsc[1]             = 4.3117*binary_x[1]+4.7429*(1-binary_x[1]);
        Vpppc[1]             = -1.2842*binary_x[1]-1.4302*(1-binary_x[1]);
        Vpdsc[1]             = -1.2470*binary_x[1]-1.3286*(1-binary_x[1]);
        Vpdpc[1]             = 3.1629*binary_x[1]+3.2195*(1-binary_x[1]);
        Vddsc[1]             = -1.0791*binary_x[1]-1.2016*(1-binary_x[1]);
        Vddpc[1]             = 6.4201*binary_x[1]+6.0757*(1-binary_x[1]);
        Vdddc[1]             = -4.5012*binary_x[1]-4.4436*(1-binary_x[1]);
//Strain
        Cscsa[1]             = 0.0;
        Cstarcstara[1]       = 0.0;
        Cscstara[1]          = 0.0;
        Cscpa[1]             = 0.0;
        Cstarcpa[1]          = 0.0;
        Cscda[1]             = 0.0;
        Cstarcda[1]          = 0.0;
        Cpcpa[1]             = 0.0;
        Cpcda[1]             = 0.0;
        Cdcda[1]             = 0.0;

        eta_sss[1]           = 3.55*binary_x[1]+3.65*(1-binary_x[1]);
        eta_sstars[1]        = 1.05*binary_x[1]+1.98*(1-binary_x[1]);
        eta_starstars[1]     = 0.0;
        eta_sps[1]           = 4.71*binary_x[1]+4.91*(1-binary_x[1]);
        eta_starps[1]        = 1.36*binary_x[1]+1.05*(1-binary_x[1]);
        eta_sds[1]           = 1.12*binary_x[1]+1.02*(1-binary_x[1]);
        eta_stards[1]        = 4.70*binary_x[1]+4.32*(1-binary_x[1]);
        eta_pps[1]           = 3.56*binary_x[1]+4.09*(1-binary_x[1]);
        eta_ppp[1]           = 2.25*binary_x[1]+3.35*(1-binary_x[1]);
        eta_pds[1]           = 2.02*binary_x[1]+1.24*(1-binary_x[1]);
        eta_pdp[1]           = 1.05*binary_x[1]+1.02*(1-binary_x[1]);
        eta_dds[1]           = 2.0;
        eta_ddp[1]           = 2.0;
        eta_ddd[1]           = 2.0;

        Eshift[1]            = 27;

        bond_length[1]       = 0.1903*binary_x[1]+0.1954*(1-binary_x[1]);

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;
	neighbor_table[2][3] = 11;
        neighbor_table[3][2] = 12;
	neighbor_table[0][3] = 11;
        neighbor_table[3][0] = 12;
	neighbor_table[2][1] = 11;
        neighbor_table[1][2] = 12;
	neighbor_table[4][5] = 21;
        neighbor_table[5][4] = 22;
	neighbor_table[6][7] = 21;
        neighbor_table[7][6] = 22;
	neighbor_table[4][7] = 21;
        neighbor_table[7][4] = 22;
	neighbor_table[6][5] = 21;
        neighbor_table[5][6] = 22;
	neighbor_table[0][5] = 11;
        neighbor_table[5][0] = 12;
	neighbor_table[2][7] = 11;
        neighbor_table[7][2] = 12;
	neighbor_table[0][7] = 11;
        neighbor_table[7][0] = 12;
	neighbor_table[2][5] = 11;
        neighbor_table[5][2] = 12;

	no_orb[0]            = sp3d5ss; //N
	no_orb[1]            = sp3d5ss; //Ga
	no_orb[2]            = sp3d5ss; //N
	no_orb[3]            = sp3d5ss; //Ga
	no_orb[4]            = sp3d5ss; //N
	no_orb[5]            = sp3d5ss; //AlGa
	no_orb[6]            = sp3d5ss; //N
	no_orb[7]            = sp3d5ss; //AlGa

	atomic_mass[0]       = 14.0067*amu;  //N
	atomic_mass[1]       = 69.7230*amu;  //Ga
	atomic_mass[2]       = 14.0067*amu;  //N
	atomic_mass[3]       = 69.7230*amu;  //Ga
	atomic_mass[4]       = 14.0067*amu;  //N
	atomic_mass[5]       = (26.981539*binary_x[0]+69.7230*(1-binary_x[0]))*amu;  //AlxGa1-x
	atomic_mass[6]       = 14.0067*amu;  //N
	atomic_mass[7]       = (26.981539*binary_x[0]+69.7230*(1-binary_x[0]))*amu;  //AlxGa1-x

        break;

    case 50:  // carbon_pz
      
        Eg                   = 0.1;
        ECmin                = 0.05;
        EVmax                = -0.05;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.142;

//Anion and from anion to
        Esa[0]               = 0.0;
        Epa[0]               = 10*tollim;
        Estara[0]            = 0.0;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           = 0.0;     //(Delta/3)
        Vsssa[0]             = 0.0;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 0.0;
        Vstarpsa[0]          = 0.0;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 0.0;
        Vpppa[0]             = -3.0;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = 0.0;
        Epc[0]               = 10*tollim;
        Estarc[0]            = 0.0;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.0;     //(Delta/3)
        Vsssc[0]             = 0.0;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 0.0;
        Vstarpsc[0]          = 0.0;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 0.0;
        Vpppc[0]             = -3.0;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 0.0;
        eta_sstars[0]        = 0.0;
        eta_starstars[0]     = 0.0;
        eta_sps[0]           = 0.0;
        eta_starps[0]        = 0.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 0.0;
        eta_ppp[0]           = 3.124;  //to ensure q0 = 1/t0*dt0/dr = -eta/r0 = 2.2 1/Ang (Mahan)
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.142;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 11;

	no_orb[0]            = pzorb;   //C anion
	no_orb[1]            = pzorb;   //C cation
	
	atomic_mass[0]       = 12.0107*amu;  //C
	atomic_mass[1]       = 12.0107*amu;  //C

	phiR[0][0]           = 365.0;
	phiR[0][1]           = 88.0;
	phiR[0][2]           = 30.0;
	phiR[0][3]           = -19.2;

	phiT[0][0]           = 245.0;
	phiT[0][1]           = -32.3;;
	phiT[0][2]           = -52.5;
	phiT[0][3]           = 22.9;

	phiZ[0][0]           = 98.2;
	phiZ[0][1]           = -4.0;
	phiZ[0][2]           = 1.5;
	phiZ[0][3]           = -5.8;

	phiR[1][0]           = 365.0;
	phiR[1][1]           = 88.0;
	phiR[1][2]           = 30.0;
	phiR[1][3]           = -19.2;

	phiT[1][0]           = 245.0;
	phiT[1][1]           = -32.3;;
	phiT[1][2]           = -52.5;
	phiT[1][3]           = 22.9;

	phiZ[1][0]           = 98.2;
	phiZ[1][1]           = -4.0;
	phiZ[1][2]           = 1.5;
	phiZ[1][3]           = -5.8;

	IPZ                  = 0;
	IS                   = 3;

	break;

    case 51:  // carbon_sp3d5ss
      
        Eg                   = 0.1;
        ECmin                = 0.05;
        EVmax                = -0.05;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.142;
	
//Anion and from anion to
        Esa[0]               = -1.0458-4.18;
        Epa[0]               = 7.0850-4.18;
        Estara[0]            = 38.2661-4.18;
        Eda12[0]             = 27.9267-4.18;
	Eda15[0]             = 27.9267-4.18;
        lambdaa[0]           = 0.0;     //(Delta/3)
        Vsssa[0]             = -4.3882;
        Vstarstarsa[0]       = -2.6737;
        Vsstarsa[0]          = -2.3899;
        Vspsa[0]             = 5.4951;
        Vstarpsa[0]          = 5.1709;
        Vsdsa[0]             = -2.7655;
        Vstardsa[0]          = -2.3034;
        Vppsa[0]             = 7.5480;
        Vpppa[0]             = -2.6363;
        Vpdsa[0]             = -2.1621;
        Vpdpa[0]             = 3.9281;
        Vddsa[0]             = -4.1813;
        Vddpa[0]             = 4.9779;
        Vddda[0]             = -3.9884;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = -1.0458-4.18;
        Epc[0]               = 7.0850-4.18;
        Estarc[0]            = 38.2661-4.18;
        Edc12[0]             = 27.9267-4.18;
	Edc15[0]             = 27.9267-4.18;
        lambdac[0]           = 0.0;     //(Delta/3)
        Vsssc[0]             = -4.3882;
        Vstarstarsc[0]       = -2.6737;
        Vsstarsc[0]          = -2.3899;
        Vspsc[0]             = 5.4951;
        Vstarpsc[0]          = 5.1709;
        Vsdsc[0]             = -2.7655;
        Vstardsc[0]          = -2.3034;
        Vppsc[0]             = 7.5480;
        Vpppc[0]             = -2.6363;
        Vpdsc[0]             = -2.1621;
        Vpdpc[0]             = 3.9281;
        Vddsc[0]             = -4.1813;
        Vddpc[0]             = 4.9779;
        Vdddc[0]             = -3.9884;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 2.0;
        eta_stards[0]        = 2.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 2.0;
        eta_pdp[0]           = 2.0;
        eta_dds[0]           = 2.0;
        eta_ddp[0]           = 2.0;
        eta_ddd[0]           = 2.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.142;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;

	no_orb[0]            = sp3d5ss;   //C anion
	no_orb[1]            = sp3d5ss;   //C cation
	
	atomic_mass[0]       = 12.0107*amu;  //C
	atomic_mass[1]       = 12.0107*amu;  //C

	phiR[0][0]           = 365.0;
	phiR[0][1]           = 88.0;
	phiR[0][2]           = 30.0;
	phiR[0][3]           = -19.2;

	phiT[0][0]           = 245.0;
	phiT[0][1]           = -32.3;;
	phiT[0][2]           = -52.5;
	phiT[0][3]           = 22.9;

	phiZ[0][0]           = 98.2;
	phiZ[0][1]           = -4.0;
	phiZ[0][2]           = 1.5;
	phiZ[0][3]           = -5.8;

	phiR[1][0]           = 365.0;
	phiR[1][1]           = 88.0;
	phiR[1][2]           = 30.0;
	phiR[1][3]           = -19.2;

	phiT[1][0]           = 245.0;
	phiT[1][1]           = -32.3;;
	phiT[1][2]           = -52.5;
	phiT[1][3]           = 22.9;

	phiZ[1][0]           = 98.2;
	phiZ[1][1]           = -4.0;
	phiZ[1][2]           = 1.5;
	phiZ[1][3]           = -5.8;

	break;

    case 52:  // carbon_sp3
      
        Eg                   = 0.1;
        ECmin                = 0.05;
        EVmax                = -0.05;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.142;
	
//Anion and from anion to
	Esa[0]               = -7.3;
        Epa[0]               = 10*tollim;
        Estara[0]            = 0.0;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           = 0.0;     //(Delta/3)
        Vsssa[0]             = -4.3;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 4.98;
        Vstarpsa[0]          = 0.0;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 6.38;
        Vpppa[0]             = -2.66;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
	Esc[0]               = -7.3;
        Epc[0]               = 10*tollim;
        Estarc[0]            = 0.0;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.0;     //(Delta/3)
        Vsssc[0]             = -4.3;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 4.98;
        Vstarpsc[0]          = 0.0;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 6.38;
        Vpppc[0]             = -2.66;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 2.0;
        eta_stards[0]        = 2.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 2.0;
        eta_pdp[0]           = 2.0;
        eta_dds[0]           = 2.0;
        eta_ddp[0]           = 2.0;
        eta_ddd[0]           = 2.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.142;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;

	no_orb[0]            = sp3;   //C anion
	no_orb[1]            = sp3;   //C cation

	atomic_mass[0]       = 12.0107*amu;  //C
	atomic_mass[1]       = 12.0107*amu;  //C

	phiR[0][0]           = 365.0;
	phiR[0][1]           = 88.0;
	phiR[0][2]           = 30.0;
	phiR[0][3]           = -19.2;

	phiT[0][0]           = 245.0;
	phiT[0][1]           = -32.3;;
	phiT[0][2]           = -52.5;
	phiT[0][3]           = 22.9;

	phiZ[0][0]           = 98.2;
	phiZ[0][1]           = -4.0;
	phiZ[0][2]           = 1.5;
	phiZ[0][3]           = -5.8;

	phiR[1][0]           = 365.0;
	phiR[1][1]           = 88.0;
	phiR[1][2]           = 30.0;
	phiR[1][3]           = -19.2;

	phiT[1][0]           = 245.0;
	phiT[1][1]           = -32.3;;
	phiT[1][2]           = -52.5;
	phiT[1][3]           = 22.9;

	phiZ[1][0]           = 98.2;
	phiZ[1][1]           = -4.0;
	phiZ[1][2]           = 1.5;
	phiZ[1][3]           = -5.8;
	
	break;

    case 53:  // carbon_pzdxzdyz
      
        Eg                   = 0.1;
        ECmin                = 0.05;
        EVmax                = -0.05;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.142;
	
//Anion and from anion to
	Esa[0]               = 0.0;
        Epa[0]               = 1.0783;
        Estara[0]            = 0.0;
        Eda12[0]             = 24.0383;
	Eda15[0]             = 24.0383;
        lambdaa[0]           = 0.0;     //(Delta/3)
	Vsssa[0]             = 0.0;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 0.0;
        Vstarpsa[0]          = 0.0;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 0.0;
        Vpppa[0]             = -3.26;
	Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 2.4;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 3.6;
        Vddda[0]             = -7.4;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
	Esc[0]               = 0.0;
        Epc[0]               = 1.0783;
        Estarc[0]            = 0.0;
        Edc12[0]             = 24.0383;
	Edc15[0]             = 24.0383;
        lambdac[0]           = 0.0;     //(Delta/3)
        Vsssc[0]             = 0.0;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 0.0;
        Vstarpsc[0]          = 0.0;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 0.0;
        Vpppc[0]             = -3.26;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 2.4;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 3.6;
        Vdddc[0]             = -7.4;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 2.0;
        eta_stards[0]        = 2.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 2.0;
        eta_pdp[0]           = 2.0;
        eta_dds[0]           = 2.0;
        eta_ddp[0]           = 2.0;
        eta_ddd[0]           = 2.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.142;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 12;

	no_orb[0]            = pzdxzdyz;   //C anion
	no_orb[1]            = pzdxzdyz;   //C cation

	atomic_mass[0]       = 12.0107*amu;  //C
	atomic_mass[1]       = 12.0107*amu;  //C

	phiR[0][0]           = 365.0;
	phiR[0][1]           = 88.0;
	phiR[0][2]           = 30.0;
	phiR[0][3]           = -19.2;

	phiT[0][0]           = 245.0;
	phiT[0][1]           = -32.3;;
	phiT[0][2]           = -52.5;
	phiT[0][3]           = 22.9;

	phiZ[0][0]           = 98.2;
	phiZ[0][1]           = -4.0;
	phiZ[0][2]           = 1.5;
	phiZ[0][3]           = -5.8;

	phiR[1][0]           = 365.0;
	phiR[1][1]           = 88.0;
	phiR[1][2]           = 30.0;
	phiR[1][3]           = -19.2;

	phiT[1][0]           = 245.0;
	phiT[1][1]           = -32.3;;
	phiT[1][2]           = -52.5;
	phiT[1][3]           = 22.9;

	phiZ[1][0]           = 98.2;
	phiZ[1][1]           = -4.0;
	phiZ[1][2]           = 1.5;
	phiZ[1][3]           = -5.8;

	IPZ                  = 0;
	ID2                  = 1;
	ID3                  = 2;
	IS                   = 3;
	IPX                  = 4;
	IPY                  = 5;
	ISS                  = 6;
	ID1                  = 7;
	
	break;

    case 54:  // cnt_pz  (problem because x,y,z basis instead of x,r,phi) 
      
        Eg                   = 0.1;
        ECmin                = 0.05;
        EVmax                = -0.05;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.142;

//Anion and from anion to
        Esa[0]               = 10*tollim;
        Epa[0]               = 0.0;
        Estara[0]            = 0.0;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           = 0.0;     //(Delta/3)
        Vsssa[0]             = -3.0;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 0.0;
        Vstarpsa[0]          = 0.0;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 0.0;
        Vpppa[0]             = 0.0;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = 10*tollim;
        Epc[0]               = 0.0;
        Estarc[0]            = 0.0;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.0;     //(Delta/3)
        Vsssc[0]             = -3.0;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 0.0;
        Vstarpsc[0]          = 0.0;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 0.0;
        Vpppc[0]             = 0.0;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 3.124;  //to ensure q0 = 1/t0*dt0/dr = -eta/r0 = 2.2 1/Ang (Mahan)
        eta_sstars[0]        = 0.0;
        eta_starstars[0]     = 0.0;
        eta_sps[0]           = 0.0;
        eta_starps[0]        = 0.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 0.0;
        eta_ppp[0]           = 0.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.142;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 11;

	no_orb[0]            = pzorb;   //C anion
	no_orb[1]            = pzorb;   //C cation
	
	atomic_mass[0]       = 12.0107*amu;  //C
	atomic_mass[1]       = 12.0107*amu;  //C

	phiR[0][0]           = 365.0;
	phiR[0][1]           = 88.0;
	phiR[0][2]           = 30.0;
	phiR[0][3]           = -19.2;

	phiT[0][0]           = 245.0;
	phiT[0][1]           = -32.3;;
	phiT[0][2]           = -52.5;
	phiT[0][3]           = 22.9;

	phiZ[0][0]           = 98.2;
	phiZ[0][1]           = -4.0;
	phiZ[0][2]           = 1.5;
	phiZ[0][3]           = -5.8;

	phiR[1][0]           = 365.0;
	phiR[1][1]           = 88.0;
	phiR[1][2]           = 30.0;
	phiR[1][3]           = -19.2;

	phiT[1][0]           = 245.0;
	phiT[1][1]           = -32.3;;
	phiT[1][2]           = -52.5;
	phiT[1][3]           = 22.9;

	phiZ[1][0]           = 98.2;
	phiZ[1][1]           = -4.0;
	phiZ[1][2]           = 1.5;
	phiZ[1][3]           = -5.8;

	break;

    case 55:  // carbon_pzdxzdyz_hydro_pz
      
        Eg                   = 0.1;
        ECmin                = 0.05;//0.042;//0.05;
        EVmax                = -0.05;//-0.495;//-0.05;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.142;
	
//Anion C and from anion C to cation C
	Esa[0]               = 0.0;
        Epa[0]               = 1.0783;
        Estara[0]            = 0.0;
        Eda12[0]             = 24.0383;
	Eda15[0]             = 24.0383;
        lambdaa[0]           = 0.0;     //(Delta/3)
	Vsssa[0]             = 0.0;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 0.0;
        Vstarpsa[0]          = 0.0;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 0.0;
        Vpppa[0]             = -3.26;
	Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 2.4;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 3.6;
        Vddda[0]             = -7.4;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation C and from cation C to anion C
	Esc[0]               = 0.0;
        Epc[0]               = 1.0783;
        Estarc[0]            = 0.0;
        Edc12[0]             = 24.0383;
	Edc15[0]             = 24.0383;
        lambdac[0]           = 0.0;     //(Delta/3)
        Vsssc[0]             = 0.0;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 0.0;
        Vstarpsc[0]          = 0.0;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 0.0;
        Vpppc[0]             = -3.26;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 2.4;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 3.6;
        Vdddc[0]             = -7.4;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 0.0;
        eta_sstars[0]        = 0.0;
        eta_starstars[0]     = 0.0;
        eta_sps[0]           = 0.0;
        eta_starps[0]        = 0.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 0.0;
        eta_ppp[0]           = 0.0;
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.142;

//Anion C and from anion C to cation H
	Esa[1]               = 0.0;
        Epa[1]               = 1.0783;
        Estara[1]            = 0.0;
        Eda12[1]             = 24.0383;
	Eda15[1]             = 24.0383;
        lambdaa[1]           = 0.0;     //(Delta/3)
	Vsssa[1]             = 0.0;
        Vstarstarsa[1]       = 0.0;
        Vsstarsa[1]          = 0.0;
        Vspsa[1]             = 0.0;
        Vstarpsa[1]          = 0.0;
        Vsdsa[1]             = 0.0;
        Vstardsa[1]          = 0.0;
        Vppsa[1]             = 0.0;
        Vpppa[1]             = -0.61754;
	Vpdsa[1]             = 0.0;
        Vpdpa[1]             = 3.4117;
        Vddsa[1]             = 0.0;
        Vddpa[1]             = 10.4466;
        Vddda[1]             = -13.9634;
//Strain
        Csasc[1]             = 0.0;
        Cstarastarc[1]       = 0.0;
        Csastarc[1]          = 0.0;
        Csapc[1]             = 0.0;
        Cstarapc[1]          = 0.0;
        Csadc[1]             = 0.0;
        Cstaradc[1]          = 0.0;
        Cpapc[1]             = 0.0;
        Cpadc[1]             = 0.0;
        Cdadc[1]             = 0.0;

//Cation H and from cation H to anion C
	Esc[1]               = 0.0;
        Epc[1]               = 13.0402;
        Estarc[1]            = 0.0;
        Edc12[1]             = 20.902;
	Edc15[1]             = 20.902;
        lambdac[1]           = 0.0;     //(Delta/3)
        Vsssc[1]             = 0.0;
        Vstarstarsc[1]       = 0.0;
        Vsstarsc[1]          = 0.0;
        Vspsc[1]             = 0.0;
        Vstarpsc[1]          = 0.0;
        Vsdsc[1]             = 0.0;
        Vstardsc[1]          = 0.0;
        Vppsc[1]             = 0.0;
        Vpppc[1]             = -0.61754;
        Vpdsc[1]             = 0.0;
        Vpdpc[1]             = 3.4117;
        Vddsc[1]             = 0.0;
        Vddpc[1]             = 10.4466;
        Vdddc[1]             = -13.9634;
//Strain
        Cscsa[1]             = 0.0;
        Cstarcstara[1]       = 0.0;
        Cscstara[1]          = 0.0;
        Cscpa[1]             = 0.0;
        Cstarcpa[1]          = 0.0;
        Cscda[1]             = 0.0;
        Cstarcda[1]          = 0.0;
        Cpcpa[1]             = 0.0;
        Cpcda[1]             = 0.0;
        Cdcda[1]             = 0.0;

        eta_sss[1]           = 0.0;
        eta_sstars[1]        = 0.0;
        eta_starstars[1]     = 0.0;
        eta_sps[1]           = 0.0;
        eta_starps[1]        = 0.0;
        eta_sds[1]           = 0.0;
        eta_stards[1]        = 0.0;
        eta_pps[1]           = 0.0;
        eta_ppp[1]           = 0.0;
        eta_pds[1]           = 0.0;
        eta_pdp[1]           = 0.0;
        eta_dds[1]           = 0.0;
        eta_ddp[1]           = 0.0;
        eta_ddd[1]           = 0.0;

        Eshift[1]            = 27;

        bond_length[1]       = 0.142;

//Anion H and from anion H to cation C
	Esa[2]               = 0.0;
        Epa[2]               = 13.0402;
        Estara[2]            = 0.0;
        Eda12[2]             = 20.902;
	Eda15[2]             = 20.902;
        lambdaa[2]           = 0.0;     //(Delta/3)
	Vsssa[2]             = 0.0;
        Vstarstarsa[2]       = 0.0;
        Vsstarsa[2]          = 0.0;
        Vspsa[2]             = 0.0;
        Vstarpsa[2]          = 0.0;
        Vsdsa[2]             = 0.0;
        Vstardsa[2]          = 0.0;
        Vppsa[2]             = 0.0;
        Vpppa[2]             = -0.61754;
	Vpdsa[2]             = 0.0;
        Vpdpa[2]             = 3.4117;
        Vddsa[2]             = 0.0;
        Vddpa[2]             = 10.4466;
        Vddda[2]             = -13.9634;
//Strain
        Csasc[2]             = 0.0;
        Cstarastarc[2]       = 0.0;
        Csastarc[2]          = 0.0;
        Csapc[2]             = 0.0;
        Cstarapc[2]          = 0.0;
        Csadc[2]             = 0.0;
        Cstaradc[2]          = 0.0;
        Cpapc[2]             = 0.0;
        Cpadc[2]             = 0.0;
        Cdadc[2]             = 0.0;

//Cation C and from cation C to anion H
	Esc[2]               = 0.0;
        Epc[2]               = 1.0783;
        Estarc[2]            = 0.0;
        Edc12[2]             = 24.0383;
	Edc15[2]             = 24.0383;
        lambdac[2]           = 0.0;     //(Delta/3)
        Vsssc[2]             = 0.0;
        Vstarstarsc[2]       = 0.0;
        Vsstarsc[2]          = 0.0;
        Vspsc[2]             = 0.0;
        Vstarpsc[2]          = 0.0;
        Vsdsc[2]             = 0.0;
        Vstardsc[2]          = 0.0;
        Vppsc[2]             = 0.0;
        Vpppc[2]             = -0.61754;
        Vpdsc[2]             = 0.0;
        Vpdpc[2]             = 3.4117;
        Vddsc[2]             = 0.0;
        Vddpc[2]             = 10.4466;
        Vdddc[2]             = -13.9634;
//Strain
        Cscsa[2]             = 0.0;
        Cstarcstara[2]       = 0.0;
        Cscstara[2]          = 0.0;
        Cscpa[2]             = 0.0;
        Cstarcpa[2]          = 0.0;
        Cscda[2]             = 0.0;
        Cstarcda[2]          = 0.0;
        Cpcpa[2]             = 0.0;
        Cpcda[2]             = 0.0;
        Cdcda[2]             = 0.0;

        eta_sss[2]           = 0.0;
        eta_sstars[2]        = 0.0;
        eta_starstars[2]     = 0.0;
        eta_sps[2]           = 0.0;
        eta_starps[2]        = 0.0;
        eta_sds[2]           = 0.0;
        eta_stards[2]        = 0.0;
        eta_pps[2]           = 0.0;
        eta_ppp[2]           = 0.0;
        eta_pds[2]           = 0.0;
        eta_pdp[2]           = 0.0;
        eta_dds[2]           = 0.0;
        eta_ddp[2]           = 0.0;
        eta_ddd[2]           = 0.0;

        Eshift[2]            = 27;

        bond_length[2]       = 0.142;

//Anion H and from anion H to cation H
	Esa[3]               = 0.0;
        Epa[3]               = 13.0402;
        Estara[3]            = 0.0;
        Eda12[3]             = 20.902;
	Eda15[3]             = 20.902;
        lambdaa[3]           = 0.0;     //(Delta/3)
	Vsssa[3]             = 0.0;
        Vstarstarsa[3]       = 0.0;
        Vsstarsa[3]          = 0.0;
        Vspsa[3]             = 0.0;
        Vstarpsa[3]          = 0.0;
        Vsdsa[3]             = 0.0;
        Vstardsa[3]          = 0.0;
        Vppsa[3]             = 0.0;
        Vpppa[3]             = 0.0;
	Vpdsa[3]             = 0.0;
        Vpdpa[3]             = 0.0;
        Vddsa[3]             = 0.0;
        Vddpa[3]             = 0.0;
        Vddda[3]             = 0.0;
//Strain
        Csasc[3]             = 0.0;
        Cstarastarc[3]       = 0.0;
        Csastarc[3]          = 0.0;
        Csapc[3]             = 0.0;
        Cstarapc[3]          = 0.0;
        Csadc[3]             = 0.0;
        Cstaradc[3]          = 0.0;
        Cpapc[3]             = 0.0;
        Cpadc[3]             = 0.0;
        Cdadc[3]             = 0.0;

//Cation H and from cation H to anion H
	Esc[3]               = 0.0;
        Epc[3]               = 13.0402;
        Estarc[3]            = 0.0;
        Edc12[3]             = 20.902;
	Edc15[3]             = 20.902;
        lambdac[3]           = 0.0;     //(Delta/3)
        Vsssc[3]             = 0.0;
        Vstarstarsc[3]       = 0.0;
        Vsstarsc[3]          = 0.0;
        Vspsc[3]             = 0.0;
        Vstarpsc[3]          = 0.0;
        Vsdsc[3]             = 0.0;
        Vstardsc[3]          = 0.0;
        Vppsc[3]             = 0.0;
        Vpppc[3]             = 0.0;
        Vpdsc[3]             = 0.0;
        Vpdpc[3]             = 0.0;
        Vddsc[3]             = 0.0;
        Vddpc[3]             = 0.0;
        Vdddc[3]             = 0.0;
//Strain
        Cscsa[3]             = 0.0;
        Cstarcstara[3]       = 0.0;
        Cscstara[3]          = 0.0;
        Cscpa[3]             = 0.0;
        Cstarcpa[3]          = 0.0;
        Cscda[3]             = 0.0;
        Cstarcda[3]          = 0.0;
        Cpcpa[3]             = 0.0;
        Cpcda[3]             = 0.0;
        Cdcda[3]             = 0.0;

        eta_sss[3]           = 0.0;
        eta_sstars[3]        = 0.0;
        eta_starstars[3]     = 0.0;
        eta_sps[3]           = 0.0;
        eta_starps[3]        = 0.0;
        eta_sds[3]           = 0.0;
        eta_stards[3]        = 0.0;
        eta_pps[3]           = 0.0;
        eta_ppp[3]           = 0.0;
        eta_pds[3]           = 0.0;
        eta_pdp[3]           = 0.0;
        eta_dds[3]           = 0.0;
        eta_ddp[3]           = 0.0;
        eta_ddd[3]           = 0.0;

        Eshift[3]            = 27;

        bond_length[3]       = 0.142;

 	neighbor_table[0][1] = 11;
        neighbor_table[0][3] = 21;
        neighbor_table[1][0] = 12;
        neighbor_table[1][2] = 32;
        neighbor_table[2][1] = 31;
        neighbor_table[2][3] = 41;
        neighbor_table[3][0] = 22;
        neighbor_table[3][2] = 42;

	no_orb[0]            = pzdxzdyz;   //C anion
	no_orb[1]            = pzdxzdyz;   //C cation
	no_orb[2]            = pzdxzdyz;   //H anion
	no_orb[3]            = pzdxzdyz;   //H cation

	atomic_mass[0]       = 12.0107*amu;  //C
	atomic_mass[1]       = 12.0107*amu;  //C
	atomic_mass[2]       = 1.0*amu;      //H
	atomic_mass[3]       = 1.0*amu;      //H

	phiR[0][0]           = 365.0;
	phiR[0][1]           = 88.0;
	phiR[0][2]           = 30.0;
	phiR[0][3]           = -19.2;

	phiT[0][0]           = 245.0;
	phiT[0][1]           = -32.3;;
	phiT[0][2]           = -52.5;
	phiT[0][3]           = 22.9;

	phiZ[0][0]           = 98.2;
	phiZ[0][1]           = -4.0;
	phiZ[0][2]           = 1.5;
	phiZ[0][3]           = -5.8;

	phiR[1][0]           = 365.0;
	phiR[1][1]           = 88.0;
	phiR[1][2]           = 30.0;
	phiR[1][3]           = -19.2;

	phiT[1][0]           = 245.0;
	phiT[1][1]           = -32.3;;
	phiT[1][2]           = -52.5;
	phiT[1][3]           = 22.9;

	phiZ[1][0]           = 98.2;
	phiZ[1][1]           = -4.0;
	phiZ[1][2]           = 1.5;
	phiZ[1][3]           = -5.8;

	phiR[2][0]           = 365.0;
	phiR[2][1]           = 88.0;
	phiR[2][2]           = 30.0;
	phiR[2][3]           = -19.2;

	phiT[2][0]           = 245.0;
	phiT[2][1]           = -32.3;;
	phiT[2][2]           = -52.5;
	phiT[2][3]           = 22.9;

	phiZ[2][0]           = 98.2;
	phiZ[2][1]           = -4.0;
	phiZ[2][2]           = 1.5;
	phiZ[2][3]           = -5.8;

	phiR[3][0]           = 365.0;
	phiR[3][1]           = 88.0;
	phiR[3][2]           = 30.0;
	phiR[3][3]           = -19.2;

	phiT[3][0]           = 245.0;
	phiT[3][1]           = -32.3;;
	phiT[3][2]           = -52.5;
	phiT[3][3]           = 22.9;

	phiZ[3][0]           = 98.2;
	phiZ[3][1]           = -4.0;
	phiZ[3][2]           = 1.5;
	phiZ[3][3]           = -5.8;

	IPZ                  = 0;
	ID2                  = 1;
	ID3                  = 2;
	IS                   = 3;
	IPX                  = 4;
	IPY                  = 5;
	ISS                  = 6;
	ID1                  = 7;
	
	break;

    case 56:  // carbon_py
      
        Eg                   = 0.1;
        ECmin                = 0.05;
        EVmax                = -0.05;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.142;

//Anion and from anion to
        Esa[0]               = 0.0;
        Epa[0]               = 10*tollim;
        Estara[0]            = 0.0;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           = 0.0;     //(Delta/3)
        Vsssa[0]             = 0.0;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 0.0;
        Vstarpsa[0]          = 0.0;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 0.0;
        Vpppa[0]             = -3.0;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = 0.0;
        Epc[0]               = 10*tollim;
        Estarc[0]            = 0.0;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.0;     //(Delta/3)
        Vsssc[0]             = 0.0;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 0.0;
        Vstarpsc[0]          = 0.0;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 0.0;
        Vpppc[0]             = -3.0;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 0.0;
        eta_sstars[0]        = 0.0;
        eta_starstars[0]     = 0.0;
        eta_sps[0]           = 0.0;
        eta_starps[0]        = 0.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 0.0;
        eta_ppp[0]           = 3.124;  //to ensure q0 = 1/t0*dt0/dr = -eta/r0 = 2.2 1/Ang (Mahan)
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.142;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 11;

	no_orb[0]            = pzorb;   //C anion
	no_orb[1]            = pzorb;   //C cation
	
	atomic_mass[0]       = 12.0107*amu;  //C
	atomic_mass[1]       = 12.0107*amu;  //C

	phiR[0][0]           = 365.0;
	phiR[0][1]           = 88.0;
	phiR[0][2]           = 30.0;
	phiR[0][3]           = -19.2;

	phiT[0][0]           = 245.0;
	phiT[0][1]           = -32.3;;
	phiT[0][2]           = -52.5;
	phiT[0][3]           = 22.9;

	phiZ[0][0]           = 98.2;
	phiZ[0][1]           = -4.0;
	phiZ[0][2]           = 1.5;
	phiZ[0][3]           = -5.8;

	phiR[1][0]           = 365.0;
	phiR[1][1]           = 88.0;
	phiR[1][2]           = 30.0;
	phiR[1][3]           = -19.2;

	phiT[1][0]           = 245.0;
	phiT[1][1]           = -32.3;;
	phiT[1][2]           = -52.5;
	phiT[1][3]           = 22.9;

	phiZ[1][0]           = 98.2;
	phiZ[1][1]           = -4.0;
	phiZ[1][2]           = 1.5;
	phiZ[1][3]           = -5.8;

	IPY                  = 0;
	IS                   = 3;

	break;

    case 57:  // carbon_sigma_pi
      
        Eg                   = 0.1;
        ECmin                = 0.05;//0.042;//0.05;
        EVmax                = -0.05;//-0.495;//-0.05;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.142;
	//sigma              = 0.165; //Eps_yy=-sigma*Eps_xx
	
//Anion C and from anion C to cation C
	Esa[0]               = -4.0-0.1274;
	Epa1[0]              = 6.0-0.1274;
	Epa2[0]              = 6.0-0.1274;
        Epa3[0]              = 1.2057-0.1274;
        Estara[0]            = 0.0;
	Eda1[0]              = 30.5-0.1274;
        Eda2[0]              = 24.1657-0.1274;
	Eda3[0]              = 24.1657-0.1274;
	Eda4[0]              = 30.5-0.1274;
	Eda5[0]              = 30.5-0.1274;
        lambdaa[0]           = 0.0;     //(Delta/3)
	Vsssa[0]             = -5.1;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 6.0;
        Vstarpsa[0]          = 0.0;
        Vsdsa[0]             = -1.7;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 8.2;
        Vpppa[0]             = -3.26;
	Vpdsa[0]             = -4.0;
        Vpdpa[0]             = 2.4;
        Vddsa[0]             = -2.5;
        Vddpa[0]             = 3.6;
        Vddda[0]             = -7.4;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation C and from cation C to anion C
	Esc[0]               = -4.0-0.1274;
	Epc1[0]              = 6.0-0.1274;
	Epc2[0]              = 6.0-0.1274;
        Epc3[0]              = 1.2057-0.1274;
        Estarc[0]            = 0.0;
	Edc1[0]              = 30.5-0.1274;
        Edc2[0]              = 24.1657-0.1274;
	Edc3[0]              = 24.1657-0.1274;
	Edc4[0]              = 30.5-0.1274;
	Edc5[0]              = 30.5-0.1274;
        lambdac[0]           = 0.0;     //(Delta/3)
        Vsssc[0]             = -5.1;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 6.0;
        Vstarpsc[0]          = 0.0;
        Vsdsc[0]             = -1.7;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 8.2;
        Vpppc[0]             = -3.26;
        Vpdsc[0]             = -4.0;
        Vpdpc[0]             = 2.4;
        Vddsc[0]             = -2.5;
        Vddpc[0]             = 3.6;
        Vdddc[0]             = -7.4;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 1.87;
        eta_sstars[0]        = 0.0;
        eta_starstars[0]     = 0.0;
        eta_sps[0]           = 1.5;
        eta_starps[0]        = 0.0;
        eta_sds[0]           = 0.7;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 0.62;
        eta_ppp[0]           = 2.6;
        eta_pds[0]           = 1.72;
        eta_pdp[0]           = 1.6;
        eta_dds[0]           = 0.3;
        eta_ddp[0]           = 0.3;
        eta_ddd[0]           = 0.3;

        Eshift[0]            = 27;

        bond_length[0]       = 0.142;

        mod_scaling_fcn      = 1;
        mu_scal              = 4.6;

 	neighbor_table[0][1] = 11;
        neighbor_table[0][3] = 11;
        neighbor_table[1][0] = 12;
        neighbor_table[1][2] = 12;
        neighbor_table[2][1] = 11;
        neighbor_table[2][3] = 11;
        neighbor_table[3][0] = 12;
        neighbor_table[3][2] = 12;

	no_orb[0]            = sp3d5;   //C anion
	no_orb[1]            = sp3d5;   //C cation
	no_orb[2]            = sp3d5;   //C anion
	no_orb[3]            = sp3d5;   //C cation

	atomic_mass[0]       = 12.0107*amu;  //C
	atomic_mass[1]       = 12.0107*amu;  //C
	atomic_mass[2]       = 12.0107*amu;  //C
	atomic_mass[3]       = 12.0107*amu;  //C

	phiR[0][0]           = 365.0;
	phiR[0][1]           = 88.0;
	phiR[0][2]           = 30.0;
	phiR[0][3]           = -19.2;

	phiT[0][0]           = 245.0;
	phiT[0][1]           = -32.3;;
	phiT[0][2]           = -52.5;
	phiT[0][3]           = 22.9;

	phiZ[0][0]           = 98.2;
	phiZ[0][1]           = -4.0;
	phiZ[0][2]           = 1.5;
	phiZ[0][3]           = -5.8;

	phiR[1][0]           = 365.0;
	phiR[1][1]           = 88.0;
	phiR[1][2]           = 30.0;
	phiR[1][3]           = -19.2;

	phiT[1][0]           = 245.0;
	phiT[1][1]           = -32.3;;
	phiT[1][2]           = -52.5;
	phiT[1][3]           = 22.9;

	phiZ[1][0]           = 98.2;
	phiZ[1][1]           = -4.0;
	phiZ[1][2]           = 1.5;
	phiZ[1][3]           = -5.8;

	phiR[2][0]           = 365.0;
	phiR[2][1]           = 88.0;
	phiR[2][2]           = 30.0;
	phiR[2][3]           = -19.2;

	phiT[2][0]           = 245.0;
	phiT[2][1]           = -32.3;;
	phiT[2][2]           = -52.5;
	phiT[2][3]           = 22.9;

	phiZ[2][0]           = 98.2;
	phiZ[2][1]           = -4.0;
	phiZ[2][2]           = 1.5;
	phiZ[2][3]           = -5.8;

	phiR[3][0]           = 365.0;
	phiR[3][1]           = 88.0;
	phiR[3][2]           = 30.0;
	phiR[3][3]           = -19.2;

	phiT[3][0]           = 245.0;
	phiT[3][1]           = -32.3;;
	phiT[3][2]           = -52.5;
	phiT[3][3]           = 22.9;

	phiZ[3][0]           = 98.2;
	phiZ[3][1]           = -4.0;
	phiZ[3][2]           = 1.5;
	phiZ[3][3]           = -5.8;

	IS                   = 0;
	IPX                  = 1;
	IPY                  = 2;
	IPZ                  = 3;
	ID1                  = 4;
	ID2                  = 5;
	ID3                  = 6;
	ID4                  = 7;
	ID5                  = 8;
	ISS                  = 9;
	
	break;

    case 58:  // bilayer_pz
      
        Eg                   = 0.1;
        ECmin                = 0.05;//0.042;//0.05;
        EVmax                = -0.05;//-0.495;//-0.05;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.142;

//Anion and from anion to
        Esa[0]               = 0.0;
        Epa[0]               = 10*tollim;
        Estara[0]            = 0.0;
        Eda12[0]             = 0.0;
	Eda15[0]             = 0.0;
        lambdaa[0]           = 0.0;     //(Delta/3)
        Vsssa[0]             = 0.0;
        Vstarstarsa[0]       = 0.0;
        Vsstarsa[0]          = 0.0;
        Vspsa[0]             = 0.0;
        Vstarpsa[0]          = 0.0;
        Vsdsa[0]             = 0.0;
        Vstardsa[0]          = 0.0;
        Vppsa[0]             = 0.35;
        Vpppa[0]             = -3.0;
        Vpdsa[0]             = 0.0;
        Vpdpa[0]             = 0.0;
        Vddsa[0]             = 0.0;
        Vddpa[0]             = 0.0;
        Vddda[0]             = 0.0;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//Cation and from cation to anion
        Esc[0]               = 0.0;
        Epc[0]               = 10*tollim;
        Estarc[0]            = 0.0;
        Edc12[0]             = 0.0;
	Edc15[0]             = 0.0;
        lambdac[0]           = 0.0;     //(Delta/3)
        Vsssc[0]             = 0.0;
        Vstarstarsc[0]       = 0.0;
        Vsstarsc[0]          = 0.0;
        Vspsc[0]             = 0.0;
        Vstarpsc[0]          = 0.0;
        Vsdsc[0]             = 0.0;
        Vstardsc[0]          = 0.0;
        Vppsc[0]             = 0.35;
        Vpppc[0]             = -3.0;
        Vpdsc[0]             = 0.0;
        Vpdpc[0]             = 0.0;
        Vddsc[0]             = 0.0;
        Vddpc[0]             = 0.0;
        Vdddc[0]             = 0.0;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 0.0;
        eta_sstars[0]        = 0.0;
        eta_starstars[0]     = 0.0;
        eta_sps[0]           = 0.0;
        eta_starps[0]        = 0.0;
        eta_sds[0]           = 0.0;
        eta_stards[0]        = 0.0;
        eta_pps[0]           = 0.0;
        eta_ppp[0]           = 3.124;  //to ensure q0 = 1/t0*dt0/dr = -eta/r0 = 2.2 1/Ang (Mahan)
        eta_pds[0]           = 0.0;
        eta_pdp[0]           = 0.0;
        eta_dds[0]           = 0.0;
        eta_ddp[0]           = 0.0;
        eta_ddd[0]           = 0.0;        
	
        Eshift[0]            = 27;

        bond_length[0]       = 0.142;

 	neighbor_table[0][1] = 11;
        neighbor_table[0][3] = 11;
        neighbor_table[1][0] = 12;
        neighbor_table[1][2] = 12;
        neighbor_table[2][1] = 11;
        neighbor_table[2][3] = 11;
        neighbor_table[3][0] = 12;
        neighbor_table[3][2] = 12;

	no_orb[0]            = pzorb;   //C anion
	no_orb[1]            = pzorb;   //C cation
	no_orb[2]            = pzorb;   //C anion
	no_orb[3]            = pzorb;   //C cation

	atomic_mass[0]       = 12.0107*amu;  //C
	atomic_mass[1]       = 12.0107*amu;  //C
	atomic_mass[2]       = 12.0107*amu;  //C
	atomic_mass[3]       = 12.0107*amu;  //C

	phiR[0][0]           = 365.0;
	phiR[0][1]           = 88.0;
	phiR[0][2]           = 30.0;
	phiR[0][3]           = -19.2;

	phiT[0][0]           = 245.0;
	phiT[0][1]           = -32.3;;
	phiT[0][2]           = -52.5;
	phiT[0][3]           = 22.9;

	phiZ[0][0]           = 98.2;
	phiZ[0][1]           = -4.0;
	phiZ[0][2]           = 1.5;
	phiZ[0][3]           = -5.8;

	phiR[1][0]           = 365.0;
	phiR[1][1]           = 88.0;
	phiR[1][2]           = 30.0;
	phiR[1][3]           = -19.2;

	phiT[1][0]           = 245.0;
	phiT[1][1]           = -32.3;;
	phiT[1][2]           = -52.5;
	phiT[1][3]           = 22.9;

	phiZ[1][0]           = 98.2;
	phiZ[1][1]           = -4.0;
	phiZ[1][2]           = 1.5;
	phiZ[1][3]           = -5.8;

	phiR[2][0]           = 365.0;
	phiR[2][1]           = 88.0;
	phiR[2][2]           = 30.0;
	phiR[2][3]           = -19.2;

	phiT[2][0]           = 245.0;
	phiT[2][1]           = -32.3;;
	phiT[2][2]           = -52.5;
	phiT[2][3]           = 22.9;

	phiZ[2][0]           = 98.2;
	phiZ[2][1]           = -4.0;
	phiZ[2][2]           = 1.5;
	phiZ[2][3]           = -5.8;

	phiR[3][0]           = 365.0;
	phiR[3][1]           = 88.0;
	phiR[3][2]           = 30.0;
	phiR[3][3]           = -19.2;

	phiT[3][0]           = 245.0;
	phiT[3][1]           = -32.3;;
	phiT[3][2]           = -52.5;
	phiT[3][3]           = 22.9;

	phiZ[3][0]           = 98.2;
	phiZ[3][1]           = -4.0;
	phiZ[3][2]           = 1.5;
	phiZ[3][3]           = -5.8;

	IPZ                  = 0;
        IS                   = 3;

	break;        

    case 100:  // Bi2Te3
      
        Eg                   = 0.158;
        ECmin                = Eg;
        EVmax                = 0.0;

	// Strain constants
        alpha                = 0.0;
        beta                 = 0.0;
        ideal_a0             = 0.0;

//Bi and from Bi to TeI
        Esa[0]               = -9.0657;
        Epa[0]               = 2.0431;
        Estara[0]            = 8.2718;
        Eda12[0]             = 18.0680;
	Eda15[0]             = 18.0680;
        lambdaa[0]           = 0.8982;     //(Delta/3)
        Vsssa[0]             = -0.3735;
        Vstarstarsa[0]       = -0.4424;
        Vsstarsa[0]          = -0.6194;
        Vspsa[0]             = 1.8234;
        Vstarpsa[0]          = 1.2503;
        Vsdsa[0]             = -0.7240;
        Vstardsa[0]          = -0.1250;
        Vppsa[0]             = 2.2860;
        Vpppa[0]             = -0.6192;
        Vpdsa[0]             = -1.8266;
        Vpdpa[0]             = 0.4966;
        Vddsa[0]             = -1.2178;
        Vddpa[0]             = 2.1665;
        Vddda[0]             = -0.0953;
//Strain
        Csasc[0]             = 0.0;
        Cstarastarc[0]       = 0.0;
        Csastarc[0]          = 0.0;
        Csapc[0]             = 0.0;
        Cstarapc[0]          = 0.0;
        Csadc[0]             = 0.0;
        Cstaradc[0]          = 0.0;
        Cpapc[0]             = 0.0;
        Cpadc[0]             = 0.0;
        Cdadc[0]             = 0.0;

//TeI and from TeI to Bi
        Esc[0]               = -10.2050;
        Epc[0]               = -0.5410;
        Estarc[0]            = 14.2024;
        Edc12[0]             = 12.5410;
	Edc15[0]             = 12.5410;
        lambdac[0]           = 0.2362;     //(Delta/3)
        Vsssc[0]             = -0.3735;
        Vstarstarsc[0]       = -0.4424;
        Vsstarsc[0]          = -0.0010;
        Vspsc[0]             = 1.0165;
        Vstarpsc[0]          = 0.6406;
        Vsdsc[0]             = -1.1863;
        Vstardsc[0]          = -0.0010;
        Vppsc[0]             = 2.2860;
        Vpppc[0]             = -0.6192;
        Vpdsc[0]             = -1.4842;
        Vpdpc[0]             = 1.4372;
        Vddsc[0]             = -1.2178;
        Vddpc[0]             = 2.1665;
        Vdddc[0]             = -0.0953;
//Strain
        Cscsa[0]             = 0.0;
        Cstarcstara[0]       = 0.0;
        Cscstara[0]          = 0.0;
        Cscpa[0]             = 0.0;
        Cstarcpa[0]          = 0.0;
        Cscda[0]             = 0.0;
        Cstarcda[0]          = 0.0;
        Cpcpa[0]             = 0.0;
        Cpcda[0]             = 0.0;
        Cdcda[0]             = 0.0;

        eta_sss[0]           = 2.0;
        eta_sstars[0]        = 2.0;
        eta_starstars[0]     = 2.0;
        eta_sps[0]           = 2.0;
        eta_starps[0]        = 2.0;
        eta_sds[0]           = 2.0;
        eta_stards[0]        = 2.0;
        eta_pps[0]           = 2.0;
        eta_ppp[0]           = 2.0;
        eta_pds[0]           = 2.0;
        eta_pdp[0]           = 2.0;
        eta_dds[0]           = 2.0;
        eta_ddp[0]           = 2.0;
        eta_ddd[0]           = 2.0;

        Eshift[0]            = 27;

        bond_length[0]       = 0.30746;
	
//Bi and from Bi to TeII
        Esa[1]               = -9.0657;
        Epa[1]               = 2.0431;
        Estara[1]            = 8.2718;
        Eda12[1]             = 18.0680;
	Eda15[1]             = 18.0680;
        lambdaa[1]           = 0.8982;     //(Delta/3)
        Vsssa[1]             = -0.6797;
        Vstarstarsa[1]       = -0.6857;
        Vsstarsa[1]          = -0.0010;
        Vspsa[1]             = 0.7972;
        Vstarpsa[1]          = 0.0011;
        Vsdsa[1]             = -1.0984;
        Vstardsa[1]          = -0.0378;
        Vppsa[1]             = 1.6122;
        Vpppa[1]             = -0.4125;
        Vpdsa[1]             = -1.8716;
        Vpdpa[1]             = 1.1729;
        Vddsa[1]             = -1.3396;
        Vddpa[1]             = 2.6219;
        Vddda[1]             = -1.2809;
//Strain
        Csasc[1]             = 0.0;
        Cstarastarc[1]       = 0.0;
        Csastarc[1]          = 0.0;
        Csapc[1]             = 0.0;
        Cstarapc[1]          = 0.0;
        Csadc[1]             = 0.0;
        Cstaradc[1]          = 0.0;
        Cpapc[1]             = 0.0;
        Cpadc[1]             = 0.0;
        Cdadc[1]             = 0.0;

//TeII and from TeII to Bi
        Esc[1]               = -10.9796;
        Epc[1]               = -1.3313;
        Estarc[1]            = 14.1324;
        Edc12[1]             = 10.2980;
	Edc15[1]             = 10.2980;
        lambdac[1]           = 0.4326;     //(Delta/3)
        Vsssc[1]             = -0.6797;
        Vstarstarsc[1]       = -0.6857;
        Vsstarsc[1]          = -0.1268;
        Vspsc[1]             = 0.4921;
        Vstarpsc[1]          = 0.4653;
        Vsdsc[1]             = -0.9489;
        Vstardsc[1]          = -0.6360;
        Vppsc[1]             = 1.6122;
        Vpppc[1]             = -0.4125;
        Vpdsc[1]             = -0.6893;
        Vpdpc[1]             = 1.1182;
        Vddsc[1]             = -1.3396;
        Vddpc[1]             = 2.6219;
        Vdddc[1]             = -1.2809;
//Strain
        Cscsa[1]             = 0.0;
        Cstarcstara[1]       = 0.0;
        Cscstara[1]          = 0.0;
        Cscpa[1]             = 0.0;
        Cstarcpa[1]          = 0.0;
        Cscda[1]             = 0.0;
        Cstarcda[1]          = 0.0;
        Cpcpa[1]             = 0.0;
        Cpcda[1]             = 0.0;
        Cdcda[1]             = 0.0;

        eta_sss[1]           = 2.0;
        eta_sstars[1]        = 2.0;
        eta_starstars[1]     = 2.0;
        eta_sps[1]           = 2.0;
        eta_starps[1]        = 2.0;
        eta_sds[1]           = 2.0;
        eta_stards[1]        = 2.0;
        eta_pps[1]           = 2.0;
        eta_ppp[1]           = 2.0;
        eta_pds[1]           = 2.0;
        eta_pdp[1]           = 2.0;
        eta_dds[1]           = 2.0;
        eta_ddp[1]           = 2.0;
        eta_ddd[1]           = 2.0;

        Eshift[1]            = 27;

        bond_length[1]       = 0.324591;

//TeI and from TeI to TeI
        Esa[2]               = -10.2050;
        Epa[2]               = -0.5410;
        Estara[2]            = 14.2024;
        Eda12[2]             = 12.5410;
	Eda15[2]             = 12.5410;
        lambdaa[2]           = 0.2362;     //(Delta/3)
        Vsssa[2]             = -0.3402;
        Vstarstarsa[2]       = -1.3463;
        Vsstarsa[2]          = -0.0011;
        Vspsa[2]             = 0.5896;
        Vstarpsa[2]          = 0.7384;
        Vsdsa[2]             = -0.5926;
        Vstardsa[2]          = -1.0654;
        Vppsa[2]             = 1.1292;
        Vpppa[2]             = -0.0010;
        Vpdsa[2]             = -1.5154;
        Vpdpa[2]             = 1.9631;
        Vddsa[2]             = -1.4381;
        Vddpa[2]             = 2.5731;
        Vddda[2]             = -1.4284;
//Strain
        Csasc[2]             = 0.0;
        Cstarastarc[2]       = 0.0;
        Csastarc[2]          = 0.0;
        Csapc[2]             = 0.0;
        Cstarapc[2]          = 0.0;
        Csadc[2]             = 0.0;
        Cstaradc[2]          = 0.0;
        Cpapc[2]             = 0.0;
        Cpadc[2]             = 0.0;
        Cdadc[2]             = 0.0;

//TeI and from TeI to TeI
        Esc[2]               = -10.2050;
        Epc[2]               = -0.5410;
        Estarc[2]            = 14.2024;
        Edc12[2]             = 12.5410;
	Edc15[2]             = 12.5410;
        lambdac[2]           = 0.2362;     //(Delta/3)
        Vsssc[2]             = -0.3402;
        Vstarstarsc[2]       = -1.3463;
        Vsstarsc[2]          = -0.0011;
        Vspsc[2]             = 0.5896;
        Vstarpsc[2]          = 0.7384;
        Vsdsc[2]             = -0.5926;
        Vstardsc[2]          = -1.0654;
        Vppsc[2]             = 1.1292;
        Vpppc[2]             = -0.0010;
        Vpdsc[2]             = -1.5154;
        Vpdpc[2]             = 1.9631;
        Vddsc[2]             = -1.4381;
        Vddpc[2]             = 2.5731;
        Vdddc[2]             = -1.4284;
//Strain
        Cscsa[2]             = 0.0;
        Cstarcstara[2]       = 0.0;
        Cscstara[2]          = 0.0;
        Cscpa[2]             = 0.0;
        Cstarcpa[2]          = 0.0;
        Cscda[2]             = 0.0;
        Cstarcda[2]          = 0.0;
        Cpcpa[2]             = 0.0;
        Cpcda[2]             = 0.0;
        Cdcda[2]             = 0.0;

        eta_sss[2]           = 2.0;
        eta_sstars[2]        = 2.0;
        eta_starstars[2]     = 2.0;
        eta_sps[2]           = 2.0;
        eta_starps[2]        = 2.0;
        eta_sds[2]           = 2.0;
        eta_stards[2]        = 2.0;
        eta_pps[2]           = 2.0;
        eta_ppp[2]           = 2.0;
        eta_pds[2]           = 2.0;
        eta_pdp[2]           = 2.0;
        eta_dds[2]           = 2.0;
        eta_ddp[2]           = 2.0;
        eta_ddd[2]           = 2.0;

        Eshift[2]            = 27;

        bond_length[2]       = 0.363239;

	neighbor_table[0][2] = 22;
        neighbor_table[0][3] = 22;
	neighbor_table[1][3] = 12;
        neighbor_table[1][4] = 31;
	neighbor_table[2][0] = 21;
        neighbor_table[2][4] = 11;
	neighbor_table[3][0] = 21;
        neighbor_table[3][1] = 11;
	neighbor_table[4][1] = 31;
        neighbor_table[4][2] = 12;

	no_orb[0]            = sp3d5ss;   //TeII anion
	no_orb[1]            = sp3d5ss;   //TeI cation
	no_orb[2]            = sp3d5ss;   //Bi anion
	no_orb[3]            = sp3d5ss;   //Bi cation
	no_orb[4]            = sp3d5ss;   //TeI anion

	atomic_mass[0]       = 127.60000*amu;  //TeII
	atomic_mass[1]       = 127.60000*amu;  //TeI
	atomic_mass[0]       = 208.98038*amu;  //Bi
	atomic_mass[1]       = 208.98038*amu;  //Bi
	atomic_mass[0]       = 127.60000*amu;  //TeI

	break;

    case 200:  // Si new strain model

        Eg                   = 1.13;
        ECmin                = Eg;
        EVmax                = 0;

	// Strain constants
        alpha                = 48.5;
        beta                 = 13.8;
        ideal_a0             = 0.543095;

//Anion and from anion to
        Esa[0]               = -2.15168;
        Epa[0]               = 4.22925;
        Estara[0]            = 19.11650;
        Eda12[0]             = 13.78950;
	Eda15[0]             = 13.78950;
        lambdaa[0]           = 0.01989;     //(Delta/3)
        Vsssa[0]             = -1.95933;
        Vstarstarsa[0]       = -4.24135;
        Vsstarsa[0]          = -1.52230;
        Vspsa[0]             = 3.02562;
        Vstarpsa[0]          = 3.15565;
        Vsdsa[0]             = -2.28485;
        Vstardsa[0]          = -0.80993;
        Vppsa[0]             = 4.10364;
        Vpppa[0]             = -1.51801;
        Vpdsa[0]             = -1.35554;
        Vpdpa[0]             = 2.38479;
        Vddsa[0]             = -1.68136;
        Vddpa[0]             = 2.58880;
        Vddda[0]             = -1.81400;
//Strain
        Csasc[0]             = 1.6788;
        Cstarastarc[0]       = 0.7777;
        Csastarc[0]          = 1.7843;
        Csapc[0]             = 0.4801;
        Cstarapc[0]          = 3.5888;
        Csadc[0]             = 0;
        Cstaradc[0]          = 0.3421;
        Cpapc[0]             = 0.0;//
        Cpadc[0]             = 0.06;//
        Cdadc[0]             = 4.33;//

//Cation and from cation to anion
        Esc[0]               = -2.15168;
        Epc[0]               = 4.22925;
        Estarc[0]            = 19.11650;
        Edc12[0]             = 13.78950;
	Edc15[0]             = 13.78950;
        lambdac[0]           = 0.01989;     //(Delta/3)
        Vsssc[0]             = -1.95933;
        Vstarstarsc[0]       = -4.24135;
        Vsstarsc[0]          = -1.52230;
        Vspsc[0]             = 3.02562;
        Vstarpsc[0]          = 3.15565;
        Vsdsc[0]             = -2.28485;
        Vstardsc[0]          = -0.80993;
        Vppsc[0]             = 4.10364;
        Vpppc[0]             = -1.51801;
        Vpdsc[0]             = -1.35554;
        Vpdpc[0]             = 2.38479;
        Vddsc[0]             = -1.68136;
        Vddpc[0]             = 2.58880;
        Vdddc[0]             = -1.81400;
//Strain
        Cscsa[0]             = 1.6788;
        Cstarcstara[0]       = 0.7777;
        Cscstara[0]          = 1.7843;
        Cscpa[0]             = 0.4801;
        Cstarcpa[0]          = 3.5888;
        Cscda[0]             = 0;
        Cstarcda[0]          = 0.3421;
        Cpcpa[0]             = 0.0;//
        Cpcda[0]             = 0.06;//
        Cdcda[0]             = 4.33;//

        eta_sss[0]           = 0.562469;
        eta_sstars[0]        = 0.132030;
        eta_starstars[0]     = 0.192369;
        eta_sps[0]           = 2.365484;
        eta_starps[0]        = 0.344918;
        eta_sds[0]           = 2.567199;
        eta_stards[0]        = 1.086006;
        eta_pps[0]           = 0.2;//
        eta_ppp[0]           = 1.68;//
        eta_pds[0]           = 0.2;//
        eta_pdp[0]           = 4.43;//
        eta_dds[0]           = 0.1;//
        eta_ddp[0]           = 6.0;//
        eta_ddd[0]           = 6.0;//

        Eshift[0]            = 27;

	r2av_pp[0]           = 0.12;
	r2av_dd[0]           = 0.51;
	r4av_dd[0]           = 0.77;
	r_pd[0]              = 0.12;
	r3_pd[0]             = 0.0;
	Zeff[0]              = 3.0;

        bond_length[0]       = 0.543*sqrt(3.0)/4.0;

        neighbor_table[0][1] = 11;
        neighbor_table[1][0] = 11;

	no_orb[0]            = sp3d5ss; //Si anion
	no_orb[1]            = sp3d5ss; //Si cation

	atomic_mass[0]       = 28.0855*amu;
	atomic_mass[1]       = 28.0855*amu;

	alpha_ph[0]          = 49.4;
	beta_ph[0]           = 4.79;
	kappa_ph[0]          = 6.99;
	tau_ph[0]            = 5.2;
	gamma_ph[0]          = 0.0;

	alpha_ph[1]          = 49.4;
	beta_ph[1]           = 4.79;
	kappa_ph[1]          = 6.99;
	tau_ph[1]            = 5.2;
	gamma_ph[1]          = 0.0;

	alphap_ph[0]         = -2.3241e12/64.0;
	betap_ph[0]          = -1.0388e12/64.0;
	taup_ph[0]           = -2.2447e11/64.0;
	delta1p_ph[0]        = -2.2300e11/64.0;
	delta3p_ph[0]        = 6.3894e11/64.0;

	alphap_ph[1]         = -2.3241e12/64.0;
	betap_ph[1]          = -1.0388e12/64.0;
	taup_ph[1]           = -2.2447e11/64.0;
	delta1p_ph[1]        = -2.2300e11/64.0;
	delta3p_ph[1]        = 6.3894e11/64.0;

        break;

    default:
        read_material_from_file(F);
    }

    if(init_gap_tables){
        for(int i=0;i<table_dim;i++){
	    for(int j=0;j<table_dim;j++){
	        mid_gap_energy[i][j] = (ECmin+EVmax)/2.0;
		band_gap_table[i][j] = Eg;
	    }
	}
    }
    
}

/************************************************************************************************/

void Material::read_material_from_file(FILE *F)
{
    int condition;
    int IA,IC,IT,ind1,ind2,NT,NCA;
    
    condition = (mat_name[0]=='p')&&(mat_name[1]=='h');

    if(!condition){

	fscanf(F,"%lg",&Eg);
	fscanf(F,"%lg",&ECmin);
	fscanf(F,"%lg",&EVmax);
    
	for(IA=0;IA<no_anion;IA++){

	    for(IC=0;IC<no_cation;IC++){
	  
	        //Anion IA and from anion IA to cation IC
	        fscanf(F,"%lg",&Esa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Epa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Estara[IA*no_cation+IC]);
		fscanf(F,"%lg",&Eda12[IA*no_cation+IC]);
		Eda15[IA*no_cation+IC] = Eda12[IA*no_cation+IC];
		fscanf(F,"%lg",&lambdaa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vsssa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vstarstarsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vsstarsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vspsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vstarpsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vsdsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vstardsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vppsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vpppa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vpdsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vpdpa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vddsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vddpa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vddda[IA*no_cation+IC]);

		//Strain
		fscanf(F,"%lg",&Csasc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cstarastarc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Csastarc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Csapc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cstarapc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Csadc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cstaradc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cpapc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cpadc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cdadc[IA*no_cation+IC]);

		//Cation IC and from cation IC to anion IA
		fscanf(F,"%lg",&Esc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Epc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Estarc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Edc12[IA*no_cation+IC]);
		Edc15[IA*no_cation+IC] = Edc12[IA*no_cation+IC];
		fscanf(F,"%lg",&lambdac[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vsssc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vstarstarsc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vsstarsc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vspsc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vstarpsc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vsdsc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vstardsc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vppsc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vpppc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vpdsc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vpdpc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vddsc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vddpc[IA*no_cation+IC]);
		fscanf(F,"%lg",&Vdddc[IA*no_cation+IC]);

		//Strain
		fscanf(F,"%lg",&Cscsa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cstarcstara[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cscstara[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cscpa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cstarcpa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cscda[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cstarcda[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cpcpa[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cpcda[IA*no_cation+IC]);
		fscanf(F,"%lg",&Cdcda[IA*no_cation+IC]);

		fscanf(F,"%lg",&eta_sss[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_sstars[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_starstars[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_sps[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_starps[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_sds[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_stards[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_pps[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_ppp[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_pds[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_pdp[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_dds[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_ddp[IA*no_cation+IC]);
		fscanf(F,"%lg",&eta_ddd[IA*no_cation+IC]);

		fscanf(F,"%lg",&Eshift[IA*no_cation+IC]);

		fscanf(F,"%lg",&bond_length[IA*no_cation+IC]);

	    }   
	}
    }else{

        for(IA=0;IA<no_anion;IA++){
	    for(IC=0;IC<no_cation;IC++){
	        fscanf(F,"%lg",&bond_length[IA*no_cation+IC]);
	    }
	}
    }

    fscanf(F,"%i",&NT);
    
    for(IT=0;IT<NT;IT++){
        fscanf(F,"%i",&ind1);
	fscanf(F,"%i",&ind2);
	fscanf(F,"%i",&neighbor_table[ind1-1][ind2-1]);
    }
    
    fscanf(F,"%i",&NCA);
    for(IA=0;IA<NCA;IA++){
        fscanf(F,"%i",&no_orb[IA]);
    }
    
    for(IA=0;IA<NCA;IA++){
        if(condition){
	    fscanf(F,"%lg",&atomic_mass[IA]);
	    fscanf(F,"%lg",&alpha_ph[IA]);
	    fscanf(F,"%lg",&beta_ph[IA]);
	    fscanf(F,"%lg",&kappa_ph[IA]);
	    fscanf(F,"%lg",&tau_ph[IA]);
	    fscanf(F,"%lg",&gamma_ph[IA]);
	}
	atomic_mass[IA] = atomic_mass[IA]*amu;
    }
    
    fclose(F);

}

/************************************************************************************************/

void Material::get_band_edge(double *pECmin,double *pEVmax)
{
    *pECmin = ECmin;
    *pEVmax = EVmax;
}

/************************************************************************************************/

void Material::calc_matrix_element(double *V, double l, double m, double n, double d0_d, \
                                   int fst_entry, int snd_entry)
{

    double vsssa,vstarstarsa,vsstarsa;
    double vspsa,vstarpsa,vsdsa;
    double vstardsa,vppsa,vpppa;
    double vpdsa,vpdpa,vddsa;
    double vddpa,vddda;
    double vsssc,vstarstarsc,vsstarsc;
    double vspsc,vstarpsc,vsdsc;
    double vstardsc,vppsc,vpppc;
    double vpdsc,vpdpc,vddsc;
    double vddpc,vdddc;
    
    if(!mod_scaling_fcn){
        
        vsssa       = Vsssa[fst_entry]*pow(d0_d,eta_sss[fst_entry]);
        vstarstarsa = Vstarstarsa[fst_entry]*pow(d0_d,eta_starstars[fst_entry]);
        vsstarsa    = Vsstarsa[fst_entry]*pow(d0_d,eta_sstars[fst_entry]);
        vspsa       = Vspsa[fst_entry]*pow(d0_d,eta_sps[fst_entry]);
        vstarpsa    = Vstarpsa[fst_entry]*pow(d0_d,eta_starps[fst_entry]);
        vsdsa       = Vsdsa[fst_entry]*pow(d0_d,eta_sds[fst_entry]);
        vstardsa    = Vstardsa[fst_entry]*pow(d0_d,eta_stards[fst_entry]);
        vppsa       = Vppsa[fst_entry]*pow(d0_d,eta_pps[fst_entry]);
        vpppa       = Vpppa[fst_entry]*pow(d0_d,eta_ppp[fst_entry]);
        vpdsa       = Vpdsa[fst_entry]*pow(d0_d,eta_pds[fst_entry]);
        vpdpa       = Vpdpa[fst_entry]*pow(d0_d,eta_pdp[fst_entry]);
        vddsa       = Vddsa[fst_entry]*pow(d0_d,eta_dds[fst_entry]);
        vddpa       = Vddpa[fst_entry]*pow(d0_d,eta_ddp[fst_entry]);
        vddda       = Vddda[fst_entry]*pow(d0_d,eta_ddd[fst_entry]);

        vsssc       = Vsssc[fst_entry]*pow(d0_d,eta_sss[fst_entry]);
        vstarstarsc = Vstarstarsc[fst_entry]*pow(d0_d,eta_starstars[fst_entry]);
        vsstarsc    = Vsstarsc[fst_entry]*pow(d0_d,eta_sstars[fst_entry]);
        vspsc       = Vspsc[fst_entry]*pow(d0_d,eta_sps[fst_entry]);
        vstarpsc    = Vstarpsc[fst_entry]*pow(d0_d,eta_starps[fst_entry]);
        vsdsc       = Vsdsc[fst_entry]*pow(d0_d,eta_sds[fst_entry]);
        vstardsc    = Vstardsc[fst_entry]*pow(d0_d,eta_stards[fst_entry]);
        vppsc       = Vppsc[fst_entry]*pow(d0_d,eta_pps[fst_entry]);
        vpppc       = Vpppc[fst_entry]*pow(d0_d,eta_ppp[fst_entry]);
        vpdsc       = Vpdsc[fst_entry]*pow(d0_d,eta_pds[fst_entry]);
        vpdpc       = Vpdpc[fst_entry]*pow(d0_d,eta_pdp[fst_entry]);
        vddsc       = Vddsc[fst_entry]*pow(d0_d,eta_dds[fst_entry]);
        vddpc       = Vddpc[fst_entry]*pow(d0_d,eta_ddp[fst_entry]);
        vdddc       = Vdddc[fst_entry]*pow(d0_d,eta_ddd[fst_entry]);
        
    }else{

        vsssa       = Vsssa[fst_entry]*pow(d0_d,eta_sss[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vstarstarsa = Vstarstarsa[fst_entry]*pow(d0_d,eta_starstars[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vsstarsa    = Vsstarsa[fst_entry]*pow(d0_d,eta_sstars[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vspsa       = Vspsa[fst_entry]*pow(d0_d,eta_sps[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vstarpsa    = Vstarpsa[fst_entry]*pow(d0_d,eta_starps[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vsdsa       = Vsdsa[fst_entry]*pow(d0_d,eta_sds[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vstardsa    = Vstardsa[fst_entry]*pow(d0_d,eta_stards[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vppsa       = Vppsa[fst_entry]*pow(d0_d,eta_pps[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vpppa       = Vpppa[fst_entry]*pow(d0_d,eta_ppp[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vpdsa       = Vpdsa[fst_entry]*pow(d0_d,eta_pds[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vpdpa       = Vpdpa[fst_entry]*pow(d0_d,eta_pdp[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vddsa       = Vddsa[fst_entry]*pow(d0_d,eta_dds[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vddpa       = Vddpa[fst_entry]*pow(d0_d,eta_ddp[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vddda       = Vddda[fst_entry]*pow(d0_d,eta_ddd[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));

        vsssc       = Vsssc[fst_entry]*pow(d0_d,eta_sss[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vstarstarsc = Vstarstarsc[fst_entry]*pow(d0_d,eta_starstars[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vsstarsc    = Vsstarsc[fst_entry]*pow(d0_d,eta_sstars[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vspsc       = Vspsc[fst_entry]*pow(d0_d,eta_sps[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vstarpsc    = Vstarpsc[fst_entry]*pow(d0_d,eta_starps[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vsdsc       = Vsdsc[fst_entry]*pow(d0_d,eta_sds[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vstardsc    = Vstardsc[fst_entry]*pow(d0_d,eta_stards[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vppsc       = Vppsc[fst_entry]*pow(d0_d,eta_pps[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vpppc       = Vpppc[fst_entry]*pow(d0_d,eta_ppp[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vpdsc       = Vpdsc[fst_entry]*pow(d0_d,eta_pds[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vpdpc       = Vpdpc[fst_entry]*pow(d0_d,eta_pdp[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vddsc       = Vddsc[fst_entry]*pow(d0_d,eta_dds[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vddpc       = Vddpc[fst_entry]*pow(d0_d,eta_ddp[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));
        vdddc       = Vdddc[fst_entry]*pow(d0_d,eta_ddd[fst_entry])*pow(d0_d,-mu_scal)*exp(-mu_scal*(1.0/d0_d-1.0));    
            
    }
    
    
//Orbitals
//s
//px
//py
//pz
//st
//d1=dxy
//d2=dyz
//d3=dzx
//d4=dx^2-y^2
//d5=d3z^2-r^2
    
    switch(snd_entry){
    case 0:

        //Coupling with <sa|
        V[IS*TB+IS]   = vsssa;
        V[IS*TB+IPX]  = l*vspsa;
        V[IS*TB+IPY]  = m*vspsa;
        V[IS*TB+IPZ]  = n*vspsa;
        V[IS*TB+ISS]  = vsstarsa;
        V[IS*TB+ID1]  = sqrt(3.0)*l*m*vsdsa;
        V[IS*TB+ID2]  = sqrt(3.0)*m*n*vsdsa;
        V[IS*TB+ID3]  = sqrt(3.0)*n*l*vsdsa;
        V[IS*TB+ID4]  = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*vsdsa;
        V[IS*TB+ID5]  = (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vsdsa;

        //Coupling with <pxa|
        V[IPX*TB+IS]  = -l*vspsc;
        V[IPX*TB+IPX] = pow(l,2.0)*vppsa+(1-pow(l,2.0))*vpppa;
        V[IPX*TB+IPY] = l*m*(vppsa-vpppa);
        V[IPX*TB+IPZ] = l*n*(vppsa-vpppa);
        V[IPX*TB+ISS] = -l*vstarpsc;
        V[IPX*TB+ID1] = sqrt(3.0)*pow(l,2.0)*m*vpdsa+m*(1-2.0*pow(l,2.0))*vpdpa;
        V[IPX*TB+ID2] = sqrt(3.0)*l*m*n*vpdsa-2.0*l*m*n*vpdpa;
        V[IPX*TB+ID3] = sqrt(3.0)*pow(l,2.0)*n*vpdsa+n*(1-2.0*pow(l,2.0))*vpdpa;
        V[IPX*TB+ID4] = sqrt(3.0)/2.0*l*(pow(l,2.0)-pow(m,2.0))*vpdsa+l*(1-pow(l,2.0)+pow(m,2.0))*vpdpa;
        V[IPX*TB+ID5] = l*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsa-sqrt(3.0)*l*pow(n,2.0)*vpdpa;

        //Coupling with <pya|
        V[IPY*TB+IS]  = -m*vspsc;
        V[IPY*TB+IPX] = m*l*(vppsa-vpppa);
        V[IPY*TB+IPY] = pow(m,2.0)*vppsa+(1-pow(m,2.0))*vpppa;
        V[IPY*TB+IPZ] = m*n*(vppsa-vpppa);
        V[IPY*TB+ISS] = -m*vstarpsc;
        V[IPY*TB+ID1] = sqrt(3.0)*pow(m,2.0)*l*vpdsa+l*(1-2.0*pow(m,2.0))*vpdpa;
        V[IPY*TB+ID2] = sqrt(3.0)*pow(m,2.0)*n*vpdsa+n*(1-2.0*pow(m,2.0))*vpdpa;
        V[IPY*TB+ID3] = sqrt(3.0)*l*m*n*vpdsa-2.0*l*m*n*vpdpa;
        V[IPY*TB+ID4] = sqrt(3.0)/2.0*m*(pow(l,2.0)-pow(m,2.0))*vpdsa-m*(1-pow(m,2.0)+pow(l,2.0))*vpdpa;
        V[IPY*TB+ID5] = m*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsa-sqrt(3.0)*m*pow(n,2.0)*vpdpa;

        //Coupling with <pza|
        V[IPZ*TB+IS]  = -n*vspsc;
        V[IPZ*TB+IPX] = l*n*(vppsa-vpppa);
        V[IPZ*TB+IPY] = n*m*(vppsa-vpppa);
        V[IPZ*TB+IPZ] = pow(n,2.0)*vppsa+(1-pow(n,2.0))*vpppa;
        V[IPZ*TB+ISS] = -n*vstarpsc;
        V[IPZ*TB+ID1] = sqrt(3.0)*l*m*n*vpdsa-2.0*l*m*n*vpdpa;
        V[IPZ*TB+ID2] = sqrt(3.0)*pow(n,2.0)*m*vpdsa+m*(1-2.0*pow(n,2.0))*vpdpa;
        V[IPZ*TB+ID3] = sqrt(3.0)*pow(n,2.0)*l*vpdsa+l*(1-2.0*pow(n,2.0))*vpdpa;
        V[IPZ*TB+ID4] = sqrt(3.0)/2.0*n*(pow(l,2.0)-pow(m,2.0))*vpdsa-n*(pow(l,2.0)-pow(m,2.0))*vpdpa;
        V[IPZ*TB+ID5] = n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsa+
	  sqrt(3.0)*n*(pow(l,2.0)+pow(m,2.0))*vpdpa;

        //Coupling with <sta|
        V[ISS*TB+IS]  = vsstarsc;
        V[ISS*TB+IPX] = l*vstarpsa;
        V[ISS*TB+IPY] = m*vstarpsa;
        V[ISS*TB+IPZ] = n*vstarpsa;
        V[ISS*TB+ISS] = vstarstarsa;
        V[ISS*TB+ID1] = sqrt(3.0)*l*m*vstardsa;
        V[ISS*TB+ID2] = sqrt(3.0)*m*n*vstardsa;
        V[ISS*TB+ID3] = sqrt(3.0)*n*l*vstardsa;
        V[ISS*TB+ID4] = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*vstardsa;
        V[ISS*TB+ID5] = (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vstardsa;

        //Coupling with <d1a|
        V[ID1*TB+IS]  = sqrt(3.0)*l*m*vsdsc;
        V[ID1*TB+IPX] = -sqrt(3.0)*pow(l,2.0)*m*vpdsc-m*(1-2.0*pow(l,2.0))*vpdpc;
        V[ID1*TB+IPY] = -sqrt(3.0)*pow(m,2.0)*l*vpdsc-l*(1-2.0*pow(m,2.0))*vpdpc;
        V[ID1*TB+IPZ] = -sqrt(3.0)*l*m*n*vpdsc+2.0*l*m*n*vpdpc;
        V[ID1*TB+ISS] = sqrt(3.0)*l*m*vstardsc;
        V[ID1*TB+ID1] = 3.0*pow(l,2.0)*pow(m,2.0)*vddsa+\
	  (pow(l,2.0)+pow(m,2.0)-4.0*pow(l,2.0)*pow(m,2.0))*vddpa+\
	  (pow(n,2.0)+pow(l,2.0)*pow(m,2.0))*vddda;
        V[ID1*TB+ID2] = 3.0*l*pow(m,2.0)*n*vddsa+l*n*(1-4.0*pow(m,2.0))*vddpa+l*n*(pow(m,2.0)-1)*vddda;
        V[ID1*TB+ID3] = 3.0*pow(l,2.0)*m*n*vddsa+m*n*(1-4.0*pow(l,2.0))*vddpa+m*n*(pow(l,2.0)-1)*vddda;
        V[ID1*TB+ID4] = 3.0/2.0*l*m*(pow(l,2.0)-pow(m,2.0))*vddsa+\
	  2.0*l*m*(pow(m,2.0)-pow(l,2.0))*vddpa+l*m*(pow(l,2.0)-pow(m,2.0))/2.0*vddda;
        V[ID1*TB+ID5] = sqrt(3.0)*l*m*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsa-\
	  2.0*sqrt(3.0)*l*m*pow(n,2.0)*vddpa+sqrt(3.0)/2.0*m*l*(1+pow(n,2.0))*vddda;

        //Coupling with <d2a|
        V[ID2*TB+IS]  = sqrt(3.0)*m*n*vsdsc;
        V[ID2*TB+IPX] = -sqrt(3.0)*l*m*n*vpdsc+2.0*l*m*n*vpdpc;
        V[ID2*TB+IPY] = -sqrt(3.0)*pow(m,2.0)*n*vpdsc-n*(1-2.0*pow(m,2.0))*vpdpc;
        V[ID2*TB+IPZ] = -sqrt(3.0)*pow(n,2.0)*m*vpdsc-m*(1-2.0*pow(n,2.0))*vpdpc;
        V[ID2*TB+ISS] = sqrt(3.0)*m*n*vstardsc;
        V[ID2*TB+ID1] = 3.0*l*pow(m,2.0)*n*vddsa+l*n*(1-4.0*pow(m,2.0))*vddpa+l*n*(pow(m,2.0)-1)*vddda;
        V[ID2*TB+ID2] = 3.0*pow(n,2.0)*pow(m,2.0)*vddsa+\
	  (pow(n,2.0)+pow(m,2.0)-4.0*pow(n,2.0)*pow(m,2.0))*vddpa+\
	  (pow(l,2.0)+pow(n,2.0)*pow(m,2.0))*vddda;
        V[ID2*TB+ID3] = 3.0*pow(n,2.0)*l*m*vddsa+m*l*(1-4.0*pow(n,2.0))*vddpa+m*l*(pow(n,2.0)-1)*vddda;
        V[ID2*TB+ID4] = 3.0/2.0*m*n*(pow(l,2.0)-pow(m,2.0))*vddsa-\
	  m*n*(1+2.0*(pow(l,2.0)-pow(m,2.0)))*vddpa+\
	  m*n*(1+(pow(l,2.0)-pow(m,2.0))/2.0)*vddda;
        V[ID2*TB+ID5] = sqrt(3.0)*m*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsa+\
	  sqrt(3.0)*m*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*vddpa-\
	  sqrt(3.0)/2.0*m*n*(pow(l,2.0)+pow(m,2.0))*vddda;

        //Coupling with <d3a|
        V[ID3*TB+IS]  = sqrt(3.0)*l*n*vsdsc;
        V[ID3*TB+IPX] = -sqrt(3.0)*pow(l,2.0)*n*vpdsc-n*(1-2.0*pow(l,2.0))*vpdpc;
        V[ID3*TB+IPY] = -sqrt(3.0)*l*m*n*vpdsc+2.0*l*m*n*vpdpc;
        V[ID3*TB+IPZ] = -sqrt(3.0)*pow(n,2.0)*l*vpdsc-l*(1-2.0*pow(n,2.0))*vpdpc;
        V[ID3*TB+ISS] = sqrt(3.0)*l*n*vstardsc;
        V[ID3*TB+ID1] = 3.0*pow(l,2.0)*m*n*vddsa+m*n*(1-4.0*pow(l,2.0))*vddpa+m*n*(pow(l,2.0)-1)*vddda;
        V[ID3*TB+ID2] = 3.0*pow(n,2.0)*l*m*vddsa+m*l*(1-4.0*pow(n,2.0))*vddpa+m*l*(pow(n,2.0)-1)*vddda;
        V[ID3*TB+ID3] = 3.0*pow(n,2.0)*pow(l,2.0)*vddsa+\
	  (pow(n,2.0)+pow(l,2.0)-4.0*pow(n,2.0)*pow(l,2.0))*vddpa+\
	  (pow(m,2.0)+pow(n,2.0)*pow(l,2.0))*vddda;
        V[ID3*TB+ID4] = 3.0/2.0*n*l*(pow(l,2.0)-pow(m,2.0))*vddsa+\
	  n*l*(1-2.0*(pow(l,2.0)-pow(m,2.0)))*vddpa-\
	  n*l*(1-(pow(l,2.0)-pow(m,2.0))/2.0)*vddda;
        V[ID3*TB+ID5] = sqrt(3.0)*l*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsa+\
	  sqrt(3.0)*l*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*vddpa-\
	  sqrt(3.0)/2.0*l*n*(pow(l,2.0)+pow(m,2.0))*vddda;

        //Coupling with <d4a|
        V[ID4*TB+IS]  = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*vsdsc;
        V[ID4*TB+IPX] = -sqrt(3.0)/2.0*l*(pow(l,2.0)-pow(m,2.0))*vpdsc-l*(1-pow(l,2.0)+pow(m,2.0))*vpdpc;
        V[ID4*TB+IPY] = -sqrt(3.0)/2.0*m*(pow(l,2.0)-pow(m,2.0))*vpdsc+m*(1-pow(m,2.0)+pow(l,2.0))*vpdpc;
        V[ID4*TB+IPZ] = -sqrt(3.0)/2.0*n*(pow(l,2.0)-pow(m,2.0))*vpdsc+n*(pow(l,2.0)-pow(m,2.0))*vpdpc;
        V[ID4*TB+ISS] = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*vstardsc;
        V[ID4*TB+ID1] = 3.0/2.0*l*m*(pow(l,2.0)-pow(m,2.0))*vddsa+\
	  2.0*l*m*(pow(m,2.0)-pow(l,2.0))*vddpa+l*m*(pow(l,2.0)-pow(m,2.0))/2.0*vddda;
        V[ID4*TB+ID2] = 3.0/2.0*m*n*(pow(l,2.0)-pow(m,2.0))*vddsa-\
	  m*n*(1+2.0*(pow(l,2.0)-pow(m,2.0)))*vddpa+m*n*(1+(pow(l,2.0)-pow(m,2.0))/2.0)*vddda;
        V[ID4*TB+ID3] = 3.0/2.0*n*l*(pow(l,2.0)-pow(m,2.0))*vddsa+\
	  n*l*(1-2.0*(pow(l,2.0)-pow(m,2.0)))*vddpa-\
	  n*l*(1-(pow(l,2.0)-pow(m,2.0))/2.0)*vddda;
        V[ID4*TB+ID4] = 3.0/4.0*pow(pow(l,2.0)-pow(m,2.0),2.0)*vddsa+\
	  (pow(l,2.0)+pow(m,2.0)-pow(pow(l,2.0)-pow(m,2.0),2.0))*vddpa+\
	  (pow(n,2.0)+pow(pow(l,2.0)-pow(m,2.0),2.0)/4.0)*vddda;
        V[ID4*TB+ID5] = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*\
	  (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsa+\
	  sqrt(3.0)*pow(n,2.0)*(pow(m,2.0)-pow(l,2.0))*vddpa+\
	  sqrt(3.0)/4.0*(1+pow(n,2.0))*(pow(l,2.0)-pow(m,2.0))*vddda;

        //Coupling with <d5a|
        V[ID5*TB+IS]  = (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vsdsc;
        V[ID5*TB+IPX] = -l*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsc+sqrt(3.0)*l*pow(n,2.0)*vpdpc;
        V[ID5*TB+IPY] = -m*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsc+sqrt(3.0)*m*pow(n,2.0)*vpdpc;
        V[ID5*TB+IPZ] = -n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsc-\
	  sqrt(3.0)*n*(pow(l,2.0)+pow(m,2.0))*vpdpc;
        V[ID5*TB+ISS] = (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vstardsc;
        V[ID5*TB+ID1] = sqrt(3.0)*l*m*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsa-\
	  2.0*sqrt(3.0)*l*m*pow(n,2.0)*vddpa+sqrt(3.0)/2.0*m*l*(1+pow(n,2.0))*vddda;
        V[ID5*TB+ID2] = sqrt(3.0)*m*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsa+\
	  sqrt(3.0)*m*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*vddpa-
	  sqrt(3.0)/2.0*m*n*(pow(l,2.0)+pow(m,2.0))*vddda;
        V[ID5*TB+ID3] = sqrt(3.0)*l*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsa+\
	  sqrt(3.0)*l*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*vddpa-
	  sqrt(3.0)/2.0*l*n*(pow(l,2.0)+pow(m,2.0))*vddda;
        V[ID5*TB+ID4] = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*\
	  (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsa+\
	  sqrt(3.0)*pow(n,2.0)*(pow(m,2.0)-pow(l,2.0))*vddpa+\
	  sqrt(3.0)/4.0*(1+pow(n,2.0))*(pow(l,2.0)-pow(m,2.0))*vddda;
        V[ID5*TB+ID5] = pow(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0,2.0)*vddsa+\
	  3.0*pow(n,2.0)*(pow(l,2.0)+pow(m,2.0))*vddpa+\
	  3.0/4.0*pow(pow(l,2.0)+pow(m,2.0),2.0)*vddda;
        break;
        
    case 1:

        //Coupling with <sc|
        V[IS*TB+IS]   = vsssc;
        V[IS*TB+IPX]  = l*vspsc;
        V[IS*TB+IPY]  = m*vspsc;
        V[IS*TB+IPZ]  = n*vspsc;
        V[IS*TB+ISS]  = vsstarsc;
        V[IS*TB+ID1]  = sqrt(3.0)*l*m*vsdsc;
        V[IS*TB+ID2]  = sqrt(3.0)*m*n*vsdsc;
        V[IS*TB+ID3]  = sqrt(3.0)*n*l*vsdsc;
        V[IS*TB+ID4]  = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*vsdsc;
        V[IS*TB+ID5]  = (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vsdsc;

        //Coupling with <pxc|
        V[IPX*TB+IS]  = -l*vspsa;
        V[IPX*TB+IPX] = pow(l,2.0)*vppsc+(1-pow(l,2.0))*vpppc;
        V[IPX*TB+IPY] = l*m*(vppsc-vpppc);
        V[IPX*TB+IPZ] = l*n*(vppsc-vpppc);
        V[IPX*TB+ISS] = -l*vstarpsa;
        V[IPX*TB+ID1] = sqrt(3.0)*pow(l,2.0)*m*vpdsc+m*(1-2.0*pow(l,2.0))*vpdpc;
        V[IPX*TB+ID2] = sqrt(3.0)*l*m*n*vpdsc-2.0*l*m*n*vpdpc;
        V[IPX*TB+ID3] = sqrt(3.0)*pow(l,2.0)*n*vpdsc+n*(1-2.0*pow(l,2.0))*vpdpc;
        V[IPX*TB+ID4] = sqrt(3.0)/2.0*l*(pow(l,2.0)-pow(m,2.0))*vpdsc+l*(1-pow(l,2.0)+pow(m,2.0))*vpdpc;
        V[IPX*TB+ID5] = l*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsc-sqrt(3.0)*l*pow(n,2.0)*vpdpc;

        //Coupling with <pyc|
        V[IPY*TB+IS]  = -m*vspsa;
        V[IPY*TB+IPX] = m*l*(vppsc-vpppc);
        V[IPY*TB+IPY] = pow(m,2.0)*vppsc+(1-pow(m,2.0))*vpppc;
        V[IPY*TB+IPZ] = m*n*(vppsc-vpppc);
        V[IPY*TB+ISS] = -m*vstarpsa;
        V[IPY*TB+ID1] = sqrt(3.0)*pow(m,2.0)*l*vpdsc+l*(1-2.0*pow(m,2.0))*vpdpc;
        V[IPY*TB+ID2] = sqrt(3.0)*pow(m,2.0)*n*vpdsc+n*(1-2.0*pow(m,2.0))*vpdpc;
        V[IPY*TB+ID3] = sqrt(3.0)*l*m*n*vpdsc-2.0*l*m*n*vpdpc;
        V[IPY*TB+ID4] = sqrt(3.0)/2.0*m*(pow(l,2.0)-pow(m,2.0))*vpdsc-m*(1-pow(m,2.0)+pow(l,2.0))*vpdpc;
        V[IPY*TB+ID5] = m*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsc-sqrt(3.0)*m*pow(n,2.0)*vpdpc;

        //Coupling with <pzc|
        V[IPZ*TB+IS]  = -n*vspsa;
        V[IPZ*TB+IPX] = l*n*(vppsc-vpppc);
        V[IPZ*TB+IPY] = n*m*(vppsc-vpppc);
        V[IPZ*TB+IPZ] = pow(n,2.0)*vppsc+(1-pow(n,2.0))*vpppc;
        V[IPZ*TB+ISS] = -n*vstarpsa;
        V[IPZ*TB+ID1] = sqrt(3.0)*l*m*n*vpdsc-2.0*l*m*n*vpdpc;
        V[IPZ*TB+ID2] = sqrt(3.0)*pow(n,2.0)*m*vpdsc+m*(1-2.0*pow(n,2.0))*vpdpc;
        V[IPZ*TB+ID3] = sqrt(3.0)*pow(n,2.0)*l*vpdsc+l*(1-2.0*pow(n,2.0))*vpdpc;
        V[IPZ*TB+ID4] = sqrt(3.0)/2.0*n*(pow(l,2.0)-pow(m,2.0))*vpdsc-n*(pow(l,2.0)-pow(m,2.0))*vpdpc;
        V[IPZ*TB+ID5] = n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsc+\
	  sqrt(3.0)*n*(pow(l,2.0)+pow(m,2.0))*vpdpc;

        //Coupling with <stc|
        V[ISS*TB+IS]  = vsstarsa;
        V[ISS*TB+IPX] = l*vstarpsc;
        V[ISS*TB+IPY] = m*vstarpsc;
        V[ISS*TB+IPZ] = n*vstarpsc;
        V[ISS*TB+ISS] = vstarstarsc;
        V[ISS*TB+ID1] = sqrt(3.0)*l*m*vstardsc;
        V[ISS*TB+ID2] = sqrt(3.0)*m*n*vstardsc;
        V[ISS*TB+ID3] = sqrt(3.0)*n*l*vstardsc;
        V[ISS*TB+ID4] = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*vstardsc;
        V[ISS*TB+ID5] = (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vstardsc;

        //Coupling with <d1c|
        V[ID1*TB+IS]  = sqrt(3.0)*l*m*vsdsa;
        V[ID1*TB+IPX] = -sqrt(3.0)*pow(l,2.0)*m*vpdsa-m*(1-2.0*pow(l,2.0))*vpdpa;
        V[ID1*TB+IPY] = -sqrt(3.0)*pow(m,2.0)*l*vpdsa-l*(1-2.0*pow(m,2.0))*vpdpa;
        V[ID1*TB+IPZ] = -sqrt(3.0)*l*m*n*vpdsa+2.0*l*m*n*vpdpa;
        V[ID1*TB+ISS] = sqrt(3.0)*l*m*vstardsa;
        V[ID1*TB+ID1] = 3.0*pow(l,2.0)*pow(m,2.0)*vddsc+\
	  (pow(l,2.0)+pow(m,2.0)-4.0*pow(l,2.0)*pow(m,2.0))*vddpc+\
	  (pow(n,2.0)+pow(l,2.0)*pow(m,2.0))*vdddc;
        V[ID1*TB+ID2] = 3.0*l*pow(m,2.0)*n*vddsc+l*n*(1-4.0*pow(m,2.0))*vddpc+l*n*(pow(m,2.0)-1)*vdddc;
        V[ID1*TB+ID3] = 3.0*pow(l,2.0)*m*n*vddsc+m*n*(1-4.0*pow(l,2.0))*vddpc+m*n*(pow(l,2.0)-1)*vdddc;
        V[ID1*TB+ID4] = 3.0/2.0*l*m*(pow(l,2.0)-pow(m,2.0))*vddsc+\
	  2.0*l*m*(pow(m,2.0)-pow(l,2.0))*vddpc+l*m*(pow(l,2.0)-pow(m,2.0))/2.0*vdddc;
        V[ID1*TB+ID5] = sqrt(3.0)*l*m*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsc-\
	  2.0*sqrt(3.0)*l*m*pow(n,2.0)*vddpc+sqrt(3.0)/2.0*m*l*(1+pow(n,2.0))*vdddc;

        //Coupling with <d2c|
        V[ID2*TB+IS]  = sqrt(3.0)*m*n*vsdsa;
        V[ID2*TB+IPX] = -sqrt(3.0)*l*m*n*vpdsa+2.0*l*m*n*vpdpa;
        V[ID2*TB+IPY] = -sqrt(3.0)*pow(m,2.0)*n*vpdsa-n*(1-2.0*pow(m,2.0))*vpdpa;
        V[ID2*TB+IPZ] = -sqrt(3.0)*pow(n,2.0)*m*vpdsa-m*(1-2.0*pow(n,2.0))*vpdpa;
        V[ID2*TB+ISS] = sqrt(3.0)*m*n*vstardsa;
        V[ID2*TB+ID1] = 3.0*l*pow(m,2.0)*n*vddsc+l*n*(1-4.0*pow(m,2.0))*vddpc+l*n*(pow(m,2.0)-1)*vdddc;
        V[ID2*TB+ID2] = 3.0*pow(n,2.0)*pow(m,2.0)*vddsc+\
	  (pow(n,2.0)+pow(m,2.0)-4.0*pow(n,2.0)*pow(m,2.0))*vddpc+\
	  (pow(l,2.0)+pow(n,2.0)*pow(m,2.0))*vdddc;
        V[ID2*TB+ID3] = 3.0*pow(n,2.0)*l*m*vddsc+m*l*(1-4.0*pow(n,2.0))*vddpc+m*l*(pow(n,2.0)-1)*vdddc;
        V[ID2*TB+ID4] = 3.0/2.0*m*n*(pow(l,2.0)-pow(m,2.0))*vddsc-\
	  m*n*(1+2.0*(pow(l,2.0)-pow(m,2.0)))*vddpc+\
	  m*n*(1+(pow(l,2.0)-pow(m,2.0))/2.0)*vdddc;
        V[ID2*TB+ID5] = sqrt(3.0)*m*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsc+\
	  sqrt(3.0)*m*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*vddpc-\
	  sqrt(3.0)/2.0*m*n*(pow(l,2.0)+pow(m,2.0))*vdddc;

        //Coupling with <d3c|
        V[ID3*TB+IS]  = sqrt(3.0)*l*n*vsdsa;
        V[ID3*TB+IPX] = -sqrt(3.0)*pow(l,2.0)*n*vpdsa-n*(1-2.0*pow(l,2.0))*vpdpa;
        V[ID3*TB+IPY] = -sqrt(3.0)*l*m*n*vpdsa+2.0*l*m*n*vpdpa;
        V[ID3*TB+IPZ] = -sqrt(3.0)*pow(n,2.0)*l*vpdsa-l*(1-2.0*pow(n,2.0))*vpdpa;
        V[ID3*TB+ISS] = sqrt(3.0)*l*n*vstardsa;
        V[ID3*TB+ID1] = 3.0*pow(l,2.0)*m*n*vddsc+m*n*(1-4.0*pow(l,2.0))*vddpc+m*n*(pow(l,2.0)-1)*vdddc;
        V[ID3*TB+ID2] = 3.0*pow(n,2.0)*l*m*vddsc+m*l*(1-4.0*pow(n,2.0))*vddpc+m*l*(pow(n,2.0)-1)*vdddc;
        V[ID3*TB+ID3] = 3.0*pow(n,2.0)*pow(l,2.0)*vddsc+\
	  (pow(n,2.0)+pow(l,2.0)-4.0*pow(n,2.0)*pow(l,2.0))*vddpc+\
	  (pow(m,2.0)+pow(n,2.0)*pow(l,2.0))*vdddc;
        V[ID3*TB+ID4] = 3.0/2.0*n*l*(pow(l,2.0)-pow(m,2.0))*vddsc+\
	  n*l*(1-2.0*(pow(l,2.0)-pow(m,2.0)))*vddpc-\
	  n*l*(1-(pow(l,2.0)-pow(m,2.0))/2.0)*vdddc;
        V[ID3*TB+ID5] = sqrt(3.0)*l*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsc+\
	  sqrt(3.0)*l*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*vddpc-\
	  sqrt(3.0)/2.0*l*n*(pow(l,2.0)+pow(m,2.0))*vdddc;

        //Coupling with <d4c|
        V[ID4*TB+IS]  = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*vsdsa;
        V[ID4*TB+IPX] = -sqrt(3.0)/2.0*l*(pow(l,2.0)-pow(m,2.0))*vpdsa-l*(1-pow(l,2.0)+pow(m,2.0))*vpdpa;
        V[ID4*TB+IPY] = -sqrt(3.0)/2.0*m*(pow(l,2.0)-pow(m,2.0))*vpdsa+m*(1-pow(m,2.0)+pow(l,2.0))*vpdpa;
        V[ID4*TB+IPZ] = -sqrt(3.0)/2.0*n*(pow(l,2.0)-pow(m,2.0))*vpdsa+n*(pow(l,2.0)-pow(m,2.0))*vpdpa;
        V[ID4*TB+ISS] = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*vstardsa;
        V[ID4*TB+ID1] = 3.0/2.0*l*m*(pow(l,2.0)-pow(m,2.0))*vddsc+\
	  2.0*l*m*(pow(m,2.0)-pow(l,2.0))*vddpc+l*m*(pow(l,2.0)-pow(m,2.0))/2.0*vdddc;
        V[ID4*TB+ID2] = 3.0/2.0*m*n*(pow(l,2.0)-pow(m,2.0))*vddsc-\
	  m*n*(1+2.0*(pow(l,2.0)-pow(m,2.0)))*vddpc+m*n*(1+(pow(l,2.0)-pow(m,2.0))/2.0)*vdddc;
        V[ID4*TB+ID3] = 3.0/2.0*n*l*(pow(l,2.0)-pow(m,2.0))*vddsc+\
	  n*l*(1-2.0*(pow(l,2.0)-pow(m,2.0)))*vddpc-\
	  n*l*(1-(pow(l,2.0)-pow(m,2.0))/2.0)*vdddc;
        V[ID4*TB+ID4] = 3.0/4.0*pow(pow(l,2.0)-pow(m,2.0),2.0)*vddsc+\
	  (pow(l,2.0)+pow(m,2.0)-pow(pow(l,2.0)-pow(m,2.0),2.0))*vddpc+\
	  (pow(n,2.0)+pow(pow(l,2.0)-pow(m,2.0),2.0)/4.0)*vdddc;
        V[ID4*TB+ID5] = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*\
	  (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsc+\
	  sqrt(3.0)*pow(n,2.0)*(pow(m,2.0)-pow(l,2.0))*vddpc+\
	  sqrt(3.0)/4.0*(1+pow(n,2.0))*(pow(l,2.0)-pow(m,2.0))*vdddc;

        //Coupling with <d5c|
        V[ID5*TB+IS]  = (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vsdsa;
        V[ID5*TB+IPX] = -l*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsa+sqrt(3.0)*l*pow(n,2.0)*vpdpa;
        V[ID5*TB+IPY] = -m*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsa+sqrt(3.0)*m*pow(n,2.0)*vpdpa;
        V[ID5*TB+IPZ] = -n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vpdsa-\
	  sqrt(3.0)*n*(pow(l,2.0)+pow(m,2.0))*vpdpa;
        V[ID5*TB+ISS] = (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vstardsa;
        V[ID5*TB+ID1] = sqrt(3.0)*l*m*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsc-\
	  2.0*sqrt(3.0)*l*m*pow(n,2.0)*vddpc+sqrt(3.0)/2.0*m*l*(1+pow(n,2.0))*vdddc;
        V[ID5*TB+ID2] = sqrt(3.0)*m*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsc+\
	  sqrt(3.0)*m*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*vddpc-
	  sqrt(3.0)/2.0*m*n*(pow(l,2.0)+pow(m,2.0))*vdddc;
        V[ID5*TB+ID3] = sqrt(3.0)*l*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsc+\
	  sqrt(3.0)*l*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*vddpc-
	  sqrt(3.0)/2.0*l*n*(pow(l,2.0)+pow(m,2.0))*vdddc;
        V[ID5*TB+ID4] = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*\
	  (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*vddsc+\
	  sqrt(3.0)*pow(n,2.0)*(pow(m,2.0)-pow(l,2.0))*vddpc+\
	  sqrt(3.0)/4.0*(1+pow(n,2.0))*(pow(l,2.0)-pow(m,2.0))*vdddc;
        V[ID5*TB+ID5] = pow(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0,2.0)*vddsc+\
	  3.0*pow(n,2.0)*(pow(l,2.0)+pow(m,2.0))*vddpc+\
	  3.0/4.0*pow(pow(l,2.0)+pow(m,2.0),2.0)*vdddc;

        break;
    }
}

/************************************************************************************************/

void Material::get_d_matrix(CSR *H00, int fst_entry, int snd_entry, int tb)
{	
    
    int IR,IC;
    int ip      = snd_entry*TB;
    int il      = snd_entry;
    double L[2] = {lambdaa[fst_entry],lambdac[fst_entry]};
    double P[20];

    if(mat_code==57){

        P[IS]       = Esa[fst_entry];
	P[IPX]      = Epa1[fst_entry];
	P[IPY]      = Epa2[fst_entry];
	P[IPZ]      = Epa3[fst_entry];
	P[ISS]      = Estara[fst_entry];
	P[ID1]      = Eda1[fst_entry];
	P[ID2]      = Eda2[fst_entry];
	P[ID3]      = Eda3[fst_entry];
	P[ID4]      = Eda4[fst_entry];
	P[ID5]      = Eda5[fst_entry];

	P[TB+IS]    = Esc[fst_entry];
	P[TB+IPX]   = Epc1[fst_entry];
	P[TB+IPY]   = Epc2[fst_entry];
	P[TB+IPZ]   = Epc3[fst_entry];
	P[TB+ISS]   = Estarc[fst_entry];
	P[TB+ID1]   = Edc1[fst_entry];
	P[TB+ID2]   = Edc2[fst_entry];
	P[TB+ID3]   = Edc3[fst_entry];
	P[TB+ID4]   = Edc4[fst_entry];
	P[TB+ID5]   = Edc5[fst_entry];

    }else{

        P[IS]       = Esa[fst_entry];
	P[IPX]      = Epa[fst_entry];
	P[IPY]      = Epa[fst_entry];
	P[IPZ]      = Epa[fst_entry];
	P[ISS]      = Estara[fst_entry];
	P[ID1]      = (2.0*Eda12[fst_entry]+Eda15[fst_entry])/3.0;
	P[ID2]      = (Eda12[fst_entry]+2.0*Eda15[fst_entry])/3.0;
	P[ID3]      = (Eda12[fst_entry]+2.0*Eda15[fst_entry])/3.0;
	P[ID4]      = (2.0*Eda12[fst_entry]+Eda15[fst_entry])/3.0;
	P[ID5]      = Eda15[fst_entry];

	P[TB+IS]    = Esc[fst_entry];
	P[TB+IPX]   = Epc[fst_entry];
	P[TB+IPY]   = Epc[fst_entry];
	P[TB+IPZ]   = Epc[fst_entry];
	P[TB+ISS]   = Estarc[fst_entry];
	P[TB+ID1]   = (2.0*Edc12[fst_entry]+Edc15[fst_entry])/3.0;
	P[TB+ID2]   = (Edc12[fst_entry]+2.0*Edc15[fst_entry])/3.0;
	P[TB+ID3]   = (Edc12[fst_entry]+2.0*Edc15[fst_entry])/3.0;
	P[TB+ID4]   = (2.0*Edc12[fst_entry]+Edc15[fst_entry])/3.0;
	P[TB+ID5]   = Edc15[fst_entry];
    }

    init_var(H00->r_nnz,tb*tb);
    init_var(H00->i_nnz,tb*tb);

    for(IR=0;IR<tb;IR++){
        H00->r_nnz[IR*tb+IR] = P[ip+IR%TB];
	H00->index_i[IR]     = tb;
	H00->diag_pos[IR]    = (tb+1)*IR;
	for(IC=0;IC<tb;IC++){
	    H00->index_j[IR*tb+IC] = IC;
	}	    
    }

    if(tb==20){
        H00->r_nnz[tb+TB+3]         = L[il];
	H00->r_nnz[3*tb+TB+1]       = -L[il];
	H00->r_nnz[TB*tb+3*tb+1]    = L[il];
	H00->r_nnz[TB*tb+tb+3]      = -L[il];
	H00->i_nnz[tb+2]            = -L[il];
	H00->i_nnz[2*tb+1]          = L[il];
	H00->i_nnz[TB*tb+tb+TB+2]   = L[il];
	H00->i_nnz[TB*tb+2*tb+TB+1] = -L[il];
	H00->i_nnz[2*tb+TB+3]       = -L[il];
	H00->i_nnz[3*tb+TB+2]       = L[il];
	H00->i_nnz[TB*tb+3*tb+2]    = L[il];
	H00->i_nnz[TB*tb+2*tb+3]    = -L[il];
    }

    H00->size       = tb;
    H00->n_nonzeros = tb*tb;

    H00->get_row_edge();

    /*
    switch(tb){
        
    case 10:{
      
        double elements_ns[10] = {P[ip],P[ip+1],P[ip+2],P[ip+3],P[ip+4],P[ip+5],P[ip+6],\
				  P[ip+7],P[ip+8],P[ip+9]};
        int index_i[10]        = {1,1,1,1,1,1,1,1,1,1};
        int index_j[10]        = {0,1,2,3,4,5,6,7,8,9};
        int edge_i[10]         = {0,1,2,3,4,5,6,7,8,9};
        int diag_pos[10]       = {0,1,2,3,4,5,6,7,8,9};

        c_dcopy(TB,elements_ns,1,H00->r_nnz,1);
        icopy(TB,index_i,H00->index_i);
        icopy(TB,index_j,H00->index_j);
        icopy(TB,edge_i,H00->edge_i);
        icopy(TB,diag_pos,H00->diag_pos);
	
	H00->size       = tb;
	H00->n_nonzeros = TB;

        break;}
    
    case 20:{
      
        double r_elements_s[32] = {P[ip],P[ip+1],0.0,L[il],0.0,P[ip+2],0.0,P[ip+3],-L[il],0.0,\
				   P[ip+4],P[ip+5],P[ip+6],P[ip+7],P[ip+8],P[ip+9],P[ip],-L[il],\
				   P[ip+1],0.0,0.0,0.0,P[ip+2],L[il],0.0,P[ip+3],P[ip+4],\
				   P[ip+5],P[ip+6],P[ip+7],P[ip+8],P[ip+9]};
        double i_elements_s[32] = {0.0,0.0,-L[il],0.0,L[il],0.0,-L[il],0.0,0.0,L[il],0.0,0.0,\
				   0.0,0.0,0.0,0.0,0.0,0.0,0.0,L[il],-L[il],-L[il],0.0,0.0,\
				   L[il],0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        int index_i[20]         = {1,3,3,3,1,1,1,1,1,1,1,3,3,3,1,1,1,1,1,1};
        int index_j[32]         = {0,1,2,13,1,2,13,3,11,12,4,5,6,7,8,9,10,3,11,12,3,11,12,1,2,13,14,\
				   15,16,17,18,19};
        int edge_i[20]          = {0,1,4,7,10,11,12,13,14,15,16,17,20,23,26,27,28,29,30,31};
        int diag_pos[20]        = {0,1,5,7,10,11,12,13,14,15,16,18,22,25,26,27,28,29,30,31};

        c_dcopy(32,r_elements_s,1,H00->r_nnz,1);
        c_dcopy(32,i_elements_s,1,H00->i_nnz,1);
        icopy(tb,index_i,H00->index_i);
        icopy(32,index_j,H00->index_j);
        icopy(tb,edge_i,H00->edge_i);
        icopy(tb,diag_pos,H00->diag_pos);

	H00->size       = tb;
	H00->n_nonzeros = tb+12;

        break;}
    }
    */
}

/************************************************************************************************/

void Material::get_nd_matrix(CSR *H01,double l,double m,double n,double d,\
                             int fst_entry,int snd_entry,int tb)
{
    int i,j,k,ind=0;
    double d0_d,V0[100];

    d0_d = bond_length[fst_entry]/d;
    
    calc_matrix_element(V0,l,m,n,d0_d,fst_entry,snd_entry);

    for(k=0;k<tb/TB;k++){
        for(i=0;i<TB;i++){
            H01->index_i[k*TB+i] = 0;
        }
        for(i=0;i<TB;i++){
            for(j=0;j<TB;j++){
                H01->index_i[k*TB+i] = H01->index_i[k*TB+i]+1;
                H01->index_j[ind]    = k*TB+j;
                H01->edge_i[k*TB+i]  = TB*(k*TB+i);
                H01->r_nnz[ind]      = V0[i*TB+j];
		H01->i_nnz[ind]      = 0.0;
                ind++;
            }
        }
    }

    H01->size       = tb;
    H01->n_nonzeros = ind;
    H01->edge_i[tb] = ind;
}

/************************************************************************************************/

void Material::get_n_and_nd_matrix_derivative(double *dH00,CSR *dH01,CSR *H01,double l,double m,\
					      double n,double d,int fst_entry,int snd_entry,\
					      int direction,int tb)
{

    double xyz[3],xyz0[3];
    double dxyz    = 1e-8;
    double d0      = d;

    xyz[0]         = l*d;
    xyz[1]         = m*d;
    xyz[2]         = n*d;
   
    c_dcopy(N3D,xyz,1,xyz0,1);

    xyz[direction] = xyz[direction]+dxyz;
    
    d              = sqrt(pow(xyz[0],2.0)+pow(xyz[1],2.0)+pow(xyz[2],2.0));
    l              = xyz[0]/d;
    m              = xyz[1]/d;
    n              = xyz[2]/d;
    
    if(sc_dist_dep) d0 = d;

    //Harrisson's Scaling Rule of Matrix Element taken into Account if sc_dist_dep=1
    get_nd_matrix(dH01,l,m,n,d0,fst_entry,snd_entry,tb);

    //Variation of the On-Site Energies taken into Account if sc_diag_def=1
    get_strain_shift(dH00,H01->r_nnz,dH01->r_nnz,xyz0,xyz,fst_entry,snd_entry,tb);
    c_dscal(TB*TB,1.0/dxyz,dH00,1);
    
    c_daxpy(dH01->n_nonzeros,-1.0,H01->r_nnz,1,dH01->r_nnz,1);
    c_dscal(dH01->n_nonzeros,1.0/dxyz,dH01->r_nnz,1);
}

/************************************************************************************************/

void Material::prepare_strain_constants(double *strain_const,int fst_entry,int snd_entry)
{

    double csasc       = Csasc[fst_entry];
    double cstarastarc = Cstarastarc[fst_entry];
    double csastarc    = Csastarc[fst_entry];
    double csapc       = Csapc[fst_entry];
    double cstarapc    = Cstarapc[fst_entry];
    double csadc       = Csadc[fst_entry];
    double cstaradc    = Cstaradc[fst_entry];
    double cpapc       = Cpapc[fst_entry];
    double cpadc       = Cpadc[fst_entry];
    double cdadc       = Cdadc[fst_entry];

    double cscsa       = Cscsa[fst_entry];
    double cstarcstara = Cstarcstara[fst_entry];
    double cscstara    = Cscstara[fst_entry];
    double cscpa       = Cscpa[fst_entry];
    double cstarcpa    = Cstarcpa[fst_entry];
    double cscda       = Cscda[fst_entry];
    double cstarcda    = Cstarcda[fst_entry];
    double cpcpa       = Cpcpa[fst_entry];
    double cpcda       = Cpcda[fst_entry];
    double cdcda       = Cdcda[fst_entry];
    
    switch(snd_entry){
    case 0:{
        double SC0[100]={csasc,csapc,csapc,csapc,csastarc,csadc,csadc,csadc,csadc,csadc,\
                         cscpa,cpapc,cpapc,cpapc,cstarcpa,cpadc,cpadc,cpadc,cpadc,cpadc,\
                         cscpa,cpapc,cpapc,cpapc,cstarcpa,cpadc,cpadc,cpadc,cpadc,cpadc,\
                         cscpa,cpapc,cpapc,cpapc,cstarcpa,cpadc,cpadc,cpadc,cpadc,cpadc,\
                         cscstara,cstarapc,cstarapc,cstarapc,cstarastarc,cstaradc,cstaradc,\
                         cstaradc,cstaradc,cstaradc,\
                         cscda,cpcda,cpcda,cpcda,cstarcda,cdadc,cdadc,cdadc,cdadc,cdadc,\
                         cscda,cpcda,cpcda,cpcda,cstarcda,cdadc,cdadc,cdadc,cdadc,cdadc,\
                         cscda,cpcda,cpcda,cpcda,cstarcda,cdadc,cdadc,cdadc,cdadc,cdadc,\
                         cscda,cpcda,cpcda,cpcda,cstarcda,cdadc,cdadc,cdadc,cdadc,cdadc,\
                         cscda,cpcda,cpcda,cpcda,cstarcda,cdadc,cdadc,cdadc,cdadc,cdadc};
        c_dcopy(100,SC0,1,strain_const,1);
        break;}
    case 1:{
        double SC1[100]={cscsa,cscpa,cscpa,cscpa,cscstara,cscda,cscda,cscda,cscda,cscda,\
                         csapc,cpcpa,cpcpa,cpcpa,cstarapc,cpcda,cpcda,cpcda,cpcda,cpcda,\
                         csapc,cpcpa,cpcpa,cpcpa,cstarapc,cpcda,cpcda,cpcda,cpcda,cpcda,\
                         csapc,cpcpa,cpcpa,cpcpa,cstarapc,cpcda,cpcda,cpcda,cpcda,cpcda,\
                         csastarc,cstarcpa,cstarcpa,cstarcpa,cstarcstara,cstarcda,cstarcda,\
                         cstarcda,cstarcda,cstarcda,\
                         csadc,cpadc,cpadc,cpadc,cstaradc,cdcda,cdcda,cdcda,cdcda,cdcda,\
                         csadc,cpadc,cpadc,cpadc,cstaradc,cdcda,cdcda,cdcda,cdcda,cdcda,\
                         csadc,cpadc,cpadc,cpadc,cstaradc,cdcda,cdcda,cdcda,cdcda,cdcda,\
                         csadc,cpadc,cpadc,cpadc,cstaradc,cdcda,cdcda,cdcda,cdcda,cdcda,\
                         csadc,cpadc,cpadc,cpadc,cstaradc,cdcda,cdcda,cdcda,cdcda,cdcda};
        c_dcopy(100,SC1,1,strain_const,1);
        break;}
    }
}

/************************************************************************************************/

void Material::get_strain_shift(double *strain_shift,double *Vns,double *Vs,double *bond0,\
				double *bond,int fst_entry,int snd_entry,int tb)
{
    int i,j,k,trd_entry=(snd_entry+1)%2;
    double strain_const[100],s_act,K1,K2;
    double Eda1   = (Eda12[fst_entry]+2.0*Eda15[fst_entry])/3.0;
    double Eda2   = (2.0*Eda12[fst_entry]+Eda15[fst_entry])/3.0;
    double Edc1   = (Edc12[fst_entry]+2.0*Edc15[fst_entry])/3.0;
    double Edc2   = (2.0*Edc12[fst_entry]+Edc15[fst_entry])/3.0;
    double P[20]  = {Esa[fst_entry],Epa[fst_entry],Epa[fst_entry],Epa[fst_entry],\
		     Estara[fst_entry],Eda2,Eda1,Eda1,Eda2,Eda15[fst_entry],\
		     Esc[fst_entry],Epc[fst_entry],Epc[fst_entry],Epc[fst_entry],\
		     Estarc[fst_entry],Edc2,Edc1,Edc1,Edc2,Edc15[fst_entry]};
    double eshift = Eshift[fst_entry];
    double V0     = 6.1244087;

    prepare_strain_constants(strain_const,fst_entry,snd_entry);
    
    init_var(strain_shift,TB*TB);
    
    if(!strain_model){

        for(i=0;i<TB;i++){
    
	    s_act = 0.0;

	    for(k=0;k<TB;k++){

	        K1    = -1.0+sqrt(1.0+2.0*strain_const[i*TB+k]);
		K2    = -1.0+sqrt(1.0+2.0*strain_const[i*TB+k]);

		s_act = s_act+(Vns[i*TB+k]*Vns[i*TB+k]-Vs[i*TB+k]*Vs[i*TB+k])/2.0*\
		  (K1/(P[snd_entry*TB+i]+P[trd_entry*TB+k]-2*eshift)+\
		   K2/(P[snd_entry*TB+i]+P[trd_entry*TB+k]-2*eshift)+K1*K2/2.0*\
		   (1/(P[snd_entry*TB+i]+P[trd_entry*TB+k]-2*eshift)+
		   1/(P[snd_entry*TB+i]+P[trd_entry*TB+k]-2*eshift)));							       

	    }

	    strain_shift[i*TB+i] = s_act;
	}

    }else{

        for(i=0;i<TB;i++){

	    for(j=0;j<TB;j++){
    
	        s_act = 0.0;

		for(k=0;k<TB;k++){

		    K1    = -1.0+sqrt(1.0+2.0*strain_const[i*TB+k]);
		    K2    = -1.0+sqrt(1.0+2.0*strain_const[j*TB+k]);

		    s_act = s_act+(Vns[i*TB+k]*Vns[j*TB+k]-Vs[i*TB+k]*Vs[j*TB+k])/2.0* \
		      (K1/(P[snd_entry*TB+i]+P[trd_entry*TB+k]-2*eshift)+ \
		       K2/(P[snd_entry*TB+j]+P[trd_entry*TB+k]-2*eshift)+K1*K2/2.0* \
		       (1/(P[snd_entry*TB+i]+P[trd_entry*TB+k]-2*eshift)+
			1/(P[snd_entry*TB+j]+P[trd_entry*TB+k]-2*eshift)));
									       
		}

		strain_shift[i*TB+j] = s_act;
	    }
	}
    }
    
    get_Hdd_shift(strain_shift,bond0,-1.0,c_dnrm2(N3D,bond0,1),r2av_dd[fst_entry],\
		  r4av_dd[fst_entry],Zeff[fst_entry]*V0);
    get_Hdd_shift(strain_shift,bond,1.0,c_dnrm2(N3D,bond0,1),r2av_dd[fst_entry],\
		  r4av_dd[fst_entry],Zeff[fst_entry]*V0);
    
    get_Hpp_shift(strain_shift,bond0,-1.0,c_dnrm2(N3D,bond0,1),r2av_pp[fst_entry],\
		  Zeff[fst_entry]*V0);
    get_Hpp_shift(strain_shift,bond,1.0,c_dnrm2(N3D,bond0,1),r2av_pp[fst_entry],\
		  Zeff[fst_entry]*V0);
    
    get_Hpd_shift(strain_shift,bond0,-1.0,c_dnrm2(N3D,bond0,1),r_pd[fst_entry],\
		  r3_pd[fst_entry],Zeff[fst_entry]*V0);
    get_Hpd_shift(strain_shift,bond,1.0,c_dnrm2(N3D,bond0,1),r_pd[fst_entry],\
		  r3_pd[fst_entry],Zeff[fst_entry]*V0);

}

/************************************************************************************************/

void Material::get_Hdd_shift(double *strain_shift,double *bond,double ssign,double d0,\
			     double r2av0,double r4av0,double Zeff_V0)
{
  
    double x,y,z,x2,y2,z2,r2,rabs,l,m,n;
    double d,r2av,r4av,V0;
    double V_dd_s,V_dd_p,V_dd_d;
    double H1,H2,H3,H4,H5;
    	
    //Get building blocks
    x       = bond[0];
    y       = bond[1];
    z       = bond[2];
    x2      = x*x;
    y2      = y*y;
    z2      = z*z;
    r2      = x2 + y2 + z2;
    rabs    = sqrt(r2);
    l       = x/rabs;
    m       = y/rabs;
    n       = z/rabs; 
    d       = c_dnrm2(N3D,bond,1);
 
    //Scale r2av0, r4av0, V0u
    r2av  = r2av0*pow(d0/d,2.0);
    r4av  = r4av0*pow(d0/d,4.0);
    V0    = Zeff_V0*d0/d;
			
    //Effective Slater-Koster parameters
    V_dd_s  = -V0*(r2av+r4av)*2.0/7.0;
    V_dd_p  = -V0*(r2av-4.0*r4av/3.0)/7.0;
    V_dd_d  = -V0*(-2.0*r2av+r4av/3.0)/7.0;

    //XY
    H1      = 3.0*pow(l,2.0)*pow(m,2.0)*V_dd_s+\
      (pow(l,2.0)+pow(m,2.0)-4.0*pow(l,2.0)*pow(m,2.0))*V_dd_p+\
      (pow(n,2.0)+pow(l,2.0)*pow(m,2.0))*V_dd_d;
    H2      = 3.0*l*pow(m,2.0)*n*V_dd_s+l*n*(1-4.0*pow(m,2.0))*V_dd_p+l*n*(pow(m,2.0)-1)*V_dd_d;
    H3      = 3.0*pow(l,2.0)*m*n*V_dd_s+m*n*(1-4.0*pow(l,2.0))*V_dd_p+m*n*(pow(l,2.0)-1)*V_dd_d;
    H4      = 3.0/2.0*l*m*(pow(l,2.0)-pow(m,2.0))*V_dd_s+\
      2.0*l*m*(pow(m,2.0)-pow(l,2.0))*V_dd_p+l*m*(pow(l,2.0)-pow(m,2.0))/2.0*V_dd_d;
    H5      = sqrt(3.0)*l*m*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*V_dd_s-\
      2.0*sqrt(3.0)*l*m*pow(n,2.0)*V_dd_p+sqrt(3.0)/2.0*m*l*(1+pow(n,2.0))*V_dd_d;

    strain_shift[ID1*TB+ID1] = strain_shift[ID1*TB+ID1]+ssign*H1;
    strain_shift[ID1*TB+ID2] = strain_shift[ID1*TB+ID2]+ssign*H2;
    strain_shift[ID1*TB+ID3] = strain_shift[ID1*TB+ID3]+ssign*H3;
    strain_shift[ID1*TB+ID4] = strain_shift[ID1*TB+ID4]+ssign*H4;
    strain_shift[ID1*TB+ID5] = strain_shift[ID1*TB+ID5]+ssign*H5;
    strain_shift[ID2*TB+ID1] = strain_shift[ID2*TB+ID1]+ssign*H2;
    strain_shift[ID3*TB+ID1] = strain_shift[ID3*TB+ID1]+ssign*H3;
    strain_shift[ID4*TB+ID1] = strain_shift[ID4*TB+ID1]+ssign*H4;
    strain_shift[ID5*TB+ID1] = strain_shift[ID5*TB+ID1]+ssign*H5;

    //YZ
    H2      = 3.0*pow(n,2.0)*pow(m,2.0)*V_dd_s+\
      (pow(n,2.0)+pow(m,2.0)-4.0*pow(n,2.0)*pow(m,2.0))*V_dd_p+\
      (pow(l,2.0)+pow(n,2.0)*pow(m,2.0))*V_dd_d;
    H3      = 3.0*pow(n,2.0)*l*m*V_dd_s+m*l*(1-4.0*pow(n,2.0))*V_dd_p+m*l*(pow(n,2.0)-1)*V_dd_d;
    H4      = 3.0/2.0*m*n*(pow(l,2.0)-pow(m,2.0))*V_dd_s-\
      m*n*(1+2.0*(pow(l,2.0)-pow(m,2.0)))*V_dd_p+\
      m*n*(1+(pow(l,2.0)-pow(m,2.0))/2.0)*V_dd_d;
    H5      = sqrt(3.0)*m*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*V_dd_s+\
      sqrt(3.0)*m*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*V_dd_p-\
      sqrt(3.0)/2.0*m*n*(pow(l,2.0)+pow(m,2.0))*V_dd_d;

    strain_shift[ID2*TB+ID2] = strain_shift[ID2*TB+ID2]+ssign*H2;
    strain_shift[ID2*TB+ID3] = strain_shift[ID2*TB+ID3]+ssign*H3;
    strain_shift[ID2*TB+ID4] = strain_shift[ID2*TB+ID4]+ssign*H4;
    strain_shift[ID2*TB+ID5] = strain_shift[ID2*TB+ID5]+ssign*H5;
    strain_shift[ID3*TB+ID2] = strain_shift[ID3*TB+ID2]+ssign*H3;
    strain_shift[ID4*TB+ID2] = strain_shift[ID4*TB+ID2]+ssign*H4;
    strain_shift[ID5*TB+ID2] = strain_shift[ID5*TB+ID2]+ssign*H5;

    //ZX
    H3      = 3.0*pow(n,2.0)*pow(l,2.0)*V_dd_s+\
      (pow(n,2.0)+pow(l,2.0)-4.0*pow(n,2.0)*pow(l,2.0))*V_dd_p+\
      (pow(m,2.0)+pow(n,2.0)*pow(l,2.0))*V_dd_d;
    H4      = 3.0/2.0*n*l*(pow(l,2.0)-pow(m,2.0))*V_dd_s+\
      n*l*(1-2.0*(pow(l,2.0)-pow(m,2.0)))*V_dd_p-\
      n*l*(1-(pow(l,2.0)-pow(m,2.0))/2.0)*V_dd_d;
    H5      = sqrt(3.0)*l*n*(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*V_dd_s+\
      sqrt(3.0)*l*n*(pow(l,2.0)+pow(m,2.0)-pow(n,2.0))*V_dd_p-\
      sqrt(3.0)/2.0*l*n*(pow(l,2.0)+pow(m,2.0))*V_dd_d;

    strain_shift[ID3*TB+ID3] = strain_shift[ID3*TB+ID3]+ssign*H3;
    strain_shift[ID3*TB+ID4] = strain_shift[ID3*TB+ID4]+ssign*H4;
    strain_shift[ID3*TB+ID5] = strain_shift[ID3*TB+ID5]+ssign*H5;
    strain_shift[ID4*TB+ID3] = strain_shift[ID4*TB+ID3]+ssign*H4;
    strain_shift[ID5*TB+ID3] = strain_shift[ID5*TB+ID3]+ssign*H5;

    //X2_Y2
    H4      = 3.0/4.0*pow(pow(l,2.0)-pow(m,2.0),2.0)*V_dd_s+\
      (pow(l,2.0)+pow(m,2.0)-pow(pow(l,2.0)-pow(m,2.0),2.0))*V_dd_p+\
      (pow(n,2.0)+pow(pow(l,2.0)-pow(m,2.0),2.0)/4.0)*V_dd_d;
    H5      = sqrt(3.0)/2.0*(pow(l,2.0)-pow(m,2.0))*\
      (pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0)*V_dd_s+\
      sqrt(3.0)*pow(n,2.0)*(pow(m,2.0)-pow(l,2.0))*V_dd_p+\
      sqrt(3.0)/4.0*(1+pow(n,2.0))*(pow(l,2.0)-pow(m,2.0))*V_dd_d;

    strain_shift[ID4*TB+ID4] = strain_shift[ID4*TB+ID4]+ssign*H4;
    strain_shift[ID4*TB+ID5] = strain_shift[ID4*TB+ID5]+ssign*H5;
    strain_shift[ID5*TB+ID4] = strain_shift[ID5*TB+ID4]+ssign*H5;

    //3Z2_R2
    H5      = pow(pow(n,2.0)-(pow(l,2.0)+pow(m,2.0))/2.0,2.0)*V_dd_s+\
      3.0*pow(n,2.0)*(pow(l,2.0)+pow(m,2.0))*V_dd_p+\
	  3.0/4.0*pow(pow(l,2.0)+pow(m,2.0),2.0)*V_dd_d;

    strain_shift[ID5*TB+ID5] = strain_shift[ID5*TB+ID5]+ssign*H5;
  
  /*
    double Y20,Y21_s,Y21_d,Y22_s,Y22_d,Y40,Y41_s,Y41_d,Y42_s;
    double Y42_d,Y43_s,Y43_d,Y44_s,Y44_d,Y2v[5],Y4v[9];
    double d,Y00,r2av,r4av,V0;
    double H1,H2,H3,H4,H5;
	
    Y00   = 0.0;
	
    d     = c_dnrm2(N3D,bond,1);
			
    //Get spherical harmonics for this neighbor
    Y_2_4_combos(bond,Y2v,Y4v);
		
    //Assign spherical harmonic functions
    Y20   = Y2v[0];
    Y21_s = Y2v[1];
    Y21_d = Y2v[2];
    Y22_s = Y2v[3];
    Y22_d = Y2v[4];
		
    Y40   = Y4v[0];
    Y41_s = Y4v[1];
    Y41_d = Y4v[2];
    Y42_s = Y4v[3];
    Y42_d = Y4v[4];
    Y43_s = Y4v[5];
    Y43_d = Y4v[6];
    Y44_s = Y4v[7];
    Y44_d = Y4v[8];
		
    //Scale r2av0, r4av0, V0u
    r2av  = r2av0*pow(d0/d,2.0);
    r4av  = r4av0*pow(d0/d,4.0);
    V0    = Zeff_V0*d0/d;

    //XY Row
    H1    = V0*sqrt(PI)*(-2.0*Y00 + 4.0*r2av*Y20*sqrt(5.0)/35.0 +
			 (r4av/63.0)*(-2.0*Y40 + sqrt(70.0)*Y44_s));
    H2    = V0*sqrt(PI)*(sqrt(30.0)*r2av*Y21_d/35.0 - (r4av/63.0)*
			 (sqrt(5.0)*Y41_d + sqrt(35.0)*Y43_d));
    H3    = V0*sqrt(PI)*(sqrt(30.0)*r2av*Y21_s/35.0 - (r4av/63.0)*
			 (sqrt(5.0)*Y41_s - sqrt(35.0)*Y43_s));
    H4    = -V0*sqrt(70.0)*sqrt(PI)*r4av*Y44_d/63.0;
    H5    = V0*sqrt(PI)*(2.0*sqrt(10.0)*r2av*Y22_d/35.0 - 
			 sqrt(30.0)*r4av*Y42_d/63.0);

    strain_shift[ID1*TB+ID1] = strain_shift[ID1*TB+ID1]+ssign*H1;
    strain_shift[ID1*TB+ID2] = strain_shift[ID1*TB+ID2]+ssign*H2;
    strain_shift[ID1*TB+ID3] = strain_shift[ID1*TB+ID3]+ssign*H3;
    strain_shift[ID1*TB+ID4] = strain_shift[ID1*TB+ID4]+ssign*H4;
    strain_shift[ID1*TB+ID5] = strain_shift[ID1*TB+ID5]+ssign*H5;
    strain_shift[ID2*TB+ID1] = strain_shift[ID2*TB+ID1]+ssign*H2;
    strain_shift[ID3*TB+ID1] = strain_shift[ID3*TB+ID1]+ssign*H3;
    strain_shift[ID4*TB+ID1] = strain_shift[ID4*TB+ID1]+ssign*H4;
    strain_shift[ID5*TB+ID1] = strain_shift[ID5*TB+ID1]+ssign*H5;

    H2    = V0*sqrt(PI)*(-2.0*Y00 + r2av*(-2.0*Y20*sqrt(5.0) +
				       sqrt(30.0)*Y22_s)/35.0 + r4av*(8.0*Y40 + 
								      sqrt(40.0)*Y42_s)/63.0);
    H3    = V0*sqrt(PI)*(-sqrt(30.0)*r2av*Y22_d/35.0 - 
			 r4av*sqrt(40.0)*Y42_d/63.0);
    H4    = V0*sqrt(PI)*(-sqrt(30.0)*r2av*Y21_s/35.0 +
		      r4av*(sqrt(5.0)*Y41_s + sqrt(35.0)*Y43_s)/63.0);
    H5    = V0*sqrt(PI)*(sqrt(10.0)*r2av*Y21_s/35.0 + 
		      2.0*sqrt(15.0)*r4av*Y41_s/63.0);

    strain_shift[ID2*TB+ID2] = strain_shift[ID2*TB+ID2]+ssign*H2;
    strain_shift[ID2*TB+ID3] = strain_shift[ID2*TB+ID3]+ssign*H3;
    strain_shift[ID2*TB+ID4] = strain_shift[ID2*TB+ID4]+ssign*H4;
    strain_shift[ID2*TB+ID5] = strain_shift[ID2*TB+ID5]+ssign*H5;
    strain_shift[ID3*TB+ID2] = strain_shift[ID3*TB+ID2]+ssign*H3;
    strain_shift[ID4*TB+ID2] = strain_shift[ID4*TB+ID2]+ssign*H4;
    strain_shift[ID5*TB+ID2] = strain_shift[ID5*TB+ID2]+ssign*H5;


    H3    = V0*sqrt(PI)*(-2.0*Y00 + r2av*(-2.0*Y20*sqrt(5.0) -
					  sqrt(30.0)*Y22_s)/35.0 + r4av*(8.0*Y40 - 
									 sqrt(40.0)*Y42_s)/63.0);
    H4    = V0*sqrt(PI)*(sqrt(30.0)*r2av*Y21_d/35.0 +
			 r4av*(-sqrt(5.0)*Y41_d + sqrt(35.0)*Y43_d)/63.0);
    H5    = V0*sqrt(PI)*(sqrt(10.0)*r2av*Y21_d/35.0 + 
			 2.0*sqrt(15.0)*r4av*Y41_d/63.0);
  
    strain_shift[ID3*TB+ID3] = strain_shift[ID3*TB+ID3]+ssign*H3;
    strain_shift[ID3*TB+ID4] = strain_shift[ID3*TB+ID4]+ssign*H4;
    strain_shift[ID3*TB+ID5] = strain_shift[ID3*TB+ID5]+ssign*H5;
    strain_shift[ID4*TB+ID3] = strain_shift[ID4*TB+ID3]+ssign*H4;
    strain_shift[ID5*TB+ID3] = strain_shift[ID5*TB+ID3]+ssign*H5;

    H4    = V0*sqrt(PI)*(-2.0*Y00 + 4.0*r2av*sqrt(5.0)*Y20/35.0
			 -r4av*(2.0*Y40 + sqrt(70.0)*Y44_s)/63.0);
    H5    = V0*sqrt(PI)*(2.0*sqrt(10.0)*r2av*Y22_s/35.0 - 
			 sqrt(30.0)*r4av*Y42_s/63.0);

    strain_shift[ID4*TB+ID4] = strain_shift[ID4*TB+ID4]+ssign*H4;
    strain_shift[ID4*TB+ID5] = strain_shift[ID4*TB+ID5]+ssign*H5;
    strain_shift[ID5*TB+ID4] = strain_shift[ID5*TB+ID4]+ssign*H5;
	

    H5    = -V0*sqrt(PI)*(2.0*Y00 + 4.0*sqrt(5.0)*Y20*r2av/35.0 +
			  4.0*Y40*r4av/21.0);

    strain_shift[ID5*TB+ID5] = strain_shift[ID5*TB+ID5]+ssign*H5;    
  */
}

/************************************************************************************************/

void Material::get_Hpp_shift(double *strain_shift,double *bond,double ssign,double d0,\
			     double r2av0,double Zeff_V0)
{

    double Y20,Y21_s,Y21_d,Y22_s,Y22_d,Y2v[5],Y4v[9];
    double Y00,d,r2av,V0;
    double H1,H2,H3;
	
    //Y00 = 0.5/srtpi;
    //Disable universal diagonal shift
    Y00   = 0.0;

    d     = c_dnrm2(N3D,bond,1);
			
    //Get spherical harmonics for this neighbor
    Y_2_4_combos(bond,Y2v,Y4v);
	
    /* Assign spherical harmonic functions; Y4n_s, Y4_d unused */
    Y20   = Y2v[0];
    Y21_s = Y2v[1];
    Y21_d = Y2v[2];
    Y22_s = Y2v[3];
    Y22_d = Y2v[4];
		
    /* Scale r2av0, V0u; unused */
    r2av  = r2av0*pow(d0/d,2.0);
    V0    = Zeff_V0*d0/d;
	
    H1    = V0*sqrt(PI)*(-2.0*Y00 + (2.0*Y20*sqrt(5.0) -
				     sqrt(30.0)*Y22_s)*r2av/25.0);
    H2    = -V0*sqrt(PI)*sqrt(30.0)*r2av*Y22_d/25.0;
    H3    = V0*sqrt(PI)*sqrt(30.0)*r2av*Y21_d/25.0;
    
    strain_shift[IPX*TB+IPX] = strain_shift[IPX*TB+IPX]+ssign*H1;
    strain_shift[IPX*TB+IPY] = strain_shift[IPX*TB+IPY]+ssign*H2;
    strain_shift[IPX*TB+IPZ] = strain_shift[IPX*TB+IPZ]+ssign*H3;
    strain_shift[IPY*TB+IPX] = strain_shift[IPY*TB+IPX]+ssign*H2;
    strain_shift[IPZ*TB+IPX] = strain_shift[IPZ*TB+IPX]+ssign*H3;

    H2    = V0*sqrt(PI)*(-2.0*Y00 + (2.0*Y20*sqrt(5.0) +
				     sqrt(30.0)*Y22_s)*r2av/25.0);
    H3    = V0*sqrt(PI)*sqrt(30.0)*r2av*Y21_s/25.0;
	
    strain_shift[IPY*TB+IPY] = strain_shift[IPY*TB+IPY]+ssign*H2;
    strain_shift[IPY*TB+IPZ] = strain_shift[IPY*TB+IPZ]+ssign*H3;
    strain_shift[IPZ*TB+IPY] = strain_shift[IPZ*TB+IPY]+ssign*H3;
    

    H3    = -V0*sqrt(PI)*(2.0*Y00 + 4.0*r2av*Y20*sqrt(5.0)/25.0);	

    strain_shift[IPZ*TB+IPZ] = strain_shift[IPZ*TB+IPZ]+ssign*H3;
}

/************************************************************************************************/

void Material::get_Hpd_shift(double *strain_shift,double *bond,double ssign,double d0,\
			     double r_pd0,double r3_pd0,double Zeff_V0)
{

    double x,y,z,x2,y2,z2,r2,rabs,r3,l,m,n,l2,m2,n2,lmn;
    double r_pd,r3_pd,V0,temp,temp2;
    double l2_m_m2,l2_p_m2,V_pd_s,V_pd_p;
    double H1,H2,H3;
	
    //Get building blocks
    x       = bond[0];
    y       = bond[1];
    z       = bond[2];
    x2      = x*x;
    y2      = y*y;
    z2      = z*z;
    r2      = x2 + y2 + z2;
    rabs    = sqrt(r2);
    r3      = rabs*r2;
    l       = x/rabs;
    m       = y/rabs;
    n       = z/rabs;
    l2      = l*l;
    m2      = m*m;
    n2      = n*n;
    l2_p_m2 = l2 + m2;
    l2_m_m2 = l2 - m2;
    lmn     = l*m*n;
 
    //Scale parameters r_pd, r3_pd 
    r_pd    = r_pd0*d0/rabs;
    r3_pd   = r3_pd0*(d0*d0*d0)/r3;
    V0      = Zeff_V0*d0/rabs;
			
    //Effective Slater-Koster parameters
    V_pd_s  = -V0*(2.0*r_pd + 9.0*r3_pd/7.0)/sqrt(15.0);
    V_pd_p  = -V0*(r_pd - 3.0*r3_pd/7.0)/sqrt(5.0);
		
    /* Add contributions into Hamiltonian */
    temp    = 1.0 - 2.0*l2;
    H1      = sqrt(3.0)*l2*m*V_pd_s + m*temp*V_pd_p;
    H2      = lmn*(sqrt(3.0)*V_pd_s - 2.0*V_pd_p);
    H3      = sqrt(3.0)*l2*n*V_pd_s + n*temp*V_pd_p;

    strain_shift[IPX*TB+ID1] = strain_shift[IPX*TB+ID1]+ssign*H1;
    strain_shift[IPX*TB+ID2] = strain_shift[IPX*TB+ID2]+ssign*H2;
    strain_shift[IPX*TB+ID3] = strain_shift[IPX*TB+ID3]+ssign*H3;
    strain_shift[ID1*TB+IPX] = strain_shift[ID1*TB+IPX]+ssign*H1;
    strain_shift[ID2*TB+IPX] = strain_shift[ID2*TB+IPX]+ssign*H2;
    strain_shift[ID3*TB+IPX] = strain_shift[ID3*TB+IPX]+ssign*H3;

    temp    = 1.0 - 2.0*m2;
    H1      = sqrt(3.0)*m2*l*V_pd_s + l*temp*V_pd_p;
    H2      = sqrt(3.0)*m2*n*V_pd_s + n*temp*V_pd_p;
    H3      = lmn*(sqrt(3.0)*V_pd_s - 2.0*V_pd_p);

    strain_shift[IPY*TB+ID1] = strain_shift[IPY*TB+ID1]+ssign*H1;
    strain_shift[IPY*TB+ID2] = strain_shift[IPY*TB+ID2]+ssign*H2;
    strain_shift[IPY*TB+ID3] = strain_shift[IPY*TB+ID3]+ssign*H3;
    strain_shift[ID1*TB+IPY] = strain_shift[ID1*TB+IPY]+ssign*H1;
    strain_shift[ID2*TB+IPY] = strain_shift[ID2*TB+IPY]+ssign*H2;
    strain_shift[ID3*TB+IPY] = strain_shift[ID3*TB+IPY]+ssign*H3;

    temp    = 1.0 - 2.0*n2;
    H1      = lmn*(sqrt(3.0)*V_pd_s - 2.0*V_pd_p);
    H2      = sqrt(3.0)*n2*m*V_pd_s + m*temp*V_pd_p;
    H3      = sqrt(3.0)*n2*l*V_pd_s + l*temp*V_pd_p;

    strain_shift[IPZ*TB+ID1] = strain_shift[IPZ*TB+ID1]+ssign*H1;
    strain_shift[IPZ*TB+ID2] = strain_shift[IPZ*TB+ID2]+ssign*H2;
    strain_shift[IPZ*TB+ID3] = strain_shift[IPZ*TB+ID3]+ssign*H3;
    strain_shift[ID1*TB+IPZ] = strain_shift[ID1*TB+IPZ]+ssign*H1;
    strain_shift[ID2*TB+IPZ] = strain_shift[ID2*TB+IPZ]+ssign*H2;
    strain_shift[ID3*TB+IPZ] = strain_shift[ID3*TB+IPZ]+ssign*H3;

    temp    = 0.5*sqrt(3.0)*l2_m_m2*V_pd_s;
    H1      = l*(temp + (1.0 - l2_m_m2)*V_pd_p);
    H2      = m*(temp - (1.0 + l2_m_m2)*V_pd_p);
    H3      = n*(temp - l2_m_m2*V_pd_p);

    strain_shift[IPX*TB+ID4] = strain_shift[IPX*TB+ID4]+ssign*H1;
    strain_shift[IPY*TB+ID4] = strain_shift[IPY*TB+ID4]+ssign*H2;
    strain_shift[IPZ*TB+ID4] = strain_shift[IPZ*TB+ID4]+ssign*H3;
    strain_shift[ID4*TB+IPX] = strain_shift[ID4*TB+IPX]+ssign*H1;
    strain_shift[ID4*TB+IPY] = strain_shift[ID4*TB+IPY]+ssign*H2;
    strain_shift[ID4*TB+IPZ] = strain_shift[ID4*TB+IPZ]+ssign*H3;
    
    temp    = (n2 - 0.5*l2_p_m2)*V_pd_s;
    temp2   = temp - n2*sqrt(3.0)*V_pd_p;
    H1      = l*temp2;
    H2      = m*temp2;
    H3      = n*(temp + sqrt(3.0)*l2_p_m2*V_pd_p);

    strain_shift[IPX*TB+ID5] = strain_shift[IPX*TB+ID5]+ssign*H1;
    strain_shift[IPY*TB+ID5] = strain_shift[IPY*TB+ID5]+ssign*H2;
    strain_shift[IPZ*TB+ID5] = strain_shift[IPZ*TB+ID5]+ssign*H3;
    strain_shift[ID5*TB+IPX] = strain_shift[ID5*TB+IPX]+ssign*H1;
    strain_shift[ID5*TB+IPY] = strain_shift[ID5*TB+IPY]+ssign*H2;
    strain_shift[ID5*TB+IPZ] = strain_shift[ID5*TB+IPZ]+ssign*H3;

}

/************************************************************************************************/

void Material::Y_2_4_combos(double *bond,double *Y2,double *Y4)
{

    double x,y,z,xy,yz,zx,x2_y2,r2,r4,pr2,pr4,x2,y2,z2;
    double sz2_r2,sz2_3r2,s5dp,temp;

    /* Get building blocks */
    x       = bond[0];
    y       = bond[1];
    z       = bond[2];
    x2      = x*x;
    y2      = y*y;
    z2      = z*z;
    r2      = x2+y2+z2;
    r4      = r2*r2;
    xy      = x*y;
    yz      = y*z;
    zx      = z*x;
    x2_y2   = x2-y2;

    /* Compute Y2 quantities */
    pr2     = sqrt(0.5*15.0/PI);
    Y2[0]   = 0.25*sqrt(5.0/PI)*(3.0*z2 - r2)/r2;
    Y2[1]   = -pr2*yz/r2;
    Y2[2]   = -pr2*zx/r2;
    Y2[3]   = 0.5*pr2*x2_y2/r2;
    Y2[4]   = pr2*xy/r2;

    /* Compute Y4 quantities */
    temp    = 0.1875/sqrt(PI);
    Y4[0]   = temp*(35.0*z2*z2 - 30.0*z2*r2 + 3.0*r4)/r4;
    s5dp    = sqrt(5.0/PI);
    pr4     = 0.75*s5dp;
    sz2_r2  = 7.0*z2 - r2;
    sz2_3r2 = sz2_r2 - 2.0*r2;

    Y4[1]   = -pr4*sz2_3r2*yz/r4;
    Y4[2]   = -pr4*sz2_3r2*zx/r4;

    /* 3*sqrt(5/8PI) */
    pr4     = 1.5*s5dp/sqrt(2.0);
    Y4[3]   = 0.5*pr4*sz2_r2*x2_y2/r4;
    Y4[4]   = pr4*sz2_r2*xy/r4;

    pr4     = 1.5*sqrt(35.0/PI);
    Y4[5]   = -0.5*pr4*(3.0*x2 - y2)*yz/r4;
    Y4[6]   = -0.5*pr4*(x2 - 3.0*y2)*zx/r4;

    pr4     = pr4/sqrt(2.0);
    Y4[7]   = 0.25*pr4*(x2_y2*x2_y2 - 4.0*x2*y2)/r4;
    Y4[8]   = pr4*x2_y2*xy/r4;

}

/************************************************************************************************/
