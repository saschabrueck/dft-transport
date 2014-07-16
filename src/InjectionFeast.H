#ifndef __INJECTIONFEAST
#define __INJECTIONFEAST

#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Types.H"
#include "Utilities.H"
#include "Blas.H"
#include "Injection.H"

template<class T>
class InjectionFeast : public Injection<T>{
public:
    
    InjectionFeast(int,int,double);

    virtual ~InjectionFeast();

    virtual void initialize(int,int,int);
    virtual void calc_kphase(TCSR<T>*,TCSR<T>*,TCSR<CPX>*,int,int,CPX*,CPX*,double*,\
			     int*,int*,int*,int*,int,int,int,MPI_Comm,int*);
    virtual void calc_kphase_sym(TCSR<T>*,TCSR<T>*,TCSR<CPX>*,int,int,CPX*,CPX*,CPX*,\
				 CPX*,double*,double*,int*,int*,int*,int*,int*,int*,\
				 int*,int*,int,int,MPI_Comm,int*);
    virtual void calc_kphase(TCSR<T>*,T*,TCSR<CPX>*,int,int,CPX*,CPX*,double*,int*,\
			     int*,int*,int*,int,int,int,MPI_Comm,int*);
    virtual void calc_kphase_sym(TCSR<T>*,T*,TCSR<CPX>*,int,int,CPX*,CPX*,CPX*,CPX*,\
				 double*,double*,int*,int*,int*,int*,int*,int*,int*,\
				 int*,int,int,MPI_Comm,int*);
    
private:

    double error_criterion;
    double mult_factor;
    double feast_radius;
    double imag_limit;
    double *alphaF,*dalphaF;
    double *theta,*dtheta;
    double *rwork;
    int mpi_size,mpi_rank;
    int loc_size,loc_rank;
    int root_rank;
    int factor_lwork;
    int contour_factor;
    int init_f_matrix;
    int NM_for_feast;
    int max_iteration;
    int slab_per_bc;
    int NQDR;
    int *index_to_keep;
    CPX *AR,*BR,*CR;
    CPX *z_out,*dz_dtheta_out;
    CPX *z_in,*dz_dtheta_in;
    CPX *lambda_up,*lambda_down;
    CPX *lambdaR,*VR;
    CPX *work;

    void sort_eig(CPX*,CPX*,double*,int*,int*,int*,int*,TCSR<CPX>*,int,int,int);
    int get_type(TCSR<CPX>*,CPX*,CPX,double*,int,int);
    double calc_derivative(TCSR<CPX>*,CPX*,CPX);
    void init_matrix_single(TCSR<CPX>*,TCSR<CPX>*,int*,TCSR<T>*,TCSR<T>*,int,int);
    void init_matrix_double(TCSR<CPX>*,TCSR<CPX>*,int*,TCSR<T>*,TCSR<T>*,int,int);
    void GaussLegendre(double*,double*,int);
    void eig(double*,double*,int);
    int gen_eig(CPX*,CPX*,int,CPX*,CPX*);
    int gen_eig_real(CPX*,CPX*,int,CPX*,CPX*);
    void get_Q_matrix(CPX*,CPX*,CPX*,TCSR<CPX>*,TCSR<CPX>*,int*,int,int,MPI_Comm);
    void solve_Q_sparse(CPX*,CPX*,TCSR<CPX>*,TCSR<CPX>*,int*,int,int,CPX,CPX);
    void solve_Q_full(CPX*,CPX*,TCSR<CPX>*,TCSR<CPX>*,int*,int,int,CPX,CPX);
    void init_Y_matrix(CPX*,int,int);
    void add_sparse(TCSR<CPX>*,TCSR<CPX>*,int*,CPX factor);
    void transpose(CPX*,CPX*,int,int);
    void merge_A_and_B(TCSR<CPX>*,TCSR<CPX>*,CPX);
    void map_A_and_B(TCSR<CPX>*,TCSR<CPX>*,int*);
    void add_sparse_block_to_full(CPX*,TCSR<CPX>*,CPX,int);
    void GetIntegrationPoints(double,double);
};

template<class T>
InjectionFeast<T>::InjectionFeast(int cpxpoint,int pslab_per_bc,double pfeast_radius)
{
    error_criterion = 1.0e-3;
    imag_limit      = 5.0e-5;
    mult_factor     = 1.2;
    max_iteration   = 2;
    NQDR            = 16;
    factor_lwork    = 20;
    contour_factor  = 2-cpxpoint;
    init_f_matrix   = 0;

    slab_per_bc     = pslab_per_bc;
    feast_radius    = pfeast_radius;
    imag_limit      = slab_per_bc*imag_limit;

    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
}

/************************************************************************************************/

template <class T>
InjectionFeast<T>::~InjectionFeast()
{
    delete[] AR;
    delete[] BR;
    delete[] CR;
    delete[] work;
    delete[] rwork;
    delete[] lambda_up;
    delete[] lambda_down;
    delete[] lambdaR;
    delete[] VR;
    delete[] index_to_keep;
    delete[] alphaF;
    delete[] dalphaF;
    delete[] theta;
    delete[] dtheta;
    delete[] z_out;
    delete[] dz_dtheta_out;
    delete[] z_in;
    delete[] dz_dtheta_in;
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::initialize(int ND,int NP,int NM)
{

    NM_for_feast  = Round(1.5*NM);
 
    AR            = new CPX[NM_for_feast*NM_for_feast];
    BR            = new CPX[NM_for_feast*NM_for_feast];
    CR            = new CPX[NM_for_feast*ND];

    work          = new CPX[factor_lwork*NM_for_feast];
    rwork         = new double[8*NM_for_feast];

    lambda_up     = new CPX[NM_for_feast];
    lambda_down   = new CPX[NM_for_feast];
    lambdaR       = new CPX[NM_for_feast];
    VR            = new CPX[NM_for_feast*NM_for_feast];

    index_to_keep = new int[NM_for_feast];

    alphaF        = new double[NQDR];
    dalphaF       = new double[NQDR];

    theta         = new double[NQDR];
    dtheta        = new double[NQDR];

    z_out         = new CPX[NQDR];
    dz_dtheta_out = new CPX[NQDR];

    z_in          = new CPX[NQDR];
    dz_dtheta_in  = new CPX[NQDR];
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::sort_eig(CPX *kphase,CPX *Vsurf,double *dEk_dk,int *Ntr,int *ind_Ntr,\
				 int *Nref,int *ind_Nref,TCSR<CPX> *T01,int ND,int NM,int sign)
{

    int IP,IM;
    int iter;
    int type;

    *Ntr  = 0;
    *Nref = 0;    
    IP    = 0;

    for(IM=0;IM<NM;IM++){

        if(index_to_keep[IM]){

	    c_zcopy(ND,&CR[IM*ND],1,&Vsurf[IP*ND],1);

	    kphase[IP] = CPX(0.0,1.0)*log(lambdaR[IM]);

	    for(iter=0;iter<floor(slab_per_bc/2.0);iter++){
                kphase[IP] = 2.0*kphase[IP];
                if(real(kphase[IP])>PI){
                    kphase[IP] = kphase[IP]-CPX(2.0*PI,0.0);
                }
		if(real(kphase[IP])<-PI){
                    kphase[IP] = kphase[IP]+CPX(2.0*PI,0.0);
                }
            }
	  
	    type = get_type(T01,&Vsurf[IP*ND],kphase[IP],&dEk_dk[(*Ntr)],index_to_keep[IM],\
			    sign);

	    if((type>0)&&(dEk_dk[(*Ntr)]>0)){
	        ind_Ntr[(*Ntr)] = IP;
		(*Ntr)          = (*Ntr)+1;
	    }else if (type<0){
	        ind_Nref[(*Nref)] = IP;
		(*Nref)           = (*Nref)+1;
	    }
                
	    IP++;
	}
    }
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::calc_kphase(TCSR<T> *D,TCSR<T> *P,TCSR<CPX> *T01,int ND,int NP,CPX *kphase,\
				    CPX *Vsurf,double *dEk_dk,int *Ntr,int *ind_Ntr,int *Nref, \
				    int *ind_Nref,int sign,int NPROW,int NPCOL,MPI_Comm comm,\
				    int *info)
{
    int iteration;
    int IM,IN;
    int no_ev;
    int no_search;
    int cond1,cond2,cond3;//cond4;
    int m_size,n_nonzeros_A,n_nonzeros_B;
    int no_prop;
    int *map_index;
    double LBound,UBound;
    CPX kM,kN;
    TCSR<CPX> *AF, *BF;
    int NM            = NM_for_feast;
    int act_max_iter  = max_iteration;
    double act_radius = feast_radius;
    int *reord_index  = new int[NM];
    int *prop_index   = new int[NM];
    CPX *QF           = new CPX[ND*NM];
    CPX *YF           = new CPX[ND*NM];
    CPX *TrQF         = new CPX[NM*ND];

    MPI_Comm_size(comm,&loc_size);
    MPI_Comm_rank(comm,&loc_rank);

    root_rank = (init_f_matrix==1)*(sign==-1)*(loc_size-1);

    //GetIntegrationPoints(mult_factor*feast_radius,1.0/(mult_factor*feast_radius));
    
    if(sign>0){
        GetIntegrationPoints(1.0+mult_factor/feast_radius,1.0/(mult_factor*feast_radius));
	LBound = 1.0/feast_radius;
	UBound = 1.0+1.0/feast_radius;
    }else{
        GetIntegrationPoints(mult_factor*feast_radius,1.0/(1.0+mult_factor/feast_radius));
	LBound = 1.0/(1.0+1.0/feast_radius);
	UBound = feast_radius;
    }

    if(loc_rank==root_rank){
        if(init_f_matrix){
            AF        = new TCSR<CPX>(D->size,Round(D->n_nonzeros/2),D->findx);
            BF        = new TCSR<CPX>(D->size,Round(D->n_nonzeros/8),D->findx);
            map_index = new int[D->n_nonzeros];

            if(slab_per_bc==1){
                init_matrix_single(AF,BF,map_index,D,P,ND,NP);
            }else{
                init_matrix_double(AF,BF,map_index,D,P,ND,NP);
            }
	}else{
            AF        = new TCSR<CPX>(D->size,D->n_nonzeros,D->findx);
            BF        = new TCSR<CPX>(P->size,P->n_nonzeros,P->findx);
            map_index = new int[D->n_nonzeros];
            AF->copy_contain(D,1.0);
            BF->copy_contain(P,1.0);
        }
	init_Y_matrix(YF,NM,ND);
    }

    double t0 = get_time(0.0);

    if(loc_size>1){

        if(loc_rank==root_rank){
	    m_size       = AF->size;
	    n_nonzeros_A = AF->n_nonzeros;
	    n_nonzeros_B = BF->n_nonzeros;
	}

	MPI_Bcast(&m_size,1,MPI_INT,root_rank,comm);
	MPI_Bcast(&n_nonzeros_A,1,MPI_INT,root_rank,comm);
	MPI_Bcast(&n_nonzeros_B,1,MPI_INT,root_rank,comm);

	if(loc_rank!=root_rank){
	    AF        = new TCSR<CPX>(m_size,n_nonzeros_A,D->findx);
	    BF        = new TCSR<CPX>(m_size,n_nonzeros_B,D->findx);
	    map_index = new int[n_nonzeros_B];
	}
       
	MPI_Bcast_Sparse(AF,root_rank,comm);
	MPI_Bcast_Sparse(BF,root_rank,comm);
	MPI_Bcast(map_index,n_nonzeros_B,MPI_INT,root_rank,comm);
	MPI_Bcast(YF,NM*ND,MPI_DOUBLE_COMPLEX,root_rank,comm);
    }

    if(!mpi_rank) printf("Time to generate A and B: %e\n",get_time(t0));

    for(iteration=0;iteration<act_max_iter;iteration++){

        t0 = get_time(0.0);
    
        get_Q_matrix(QF,TrQF,YF,AF,BF,map_index,NM,ND,comm);

	if(!mpi_rank) printf("Time to generate Q: %e\n",get_time(t0));

	t0 = get_time(0.0);
	if(loc_rank==root_rank){

            omp_set_num_threads(omp_get_max_threads()*loc_size);

	    transpose(TrQF,QF,ND,NM);

	    AF->trans_mat_vec_mult_omp(TrQF,CR,NM,ND);
	    c_zgemm('C','N',NM,NM,ND,CPX(1.0,0.0),QF,ND,CR,ND,CPX(0.0,0.0),AR,NM);

	    BF->trans_mat_vec_mult_omp(TrQF,CR,NM,ND);
	    c_zgemm('C','N',NM,NM,ND,CPX(1.0,0.0),QF,ND,CR,ND,CPX(0.0,0.0),BR,NM);

            omp_set_num_threads(omp_get_max_threads()/loc_size);

            if(contour_factor==2){
                gen_eig_real(AR,BR,NM,VR,lambdaR);
            }else{
                gen_eig(AR,BR,NM,VR,lambdaR);
            }

            omp_set_num_threads(omp_get_max_threads()*loc_size);

	    sort_abs_vec(lambdaR,reord_index,NM);

	    init_var(index_to_keep,NM);
	    no_ev   = 0;
	    no_prop = 0;

            c_zgemm('N','N',ND,NM,NM,CPX(1.0,0.0),QF,ND,VR,NM,\
                    CPX(0.0,0.0),TrQF,ND);

	    for(IM=0;IM<NM;IM++){

                c_zcopy(ND,&TrQF[IM*ND],1,&CR[reord_index[IM]*ND],1);

//	        c_zgemm('N','N',ND,1,NM,CPX(1.0,0.0),QF,ND,&VR[IM*NM],NM,\
//			CPX(0.0,0.0),&CR[reord_index[IM]*ND],ND);
	    
		c_zscal(ND,CPX(1.0/c_dznrm2(ND,&CR[reord_index[IM]*ND],1),0.0),\
			&CR[reord_index[IM]*ND],1);

	    }

	    transpose(TrQF,CR,ND,NM);

	    AF->trans_mat_vec_mult_omp(TrQF,QF,NM,ND);
//	    c_zgemm('C','N',NM,NM,ND,CPX(1.0,0.0),CR,ND,QF,ND,CPX(0.0,0.0),AR,NM);
            for(IM=0;IM<NM;IM++){
                AR[IM*(1+NM)]=c_zdotc(ND,&CR[IM*ND],1,&QF[IM*ND],1);
            }

	    BF->trans_mat_vec_mult_omp(TrQF,QF,NM,ND);
//	    c_zgemm('C','N',NM,NM,ND,CPX(1.0,0.0),CR,ND,QF,ND,CPX(0.0,0.0),BR,NM);
            for(IM=0;IM<NM;IM++){
                BR[IM*(1+NM)]=c_zdotc(ND,&CR[IM*ND],1,&QF[IM*ND],1);
            }

            omp_set_num_threads(omp_get_max_threads()/loc_size);

	    for(IM=0;IM<NM;IM++){

		if((abs(lambdaR[IM])>=LBound)&&(abs(lambdaR[IM])<=UBound)){
		  
		    if(abs(AR[IM*(1+NM)]/BR[IM*(1+NM)]-lambdaR[IM])/abs(lambdaR[IM])*100<error_criterion){

		        if(abs(lambdaR[IM]*conj(lambdaR[IM])-CPX(1.0,0.0))<error_criterion){
			    index_to_keep[IM]   = 2;
			    prop_index[no_prop] = IM;
			    no_prop++;
			}else{
			    index_to_keep[IM] = 1;
			}
			no_ev++;
		    }
		}
	    }

	    /*
	    for(IM=0;IM<NM;IM++){

	        c_zgemm('N','N',ND,1,NM,CPX(1.0,0.0),QF,ND,&VR[IM*NM],NM,\
			CPX(0.0,0.0),&CR[reord_index[IM]*ND],ND);
	    
		c_zscal(ND,CPX(1.0/c_dznrm2(ND,&CR[reord_index[IM]*ND],1),0.0),\
			&CR[reord_index[IM]*ND],1);
	   
		cond1 = 0;
		cond2 = 0;
		cond3 = 0;
		cond4 = 0;

		if((abs(lambdaR[IM])>=1.0/feast_radius)&&(abs(lambdaR[IM])<=feast_radius)){
		    cond1 = 1;
		}

		for(IN=0;IN<NM;IN++){

		    if(abs(imag(lambdaR[IM])+imag(lambdaR[IN]))<error_criterion){
		        cond2     = 1;
		    }

		    if(abs(real(lambdaR[IM])-real(lambdaR[IN]))<error_criterion){
		        cond3     = 1;
		    }
	      
		    if(abs(lambdaR[IM]*lambdaR[IN]-CPX(1.0,0.0))<error_criterion){
		        cond4     = 1;
		    }
		}

		if(cond1&&cond2&&cond3&&cond4){
		    if(abs(lambdaR[IM]*conj(lambdaR[IM])-CPX(1.0,0.0))<error_criterion){
		        index_to_keep[IM]   = 2;
			prop_index[no_prop] = IM;
			no_prop++;
		    }else{
		        index_to_keep[IM] = 1;
		    }
		    no_ev++;
		}
	    }
	    */

	    for(IM=0;IM<no_prop;IM++){

	        kM = CPX(0.0,1.0)*log(lambdaR[prop_index[IM]]);

	        for(IN=0;IN<no_prop;IN++){

		    kN = CPX(0.0,1.0)*log(lambdaR[prop_index[IN]]);

		    cond1 = abs((real(kM)-real(kN))/real(kM))<0.02;
		    cond2 = abs((imag(kM)+imag(kN))/imag(kM))<0.02;
		    cond3 = abs(imag(kM))>1.0e-6;

		    if(cond1&&cond2&&cond3){
		        index_to_keep[prop_index[IM]] = 1;
			break;
		    }
		}
	    }

	    IN = 0;
	    for(IM=0;IM<NM;IM++){
	  
	        if((abs(lambdaR[IM])>=1.0/(mult_factor*feast_radius))&&\
		   (abs(lambdaR[IM])<=mult_factor*feast_radius)){
		    BF->mat_vec_mult(&CR[IM*ND],&YF[IN*ND],1);
		    IN++;
		}
	    }

	    no_search = NM;
	    NM        = IN;
	}

	if(!mpi_rank) printf("Sequential time: %e\n",get_time(t0));

	if(loc_size>1){
	    MPI_Bcast(&no_ev,1,MPI_INT,root_rank,comm);
	    MPI_Bcast(&NM,1,MPI_INT,root_rank,comm);
	    if(iteration<act_max_iter-1){
		MPI_Bcast(YF,NM*ND,MPI_DOUBLE_COMPLEX,root_rank,comm);
	    }
	}

	if((iteration==(act_max_iter-1))&&(no_ev<Round(NM_for_feast/4.0))){

	    act_max_iter = min(act_max_iter+max_iteration,2*max_iteration);
	    act_radius   = 1.5*act_radius;
	    NM           = NM_for_feast;

	    //GetIntegrationPoints(mult_factor*act_radius,1.0/(mult_factor*act_radius));
	    
	    if(sign>0){
	        GetIntegrationPoints(1.0+mult_factor/act_radius,1.0/(mult_factor*act_radius));
	    }else{
		GetIntegrationPoints(mult_factor*act_radius,1.0-mult_factor/act_radius);
	    }
	    
	    if(loc_rank==root_rank){
	        init_Y_matrix(YF,NM,ND);
	    }

	    if(loc_size>1){
	        MPI_Bcast(YF,NM*ND,MPI_DOUBLE_COMPLEX,root_rank,comm);
	    }
	}
	/*
	if(mpi_rank==1999){
	  cout<<no_ev<<" "<<no_search<<" "<<NM<<endl;
	  for(int i=0;i<no_search;i++){
	    cout<<"i: "<<i<<" "<<index_to_keep[i]<<" "<<lambdaR[i]<<" "<<abs(lambdaR[i])<<" "<<CPX(0.0,1.0)*log(lambdaR[i])<<endl;
	  }
	}
	*/
    }

    if(loc_rank==root_rank){
        sort_eig(kphase,Vsurf,dEk_dk,Ntr,ind_Ntr,Nref,ind_Nref,T01,ND,no_search,sign);
    }

    if(!mpi_rank || mpi_rank==2*loc_size-1){
      for(int i=0;i<*Nref;i++){
	cout<<"Result: "<<i<<" "<<kphase[ind_Nref[i]]<<" with sign "<<sign<<endl;
      }
    }

    delete[] reord_index;
    delete[] map_index;
    delete[] prop_index;
    delete[] QF;
    delete[] YF;
    delete[] TrQF;
    delete AF;
    delete BF;
 
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::calc_kphase_sym(TCSR<T> *D,TCSR<T> *P,TCSR<CPX> *T01,int ND,int NP,\
					CPX *kphasel,CPX *kphaser,CPX *Vsurfl,CPX *Vsurfr, \
					double *dEkl_dk, double *dEkr_dk, int *Ntrl, int *Ntrr,	\
					int *ind_Ntrl,int *ind_Ntrr,int *Nrefl,int *Nrefr, \
					int *ind_Nrefl,int *ind_Nrefr,int NPROW,int NPCOL, \
					MPI_Comm comm,int *info)
{
  /*
 
    
    sort_eig(kphasel,Vsurfl,dEkl_dk,Ntrl,ind_Ntrl,Nrefl,ind_Nrefl,phase,VUp,VDown,T01,ND,NP,1);
    sort_eig(kphaser,Vsurfr,dEkr_dk,Ntrr,ind_Ntrr,Nrefr,ind_Nrefr,phase,VUp,VDown,T01,ND,NP,-1);


  */
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::calc_kphase(TCSR<T> *D,T *P,TCSR<CPX> *T01,int ND,int NP,CPX *kphase,\
				    CPX *Vsurf,double *dEk_dk,int *Ntr,int *ind_Ntr,int *Nref, \
				    int *ind_Nref,int sign,int NPROW,int NPCOL,MPI_Comm comm,\
				    int *info)
{
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::calc_kphase_sym(TCSR<T> *D,T *P,TCSR<CPX> *T01,int ND,int NP,\
					CPX *kphasel,CPX *kphaser,CPX *Vsurfl,CPX *Vsurfr, \
					double *dEkl_dk, double *dEkr_dk, int *Ntrl, int *Ntrr,	\
					int *ind_Ntrl,int *ind_Ntrr,int *Nrefl,int *Nrefr, \
					int *ind_Nrefl,int *ind_Nrefr,int NPROW,int NPCOL, \
					MPI_Comm comm,int *info)
{
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::get_Q_matrix(CPX *Q,CPX *QLoc,CPX *Y,TCSR<CPX> *A,TCSR<CPX> *B,\
				     int *map_index,int NM,int N,MPI_Comm comm)
{

    int IR;

    init_var(QLoc,N*NM);

    if(contour_factor==2){
        c_dscal(N*NM,0.0,((double*)Y)+1,2);
    }

    for(IR=loc_rank;IR<2*NQDR;IR=IR+loc_size){

        if(IR<NQDR/contour_factor){

	    if(slab_per_bc==1){
	        solve_Q_sparse(QLoc,Y,A,B,map_index,NM,N,z_out[IR],dz_dtheta_out[IR]*dtheta[IR]);
	    }else{
	        solve_Q_full(QLoc,Y,A,B,map_index,NM,N,z_out[IR],dz_dtheta_out[IR]*dtheta[IR]);
	    }

	}else if(IR>=NQDR && IR<NQDR+NQDR/contour_factor){

	    if(slab_per_bc==1){
	        solve_Q_sparse(QLoc,Y,A,B,map_index,NM,N,z_in[IR-NQDR],-dz_dtheta_in[IR-NQDR]*\
			       dtheta[IR-NQDR]);
	    }else{
	        solve_Q_full(QLoc,Y,A,B,map_index,NM,N,z_in[IR-NQDR],-dz_dtheta_in[IR-NQDR]*\
			     dtheta[IR-NQDR]);
	    }

	}
    }

    if(loc_size>1){
        MPI_Reduce(QLoc,Q,N*NM,MPI_DOUBLE_COMPLEX,MPI_SUM,root_rank,comm);
    }else{
        c_zcopy(N*NM,QLoc,1,Q,1);
    }

    if(loc_rank==root_rank && contour_factor==2){
        c_dscal(N*NM,2.0, (double*)Q   ,2);
        c_dscal(N*NM,0.0,((double*)Q)+1,2);
    }
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::solve_Q_sparse(CPX *Q,CPX *Y,TCSR<CPX> *A,TCSR<CPX> *B,int *map_index,\
				       int NM,int N,CPX z,CPX dz_dtheta_dtheta)
{
    int IC;
    CPX *QCol = new CPX[N];
    LinearSolver<CPX> *solver;

    add_sparse(A,B,map_index,-z);

    solver = new Umfpack<CPX>(A,MPI_COMM_SELF);
 
    solver->prepare();
 
    for(IC=0;IC<NM;IC++){
        solver->solve_equation(QCol,&Y[N*IC],1);
	c_zaxpy(N,CPX(0.0,1.0/(2.0*PI))*dz_dtheta_dtheta,QCol,1,&Q[N*IC],1);
    }
    delete solver;

    add_sparse(A,B,map_index,z);

    delete[] QCol;
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::solve_Q_full(CPX *Q,CPX *Y,TCSR<CPX> *A,TCSR<CPX> *B,int *map_index,\
				     int NM,int N,CPX z,CPX dz_dtheta_dtheta)
{

    int iter;
    int info;
    int IC;
    int NR    = A->size/slab_per_bc;
    int *ipiv = new int[NR];
    CPX *M    = new CPX[NR*NR];
    CPX *rhs  = new CPX[NR*NM];
    CPX *Yact = new CPX[NR*NM];

    double *Mr;
    double *rhsr;
    double *Yr;
    if(contour_factor==2){
        Mr   = new double[NR*NR];
        rhsr = new double[NR*NM];
        Yr   = new double[NR*NM];
    }


    init_var(M,NR*NR);

    add_sparse_block_to_full(M,B,CPX(1.0,0.0),0);

    for(IC=0;IC<NM;IC++){
        c_zcopy(NR,&Y[IC*A->size],1,&rhs[IC*NR],1);
    }

    for(IC=0;IC<NM;IC++){
        c_zcopy(NR,&Y[IC*A->size],1,&rhs[IC*NR],1);
    }

    for(iter=1;iter<slab_per_bc;iter++){

        c_zscal(NR*NR,z,M,1);
	add_sparse_block_to_full(M,A,CPX(1.0,0.0),iter);
	
	for(IC=0;IC<NM;IC++){
	    c_zcopy(NR,&Y[iter*NR+IC*A->size],1,&Yact[IC*NR],1);
	}
        if(contour_factor==2){
            c_dcopy(NR*NM,(double*)Yact,2,Yr,1);
            c_dcopy(NR*NR,(double*)M,2,Mr,1);
            c_dgemm('N','N',NR,NM,NR,-1.0,Mr,NR,Yr,NR,0.0,rhsr,NR);
            c_daxpy(NR*NM,1.0,rhsr,1,(double*)rhs,2);
            c_dcopy(NR*NR,((double*)M)+1,2,Mr,1);
            c_dgemm('N','N',NR,NM,NR,-1.0,Mr,NR,Yr,NR,0.0,rhsr,NR);
            c_daxpy(NR*NM,1.0,rhsr,1,((double*)rhs)+1,2);
        }else{
            c_zgemm('N','N',NR,NM,NR,CPX(-1.0,0.0),M,NR,Yact,NR,CPX(1.0,0.0),rhs,NR);
        }

    }

    c_zscal(NR*NR,z,M,1);
    add_sparse_block_to_full(M,A,CPX(1.0,0.0),slab_per_bc);
    
    c_zgetrf(NR,NR,M,NR,ipiv,&info);
    c_zgetrs('N',NR,NM,M,NR,ipiv,rhs,NR,&info);

    for(IC=0;IC<NM;IC++){
	c_zaxpy(NR,CPX(0.0,1.0/(2.0*PI))*dz_dtheta_dtheta,&rhs[IC*NR],1,\
		    &Q[(slab_per_bc-1)*NR+IC*A->size],1);
    }

    for(iter=slab_per_bc-2;iter>=0;iter--){

        c_zscal(NR*NM,z,rhs,1);

	for(IC=0;IC<NM;IC++){
	    c_zaxpy(NR,CPX(1.0,0.0),&Y[(iter+1)*NR+IC*A->size],1,&rhs[IC*NR],1);
	    c_zaxpy(NR,CPX(0.0,1.0/(2.0*PI))*dz_dtheta_dtheta,&rhs[IC*NR],1,\
		    &Q[iter*NR+IC*A->size],1);
	}
    }

    delete[] Yact;
    delete[] rhs;
    delete[] M;
    delete[] ipiv;

    if(contour_factor==2){
        delete[] Mr;
        delete[] rhsr;
        delete[] Yr;
    }
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::init_Y_matrix(CPX *Y,int NM,int N)
{
    //all the CPUs should see the same Y
    int IR,IC;

    //srand(time(0));
    srand(0);

    for(IC=0;IC<NM;IC++){
        for(IR=0;IR<N;IR++){
	    Y[IR+IC*N] = CPX(randn(),randn());
	}
	c_zscal(N,CPX(1.0/c_dznrm2(N,&Y[IC*N],1),0.0),&Y[IC*N],1);
    }

}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::init_matrix_single(TCSR<CPX> *A,TCSR<CPX> *B,int *map_index,TCSR<T> *D,\
					   TCSR<T> *P,int ND,int NP)
{

    int i,j;
    
    A->copy_contain(D,1.0);

    B->n_nonzeros = 0;

    for(i=0;i<P->size;i++){

        B->index_i[i] = 0;

        for(j=P->edge_i[i]-P->findx;j<P->edge_i[i+1]-P->findx;j++){

	    B->index_j[B->n_nonzeros] = P->index_j[j]-P->findx+ND-NP+B->findx;
	    B->nnz[B->n_nonzeros]     = convert<CPX,T>(P->nnz[j]);

	    B->index_i[i]++;
	    B->n_nonzeros++;
	}
    }
    B->get_row_edge();

    merge_A_and_B(A,B,CPX(1.0,0.0));
    map_A_and_B(A,B,map_index);

}
/************************************************************************************************/

template <class T>
void InjectionFeast<T>::init_matrix_double(TCSR<CPX> *A,TCSR<CPX> *B,int *map_index,TCSR<T> *D,\
					   TCSR<T> *P,int ND,int NP)
{

    int i,j;
    int row_ind;
    int col_ind;
    CPX *line = new CPX[D->size];

    A->n_nonzeros = 0;
    B->n_nonzeros = 0;
  
    for(i=0;i<ND;i++){ 

        A->index_i[i] = 0;
        B->index_i[i] = 0;

	if(i<ND/slab_per_bc){

	    init_var(line,D->size);

	    row_ind = i+ND/2-ND/slab_per_bc;

	    for(j=D->edge_i[row_ind]-D->findx;j<D->edge_i[row_ind+1]-D->findx;j++){

	        col_ind           = D->index_j[j]-D->findx;

		if((col_ind<ND/2)||(col_ind>=ND-ND/slab_per_bc)){
		    line[col_ind] = convert<CPX,T>(D->nnz[j]);
		}
	    }

	    for(j=P->edge_i[row_ind]-P->findx;j<P->edge_i[row_ind+1]-P->findx;j++){

	        col_ind                   = P->index_j[j]-P->findx+ND-NP;

		if((col_ind<ND/2)||(col_ind>=ND-ND/slab_per_bc)){
		    line[col_ind]         = line[col_ind]+convert<CPX,T>(P->nnz[j]);
		}
		  
		B->index_j[B->n_nonzeros] = col_ind-ND+ND/slab_per_bc+B->findx;
		B->nnz[B->n_nonzeros]     = convert<CPX,T>(P->nnz[j]);

		B->index_i[i]++;
		B->n_nonzeros++;
	    }

	    row_ind = i+ND/2-2*ND/slab_per_bc;

	    for(j=D->edge_i[row_ind]-D->findx;j<D->edge_i[row_ind+1]-D->findx;j++){

	        col_ind           = D->index_j[j]-D->findx+ND/slab_per_bc;

		if((col_ind>=ND/2)&&(col_ind<ND-ND/slab_per_bc)){
		    line[col_ind] = convert<CPX,T>(D->nnz[j]);
		}
	    }

	    for(j=P->edge_i[row_ind]-P->findx;j<P->edge_i[row_ind+1]-P->findx;j++){

	        col_ind           = P->index_j[j]-P->findx+ND-NP+ND/slab_per_bc;

		if((col_ind>=ND/2)&&(col_ind<ND-ND/slab_per_bc)){
		    line[col_ind] = line[col_ind]+convert<CPX,T>(P->nnz[j]);
		}
	    }

	    for(j=0;j<D->size;j++){

	        if(abs(line[j])){

		    A->nnz[A->n_nonzeros]     = line[j];
		    A->index_j[A->n_nonzeros] = j+A->findx;

		    A->index_i[i]++;
		    A->n_nonzeros++;
		}
	    }

	}else{

	    A->index_j[A->n_nonzeros] = i-ND/slab_per_bc+A->findx;
	    A->nnz[A->n_nonzeros]     = CPX(1.0,0.0);

	    A->index_i[i]++;
	    A->n_nonzeros++;

	    B->index_j[B->n_nonzeros] = i+B->findx;
	    B->nnz[B->n_nonzeros]     = CPX(1.0,0.0);

	    B->index_i[i]++;
	    B->n_nonzeros++;
	}
    }
    A->get_row_edge();
    B->get_row_edge();

    //merge_A_and_B(A,B,CPX(0.0,0.0));
    //map_A_and_B(A,B,map_index);

    delete[] line;
}
/************************************************************************************************/

template <class T>
int InjectionFeast<T>::get_type(TCSR<CPX> *T01,CPX* Vin,CPX kin,double *dEk_dk,int ind_to_keep,\
				int sign)
{
    int type;
    //double dEk;

    *dEk_dk = 0.0;
    
    if(((imag(kin)<imag_limit)&&(imag(kin)>-imag_limit))&&(ind_to_keep==2)){
        
        *dEk_dk = sign*calc_derivative(T01,Vin,kin);
 
        if(*dEk_dk>0){
	  type = 1;
	}else{
	  type = -1;
	}
        
    }else{

        if(imag(kin)>imag_limit){
	    type = sign;
	}else{
	    type = -sign;
	}
	/*
        if(ind_to_keep==2){

	    dEk = sign*calc_derivative(T01,Vin,kin);

	    if(dEk>0){
	        type = 1;
	    }else{
	        type = -1;
	    }

	}else{

	    if(imag(kin)>imag_limit){
	        type = sign;
	    }else{
	        type = -sign;
	    }
	}
	*/
    }

    return type;
    
}

/************************************************************************************************/

template <class T>
double InjectionFeast<T>::calc_derivative(TCSR<CPX> *T01,CPX* Vin,CPX kin)
{
    double dEk_dk;
    CPX *Vout = new CPX[T01->size];
    
    T01->mat_vec_mult(Vin,Vout,1);
    //no minus sign in dEk_dk because T01=-H01 used instead of H01
    dEk_dk = 2.0*imag(c_zdotc(T01->size,Vin,1,Vout,1)*exp(CPX(0.0,1.0)*kin)); 

    delete[] Vout;
    
    return dEk_dk;
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::GaussLegendre(double *alpha,double *dalpha,int NR)
{

    int IR;
    int *sort_index = new int[NR];
    double *CM      = new double[NR*NR];

    init_var(CM,NR*NR);

    for(IR=0;IR<NR-1;IR++){
	CM[IR+(IR+1)*NR] = (IR+1.0)/sqrt(4.0*(IR+1.0)*(IR+1.0)-1.0);
    }

    eig(CM,alpha,NR);

    sort_vec(alpha,sort_index,NR);

    for(IR=0;IR<NR;IR++){
        dalpha[sort_index[IR]] = 2.0*pow(CM[IR*NR],2.0);
    }

    delete[] sort_index;
    delete[] CM;
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::eig(double *A,double *lambda,int N)
{
    int info;
    char jobz        = 'V';
    char uplo        = 'U';
    int lwork        = 20*N;
    double *loc_work = new double[lwork];

    c_dsyev(jobz,uplo,N,A,N,lambda,loc_work,lwork,&info);

    delete[] loc_work;
}

/************************************************************************************************/

template <class T>
int InjectionFeast<T>::gen_eig(CPX *a,CPX *b,int N,CPX *vr,CPX *lambda)
{
    
    char jobvl   = 'N';
    char jobvr   = 'V';
    int n        = N;
    int lda      = n;
    int ldb      = n;
    CPX *vl      = NULL;
    int ldvl     = 1;
    int ldvr     = n;
    int lwork    = factor_lwork*n;
    int info,IN;

    c_zggev(jobvl,jobvr,n,a,lda,b,ldb,lambda_up,lambda_down,vl,ldvl,vr,ldvr,work,lwork,\
	    rwork,&info);

    for(IN=0;IN<n;IN++){
        if(abs(lambda_down[IN])>tollim){
            lambda[IN] = lambda_up[IN]/lambda_down[IN];
        }else{
            lambda[IN] = CPX(INF,INF);
        }
    }

    return info;
}

/************************************************************************************************/

template <class T>
int InjectionFeast<T>::gen_eig_real(CPX *a,CPX *b,int N,CPX *vr,CPX *lambda)
{
    
    char jobvl   = 'N';
    char jobvr   = 'V';
    int n        = N;
    int lda      = n;
    int ldb      = n;
    CPX *vl      = NULL;
    int ldvl     = 1;
    int ldvr     = n;
    int lwork    = factor_lwork*n;
    int info,IN;

    double* lambda_up_real = new double[N];
    double* lambda_up_imag = new double[N];
    double* lambda_down_real = new double[N];
    double* vl_real = new double[1];
    double* vr_real = new double[N*N];

    double* areal = new double[N*N];
    double* breal = new double[N*N];
    c_dcopy(N*N,(double*)a,2,areal,1);
    c_dcopy(N*N,(double*)b,2,breal,1);

    double work_test;
    c_dggev(jobvl,jobvr,n,areal,lda,breal,ldb,lambda_up_real,lambda_up_imag,lambda_down_real,\
            vl_real,ldvl,vr_real,ldvr,&work_test,-1,&info);
    int lwork_real=int(work_test);
    double* work_real = new double[lwork_real];

    c_dggev(jobvl,jobvr,n,areal,lda,breal,ldb,lambda_up_real,lambda_up_imag,lambda_down_real,\
            vl_real,ldvl,vr_real,ldvr,work_real,lwork_real,&info);

    delete[] work_real;
    delete[] areal;
    delete[] breal;

    for(IN=0;IN<n;IN++){
        if(abs(lambda_down_real[IN])>tollim){
            lambda[IN] = CPX(lambda_up_real[IN],lambda_up_imag[IN])/lambda_down_real[IN];
        }else{
            lambda[IN] = CPX(INF,INF);
        }
    }

    c_zscal(N*N,CPX(0.0,0.0),vr,1);
    int iwhile=0;
    while (iwhile<N) {
        if (!lambda_up_imag[iwhile]) {
            c_dcopy(N,&vr_real[iwhile*N],1,(double*)&vr[iwhile*N],2);
            iwhile++;
        }
        else {
            c_dcopy(N,&vr_real[ iwhile   *N],1, (double*)&vr[ iwhile   *N]   ,2);
            c_dcopy(N,&vr_real[(iwhile+1)*N],1,((double*)&vr[ iwhile   *N])+1,2);
            c_dcopy(N,&vr_real[ iwhile   *N],1,((double*)&vr[(iwhile+1)*N])+1,2);
            c_dcopy(N,&vr_real[(iwhile+1)*N],1, (double*)&vr[(iwhile+1)*N]   ,2);
            iwhile+=2;
        }
    } // END WHILE

    delete[] lambda_up_real;
    delete[] lambda_up_imag;
    delete[] lambda_down_real;
    delete[] vl_real;
    delete[] vr_real;

    return info;
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::add_sparse(TCSR<CPX> *A,TCSR<CPX> *B,int *map_index,CPX factor)
{

    int i;

    for(i=0;i<B->n_nonzeros;i++){
        A->nnz[map_index[i]] = A->nnz[map_index[i]]+factor*B->nnz[i];
    }

}
/************************************************************************************************/

template <class T>
void InjectionFeast<T>::transpose(CPX *A,CPX *B,int nrow,int ncol)
{
    int i;

    for(i=0;i<ncol;i++) c_zcopy(nrow,&B[i*nrow],1,&A[i],ncol);
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::merge_A_and_B(TCSR<CPX> *A,TCSR<CPX> *B,CPX factor)
{

    int IR,IC;
    int n_nonzeros;
    TCSR<CPX> *sparse_patt = new TCSR<CPX>(A->size,A->n_nonzeros+B->n_nonzeros,A->findx);
    CPX *line              = new CPX[A->size];

    n_nonzeros             = 0;

    for(IR=0;IR<A->size;IR++){

        init_var(line,A->size);

        sparse_patt->index_i[IR] = 0;        

	if(A->n_nonzeros>0){

	    for(IC=A->edge_i[IR]-A->findx;IC<A->edge_i[IR+1]-A->findx;IC++){

		line[A->index_j[IC]-A->findx] = A->nnz[IC];
	    }	
	}

	if(B->n_nonzeros>0){

	    for(IC=B->edge_i[IR]-B->findx;IC<B->edge_i[IR+1]-B->findx;IC++){

		  line[B->index_j[IC]-B->findx] = line[B->index_j[IC]-B->findx]+\
		    factor*B->nnz[IC]+CPX(2*tollim,0.0);
	    }	
	}

	for(IC=0;IC<A->size;IC++){

	    if(abs(line[IC])){
	        sparse_patt->nnz[n_nonzeros]     = line[IC];
		sparse_patt->index_j[n_nonzeros] = IC+A->findx;
		sparse_patt->index_i[IR]++;
		n_nonzeros++;
	    }
	}
    }

    sparse_patt->n_nonzeros = n_nonzeros;
    sparse_patt->first_row  = A->first_row;
    sparse_patt->size_tot   = A->size_tot;
    sparse_patt->get_row_edge();
    sparse_patt->get_diag_pos();

    delete A;

    A = new TCSR<CPX>(sparse_patt);

    delete sparse_patt;
    delete[] line;
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::map_A_and_B(TCSR<CPX> *A,TCSR<CPX> *B,int *index)
{

    int IR,IC,icc;
    int ICmin;
    int index_ok;

    for(IR=0;IR<A->size;IR++){

        if(B->n_nonzeros>0){

	    if(B->index_i[IR]>0){

	        ICmin = A->edge_i[IR]-A->findx;

		for(icc=B->edge_i[IR]-B->findx;icc<B->edge_i[IR+1]-B->findx;icc++){

		    index_ok = 0;

		    for(IC=ICmin;IC<A->edge_i[IR+1]-A->findx;IC++){

		        if(A->index_j[IC]==B->index_j[icc]){
			    index[icc] = IC;
			    ICmin      = IC+1;
			    index_ok   = 1;
			    break;
			}
		    }

		    if(!index_ok){
		        printf("There was a problem with the A and B index mapping\n");
			abort();
		    }
		}
	    }
	}
    }
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::add_sparse_block_to_full(CPX *M,TCSR<CPX> *A,CPX factor,int iter)
{

    int IR,IC;
    int col_ind;
    int NR = A->size/slab_per_bc;

    if(!iter){

        for(IR=0;IR<NR;IR++){

	    for(IC=A->edge_i[IR]-A->findx;IC<A->edge_i[IR+1]-A->findx;IC++){

	        col_ind          = A->index_j[IC]-A->findx;
	        M[IR+col_ind*NR] = M[IR+col_ind*NR]-factor*A->nnz[IC];
	    }
	}

    }else{

        for(IR=0;IR<NR;IR++){

	    for(IC=A->edge_i[IR]-A->findx;IC<A->edge_i[IR+1]-A->findx;IC++){

	        col_ind = A->index_j[IC]-A->findx;

		if((col_ind>=(iter-1)*NR)&&(col_ind<iter*NR)){

		    M[IR+(col_ind-(iter-1)*NR)*NR] = M[IR+(col_ind-(iter-1)*NR)*NR]+\
		      factor*A->nnz[IC];
		}
	    }
	}
    }
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::GetIntegrationPoints(double R_out,double R_in)
{

    int IR;

    for(IR=0;IR<NQDR;IR++){
        dtheta[IR]        = 2.0*M_PI/NQDR;
        z_out[IR]         = R_out*exp(CPX(0.0,2.0*(IR+0.5)*M_PI/NQDR));
        dz_dtheta_out[IR] = CPX(0.0,1.0)*z_out[IR];
        z_in[IR]          = R_in*exp(CPX(0.0,2.0*(IR+0.5)*M_PI/NQDR));
        dz_dtheta_in[IR]  = CPX(0.0,1.0)*z_in[IR];
    }
}

/************************************************************************************************/

#endif
