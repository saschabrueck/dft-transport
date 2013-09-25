#ifndef __INJECTIONFEAST
#define __INJECTIONFEAST

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
    
    InjectionFeast();

    virtual ~InjectionFeast();

    virtual void initialize(int,int,int);
    virtual void calc_kphase(TCSR<T>*,T*,TCSR<CPX>*,int,int,CPX*,CPX*,double*,int*,\
			     int*,int*,int*,int,int,int,MPI_Comm,int*);
    virtual void calc_kphase_sym(TCSR<T>*,T*,TCSR<CPX>*,int,int,CPX*,CPX*,CPX*,CPX*,\
				 double*,double*,int*,int*,int*,int*,int*,int*,int*,\
				 int*,int,int,MPI_Comm,int*);
    
private:

    double error_criterion,mult_factor,radius;
    double imag_limit;
    double *alphaF,*dalphaF;
    double *theta,*dtheta;
    double *rwork;
    int mpi_size,mpi_rank;
    int factor_lwork;
    int NM_for_feast;
    int max_iteration;
    int NQDR;
    int *index_to_keep;
    CPX *AR,*BR,*CR;
    CPX *z_out,*dz_dtheta_out;
    CPX *z_in,*dz_dtheta_in;
    CPX *lambda_up,*lambda_down;
    CPX *lambdaR,*VR;
    CPX *work;

    void sort_eig(CPX*,CPX*,double*,int*,int*,int*,int*,TCSR<CPX>*,int,int,int);
    int get_type(TCSR<CPX>*,CPX*,CPX,double*,int);
    double calc_derivative(TCSR<CPX>*,CPX*,CPX);
    void init_matrix(TCSR<CPX>*,TCSR<CPX>*,int*,TCSR<T>*,T*,int,int);
    void GaussLegendre(double*,double*,int);
    void eig(double*,double*,int);
    int gen_eig(CPX*,CPX*,int,CPX*,CPX*);
    void get_Q_matrix(CPX*,CPX*,TCSR<CPX>*,TCSR<CPX>*,int*,int,int,MPI_Comm);
    void init_Y_matrix(CPX*,int,int);
    void add_sparse(TCSR<CPX>*,TCSR<CPX>*,int*,CPX factor);
    void transpose(CPX*,CPX*,int,int);
};

template<class T>
InjectionFeast<T>::InjectionFeast()
{
    error_criterion = 1.0e-3;
    imag_limit      = 1.0e-5;
    mult_factor     = 1.25;
    radius          = 100.0;
    max_iteration   = 2;
    NQDR            = 16;
    factor_lwork    = 20;

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

    int IR;

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

    GaussLegendre(alphaF,dalphaF,NQDR);

    for(IR=0;IR<NQDR;IR++){
        theta[IR]         = PI*(alphaF[IR]+1.0);
	dtheta[IR]        = PI*dalphaF[IR];
	z_out[IR]         = mult_factor*radius*exp(CPX(0.0,theta[IR]));
	dz_dtheta_out[IR] = mult_factor*radius*CPX(0.0,1.0)*exp(CPX(0.0,theta[IR])); 
	z_in[IR]          = 1.0/(mult_factor*radius)*exp(CPX(0.0,theta[IR]));
	dz_dtheta_in[IR]  = 1.0/(mult_factor*radius)*CPX(0.0,1.0)*exp(CPX(0.0,theta[IR])); 
    }

    if(!mpi_rank){
      cout<<"Feast parameter: "<<NM_for_feast<<endl;
    }
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::sort_eig(CPX *kphase,CPX *Vsurf,double *dEk_dk,int *Ntr,int *ind_Ntr,\
				 int *Nref,int *ind_Nref,TCSR<CPX> *T01,int ND,int NM,int sign)
{

    int IP,IM;
    int type;

    *Ntr=0;
    *Nref=0;    
    
    IP = 0;
    for(IM=0;IM<NM;IM++){

        if(index_to_keep[IM]){

	    c_zcopy(ND,&CR[IM*ND],1,&Vsurf[IP*ND],1);

	    kphase[IP] = CPX(0.0,1.0)*log(lambdaR[IM]);
	    
	    type       = get_type(T01,&Vsurf[IP*ND],kphase[IP],&dEk_dk[(*Ntr)],sign);

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
void InjectionFeast<T>::calc_kphase(TCSR<T> *D,T *P,TCSR<CPX> *T01,int ND,int NP,CPX *kphase, \
				    CPX *Vsurf,double *dEk_dk,int *Ntr,int *ind_Ntr,int *Nref, \
				    int *ind_Nref,int sign,int NPROW,int NPCOL,MPI_Comm comm,\
				    int *info)
{
    int iteration;
    int IM,IN;
    int no_ev;
    int cond1,cond2,cond3,cond4;
    int NM           = NM_for_feast;
    int *reord_index = new int[NM];
    int *map_index   = new int[D->n_nonzeros];
    CPX *QF          = new CPX[ND*NM];
    CPX *YF          = new CPX[ND*NM];
    CPX *TrQF        = new CPX[NM*ND];
    TCSR<CPX> *AF    = new TCSR<CPX>(D->size,D->n_nonzeros,D->findx);
    TCSR<CPX> *BF    = new TCSR<CPX>(D->size,D->n_nonzeros,D->findx);
    
    init_matrix(AF,BF,map_index,D,P,ND,NP);
    init_Y_matrix(YF,NM,ND);

    for(iteration=0;iteration<max_iteration;iteration++){

        get_Q_matrix(QF,YF,AF,BF,map_index,NM,ND,comm);

	transpose(TrQF,QF,ND,NM);

	AF->trans_mat_vec_mult(TrQF,CR,NM,ND);
	c_zgemm('C','N',NM,NM,ND,CPX(1.0,0.0),QF,ND,CR,ND,CPX(0.0,0.0),AR,NM);

	BF->trans_mat_vec_mult(TrQF,CR,NM,ND);
	c_zgemm('C','N',NM,NM,ND,CPX(1.0,0.0),QF,ND,CR,ND,CPX(0.0,0.0),BR,NM);

	gen_eig(AR,BR,NM,VR,lambdaR);

	sort_abs_vec(lambdaR,reord_index,NM);

	init_var(index_to_keep,NM);
	no_ev = 0;

	for(IM=0;IM<NM;IM++){

	    c_zgemm('N','N',ND,1,NM,CPX(1.0,0.0),QF,ND,&VR[IM*NM],NM,\
		    CPX(0.0,0.0),&CR[reord_index[IM]*ND],ND);
	    
	    c_zscal(ND,CPX(1.0/c_dznrm2(ND,&CR[reord_index[IM]*ND],1),0.0),\
		    &CR[reord_index[IM]*ND],1);
	   
	    cond1 = 0;
	    cond2 = 0;
	    cond3 = 0;
	    cond4 = 0;

	    if((abs(lambdaR[IM])>=1.0/radius)&&(abs(lambdaR[IM])<=radius)){
	        cond1 = 1;
	    }

	    for(IN=0;IN<NM;IN++){

	        if(abs(imag(lambdaR[IM])+imag(lambdaR[IN]))<error_criterion){
		    cond2 = 1;
		}

		if(abs(real(lambdaR[IM])-real(lambdaR[IN]))<error_criterion){
		    cond3 = 1;
		}
	      
		if(abs(lambdaR[IM]*lambdaR[IN]-CPX(1.0,0.0))<error_criterion){
		    cond4 = 1;
		}
	    }

	    if(cond1&&cond2&&cond3&&cond4){
	        index_to_keep[IM] = 1;
		no_ev++;
	    }
	}

	IN = 0;
	for(IM=0;IM<NM;IM++){
	  
	    if((abs(lambdaR[IM])>=1.0/(mult_factor*radius))&&\
		(abs(lambdaR[IM])<=mult_factor*radius)){
	       BF->mat_vec_mult(&CR[IM*ND],&YF[IN*ND],1);
	       IN++;
	    }
	}

	NM = IN;
    }

    sort_eig(kphase,Vsurf,dEk_dk,Ntr,ind_Ntr,Nref,ind_Nref,T01,ND,NM,sign);
    
    /*
    if(!mpi_rank){
        cout<<"number of propagating states: "<<*Ntr<<" "<<*Nref<<" "<<no_ev<<endl;
	for(int i=0;i<*Ntr;i++){
	    cout<<"i: "<<i<<" "<<kphase[ind_Ntr[i]]<<" "<<dEk_dk[i]<<endl;
	}
    }
    */
    
    delete[] reord_index;
    delete[] map_index;
    delete[] QF;
    delete[] YF;
    delete[] TrQF;
    delete AF;
    delete BF;
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
  /*
 
    
    sort_eig(kphasel,Vsurfl,dEkl_dk,Ntrl,ind_Ntrl,Nrefl,ind_Nrefl,phase,VUp,VDown,T01,ND,NP,1);
    sort_eig(kphaser,Vsurfr,dEkr_dk,Ntrr,ind_Ntrr,Nrefr,ind_Nrefr,phase,VUp,VDown,T01,ND,NP,-1);


  */
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::get_Q_matrix(CPX *Q,CPX *Y,TCSR<CPX> *A,TCSR<CPX> *B,int *map_index,\
				     int NM,int N,MPI_Comm comm)
{

    int IR,IC;
    CPX *QCol = new CPX[N];
    LinearSolver<CPX> *solver;

    init_var(Q,N*NM);

    for(IR=0;IR<NQDR;IR++){

	add_sparse(A,B,map_index,-z_out[IR]);

	solver = new Umfpack<CPX>(A,comm);
 
	solver->prepare();
 
	for(IC=0;IC<NM;IC++){
	    solver->solve_equation(QCol,&Y[N*IC],1);
	    c_zaxpy(N,CPX(0.0,1.0/(2.0*PI))*dz_dtheta_out[IR]*dtheta[IR],QCol,1,&Q[N*IC],1);
	}
	delete solver;

	add_sparse(A,B,map_index,z_out[IR]-z_in[IR]);

	solver = new Umfpack<CPX>(A,comm);

	solver->prepare();
	for(IC=0;IC<NM;IC++){
	    solver->solve_equation(QCol,&Y[N*IC],1);
	    c_zaxpy(N,CPX(0.0,-1.0/(2.0*PI))*dz_dtheta_in[IR]*dtheta[IR],QCol,1,&Q[N*IC],1);
	}
	delete solver;

	add_sparse(A,B,map_index,z_in[IR]);
    }

    delete[] QCol;
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
void InjectionFeast<T>::init_matrix(TCSR<CPX> *A,TCSR<CPX> *B,int *map_index,TCSR<T> *D,T *P,\
				    int ND,int NP)
{

    int i,j,IC;
    
    A->copy_contain(D,1.0);

    B->n_nonzeros = 0;

    for(i=0;i<ND;i++){ 

        B->index_i[i] = 0;

	for(j=0;j<NP;j++){

	    if(abs(P[i+j*ND])>tollim){

	        B->index_j[B->n_nonzeros] = j+ND-NP+B->findx;
                B->nnz[B->n_nonzeros]     = convert<CPX,T>(P[i+j*ND]);

		for(IC=A->edge_i[i]-A->findx;IC<A->edge_i[i+1]-A->findx;IC++){

		    if(A->index_j[IC]==B->index_j[B->n_nonzeros]){

		        A->nnz[IC]               = A->nnz[IC]+B->nnz[B->n_nonzeros];
			map_index[B->n_nonzeros] = IC;
			break;
		    }
		}

                B->index_i[i]++;
                B->n_nonzeros++;
	    }
	}
    }
    B->get_row_edge();

}
/************************************************************************************************/

template <class T>
int InjectionFeast<T>::get_type(TCSR<CPX> *T01,CPX* Vin,CPX kin,double *dEk_dk,int sign)
{
    int type;
    *dEk_dk = 0.0;
    
    if((imag(kin)<imag_limit)&&(imag(kin)>-imag_limit)){
        
        *dEk_dk = sign*calc_derivative(T01,Vin,kin);
 
        if(*dEk_dk>0){type=1;}
        else{type=-1;}
        
    }else{
        
        if(imag(kin)>imag_limit){
            type = sign;
        }else{
            type = -sign;
        }
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

#endif
