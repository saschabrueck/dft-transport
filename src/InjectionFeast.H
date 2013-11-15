#ifndef __INJECTIONFEAST
#define __INJECTIONFEAST

#include <mpi.h>
#include <omp.h>
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
#include "LinearSolver.H"
#include "Umfpack.H"

template<class T>
class InjectionFeast{
public:
    
    InjectionFeast();

    virtual ~InjectionFeast();

    virtual void initialize(int,int);
    virtual void calc_kphase(TCSR<CPX>*,TCSR<CPX>*,int,CPX*,CPX*,int*,MPI_Comm,int*);
    
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
    void get_Q_matrix(CPX*,CPX*,TCSR<CPX>*,TCSR<CPX>*,int,int,MPI_Comm);
    void get_Q_matrix_mpi(CPX*,CPX*,TCSR<CPX>*,TCSR<CPX>*,int,int,MPI_Comm);
    void get_Q_matrix_omp(CPX*,CPX*,TCSR<CPX>*,TCSR<CPX>*,int,int,MPI_Comm);
    void get_Q_matrix_omp_solver(CPX*,CPX*,LinearSolver<CPX>**,LinearSolver<CPX>**,int,int);
    void get_Q_matrix_omp_prepare(TCSR<CPX>*,TCSR<CPX>*,LinearSolver<CPX>**,\
                                  LinearSolver<CPX>**,int,int,MPI_Comm);
    void get_Q_matrix_omp_delete(LinearSolver<CPX>**,LinearSolver<CPX>**);
    void init_Y_matrix(CPX*,int,int);
    void add_sparse(TCSR<CPX>*,TCSR<CPX>*,int*,CPX factor);
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
void InjectionFeast<T>::initialize(int ND,int NM)
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

}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::calc_kphase(TCSR<CPX> *AF,TCSR<CPX> *BF,int ND,CPX *kphase, \
				    CPX *Vsurf,int *IPout,MPI_Comm comm,int *info)
{
    int iteration;
    int IM,IN;
    int cond1,cond2,cond3,cond4;
    int NM           = NM_for_feast;
    int *reord_index = new int[NM];
    CPX *QF          = NULL;
    CPX *TrQF        = NULL;
    CPX *YF          = NULL;

    int iam,nprocs;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&iam);

    if (!iam){
        QF           = new CPX[ND*NM];
        TrQF         = new CPX[NM*ND];
        YF           = new CPX[ND*NM];
        init_Y_matrix(YF,NM,ND);
    }

//    LinearSolver<CPX> **solver_out = new LinearSolver<CPX>*[NQDR];
//    LinearSolver<CPX> **solver_in = new LinearSolver<CPX>*[NQDR];
//    get_Q_matrix_omp_prepare(AF,BF,solver_out,solver_in,NM,ND,comm);

    for(iteration=0;iteration<max_iteration;iteration++){

//        get_Q_matrix(QF,YF,AF,BF,NM,ND,comm);

//        get_Q_matrix_omp(QF,YF,AF,BF,NM,ND,comm);

        get_Q_matrix_mpi(QF,YF,AF,BF,NM,ND,comm);

//        get_Q_matrix_omp_solver(QF,YF,solver_out,solver_in,NM,ND);

        if (!iam){
            full_transpose(NM,ND,QF,TrQF);
            AF->trans_mat_vec_mult(TrQF,CR,NM,ND);
            c_zgemm('C','N',NM,NM,ND,CPX(1.0,0.0),QF,ND,CR,ND,CPX(0.0,0.0),AR,NM);
            BF->trans_mat_vec_mult(TrQF,CR,NM,ND);
            c_zgemm('C','N',NM,NM,ND,CPX(1.0,0.0),QF,ND,CR,ND,CPX(0.0,0.0),BR,NM);
        }
        MPI_Bcast(AR,NM*NM,MPI_DOUBLE_COMPLEX,0,comm);
        MPI_Bcast(BR,NM*NM,MPI_DOUBLE_COMPLEX,0,comm);

	gen_eig(AR,BR,NM,VR,lambdaR);

	sort_abs_vec(lambdaR,reord_index,NM);

	init_var(index_to_keep,NM);

	if (!iam) for(IM=0;IM<NM;IM++){

	    c_zgemm('N','N',ND,1,NM,CPX(1.0,0.0),QF,ND,&VR[IM*NM],NM,\
		    CPX(0.0,0.0),&CR[reord_index[IM]*ND],ND);
	    
	    c_zscal(ND,CPX(1.0/c_dznrm2(ND,&CR[reord_index[IM]*ND],1),0.0),\
		    &CR[reord_index[IM]*ND],1);

        }
        MPI_Bcast(CR,NM*ND,MPI_DOUBLE_COMPLEX,0,comm);


	for(IM=0;IM<NM;IM++){
	   
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
	    }
	}

        if (!iam){
            IN = 0;
            for(IM=0;IM<NM;IM++){
                double abslambdaR=abs(lambdaR[IM]);
                if((abslambdaR>=1.0/(mult_factor*radius))&&(abslambdaR<=mult_factor*radius)){
                    BF->mat_vec_mult(&CR[IM*ND],&YF[IN*ND],1);
                    IN++;
                }
            }
            NM = IN;
        }
        MPI_Bcast(&NM,1,MPI_INT,0,comm);
    }

//    get_Q_matrix_omp_delete(solver_out,solver_in);
//    delete[] solver_out;
//    delete[] solver_in;

    int IP = 0;
    for(IM=0;IM<NM;IM++){

        if(index_to_keep[IM]){

            c_zcopy(ND,&CR[IM*ND],1,&Vsurf[IP*ND],1);

            kphase[IP] = CPX(0.0,1.0)*log(lambdaR[IM]);

            IP++;
        }
    }

    IPout[0]=IP;

    delete[] reord_index;

    if (!iam){
        delete[] QF;
        delete[] TrQF;
        delete[] YF;
    }

}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::get_Q_matrix_omp_prepare(TCSR<CPX> *A,TCSR<CPX> *B,\
                                                 LinearSolver<CPX> **solver_out,\
                                                 LinearSolver<CPX> **solver_in,\
                                                 int NM,int N,MPI_Comm comm)
{

    int IR;
    TCSR<CPX> *C;

    for(IR=0;IR<NQDR;IR++){
        C = new TCSR<CPX>(CPX(1.0,0.0),A,-z_out[IR],B);
	solver_out[IR] = new Umfpack<CPX>(C,comm);
        delete C;
    }

#pragma omp parallel for
    for(IR=0;IR<NQDR;IR++){
	solver_out[IR]->prepare();
    }

    for(IR=0;IR<NQDR;IR++){
        C = new TCSR<CPX>(CPX(1.0,0.0),A,-z_in[IR],B);
	solver_in[IR] = new Umfpack<CPX>(C,comm);
        delete C;
    }

#pragma omp parallel for
    for(IR=0;IR<NQDR;IR++){
	solver_in[IR]->prepare();
    }

}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::get_Q_matrix_omp_delete(LinearSolver<CPX> **solver_out,\
                                                LinearSolver<CPX> **solver_in)
{

    for(int IR=0;IR<NQDR;IR++){
        delete solver_out[IR];
        delete solver_in[IR];
    }

}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::get_Q_matrix_omp_solver(CPX *Q,CPX *Y,LinearSolver<CPX> **solver_out,\
                                                LinearSolver<CPX> **solver_in,int NM,int N)
{

    int IR;
    CPX **QCol = new CPX*[NQDR];

    init_var(Q,N*NM);

    for(IR=0;IR<NQDR;IR++){
        QCol[IR] = new CPX[N*NM];
    }

#pragma omp parallel for
    for(IR=0;IR<NQDR;IR++){
        solver_out[IR]->solve_equation(QCol[IR],Y,NM);
    }

    for(IR=0;IR<NQDR;IR++){
        c_zaxpy(N*NM,CPX(0.0,1.0/(2.0*PI))*dz_dtheta_out[IR]*dtheta[IR],QCol[IR],1,Q,1);
    }

#pragma omp parallel for
    for(IR=0;IR<NQDR;IR++){
        solver_in[IR]->solve_equation(QCol[IR],Y,NM);
    }

    for(IR=0;IR<NQDR;IR++){
        c_zaxpy(N*NM,CPX(0.0,-1.0/(2.0*PI))*dz_dtheta_in[IR]*dtheta[IR],QCol[IR],1,Q,1);
    }

    for(IR=0;IR<NQDR;IR++){
        delete[] QCol[IR];
    }

    delete[] QCol;

}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::get_Q_matrix_omp(CPX *Q,CPX *Y,TCSR<CPX> *A,TCSR<CPX> *B,\
                                         int NM,int N,MPI_Comm comm)
{
    int n_threads=omp_get_num_threads();
    int IT,IR,IP;

    CPX **QCol = new CPX*[n_threads];
    LinearSolver<CPX> **solver = new LinearSolver<CPX>*[n_threads];
    TCSR<CPX> *C;

    init_var(Q,N*NM);

    for(IT=0;IT<int(ceil(double(NQDR)/n_threads));IT++){

        for(IR=0;IR<n_threads;IR++){
            if ((IP=IR+IT*n_threads)<NQDR){
                QCol[IR] = new CPX[N*NM];
                C = new TCSR<CPX>(CPX(1.0,0.0),A,-z_out[IP],B);
                solver[IR] = new Umfpack<CPX>(C,comm);
                delete C;
            }
        }

#pragma omp parallel for
        for(IR=0;IR<n_threads;IR++){
            if ((IP=IR+IT*n_threads)<NQDR){
                solver[IR]->prepare();
                solver[IR]->solve_equation(QCol[IR],Y,NM);
            }
        }
 
        for(IR=0;IR<n_threads;IR++){
            if ((IP=IR+IT*n_threads)<NQDR){
                delete solver[IR];
                c_zaxpy(N*NM,CPX(0.0,1.0/(2.0*PI))*dz_dtheta_out[IP]*dtheta[IP],QCol[IR],1,Q,1);
                delete[] QCol[IR];
            }
        }
 
        for(IR=0;IR<n_threads;IR++){
            if ((IP=IR+IT*n_threads)<NQDR){
                QCol[IR] = new CPX[N*NM];
                C = new TCSR<CPX>(CPX(1.0,0.0),A,-z_in[IP],B);
                solver[IR] = new Umfpack<CPX>(C,comm);
                delete C;
            }
        }

#pragma omp parallel for
        for(IR=0;IR<n_threads;IR++){
            if ((IP=IR+IT*n_threads)<NQDR){
                solver[IR]->prepare();
                solver[IR]->solve_equation(QCol[IR],Y,NM);
            }
        }
 
        for(IR=0;IR<n_threads;IR++){
            if ((IP=IR+IT*n_threads)<NQDR){
                delete solver[IR];
                c_zaxpy(N*NM,CPX(0.0,-1.0/(2.0*PI))*dz_dtheta_in[IP]*dtheta[IP],QCol[IR],1,Q,1);
                delete[] QCol[IR];
            }
        }

    }

    delete[] QCol;
    delete[] solver;

}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::get_Q_matrix_mpi(CPX *Q,CPX *Y,TCSR<CPX> *A,TCSR<CPX> *B,\
				     int NM,int N,MPI_Comm comm)
{
    int iam,nprocs;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&iam);
    int NT=int(ceil(double(NQDR)/nprocs));

    CPX *QCol = new CPX[N*NM];
    CPX *Qtmp;
    LinearSolver<CPX> *solver;
    TCSR<CPX> *C;

    if (!iam){
        init_var(Q,N*NM);
        Qtmp = new CPX[N*NM];
    }

    int A_n_nonzeros;
    int B_n_nonzeros;
    int AB_findx;
    if (!iam) {
        A_n_nonzeros=A->n_nonzeros;
        B_n_nonzeros=B->n_nonzeros;
        AB_findx=A->findx;
    }
    MPI_Bcast(&A_n_nonzeros,1,MPI_INT,0,comm);
    MPI_Bcast(&B_n_nonzeros,1,MPI_INT,0,comm);
    MPI_Bcast(&AB_findx,1,MPI_INT,0,comm);

    for(int IT=0;IT<NT;IT++){
        int IR=iam+IT*nprocs;
        if(iam){
            A = new TCSR<CPX>(N,A_n_nonzeros,AB_findx);
            B = new TCSR<CPX>(N,B_n_nonzeros,AB_findx);
        }
        A->Bcast(0,comm);
        B->Bcast(0,comm);
        if(IR<NQDR){
            C = new TCSR<CPX>(CPX(1.0,0.0),A,-z_out[IR],B);
        }
        if(iam){
            delete A;
            delete B;
        }
        if(IR<NQDR){
 	    solver = new Umfpack<CPX>(C,comm);
            delete C;
 	    solver->prepare();
        }
        if(iam){
            Y = new CPX[NM*N];
        }
        MPI_Bcast(Y,NM*N,MPI_DOUBLE_COMPLEX,0,comm);
        if(IR<NQDR){
            solver->solve_equation(QCol,Y,NM);
 	    delete solver;
            c_zscal(N*NM,CPX(0.0,1.0/(2.0*PI))*dz_dtheta_out[IR]*dtheta[IR],QCol,1);
        }else{
            c_zscal(N*NM,CPX(0.0,0.0),QCol,1);
        }
        if(iam){
            delete[] Y;
        }
        MPI_Reduce(QCol,Qtmp,N*NM,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm);
        if(!iam){
            c_zaxpy(N*NM,CPX(1.0,0.0),Qtmp,1,Q,1);
        }
        if(iam){
            A = new TCSR<CPX>(N,A_n_nonzeros,AB_findx);
            B = new TCSR<CPX>(N,B_n_nonzeros,AB_findx);
        }
        A->Bcast(0,comm);
        B->Bcast(0,comm);
        if(IR<NQDR){
            C = new TCSR<CPX>(CPX(1.0,0.0),A,-z_in[IR],B);
        }
        if(iam){
            delete A;
            delete B;
        }
        if(IR<NQDR){
 	    solver = new Umfpack<CPX>(C,comm);
            delete C;
 	    solver->prepare();
        }
        if(iam){
            Y = new CPX[NM*N];
        }
        MPI_Bcast(Y,NM*N,MPI_DOUBLE_COMPLEX,0,comm);
        if(IR<NQDR){
            solver->solve_equation(QCol,Y,NM);
 	    delete solver;
            c_zscal(N*NM,CPX(0.0,-1.0/(2.0*PI))*dz_dtheta_in[IR]*dtheta[IR],QCol,1);
        }else{
            c_zscal(N*NM,CPX(0.0,0.0),QCol,1);
        }
        if(iam){
            delete[] Y;
        }
        MPI_Reduce(QCol,Qtmp,N*NM,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm);
        if(!iam){
            c_zaxpy(N*NM,CPX(1.0,0.0),Qtmp,1,Q,1);
        }
    }

    delete[] QCol;
    if (!iam){
        delete[] Qtmp;
    }
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::get_Q_matrix(CPX *Q,CPX *Y,TCSR<CPX> *A,TCSR<CPX> *B,\
				     int NM,int N,MPI_Comm comm)
{

    int IR,IC;
    CPX *QCol = new CPX[N];
    LinearSolver<CPX> *solver;
    TCSR<CPX> *C;

    init_var(Q,N*NM);

    for(IR=0;IR<NQDR;IR++){

        C = new TCSR<CPX>(CPX(1.0,0.0),A,-z_out[IR],B);
	solver = new Umfpack<CPX>(C,comm);
        delete C;
 
	solver->prepare();
 
	for(IC=0;IC<NM;IC++){
	    solver->solve_equation(QCol,&Y[N*IC],1);
	    c_zaxpy(N,CPX(0.0,1.0/(2.0*PI))*dz_dtheta_out[IR]*dtheta[IR],QCol,1,&Q[N*IC],1);
	}
	delete solver;

        C = new TCSR<CPX>(CPX(1.0,0.0),A,-z_in[IR],B);
	solver = new Umfpack<CPX>(C,comm);
        delete C;

	solver->prepare();

	for(IC=0;IC<NM;IC++){
	    solver->solve_equation(QCol,&Y[N*IC],1);
	    c_zaxpy(N,CPX(0.0,-1.0/(2.0*PI))*dz_dtheta_in[IR]*dtheta[IR],QCol,1,&Q[N*IC],1);
	}
	delete solver;

    }

    delete[] QCol;
}

/************************************************************************************************/

template <class T>
void InjectionFeast<T>::init_Y_matrix(CPX *Y,int NM,int N)
{
    int IR,IC;
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

#endif
