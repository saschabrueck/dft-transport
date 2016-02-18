#ifndef __INJECTIONBEYN
#define __INJECTIONBEYN

#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "Types.H"
#include "Utilities.H"
#include "Blas.H"
#include "Connection.H"
#include "LinearSolver.H"
#include "Umfpack.H"
#include "Injection.H"
#include "p_eig.H"

template<class T>
class InjectionBeyn : public Injection<T>{
public:
    
    InjectionBeyn(int,double);

    virtual ~InjectionBeyn();

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
    int execute(TCSR<CPX>**,int,int,CPX*,CPX*,int,MPI_Comm);
    
private:

    double mult_factor;
    double svd_fac;
    double contour_radius;
    double imag_limit;
    int mpi_size,mpi_rank;
    int loc_size,loc_rank;
    int NM_input;
    int slab_per_bc;
    int NQDR;
    int NCRC;
    int parallel_svd;

    void sort_eig(CPX*,CPX*,double*,int*,int*,int*,int*,TCSR<CPX>*,int,int,int);
    int get_type(TCSR<CPX>*,CPX*,CPX,double*,int,int);
    double calc_derivative(TCSR<CPX>*,CPX*,CPX);
    void init_matrix_single(TCSR<CPX>*,TCSR<CPX>*,TCSR<T>*,TCSR<T>*,int,int);
    void init_matrix_double(TCSR<CPX>*,TCSR<CPX>*,TCSR<T>*,TCSR<T>*,int,int);
    int parallel_decompose(CPX*,CPX*,int,int,MPI_Comm);
    int decompose(CPX*,CPX*,int,int);
    int svd(CPX*,CPX*,double*,int,int);
    int eig(CPX*,int,CPX*);
    void allocate_bcast_H(TCSR<CPX>**&,int,MPI_Comm);
    void get_Q_matrix_int(CPX*,CPX*,CPX*,TCSR<CPX>**,int,int,int,MPI_Comm);
    void get_Q_matrix(CPX*,CPX*,CPX*,TCSR<CPX>**,int,int,int,int,MPI_Comm);
    void create_M_matrix(CPX*,TCSR<CPX>**,int,CPX);
    void solve_Q_sparse(CPX*,CPX*,int,int,CPX,MPI_Comm);
    void solve_Q_full(CPX*,CPX*,int,int,CPX,MPI_Comm);
    void inv_Q_full(CPX*,CPX*,int,CPX,MPI_Comm);
    void init_Y_matrix(CPX*,int,int);
};

template<class T>
InjectionBeyn<T>::InjectionBeyn(int pslab_per_bc,double pcontour_radius)
{
    imag_limit      = 5.0e-5; //eps_decay for OMEN part
    mult_factor     = 1.2;
    svd_fac         = 1.0;
    NQDR            = 32;
    NCRC            = 2;
    parallel_svd    = 1;

    slab_per_bc    = pslab_per_bc;
    contour_radius = pcontour_radius;
    imag_limit     = slab_per_bc*imag_limit;

    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
}

/************************************************************************************************/

template <class T>
InjectionBeyn<T>::~InjectionBeyn()
{
}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::initialize(int ND,int NP,int NM)
{

    NM_input = min(2*NM,ND);
 
}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::sort_eig(CPX *kphase,CPX *Vsurf,double *dEk_dk,int *Ntr,int *ind_Ntr,\
				 int *Nref,int *ind_Nref,TCSR<CPX> *T01,int ND,int NM,int sign)
{

    int IP,IM;
    int type;

    int NR = T01->size/slab_per_bc;

    *Ntr   = 0;
    *Nref  = 0;    

    for(IM=0;IM<NM;IM++){

        for (IP=1;IP<slab_per_bc;IP++){
            c_zcopy(NR,&Vsurf[IM*ND],1,&Vsurf[IM*ND+IP*NR],1);
            c_zscal(NR,pow(kphase[IM],-IP),&Vsurf[IM*ND+IP*NR],1);
        }

        kphase[IM] = CPX(0.0,1.0)*log(pow(kphase[IM],slab_per_bc));
	  
        type = get_type(T01,&Vsurf[IM*ND],kphase[IM],&dEk_dk[(*Ntr)],0,sign);

        if((type>0)&&(dEk_dk[(*Ntr)]>0)){
            ind_Ntr[(*Ntr)] = IM;
            (*Ntr)          = (*Ntr)+1;
        }else if (type<0){
            ind_Nref[(*Nref)] = IM;
            (*Nref)           = (*Nref)+1;
        }
    }
}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::calc_kphase(TCSR<T> *D,TCSR<T> *P,TCSR<CPX> *T01,int ND,int NP,CPX *kphase,\
				    CPX *Vsurf,double *dEk_dk,int *Ntr,int *ind_Ntr,int *Nref, \
				    int *ind_Nref,int sign,int NPROW,int NPCOL,MPI_Comm comm,\
				    int *info)
{
    int IM;
    int NR;
    int NM = NM_input;

    int inp_rank;
    MPI_Comm_rank(comm,&inp_rank);
    MPI_Comm loc_comm;
    MPI_Comm_split(comm,0,sign*inp_rank,&loc_comm);
    MPI_Comm_size(loc_comm,&loc_size);
    MPI_Comm_rank(loc_comm,&loc_rank);

    if (!loc_rank) NR=D->size/slab_per_bc;
    MPI_Bcast(&NR,1,MPI_INT,0,loc_comm);
    CPX *YF      = NULL;
    CPX *lambdaR = NULL;
    if(!loc_rank){
        YF      = new CPX[NR*NM];
        lambdaR = new CPX[NR];
    }
    TCSR<CPX> **H = NULL;

    double t0 = get_time(0.0);

    if(!loc_rank){

        TCSR<CPX> *AF = new TCSR<CPX>(D->size,D->n_nonzeros,D->findx);
        TCSR<CPX> *BF = new TCSR<CPX>(D->size,D->n_nonzeros,D->findx);

        if(slab_per_bc==1){
            init_matrix_single(AF,BF,D,P,ND,NP);
        }else{
            init_matrix_double(AF,BF,D,P,ND,NP);
        }

        H = new TCSR<CPX>*[slab_per_bc+1];
        H[0] = new TCSR<CPX>(NR,BF->n_nonzeros,BF->findx);
        H[0]->n_nonzeros=0;
        init_var(H[0]->index_i,NR);
        for(int iter=1;iter<=slab_per_bc;iter++){
            H[iter] = new TCSR<CPX>(NR,AF->n_nonzeros,AF->findx);
            H[iter]->n_nonzeros=0;
            init_var(H[iter]->index_i,NR);
        }
        for(int IR=0;IR<NR;IR++){
            for(int IC=BF->edge_i[IR]-BF->findx;IC<BF->edge_i[IR+1]-BF->findx;IC++){
                H[0]->index_i[IR]++;
                H[0]->index_j[H[0]->n_nonzeros]=BF->index_j[IC]-BF->findx+H[0]->findx;
                H[0]->nnz[H[0]->n_nonzeros++]=-BF->nnz[IC];
            }
        }
        for(int IR=0;IR<NR;IR++){
            for(int IC=AF->edge_i[IR]-AF->findx;IC<AF->edge_i[IR+1]-AF->findx;IC++){
                int col_ind = AF->index_j[IC]-AF->findx;
                int mat_ind = col_ind/NR+1;
                H[mat_ind]->index_i[IR]++;
                H[mat_ind]->index_j[H[mat_ind]->n_nonzeros]=col_ind%NR+H[mat_ind]->findx;
                H[mat_ind]->nnz[H[mat_ind]->n_nonzeros++]=AF->nnz[IC];
            }
        }
        for(int iter=0;iter<=slab_per_bc;iter++){
            H[iter]->get_row_edge();
            H[iter]->get_diag_pos();
        }
 
        delete AF;
        delete BF;

    }

    if(!mpi_rank) printf("Time to generate matrices: %e\n",get_time(t0));

    NM=execute(H,NR,NM,lambdaR,YF,sign,loc_comm);

    t0 = get_time(0.0);

    if(!loc_rank){

        int *reord_index = new int[NM];

        sort_abs_vec(lambdaR,reord_index,NM);

        for(IM=0;IM<NM;IM++){
            c_zcopy(NR,&YF[IM*NR],1,&Vsurf[reord_index[IM]*ND],1);
        }
        c_zcopy(NM,lambdaR,1,kphase,1);

        delete[] reord_index;

        sort_eig(kphase,Vsurf,dEk_dk,Ntr,ind_Ntr,Nref,ind_Nref,T01,ND,NM,sign);

    }

    if(!mpi_rank) printf("Sort time: %e\n",get_time(t0));

    if (!loc_rank){
        delete[] YF;
        delete[] lambdaR;
        for(int IP=0;IP<slab_per_bc+1;IP++){
            delete H[IP];
        }
        delete[] H;
    }
    MPI_Comm_free(&loc_comm);
 
}

/************************************************************************************************/

template <class T>
int InjectionBeyn<T>::execute(TCSR<CPX> **A,int NR,int NM,CPX *lambdaR,CPX *YF, \
                               int sign,MPI_Comm comm)
{

    int IP,IM;
    double LBound,UBound;

    MPI_Comm_size(comm,&loc_size);
    MPI_Comm_rank(comm,&loc_rank);

    if(sign>0){
        UBound = 1.0+1.0/contour_radius;
        LBound = 1.0/contour_radius;
    }else{
        UBound = contour_radius;
        LBound = 1.0/(1.0+1.0/contour_radius);
    }

    CPX *QF1 = NULL;
    CPX *QF2 = NULL;

    if(!loc_rank){
        init_Y_matrix(YF,NM,NR);
        QF1 = new CPX[NR*NM];
        QF2 = new CPX[NR*NM];
    }

    double t0 = get_time(0.0);

    get_Q_matrix_int(QF1,QF2,YF,A,sign,NR,NM,comm);

    if(!mpi_rank) printf("Time to generate Q: %e\n",get_time(t0));

    t0 = get_time(0.0);

    if(parallel_svd){

        NM=parallel_decompose(QF1,QF2,NR,NM,comm);

    }else if(!loc_rank){

        NM=decompose(QF1,QF2,NR,NM);

    }

    if(!mpi_rank) printf("Time for SVD: %e\n",get_time(t0));

    t0 = get_time(0.0);

    if(!loc_rank){

        eig(QF2,NM,lambdaR);

        c_zgemm('N','N',NR,NM,NM,CPX(1.0,0.0),QF1,NR,QF2,NM,CPX(0.0,0.0),YF,NR);

        IP = 0;
        for(IM=0;IM<NM;IM++){
            if((abs(lambdaR[IM])>=LBound)&&(abs(lambdaR[IM])<=UBound)){
                if(IP<IM){
                    c_zcopy(NR,&YF[IM*NR],1,&YF[IP*NR],1);
                    lambdaR[IP]=lambdaR[IM];
                }
                IP++;
            }
        }
        NM=IP;

    }

    if(!mpi_rank) printf("Sequential time: %e\n",get_time(t0));

    if(!loc_rank){
        delete[] QF1;
        delete[] QF2;
    }

    return NM;

}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::calc_kphase_sym(TCSR<T> *D,TCSR<T> *P,TCSR<CPX> *T01,int ND,int NP,\
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
void InjectionBeyn<T>::calc_kphase(TCSR<T> *D,T *P,TCSR<CPX> *T01,int ND,int NP,CPX *kphase,\
				    CPX *Vsurf,double *dEk_dk,int *Ntr,int *ind_Ntr,int *Nref, \
				    int *ind_Nref,int sign,int NPROW,int NPCOL,MPI_Comm comm,\
				    int *info)
{
}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::calc_kphase_sym(TCSR<T> *D,T *P,TCSR<CPX> *T01,int ND,int NP,\
					CPX *kphasel,CPX *kphaser,CPX *Vsurfl,CPX *Vsurfr, \
					double *dEkl_dk, double *dEkr_dk, int *Ntrl, int *Ntrr,	\
					int *ind_Ntrl,int *ind_Ntrr,int *Nrefl,int *Nrefr, \
					int *ind_Nrefl,int *ind_Nrefr,int NPROW,int NPCOL, \
					MPI_Comm comm,int *info)
{
}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::allocate_bcast_H(TCSR<CPX> **&H,int broot_rank,MPI_Comm bcast_comm)
{
    int bcast_rank;
    MPI_Comm_rank(bcast_comm,&bcast_rank);
    int *nonzH = new int[slab_per_bc+1];
    int findxH;
    int sizeH;
    if(!bcast_rank){
        for (int ibw=0;ibw<slab_per_bc+1;ibw++) {
            nonzH[ibw] = H[ibw]->n_nonzeros;
        }
        findxH = H[0]->findx;
        sizeH  = H[0]->size;
    }
    MPI_Bcast(nonzH,slab_per_bc+1,MPI_INT,broot_rank,bcast_comm);
    MPI_Bcast(&findxH,1,MPI_INT,broot_rank,bcast_comm);
    MPI_Bcast(&sizeH,1,MPI_INT,broot_rank,bcast_comm);
    if(bcast_rank){
        H = new TCSR<CPX>*[slab_per_bc+1];
        for(int ibw=0;ibw<slab_per_bc+1;ibw++){
            H[ibw] = new TCSR<CPX>(sizeH,nonzH[ibw],findxH);
        }
    }
    delete[] nonzH;
    for(int ibw=0;ibw<slab_per_bc+1;ibw++){
        MPI_Bcast_Sparse(H[ibw],broot_rank,bcast_comm);
    }
}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::get_Q_matrix(CPX *Q,CPX *Q2,CPX *Y,TCSR<CPX> **A,int sign,int NR,int NM,\
                                     int NPOC,MPI_Comm comm)
{

    int IR;
    CPX z;
    double R_in,R_out;

    if(slab_per_bc==1) NCRC=2; //NOT IMPLEMENTED FOR SPARSE SOLVER

    int size_lin_comm = max(1,loc_size/(NCRC*NPOC));
if (!mpi_rank) cout << "SIZE OF COMM FOR LINEAR SYSTEM " << size_lin_comm << endl;
if (!mpi_rank) cout << "NUMBER OF INTEGRATION POINTS " << NCRC*NPOC << endl;
if (!mpi_rank) cout << "SIZE OF COMM PER BOUNDARY CONDITION " << loc_size << endl;

    R_out = mult_factor*contour_radius;
    if(NCRC==2){
        if(sign>0){
            R_out = 1.0+mult_factor/contour_radius;
            R_in  = 1.0/(mult_factor*contour_radius);
        }else{
            R_out = mult_factor*contour_radius;
            R_in  = 1.0/(1.0+mult_factor/contour_radius);
        }
    }

    if(size_lin_comm>1){

        IR = loc_rank / size_lin_comm;
        MPI_Comm lin_comm;
        MPI_Comm_split(comm,IR,loc_rank,&lin_comm);
        int rank_lin_comm;
        MPI_Comm_rank(lin_comm,&rank_lin_comm);
        MPI_Comm master_comm;
        int color = rank_lin_comm;
        if (loc_rank>=size_lin_comm*NCRC*NPOC) {
            color = MPI_UNDEFINED;
            rank_lin_comm = -1;
        }
        MPI_Comm_split(comm,color,loc_rank,&master_comm);

        CPX *QLoc  = NULL;
        CPX *QLoc2 = NULL;
        if(!rank_lin_comm){
            QLoc  = new CPX[NR*NR];
            QLoc2 = new CPX[NR*NM];
        }

        if(!loc_rank){
            c_zcopy(NM*NR,Y,1,QLoc2,1);
        }
        if(!rank_lin_comm){
            allocate_bcast_H(A,0,master_comm);
            MPI_Bcast(QLoc2,NM*NR,MPI_DOUBLE_COMPLEX,0,master_comm);
        }

        if(IR<NPOC){
            z = R_out*exp(CPX(0.0,2.0*(IR+0.5)*M_PI/NQDR));
	}else if(IR>=NPOC && IR<NCRC*NPOC){
            z = R_in*exp(CPX(0.0,2.0*((IR-NPOC)+0.5)*M_PI/NQDR));
	}

        if(IR<NCRC*NPOC){

            if(!rank_lin_comm){
                create_M_matrix(QLoc,A,NR,z);
            }

	    if(slab_per_bc==1){
                solve_Q_sparse(QLoc,QLoc2,NR,NM,z,lin_comm);
  	    }else if(NR==NM){
  	        inv_Q_full(QLoc,QLoc2,NM,z,lin_comm);
	    }else{
	        solve_Q_full(QLoc,QLoc2,NR,NM,z,lin_comm);
	    }

        }

        if(!rank_lin_comm){
            MPI_Reduce(QLoc ,Q ,NR*NM,MPI_DOUBLE_COMPLEX,MPI_SUM,0,master_comm);
            MPI_Reduce(QLoc2,Q2,NR*NM,MPI_DOUBLE_COMPLEX,MPI_SUM,0,master_comm);
            delete[] QLoc;
            delete[] QLoc2;
            if(loc_rank){
                for(int IP=0;IP<slab_per_bc+1;IP++){
                    delete A[IP];
                }
                delete[] A;
            }
        }

        MPI_Comm_free(&lin_comm);
        MPI_Comm_free(&master_comm);

    }else{

        CPX *QLoc  = new CPX[NR*NR];
        CPX *QLoc2 = new CPX[NR*NM];
        CPX *QRed  = NULL;
        if(!loc_rank){
            init_var(Q ,NR*NM);
            init_var(Q2,NR*NM);
            QRed  = new CPX[NR*NM];
        }
        if(loc_size>1){
            allocate_bcast_H(A,0,comm);
        }

        for(IR=loc_rank;IR<NCRC*NPOC;IR=IR+loc_size){

            if(!loc_rank){
                c_zcopy(NM*NR,Y,1,QLoc2,1);
            }
            if(loc_size>1){
                MPI_Bcast(QLoc2,NM*NR,MPI_DOUBLE_COMPLEX,0,comm);
            }
 
            if(IR<NPOC){
                z = R_out*exp(CPX(0.0,2.0*(IR+0.5)*M_PI/NQDR));
            }else if(IR>=NPOC && IR<NCRC*NPOC){
                z = R_in*exp(CPX(0.0,2.0*((IR-NPOC)+0.5)*M_PI/NQDR));
            }

            if(IR<NCRC*NPOC){
 
                create_M_matrix(QLoc,A,NR,z);
                if(slab_per_bc==1){
                    solve_Q_sparse(QLoc,QLoc2,NR,NM,z,MPI_COMM_SELF);
                }else{
                    solve_Q_full(QLoc,QLoc2,NR,NM,z,MPI_COMM_SELF);
                }
 
	    }else{

                init_var(QLoc ,NR*NM);
                init_var(QLoc2,NR*NM);

            }

            MPI_Reduce(QLoc ,QRed,NR*NM,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm);
            if(!loc_rank){
                c_zaxpy(NR*NM,CPX(1.0,0.0),QRed,1,Q ,1);
            }
            MPI_Reduce(QLoc2,QRed,NR*NM,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm);
            if(!loc_rank){
                c_zaxpy(NR*NM,CPX(1.0,0.0),QRed,1,Q2,1);
            }
 
        }

        delete[] QLoc;
        delete[] QLoc2;
        if(!loc_rank){
            delete[] QRed;
        }
        if(loc_rank){
            for(int IP=0;IP<slab_per_bc+1;IP++){
                delete A[IP];
            }
            delete[] A;
        }

    }

}

/************************************************************************************************/

template <>
void InjectionBeyn<CPX>::get_Q_matrix_int(CPX *Q,CPX *Q2,CPX *Y,TCSR<CPX> **A,int sign,\
                                           int NR,int NM,MPI_Comm comm)
{

    get_Q_matrix(Q,Q2,Y,A,sign,NR,NM,NQDR,comm);

}

/************************************************************************************************/

template <>
void InjectionBeyn<double>::get_Q_matrix_int(CPX *Q,CPX *Q2,CPX *Y,TCSR<CPX> **A,int sign,\
                                              int NR,int NM,MPI_Comm comm)
{

    if(!loc_rank){
        c_dscal(NR*NM,0.0,((double*)Y)+1,2);
    }

    get_Q_matrix(Q,Q2,Y,A,sign,NR,NM,NQDR/2,comm);

    if(!loc_rank){
        c_dscal(NR*NM,2.0, (double*)Q    ,2);
        c_dscal(NR*NM,0.0,((double*)Q )+1,2);
        c_dscal(NR*NM,2.0, (double*)Q2   ,2);
        c_dscal(NR*NM,0.0,((double*)Q2)+1,2);
    }

}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::create_M_matrix(CPX *M,TCSR<CPX> **A,int NR,CPX z)
{
    init_var(M,NR*NR);
 
    for(int IP=0;IP<slab_per_bc+1;IP++){
        TCSR<CPX>* AA = A[IP];
        int ffindx = AA->findx;
        CPX zz = pow(z,slab_per_bc/2-IP);
        for(int IR=0;IR<NR;IR++){
            for(int IJ=AA->edge_i[IR]-ffindx;IJ<AA->edge_i[IR+1]-ffindx;IJ++){
                M[IR+(AA->index_j[IJ]-ffindx)*NR] += zz*AA->nnz[IJ];
            }
        }
    }
}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::solve_Q_sparse(CPX *QM,CPX *QY,int NR,int NM,CPX z,MPI_Comm lin_comm)
{
    int lin_size,lin_rank;
    MPI_Comm_size(lin_comm,&lin_size);
    MPI_Comm_rank(lin_comm,&lin_rank);

    if(!lin_rank){
 
        TCSR<CPX> *M_sp = new TCSR<CPX>(NR,NR*NR,1);
        M_sp->full_to_sparse(QM,NR,NR);
 
        LinearSolver<CPX> *M_solver = new Umfpack<CPX>(M_sp,MPI_COMM_SELF);
        M_solver->prepare();
        M_solver->solve_equation(QM,QY,NM);
 
        delete M_solver;
        delete M_sp;
 
        if(abs(z)<1.0) c_zscal(NR*NM,-1.0,QM,1);
        c_zcopy(NR*NM,QM,1,QY,1);
        c_zscal(NR*NM,z,QY,1);

    }

}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::solve_Q_full(CPX *QM,CPX *QY,int NR,int NM,CPX z,MPI_Comm lin_comm)
{
    int lin_size,lin_rank;
    MPI_Comm_size(lin_comm,&lin_size);
    MPI_Comm_rank(lin_comm,&lin_rank);

    if(NCRC==1){
        p_lin(QM,QY,NR,NM,z,lin_comm);
    }else{
        if(lin_size==1){
            int info;
            int *ipiv = new int[NR];
            c_zgetrf(NR,NR,QM,NR,ipiv,&info);
            c_zgetrs('N',NR,NM,QM,NR,ipiv,QY,NR,&info);
            delete[] ipiv;
            c_zcopy(NR*NM,QY,1,QM,1);
        }else{
            p_lin(QM,QY,QM,NR,NM,lin_comm);
        }
        if(!lin_rank){
            if(abs(z)<1.0) c_zscal(NR*NM,-1.0,QM,1);
            c_zcopy(NR*NM,QM,1,QY,1);
            c_zscal(NR*NM,z,QY,1);
        }
    }

}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::inv_Q_full(CPX *QM,CPX *QY,int NR,CPX z,MPI_Comm lin_comm)
{
    int lin_rank;
    MPI_Comm_rank(lin_comm,&lin_rank);

    if(NCRC==1){
        p_inv(QM,QY,NR,z,lin_comm);
    }else{
        p_inv(QM,NR,lin_comm);
        if(!lin_rank){
            if(abs(z)<1.0) c_zscal(NR*NR,-1.0,QM,1);
            c_zcopy(NR*NR,QM,1,QY,1);
            c_zscal(NR*NR,z,QY,1);
        }
    }

}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::init_Y_matrix(CPX *Y,int NM,int N)
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
void InjectionBeyn<T>::init_matrix_single(TCSR<CPX> *A,TCSR<CPX> *B,TCSR<T> *D,\
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

}

/************************************************************************************************/

template <class T>
void InjectionBeyn<T>::init_matrix_double(TCSR<CPX> *A,TCSR<CPX> *B,TCSR<T> *D,\
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

    delete[] line;
}

/************************************************************************************************/

template <class T>
int InjectionBeyn<T>::get_type(TCSR<CPX> *T01,CPX* Vin,CPX kin,double *dEk_dk,int ind_to_keep,\
				int sign)
{
    int type;
    //double dEk;

    *dEk_dk = 0.0;
    
    if(((imag(kin)<imag_limit)&&(imag(kin)>-imag_limit))){
        
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
double InjectionBeyn<T>::calc_derivative(TCSR<CPX> *T01,CPX* Vin,CPX kin)
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
int InjectionBeyn<T>::parallel_decompose(CPX *Q1,CPX *Q2,int ND,int NM, MPI_Comm comm)
{

    double *Sval = new double[NM];

    int icontxt=MPI_Comm_c2f(comm);
    int nprowcol[2]={0,0};
    MPI_Dims_create(loc_size,2,nprowcol);
    int nprow = nprowcol[0];
    int npcol = nprowcol[1];
    int myrow, mycol;
    char gridr[1] = {'R'};
    Cblacs_gridinit(&icontxt,gridr,nprow,npcol);
    Cblacs_gridinfo(icontxt,&nprow,&npcol,&myrow,&mycol);
    int row_per_processor    = ND/nprow;
    int block_per_rprocessor = int(ceil(double(row_per_processor)/64));
    int mb                   = row_per_processor/block_per_rprocessor;

    int col_per_processor    = NM/npcol;
    int block_per_cprocessor = int(ceil(double(col_per_processor)/64));
    int nb                   = col_per_processor/block_per_cprocessor;

    nb=(nb+mb)/2;
    mb=nb;

    int rloc      = max(1,c_numroc(ND,mb,myrow,0,nprow));
    int cloc      = c_numroc(NM,nb,mycol,0,npcol);

    int info;
    int descQ[9];
    c_descinit(descQ,ND,NM,ND,NM,0,0,icontxt,ND,&info);
    int descq[9]; //LOCAL MATRICES ALL (ND,NM) EVEN IF LESS SPACE NEEDED
    c_descinit(descq,ND,NM,mb,nb,0,0,icontxt,rloc,&info);

    CPX *q_dist      = new CPX[2*rloc*cloc];
    CPX *u_dist      = new CPX[2*rloc*cloc];
    CPX *vtrans_dist = new CPX[2*rloc*cloc];

    c_pzgeadd('N',ND,NM,1.0,Q1,1,1,descQ,0.0,q_dist,1,1,descq);
    
    int lwork     = -1;
    CPX workq;
    double rworkq;
    CPX *work     = &workq;
    double *rwork = &rworkq;

    c_pzgesvd('V','V',ND,NM,q_dist,1,1,descq,Sval,u_dist,1,1,descq,vtrans_dist,1,1,descq,work,lwork,rwork,&info);

    lwork = 2*int(real(workq));
    work  = new CPX[lwork];
    rwork = new double[2*int(rworkq)];

    c_pzgesvd('V','V',ND,NM,q_dist,1,1,descq,Sval,u_dist,1,1,descq,vtrans_dist,1,1,descq,work,lwork,rwork,&info);

    delete[] work;
    delete[] rwork;

    c_pzgeadd('N',ND,NM,1.0,u_dist,1,1,descq,0.0,Q1,1,1,descQ);
    c_pzgeadd('N',ND,NM,1.0,Q2,1,1,descQ,0.0,q_dist,1,1,descq);

    int NM_new = 0;
    while(Sval[NM_new+1]>svd_fac*Sval[0]*ND*(numeric_limits<double>::epsilon)()) NM_new++;
if (NM_new>=NM) cout << "BAD ERROR IN BEYN" << endl;

    CPX *UtxQ2_dist = new CPX[rloc*cloc];
    c_pzgemm('C','N',NM_new,NM,ND,CPX(1.0,0.0),u_dist,1,1,descq,q_dist,1,1,descq,CPX(0.0,0.0),UtxQ2_dist,1,1,descq);
    delete[] u_dist;
    c_pzgemm('N','C',NM_new,NM_new,NM,CPX(1.0,0.0),UtxQ2_dist,1,1,descq,vtrans_dist,1,1,descq,CPX(0.0,0.0),q_dist,1,1,descq);
    delete[] UtxQ2_dist;
    delete[] vtrans_dist;

    int descQ_out[9];
    c_descinit(descQ_out,NM_new,NM_new,NM_new,NM_new,0,0,icontxt,NM_new,&info);
    c_pzgeadd('N',NM_new,NM_new,1.0,q_dist,1,1,descq,0.0,Q2,1,1,descQ_out);
    delete[] q_dist;
    if(!loc_rank) for (int inm=0;inm<NM_new;inm++) c_zscal(NM_new,CPX(1.0/Sval[inm],0.0),&Q2[inm*NM_new],1);
    delete[] Sval;

    Cblacs_gridexit(icontxt);

    return NM_new;
}

/************************************************************************************************/

template <class T>
int InjectionBeyn<T>::decompose(CPX *Q1,CPX *Q2,int ND,int NM)
{

    double *Sval = new double[NM];
    CPX *vtrans = new CPX[NM*NM];
    if (svd(Q1,vtrans,Sval,ND,NM)) return 0;
    int NM_new = 0;
    while(Sval[NM_new+1]>svd_fac*Sval[0]*ND*(numeric_limits<double>::epsilon)()) NM_new++;
    CPX *UtxQ2 = new CPX[NM_new*NM];
    c_zgemm('C','N',NM_new,NM,ND,CPX(1.0,0.0),Q1,ND,Q2,ND,CPX(0.0,0.0),UtxQ2,NM_new);
    c_zgemm('N','C',NM_new,NM_new,NM,CPX(1.0,0.0),UtxQ2,NM_new,vtrans,NM,CPX(0.0,0.0),Q2,NM_new);
    delete[] UtxQ2;
    delete[] vtrans;
    for (int inm=0;inm<NM_new;inm++) c_zscal(NM_new,CPX(1.0/Sval[inm],0.0),&Q2[inm*NM_new],1);
    delete[] Sval;
    return NM_new;
}

/************************************************************************************************/

template <class T>
int InjectionBeyn<T>::svd(CPX *u,CPX *vtrans,double *Sval,int ND,int NM)
{
    
    int info;
    int lwork     = -1;
    double *rwork = new double[NM*max(5*NM+7,2*ND+2*NM+1)];
    int *iwork    = new int[8*NM];
    CPX workq;
    CPX *work     = &workq;

    c_zgesdd('O',ND,NM,u,ND,Sval,NULL,1,vtrans,NM,work,lwork,rwork,iwork,&info);

    lwork = int(real(workq));
    work  = new CPX[lwork];

    c_zgesdd('O',ND,NM,u,ND,Sval,NULL,1,vtrans,NM,work,lwork,rwork,iwork,&info);

    delete[] work;
    delete[] rwork;
    delete[] iwork;

    return info;
}

/************************************************************************************************/

template <class T>
int InjectionBeyn<T>::eig(CPX *a,int n,CPX *lambda)
{
    
    char jobvl    = 'N';
    char jobvr    = 'V';
    int lda       = n;
    CPX *vl       = NULL;
    int ldvl      = 1;
    CPX *vr       = new CPX[n*n];
    int ldvr      = n;
    int lwork     = -1;
    int info;
    CPX worktest;
    CPX *work     = &worktest;
    double *rwork = new double[2*n];

    c_zgeev(jobvl,jobvr,n,a,lda,lambda,vl,ldvl,vr,ldvr,work,lwork,rwork,&info);

    lwork = int(real(worktest));
    work  = new CPX[lwork];

    c_zgeev(jobvl,jobvr,n,a,lda,lambda,vl,ldvl,vr,ldvr,work,lwork,rwork,&info);

    delete[] work;
    delete[] rwork;

    c_zcopy(n*n,vr,1,a,1);

    delete[] vr;

    return info;
}

/************************************************************************************************/

#endif
