#ifndef __HYPRE
#define __HYPRE

#include "CSR.H"
#include "Types.H"
#include "LinearSolver.H"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

struct H_Exception{
    H_Exception(const int line,const char* file,const int rank,const int info) {
        char descr[128]="";
        HYPRE_DescribeError(info,descr);
        std::cerr<<"Error "<<descr<<" in line "<<line<<" of file "<<file<<" on rank "<<rank<<std::endl;
    }
};

template <class T>
class Hypre : public LinearSolver<T>{
	
public:

    Hypre(TCSR<T>*,MPI_Comm);

    virtual ~Hypre();
    virtual void prepare();
    virtual void prepare(int*,int*,int,int,int*,int);
    virtual void prepare_corner(CPX*,CPX*,int*,int*,int*,int,int,int*,int);
    virtual void solve_equation(T* res, T* rhs, int no_rhs);
				
private:

    inline MPI_Datatype give_me_MPI_datatype (double*)
    {
        return MPI_DOUBLE;
    }
    inline MPI_Datatype give_me_MPI_datatype (CPX*)
    {
        return MPI_DOUBLE_COMPLEX;
    }

    TCSR<T>* P;
    MPI_Comm internal_comm;

};

/************************************************************************************************/

template <class T>
Hypre<T>::Hypre(TCSR<T>* mat,MPI_Comm solver_comm)
{   

    P = mat;
    internal_comm = solver_comm;

}

/************************************************************************************************/

template <class T>
Hypre<T>::~Hypre()
{
}

/************************************************************************************************/

template <class T>
void Hypre<T>::prepare()
{
}

/************************************************************************************************/

template <class T>
void Hypre<T>::prepare(int *Bmin,int *Bmax,int NBlock,int Bsize,int *n_of_el,int tb)
{
}

/************************************************************************************************/

template <class T>
void Hypre<T>::prepare_corner(CPX *ML_ref,CPX *MR_ref,int *NBC,int *Bmin,int *Bmax,int NBlock,\
			      int Bsize,int *orb_per_at,int tb)
{
}

/************************************************************************************************/

template <class T>
void Hypre<T>::solve_equation(T* res, T* arg_rhs, int no_rhs)
{

    int maxiter=100;
    double solvetol=1e-6;

    int rank,commsize,start_row,start_ele;
    MPI_Comm_rank(internal_comm,&rank);
    MPI_Comm_size(internal_comm,&commsize);
    int *rowindx=new int[P->size];
    for (int i=0;i<P->size;i++) rowindx[i]=i+P->findx+P->first_row;
    int *dist = new int[commsize];
    int *disp = new int[commsize+1]();
    if (P->size!=P->size_tot) {
        MPI_Allgather(&P->size,1,MPI_INT,dist,1,MPI_INT,internal_comm);
        for (int i=0;i<commsize;i++) {
            disp[i+1]=disp[i]+dist[i];
        }
        start_row=0;
        start_ele=0;
    } else {
        int loc_size_max=int(ceil(double(P->size_tot)/double(commsize)));
        int sizepool=P->size_tot;
        for (int i=0;i<commsize;i++) {
            int loc_size_act=loc_size_max;
            if (sizepool<loc_size_max) {
                loc_size_act=sizepool;
            }    
            dist[i]=loc_size_act;
            sizepool-=loc_size_act;
            disp[i+1]=disp[i]+dist[i];
        }
        start_row=disp[rank];
        start_ele=P->edge_i[start_row]-P->findx;
    }
    int ilower=disp[rank]+P->findx;
    int iupper=disp[rank+1]-1+P->findx;

    HYPRE_IJMatrix A;
    HYPRE_IJVector b;
    HYPRE_IJVector x;
    HYPRE_Solver solver;
    HYPRE_Solver precond;
    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;
    int info;

    if ( (info=HYPRE_IJMatrixCreate(internal_comm, ilower, iupper, ilower, iupper, &A)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJMatrixInitialize(A)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJMatrixSetValues(A, dist[rank], &P->index_i[start_row], &rowindx[start_row], &P->index_j[start_ele], &P->nnz[start_ele])) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    delete[] rowindx;
    if ( (info=HYPRE_IJMatrixAssemble(A)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A)) ) throw H_Exception(__LINE__,__FILE__,rank,info);

    if ( (info=HYPRE_IJVectorCreate(internal_comm, ilower, iupper,&x)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorInitialize(x)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorAssemble(x)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorGetObject(x, (void **) &par_x)) ) throw H_Exception(__LINE__,__FILE__,rank,info);

    if ( (info=HYPRE_IJVectorCreate(internal_comm, ilower, iupper,&b)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorInitialize(b)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorSetValues(b, dist[rank], NULL, &arg_rhs[disp[rank]])) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorAssemble(b)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorGetObject(b, (void **) &par_b)) ) throw H_Exception(__LINE__,__FILE__,rank,info);

    if ( (info=HYPRE_ParCSRBiCGSTABCreate(internal_comm, &solver)) ) throw H_Exception(__LINE__,__FILE__,rank,info);

    if ( (info=HYPRE_BiCGSTABSetMaxIter(solver, maxiter)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_BiCGSTABSetTol(solver, solvetol)) ) throw H_Exception(__LINE__,__FILE__,rank,info);

    if ( (info=HYPRE_EuclidCreate(internal_comm, &precond)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_BiCGSTABSetPrecond(solver,(HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,(HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup,precond)) ) throw H_Exception(__LINE__,__FILE__,rank,info);

    if ( (info=HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_b, par_x)) ) throw H_Exception(__LINE__,__FILE__,rank,info);

    if ( (info=HYPRE_ParCSRBiCGSTABDestroy(solver)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJMatrixDestroy(A)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorDestroy(b)) ) throw H_Exception(__LINE__,__FILE__,rank,info);

    T *solution=new T[dist[rank]];
    if ( (info=HYPRE_IJVectorGetValues(x,dist[rank],NULL,solution)) ) throw H_Exception(__LINE__,__FILE__,rank,info);
    if ( (info=HYPRE_IJVectorDestroy(x)) ) throw H_Exception(__LINE__,__FILE__,rank,info);

    MPI_Allgatherv(solution,dist[rank],give_me_MPI_datatype(solution),res,dist,disp,give_me_MPI_datatype(res),internal_comm);
    delete[] dist;
    delete[] disp;
    delete[] solution;

}

/************************************************************************************************/

#endif
