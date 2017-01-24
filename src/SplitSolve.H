#ifndef _SPLITSOLVE
#define _SPLITSOLVE

#include "cublas_v2.h"
#include "CSR.H"
#include "CSC.H"
#include "Types.H"
#include "Blas.H"
#include "LinearSolver.H"

template <class T>
class SplitSolve : public LinearSolver<CPX>{
	
public:

    SplitSolve(TCSR<CPX>* mat,double*,int,int,MPI_Comm,MPI_Comm);

    virtual ~SplitSolve();
    virtual void prepare();
    virtual void prepare(int*,int*,int,int,int*,int);
    virtual void prepare_corner(CPX*,CPX*,int*,int*,int*,int,int,int*,int);
    virtual void solve_equation(CPX* res, CPX* rhs, int no_rhs);    
				
private:

    int glob_rank,glob_size;
    int left_rank,right_rank;
    int mpi_size,mpi_rank;
    int CPU_per_bc;
    int CPU_per_solver;
    int interleave;
    int cpu_color;
    int glob_cond;
    int TB;
    int NBCL,NBCR;
    int n_block;
    int b_size;
    int pivot;
    int neigh_pivot;
    int min_rank;
    int max_rank;
    size_t sigma_size;
    size_t bL_size;
    size_t bR_size;
    size_t y_size;
    size_t z_size;
    MPI_Comm wco_comm,glob_wco_comm;
    TCSR<CPX>* matrix;

    cudaStream_t stream_c,stream_m;
    void *cublas_handle;
    int *IL_first;
    T *HostMem;
    T *M_host;
    T *M1_host;
    T *M_dev;
    T *M1_dev;
    T *M2_dev;
    T *C_dev;
    T *xLU_dev;
    CPX *SigmaB_dev;
    CPX *bL_dev;
    CPX *bR_dev;
    CPX *M3_dev;
    CPX *y_dev;
    CPX *z_dev;
    T *nnz_dev;
    int *edge_i_dev;
    int *index_j_dev;
    int *b_min,*b_max;

    void extract_diag(int*,int*,T*,T*,int*,int*,int,cudaStream_t);
    void extract_not_diag(int*,int*,T*,T*,int*,int*,int,int,cudaStream_t);
    void change_var_type(T*,CPX*,int,cudaStream_t);
    void symmetrize_matrix(T*,int,cudaStream_t);
    void create_blocks(int,int,int*,int*,int*,int);
    void FirstStage();
    void SecondStage();
    void ThirdStage();
    void FourthStage();
    void RecursiveStage();
    void FifthStage(int,int,int,int);
    void SixthStage(int,int,int,int);
    void InitSolve(CPX*,int);
    void FirstSolve(int);
    void SecondSolve(int);
    void ThirdSolve(int);
    void FourthSolve(CPX*,int);
    void write_dev_matrix(char*,double*,int,int);
    void write_dev_matrix(char*,CPX*,int,int);
    void write_matrix(char*,double*,int,int);
    void write_matrix(char*,CPX*,int,int);
    void write_matrix(char*,int*,int,int);
};

/************************************************************************************************/

template <class T>
SplitSolve<T>::SplitSolve(TCSR<CPX>* mat,double *pHostMem,int pCPU_per_bc,int pinterleave,\
			  MPI_Comm solver_comm,MPI_Comm glob_solver_comm)
{
    int cond1,cond2,cond3;
    int cond11,cond12,cond13,cond14;
    int cond21,cond22,cond23;

    TB              = 10;

    matrix          = mat;

    HostMem         = (T*)pHostMem;

    wco_comm        = solver_comm;
    glob_wco_comm   = glob_solver_comm;

    CPU_per_bc      = pCPU_per_bc;
    interleave      = pinterleave;

    MPI_Comm_size(glob_wco_comm,&glob_size);
    MPI_Comm_rank(glob_wco_comm,&glob_rank);

    CPU_per_solver  = glob_size;

    cpu_color      = get_cpu_color(glob_rank,interleave,glob_size,glob_size-CPU_per_bc);

    cond11          = glob_rank>0;
    cond12          = (glob_rank<glob_size-1)&&(CPU_per_bc==2);
    cond13          = (glob_rank<glob_size)&&(CPU_per_bc==1);
    cond14          = CPU_per_bc<=CPU_per_solver;
    cond1           = ((cond11&&(cond12||cond13))||(!interleave))&&cond14;
    cond21          = glob_rank<ceil(CPU_per_solver/2.0);
    cond22          = glob_rank>=(glob_size-floor(CPU_per_solver/2.0));
    cond23          = CPU_per_bc>CPU_per_solver;
    cond2           = (cond21||cond22)&&cond23&&(!interleave);
    cond3           = (!cpu_color)&&(CPU_per_bc>2)&&interleave;
    glob_cond       = cond1||cond2||cond3;

    left_rank       = ceil((double)glob_size/(CPU_per_solver-CPU_per_bc))-1;
    right_rank      = glob_size-ceil((double)glob_size/(CPU_per_solver-CPU_per_bc));

    if(glob_cond){

        MPI_Comm_size(wco_comm,&mpi_size);
	MPI_Comm_rank(wco_comm,&mpi_rank);

	cudaStreamCreate(&stream_c);
	cudaStreamCreate(&stream_m);

	magma_init();
	cublas_init(&cublas_handle);
    }

}

/************************************************************************************************/

template <class T>
SplitSolve<T>::~SplitSolve()
{
    if(glob_cond){

        magma_finalize();
        cublas_finalize(cublas_handle);

	cudaStreamDestroy(stream_c);
        cudaStreamDestroy(stream_m);

	deallocate_data_on_dev(M_dev,matrix->size*b_size*sizeof(T));
	deallocate_data_on_dev(C_dev,b_size*b_size*sizeof(T));
	deallocate_data_on_dev(bL_dev,bL_size);
	deallocate_data_on_dev(bR_dev,bR_size);
	deallocate_data_on_dev(M3_dev,3*b_size*b_size*sizeof(CPX));
	deallocate_data_on_dev(y_dev,y_size);
	deallocate_data_on_dev(z_dev,z_size);

	delete[] IL_first;
	delete[] b_min;
	delete[] b_max;
    }
}

/************************************************************************************************/

template <>
void SplitSolve<double>::change_var_type(double *D1,CPX *D2,int N,cudaStream_t stream)
{

    change_var_type_on_dev(D1,D2,N,stream);

}

/************************************************************************************************/

template <>
void SplitSolve<CPX>::change_var_type(CPX *D1,CPX *D2,int N,cudaStream_t stream)
{

    //cudaMemcpy(D2,D1,N*sizeof(CPX),cudaMemcpyDeviceToDevice);
    cudaMemcpyAsync(D2,D1,N*sizeof(CPX),cudaMemcpyDeviceToDevice,stream);
}

/************************************************************************************************/

template <>
void SplitSolve<double>::symmetrize_matrix(double *M,int N,cudaStream_t stream)
{
    d_symmetrize_matrix(M,N,stream);
}

/************************************************************************************************/

template <>
void SplitSolve<CPX>::symmetrize_matrix(CPX *M,int N,cudaStream_t stream)
{
    z_symmetrize_matrix(M,N,stream);
}

/************************************************************************************************/

template <>
void SplitSolve<double>::extract_diag(int *edge_i,int *index_j,double *nnz,double *D,int *Bmin,\
				      int *Bmax,int index,cudaStream_t stream)
{

    int NR;
    int imin,imax;

    NR   = Bmax[index]-Bmin[index];
    imin = Bmin[index]-matrix->first_row;
    imax = Bmax[index]-matrix->first_row;

    d_init_var_on_dev(D,NR*NR,stream);

    d_extract_diag_on_dev(D,edge_i,index_j,nnz,NR,imin,imax,matrix->first_row,matrix->findx,stream);
}

/************************************************************************************************/

template <>
void SplitSolve<CPX>::extract_diag(int *edge_i,int *index_j,CPX *nnz,CPX *D,int *Bmin,\
				   int *Bmax,int index,cudaStream_t stream)
{

    int NR;
    int imin,imax;

    NR   = Bmax[index]-Bmin[index];
    imin = Bmin[index]-matrix->first_row;
    imax = Bmax[index]-matrix->first_row;

    z_init_var_on_dev(D,NR*NR,stream);

    z_extract_diag_on_dev(D,edge_i,index_j,nnz,NR,imin,imax,matrix->first_row,matrix->findx,\
			  stream);

}

/************************************************************************************************/

template <>
void SplitSolve<double>::extract_not_diag(int *edge_i,int *index_j,double *nnz,double *D,int *Bmin,\
					  int *Bmax,int index,int side,cudaStream_t stream)
{

    int NR,NC;
    int imin,imax;
    int jmin;

    NR   = Bmax[index]-Bmin[index];
    NC   = Bmax[index+side]-Bmin[index+side];
    imin = Bmin[index]-matrix->first_row;
    imax = Bmax[index]-matrix->first_row;
    jmin = Bmin[index+side];

    d_init_var_on_dev(D,NR*NC,stream);

    d_extract_not_diag_on_dev(D,edge_i,index_j,nnz,NR,imin,imax,jmin,side,matrix->first_row,\
			      matrix->findx,stream);
}

/************************************************************************************************/

template <>
void SplitSolve<CPX>::extract_not_diag(int *edge_i,int *index_j,CPX *nnz,CPX *D,int *Bmin,\
				       int *Bmax,int index,int side,cudaStream_t stream)
{

    int NR,NC;
    int imin,imax;
    int jmin;

    NR   = Bmax[index]-Bmin[index];
    NC   = Bmax[index+side]-Bmin[index+side];
    imin = Bmin[index]-matrix->first_row;
    imax = Bmax[index]-matrix->first_row;
    jmin = Bmin[index+side];

    z_init_var_on_dev(D,NR*NC,stream);

    z_extract_not_diag_on_dev(D,edge_i,index_j,nnz,NR,imin,imax,jmin,side,matrix->first_row,\
			      matrix->findx,stream);

}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::FirstStage()
{

    int info;
    int IB;
    int NR,NM,NP;
    int NS1,NS2;
    T ONE     = convert<T,double>(1.0);
    int *ipiv = new int[b_size];

    cublasSetStream((cublasHandle_t)cublas_handle,stream_c);

    if(!(mpi_rank%2)){
        
        NR = b_max[IL_first[mpi_rank]]-b_min[IL_first[mpi_rank]];
	NP = b_max[IL_first[mpi_rank]+1]-b_min[IL_first[mpi_rank]+1];
	
        //extract E-H11
        extract_diag(edge_i_dev,index_j_dev,nnz_dev,M1_dev,b_min,b_max,IL_first[mpi_rank],0);
	//extract -H12
	extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,M_dev,b_min,b_max,IL_first[mpi_rank],1,0);
	//compute (E-H11)^{-1}*(-H12)
	tgesv_dev(NR,NP,M1_dev,NR,ipiv,M_dev,NR,0,&info);
    
	for(IB=(IL_first[mpi_rank]+1);IB<IL_first[mpi_rank+1];IB++){

	    NR  = b_max[IB]-b_min[IB];
	    NM  = b_max[IB-1]-b_min[IB-1];
	    NP  = b_max[IB+1]-b_min[IB+1];
	    NS1 = b_min[IB]-matrix->first_row;
	    NS2 = b_min[IB-1]-matrix->first_row;

	    //extract -Hii-1=M2_dev
	    extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,M2_dev,b_min,b_max,IB,-1,stream_c);
	    //extract E-Hii=M1_dev
	    extract_diag(edge_i_dev,index_j_dev,nnz_dev,M1_dev,b_min,b_max,IB,stream_c);
	    //extract -Hii+1=M_dev[add1]
	    extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,&M_dev[NS1*b_size],b_min,b_max,IB,1,stream_c);

	    //compute E-Hii-Hii-1*gi-1*Hi-1i=M1_dev-M2_dev*M_dev[add2]
	    tgemm_dev(cublas_handle,'N','N',NR,NR,NM,-ONE,M2_dev,NR,&M_dev[NS2*b_size],NM,\
		      ONE,M1_dev,NR);
	    symmetrize_matrix(M1_dev,NR,stream_c);

	    cudaMemcpyAsync(&M_host[NS2*b_size],&M_dev[NS2*b_size],NM*NR*sizeof(T), \
			    cudaMemcpyDeviceToHost,stream_m);

	    cudaStreamSynchronize(stream_c);
	    	    	    
	    //compute (E-Hii-Hii-1*gi-1*Hi-1i)^{-1}*(-Hii+1)=M1_dev^{-1}*M_dev[add1]
	    tgesv_dev(NR,NP,M1_dev,NR,ipiv,&M_dev[NS1*b_size],NR,0,&info);
	}
	
	IB  = IL_first[mpi_rank+1]-1;
	NR  = b_max[IB]-b_min[IB];
	NP  = b_max[IB+1]-b_min[IB+1];
	NS1 = b_min[IB]-matrix->first_row;

	cudaMemcpyAsync(&M_host[NS1*b_size],&M_dev[NS1*b_size],NR*NP*sizeof(T), \
			cudaMemcpyDeviceToHost,stream_m);
	
	//memcpy_to_host(M_host,M_dev,matrix->size*b_size*sizeof(T));

    }else{

        NR  = b_max[IL_first[mpi_rank+1]-1]-b_min[IL_first[mpi_rank+1]-1];
	NM  = b_max[IL_first[mpi_rank+1]-2]-b_min[IL_first[mpi_rank+1]-2];
	NS1 = b_min[IL_first[mpi_rank+1]-1]-matrix->first_row;

        //extract E-HNN
        extract_diag(edge_i_dev,index_j_dev,nnz_dev,M1_dev,b_min,b_max,IL_first[mpi_rank+1]-1,0);
	//extract -HNN-1
	extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,&M_dev[NS1*b_size],b_min,b_max,\
			 IL_first[mpi_rank+1]-1,-1,0);
	//compute (E-HNN)^{-1}*(-HNN-1)
	tgesv_dev(NR,NM,M1_dev,NR,ipiv,&M_dev[NS1*b_size],NR,0,&info);

        for(IB=(IL_first[mpi_rank+1]-2);IB>=IL_first[mpi_rank];IB--){

	    NR  = b_max[IB]-b_min[IB];
	    NM  = b_max[IB-1]-b_min[IB-1];
	    NP  = b_max[IB+1]-b_min[IB+1];
	    NS1 = b_min[IB]-matrix->first_row;
	    NS2 = b_min[IB+1]-matrix->first_row;

	    //extract -Hii+1=M2_dev
	    extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,M2_dev,b_min,b_max,IB,1,stream_c);
	    //extract E-Hii=M1_dev
	    extract_diag(edge_i_dev,index_j_dev,nnz_dev,M1_dev,b_min,b_max,IB,stream_c);
	    //extract -Hii-1=M_dev[add1]
	    extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,&M_dev[NS1*b_size],b_min,\
			     b_max,IB,-1,stream_c);

	    //compute E-Hii-Hii+1*gi+1*Hi+1i=M1_dev-M2_dev*M_dev[add2]
	    tgemm_dev(cublas_handle,'N','N',NR,NR,NP,-ONE,M2_dev,NR,&M_dev[NS2*b_size],NP,\
		      ONE,M1_dev,NR);
	    symmetrize_matrix(M1_dev,NR,stream_c);

	    cudaMemcpyAsync(&M_host[NS2*b_size],&M_dev[NS2*b_size],NP*NR*sizeof(T), \
			    cudaMemcpyDeviceToHost,stream_m);

	    cudaStreamSynchronize(stream_c);

	    //compute (E-Hii-Hii+1*gi+1*Hi+1i)^{-1}*(-Hii-1)=M1_dev^{-1}*M_dev[add1]
	    tgesv_dev(NR,NM,M1_dev,NR,ipiv,&M_dev[NS1*b_size],NR,0,&info);
	}

	IB  = IL_first[mpi_rank];
	NR  = b_max[IB]-b_min[IB];
	NM  = b_max[IB-1]-b_min[IB-1];
	NS1 = b_min[IB]-matrix->first_row;

	cudaMemcpyAsync(&M_host[NS1*b_size],&M_dev[NS1*b_size],NR*NM*sizeof(T), \
			cudaMemcpyDeviceToHost,stream_m);

	//memcpy_to_host(M_host,M_dev,matrix->size*b_size*sizeof(T));
    }

    cudaStreamSynchronize(stream_c);
    cudaStreamSynchronize(stream_m);
    
    cublasSetStream((cublasHandle_t)cublas_handle,0);

    delete[] ipiv;
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::SecondStage()
{

    int info;
    int IB;
    int NR,NM,NP;
    int NS1,NS2;
    int N1,NN;
    T ONE     = convert<T,double>(1.0);
    int *ipiv = new int[b_size];

    if(!(mpi_rank%2)){

        NR  = b_max[IL_first[mpi_rank+1]-1]-b_min[IL_first[mpi_rank+1]-1];
	NP  = b_max[IL_first[mpi_rank+1]]-b_min[IL_first[mpi_rank+1]];
	NS1 = b_min[IL_first[mpi_rank+1]-1]-matrix->first_row;

	//MPI_Send_Full(&M_host[NS1*b_size],NR*NP,mpi_rank+1,0,wco_comm);

	//MPI_Recv_Full(M1_host,NR*NP,mpi_rank+1,1,wco_comm);

	MPI_Sendrecv((double*)&M_host[NS1*b_size],NR*NP*sizeof(T)/sizeof(double),MPI_DOUBLE,
		     mpi_rank+1,0,(double*)M1_host,NR*NP*sizeof(T)/sizeof(double),MPI_DOUBLE,
		     mpi_rank+1,1,wco_comm,MPI_STATUS_IGNORE);

	//memcpy_to_device(M1_host,M_dev,NP*NR*sizeof(T));
	cudaMemcpy(M_dev,M1_host,NP*NR*sizeof(T),cudaMemcpyHostToDevice);

        for(IB=(IL_first[mpi_rank+1]-1);IB>=IL_first[mpi_rank];IB--){

	    NR  = b_max[IB]-b_min[IB];
	    NP  = b_max[IB+1]-b_min[IB+1];
	    NS1 = b_min[IB]-matrix->first_row;
	    NS2 = b_min[IB+1]-matrix->first_row;

	    //extract -Hii+1=M2_dev
	    extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,M2_dev,b_min,b_max,IB,1,0);
	    //extract E-Hii=M1_dev
	    extract_diag(edge_i_dev,index_j_dev,nnz_dev,M1_dev,b_min,b_max,IB,0);
	    //compute E-Hii-Hii+1*gi+1*Hi+1i=M1_dev-M2_dev*M_dev[add2]
	    if(IB==(IL_first[mpi_rank+1]-1)){
	        tgemm_dev(cublas_handle,'N','N',NR,NR,NP,-ONE,M2_dev,NR,M_dev,NP,ONE,M1_dev,NR);
	    }else{
	        tgemm_dev(cublas_handle,'N','N',NR,NR,NP,-ONE,M2_dev,NR,&M_dev[NS2*b_size],NP,\
			  ONE,M1_dev,NR);
	    }
	    symmetrize_matrix(M1_dev,NR,0);
 
	    if(IB>IL_first[mpi_rank]){

	        NM  = b_max[IB-1]-b_min[IB-1];

	        //extract -Hii-1=M_dev[add1]
	        extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,&M_dev[NS1*b_size],b_min,b_max,\
				 IB,-1,0);

		//compute (E-Hii-Hii+1*gi+1*Hi+1i)^{-1}*(-Hii-1)=M1_dev^{-1}*M_dev[add1]
		tgesv_dev(NR,NM,M1_dev,NR,ipiv,&M_dev[NS1*b_size],NR,0,&info);

	    }else{

	        //cudaMemcpyAsync(&M_dev[NS1*b_size],M1_dev,NR*NR*sizeof(T),cudaMemcpyDeviceToDevice,NULL);

	        //compute (E-H11-H12*g2*H21)^{-1}
		//tgetrf_dev(NR,NR,&M_dev[NS1*b_size],NR,ipiv,&info);
		//tgetri_dev(NR,&M_dev[NS1*b_size],NR,ipiv,M1_dev,NR*NR,&info);

	        //compute (E-H11-H12*g2*H21)^{-1}
	        if(mpi_rank){
		    N1 = b_max[IB-1]-b_min[IB-1];
		    extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,&M_dev[NS1*b_size],b_min,b_max,IB,-1,0);
		}else{
		    N1 = NR;
		    init_eye_on_dev(&M_dev[NS1*b_size],NR,0);
		}
		tgesv_dev(NR,N1,M1_dev,NR,ipiv,&M_dev[NS1*b_size],NR,0,&info);

		cudaMemcpyAsync(C_dev,&M_dev[NS1*b_size],NR*N1*sizeof(T),cudaMemcpyDeviceToDevice,NULL);
	    }
	}

    }else{

        NR = b_max[IL_first[mpi_rank]]-b_min[IL_first[mpi_rank]];
	NM = b_max[IL_first[mpi_rank]-1]-b_min[IL_first[mpi_rank]-1];

	//MPI_Recv_Full(M1_host,NR*NM,mpi_rank-1,0,wco_comm);

	//MPI_Send_Full(M_host,NM*NR,mpi_rank-1,1,wco_comm);

	MPI_Sendrecv((double*)M_host,NM*NR*sizeof(T)/sizeof(double),MPI_DOUBLE,mpi_rank-1,1,\
		     (double*)M1_host,NR*NM*sizeof(T)/sizeof(double),MPI_DOUBLE,
		     mpi_rank-1,0,wco_comm,MPI_STATUS_IGNORE);

	//memcpy_to_device(M1_host,M_dev,NM*NR*sizeof(T));
	cudaMemcpy(M_dev,M1_host,NM*NR*sizeof(T),cudaMemcpyHostToDevice);

        for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){
	  
	    NR  = b_max[IB]-b_min[IB];
	    NM  = b_max[IB-1]-b_min[IB-1];
	    NS1 = b_min[IB]-matrix->first_row;
	    NS2 = b_min[IB-1]-matrix->first_row;

	    //extract -Hii-1=M2_dev
	    extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,M2_dev,b_min,b_max,IB,-1,0);
	    //extract E-Hii=M1_dev
	    extract_diag(edge_i_dev,index_j_dev,nnz_dev,M1_dev,b_min,b_max,IB,0);

	    //compute E-Hii-Hii-1*gi-1*Hi-1i=M1_dev-M2_dev*M_dev[add2]
	    if(IB==IL_first[mpi_rank]){
	        tgemm_dev(cublas_handle,'N','N',NR,NR,NM,-ONE,M2_dev,NR,M_dev,NM,ONE,M1_dev,NR);
	    }else{
	        tgemm_dev(cublas_handle,'N','N',NR,NR,NM,-ONE,M2_dev,NR,&M_dev[NS2*b_size],NM,\
			  ONE,M1_dev,NR);
	    }
	    symmetrize_matrix(M1_dev,NR,0);

	    if(IB<(IL_first[mpi_rank+1]-1)){

	        NP = b_max[IB+1]-b_min[IB+1];

	        //extract -Hii+1=M_dev[add1]
	        extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,&M_dev[NS1*b_size],b_min,b_max,\
				 IB,1,0);

		//compute (E-Hii-Hii-1*gi-1*Hi-1i)^{-1}*(-Hii+1)=M1_dev^{-1}*M_dev[add1]
		tgesv_dev(NR,NP,M1_dev,NR,ipiv,&M_dev[NS1*b_size],NR,0,&info);

	    }else{

	        //cudaMemcpyAsync(&M_dev[NS1*b_size],M1_dev,NR*NR*sizeof(T),cudaMemcpyDeviceToDevice,NULL);

	        //compute (E-HNN-HNN-1*gN-1*HN-1N)^{-1}
		//tgetrf_dev(NR,NR,&M_dev[NS1*b_size],NR,ipiv,&info);
	        //tgetri_dev(NR,&M_dev[NS1*b_size],NR,ipiv,M1_dev,NR*NR,&info);

	        //compute (E-HNN-HNN-1*gN-1*HN-1N)^{-1}
	        if(mpi_rank<(mpi_size-1)){
		    NN = b_max[IB+1]-b_min[IB+1];
		    extract_not_diag(edge_i_dev,index_j_dev,nnz_dev,&M_dev[NS1*b_size],b_min,b_max,IB,1,0);
		}else{
		    NN = NR;
		    init_eye_on_dev(&M_dev[NS1*b_size],NR,0);
		}
		tgesv_dev(NR,NN,M1_dev,NR,ipiv,&M_dev[NS1*b_size],NR,0,&info);

		cudaMemcpyAsync(C_dev,&M_dev[NS1*b_size],NR*NN*sizeof(T),cudaMemcpyDeviceToDevice,NULL);
	    }
	}
    }

    deallocate_data_on_dev(edge_i_dev,(matrix->size+1)*sizeof(int));
    deallocate_data_on_dev(index_j_dev,matrix->n_nonzeros*sizeof(int));
    deallocate_data_on_dev(nnz_dev,matrix->n_nonzeros*sizeof(T));

    delete[] ipiv;
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::ThirdStage()
{

    int IB;
    int N1,NN;
    int NM,NR,NP;
    int NS1,NS2;
    T ONE  = convert<T,double>(1.0);
    T ZERO = convert<T,double>(0.0);

    cublasSetStream((cublasHandle_t)cublas_handle,stream_c);

    if(!(mpi_rank%2)){

        if(mpi_rank){
	    N1 = b_max[IL_first[mpi_rank]-1]-b_min[IL_first[mpi_rank]-1];
	}else{
	    N1 = b_max[IL_first[mpi_rank]]-b_min[IL_first[mpi_rank]];
	}

	for(IB=(IL_first[mpi_rank]+1);IB<IL_first[mpi_rank+1];IB++){

	    NR  = b_max[IB]-b_min[IB];
	    NM  = b_max[IB-1]-b_min[IB-1];
	    NS1 = b_min[IB]-matrix->first_row;
	    NS2 = b_min[IB-1]-matrix->first_row;

	    tgemm_dev(cublas_handle,'N','N',NR,N1,NM,-ONE,&M_dev[NS1*b_size],NR,&M_dev[NS2*b_size],\
		      NM,ZERO,M1_dev,NR);

	    cudaMemcpyAsync(&M_dev[NS1*b_size],M1_dev,NR*N1*sizeof(T),cudaMemcpyDeviceToDevice,stream_c);

	    if(IB>(IL_first[mpi_rank]+1)){

	        NR  = b_max[IB-2]-b_min[IB-2];
		NP  = b_max[IB-1]-b_min[IB-1];
		NS1 = b_min[IB-2]-matrix->first_row;

	        cudaMemcpyAsync(M1_host,&M_dev[NS1*b_size],NR*N1*sizeof(T),cudaMemcpyDeviceToHost,stream_m);
		cudaMemcpyAsync(&M_dev[NS1*b_size],&M_host[NS1*b_size],NR*NP*sizeof(T),cudaMemcpyHostToDevice, \
				stream_m);

		cudaStreamSynchronize(stream_m);

		c_tcopy(NR*N1,M1_host,1,&M_host[NS1*b_size],1);
	    }

	    cudaStreamSynchronize(stream_c);
	}

	for(IB=(IL_first[mpi_rank+1]-2);IB<IL_first[mpi_rank+1];IB++){

	    NR  = b_max[IB]-b_min[IB];
	    NP  = b_max[IB+1]-b_min[IB+1];
	    NS1 = b_min[IB]-matrix->first_row;
 
	    cudaMemcpyAsync(M1_host,&M_dev[NS1*b_size],NR*N1*sizeof(T),cudaMemcpyDeviceToHost,stream_m);
	    cudaMemcpyAsync(&M_dev[NS1*b_size],&M_host[NS1*b_size],NR*NP*sizeof(T),cudaMemcpyHostToDevice, \
			    stream_m);

	    cudaStreamSynchronize(stream_m);

	    c_tcopy(NR*N1,M1_host,1,&M_host[NS1*b_size],1);
	}

	cudaStreamSynchronize(stream_c);

	/*
	for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

	    NR  = b_max[IB]-b_min[IB];
	    NP  = b_max[IB+1]-b_min[IB+1];
	    NS1 = b_min[IB]-matrix->first_row;

	    memcpy_to_host(M1_host,&M_dev[NS1*b_size],NR*N1*sizeof(T));
	    memcpy_to_device(&M_host[NS1*b_size],&M_dev[NS1*b_size],NR*NP*sizeof(T));
	    c_tcopy(NR*N1,M1_host,1,&M_host[NS1*b_size],1);
	}
	*/

    }else{

        if(mpi_rank<(mpi_size-1)){
	    NN = b_max[IL_first[mpi_rank+1]]-b_min[IL_first[mpi_rank+1]];
	}else{
	    NN = b_max[IL_first[mpi_rank+1]-1]-b_min[IL_first[mpi_rank+1]-1];
	}

        for(IB=(IL_first[mpi_rank+1]-2);IB>=IL_first[mpi_rank];IB--){

	    NR  = b_max[IB]-b_min[IB];
	    NP  = b_max[IB+1]-b_min[IB+1];
	    NS1 = b_min[IB]-matrix->first_row;
	    NS2 = b_min[IB+1]-matrix->first_row;

	    tgemm_dev(cublas_handle,'N','N',NR,NN,NP,-ONE,&M_dev[NS1*b_size],NR,&M_dev[NS2*b_size],\
		      NP,ZERO,M1_dev,NR);

	    cudaMemcpyAsync(&M_dev[NS1*b_size],M1_dev,NR*NN*sizeof(T),cudaMemcpyDeviceToDevice,stream_c);

	    if(IB<(IL_first[mpi_rank+1]-2)){

	        NR  = b_max[IB+2]-b_min[IB+2];
		NM  = b_max[IB+1]-b_min[IB+1];
		NS1 = b_min[IB+2]-matrix->first_row;

	        cudaMemcpyAsync(M1_host,&M_dev[NS1*b_size],NR*NN*sizeof(T),cudaMemcpyDeviceToHost,stream_m);
		cudaMemcpyAsync(&M_dev[NS1*b_size],&M_host[NS1*b_size],NR*NM*sizeof(T),cudaMemcpyHostToDevice, \
				stream_m);

		cudaStreamSynchronize(stream_m);

		c_tcopy(NR*NN,M1_host,1,&M_host[NS1*b_size],1);
	    }

	    cudaStreamSynchronize(stream_c);
	}

	for(IB=IL_first[mpi_rank]+1;IB>=IL_first[mpi_rank];IB--){

	    NR  = b_max[IB]-b_min[IB];
	    NM  = b_max[IB-1]-b_min[IB-1];
	    NS1 = b_min[IB]-matrix->first_row;

	    cudaMemcpyAsync(M1_host,&M_dev[NS1*b_size],NR*NN*sizeof(T),cudaMemcpyDeviceToHost,stream_m);
	    cudaMemcpyAsync(&M_dev[NS1*b_size],&M_host[NS1*b_size],NR*NM*sizeof(T),cudaMemcpyHostToDevice, \
			    stream_m);

	    cudaStreamSynchronize(stream_m);

	    c_tcopy(NR*NN,M1_host,1,&M_host[NS1*b_size],1);
	}

	cudaStreamSynchronize(stream_c);

	/*
	for(IB=(IL_first[mpi_rank+1]-1);IB>=IL_first[mpi_rank];IB--){

	    NR  = b_max[IB]-b_min[IB];
	    NM  = b_max[IB-1]-b_min[IB-1];
	    NS1 = b_min[IB]-matrix->first_row;

	    memcpy_to_host(M1_host,&M_dev[NS1*b_size],NR*NN*sizeof(T));
	    memcpy_to_device(&M_host[NS1*b_size],&M_dev[NS1*b_size],NR*NM*sizeof(T));
	    c_tcopy(NR*NN,M1_host,1,&M_host[NS1*b_size],1);
	}
	*/
    }

    cublasSetStream((cublasHandle_t)cublas_handle,0);

    cudaStreamSynchronize(stream_c);
    cudaStreamSynchronize(stream_m);
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::FourthStage()
{

    int IB;
    int N1,NN;
    int NR,NP,NM;
    int NS1,NS2;
    T ONE  = convert<T,double>(1.0);
    T ZERO = convert<T,double>(0.0);

    if(!(mpi_rank%2)){

        if(mpi_rank){
	    N1  = b_max[IL_first[mpi_rank]-1]-b_min[IL_first[mpi_rank]-1];
	}else{
	    N1  = b_max[IL_first[mpi_rank]]-b_min[IL_first[mpi_rank]];
	}
	if(mpi_rank<(mpi_size-2)){
	    NN  = b_max[IL_first[mpi_rank+2]]-b_min[IL_first[mpi_rank+2]];
	}else{
	    NN  = b_max[IL_first[mpi_rank+2]-1]-b_min[IL_first[mpi_rank+2]-1];
	}

        NR  = b_max[IL_first[mpi_rank+1]-1]-b_min[IL_first[mpi_rank+1]-1];
	NP  = b_max[IL_first[mpi_rank+1]]-b_min[IL_first[mpi_rank+1]];
	NS1 = b_min[IL_first[mpi_rank+1]-1]-matrix->first_row;

	//MPI_Send_Full(&M_host[NS1*b_size],NR*N1,mpi_rank+1,0,wco_comm);

	//MPI_Recv_Full(M1_host,NP*NN,mpi_rank+1,1,wco_comm);

	MPI_Sendrecv((double*)&M_host[NS1*b_size],NR*N1*sizeof(T)/sizeof(double),MPI_DOUBLE,
		     mpi_rank+1,0,(double*)M1_host,NP*NN*sizeof(T)/sizeof(double),MPI_DOUBLE,
		     mpi_rank+1,1,wco_comm,MPI_STATUS_IGNORE);

	//memcpy_to_device(M1_host,M2_dev,NP*NN*sizeof(T));
	cudaMemcpy(M2_dev,M1_host,NP*NN*sizeof(T),cudaMemcpyHostToDevice);

	for(IB=(IL_first[mpi_rank+1]-1);IB>=IL_first[mpi_rank];IB--){

	    NR  = b_max[IB]-b_min[IB];
	    NP  = b_max[IB+1]-b_min[IB+1];
	    NS1 = b_min[IB]-matrix->first_row;
	    NS2 = b_min[IB+1]-matrix->first_row;

	    if(IB==(IL_first[mpi_rank+1]-1)){
	        tgemm_dev(cublas_handle,'N','N',NR,NN,NP,-ONE,&M_dev[NS1*b_size],NR,M2_dev,\
			  NP,ZERO,M1_dev,NR);
	    }else{
	        tgemm_dev(cublas_handle,'N','N',NR,NN,NP,-ONE,&M_dev[NS1*b_size],NR,&M_dev[NS2*b_size],\
			  NP,ZERO,M1_dev,NR);
	    }

	    cudaMemcpyAsync(&M_dev[NS1*b_size],M1_dev,NR*NN*sizeof(T),cudaMemcpyDeviceToDevice,NULL);
	}

    }else{

        if(mpi_rank>1){
	    N1  = b_max[IL_first[mpi_rank-1]-1]-b_min[IL_first[mpi_rank-1]-1];
	}else{
	    N1  = b_max[IL_first[mpi_rank-1]]-b_min[IL_first[mpi_rank-1]];
	}
	if(mpi_rank<(mpi_size-1)){
	    NN  = b_max[IL_first[mpi_rank+1]]-b_min[IL_first[mpi_rank+1]];
	}else{
	    NN  = b_max[IL_first[mpi_rank+1]-1]-b_min[IL_first[mpi_rank+1]-1];
	}

        NR  = b_max[IL_first[mpi_rank]]-b_min[IL_first[mpi_rank]];
	NM  = b_max[IL_first[mpi_rank]-1]-b_min[IL_first[mpi_rank]-1];

	//MPI_Recv_Full(M1_host,NM*N1,mpi_rank-1,0,wco_comm);

	//MPI_Send_Full(M_host,NR*NN,mpi_rank-1,1,wco_comm);

	MPI_Sendrecv((double*)M_host,NR*NN*sizeof(T)/sizeof(double),MPI_DOUBLE,
		     mpi_rank-1,1,(double*)M1_host,NM*N1*sizeof(T)/sizeof(double),MPI_DOUBLE,
		     mpi_rank-1,0,wco_comm,MPI_STATUS_IGNORE);

	//memcpy_to_device(M1_host,M2_dev,NM*N1*sizeof(T));
	cudaMemcpy(M2_dev,M1_host,NM*N1*sizeof(T),cudaMemcpyHostToDevice);

	for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

	    NR  = b_max[IB]-b_min[IB];
	    NM  = b_max[IB-1]-b_min[IB-1];
	    NS1 = b_min[IB]-matrix->first_row;
	    NS2 = b_min[IB-1]-matrix->first_row;

	    if(IB==IL_first[mpi_rank]){
	        tgemm_dev(cublas_handle,'N','N',NR,N1,NM,-ONE,&M_dev[NS1*b_size],NR,M2_dev,\
			  NM,ZERO,M1_dev,NR);
	    }else{
	        tgemm_dev(cublas_handle,'N','N',NR,N1,NM,-ONE,&M_dev[NS1*b_size],NR,&M_dev[NS2*b_size],\
			  NM,ZERO,M1_dev,NR);
	    }

	    cudaMemcpyAsync(&M_dev[NS1*b_size],M1_dev,NR*N1*sizeof(T),cudaMemcpyDeviceToDevice,NULL);
	}
    }
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::RecursiveStage()
{

    int IS;
    int ind;
    int group_index;
    int N1,NN;
    int factor = 2;

    for(IS=1;IS<log2(mpi_size);IS++){

        min_rank    = mpi_rank-mpi_rank%factor;
	max_rank    = min_rank+factor-1;
	group_index = min_rank/factor;

	N1          = b_max[0]-b_min[0];
	NN          = b_max[n_block-1]-b_min[n_block-1];

	if(!(group_index%2)){

	    pivot       = max_rank;
	    neigh_pivot = pivot+1;

	    if(min_rank){
	        ind = IL_first[min_rank]-1;
		N1  = b_max[ind]-b_min[ind];
	    }

	    if((min_rank+2*factor)<mpi_size){
	        ind = IL_first[min_rank+2*factor];
		NN  = b_max[ind]-b_min[ind];
	    }

	}else{

	    pivot       = min_rank;
	    neigh_pivot = pivot-1;

	    if((min_rank-factor)){
	        ind = IL_first[min_rank-factor]-1;
		N1  = b_max[ind]-b_min[ind];
	    }

	    if(max_rank<(mpi_size-1)){
	        ind = IL_first[max_rank+1];
		NN  = b_max[ind]-b_min[ind];
	    }

	}

	//printf("IS: %i N1: %i NN: %i pivot: %i neigh_pivot: %i group_index %i (%i)\n",IS,N1,NN,pivot,neigh_pivot,group_index,mpi_rank);

	if(IS==1){
	    allocate_data_on_device((void**)&xLU_dev,2*b_size*b_size*sizeof(T));
	}

	FifthStage(group_index,factor,N1,NN);
	SixthStage(group_index,factor,N1,NN);

	if(IS==(log2(mpi_size)-1)){
	    deallocate_data_on_dev(xLU_dev,2*b_size*b_size*sizeof(T));
	}

	factor = 2*factor;
    }

    deallocate_data_on_dev(M1_dev,2*b_size*b_size*sizeof(T));
    deallocate_data_on_dev(M2_dev,b_size*b_size*sizeof(T));
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::FifthStage(int group_index,int group_size,int N1,int NN)
{

    int info;
    int IB,IP;
    int NR,NM,NP;
    int NS1;
    int ind;
    T ONE      = convert<T,double>(1.0);
    T ZERO     = convert<T,double>(0.0);
    //T *M2_host = new T[2*b_size*b_size];
    T *M2_host = &HostMem[b_size*(matrix->size+b_size)];
    int *ipiv  = new int[b_size];
  
    if(!(group_index%2)){

        if(mpi_rank!=pivot){

	    ind = IL_first[neigh_pivot];
	    NP  = b_max[ind]-b_min[ind];

	    //receiving xL                                                                                                                     
	    MPI_Recv_Full(M2_host,NP*(N1+NN),neigh_pivot,200+mpi_rank,wco_comm);
	    //memcpy_to_device(M2_host,xLU_dev,NP*(N1+NN)*sizeof(T));
	    cudaMemcpy(xLU_dev,M2_host,NP*(N1+NN)*sizeof(T),cudaMemcpyHostToDevice);

	}else{

	    ind = IL_first[mpi_rank+1]-1;

	    NR  = b_max[ind]-b_min[ind];
	    NP  = b_max[ind+1]-b_min[ind+1];
	    NS1 = b_min[ind]-matrix->first_row;

	    //M1_host = W^U_N
	    //memcpy_to_host(M1_host,C_dev,NR*NP*sizeof(T));

	    //sending W^U_N
	    //MPI_Send_Full(M1_host,NR*NP,mpi_rank+1,0,wco_comm);
	    //MPI_Send_Full(&M_host[NS1*b_size],NR*NP,mpi_rank+1,0,wco_comm);
	    
	    //receiving W^L_1
	    //MPI_Recv_Full(M1_host,NR*NP,mpi_rank+1,1,wco_comm);

	    //sending W^U_N and receiving W^L_1
	    MPI_Sendrecv((double*)&M_host[NS1*b_size],NR*NP*sizeof(T)/sizeof(double),MPI_DOUBLE,
			 mpi_rank+1,0,(double*)M1_host,NR*NP*sizeof(T)/sizeof(double),MPI_DOUBLE,
			 mpi_rank+1,1,wco_comm,MPI_STATUS_IGNORE);

	    //M1_dev = W^L_1
	    //memcpy_to_device(M1_host,M1_dev,NP*NR*sizeof(T));
	    cudaMemcpy(M1_dev,M1_host,NP*NR*sizeof(T),cudaMemcpyHostToDevice);

	    //M1_host = G^U_{N1}
	    //memcpy_to_host(M1_host,&M_dev[NS1*b_size],NR*N1*sizeof(T));
	    cudaMemcpy(M1_host,&M_dev[NS1*b_size],NR*N1*sizeof(T),cudaMemcpyDeviceToHost);

	    //sending G^U_{N1}
	    MPI_Send_Full(M1_host,NR*N1,mpi_rank+1,2,wco_comm);

	    //receiving G^L_{1N}
	    MPI_Recv_Full(M1_host,NP*NN,mpi_rank+1,3,wco_comm);

	    //M2_dev = G^L_{1N}
	    //memcpy_to_device(M1_host,M2_dev,NP*NN*sizeof(T));
	    cudaMemcpy(M2_dev,M1_host,NP*NN*sizeof(T),cudaMemcpyHostToDevice);

	    //Start Calculation
	    //copy G^U_{N1} to left part of xLU
	    cudaMemcpyAsync(xLU_dev,&M_dev[NS1*b_size],NR*N1*sizeof(T),\
			    cudaMemcpyDeviceToDevice,NULL);
	    //put -WU_N*G^L_{1N} to right part of xLU
	    tgemm_dev(cublas_handle,'N','N',NR,NN,NP,-ONE,C_dev,NR,M2_dev,\
		      NP,ZERO,&xLU_dev[NR*N1],NR);

	    //creating I
	    init_eye_on_dev(M2_dev,NR,0);
	    //adding -W^U_N*W^L_1
	    tgemm_dev(cublas_handle,'N','N',NR,NR,NP,-ONE,C_dev,NR,M1_dev,NP,ONE,M2_dev,NR);
	    //calculating xLU = inv(I-W^U_N*W^L_1)*[G^U_{N1} -WU_N*G^L_{1N}]
	    tgesv_dev(NR,N1+NN,M2_dev,NR,ipiv,xLU_dev,NR,1,&info);
	    //copying xU to M1_dev
	    cudaMemcpyAsync(M1_dev,xLU_dev,NR*(N1+NN)*sizeof(T),cudaMemcpyDeviceToDevice,NULL);

	    //receiving xL
	    MPI_Recv_Full(M2_host,NP*(N1+NN),mpi_rank+1,200+mpi_rank,wco_comm);

	    //Copying xL to device
	    //memcpy_to_device(M2_host,xLU_dev,NP*(N1+NN)*sizeof(T));
	    //memcpy_to_host(M2_host,M1_dev,NR*(N1+NN)*sizeof(T));
	    cudaMemcpy(xLU_dev,M2_host,NP*(N1+NN)*sizeof(T),cudaMemcpyHostToDevice);
	    cudaMemcpy(M2_host,M1_dev,NR*(N1+NN)*sizeof(T),cudaMemcpyDeviceToHost);

	    //sending xU
	    for(IP=neigh_pivot;IP<(neigh_pivot+group_size);IP++){
	        MPI_Send_Full(M2_host,NR*(N1+NN),IP,100+IP,wco_comm);
	    }
	}

	if(mpi_rank>=(min_rank+group_size/2)){

	    ind = IL_first[neigh_pivot];
	    NP  = b_max[ind]-b_min[ind];

	    //swap GU_n1 (GPU) and WU (CPU)
	    for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

	        NR  = b_max[IB]-b_min[IB];
		NS1 = b_min[IB]-matrix->first_row;

		cudaMemcpyAsync(M1_dev,&M_host[NS1*b_size],NR*NP*sizeof(T),\
				cudaMemcpyHostToDevice,NULL);
		cudaMemcpyAsync(&M_host[NS1*b_size],&M_dev[NS1*b_size],NR*N1*sizeof(T),\
				cudaMemcpyDeviceToHost,NULL);
		cudaMemcpyAsync(&M_dev[NS1*b_size],M1_dev,NR*NP*sizeof(T),\
				cudaMemcpyDeviceToDevice,NULL);

		//memcpy_to_host(M1_host,&M_dev[NS1*b_size],NR*N1*sizeof(T));
		//memcpy_to_device(&M_host[NS1*b_size],&M_dev[NS1*b_size],NR*NP*sizeof(T));
		//c_tcopy(NR*N1,M1_host,1,&M_host[NS1*b_size],1);
	    }
	}

    }else{

	if(mpi_rank!=pivot){
	  
	    ind = IL_first[pivot]-1;
	    NM  = b_max[ind]-b_min[ind];

	    //receiving xU                                                                                
	    MPI_Recv_Full(M2_host,NM*(N1+NN),neigh_pivot,100+mpi_rank,wco_comm);
	    //memcpy_to_device(M2_host,xLU_dev,NM*(N1+NN)*sizeof(T));
	    cudaMemcpy(xLU_dev,M2_host,NM*(N1+NN)*sizeof(T),cudaMemcpyHostToDevice);

	}else{

	    ind = IL_first[mpi_rank];

	    NR  = b_max[ind]-b_min[ind];
	    NM  = b_max[ind-1]-b_min[ind-1];
	    NS1 = b_min[ind]-matrix->first_row;

	    //receiving W^U_N
	    //MPI_Recv_Full(M1_host,NR*NM,mpi_rank-1,0,wco_comm);

	    //sending W^L_1 and receiving W^U_N1
	    MPI_Sendrecv((double*)M_host,NR*NM*sizeof(T)/sizeof(double),MPI_DOUBLE,
			 mpi_rank-1,1,(double*)M1_host,NR*NM*sizeof(T)/sizeof(double),MPI_DOUBLE,
			 mpi_rank-1,0,wco_comm,MPI_STATUS_IGNORE);

	    //M1_dev = W^U_N
	    //memcpy_to_device(M1_host,M1_dev,NR*NM*sizeof(T));
	    cudaMemcpy(M1_dev,M1_host,NR*NM*sizeof(T),cudaMemcpyHostToDevice);

	    //M1_host = W^L_1
	    //memcpy_to_host(M1_host,C_dev,NR*NM*sizeof(T));

	    //sending W^L_1
	    //MPI_Send_Full(M1_host,NR*NM,mpi_rank-1,1,wco_comm);
	    //MPI_Send_Full(M_host,NR*NM,mpi_rank-1,1,wco_comm);

	    //receiving G^U_{N1}
	    MPI_Recv_Full(M1_host,NM*N1,mpi_rank-1,2,wco_comm);

	    //M2_dev = G^U_{N1}
	    //memcpy_to_device(M1_host,M2_dev,NM*N1*sizeof(T));
	    cudaMemcpy(M2_dev,M1_host,NM*N1*sizeof(T),cudaMemcpyHostToDevice);

	    //M1_host = G^L_{1N}
	    //memcpy_to_host(M1_host,M_dev,NR*NN*sizeof(T));
	    cudaMemcpy(M1_host,M_dev,NR*NN*sizeof(T),cudaMemcpyDeviceToHost);

	    //sending G^L_{1N}
	    MPI_Send_Full(M1_host,NR*NN,mpi_rank-1,3,wco_comm);

	    //Start Calculation
	    //copy G^L_{1N} to right part of xLU
	    cudaMemcpyAsync(&xLU_dev[NR*N1],M_dev,NR*NN*sizeof(T),\
			    cudaMemcpyDeviceToDevice,NULL);
	    //put -WL_1*G^U_{N1} to left part of xLU
	    tgemm_dev(cublas_handle,'N','N',NR,N1,NM,-ONE,C_dev,NR,M2_dev,\
		      NM,ZERO,xLU_dev,NR);

	    //create I
	    init_eye_on_dev(M2_dev,NR,0);
	    //add -W^L_1*W^U_N
	    tgemm_dev(cublas_handle,'N','N',NR,NR,NM,-ONE,C_dev,NR,M1_dev,NM,ONE,M2_dev,NR);
	    //calc xLU = inv(I-W^L_1*W^U_N)*[-WL_1*G^U_{N1} G^L_{1N}]
	    tgesv_dev(NR,N1+NN,M2_dev,NR,ipiv,xLU_dev,NR,1,&info);

	    //sending xL
	    //memcpy_to_host(M2_host,xLU_dev,NR*(N1+NN)*sizeof(T));
	    cudaMemcpy(M2_host,xLU_dev,NR*(N1+NN)*sizeof(T),cudaMemcpyDeviceToHost);
	    for(IP=neigh_pivot;IP>(neigh_pivot-group_size);IP--){
	        MPI_Send_Full(M2_host,NR*(N1+NN),IP,200+IP,wco_comm);
	    }

	    //receiving xU
	    MPI_Recv_Full(M2_host,NM*(N1+NN),neigh_pivot,100+mpi_rank,wco_comm);
	    //memcpy_to_device(M2_host,xLU_dev,NM*(N1+NN)*sizeof(T));
	    cudaMemcpy(xLU_dev,M2_host,NM*(N1+NN)*sizeof(T),cudaMemcpyHostToDevice);
	}

	if(mpi_rank<(min_rank+group_size/2)){

	    ind = IL_first[pivot]-1;
	    NM  = b_max[ind]-b_min[ind];

	    //swap GL_nN (GPU) and WL (CPU)
	    for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

	        NR  = b_max[IB]-b_min[IB];
		NS1 = b_min[IB]-matrix->first_row;

		cudaMemcpyAsync(M1_dev,&M_host[NS1*b_size],NR*NM*sizeof(T),\
				cudaMemcpyHostToDevice,NULL);
		cudaMemcpyAsync(&M_host[NS1*b_size],&M_dev[NS1*b_size],NR*NN*sizeof(T),\
				cudaMemcpyDeviceToHost,NULL);
		cudaMemcpyAsync(&M_dev[NS1*b_size],M1_dev,NR*NM*sizeof(T),\
				cudaMemcpyDeviceToDevice,NULL);

		//memcpy_to_host(M1_host,&M_dev[NS1*b_size],NR*NN*sizeof(T));
		//memcpy_to_device(&M_host[NS1*b_size],&M_dev[NS1*b_size],NR*NM*sizeof(T));
		//c_tcopy(NR*NN,M1_host,1,&M_host[NS1*b_size],1);
	    }
	}
    }

    delete[] ipiv;
    //delete[] M2_host;
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::SixthStage(int group_index,int group_size,int N1,int NN)
{
    int IB;
    int ind;
    int NR;
    int NM;
    int NLU;
    int NS1,NS2;
    int NSide;
    int shift1,shift2;
    T ONE  = convert<T,double>(1.0);
    T ZERO = convert<T,double>(0.0);

    if(!(group_index%2)){
        ind    = IL_first[neigh_pivot];
	NLU    = b_max[ind]-b_min[ind];
	NSide  = N1;
	shift1 = 0;
	shift2 = N1;
    }else{
        ind    = IL_first[pivot]-1;
	NLU    = b_max[ind]-b_min[ind];
	NSide  = NN;
	shift1 = N1;
	shift2 = 0;
    }

    cublasSetStream((cublasHandle_t)cublas_handle,stream_c);
 
    for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

        NR  = b_max[IB]-b_min[IB];
	NM  = b_max[IB-1]-b_min[IB-1];
	NS1 = b_min[IB]-matrix->first_row;
	NS2 = b_min[IB-1]-matrix->first_row;

	//calc -WLU*XLU
	tgemm_dev(cublas_handle,'N','N',NR,N1+NN,NLU,-ONE,&M_dev[NS1*b_size],NR,xLU_dev,NLU, \
		  ZERO,M1_dev,NR);

	cudaMemcpyAsync(&M_dev[NS1*b_size],&M1_dev[NR*shift2],NR*(N1+NN-NSide)*sizeof(T),\
			cudaMemcpyDeviceToDevice,stream_c);

	//calc Green's Function
	if(IB==IL_first[mpi_rank]){
	    cudaMemcpyAsync(M2_dev,&M1_dev[NR*shift1],NR*NSide*sizeof(T),\
			    cudaMemcpyDeviceToDevice,stream_c);
	}else{
	    cudaMemcpyAsync(M1_host,M2_dev,NM*NSide*sizeof(T),\
			    cudaMemcpyDeviceToHost,stream_m);
	    cudaMemcpyAsync(M2_dev,&M1_dev[NR*shift1],NR*NSide*sizeof(T), \
			    cudaMemcpyDeviceToDevice,stream_c);

	    cudaStreamSynchronize(stream_m);

	    c_taxpy(NM*NSide,ONE,M1_host,1,&M_host[NS2*b_size],1);
	}

	cudaStreamSynchronize(stream_c);

	//calc Green's Function
	if(IB==(IL_first[mpi_rank+1]-1)){
	    //memcpy_to_host(M1_host,&M1_dev[NR*shift1],NR*NSide*sizeof(T));
	    cudaMemcpy(M1_host,&M1_dev[NR*shift1],NR*NSide*sizeof(T),cudaMemcpyDeviceToHost);
	    c_taxpy(NR*NSide,ONE,M1_host,1,&M_host[NS1*b_size],1);
	}
    }

    cublasSetStream((cublasHandle_t)cublas_handle,0);

    cudaStreamSynchronize(stream_m);

    //copy G11 to device
    if((!(mpi_rank%2))&&(mpi_rank==min_rank)){
        NR = b_max[IL_first[mpi_rank]]-b_min[IL_first[mpi_rank]];
	//memcpy_to_device(M_host,C_dev,NR*N1*sizeof(T));
	cudaMemcpyAsync(C_dev,M_host,NR*N1*sizeof(T),cudaMemcpyHostToDevice,NULL);
    }

    //copy GNN to device
    if((mpi_rank%2)&&(mpi_rank==max_rank)){
        ind = IL_first[mpi_rank+1]-1;
        NR  = b_max[ind]-b_min[ind];
	NS1 = b_min[ind]-matrix->first_row;
	//memcpy_to_device(&M_host[NS1*b_size],C_dev,NR*NN*sizeof(T));
	cudaMemcpyAsync(C_dev,&M_host[NS1*b_size],NR*NN*sizeof(T),cudaMemcpyHostToDevice,NULL);
    }
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::prepare()
{
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::prepare(int *Bmin,int *Bmax,int NBlock,int Bsize,int *orb_per_at,int tb)
{

    if(glob_cond){

        double T0,t0,t1,t2,t3,t4,t5;
	int tot_rank;

	T0 = get_time(0.0);

	t0 = get_time(0.0);
	create_blocks(Bsize,NBlock,Bmin,Bmax,orb_per_at,tb);

	allocate_data_on_device((void**)&edge_i_dev,(matrix->size+1)*sizeof(int));
	allocate_data_on_device((void**)&index_j_dev,matrix->n_nonzeros*sizeof(int));
	allocate_data_on_device((void**)&nnz_dev,matrix->n_nonzeros*sizeof(T));

	memcpy_to_device(matrix->edge_i,edge_i_dev,(matrix->size+1)*sizeof(int));
	memcpy_to_device(matrix->index_j,index_j_dev,matrix->n_nonzeros*sizeof(int));
	//memcpy_to_device(matrix->nnz,nnz_dev,matrix->n_nonzeros*sizeof(CPX));
	cudaMemcpy2DAsync(nnz_dev,sizeof(T),matrix->nnz,sizeof(CPX),sizeof(T),matrix->n_nonzeros,cudaMemcpyHostToDevice,NULL);


	M_host  = HostMem;
	M1_host = &HostMem[b_size*matrix->size];

	allocate_data_on_device((void**)&M_dev,b_size*matrix->size*sizeof(T));
	allocate_data_on_device((void**)&M1_dev,2*b_size*b_size*sizeof(T));
	allocate_data_on_device((void**)&M2_dev,b_size*b_size*sizeof(T));
	allocate_data_on_device((void**)&C_dev,b_size*b_size*sizeof(T));
	t0 = get_time(t0);

	t1 = get_time(0.0);
	FirstStage();
	t1 = get_time(t1);
	t2 = get_time(0.0);
	SecondStage();
	t2 = get_time(t2);
	t3 = get_time(0.0);
	ThirdStage();
	t3 = get_time(t3);
	t4 = get_time(0.0);
	FourthStage();
	t4 = get_time(t4);
	t5 = get_time(0.0);
	RecursiveStage();
	t5 = get_time(t5);
	T0 = get_time(T0);

	//printf("Prepare time (%i): %e (s) / %e (s) / %e (s) / %e (s) / %e (s) / %e (s) / %e (s)\n",mpi_rank,t0,t1,t2,t3,t4,t5,T0);
    }
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::create_blocks(int Bsize,int NBlock,int *Bmin,int *Bmax,int *orb_per_at,int tb)
{
    int IB,IR;
    int my_nb;
    int b_size_loc;
    int *loc_part;
    int *glob_part;
    int *loc_frow;
    int *glob_frow;

    IL_first   = new int[mpi_size+1];

    loc_part   = new int[mpi_size];
    glob_part  = new int[mpi_size];

    my_nb      = floor((double)matrix->size/Bsize);
    b_size_loc = ceil((double)matrix->size/my_nb);

    init_var(loc_part,mpi_size);
    loc_part[mpi_rank] = my_nb;

    loc_frow   = new int[mpi_size];
    glob_frow  = new int[mpi_size+1];

    init_var(loc_frow,mpi_size);
    loc_frow[mpi_rank]  = matrix->first_row;
    glob_frow[mpi_size] = matrix->size_tot;

    MPI_Allreduce(&my_nb,&n_block,1,MPI_INT,MPI_SUM,wco_comm);
    MPI_Allreduce(loc_part,glob_part,mpi_size,MPI_INT,MPI_SUM,wco_comm);
    MPI_Allreduce(&b_size_loc,&b_size,1,MPI_INT,MPI_MAX,wco_comm);
    MPI_Allreduce(loc_frow,glob_frow,mpi_size,MPI_INT,MPI_SUM,wco_comm);

    b_min       = new int[max(n_block,NBlock)];
    b_max       = new int[max(n_block,NBlock)];

    IL_first[0] = 0;

    if(n_block>=NBlock){

        for(IR=0;IR<mpi_size;IR++){

	    IL_first[IR+1] = IL_first[IR]+glob_part[IR]; 

	    for(IB=IL_first[IR];IB<IL_first[IR+1];IB++){
	        b_min[IB] = glob_frow[IR]+(IB-IL_first[IR])*b_size;
		b_max[IB] = min(glob_frow[IR]+(IB-IL_first[IR]+1)*b_size,glob_frow[IR+1]);
	    }
	}

    }else{

        b_size  = Bsize;
        n_block = NBlock;
	IR      = 0;

	for(IB=0;IB<NBlock;IB++){

	    b_min[IB] = get_msize(0,Bmin[IB]-1,tb,orb_per_at);
	    b_max[IB] = get_msize(0,Bmax[IB],tb,orb_per_at);

	    if(IR<(mpi_size-1)){
	        if(glob_frow[IR+1]==b_max[IB]){
		    IL_first[IR+1] = IB+1;
		    IR++;
		}
	    }
	}

	IL_first[mpi_size] = NBlock;
    }

    delete[] loc_part;
    delete[] glob_part;
    delete[] loc_frow;
    delete[] glob_frow;
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::prepare_corner(CPX *ML_ref,CPX *MR_ref,int *NBC,int *Bmin,int *Bmax,\
				   int NBlock,int Bsize,int *orb_per_at,int tb)
{

    MPI_Status status;
    int IC;
    int N1,NN;
    int data[4];
    CPX *SigmaB;

    if(!glob_rank){

        MPI_Recv(&N1,1,MPI_INT,left_rank,10,glob_wco_comm,&status);

        SigmaB  = new CPX[N1*N1];
	NBCL    = NBC[0];

	init_var(SigmaB,N1*N1);

	for(IC=0;IC<NBC[0];IC++){
	    c_zcopy(NBC[0],&ML_ref[IC*NBC[0]],1,&SigmaB[IC*N1],1);
	}

	MPI_Send(&NBCL,1,MPI_INT,left_rank,0,glob_wco_comm);
	MPI_Send(SigmaB,N1*N1,MPI_DOUBLE_COMPLEX,left_rank,1,glob_wco_comm);
	MPI_Send(&NBCL,1,MPI_INT,right_rank,2,glob_wco_comm);

	delete[] SigmaB;
    }

    if(glob_rank==left_rank){

        N1     = b_max[0]-b_min[0];
	//SigmaB = new CPX[N1*N1];
	SigmaB = (CPX*)&HostMem[b_size*(matrix->size+b_size)];

	MPI_Send(&N1,1,MPI_INT,0,10,glob_wco_comm);

	MPI_Recv(&NBCL,1,MPI_INT,0,0,glob_wco_comm,&status);
	MPI_Recv(SigmaB,N1*N1,MPI_DOUBLE_COMPLEX,0,1,glob_wco_comm,&status);
	MPI_Recv(&NBCR,1,MPI_INT,glob_size-1,5,glob_wco_comm,&status);

	sigma_size = N1*N1*sizeof(CPX);
	allocate_data_on_device((void**)&SigmaB_dev,sigma_size);

	cudaMemcpy(SigmaB_dev,SigmaB,N1*N1*sizeof(CPX),cudaMemcpyHostToDevice);

	//delete[] SigmaB;
    }

    if(glob_rank==right_rank){

        NN     = b_max[n_block-1]-b_min[n_block-1];
        //SigmaB = new CPX[NN*NN];
	SigmaB = (CPX*)&HostMem[b_size*(matrix->size+b_size)];

	MPI_Send(&NN,1,MPI_INT,glob_size-1,11,glob_wco_comm);

	MPI_Recv(&NBCR,1,MPI_INT,glob_size-1,3,glob_wco_comm,&status);
	MPI_Recv(SigmaB,NN*NN,MPI_DOUBLE_COMPLEX,glob_size-1,4,glob_wco_comm,&status);
	MPI_Recv(&NBCL,1,MPI_INT,0,2,glob_wco_comm,&status);

	sigma_size = NN*NN*sizeof(CPX);
	allocate_data_on_device((void**)&SigmaB_dev,sigma_size);

	cudaMemcpy(SigmaB_dev,SigmaB,NN*NN*sizeof(CPX),cudaMemcpyHostToDevice);

	//delete[] SigmaB;
    }

    if(glob_rank==(glob_size-1)){

        MPI_Recv(&NN,1,MPI_INT,right_rank,11,glob_wco_comm,&status);

        SigmaB = new CPX[NN*NN];
	NBCR   = NBC[1];

	init_var(SigmaB,NN*NN);

	for(IC=0;IC<NBC[1];IC++){
	    c_zcopy(NBC[1],&MR_ref[IC*NBC[1]],1,&SigmaB[(NN-NBC[1]+IC)*NN+NN-NBC[1]],1);
	}

	MPI_Send(&NBCR,1,MPI_INT,right_rank,3,glob_wco_comm);
	MPI_Send(SigmaB,NN*NN,MPI_DOUBLE_COMPLEX,right_rank,4,glob_wco_comm);
	MPI_Send(&NBCR,1,MPI_INT,left_rank,5,glob_wco_comm);

	delete[] SigmaB;
    }
}
/************************************************************************************************/

template <class T>
void SplitSolve<T>::solve_equation(CPX* res, CPX* arg_rhs, int no_rhs)
{   

    if(glob_cond){

        double t0,t1,t2,t3,t4;
       
        t0 = get_time(0.0);
        InitSolve(arg_rhs,no_rhs);
	t0 = get_time(t0);
        t1 = get_time(0.0);
	FirstSolve(no_rhs);
	t1 = get_time(t1);
	t2 = get_time(0.0);
	SecondSolve(no_rhs);
	t2 = get_time(t2);
	t3 = get_time(0.0);
	ThirdSolve(no_rhs);
	t3 = get_time(t3);
	t4 = get_time(0.0);
	FourthSolve(res,no_rhs);
	t4 = get_time(t4);
	
	//printf("Solve time (%i): %e (s) / %e (s) / %e (s) / %e (s) / %e (s)\n",mpi_rank,t0,t1,t2,t3,t4);

    }
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::InitSolve(CPX *arg_rhs,int no_rhs)
{

    MPI_Status status;
    int IC;
    int N1,NN;
    CPX *b;

    N1      = b_max[0]-b_min[0];
    NN      = b_max[n_block-1]-b_min[n_block-1];
 
    bL_size = N1*no_rhs*sizeof(CPX);
    bR_size = NN*no_rhs*sizeof(CPX);

    allocate_data_on_device((void**)&bL_dev,bL_size);
    allocate_data_on_device((void**)&bR_dev,bR_size);

    if(!mpi_rank){

        //b = new CPX[b_size*no_rhs];
        b = (CPX*)&HostMem[b_size*(matrix->size+b_size)];

	init_var(b,b_size*no_rhs);

	for(IC=0;IC<no_rhs;IC++){
	    c_zcopy(NBCL,&arg_rhs[IC*NBCL],1,&b[IC*N1],1);
	}

	cudaMemcpyAsync(bL_dev,b,bL_size,cudaMemcpyHostToDevice,NULL);

	MPI_Send(b,N1*no_rhs,MPI_DOUBLE_COMPLEX,mpi_size-1,0,wco_comm);

	MPI_Recv(b,NN*no_rhs,MPI_DOUBLE_COMPLEX,mpi_size-1,1,wco_comm,&status);

	cudaMemcpy(bR_dev,b,bR_size,cudaMemcpyHostToDevice);

	//delete[] b;

    }

    if(mpi_rank==(mpi_size-1)){

        //b = new CPX[2*b_size*no_rhs];
        b = (CPX*)&HostMem[b_size*(matrix->size+b_size)];

	init_var(b,b_size*no_rhs);

	for(IC=0;IC<no_rhs;IC++){
	    c_zcopy(NBCR,&arg_rhs[IC*NBCR],1,&b[IC*NN+NN-NBCR],1);
	}
	
	cudaMemcpyAsync(bR_dev,b,bR_size,cudaMemcpyHostToDevice,NULL);

	MPI_Recv(&b[b_size*no_rhs],N1*no_rhs,MPI_DOUBLE_COMPLEX,0,0,wco_comm,&status);

	MPI_Send(b,NN*no_rhs,MPI_DOUBLE_COMPLEX,0,1,wco_comm);

	cudaMemcpy(bL_dev,&b[b_size*no_rhs],bL_size,cudaMemcpyHostToDevice);

	//delete[] b;
    }
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::FirstSolve(int no_rhs)
{

    int N1,NN;
    int NS;

    N1     = b_max[0]-b_min[0];
    NN     = b_max[n_block-1]-b_min[n_block-1];

    y_size = matrix->size*no_rhs*sizeof(CPX);
    
    allocate_data_on_device((void**)&M3_dev,3*b_size*b_size*sizeof(CPX));
    allocate_data_on_device((void**)&y_dev,y_size);

    if(!mpi_rank){

        //change type of G11
        change_var_type(C_dev,M3_dev,N1*N1,0);
	//change type of G1N
        change_var_type(M_dev,&M3_dev[b_size*b_size],N1*NN,0);

	//y1=G11*b1
	zgemm_on_dev(cublas_handle,'N','N',N1,no_rhs,N1,CPX(1.0,0.0),M3_dev,N1,bL_dev,N1,\
		     CPX(0.0,0.0),y_dev,N1);
	//y1=y1+G1N*bN
	zgemm_on_dev(cublas_handle,'N','N',N1,no_rhs,NN,CPX(1.0,0.0),&M3_dev[b_size*b_size],\
		     N1,bR_dev,NN,CPX(1.0,0.0),y_dev,N1);

    }

    if(mpi_rank==(mpi_size-1)){

        NS = b_min[n_block-1]-matrix->first_row;

        //change type of GNN
        change_var_type(C_dev,&M3_dev[b_size*b_size],NN*NN,0);
	//change type of GN1
        change_var_type(&M_dev[NS*b_size],M3_dev,NN*N1,0);

	//yN=GNN*bN
	zgemm_on_dev(cublas_handle,'N','N',NN,no_rhs,NN,CPX(1.0,0.0),&M3_dev[b_size*b_size],\
		     NN,bR_dev,NN,CPX(0.0,0.0),&y_dev[NS*no_rhs],NN);
	//yN=yN+GN1*b1
	zgemm_on_dev(cublas_handle,'N','N',NN,no_rhs,N1,CPX(1.0,0.0),M3_dev,NN,bL_dev,N1,\
		     CPX(1.0,0.0),&y_dev[NS*no_rhs],NN);
    }
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::SecondSolve(int no_rhs)
{

    int N1,NN;
    int NS;

    N1     = b_max[0]-b_min[0];
    NN     = b_max[n_block-1]-b_min[n_block-1];

    z_size = 2*b_size*no_rhs*sizeof(CPX);

    allocate_data_on_device((void**)&z_dev,z_size);

    if(!mpi_rank){

        //SigmaL*y1
        zgemm_on_dev(cublas_handle,'N','N',N1,no_rhs,N1,CPX(1.0,0.0),SigmaB_dev,N1,y_dev,N1,\
		     CPX(0.0,0.0),z_dev,N1);

	//I-SigmaL*G11
	z_init_eye_on_dev(&M3_dev[2*b_size*b_size],N1,0);
	zgemm_on_dev(cublas_handle,'N','N',N1,N1,N1,CPX(-1.0,0.0),SigmaB_dev,N1,M3_dev,N1,\
		     CPX(1.0,0.0),&M3_dev[2*b_size*b_size],N1);
	cudaMemcpyAsync(M3_dev,&M3_dev[2*b_size*b_size],N1*N1*sizeof(CPX),\
			cudaMemcpyDeviceToDevice,NULL);

	//-SigmaL*G1N
	zgemm_on_dev(cublas_handle,'N','N',N1,NN,N1,CPX(-1.0,0.0),SigmaB_dev,N1,&M3_dev[b_size*b_size],\
		     N1,CPX(0.0,0.0),&M3_dev[2*b_size*b_size],N1);
	cudaMemcpyAsync(&M3_dev[b_size*b_size],&M3_dev[2*b_size*b_size],N1*NN*sizeof(CPX),\
			cudaMemcpyDeviceToDevice,NULL);

	deallocate_data_on_dev(SigmaB_dev,sigma_size);
    }

    if(mpi_rank==(mpi_size-1)){

        NS = b_min[n_block-1]-matrix->first_row;

	//SigmaR*yN
	zgemm_on_dev(cublas_handle,'N','N',NN,no_rhs,NN,CPX(1.0,0.0),SigmaB_dev,NN,&y_dev[NS*no_rhs],\
		     NN,CPX(0.0,0.0),z_dev,NN);

	//-SigmaR*GN1
	zgemm_on_dev(cublas_handle,'N','N',NN,N1,NN,CPX(-1.0,0.0),SigmaB_dev,NN,M3_dev,NN,\
		     CPX(0.0,0.0),&M3_dev[2*b_size*b_size],NN);
	cudaMemcpyAsync(M3_dev,&M3_dev[2*b_size*b_size],NN*N1*sizeof(CPX),\
			cudaMemcpyDeviceToDevice,NULL);

	//I-SigmaR*GNN
	z_init_eye_on_dev(&M3_dev[2*b_size*b_size],NN,0);
	zgemm_on_dev(cublas_handle,'N','N',NN,NN,NN,CPX(-1.0,0.0),SigmaB_dev,NN,&M3_dev[b_size*b_size],\
		     NN,CPX(1.0,0.0),&M3_dev[2*b_size*b_size],NN);
	cudaMemcpyAsync(&M3_dev[b_size*b_size],&M3_dev[2*b_size*b_size],NN*NN*sizeof(CPX),\
			cudaMemcpyDeviceToDevice,NULL);

	deallocate_data_on_dev(SigmaB_dev,sigma_size);
    }
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::ThirdSolve(int no_rhs)
{

    MPI_Status status;
    int info;
    int IC;
    int N1,NN;
    int NS;
    int *ipiv;
    CPX *T_host;
    CPX *z_host;

    N1     = b_max[0]-b_min[0];
    NN     = b_max[n_block-1]-b_min[n_block-1];

    T_host = (CPX*)&HostMem[b_size*(matrix->size+b_size)];
    z_host = (CPX*)&HostMem[b_size*(matrix->size+b_size)+8*b_size*b_size*sizeof(double)/sizeof(T)];

    if((!mpi_rank)||(mpi_rank==(mpi_size-1))){
        ipiv   = new int[N1+NN];
	//T_host = new CPX[4*b_size*b_size];
	//z_host = new CPX[3*b_size*no_rhs];

	//memcpy_to_host(&T_host[2*b_size*b_size],M3_dev,2*b_size*b_size*sizeof(CPX));
	cudaMemcpy(&T_host[2*b_size*b_size],M3_dev,2*b_size*b_size*sizeof(CPX), 
		   cudaMemcpyDeviceToHost);
	//memcpy_to_host(&z_host[2*b_size*no_rhs],z_dev,b_size*no_rhs*sizeof(CPX));
	cudaMemcpy(&z_host[2*b_size*no_rhs],z_dev,b_size*no_rhs*sizeof(CPX), 
		   cudaMemcpyDeviceToHost);
    }else{
        ipiv   = new int[1];
	//T_host = new CPX[1];
	//z_host = new CPX[1];
    }

    //System to solve: [[A B];[C D]]*[x;y]=[e;f]
    //upper part stored on mpi_rank=0
    //lower part stored on mpi_rank=mpi_size-1
    //x=inv(A-B*inv(D)*C)*(e-B*inv(D)*f)
    //y=inv(D-C*inv(A)*B)*(f-C*inv(A)*e)

    if(!mpi_rank){

        //copy upper left part of matrix
        for(IC=0;IC<NBCL;IC++){
	    c_zcopy(NBCL,&T_host[2*b_size*b_size+IC*N1],1,&T_host[IC*NBCL],1);
	}

	//copy upper right part of matrix
	for(IC=0;IC<NBCR;IC++){
	    c_zcopy(NBCL,&T_host[3*b_size*b_size+(IC+NN-NBCR)*N1],1,&T_host[(NBCL+IC)*NBCL],1);
	}

	//copy upper part of z_dev
	for(IC=0;IC<no_rhs;IC++){
	    c_zcopy(NBCL,&z_host[2*b_size*no_rhs+IC*N1],1,&z_host[IC*NBCL],1);
	}

	cudaMemcpyAsync(M3_dev,T_host,NBCL*(NBCL+NBCR)*sizeof(CPX),cudaMemcpyHostToDevice,NULL);
	cudaMemcpyAsync(z_dev,z_host,NBCL*no_rhs*sizeof(CPX),cudaMemcpyHostToDevice,NULL);

	//factorize A
	zgetrf_nopiv_dev(NBCL,NBCL,M3_dev,NBCL,&info);
	//calc inv(A)*B
	zgetrs_nopiv_dev('N',NBCL,NBCR,M3_dev,NBCL,&M3_dev[NBCL*NBCL],NBCL,&info);
	//calc inv(A)*e
	zgetrs_nopiv_dev('N',NBCL,no_rhs,M3_dev,NBCL,z_dev,NBCL,&info);

	//copy to host
	cudaMemcpy(&T_host[2*b_size*b_size],&M3_dev[NBCL*NBCL],NBCL*NBCR*sizeof(CPX),\
		   cudaMemcpyDeviceToHost);
	cudaMemcpy(&z_host[b_size*no_rhs],z_dev,NBCL*no_rhs*sizeof(CPX),cudaMemcpyDeviceToHost);

	//exchange data (inv(A)*B and inv(A)*e with inv(D)*C and inv(D)*f)
	/*
	MPI_Send(&T_host[2*b_size*b_size],NBCL*NBCR,MPI_DOUBLE_COMPLEX,mpi_size-1,0,wco_comm);
	MPI_Send(&z_host[b_size*no_rhs],NBCL*no_rhs,MPI_DOUBLE_COMPLEX,mpi_size-1,1,wco_comm);

	MPI_Recv(&T_host[3*b_size*b_size],NBCR*NBCL,MPI_DOUBLE_COMPLEX,mpi_size-1,2,\
		 wco_comm,&status);
	MPI_Recv(&z_host[2*b_size*no_rhs],NBCR*no_rhs,MPI_DOUBLE_COMPLEX,mpi_size-1,3,\
		 wco_comm,&status);
	*/

	MPI_Sendrecv(&T_host[2*b_size*b_size],NBCL*NBCR,MPI_DOUBLE_COMPLEX,mpi_size-1,0,\
		     &T_host[3*b_size*b_size],NBCR*NBCL,MPI_DOUBLE_COMPLEX,mpi_size-1,1,
		     wco_comm,MPI_STATUS_IGNORE);

	MPI_Sendrecv(&z_host[b_size*no_rhs],NBCL*no_rhs,MPI_DOUBLE_COMPLEX,mpi_size-1,2,\
		     &z_host[2*b_size*no_rhs],NBCR*no_rhs,MPI_DOUBLE_COMPLEX,mpi_size-1,3,
		     wco_comm,MPI_STATUS_IGNORE);

	//copy again A and B to device
	cudaMemcpyAsync(M3_dev,T_host,NBCL*(NBCL+NBCR)*sizeof(CPX),cudaMemcpyHostToDevice,NULL);
	cudaMemcpyAsync(z_dev,z_host,NBCL*no_rhs*sizeof(CPX),cudaMemcpyHostToDevice,NULL);

	//copy inv(D)*C and inv(D)*f to device
	cudaMemcpyAsync(&M3_dev[2*b_size*b_size],&T_host[3*b_size*b_size],NBCR*NBCL*sizeof(CPX),\
			cudaMemcpyHostToDevice,NULL);
	cudaMemcpyAsync(&z_dev[b_size*no_rhs],&z_host[2*b_size*no_rhs],NBCR*no_rhs*sizeof(CPX),\
			cudaMemcpyHostToDevice,NULL);

	//calc A-B*inv(D)*C
	zgemm_on_dev(cublas_handle,'N','N',NBCL,NBCL,NBCR,CPX(-1.0,0.0),&M3_dev[NBCL*NBCL],NBCL,\
		     &M3_dev[2*b_size*b_size],NBCR,CPX(1.0,0.0),M3_dev,NBCL);
	//calc e-B*inv(D)*f
	zgemm_on_dev(cublas_handle,'N','N',NBCL,no_rhs,NBCR,CPX(-1.0,0.0),&M3_dev[NBCL*NBCL],NBCL,\
		     &z_dev[b_size*no_rhs],NBCR,CPX(1.0,0.0),z_dev,NBCL);

	//calc z_dev=inv(A-B*inv(D)*C)*(e-B*inv(D)*f)
	magma_zgesv_nopiv_gpu(NBCL,no_rhs,(magmaDoubleComplex_ptr)M3_dev,NBCL,\
			      (magmaDoubleComplex_ptr)z_dev,NBCL,&info);

	//copy z_dev to host
	cudaMemcpy(&z_host[b_size*no_rhs],z_dev,NBCL*no_rhs*sizeof(CPX),cudaMemcpyDeviceToHost);
	MPI_Recv(&z_host[2*b_size*no_rhs],NBCR*no_rhs,MPI_DOUBLE_COMPLEX,mpi_size-1,4,\
		 wco_comm,&status);

	for(IC=0;IC<no_rhs;IC++){
	    c_zcopy(NBCL,&z_host[b_size*no_rhs+IC*NBCL],1,&z_host[IC*(NBCL+NBCR)],1);
	    c_zcopy(NBCR,&z_host[2*b_size*no_rhs+IC*NBCR],1,&z_host[IC*(NBCL+NBCR)+NBCL],1);
	}

	cudaMemcpyAsync(z_dev,z_host,(NBCL+NBCR)*no_rhs*sizeof(CPX),cudaMemcpyHostToDevice,NULL);
	MPI_Send(z_host,(NBCL+NBCR)*no_rhs,MPI_DOUBLE_COMPLEX,mpi_size-1,5,wco_comm);

    }

    if(mpi_rank==(mpi_size-1)){


        //copy lower left part of matrix
	for(IC=0;IC<NBCL;IC++){
	    c_zcopy(NBCR,&T_host[2*b_size*b_size+IC*NN+NN-NBCR],1,&T_host[IC*NBCR],1);
	}

	//copy lower right part of matrix
	for(IC=0;IC<NBCR;IC++){
	    c_zcopy(NBCR,&T_host[3*b_size*b_size+(IC+NN-NBCR)*NN+NN-NBCR],1,\
		    &T_host[(NBCL+IC)*NBCR],1);
	}

	//copy lower part of z_dev
	for(IC=0;IC<no_rhs;IC++){
	    c_zcopy(NBCR,&z_host[2*b_size*no_rhs+IC*NN+NN-NBCR],1,&z_host[IC*NBCR],1);
	}

	cudaMemcpyAsync(M3_dev,T_host,NBCR*(NBCL+NBCR)*sizeof(CPX),cudaMemcpyHostToDevice,NULL);
	cudaMemcpyAsync(z_dev,z_host,NBCR*no_rhs*sizeof(CPX),cudaMemcpyHostToDevice,NULL);

	//factorize D
	zgetrf_nopiv_dev(NBCR,NBCR,&M3_dev[NBCR*NBCL],NBCR,&info);
	//calc inv(D)*C
	zgetrs_nopiv_dev('N',NBCR,NBCL,&M3_dev[NBCR*NBCL],NBCR,M3_dev,NBCR,&info);
	//calc inv(D)*f
	zgetrs_nopiv_dev('N',NBCR,no_rhs,&M3_dev[NBCR*NBCL],NBCR,z_dev,NBCR,&info);

	//copy to host
	cudaMemcpy(&T_host[2*b_size*b_size],M3_dev,NBCR*NBCL*sizeof(CPX),cudaMemcpyDeviceToHost);
	cudaMemcpy(&z_host[b_size*no_rhs],z_dev,NBCR*no_rhs*sizeof(CPX),cudaMemcpyDeviceToHost);

	//exchange data (inv(D)*C and inv(D)*f with inv(A)*B and inv(A)*e)
	/*
	MPI_Recv(&T_host[3*b_size*b_size],NBCL*NBCR,MPI_DOUBLE_COMPLEX,0,0,\
		 wco_comm,&status);
	MPI_Recv(&z_host[2*b_size*no_rhs],NBCL*no_rhs,MPI_DOUBLE_COMPLEX,0,1,\
		 wco_comm,&status);

	MPI_Send(&T_host[2*b_size*b_size],NBCR*NBCL,MPI_DOUBLE_COMPLEX,0,2,wco_comm);
	MPI_Send(&z_host[b_size*no_rhs],NBCR*no_rhs,MPI_DOUBLE_COMPLEX,0,3,wco_comm);
	*/

	MPI_Sendrecv(&T_host[2*b_size*b_size],NBCR*NBCL,MPI_DOUBLE_COMPLEX,0,1,\
		     &T_host[3*b_size*b_size],NBCL*NBCR,MPI_DOUBLE_COMPLEX,0,0,\
		     wco_comm,MPI_STATUS_IGNORE);

	MPI_Sendrecv(&z_host[b_size*no_rhs],NBCR*no_rhs,MPI_DOUBLE_COMPLEX,0,3,\
		     &z_host[2*b_size*no_rhs],NBCL*no_rhs,MPI_DOUBLE_COMPLEX,0,2,
		     wco_comm,MPI_STATUS_IGNORE);

	//copy again D and D to device
	cudaMemcpyAsync(M3_dev,T_host,NBCR*(NBCL+NBCR)*sizeof(CPX),cudaMemcpyHostToDevice,NULL);
	cudaMemcpyAsync(z_dev,z_host,NBCR*no_rhs*sizeof(CPX),cudaMemcpyHostToDevice,NULL);

	//copy inv(A)*B and inv(A)*e to device
	cudaMemcpyAsync(&M3_dev[2*b_size*b_size],&T_host[3*b_size*b_size],NBCL*NBCR*sizeof(CPX),\
			cudaMemcpyHostToDevice,NULL);
	cudaMemcpyAsync(&z_dev[b_size*no_rhs],&z_host[2*b_size*no_rhs],NBCL*no_rhs*sizeof(CPX),\
			cudaMemcpyHostToDevice,NULL);

	//calc D-C*inv(A)*B
	zgemm_on_dev(cublas_handle,'N','N',NBCR,NBCR,NBCL,CPX(-1.0,0.0),M3_dev,NBCR,\
		     &M3_dev[2*b_size*b_size],NBCL,CPX(1.0,0.0),&M3_dev[NBCR*NBCL],NBCR);
	//calc f-C*inv(A)*e
	zgemm_on_dev(cublas_handle,'N','N',NBCR,no_rhs,NBCL,CPX(-1.0,0.0),M3_dev,NBCR,\
		     &z_dev[b_size*no_rhs],NBCL,CPX(1.0,0.0),z_dev,NBCR);

	//calc z_dev=inv(D-C*inv(A)*B)*(f-C*inv(A)*e)
	magma_zgesv_nopiv_gpu(NBCR,no_rhs,(magmaDoubleComplex_ptr)&M3_dev[NBCR*NBCL],NBCR,\
			      (magmaDoubleComplex_ptr)z_dev,NBCR,&info);

	//copy z_dev to host
	cudaMemcpy(&z_host[2*b_size*no_rhs],z_dev,NBCR*no_rhs*sizeof(CPX),cudaMemcpyDeviceToHost);
	MPI_Send(&z_host[2*b_size*no_rhs],NBCR*no_rhs,MPI_DOUBLE_COMPLEX,0,4,wco_comm);

	MPI_Recv(z_host,(NBCL+NBCR)*no_rhs,MPI_DOUBLE_COMPLEX,0,5,wco_comm,&status);
	cudaMemcpyAsync(z_dev,z_host,(NBCL+NBCR)*no_rhs*sizeof(CPX),cudaMemcpyHostToDevice,NULL);
    }

    delete[] ipiv;
    //delete[] T_host;
    //delete[] z_host;
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::FourthSolve(CPX *res,int no_rhs)
{
    int IB,IC;
    int N1,NN;
    int NB,NS;
    int cond;
    CPX *Sol;

    N1   = b_max[0]-b_min[0];
    NN   = b_max[n_block-1]-b_min[n_block-1];

    cond = (matrix->size*no_rhs)<(5*b_size*b_size);

    if(cond){
        Sol = (CPX*)&HostMem[b_size*(matrix->size+b_size)];
    }else{
        cudaMallocHost((void **)&Sol,matrix->size*no_rhs*sizeof(CPX));
    }

    if((!mpi_rank)||(mpi_rank==(mpi_size-1))){
        for(IC=0;IC<no_rhs;IC++){
	    zaxpy_on_dev(cublas_handle,NBCL,CPX(1.0,0.0),&z_dev[IC*(NBCL+NBCR)],1,&bL_dev[IC*N1],1);
	    zaxpy_on_dev(cublas_handle,NBCR,CPX(1.0,0.0),&z_dev[IC*(NBCL+NBCR)+NBCL],1,	\
			 &bR_dev[IC*NN+NN-NBCR],1);
	}
    }

    cublasSetStream((cublasHandle_t)cublas_handle,stream_c);

    if(mpi_rank<(mpi_size/2)){

        if(!mpi_rank){

	    if(mpi_size>2){
	        //memcpy_to_host(Sol,bL_dev,N1*no_rhs*sizeof(CPX));
		cudaMemcpy(Sol,bL_dev,N1*no_rhs*sizeof(CPX),cudaMemcpyDeviceToHost);
		//memcpy_to_host(&Sol[N1*no_rhs],bR_dev,NN*no_rhs*sizeof(CPX));
		cudaMemcpy(&Sol[N1*no_rhs],bR_dev,NN*no_rhs*sizeof(CPX),cudaMemcpyDeviceToHost);
	    }
	  
	    for(IB=1;IB<(mpi_size/2);IB++){
	        MPI_Send_Full(Sol,(N1+NN)*no_rhs,IB,10+IB,wco_comm);
	    }

	}else{
	    MPI_Recv_Full(Sol,(N1+NN)*no_rhs,0,10+mpi_rank,wco_comm);
	    memcpy_to_device(Sol,bL_dev,N1*no_rhs*sizeof(CPX));
	    memcpy_to_device(&Sol[N1*no_rhs],bR_dev,NN*no_rhs*sizeof(CPX));
	}

        for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){
	  
	    NS = b_min[IB]-matrix->first_row;
	    NB = b_max[IB]-b_min[IB];

	    //get GnN(IB,:)
	    change_var_type(&M_dev[NS*b_size],M3_dev,NB*NN,stream_c);
	    zgemm_on_dev(cublas_handle,'N','N',NB,no_rhs,NN,CPX(1.0,0.0),M3_dev,NB,bR_dev,NN,\
			 CPX(0.0,0.0),&y_dev[NS*no_rhs],NB);

	    if(IB>(IL_first[mpi_rank])){

	        NS = b_min[IB-1]-matrix->first_row;
		NB = b_max[IB-1]-b_min[IB-1];

	        cudaMemcpyAsync(&M_dev[NS*b_size],&M_host[NS*b_size],NB*N1*sizeof(T),\
			    cudaMemcpyHostToDevice,stream_m);
	    }
	    cudaStreamSynchronize(stream_c);
	    cudaStreamSynchronize(stream_m);
	}

	IB = IL_first[mpi_rank+1]-1;
	NS = b_min[IB]-matrix->first_row;
	NB = b_max[IB]-b_min[IB];

	memcpy_to_device(&M_host[NS*b_size],&M_dev[NS*b_size],NB*N1*sizeof(T));

	//memcpy_to_device(M_host,M_dev,matrix->size*b_size*sizeof(T));

	for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

	    NS = b_min[IB]-matrix->first_row;
	    NB = b_max[IB]-b_min[IB];

	    //get Gn1(IB,:)
	    change_var_type(&M_dev[NS*b_size],M3_dev,NB*N1,stream_c);
	    zgemm_on_dev(cublas_handle,'N','N',NB,no_rhs,N1,CPX(1.0,0.0),M3_dev,NB,bL_dev,N1,\
			 CPX(1.0,0.0),&y_dev[NS*no_rhs],NB);

	    if(IB>(IL_first[mpi_rank])){

	        NS = b_min[IB-1]-matrix->first_row;
		NB = b_max[IB-1]-b_min[IB-1];

	        cudaMemcpyAsync(&Sol[NS*no_rhs],&y_dev[NS*no_rhs],NB*no_rhs*sizeof(CPX),\
				cudaMemcpyDeviceToHost,stream_m);
	    }
	    cudaStreamSynchronize(stream_c);
	    cudaStreamSynchronize(stream_m);
	}

	IB = IL_first[mpi_rank+1]-1;
	NS = b_min[IB]-matrix->first_row;
	NB = b_max[IB]-b_min[IB];

	//memcpy_to_host(&Sol[NS*no_rhs],&y_dev[NS*no_rhs],NB*no_rhs*sizeof(CPX));
	cudaMemcpy(&Sol[NS*no_rhs],&y_dev[NS*no_rhs],NB*no_rhs*sizeof(CPX), 
		   cudaMemcpyDeviceToHost);

	//memcpy_to_host(Sol,y_dev,matrix->size*no_rhs*sizeof(CPX));

	for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

	    NS = b_min[IB]-matrix->first_row;
	    NB = b_max[IB]-b_min[IB];

	    for(IC=0;IC<no_rhs;IC++){
	        c_zcopy(NB,&Sol[NS*no_rhs+IC*NB],1,&res[IC*matrix->size+NS],1);
	    }
	}

    }else{

        if(mpi_rank==(mpi_size-1)){

	    if(mpi_size>2){
	        //memcpy_to_host(Sol,bL_dev,N1*no_rhs*sizeof(CPX));
		cudaMemcpy(Sol,bL_dev,N1*no_rhs*sizeof(CPX),cudaMemcpyDeviceToHost);
		//memcpy_to_host(&Sol[N1*no_rhs],bR_dev,NN*no_rhs*sizeof(CPX));
		cudaMemcpy(&Sol[N1*no_rhs],bR_dev,NN*no_rhs*sizeof(CPX),cudaMemcpyDeviceToHost);
	    }
	  
	    for(IB=mpi_size-2;IB>=(mpi_size/2);IB--){
	        MPI_Send_Full(Sol,(N1+NN)*no_rhs,IB,20+IB,wco_comm);
	    }

	}else{
	    MPI_Recv_Full(Sol,(N1+NN)*no_rhs,mpi_size-1,20+mpi_rank,wco_comm);
	    memcpy_to_device(Sol,bL_dev,N1*no_rhs*sizeof(CPX));
	    memcpy_to_device(&Sol[N1*no_rhs],bR_dev,NN*no_rhs*sizeof(CPX));
	}

        for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

	    NS = b_min[IB]-matrix->first_row;
	    NB = b_max[IB]-b_min[IB];

	    //get Gn1(IB,:)
	    change_var_type(&M_dev[NS*b_size],M3_dev,NB*N1,stream_c);
	    zgemm_on_dev(cublas_handle,'N','N',NB,no_rhs,N1,CPX(1.0,0.0),M3_dev,NB,bL_dev,N1,\
			 CPX(0.0,0.0),&y_dev[NS*no_rhs],NB);
       
	    if(IB>(IL_first[mpi_rank])){

	        NS = b_min[IB-1]-matrix->first_row;
		NB = b_max[IB-1]-b_min[IB-1];

	        cudaMemcpyAsync(&M_dev[NS*b_size],&M_host[NS*b_size],NB*NN*sizeof(T),\
			    cudaMemcpyHostToDevice,stream_m);
	    }
	    cudaStreamSynchronize(stream_c);
	    cudaStreamSynchronize(stream_m);
	}

	IB = IL_first[mpi_rank+1]-1;
	NS = b_min[IB]-matrix->first_row;
	NB = b_max[IB]-b_min[IB];

	memcpy_to_device(&M_host[NS*b_size],&M_dev[NS*b_size],NB*NN*sizeof(T));

	//memcpy_to_device(M_host,M_dev,matrix->size*b_size*sizeof(T));

	for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

	    NS = b_min[IB]-matrix->first_row;
	    NB = b_max[IB]-b_min[IB];

	    //get GnN(IB,:)
	    change_var_type(&M_dev[NS*b_size],M3_dev,NB*NN,stream_c);
	    zgemm_on_dev(cublas_handle,'N','N',NB,no_rhs,NN,CPX(1.0,0.0),M3_dev,NB,bR_dev,NN,\
			 CPX(1.0,0.0),&y_dev[NS*no_rhs],NB);

	    if(IB>(IL_first[mpi_rank])){

	        NS = b_min[IB-1]-matrix->first_row;
		NB = b_max[IB-1]-b_min[IB-1];

	        cudaMemcpyAsync(&Sol[NS*no_rhs],&y_dev[NS*no_rhs],NB*no_rhs*sizeof(CPX),\
				cudaMemcpyDeviceToHost,stream_m);
	    }
	    cudaStreamSynchronize(stream_c);
	    cudaStreamSynchronize(stream_m);
	}

	IB = IL_first[mpi_rank+1]-1;
	NS = b_min[IB]-matrix->first_row;
	NB = b_max[IB]-b_min[IB];

	//memcpy_to_host(&Sol[NS*no_rhs],&y_dev[NS*no_rhs],NB*no_rhs*sizeof(CPX));
	cudaMemcpy(&Sol[NS*no_rhs],&y_dev[NS*no_rhs],NB*no_rhs*sizeof(CPX), 
		   cudaMemcpyDeviceToHost);

	//memcpy_to_host(Sol,y_dev,matrix->size*no_rhs*sizeof(CPX));

	for(IB=IL_first[mpi_rank];IB<IL_first[mpi_rank+1];IB++){

	    NS = b_min[IB]-matrix->first_row;
	    NB = b_max[IB]-b_min[IB];

	    for(IC=0;IC<no_rhs;IC++){
	        c_zcopy(NB,&Sol[NS*no_rhs+IC*NB],1,&res[IC*matrix->size+NS],1);
	    }
	}
    }

    cublasSetStream((cublasHandle_t)cublas_handle,0);

    cudaStreamSynchronize(stream_c);
    cudaStreamSynchronize(stream_m);

    if(!cond){
        cudaFreeHost(Sol);
    }
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::write_dev_matrix(char *filename,CPX *Mat_dev,int NR,int NC)
{

    ofstream myfile;
    int IR,IC;
    CPX *Mat_host = new CPX[NR*NC];

    cudaMemcpy(Mat_host,Mat_dev,NR*NC*sizeof(CPX),cudaMemcpyDeviceToHost);

    myfile.open(filename);
    for(IR=0;IR<NR;IR++){
        for(IC=0;IC<NC;IC++){
	  myfile<<real(Mat_host[IR+IC*NR])<<" "<<imag(Mat_host[IR+IC*NR])<<" ";
	}
	myfile<<endl;
    }
    myfile.close();

    delete[] Mat_host;
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::write_dev_matrix(char *filename,double *Mat_dev,int NR,int NC)
{

    ofstream myfile;
    int IR,IC;
    double *Mat_host = new double[NR*NC];

    cudaMemcpy(Mat_host,Mat_dev,NR*NC*sizeof(double),cudaMemcpyDeviceToHost);

    myfile.open(filename);
    for(IR=0;IR<NR;IR++){
        for(IC=0;IC<NC;IC++){
	    myfile<<Mat_host[IR+IC*NR]<<" ";
	}
	myfile<<endl;
    }
    myfile.close();

    delete[] Mat_host;
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::write_matrix(char *filename,CPX *Mat_host,int NR,int NC)
{

    ofstream myfile;
    int IR,IC;

    myfile.open(filename);
    for(IR=0;IR<NR;IR++){
        for(IC=0;IC<NC;IC++){
	  myfile<<real(Mat_host[IR+IC*NR])<<" "<<imag(Mat_host[IR+IC*NR])<<" ";
	}
	myfile<<endl;
    }
    myfile.close();
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::write_matrix(char *filename,double *Mat_host,int NR,int NC)
{

    ofstream myfile;
    int IR,IC;

    myfile.open(filename);
    for(IR=0;IR<NR;IR++){
        for(IC=0;IC<NC;IC++){
	    myfile<<Mat_host[IR+IC*NR]<<" ";
	}
	myfile<<endl;
    }
    myfile.close();
}

/************************************************************************************************/

template <class T>
void SplitSolve<T>::write_matrix(char *filename,int *Mat_host,int NR,int NC)
{

    ofstream myfile;
    int IR,IC;

    myfile.open(filename);
    for(IR=0;IR<NR;IR++){
        for(IC=0;IC<NC;IC++){
	    myfile<<Mat_host[IR+IC*NR]<<" ";
	}
	myfile<<endl;
    }
    myfile.close();
}

/************************************************************************************************/

#endif

