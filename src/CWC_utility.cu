#include <stdio.h>
#include "Types.H"
#include "cublas_v2.h"
#include "cusparse_v2.h"
#include "cuda.h"

#ifndef max_stream
#define max_stream 16
#endif

#ifndef BLOCK_DIM
#define BLOCK_DIM 16
#endif

static volatile size_t c_memory = 0;

extern "C"
void set_gpu(int dev,char *gpu_string){
     struct cudaDeviceProp dprop;
     cudaSetDevice(dev);
     cudaGetDeviceProperties(&dprop, dev);
     strcpy(gpu_string,dprop.name);	
}

extern "C"
void cublas_init(void **handle){
     cublasCreate((cublasHandle_t*)handle);
}

extern "C"
void cublas_finalize(void *handle){
     cublasDestroy((cublasHandle_t)handle);
}

extern "C"
void cusparse_init(void **handle){
     cusparseCreate((cusparseHandle_t*)handle);
}

extern "C"
void cusparse_finalize(void *handle){
     cusparseDestroy((cusparseHandle_t)handle);
}

extern "C"
size_t allocate_data_on_device(void **data,size_t size_data){

     cudaError_t mem_error;  

     mem_error = cudaMalloc(data,size_data);

     if(mem_error!=cudaSuccess){
         printf("CPU wants to allocate %e MBytes on the device, but already %e MBytes are in use\n",size_data/1e6,c_memory/1e6);
	 exit(0);
     }else{
         c_memory = c_memory+size_data;
     }

     return c_memory;
}

extern "C"
void deallocate_data_on_device(void *data){
     cudaFree(data);
}

extern "C"
size_t deallocate_data_on_dev(void *data,size_t size_data){

     cudaFree(data);

     c_memory = c_memory-size_data;

     return c_memory;
}

extern "C"
void copy_data_to_device(void *host_data,void *device_data,int N,int M,size_t size_element){
     cublasSetMatrixAsync(N,M,size_element,host_data,N,device_data,N,NULL);
}

extern "C"
void memcpy_to_device(void *host_data,void *device_data,size_t size_element){
     cudaMemcpyAsync(device_data,host_data,size_element,cudaMemcpyHostToDevice,NULL);
}

extern "C"
void copy_data_to_host(void *host_data,void *device_data,int N,int M,size_t size_element){
     cublasGetMatrixAsync(N,M,size_element,device_data,N,host_data,N,NULL);
}

extern "C"
void memcpy_to_host(void *host_data,void *device_data,size_t size_element){
     cudaMemcpyAsync(host_data,device_data,size_element,cudaMemcpyDeviceToHost,NULL);
}

extern "C"
void dgemm_on_dev(void *handle,char transa,char transb,int m,int n,int k,double alpha,\
      		  double *A,int lda,double *B,int ldb,double beta,double *C, int ldc){
     
     cublasOperation_t opA,opB;

     if(transa=='N'){
	opA = CUBLAS_OP_N;
     }
     if(transa=='C'){
	opA = CUBLAS_OP_C;
     }
     if(transa=='T'){
	opA = CUBLAS_OP_T;
     }

     if(transb=='N'){
	opB = CUBLAS_OP_N;
     }
     if(transb=='C'){
	opB = CUBLAS_OP_C;
     }
     if(transb=='T'){
	opB = CUBLAS_OP_T;
     }

     cublasDgemm((cublasHandle_t)handle,opA,opB,m,n,k,&alpha,A,lda,B,ldb,&beta,C,ldc);
}

extern "C"
void zgemm_on_dev(void *handle,char transa,char transb,int m,int n,int k,CPX alpha,\
                  CPX *A,int lda,CPX *B,int ldb,CPX beta,CPX *C, int ldc){
  
     cublasOperation_t opA,opB;

     if(transa=='N'){
	opA = CUBLAS_OP_N;
     }
     if(transa=='C'){
	opA = CUBLAS_OP_C;
     }
     if(transa=='T'){
	opA = CUBLAS_OP_T;
     }

     if(transb=='N'){
	opB = CUBLAS_OP_N;
     }
     if(transb=='C'){
	opB = CUBLAS_OP_C;
     }
     if(transb=='T'){
	opB = CUBLAS_OP_T;
     }

     cublasZgemm((cublasHandle_t)handle,opA,opB,m,n,k,(cuDoubleComplex*)&alpha,\
                 (cuDoubleComplex*)A,lda,(cuDoubleComplex*)B,ldb,(cuDoubleComplex*)&beta,\
		 (cuDoubleComplex*)C,ldc);
}

extern "C"
void zaxpy_on_dev(void *handle,int n,CPX alpha,CPX *x,int incx,CPX *y,int incy){
  
    cublasZaxpy((cublasHandle_t)handle,n,(cuDoubleComplex*)&alpha,(cuDoubleComplex*)x,\
                incx,(cuDoubleComplex*)y,incy);
}

__global__ void d_init_variable_on_dev(double *var,int N){

     int idx = blockIdx.x*blockDim.x + threadIdx.x;

     if(idx<N){
	var[idx] = 0.0;
     }	   

     __syncthreads();
}

extern "C"
void d_init_var_on_dev(double *var,int N,cudaStream_t stream){

    uint i_N = N + (BLOCK_DIM-(N%BLOCK_DIM));

    d_init_variable_on_dev<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>(var,N);
}

__global__ void d_init_eye_on_device(double *var,int N){

     int idx = blockIdx.x*blockDim.x + threadIdx.x;

     if(idx<N*N){
	var[idx] = 0.0;
	if(!(idx%(N+1))){
	    var[idx] = 1.0;
	}	   
     }

     __syncthreads();
}

extern "C"
void d_init_eye_on_dev(double *var,int N,cudaStream_t stream){

    uint i_N = N*N + (BLOCK_DIM-((N*N)%BLOCK_DIM));

    d_init_eye_on_device<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>(var,N);
}

__global__ void z_init_variable_on_dev(cuDoubleComplex *var,int N){

     int idx = blockIdx.x*blockDim.x + threadIdx.x;

     if(idx<N){
	var[idx].x = 0.0;
	var[idx].y = 0.0;
     }	   

     __syncthreads();
}

extern "C"
void z_init_var_on_dev(CPX *var,int N,cudaStream_t stream){

    uint i_N = N + (BLOCK_DIM-(N%BLOCK_DIM));

    z_init_variable_on_dev<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>((cuDoubleComplex*)var,N);
}

__global__ void z_init_eye_on_device(cuDoubleComplex *var,int N){

     int idx = blockIdx.x*blockDim.x + threadIdx.x;

     if(idx<N*N){
	var[idx].x = 0.0;
	var[idx].y = 0.0;
	if(!(idx%(N+1))){
	    var[idx].x = 1.0;
	}	   
     }

     __syncthreads();
}

extern "C"
void z_init_eye_on_dev(CPX *var,int N,cudaStream_t stream){

    uint i_N = N*N + (BLOCK_DIM-((N*N)%BLOCK_DIM));

    z_init_eye_on_device<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>((cuDoubleComplex*)var,N);
}

__global__ void correct_diag_on_device(cuDoubleComplex *var,int N){

     int idx = blockIdx.x*blockDim.x + threadIdx.x;

     if((idx<N*N)&&(!(idx%(N+1)))){
         var[idx].y = 0.0;	   
     }

     __syncthreads();
}

extern "C"
void correct_diag_on_dev(CPX *var,int N,cudaStream_t stream){

    uint i_N = N*N + (BLOCK_DIM-((N*N)%BLOCK_DIM));

    correct_diag_on_device<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>((cuDoubleComplex*)var,N);
}

__global__ void change_variable_type_on_dev(double *var1,cuDoubleComplex *var2,int N){

     int idx = blockIdx.x*blockDim.x + threadIdx.x;

     if(idx<N){
	var2[idx].x = var1[idx];
	var2[idx].y = 0.0;
     }	   

     __syncthreads();
}

extern "C"
void change_var_type_on_dev(double *var1,CPX *var2,int N,cudaStream_t stream){

    uint i_N = N + (BLOCK_DIM-(N%BLOCK_DIM));

    change_variable_type_on_dev<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>(var1,(cuDoubleComplex*)var2,N);
}

__global__ void change_sign_imaginary_part_on_dev(cuDoubleComplex *var,int N){

     int idx = blockIdx.x*blockDim.x + threadIdx.x;

     if(idx<N){
	var[idx].y = -var[idx].y;
     }	   

     __syncthreads();
}

extern "C"
void change_sign_imag_on_dev(CPX *var,int N){

    uint i_N = N + (BLOCK_DIM-(N%BLOCK_DIM));

    change_sign_imaginary_part_on_dev<<< i_N/BLOCK_DIM, BLOCK_DIM >>>((cuDoubleComplex*)var,N);
}

__global__ void d_extract_diag(double *D,int *edge_i,int *index_j,cuDoubleComplex *nnz,\
	   int NR,int imin,int imax,int shift,int findx){

     int j;
     int ind_j;	   
     int idx = blockIdx.x*blockDim.x + threadIdx.x;

     if(idx<NR){
	  for(j=edge_i[idx+imin]-findx;j<edge_i[idx+imin+1]-findx;j++){
	      ind_j = index_j[j]-findx-shift-imin;
	      if((ind_j>=0)&&(ind_j<NR)){
	          D[idx+ind_j*NR] = nnz[j].x;
	      }
	  }
     }	   

     __syncthreads();
}

extern "C"
void d_extract_diag_on_dev(double *D,int *edge_i,int *index_j,CPX *nnz,int NR,\
     int imin,int imax,int shift,int findx,cudaStream_t stream){

    uint i_N = NR + (BLOCK_DIM-(NR%BLOCK_DIM));

    d_extract_diag<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>(D,edge_i,index_j,(cuDoubleComplex*)nnz,NR,imin,imax,shift,findx);
}

__global__ void d_extract_not_diag(double *D,int *edge_i,int *index_j,cuDoubleComplex *nnz,\
	   int NR,int imin,int imax,int jmin,int side,int shift,int findx){

     int j;
     int ind_j;	   
     int limit = 0;
     int idx   = blockIdx.x*blockDim.x + threadIdx.x;

     if(side==-1){
         limit = -(imin+shift-jmin-1);
     }

     if(idx<NR){
	  for(j=edge_i[idx+imin]-findx;j<edge_i[idx+imin+1]-findx;j++){
	      ind_j = index_j[j]-findx-jmin;
	      if(side*ind_j>=limit){
	          D[idx+ind_j*NR] = nnz[j].x;
	      }
	  }
     }	   

     __syncthreads();
}

extern "C"
void d_extract_not_diag_on_dev(double *D,int *edge_i,int *index_j,CPX *nnz,int NR,\
     int imin,int imax,int jmin,int side,int shift,int findx,cudaStream_t stream){

    uint i_N = NR + (BLOCK_DIM-(NR%BLOCK_DIM));

    d_extract_not_diag<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>(D,edge_i,index_j,(cuDoubleComplex*)nnz,NR,imin,imax,jmin,side,shift,findx);
}

__global__ void z_extract_diag(cuDoubleComplex *D,int *edge_i,int *index_j,cuDoubleComplex *nnz,\
	   int NR,int imin,int imax,int shift,int findx){

     int j;
     int ind_j;	   
     int idx = blockIdx.x*blockDim.x + threadIdx.x;

     if(idx<NR){
	  for(j=edge_i[idx+imin]-findx;j<edge_i[idx+imin+1]-findx;j++){
	      ind_j = index_j[j]-findx-shift-imin;
	      if((ind_j>=0)&&(ind_j<NR)){
	          D[idx+ind_j*NR].x = nnz[j].x;
		  D[idx+ind_j*NR].y = nnz[j].y;
	      }
	  }
     }	   

     __syncthreads();
}

extern "C"
void z_extract_diag_on_dev(CPX *D,int *edge_i,int *index_j,CPX *nnz,int NR,\
     int imin,int imax,int shift,int findx,cudaStream_t stream){

    uint i_N = NR + (BLOCK_DIM-(NR%BLOCK_DIM));

    z_extract_diag<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>((cuDoubleComplex*)D,edge_i,index_j,(cuDoubleComplex*)nnz,NR,imin,imax,shift,findx);
}

__global__ void z_extract_not_diag(cuDoubleComplex *D,int *edge_i,int *index_j,cuDoubleComplex *nnz,\
	   int NR,int imin,int imax,int jmin,int side,int shift,int findx){

     int j;
     int ind_j;	   
     int limit = 0;
     int idx   = blockIdx.x*blockDim.x + threadIdx.x;

     if(side==-1){
         limit = -(imin+shift-jmin-1);
     }

     if(idx<NR){
	  for(j=edge_i[idx+imin]-findx;j<edge_i[idx+imin+1]-findx;j++){
	      ind_j = index_j[j]-findx-jmin;
	      if(side*ind_j>=limit){
	          D[idx+ind_j*NR].x = nnz[j].x;
		  D[idx+ind_j*NR].y = nnz[j].y;
	      }
	  }
     }	   

     __syncthreads();
}

extern "C"
void z_extract_not_diag_on_dev(CPX* *D,int *edge_i,int *index_j,CPX *nnz,int NR,\
     int imin,int imax,int jmin,int side,int shift,int findx,cudaStream_t stream){

    uint i_N = NR + (BLOCK_DIM-(NR%BLOCK_DIM));

    z_extract_not_diag<<< i_N/BLOCK_DIM, BLOCK_DIM, 0, stream >>>((cuDoubleComplex*)D,edge_i,index_j,(cuDoubleComplex*)nnz,NR,imin,imax,jmin,side,shift,findx);
}

extern "C"
void d_copy_csr_to_device(int size,int n_nonzeros,int *hedge_i,int *hindex_j,double *hnnz,\
          		  int *dedge_i,int *dindex_j,double *dnnz){
    
    cudaMemcpyAsync(dedge_i,hedge_i,(size+1)*sizeof(int),cudaMemcpyHostToDevice,NULL);
    cudaMemcpyAsync(dindex_j,hindex_j,n_nonzeros*sizeof(int),cudaMemcpyHostToDevice,NULL);
    cudaMemcpyAsync(dnnz,hnnz,n_nonzeros*sizeof(double),cudaMemcpyHostToDevice,NULL);
}

extern "C"
void z_copy_csr_to_device(int size,int n_nonzeros,int *hedge_i,int *hindex_j,CPX *hnnz,\
          		  int *dedge_i,int *dindex_j,CPX *dnnz){
    
    cudaMemcpyAsync(dedge_i,hedge_i,(size+1)*sizeof(int),cudaMemcpyHostToDevice,NULL);
    cudaMemcpyAsync(dindex_j,hindex_j,n_nonzeros*sizeof(int),cudaMemcpyHostToDevice,NULL);
    cudaMemcpyAsync(dnnz,hnnz,n_nonzeros*sizeof(CPX),cudaMemcpyHostToDevice,NULL);
}

extern "C"
void d_csr_mult_f(void *handle,int m,int n,int k,int n_nonzeros,int *Aedge_i,int *Aindex_j,\
                  double *Annz,double alpha,double *B,double beta,double *C){

    cusparseMatDescr_t descra;

    cusparseCreateMatDescr(&descra);
    cusparseSetMatType(descra,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descra,CUSPARSE_INDEX_BASE_ONE);

    cusparseDcsrmm((cusparseHandle_t)handle,CUSPARSE_OPERATION_NON_TRANSPOSE,m,n,k,n_nonzeros,\
                   &alpha,descra,Annz,Aedge_i,Aindex_j,B,k,&beta,C,m);

    cusparseDestroyMatDescr(descra);
}

extern "C"
void z_csr_mult_f(void *handle,int m,int n,int k,int n_nonzeros,int *Aedge_i,int *Aindex_j,\
                  CPX *Annz,CPX alpha,CPX *B,CPX beta,CPX *C){

    cusparseMatDescr_t descra;

    cusparseCreateMatDescr(&descra);
    cusparseSetMatType(descra,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descra,CUSPARSE_INDEX_BASE_ONE);

    cusparseZcsrmm((cusparseHandle_t)handle,CUSPARSE_OPERATION_NON_TRANSPOSE,m,n,k,n_nonzeros,\
                   (cuDoubleComplex*)&alpha,descra,(cuDoubleComplex*)Annz,Aedge_i,Aindex_j,\
		   (cuDoubleComplex*)B,k,(cuDoubleComplex*)&beta,(cuDoubleComplex*)C,m);

    cusparseDestroyMatDescr(descra);
}

extern "C"
void z_csr_mult_fo(void *handle,int m,int n,int k,int n_nonzeros,int *Aedge_i,int *Aindex_j,\
                   CPX *Annz,CPX alpha,CPX *B,CPX beta,CPX *C){

    cusparseMatDescr_t descra;

    cusparseCreateMatDescr(&descra);
    cusparseSetMatType(descra,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descra,CUSPARSE_INDEX_BASE_ZERO);

    cusparseZcsrmm((cusparseHandle_t)handle,CUSPARSE_OPERATION_NON_TRANSPOSE,m,n,k,n_nonzeros,\
                   (cuDoubleComplex*)&alpha,descra,(cuDoubleComplex*)Annz,Aedge_i,Aindex_j,\
		   (cuDoubleComplex*)B,k,(cuDoubleComplex*)&beta,(cuDoubleComplex*)C,m);

    cusparseDestroyMatDescr(descra);
}

// This kernel is optimized to ensure all global reads and writes are coalesced,
// and to avoid bank conflicts in shared memory.  This kernel is up to 11x faster
// than the naive kernel below.  Note that the shared memory array is sized to 
// (BLOCK_DIM+1)*BLOCK_DIM.  This pads each row of the 2D block in shared memory 
// so that bank conflicts do not occur when threads address the array column-wise.
__global__ void d_transpose(double *odata, double *idata, int width, int height)
{
	__shared__ double block[BLOCK_DIM][BLOCK_DIM+1];
	
	// read the matrix tile into shared memory
	unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
	unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
	if((xIndex < width) && (yIndex < height))
	{
		unsigned int index_in = yIndex * width + xIndex;
		block[threadIdx.y][threadIdx.x] = idata[index_in];
	}

	__syncthreads();

	// write the transposed matrix tile to global memory
	xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
	yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;
	if((xIndex < height) && (yIndex < width))
	{
		unsigned int index_out = yIndex * height + xIndex;
		odata[index_out] = block[threadIdx.x][threadIdx.y];
	}
}

extern "C"
void d_transpose_matrix(double *odata,double *idata,int size_x,int size_y){

    uint i_size_x = size_x + (BLOCK_DIM-(size_x%BLOCK_DIM));
    uint i_size_y = size_y + (BLOCK_DIM-(size_y%BLOCK_DIM));

    dim3 grid(i_size_x / BLOCK_DIM, i_size_y / BLOCK_DIM, 1);
    dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);

    d_transpose<<< grid, threads >>>(odata, idata, size_x, size_y);
}

// This kernel is optimized to ensure all global reads and writes are coalesced,
// and to avoid bank conflicts in shared memory.  This kernel is up to 11x faster
// than the naive kernel below.  Note that the shared memory array is sized to 
// (BLOCK_DIM+1)*BLOCK_DIM.  This pads each row of the 2D block in shared memory 
// so that bank conflicts do not occur when threads address the array column-wise.
__global__ void z_transpose(cuDoubleComplex *odata, cuDoubleComplex *idata, int width, int height)
{
	__shared__ cuDoubleComplex block[BLOCK_DIM][BLOCK_DIM+1];
	
	// read the matrix tile into shared memory
	unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
	unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;
	if((xIndex < width) && (yIndex < height))
	{
		unsigned int index_in = yIndex * width + xIndex;
		block[threadIdx.y][threadIdx.x] = idata[index_in];
	}

	__syncthreads();

	// write the transposed matrix tile to global memory
	xIndex = blockIdx.y * BLOCK_DIM + threadIdx.x;
	yIndex = blockIdx.x * BLOCK_DIM + threadIdx.y;
	if((xIndex < height) && (yIndex < width))
	{
		unsigned int index_out = yIndex * height + xIndex;
		odata[index_out].x = block[threadIdx.x][threadIdx.y].x;
		odata[index_out].y = -block[threadIdx.x][threadIdx.y].y;
	}
}

extern "C"
void z_transpose_matrix(CPX *odata,CPX *idata,int size_x,int size_y){

    uint i_size_x = size_x + (BLOCK_DIM-(size_x%BLOCK_DIM));
    uint i_size_y = size_y + (BLOCK_DIM-(size_y%BLOCK_DIM));

    dim3 grid(i_size_x / BLOCK_DIM, i_size_y / BLOCK_DIM, 1);
    dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);

    z_transpose<<< grid, threads >>>((cuDoubleComplex*)odata,(cuDoubleComplex*)idata,size_x,size_y);
}

__global__ void d_symmetrize(double *matrix, int N)
{

        unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
        unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;

        if((xIndex < N) && (yIndex < N) && (yIndex>=xIndex)){
            unsigned int index_1  = yIndex * N + xIndex;
            unsigned int index_2  = xIndex * N + yIndex;
            double val_1    = matrix[index_1];
            double val_2    = matrix[index_2];

            matrix[index_1] = (val_1+val_2)/2.0;
            matrix[index_2] = (val_1+val_2)/2.0;
        }

        __syncthreads();
}

extern "C"
void d_symmetrize_matrix(double *matrix,int N,cudaStream_t stream){

    uint i_size = N + (BLOCK_DIM-(N%BLOCK_DIM));

    dim3 grid(i_size / BLOCK_DIM, i_size / BLOCK_DIM, 1);
    dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);

    d_symmetrize<<< grid, threads, 0, stream >>>(matrix, N);
}

__global__ void z_symmetrize(cuDoubleComplex *matrix, int N)
{

	unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
	unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;

	if((xIndex < N) && (yIndex < N) && (yIndex>=xIndex)){
	    unsigned int index_1  = yIndex * N + xIndex;
	    unsigned int index_2  = xIndex * N + yIndex;
	    cuDoubleComplex val_1 = matrix[index_1];
	    cuDoubleComplex val_2 = matrix[index_2];
	    
	    matrix[index_1].x     = (val_1.x+val_2.x)/2.0;
	    matrix[index_1].y     = (val_1.y-val_2.y)/2.0;
	    matrix[index_2].x     = (val_1.x+val_2.x)/2.0;
	    matrix[index_2].y     = -(val_1.y-val_2.y)/2.0;
	}

	__syncthreads();
}

extern "C"
void z_symmetrize_matrix(CPX *matrix,int N,cudaStream_t stream){

    uint i_size = N + (BLOCK_DIM-(N%BLOCK_DIM));

    dim3 grid(i_size / BLOCK_DIM, i_size / BLOCK_DIM, 1);
    dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);

    z_symmetrize<<< grid, threads, 0, stream >>>((cuDoubleComplex*)matrix, N);
}

__global__ void z_symmetrize_2(cuDoubleComplex *matrix, int N)
{

	unsigned int xIndex = blockIdx.x * BLOCK_DIM + threadIdx.x;
	unsigned int yIndex = blockIdx.y * BLOCK_DIM + threadIdx.y;

	if((xIndex < N) && (yIndex < N) && (yIndex>=xIndex)){
	    unsigned int index_1  = yIndex * N + xIndex;
	    unsigned int index_2  = xIndex * N + yIndex;
	    cuDoubleComplex val_1 = matrix[index_1];

	    if(yIndex==xIndex){
	        matrix[index_1].x = 0.0;
		matrix[index_1].y = val_1.y;
	    }else{
	        matrix[index_2].x = -val_1.x;
		matrix[index_2].y = val_1.y;
	    }
	}

	__syncthreads();
}

extern "C"
void z_symmetrize_matrix_2(CPX *matrix,int N,cudaStream_t stream){

    uint i_size = N + (BLOCK_DIM-(N%BLOCK_DIM));

    dim3 grid(i_size / BLOCK_DIM, i_size / BLOCK_DIM, 1);
    dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);

    z_symmetrize_2<<< grid, threads, 0, stream >>>((cuDoubleComplex*)matrix, N);
}
