#include "Utilities.H"
#include "Blas.H"

bool sortx(const XYZPOS& a, const XYZPOS& b)
{
    return a.x < b.x;
}

/************************************************************************************************/

bool sorty(const XYZPOS& a, const XYZPOS& b)
{
    return a.y < b.y;
}

/************************************************************************************************/

bool sortz(const XYZPOS& a, const XYZPOS& b)
{
    return a.z < b.z;
}

/************************************************************************************************/

bool sortn(const CONNEC& a, const CONNEC& b)
{
    return a.neigh < b.neigh;
}

/************************************************************************************************/

bool sorti(const IJPOS& a, const IJPOS& b)
{
    return a.i < b.i;
}

/************************************************************************************************/

bool sortj(const IJPOS& a, const IJPOS& b)
{
    return a.j < b.j;
}

/************************************************************************************************/


/************************************************************************************************/

void dreshape(int N, int inc, double *v, int *index)
{
    int i;
    double *v_copy = new double[N*inc];

    c_dcopy(N*inc,v,1,v_copy,1);

    for(i=0;i<N;i++){
        c_dcopy(inc,&v_copy[inc*index[i]],1,&v[inc*i],1);
    }

    delete [] v_copy;
}

/************************************************************************************************/

void icopy(int N, int* x, int* x_copy)
{
    int i;

    for(i=0;i<N;i++){
        x_copy[i]=x[i];
    }
    
}

/************************************************************************************************/

int Round(double x)
{
    if(x<0){
        return (int)(x-0.5);
    }else{
        return (int)(x+0.5);
    }
}

/************************************************************************************************/

void init_CSR(CSR **M, int size_i, int N, int nnz)
{
    
    for(int i=0;i<size_i;i++){
        M[i] = new CSR(N,nnz);
    }
}

/************************************************************************************************/

void del_CSR(CSR **M, int size_i)
{
    for(int i=0;i<size_i;i++){
        delete M[i];
    }

    delete M;
}

/************************************************************************************************/

double min_vec(double *vec,int N,int inc)
{
    double min_val = 1.0e10;
    
    for(int i=0;i<N;i++){
        if(vec[i*inc]<min_val){
            min_val = vec[i*inc];
        }
    }
    return min_val;
}

/************************************************************************************************/

double max_sign_abs_vec(double *vec,int N,int inc)
{
    double max_val = 0.0;
    double result;
    
    for(int i=0;i<N;i++){
        if(fabs(vec[i*inc])>max_val){
            max_val = abs(vec[i*inc]);
            result  = vec[i*inc];
        }
    }
    return result;
}

/************************************************************************************************/

void change_sign(double *vec,int N,int inc,int sign)
{
    
    for(int i=0;i<N;i++){
        vec[i*inc]=sign*vec[i*inc];
    }
}

/************************************************************************************************/

void sort_vec(double *vec,int N)
{
    int ii,jj,pos;
    int* idx_order   = new int[N];
    double* vec_new  = new double[N];
    
    for(ii = 0; ii < N; ii++) {
      pos = 0;
      for(jj = 0; jj < N; jj++) {
         if(ii != jj) {
             if((vec[ii]>vec[jj])
                || (vec[ii] == vec[jj] && ii > jj)) {
                 pos++;
             }
         }
      }
      idx_order[ii] = pos;
    }
    
    for(ii = 0; ii < N; ii++) {
        vec_new[idx_order[ii]] = vec[ii];
    }
    
    c_dcopy(N,vec_new,1,vec,1);

    delete[] idx_order;
    delete[] vec_new;

}

/************************************************************************************************/

void sort_abs_imag(CPX *vec,int *index,int N)
{

    int* idx_order   = new int[N];
    int  pos,ii,jj;
    CPX* vec_new     = new CPX[N];
    
    for(ii = 0; ii < N; ii++) {
      pos = 0;
      for(jj = 0; jj < N; jj++) {
         if(ii != jj) {
             if((abs(imag(vec[index[ii]]))>abs(imag(vec[index[jj]])))
                || (abs(imag(vec[index[ii]])) == abs(imag(vec[index[jj]])) && ii > jj)) {
                 pos++;
             }
         }
      }
      idx_order[ii] = pos;
    }
    
    for(ii = 0; ii < N; ii++) {
        vec_new[idx_order[ii]] = vec[index[ii]];
    }
    
    c_zcopy(N,vec_new,1,vec,1);

    delete[] idx_order;
    delete[] vec_new;
}

/************************************************************************************************/

void check_mpi(int size,int rank,int CPU_per_vg_point,int CPU_per_vd_point,int CPU_per_kz_point,\
	       int Nky,int Nkz,int CPU_ppoint)
{   

    int kz_mpi_size;

    if(CPU_per_vg_point==-1) CPU_per_vg_point = size;
    if(CPU_per_vd_point==-1) CPU_per_vd_point = CPU_per_vg_point;
    if(CPU_per_kz_point==-1) CPU_per_kz_point = CPU_per_vd_point;

    if(size%CPU_per_vg_point){
        if(!rank){
	    printf("The number of CPUs per vg point must be a divider of the total number of CPUs\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	exit(0);
    }

    if(CPU_per_vg_point%CPU_per_vd_point){
        if(!rank){
	    printf("The number of CPUs per vd point must be a divider of the number of CPUs per vg point\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	exit(0);
    }

    if(CPU_per_vd_point%CPU_per_kz_point){
        if(!rank){
	    printf("The number of CPUs per kz point must be a divider of the number of CPUs per vd point\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	exit(0);
    }

    kz_mpi_size = CPU_per_vd_point/CPU_per_kz_point;

    if((Nkz*Nky)%kz_mpi_size){
        if(!rank){
	    printf("The number of kz points that are simultaneously treated must be a divider of Nky*Nkz\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	exit(0);
    }

    if(CPU_per_kz_point%CPU_ppoint){
        if(!rank){
	    printf("The number of CPUs per energy point must be a divider of the number of CPUs per kz point\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	exit(0);
    }
}

/************************************************************************************************/

double get_time(double t0)
{
    timeval tim;
    gettimeofday(&tim, NULL);
    
    return (double)(tim.tv_sec+(tim.tv_usec/1000000.0))-t0;
}

/************************************************************************************************/

double randn()
{
    double r1,r2;

    r1 = (double)rand()/RAND_MAX;
    r2 = (double)rand()/RAND_MAX;

    return sqrt(-2.0*log(r1))*cos(2*PI*r2);
}

/************************************************************************************************/

int get_number_of_blocks(int NLayer,int layer_per_slab,int robust_numerics,int *neighbor_layer)
{

    int IL,inc;

    inc = 1;
    for(IL=0;IL<NLayer;IL++){
        if(neighbor_layer[IL]>inc){
            inc = neighbor_layer[IL];
        }
    }
    
    if((NLayer%inc)||robust_numerics){
        inc = layer_per_slab;
    }
    
    return NLayer/inc;

}

/************************************************************************************************/

void domain_decomposition(int NLayer,int mpi_size,int NCS,int factor,int *IL_start,int *IL_stop,\
			  int spec_decomp)
{
    int IP,IC,NL_CPU,NL_local,max_rank,proc_size;

    proc_size = mpi_size/NCS;
    NL_CPU    = ceil((double)NLayer/proc_size);
    max_rank  = NLayer-(NL_CPU-1)*proc_size;

    NL_local    = NL_CPU;
    for(IP=0;IP<proc_size;IP++){
	if(IP>=max_rank) NL_local = NL_CPU-1;
	for(IC=0;IC<NCS;IC++){
	    IL_start[IP*NCS+IC] = IP*(NL_CPU-1)+min(IP,max_rank);
	    IL_stop[IP*NCS+IC]  = IL_start[IP*NCS+IC]+NL_local;
	}
    }

    if(!spec_decomp){
        for(IP=1;IP<mpi_size;IP++){
	    if(IL_start[IP]%2){
	        IL_start[IP]  = IL_start[IP]-1;
		IL_stop[IP-1] = IL_stop[IP-1]-1;
	    }
	    if(IL_start[IP]%4){
	        IL_start[IP]  = IL_start[IP]+2;
		IL_stop[IP-1] = IL_stop[IP-1]+2;
	    }
	}
	/*
	if(mpi_size>=4){
	    IL_stop[0]           = IL_stop[0]-4;
	    IL_start[1]          = IL_start[1]-4;
	    IL_stop[mpi_size-2]  = IL_stop[mpi_size-2]+4;
	    IL_start[mpi_size-1] = IL_start[mpi_size-1]+4;
	}
	*/
    }

    for(IP=0;IP<mpi_size;IP++){
        IL_start[IP] = IL_start[IP]*factor;
	IL_stop[IP]  = IL_stop[IP]*factor;
    }

}

/************************************************************************************************/

int get_msize(int start,int stop,int tb,int* n_of_el)
{ 
    /*
    msize = 0;
    for(IM=start;IM<=stop;IM++){
        msize = msize+factor*n_of_el[IM];
    }
    */

    return (tb/10)*(n_of_el[stop+1]-n_of_el[start]);
}
/************************************************************************************************/

int get_max_orb(int* n_of_el,int NA)
{ 
    int IA,max_orb;

    max_orb = 0;

    for(IA=0;IA<NA;IA++){

        if((n_of_el[IA+1]-n_of_el[IA])>max_orb){
	    max_orb = (n_of_el[IA+1]-n_of_el[IA]);
	}

    }

    return max_orb;
}
/************************************************************************************************/

int sum_vec(int N,int *vec,int *active)
{
    int i,sum;

    sum = 0;
    for(i=0;i<N;i++){
        if(!active[i]){
	    sum = sum+vec[i];
	}
    }

    return sum;
}

/************************************************************************************************/

int min_active_vec(int N,int *vec,int *active)
{
    int i,min_vec,index;

    min_vec = INF;
    for(i=0;i<N;i++){
        if((!active[i])&&(vec[i]<min_vec)){
	    min_vec = vec[i];
	    index   = i;
	}
    }
    
    return index;
}

/************************************************************************************************/

int max_active_vec(int N,int *vec,int *active)
{
    int i,max_vec,index;

    max_vec = 0;
    for(i=0;i<N;i++){
        if((!active[i])&&(vec[i]>max_vec)){
	    max_vec = vec[i];
	    index   = i;
	}
    }
    
    return index;
}

/************************************************************************************************/
