#include "RGF.H"

/************************************************************************************************/

///Constructor of the class RGF

RGF::RGF(TCSR<CPX>* mat,int pNPCS,MPI_Comm solver_comm)
{
    matrix = mat;
    findx  = 0;

    TB     = 10;

    NPCS   = pNPCS;
}

/************************************************************************************************/

///Destructor of the class RGF

RGF::~RGF()
{
}

/************************************************************************************************/

/*!
\brief Function that solves (E-H-Sigma)*GR=I using the RGF algorithm and returns the diagonal blocks and the first block column of GR

GR                 (output) Diagonal blocks of the retarded Green's Function

GRN1               (output) First block column of the retarded Green's Function

SigRl              (input) Retarded boundary self-energy of the left contact

Gammal             (input/output) Broadening function of the left contact

SigRr              (input) Retarded boundary self-energy of the left contact

Gammar             (input/output) Broadening function of the right contact

Bmin               (input) index of the first atom in each block of the Hamiltonian matrix H

Bmax               (input) index of the last atom in each block of the Hamiltonian matrix H

NBlock             (input) number of blocks

Bsize              (input) size of the largest block of H

NSlab              (input) number of unit cells in the device structure

layer_per_slab     (input) number of atomic layers per device unit cell

NBF                (input) index of the first block belonging to the quantum mechanical domain

NBL                (input) index of the last block belonging to the quantum mechanical region

atom_per_column    (input) vector containing the number of atoms with the same (x,y), but different z coordinates

orb_per_at         (input) number of orbitals per atom

tb                 (input) tight-binding parameter

*/

void RGF::solve_equation(CPX *GR,CPX *GRN1,CPX *SigRl,CPX *Gammal,CPX *SigRr,CPX *Gammar,
                         int *Bmin,int *Bmax,int NBlock,int Bsize,int NSlab,int layer_per_slab,\
			 int NBF,int NBL,int *atom_per_column,int *orb_per_at,int tb)
{
    int IB,N0,N1,NS,NT;
    int *start_index    = new int[NBlock+1];
    int *ipiv           = new int[Bsize];
    CPX *GRN1_act       = new CPX[Bsize*Bsize];
    CPX *M              = new CPX[Bsize*Bsize];
    CPX *TrM            = new CPX[Bsize*Bsize];
    CPX *inv_gR         = new CPX[Bsize*Bsize];
    CPX *SigR           = new CPX[Bsize*Bsize];
    CPX *T01_gR11       = new CPX[Bsize*Bsize];
    CPX *GRr_act        = new CPX[Bsize*Bsize];
    CPX *GRl_act        = new CPX[Bsize*Bsize];
    CPX *Zl             = new CPX[Bsize*Bsize];
    CPX *gR             = new CPX[NBlock*Bsize*Bsize];
    TCSR<CPX> *TRow     = new TCSR<CPX>(Bsize,Bsize*Bsize,findx);
    TCSC<CPX,int> *TCol = new TCSC<CPX,int>(Bsize,Bsize*Bsize,findx);

    //set all the Green's Functions to 0
    init_green(gR,GRr_act,start_index,NBlock,Bmin,Bmax,Bsize,tb,orb_per_at);

    if(NBF>0){
        calc_gR(&gR[start_index[0]],inv_gR,SigRl,Bmin,Bmax,NBlock,ipiv,tb,orb_per_at,0);

	//calculate gR_{i+1}=inv(E-H_{i+1}-T_{i+1,i}*gR_{i}*T_{i,i+1}) for 1<i<NBF
	for(IB=1;IB<NBF;IB++){
	    calc_SigR(SigR,T01_gR11,&gR[start_index[IB-1]],TrM,TRow,TCol,Bmin,Bmax,tb,\
		      orb_per_at,IB,"left");
	    calc_gR(&gR[start_index[IB]],inv_gR,SigR,Bmin,Bmax,NBlock,ipiv,tb,orb_per_at,IB);
	}
	
	//calc new boundary self-energy SigRl and broadening function Gammal for block IB=NBF
	calc_SigR(SigRl,T01_gR11,&gR[start_index[NBF-1]],TrM,TRow,TCol,Bmin,Bmax,tb,\
		  orb_per_at,NBF,"left");
	get_Gamma(Gammal,SigRl,get_msize(Bmin[NBF],Bmax[NBF],tb,orb_per_at));

    }
    
    calc_gR(&gR[start_index[NBlock-1]],inv_gR,SigRr,Bmin,Bmax,NBlock,ipiv,tb,orb_per_at,\
	    NBlock-1);

    //calculate gR_{i-1}=inv(E-H_{i-1}-T_{i-1,i}*gR_{i}*T_{i,i-1}) for NBlock-1>i>NBF
    for(IB=NBlock-2;IB>NBF;IB--){
        calc_SigR(SigR,T01_gR11,&gR[start_index[IB+1]],TrM,TRow,TCol,Bmin,Bmax,tb,\
		  orb_per_at,IB+1,"right");
        calc_gR(&gR[start_index[IB]],inv_gR,SigR,Bmin,Bmax,NBlock,ipiv,tb,orb_per_at,IB);
    }
       
    N0 = get_msize(Bmin[NBF],Bmax[NBF],tb,orb_per_at);
    
    //calculate GR_{NBF}=inv(E-H_{NBF}-T_{NBF,NBF+1}*gR_{NBF+1}*T_{NBF+1,NBF}-SigRl)
    calc_SigR(SigR,T01_gR11,&gR[start_index[NBF+1]],TrM,TRow,TCol,Bmin,Bmax,tb,orb_per_at,\
	      NBF+1,"right");
    c_zaxpy(N0*N0,CPX(1.0,0.0),SigRl,1,SigR,1);
    calc_gR(GRr_act,inv_gR,SigR,Bmin,Bmax,NBlock,ipiv,tb,orb_per_at,NBF);
    c_zcopy(N0*N0,GRr_act,1,GRN1_act,1);
    c_zcopy(N0*N0,GRr_act,1,GRl_act,1);

    NS = get_msize(0,Bmin[NBF]-1,tb,orb_per_at);
    NT = get_msize(0,Bmax[NBlock-1],tb,orb_per_at);
    
    c_zgemm('N','N',N0,N0,N0,CPX(1.0,0.0),GRN1_act,N0,Gammal,N0,CPX(0.0,0.0),M,N0);
    c_zgemm('N','C',N0,N0,N0,CPX(1.0,0.0),M,N0,GRN1_act,N0,CPX(0.0,0.0),Zl,N0);
    c_zcopy(N0,GRr_act,N0+1,&GR[NS],1);
    c_zcopy(N0,Zl,N0+1,&GR[NT+NS],1);
    
    //calculate GR_{i}=gR_{i}+gR_{i}*T_{i,i-1}*GR_{i-1}*T_{i-1,i}*gR_{i}
    for(IB=NBF+1;IB<NBlock;IB++){
        N1 = get_msize(Bmin[IB],Bmax[IB],tb,orb_per_at);
	NS = get_msize(0,Bmin[IB]-1,tb,orb_per_at);
        calc_GR(GRr_act,GRN1_act,&gR[start_index[IB]],M,TrM,TRow,TCol,Bmin,Bmax,\
		tb,orb_per_at,IB,NBL,"right");
	c_zcopy(N1,GRr_act,N1+1,&GR[NS],1);
	if(IB<=NBL){
	    //contribution from left and right contact
	    c_zgemm('N','N',N1,N0,N0,CPX(1.0,0.0),GRN1_act,N1,Gammal,N0,CPX(0.0,0.0),M,N1);
	    c_zgemm('N','C',N1,N1,N0,CPX(1.0,0.0),M,N1,GRN1_act,N1,CPX(0.0,0.0),Zl,N1);
	    c_zcopy(N1,Zl,N1+1,&GR[NT+NS],1);
	}else{
	    //region right to the quantum mechanical region
	    init_var(&GR[NT+NS],N1);
	}
    }

    //region left to the quantum mechanical region
    for(IB=NBF-1;IB>=0;IB--){
      
        N1 = get_msize(Bmin[IB],Bmax[IB],tb,orb_per_at);
	NS = get_msize(0,Bmin[IB]-1,tb,orb_per_at);
	
        calc_GR(GRl_act,NULL,&gR[start_index[IB]],M,TrM,TRow,TCol,Bmin,Bmax,\
		tb,orb_per_at,IB,NBL,"left");
	
	init_var(&GR[NT+NS],N1);
	c_daxpy(N1,-2.0,((double*)GRl_act)+1,2*(N1+1),(double*)&GR[NT+NS],2);
	c_zcopy(N1,GRl_act,N1+1,&GR[NS],1);

    }
 
    c_zcopy(Bsize*Bsize,GRN1_act,1,GRN1,1);

    if(NBL<NBlock-1){
        calc_SigR(SigRr,T01_gR11,&gR[start_index[NBL+1]],TrM,TRow,TCol,Bmin,Bmax,tb,\
		  orb_per_at,NBL+1,"right");
	get_Gamma(Gammar,SigRr,get_msize(Bmin[NBL],Bmax[NBL],tb,orb_per_at));
    }
   
    delete[] start_index;
    delete[] ipiv;
    delete[] inv_gR;
    delete[] SigR;
    delete[] T01_gR11;
    delete[] GRr_act;
    delete[] GRl_act;
    delete[] GRN1_act;
    delete[] M;
    delete[] TrM;
    delete[] Zl;
    delete[] gR;
    delete TRow;
    delete TCol;
}

/************************************************************************************************/

/*!
\brief Function to extract a diagonal block of the sparse matrix matrix and store it in the full matrix D

D          (output) diagonal block of matrix

N          (output) size of D

Bmin       (input) index of the first atom belonging to each diagonal block of matrix

Bmax       (input) index of the last atom belonging to each diagonal block of matrix

tb         (input) tight-binding parameter

orb_per_at (input) number of orbitals per atom

index      (input) index of the diagonal block to extract

*/

void RGF::extract_diag(CPX *D,int *N,int *Bmin,int *Bmax,int tb,int *orb_per_at,int index)
{
    int i,j,imin,imax,index_i,index_j;

    *N   = get_msize(Bmin[index],Bmax[index],tb,orb_per_at);
    imin = tb/TB*orb_per_at[Bmin[index]];
    imax = tb/TB*orb_per_at[Bmax[index]+1];
    init_var(D,(*N)*(*N));
    
    for(i=imin;i<imax;i++){
        index_i = i-imin;
        for(j=matrix->edge_i[i]-matrix->findx;j<matrix->edge_i[i+1]-matrix->findx;j++){
            index_j = matrix->index_j[j]-matrix->findx-imin;
            if((index_j>=0)&&(index_j<(*N))){
                D[index_i+index_j*(*N)] = matrix->nnz[j];
            }
        }
    }
}

/************************************************************************************************/

/*!
\brief Function to extract two off-diagonal blocks of the sparse matrix matrix, store it in the CSR matrix T01, and conjugate it and store it in the CSC matrix T10

T01        (output) upper off-diagonal block of matrix: matrix_{i,i+1} in CSR format

T10        (output) lower off-diagonal block of matrix: matrix_{i+1,i} in CSC format

Bmin       (input) index of the first atom belonging to each diagonal block of matrix

Bmax       (input) index of the last atom belonging to each diagonal block of matrix

tb         (input) tight-binding parameter

orb_per_at (input) number of orbitals per atom

index      (input) index of the block line i from which the 2 off-diagonal blocks should be extracted

*/

void RGF::extract_not_diag(TCSR<CPX> *T01,TCSC<CPX,int> *T10,int *Bmin,int *Bmax,int tb,\
                           int *orb_per_at,int index)
{
    int i,j,index_i,index_j,n_nonzeros,NR,imin,imax,jshift;
    
    NR        = get_msize(Bmin[index-1],Bmax[index-1],tb,orb_per_at);
    imin      = tb/TB*orb_per_at[Bmin[index-1]];
    imax      = tb/TB*orb_per_at[Bmax[index-1]+1];
    jshift    = tb/TB*orb_per_at[Bmin[index]];

    T01->size = NR;
    T10->size = NR;

    n_nonzeros = 0;
    for(i=imin;i<imax;i++){
        index_i         = i-imin;
        T01->index_i[index_i] = 0;
        T10->index_j[index_i] = 0;
        for(j=matrix->edge_i[i]-matrix->findx;j<matrix->edge_i[i+1]-matrix->findx;j++){
            index_j = matrix->index_j[j]-matrix->findx-jshift;
            if(index_j>=0){
                T01->index_j[n_nonzeros] = index_j+findx;
                T10->index_i[n_nonzeros] = index_j+findx;
                T01->nnz[n_nonzeros]     = matrix->nnz[j];
                T10->nnz[n_nonzeros]     = conj(matrix->nnz[j]);
                T01->index_i[index_i]++;
                T10->index_j[index_i]++;
                n_nonzeros++;
            }
        }
    }
    
    T01->n_nonzeros = n_nonzeros;
    T10->n_nonzeros = n_nonzeros;
    T01->get_row_edge();
    T10->get_column_edge();
}

/************************************************************************************************/

/*!
\brief Function to extract two off-diagonal blocks of the sparse matrix matrix, store it in the CSR matrix T10, and conjugate it and store it in the CSC matrix T01

T01        (output) upper off-diagonal block of matrix: matrix_{i,i+1} in CSC format

T10        (output) lower off-diagonal block of matrix: matrix_{i+1,i} in CSR format

Bmin       (input) index of the first atom belonging to each diagonal block of matrix

Bmax       (input) index of the last atom belonging to each diagonal block of matrix

tb         (input) tight-binding parameter

orb_per_at (input) number of orbitals per atom

index      (input) index of the block line i+1 from which the 2 off-diagonal blocks should be extracted

*/

void RGF::extract_not_diag(TCSC<CPX,int> *T01,TCSR<CPX> *T10,int *Bmin,int *Bmax,int tb,\
                           int *orb_per_at,int index)
{
    int i,j,index_i,index_j,n_nonzeros,NR,imin,imax,jshift;
    
    NR        = get_msize(Bmin[index],Bmax[index],tb,orb_per_at);
    imin      = tb/TB*orb_per_at[Bmin[index]];
    imax      = tb/TB*orb_per_at[Bmax[index]+1];
    jshift    = tb/TB*orb_per_at[Bmin[index-1]];

    T01->size = NR;
    T10->size = NR;

    n_nonzeros = 0;
    for(i=imin;i<imax;i++){
        index_i         = i-imin;
        T10->index_i[index_i] = 0;
        T01->index_j[index_i] = 0;
        for(j=matrix->edge_i[i]-matrix->findx;j<matrix->edge_i[i+1]-matrix->findx;j++){
            index_j = matrix->index_j[j]-matrix->findx-jshift;
            if(index_j<get_msize(Bmin[index-1],Bmax[index-1],tb,orb_per_at)){
                T01->index_i[n_nonzeros] = index_j+findx;
                T10->index_j[n_nonzeros] = index_j+findx;
                T01->nnz[n_nonzeros]     = conj(matrix->nnz[j]);
                T10->nnz[n_nonzeros]     = matrix->nnz[j];
                T01->index_j[index_i]++;
                T10->index_i[index_i]++;
                n_nonzeros++;
            }
        }
    }
    
    T01->n_nonzeros = n_nonzeros;
    T10->n_nonzeros = n_nonzeros;
    T01->get_column_edge();
    T10->get_row_edge();
}

/************************************************************************************************/

/*!
\brief Function to initialize the small retarded Green's Functions gR and the help variable GR_act

gR              (output) small retarded Green's Function

GR_act          (output) help variable to store the current GR

start_index     (output) index of the entries of gR at which each block starts

NBlock          (input) number of blocks in the block-tri-diagonal matrix matrix

Bmin            (input) index of the first atom belonging to each diagonal block of matrix

Bmax            (input) index of the last atom belonging to each diagonal block of matrix

Bsize           (input) size of the largest block of matrix

tb              (input) tight-binding parameter

orb_per_at      (input) number of orbitals per atom

*/

void RGF::init_green(CPX *gR,CPX *GR_act,int *start_index,int NBlock,int *Bmin,int *Bmax,\
                     int Bsize,int tb,int *orb_per_at)
{
    int i,N;
    CPX *vec = new CPX[Bsize];
    
    for(i=0;i<Bsize;i++) vec[i] = CPX(1.0,0.0);
    
    start_index[0] = 0;
    init_var(gR,NBlock*Bsize*Bsize);
    for(i=0;i<NBlock;i++){
        N = get_msize(Bmin[i],Bmax[i],tb,orb_per_at);
        c_zcopy(N,vec,1,&gR[start_index[i]],N+1);
        start_index[i+1] = start_index[i]+N*N;
    }
    
    init_var(GR_act,Bsize*Bsize);
    N = get_msize(Bmin[0],Bmax[0],tb,orb_per_at);
    c_zcopy(N,vec,1,GR_act,N+1);
    
    delete[] vec;

}

/************************************************************************************************/

/*!
\brief Function to calculate gR=inv(matrix_{i,i}-SigR)

gR              (output) current small retarded Green's Function to calculate

inv_gR          (input) help variable to store the inverse of gR

SigR            (input) boundary self-energy

Bmin            (input) index of the first atom belonging to each diagonal block of matrix

Bmax            (input) index of the last atom belonging to each diagonal block of matrix

NBlock          (input) number of blocks in the block-tri-diagonal matrix matrix

ipiv            (input) pivot vector used to inverse inv_gR

tb              (input) tight-binding parameter

orb_per_at      (input) number of orbitals per atom

index           (input) index of the gR that should be computed

*/

void RGF::calc_gR(CPX *gR,CPX *inv_gR,CPX *SigR,int *Bmin,int *Bmax,int NBlock,int *ipiv,\
		  int tb,int *orb_per_at,int index)
{
    int NB,info;
    
    //extract matrix_{i,i} to inv_gR
    extract_diag(inv_gR,&NB,Bmin,Bmax,tb,orb_per_at,index);
    //do inv_gR=inv_gR-SigR
    c_zaxpy(NB*NB,CPX(-1.0,0.0),SigR,1,inv_gR,1);
    //invert inv_gR and store it in gR
    c_zgetrf(NB,NB,inv_gR,NB,ipiv,&info);
    c_zgetrs('N',NB,NB,inv_gR,NB,ipiv,gR,NB,&info);
}

/************************************************************************************************/

/*!
\brief Function to calculate either SigR_{i}=T_{i,i-1}*gR_{i-1}*T_{i-1,i} (side="left") or SigR_{i}=T_{i,i+1}*gR_{i+1}*T_{i+1,i} (side="right")

SigR               (output) boundary self-energy defined above

T01_gR11           (input) help variable to store T_{i,i+-1}*gR_{i+-1}

gR                 (input) small retarded Green's Function gR_{i+-1}

TrM                (input) help variable to transpose the gR matrix 

T01                (input) T_{i,i+1} matrix if side="right", T_{i,i-1} matrix if side="left"

T10                (input) T_{i,i-1} matrix if side="right", T_{i,i+1} matrix if side="left"

Bmin               (input) index of the first atom belonging to each diagonal block of matrix

Bmax               (input) index of the last atom belonging to each diagonal block of matrix

tb                 (input) tight-binding parameter

orb_per_at         (input) number of orbitals per atom

index              (input) index of SigR that should be computed

side               (input) "left" or "right"

*/

void RGF::calc_SigR(CPX *SigR,CPX *T01_gR11,CPX *gR,CPX *TrM,TCSR<CPX> *T01,TCSC<CPX,int> *T10,\
                    int *Bmin,int *Bmax,int tb,int *orb_per_at,int index,char *side)
{
    int Nold,Nnew;
    
    //extraction of T01 and T10
    if(!strcmp(side,"right")){

        Nold = get_msize(Bmin[index],Bmax[index],tb,orb_per_at);
	Nnew = get_msize(Bmin[index-1],Bmax[index-1],tb,orb_per_at);
    
	extract_not_diag(T01,T10,Bmin,Bmax,tb,orb_per_at,index);

    }else{

        Nold = get_msize(Bmin[index-1],Bmax[index-1],tb,orb_per_at);
	Nnew = get_msize(Bmin[index],Bmax[index],tb,orb_per_at);
    
	extract_not_diag(T10,T01,Bmin,Bmax,tb,orb_per_at,index);

    }
    
    //calculate T01_gR11=T_{ii+-1}*gR_{i+-1}
    transpose(TrM,gR,Nold,Nold);
    T01->trans_mat_vec_mult(TrM,T01_gR11,Nold,Nold);
    //T01->mat_vec_mult(gR,T01_gR11,Nold,Nold);

    //calculate T01_gR11*T_{i+-1,i}
    T10->vec_mat_mult(T01_gR11,SigR,Nnew);
}

/************************************************************************************************/

///Function to write a double complex matrix matrix of size NR times NC to a file filename

void RGF::write_matrix(const char *filename,CPX *matrix,int NR,int NC)
{
    int IC,IR;
    ofstream myfile;
    
    myfile.open(filename);
    myfile.precision(8);
    for(IR=0;IR<NR;IR++){
        for(IC=0;IC<NC;IC++){
            myfile<<real(matrix[IR+IC*NR])<<" "<<imag(matrix[IR+IC*NR])<<" ";
        }
        myfile<<"\n";
    }
    myfile.close();
}

/************************************************************************************************/

/*!
\brief Function to calculate GR_{i}=gR_{i}+gR_{i}*T_{i,i+-1}*GR_{i+-1}*T_{i+-1,i}*gR_{i} and GRN1_{i}=-gR_{i}*T_{i,i-1}*GRN1_{i-1}

GR_act         (output) matrix GR_{i}

GRN1_act       (output) matrix GRN1_{i}

gR             (input) small retarded Green's Function

M              (input) matrix to store intermediate results

TrM            (input) matrix to store intermediate results

T01            (input) T_{i,i+1} matrix if side="right", T_{i,i-1} matrix if side="left"

T10            (input) T_{i,i-1} matrix if side="right", T_{i,i+1} matrix if side="left"

Bmin           (input) index of the first atom belonging to each diagonal block of matrix

Bmax           (input) index of the last atom belonging to each diagonal block of matrix

tb             (input) tight-binding parameter

orb_per_at     (input) number of orbitals per atom

index          (input) index of SigR that should be computed

stop           (input) last block for which GRN1 should be computed

side           (input) "left" or "right"

*/

void RGF::calc_GR(CPX *GR_act,CPX *GRN1_act,CPX *gR,CPX *M,CPX *TrM,TCSR<CPX> *T10,\
                  TCSC<CPX,int> *T01,int *Bmin,int *Bmax,int tb,int *orb_per_at,int index,\
		  int stop,char *side)
{

    int N0,N1;

    if(!strcmp(side,"right")){

        N1 = get_msize(Bmin[index],Bmax[index],tb,orb_per_at);
	N0 = get_msize(Bmin[index-1],Bmax[index-1],tb,orb_per_at);

	extract_not_diag(T01,T10,Bmin,Bmax,tb,orb_per_at,index);

	//calculate GRN1_act=-gR*T10*GRN1_act
	if(index<=stop){
	    int NL = get_msize(Bmin[0],Bmax[0],tb,orb_per_at);;
	    transpose(TrM,GRN1_act,N0,NL);
	    T10->trans_mat_vec_mult(TrM,M,NL,N0);
	    //T10->mat_vec_mult(GRN1_act,M,NL,N0);

	    c_zgemm('N','N',N1,NL,N1,CPX(-1.0,0.0),gR,N1,M,N1,CPX(0.0,0.0),GRN1_act,N1);
	}

    }else{

        N1 = get_msize(Bmin[index],Bmax[index],tb,orb_per_at);
	N0 = get_msize(Bmin[index+1],Bmax[index+1],tb,orb_per_at);

	extract_not_diag(T10,T01,Bmin,Bmax,tb,orb_per_at,index+1);
    }

    //calculate GR_act=T10*GR_act*T01
    transpose(TrM,GR_act,N0,N0);
    T10->trans_mat_vec_mult(TrM,M,N0,N0);
    //T10->mat_vec_mult(GR_act,M,N0,N0);
    T01->vec_mat_mult(M,GR_act,N1);

    //calculate M=GR_act*gR and then GR_act=gR*M
    c_zgemm('N','N',N1,N1,N1,CPX(1.0,0.0),GR_act,N1,gR,N1,CPX(0.0,0.0),M,N1);
    c_zgemm('N','N',N1,N1,N1,CPX(1.0,0.0),gR,N1,M,N1,CPX(0.0,0.0),GR_act,N1);

    //calculate GR_act=gR+GR_act
    c_zaxpy(N1*N1,CPX(1.0,0.0),gR,1,GR_act,1);

}

/************************************************************************************************/

///Function to transport a matrix B of size nrow times ncol into a matrix A

void RGF::transpose(CPX *A,CPX *B,int nrow,int ncol)
{
    int i;

    for(i=0;i<ncol;i++) c_zcopy(nrow,&B[i*nrow],1,&A[i],ncol);
}

/************************************************************************************************/

/*!
\brief Calculate the broadening function Gamma=SigR-SigR' of size N from the retarded self-energy SigR
*/

void RGF::get_Gamma(CPX *Gamma,CPX *SigR,int N)
{
    int IR,IC;

    for(IR=0;IR<N;IR++){
        for(IC=0;IC<N;IC++){
            Gamma[IR+IC*N] = CPX(0.0,1.0)*(SigR[IR+IC*N]-conj(SigR[IC+IR*N]));
        }
    }
}

/************************************************************************************************/
