#include "GetSingularities.H"

Singularities::Singularities(TCSR<double> *KohnSham,TCSR<double> *Overlap,c_transport_type parameter_sab)
/**  \brief Initialize array energies and fill it with n_energies energy points at which there are singularities in the DOS, in addition get integration range
 *
 *   \param KohnSham      H matrix in CSR format collected on one node
 *   \param Overlap       S matrix in CSR format collected on one node
 *   \param parameter_sab Parameters
 */
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam); 

    int n_k=parameter_sab.n_kpoint;
    double *k = new double[n_k];
    k[0]=0.0;
    if (n_k>1) for (int i=1;i<n_k;i++) k[i]=i*M_PI/(n_k-1);
    int seq_per_cpu=int(ceil(double(n_k)/nprocs));
    int ndof=Overlap->size_tot/parameter_sab.n_cells;
    double *energies_local = new double[ndof*seq_per_cpu];
    int n_energies_local=0;
    int kpos;
    double energy_gs_local=numeric_limits<double>::max();
    double energy_vbe_local=-numeric_limits<double>::max();
    double energy_gs_current;
    double energy_vbe_current;
    for (int i=0;i<seq_per_cpu;i++)
        if ( (kpos=i+iam*seq_per_cpu)<n_k ) {
            determine_velocities(KohnSham,Overlap,k[kpos],energy_gs_current,energy_vbe_current,energies_local,n_energies_local,parameter_sab);
            energy_gs_local=min(energy_gs_local,energy_gs_current);
            energy_vbe_local=max(energy_vbe_local,energy_vbe_current);
        }
    MPI_Allreduce(&energy_gs_local,&energy_gs,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&energy_vbe_local,&energy_vbe,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    int *n_energies_local_array = new int[nprocs];
    MPI_Allgather(&n_energies_local,1,MPI_INT,n_energies_local_array,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allreduce(&n_energies_local,&n_energies,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    energies = new double[n_energies];
    int* displc = new int[nprocs+1]();
    for (int i=1;i<nprocs+1;i++) displc[i]=displc[i-1]+n_energies_local_array[i-1];
    MPI_Allgatherv(energies_local,n_energies_local,MPI_DOUBLE,energies,n_energies_local_array,displc,MPI_DOUBLE,MPI_COMM_WORLD);
    delete[] displc;
    delete[] energies_local;
    delete[] n_energies_local_array;
    int info;
    c_dlasrt('I',n_energies,energies,&info);

    delete[] k;
}

Singularities::~Singularities()
{
    delete[] energies;
}

int Singularities::determine_velocities(TCSR<double> *KohnSham,TCSR<double> *Overlap,double k_in,double &energy_gs_k,double &energy_vbe_k,double *energies_k,int &n_energies_k,c_transport_type parameter_sab)
/**  \brief Determine energies where velocity vanishes and their number for a given k
 *
 *   \param KohnSham      H matrix in CSR format collected on one node
 *   \param Overlap       S matrix in CSR format collected on one node
 *   \param k_in          position on the k axis from 0 to pi
 *   \param energy_gs_k   will be the lowest energy eigenvalue
 *   \param energy_vbe_k  will be the energy eigenvalue of the highest occupied
 *   \param energies_k    will contain all energies with a vanishing derivative
 *   \param n_energies_k  will give the number of energies
 *   \param parameter_sab includes the threshold for vanishing derivative
 */
{
    double d_zer=0.0;
    CPX z_zer=CPX(d_zer,d_zer);
    CPX z_one=CPX(1.0,d_zer);
    CPX kval=exp(CPX(d_zer,k_in));

    int ndof=Overlap->size_tot/parameter_sab.n_cells;
    int ndofsq=ndof*ndof;
    int bandwidth=parameter_sab.bandwidth;
    int ndofsqbw=ndofsq*(2*bandwidth+1);
    double eps_singularities=parameter_sab.eps_singularities;
    double evoltfactor=parameter_sab.evoltfactor;
    int noccunitcell=parameter_sab.n_occ/parameter_sab.n_cells;

    double *Hreal = new double[ndofsqbw];
    double *Sreal = new double[ndofsqbw];
    CPX *H = new CPX[ndofsqbw]();
    CPX *S = new CPX[ndofsqbw]();
    KohnSham->contactunitcell(Hreal,ndof,bandwidth,1);
    Overlap->contactunitcell(Sreal,ndof,bandwidth,1);
    c_dcopy(ndofsqbw,Hreal,1,(double*)H,2);
    c_dcopy(ndofsqbw,Sreal,1,(double*)S,2);
    delete[] Hreal;
    delete[] Sreal;

    CPX *H_Sum_k = new CPX[ndofsq];
    CPX *S_Sum_k = new CPX[ndofsq];
    c_zcopy(ndofsq,&H[bandwidth*ndofsq],1,H_Sum_k,1);
    c_zcopy(ndofsq,&S[bandwidth*ndofsq],1,S_Sum_k,1);

    for (int ibandw=1;ibandw<=bandwidth;ibandw++) {
        c_zaxpy(ndofsq,CPX(ibandw,d_zer)*pow(kval,+(ibandw)),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_k,1);
        c_zaxpy(ndofsq,CPX(ibandw,d_zer)*pow(kval,+(ibandw)),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_k,1);
        c_zaxpy(ndofsq,CPX(ibandw,d_zer)*pow(kval,-(ibandw)),&H[(bandwidth-ibandw)*ndofsq],1,H_Sum_k,1);
        c_zaxpy(ndofsq,CPX(ibandw,d_zer)*pow(kval,-(ibandw)),&S[(bandwidth-ibandw)*ndofsq],1,S_Sum_k,1);
    }

    double *eigval = new double[ndof];
    double *rwork = new double[3*ndof-2];
    int iinfo;
    CPX twork;
    c_zhegv(1,'V','U',ndof,H_Sum_k,ndof,S_Sum_k,ndof,eigval,&twork,-1,rwork,&iinfo);
    int lwork=int(real(twork));
    CPX *work=new CPX[lwork];
    c_zhegv(1,'V','U',ndof,H_Sum_k,ndof,S_Sum_k,ndof,eigval,work,lwork,rwork,&iinfo);
    if (iinfo) return (cerr<<__LINE__<<endl, EXIT_FAILURE);
    delete[] work;
    delete[] rwork;
    CPX *eigvec = new CPX[ndofsq];
    c_zcopy(ndofsq,H_Sum_k,1,eigvec,1);
    energy_gs_k=evoltfactor*eigval[0];
    energy_vbe_k=evoltfactor*eigval[noccunitcell-1];

    c_zcopy(ndofsq,&H[(bandwidth+1)*ndofsq],1,H_Sum_k,1);
    c_zcopy(ndofsq,&S[(bandwidth+1)*ndofsq],1,S_Sum_k,1);
    for (int ibandw=2;ibandw<=bandwidth;ibandw++) {
        c_zaxpy(ndofsq,CPX(ibandw,d_zer)*pow(kval,+(ibandw-1)),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_k,1);
        c_zaxpy(ndofsq,CPX(ibandw,d_zer)*pow(kval,+(ibandw-1)),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_k,1);
    }
    delete[] H;
    delete[] S;

    CPX *vector   = new CPX[ndof];
    double velocity;
    for (int ivec=0;ivec<ndof;ivec++) {
        c_zgemv('N',ndof,ndof,z_one,H_Sum_k,ndof,&eigvec[ivec*ndof],1,z_zer,vector,1);
        c_zgemv('N',ndof,ndof,CPX(-eigval[ivec],d_zer),S_Sum_k,ndof,&eigvec[ivec*ndof],1,z_one,vector,1);
        velocity=-2.0*imag(z_one/kval*c_zdotc(ndof,&eigvec[ivec*ndof],1,vector,1));
        if (abs(velocity)<eps_singularities) energies_k[n_energies_k++]=evoltfactor*eigval[ivec];
    }
    delete[] vector;
    delete[] H_Sum_k;
    delete[] S_Sum_k;
    delete[] eigvec;
    delete[] eigval;

    return 0;
}
