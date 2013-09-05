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

    int n_mu=1;
    int contactvec[1]={1};

    int do_fitting=1;

    int n_k=parameter_sab.n_kpoint;
    double *k = new double[n_k];
    k[0]=0.0;
    if (n_k>1) for (int i=1;i<n_k;i++) k[i]=i*M_PI/(n_k-1);
    int seq_per_cpu=int(ceil(double(n_k)/nprocs));
    int ndof=Overlap->size_tot/parameter_sab.n_cells;
    double *energies_matrix;
    double *derivatives_matrix;
    if (!iam) {
        energies_matrix = new double[ndof*nprocs*seq_per_cpu];
        derivatives_matrix = new double[ndof*nprocs*seq_per_cpu];
        energy_gs=numeric_limits<double>::max();
        energy_vbe=-numeric_limits<double>::max();
    }
    double *energies_local = new double[ndof*seq_per_cpu]();
    double *derivatives_local = new double[ndof*seq_per_cpu]();
    int n_energies;
    int kpos;
    double extrval;
    int noccunitcell=parameter_sab.n_occ/parameter_sab.n_cells;
    for (int i_mu=0;i_mu<n_mu;i_mu++) {
        for (int iseq=0;iseq<seq_per_cpu;iseq++)
            if ( (kpos=iseq+iam*seq_per_cpu)<n_k )
                determine_velocities(KohnSham,Overlap,k[kpos],contactvec[i_mu],&energies_local[iseq*ndof],&derivatives_local[iseq*ndof],parameter_sab);
        MPI_Gather(energies_local,ndof*seq_per_cpu,MPI_DOUBLE,energies_matrix,ndof*seq_per_cpu,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(derivatives_local,ndof*seq_per_cpu,MPI_DOUBLE,derivatives_matrix,ndof*seq_per_cpu,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if (!iam && i_mu==0) {
ofstream efile;
efile.open("energiesmatrix");
efile.precision(15);
for (int iwrite=0;iwrite<ndof*nprocs*seq_per_cpu;iwrite++) efile<<energies_matrix[iwrite]<<endl;
efile.close();
efile.open("derivativesmatrix");
efile.precision(15);
for (int iwrite=0;iwrite<ndof*nprocs*seq_per_cpu;iwrite++) efile<<derivatives_matrix[iwrite]<<endl;
efile.close();
            for (int i=0;i<ndof;i++) {
                if (do_fitting) {
                    if (abs(derivatives_matrix[i])<parameter_sab.eps_singularities) energies.push_back(energies_matrix[i]);
                    for (int j=1;j<n_k-2;j++) {
                        if (derivatives_matrix[i+j*ndof]*derivatives_matrix[i+(j+1)*ndof]<0.0) {
cout<<"minpara "<<min_parabola(4,&k[j-1],&energies_matrix[i+(j-1)*ndof],ndof,extrval)<<endl;
                            if (min_parabola(4,&k[j-1],&energies_matrix[i+(j-1)*ndof],ndof,extrval)<1.0/(n_k*n_k)) energies.push_back(extrval);
                        }
                    }
                    if (n_k>1) if (abs((energies_matrix[i+(n_k-1)*ndof]-energies_matrix[i+(n_k-2)*ndof])*(n_k-1)/M_PI)<parameter_sab.eps_singularities)
                        energies.push_back(energies_matrix[i+(n_k-1)*ndof]);
                } else {
                    for (int j=0;j<n_k;j++) {
                        if (abs(derivatives_matrix[i+j*ndof])<parameter_sab.eps_singularities) energies.push_back(energies_matrix[i+j*ndof]);
                    }
                }
            }
            for (int j=0;j<n_k;j++) if (energies_matrix[j*ndof]<energy_gs) energy_gs=energies_matrix[j*ndof];
            for (int j=0;j<n_k;j++) if (energies_matrix[noccunitcell-1+j*ndof]>energy_vbe) energy_vbe=energies_matrix[noccunitcell-1+j*ndof];
        }
    }
    if (!iam) std::sort(energies.begin(),energies.end());
    if (!iam) n_energies=energies.size();
    MPI_Bcast(&n_energies,1,MPI_INT,0,MPI_COMM_WORLD);
    energies.resize(n_energies);
    MPI_Bcast(&energies[0],n_energies,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&energy_gs,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&energy_vbe,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (!iam) {
        delete[] energies_matrix;
        delete[] derivatives_matrix;
    }
    delete[] energies_local;
    delete[] derivatives_local;
    delete[] k;
}

Singularities::~Singularities()
{

}

double Singularities::min_parabola(int n,double *x,double *y,int disp,double& extr)
{
    int i,j,p;
    double a[3];
    double mat[9];
    if (n==3) {
        for (i=0;i<3;i++) {
            a[i]=y[i*disp];
            for (j=0;j<3;j++)
                mat[i+3*j]=pow(x[i],2-j);
        }
    } else if (n>3) {
        for (i=0;i<3;i++) {
            a[i]=0.0;
            for (p=0;p<n;p++)
                a[i]+=pow(x[p],2-i)*y[p*disp];
            for (j=0;j<3;j++) {
                mat[i+3*j]=0.0;
                for (p=0;p<n;p++)
                    mat[i+3*j]+=pow(x[p],4-i-j);
            }
        }
    }
    int info;
    int pivarray[3];
    c_dgetrf(3,3,mat,3,pivarray,&info);
    c_dgetrs('N',3,1,mat,3,pivarray,a,3,&info);
    extr = a[3]-a[2]*a[2]/a[1]/4.0;
    double resid = 0.0;
    for (p=0;p<n;p++)
        resid+=pow(a[1]*x[p]*x[p]+a[2]*x[p]+a[3]-y[p*disp],2);
    return sqrt(resid);
}

int Singularities::determine_velocities(TCSR<double> *KohnSham,TCSR<double> *Overlap,double k_in,int contact,double *energies_k,double *derivatives_k,c_transport_type parameter_sab)
/**  \brief Determine energies where velocity vanishes and their number for a given k
 *
 *   \param KohnSham      H matrix in CSR format collected on one node
 *   \param Overlap       S matrix in CSR format collected on one node
 *   \param k_in          position on the k axis from 0 to pi
 *   \param energies_k    will contain all eigenenergies for given k
 *   \param derivatives_k will contain all derivatives of eigenenergies for given k
 *   \param parameter_sab Parameters
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

    double *Hreal = new double[ndofsqbw];
    double *Sreal = new double[ndofsqbw];
    CPX *H = new CPX[ndofsqbw]();
    CPX *S = new CPX[ndofsqbw]();
    KohnSham->contactunitcell(Hreal,ndof,bandwidth,contact);
    Overlap->contactunitcell(Sreal,ndof,bandwidth,contact);
    c_dcopy(ndofsqbw,Hreal,1,(double*)H,2);
    c_dcopy(ndofsqbw,Sreal,1,(double*)S,2);
    delete[] Hreal;
    delete[] Sreal;

    CPX *H_Sum_k = new CPX[ndofsq]();
    CPX *S_Sum_k = new CPX[ndofsq]();
    for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
        c_zaxpy(ndofsq,pow(kval,ibandw),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_k,1);
        c_zaxpy(ndofsq,pow(kval,ibandw),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_k,1);
    }

    double *rwork = new double[3*ndof-2];
    int iinfo;
    CPX twork;
    c_zhegv(1,'V','U',ndof,H_Sum_k,ndof,S_Sum_k,ndof,energies_k,&twork,-1,rwork,&iinfo);
    int lwork=int(real(twork));
    CPX *work=new CPX[lwork];
    c_zhegv(1,'V','U',ndof,H_Sum_k,ndof,S_Sum_k,ndof,energies_k,work,lwork,rwork,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    delete[] work;
    delete[] rwork;

    CPX *H_Sum_dk = new CPX[ndofsq]();
    CPX *S_Sum_dk = new CPX[ndofsq]();
    for (int ibandw=1;ibandw<=bandwidth;ibandw++) {
        c_zaxpy(ndofsq,CPX(ibandw,d_zer)*pow(kval,ibandw),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_dk,1);
        c_zaxpy(ndofsq,CPX(ibandw,d_zer)*pow(kval,ibandw),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_dk,1);
    }
    delete[] H;
    delete[] S;

    CPX *vector = new CPX[ndof];
    for (int ivec=0;ivec<ndof;ivec++) {
        c_zgemv('N',ndof,ndof,z_one,H_Sum_dk,ndof,&H_Sum_k[ivec*ndof],1,z_zer,vector,1);
        c_zgemv('N',ndof,ndof,CPX(-energies_k[ivec],d_zer),S_Sum_dk,ndof,&H_Sum_k[ivec*ndof],1,z_one,vector,1);
        derivatives_k[ivec]=-2.0*imag(c_zdotc(ndof,&H_Sum_k[ivec*ndof],1,vector,1));
    }
    delete[] H_Sum_dk;
    delete[] S_Sum_dk;

    for (int ivec=0;ivec<ndof;ivec++) {
        c_zgemv('N',ndof,ndof,z_one,S_Sum_k,ndof,&H_Sum_k[ivec*ndof],1,z_zer,vector,1);
        derivatives_k[ivec]/=real(c_zdotc(ndof,&H_Sum_k[ivec*ndof],1,vector,1));
    }
    delete[] vector;
    delete[] H_Sum_k;
    delete[] S_Sum_k;

    return 0;
}
