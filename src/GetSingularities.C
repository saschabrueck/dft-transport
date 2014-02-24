#include "GetSingularities.H"

Singularities::Singularities(c_transport_type parameter_sab)
{
    eps_singularities=parameter_sab.eps_singularities;
    n_k=parameter_sab.n_kpoint;
    n_cells=parameter_sab.n_cells;
    bandwidth=parameter_sab.bandwidth;
    noccunitcell=parameter_sab.n_occ/parameter_sab.n_cells;

    energies_matrix    = NULL;
    derivatives_matrix = NULL;
    curvatures_matrix  = NULL;
}

int Singularities::Execute(TCSR<double> *KohnSham,TCSR<double> *Overlap,int n_mu,double *muvec,double *dopingvec,int *contactvec)
/**  \brief Initialize array energies and fill it with n_energies energy points at which there are singularities in the DOS, in addition get integration range
 *
 *   \param KohnSham      H matrix in CSR format
 *   \param Overlap       S matrix in CSR format
 */
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam); 

    int do_fitting=0;

    double *k = new double[n_k];
    k[0]=0.0;
    if (n_k>1) for (int i=1;i<n_k;i++) k[i]=i*M_PI/(n_k-1);
    int seq_per_cpu=int(ceil(double(n_k)/nprocs));
    ndof=Overlap->size_tot/n_cells;
    int ndofsq=ndof*ndof;
    int ndofsqbw=ndofsq*(2*bandwidth+1);
    if (!iam) {
        energies_matrix = new double[ndof*nprocs*seq_per_cpu];
        derivatives_matrix = new double[ndof*nprocs*seq_per_cpu];
        curvatures_matrix = new double[ndof*nprocs*seq_per_cpu];
        energy_gs=numeric_limits<double>::max();
    }
    double *energies_local = new double[ndof*seq_per_cpu]();
    double *derivatives_local = new double[ndof*seq_per_cpu]();
    double *curvatures_local = new double[ndof*seq_per_cpu]();
    int n_energies;
    int kpos;
    double extrval;
    CPX *H;
    CPX *S;
    for (int i_mu=0;i_mu<n_mu;i_mu++) {
        if (contactvec[i_mu]==1) {
            TCSR<double> *Hcut = new TCSR<double>(KohnSham,0,ndof,0,(bandwidth+1)*ndof);
            TCSR<double> *Hsp = new TCSR<double>(Hcut,MPI_COMM_WORLD);
            Hsp->shift_resize(0,ndof,0,(bandwidth+1)*ndof);
            delete Hcut;
            H = new CPX[ndofsqbw];
            Hsp->sparse_to_cmp_full(&H[bandwidth*ndofsq],ndof,(bandwidth+1)*ndof);
            delete Hsp;
            for (int ibw=1;ibw<=bandwidth;ibw++)
                full_transpose(ndof,ndof,&H[(bandwidth+ibw)*ndofsq],&H[(bandwidth-ibw)*ndofsq]);
            TCSR<double> *Scut = new TCSR<double>(Overlap,0,ndof,0,(bandwidth+1)*ndof);
            TCSR<double> *Ssp = new TCSR<double>(Scut,MPI_COMM_WORLD);
            Ssp->shift_resize(0,ndof,0,(bandwidth+1)*ndof);
            delete Scut;
            S = new CPX[ndofsqbw];
            Ssp->sparse_to_cmp_full(&S[bandwidth*ndofsq],ndof,(bandwidth+1)*ndof);
            delete Ssp;
            for (int ibw=1;ibw<=bandwidth;ibw++)
                full_transpose(ndof,ndof,&S[(bandwidth+ibw)*ndofsq],&S[(bandwidth-ibw)*ndofsq]);
        } else if (contactvec[i_mu]==2) {
            TCSR<double> *Hcut = new TCSR<double>(KohnSham,KohnSham->size_tot-ndof,ndof,KohnSham->size_tot-(bandwidth+1)*ndof,(bandwidth+1)*ndof);
            TCSR<double> *Hsp = new TCSR<double>(Hcut,MPI_COMM_WORLD);
            Hsp->shift_resize(KohnSham->size_tot-ndof,ndof,KohnSham->size_tot-(bandwidth+1)*ndof,(bandwidth+1)*ndof);
            delete Hcut;
            H = new CPX[ndofsqbw];
            Hsp->sparse_to_cmp_full(H,ndof,(bandwidth+1)*ndof);
            delete Hsp;
            for (int ibw=1;ibw<=bandwidth;ibw++)
                full_transpose(ndof,ndof,&H[(bandwidth-ibw)*ndofsq],&H[(bandwidth+ibw)*ndofsq]);
            TCSR<double> *Scut = new TCSR<double>(Overlap,Overlap->size_tot-ndof,ndof,Overlap->size_tot-(bandwidth+1)*ndof,(bandwidth+1)*ndof);
            TCSR<double> *Ssp = new TCSR<double>(Scut,MPI_COMM_WORLD);
            Ssp->shift_resize(Overlap->size_tot-ndof,ndof,Overlap->size_tot-(bandwidth+1)*ndof,(bandwidth+1)*ndof);
            delete Scut;
            S = new CPX[ndofsqbw];
            Ssp->sparse_to_cmp_full(S,ndof,(bandwidth+1)*ndof);
            delete Ssp;
            for (int ibw=1;ibw<=bandwidth;ibw++)
                full_transpose(ndof,ndof,&S[(bandwidth-ibw)*ndofsq],&S[(bandwidth+ibw)*ndofsq]);
        }
        for (int iseq=0;iseq<seq_per_cpu;iseq++)
            if ( (kpos=iseq+iam*seq_per_cpu)<n_k )
                if (determine_velocities(H,S,k[kpos],&energies_local[iseq*ndof],&derivatives_local[iseq*ndof],&curvatures_local[iseq*ndof]))
                    return (LOGCERR, EXIT_FAILURE);
        delete[] H;
        delete[] S;
        MPI_Gather(energies_local,ndof*seq_per_cpu,MPI_DOUBLE,energies_matrix,ndof*seq_per_cpu,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(derivatives_local,ndof*seq_per_cpu,MPI_DOUBLE,derivatives_matrix,ndof*seq_per_cpu,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(curvatures_local,ndof*seq_per_cpu,MPI_DOUBLE,curvatures_matrix,ndof*seq_per_cpu,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if (!iam) {
// DO I HAVE TO DETERMINE THE FERMI LEVEL WITH THE ORIGINAL OR THE BAND SORTED BANDSTRUCTURE ???
            determine_fermi(muvec[i_mu],dopingvec[i_mu]);
            follow_band();
            write_bandstructure(i_mu);
            for (int i=0;i<ndof;i++) {
                if (do_fitting) {
                    if (abs(derivatives_matrix[i])<eps_singularities) energies.push_back(energies_matrix[i]);
                    for (int j=1;j<n_k-2;j++) {
                        if (derivatives_matrix[i+j*ndof]*derivatives_matrix[i+(j+1)*ndof]<0.0) {
                            if (min_parabola(4,&k[j-1],&energies_matrix[i+(j-1)*ndof],ndof,extrval)<1.0/n_k) energies.push_back(extrval);
                        }
                    }
                    if (abs(derivatives_matrix[i+(n_k-1)*ndof])<eps_singularities) energies.push_back(energies_matrix[i+(n_k-1)*ndof]);
                } else {
                    for (int j=0;j<n_k-1;j++) {
                        double xval=-derivatives_matrix[i+j*ndof]/curvatures_matrix[i+j*ndof];
                        if (xval<M_PI/(n_k-1) && xval>=0.0) {
                            energies.push_back(energies_matrix[i+j*ndof]+derivatives_matrix[i+j*ndof]*xval/2.0);
                        }
                    }
                }
            }
            for (int j=0;j<n_k;j++) if (energies_matrix[j*ndof]<energy_gs) energy_gs=energies_matrix[j*ndof];
        }
    }
    if (!iam) std::sort(energies.begin(),energies.end());
    if (!iam) n_energies=energies.size();
    MPI_Bcast(&n_energies,1,MPI_INT,0,MPI_COMM_WORLD);
    energies.resize(n_energies);
    MPI_Bcast(&energies[0],n_energies,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&energy_gs,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&muvec[0],n_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);
    delete[] energies_local;
    delete[] derivatives_local;
    delete[] curvatures_local;
    delete[] k;

    return 0;
}

Singularities::~Singularities()
{
    if (energies_matrix!=NULL) delete[] energies_matrix;
    if (derivatives_matrix!=NULL) delete[] derivatives_matrix;
    if (curvatures_matrix!=NULL) delete[] curvatures_matrix;
}

void Singularities::get_propagating(vector< vector<double> > &prop,const vector<CPX> &evec)
{
    prop.resize(evec.size());
    for (int ie=0;ie<evec.size();ie++) {
        if (!imag(evec[ie])) {
            for (int iband=0;iband<ndof;iband++) {
                for (int i_k=0;i_k<n_k-1;i_k++) {
                   double val=-derivatives_matrix[iband+i_k*ndof]/curvatures_matrix[iband+i_k*ndof];
                   double discri=val*val+2.0/curvatures_matrix[iband+i_k*ndof]*(real(evec[ie])-energies_matrix[iband+i_k*ndof]);
                   if (discri>=0.0) {
                       double x1=(val+sqrt(discri))*(n_k-1)/M_PI;
                       if (x1<1.0 && x1>=0.0) {
                           prop[ie].push_back((i_k+x1)*M_PI/(n_k-1));
                       }
                       double x2=(val-sqrt(discri))*(n_k-1)/M_PI;
                       if (x2<1.0 && x2>=0.0) {
                           prop[ie].push_back((i_k+x2)*M_PI/(n_k-1));
                       }
                   }
                }
            }
        }
    }
}

void Singularities::delete_matrices()
{
    int iam;
    MPI_Comm_rank(MPI_COMM_WORLD,&iam); 
    if (!iam) {
        delete[] energies_matrix;
        delete[] derivatives_matrix;
        delete[] curvatures_matrix;
        energies_matrix=NULL;
        derivatives_matrix=NULL;
        curvatures_matrix=NULL;
    }
}

void Singularities::determine_fermi(double &mu,double doping)
{
    mu=(energies_matrix[noccunitcell-1]+energies_matrix[noccunitcell])/2;
    cout << "Starting Fermi Level " << mu << endl;
    int nocctry=0;
    for (int j=0;j<n_k;j++) {
        int nocc_k=-1;
        while (energies_matrix[++nocc_k+j*ndof]<mu);
        nocctry+=nocc_k;
    }
    double nocciter = (double) nocctry/n_k;
    cout << "Starting Number of Electrons " << nocciter << endl;
    double nocctol=1.0/n_k;
    while (nocciter<noccunitcell+doping-nocctol || nocciter>noccunitcell+doping+nocctol) {
        mu+=(noccunitcell+doping-nocciter)/10.0;
        cout << "New Fermi Level " << mu << endl;
        nocctry=0;
        for (int j=0;j<n_k;j++) {
            int nocc_k=-1;
            while (energies_matrix[++nocc_k+j*ndof]<mu);
            nocctry+=nocc_k;
        }
        nocciter = (double) nocctry/n_k;
        cout << "New Number of Electrons " << nocciter << endl;
    }
}

void Singularities::follow_band()
{
    valarray<double> searchvec(ndof);
    for (int iband=0;iband<ndof;iband++) {
        for (int i_k=1;i_k<n_k;i_k++) {
            c_dcopy(ndof,&energies_matrix[i_k*ndof],1,&searchvec[0],1);
            double E=energies_matrix[iband+(i_k-1)*ndof];
            double D=derivatives_matrix[iband+(i_k-1)*ndof];
            double C=curvatures_matrix[iband+(i_k-1)*ndof];
            searchvec=abs(searchvec-E-D*M_PI/(n_k-1)-C/2.0*M_PI/(n_k-1)*M_PI/(n_k-1));
            int iband2=distance(&searchvec[0],min_element(&searchvec[0]+iband,&searchvec[0]+ndof));
            swap(energies_matrix[iband+i_k*ndof],energies_matrix[iband2+i_k*ndof]);
            swap(derivatives_matrix[iband+i_k*ndof],derivatives_matrix[iband2+i_k*ndof]);
            swap(curvatures_matrix[iband+i_k*ndof],curvatures_matrix[iband2+i_k*ndof]);
        }
    }
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
    extr = a[2]-a[1]*a[1]/a[0]/4.0;
    double resid = 0.0;
    for (p=0;p<n;p++)
        resid+=pow(a[0]*x[p]*x[p]+a[1]*x[p]+a[2]-y[p*disp],2);
    return sqrt(resid);
}

void Singularities::write_bandstructure(int i_mu)
{
    ofstream myfile;
    stringstream mysstream;
    mysstream << "EnergiesWRTk" << i_mu;
    myfile.open(mysstream.str().c_str());
    myfile.precision(15);
    for (int iele=0;iele<n_k*ndof;iele++)
        myfile << energies_matrix[iele] << endl;
    myfile.close();
    mysstream.str("");
    mysstream.clear();
    mysstream << "DerivativesWRTk" << i_mu;
    myfile.open(mysstream.str().c_str());
    myfile.precision(15);
    for (int iele=0;iele<n_k*ndof;iele++)
        myfile << derivatives_matrix[iele] << endl;
    myfile.close();
    mysstream.str("");
    mysstream.clear();
    mysstream << "CurvaturesWRTk" << i_mu;
    myfile.open(mysstream.str().c_str());
    myfile.precision(15);
    for (int iele=0;iele<n_k*ndof;iele++)
        myfile << curvatures_matrix[iele] << endl;
    myfile.close();
}

int Singularities::determine_velocities(CPX *H,CPX *S,double k_in,double *energies_k,double *derivatives_k,double *curvatures_k)
/**  \brief Determine energies where velocity vanishes and their number for a given k
 *
 *   \param H             H matrix
 *   \param S             S matrix
 *   \param k_in          position on the k axis from 0 to pi
 *   \param energies_k    will contain all eigenenergies for given k
 *   \param derivatives_k will contain all derivatives of eigenenergies for given k
 *   \param curvatures_k  will contain all second derivatives of eigenenergies for given k
 */
{
    double d_zer=0.0;
    CPX z_zer=CPX(d_zer,d_zer);
    CPX z_one=CPX(1.0,d_zer);
    CPX z_img=CPX(d_zer,1.0);
    int ndofsq=ndof*ndof;

    CPX kval=exp(CPX(d_zer,k_in));

    CPX *eigvec = new CPX[ndofsq];
    CPX *ovlmat = new CPX[ndofsq];

    double sabtime=get_time(0.0);
    if (eigen(eigvec,ovlmat,H,S,kval,energies_k)) return (LOGCERR, EXIT_FAILURE);
    cout << "TIME FOR SINGULARITY DIAGONALIZATION " << get_time(sabtime) << endl;

    CPX *H_Sum_dk = new CPX[ndofsq]();
    CPX *S_Sum_dk = new CPX[ndofsq]();
    for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
        c_zaxpy(ndofsq,z_img*CPX(ibandw,d_zer)*pow(kval,ibandw),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_dk,1);
        c_zaxpy(ndofsq,z_img*CPX(ibandw,d_zer)*pow(kval,ibandw),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_dk,1);
    }

    CPX *vector = new CPX[ndof];
    for (int ivec=0;ivec<ndof;ivec++) {
        c_zgemv('N',ndof,ndof,z_one,H_Sum_dk,ndof,&eigvec[ivec*ndof],1,z_zer,vector,1);
        c_zgemv('N',ndof,ndof,CPX(-energies_k[ivec],d_zer),S_Sum_dk,ndof,&eigvec[ivec*ndof],1,z_one,vector,1);
        derivatives_k[ivec]=real(c_zdotc(ndof,&eigvec[ivec*ndof],1,vector,1));
    }

    int do_curvatures=1;
    if (do_curvatures) {

        c_zgemm('C','N',ndof,ndof,ndof,z_one,eigvec,ndof,H_Sum_dk,ndof,z_zer,ovlmat,ndof);
        c_zgemm('N','N',ndof,ndof,ndof,z_one,ovlmat,ndof,eigvec,ndof,z_zer,H_Sum_dk,ndof);
        c_zgemm('C','N',ndof,ndof,ndof,z_one,eigvec,ndof,S_Sum_dk,ndof,z_zer,ovlmat,ndof);
        c_zgemm('N','N',ndof,ndof,ndof,z_one,ovlmat,ndof,eigvec,ndof,z_zer,S_Sum_dk,ndof);
 
        for (int ivec=0;ivec<ndof;ivec++) {
            curvatures_k[ivec]=-2.0*derivatives_k[ivec]*real(S_Sum_dk[ivec+ndof*ivec]);
        }

        for (int ivec=0;ivec<ndof;ivec++) {
            c_zscal(ndof,CPX(-energies_k[ivec],d_zer),&S_Sum_dk[ivec*ndof],1);
        }
        c_zaxpy(ndofsq,z_one,S_Sum_dk,1,H_Sum_dk,1);

        for (int ii=0;ii<ndof;ii++) {
            for (int jj=0;jj<ndof;jj++) {
                double delta_e=energies_k[jj]-energies_k[ii];
                if (abs(delta_e)>5.0E-2) {
                    curvatures_k[jj]+=2.0*norm(H_Sum_dk[ii+ndof*jj])/delta_e;
                }
            }
        }

        c_zscal(ndofsq,z_zer,H_Sum_dk,1);
        c_zscal(ndofsq,z_zer,S_Sum_dk,1);
        for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
            c_zaxpy(ndofsq,-CPX(ibandw,d_zer)*CPX(ibandw,d_zer)*pow(kval,ibandw),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_dk,1);
            c_zaxpy(ndofsq,-CPX(ibandw,d_zer)*CPX(ibandw,d_zer)*pow(kval,ibandw),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_dk,1);
        }
 
        for (int ivec=0;ivec<ndof;ivec++) {
            c_zgemv('N',ndof,ndof,z_one,H_Sum_dk,ndof,&eigvec[ivec*ndof],1,z_zer,vector,1);
            c_zgemv('N',ndof,ndof,CPX(-energies_k[ivec],d_zer),S_Sum_dk,ndof,&eigvec[ivec*ndof],1,z_one,vector,1);
            curvatures_k[ivec]+=real(c_zdotc(ndof,&eigvec[ivec*ndof],1,vector,1));
        }

    }

    delete[] H_Sum_dk;
    delete[] S_Sum_dk;
    delete[] vector;
    delete[] eigvec;
    delete[] ovlmat;

    return 0;
}

int Singularities::eigen(CPX *H_Sum_k,CPX *S_Sum_k,CPX *H,CPX *S,CPX kval,double *energies_k)
{
    CPX z_zer=CPX(0.0,0.0);
    CPX z_one=CPX(1.0,0.0);
    int ndofsq=ndof*ndof;

    c_zscal(ndofsq,z_zer,H_Sum_k,1);
    c_zscal(ndofsq,z_zer,S_Sum_k,1);
    for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
        c_zaxpy(ndofsq,pow(kval,ibandw),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_k,1);
        c_zaxpy(ndofsq,pow(kval,ibandw),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_k,1);
    }

    int do_extended=1;
    if (do_extended) {
        int iinfo;
        CPX twork;
        int mfound;
        double abstol=2*c_dlamch('S');
        CPX *zmat=new CPX[ndofsq];
        double *rwork=new double[7*ndof];
        int *iwork=new int[5*ndof];
        int *ifail=new int[ndof];
        c_zhegvx(1,'V','A','U',ndof,H_Sum_k,ndof,S_Sum_k,ndof,0.0,0.0,0,0,abstol,&mfound,energies_k,zmat,ndof,&twork,-1,rwork,iwork,ifail,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        int lwork=int(real(twork));
        CPX *work = new CPX[lwork];
        c_zhegvx(1,'V','A','U',ndof,H_Sum_k,ndof,S_Sum_k,ndof,0.0,0.0,0,0,abstol,&mfound,energies_k,zmat,ndof,work,lwork,rwork,iwork,ifail,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        if (mfound!=ndof) return (LOGCERR, EXIT_FAILURE);
        delete[] work;
        delete[] rwork;
        delete[] iwork;
        delete[] ifail;
        c_zcopy(ndofsq,zmat,1,H_Sum_k,1);
        delete[] zmat;
    } else {
        int iinfo;
        CPX twork;
        double *rwork = new double[3*ndof-2];
        c_zhegv(1,'V','U',ndof,H_Sum_k,ndof,S_Sum_k,ndof,energies_k,&twork,-1,rwork,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        int lwork=int(real(twork));
        CPX *work = new CPX[lwork];
        c_zhegv(1,'V','U',ndof,H_Sum_k,ndof,S_Sum_k,ndof,energies_k,work,lwork,rwork,&iinfo);
        if (iinfo>ndof) {
            cout << "Overlap matrix not positive definite" << endl;
            c_zscal(ndofsq,z_zer,S_Sum_k,1);
            for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
                c_zaxpy(ndofsq,pow(kval,ibandw),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_k,1);
            }
            int iiinfo;
            c_zheev('V','U',ndof,S_Sum_k,ndof,energies_k,&twork,-1,rwork,&iiinfo);
            if (iiinfo) return (LOGCERR, EXIT_FAILURE);
            int lworke=int(real(twork));
            CPX *worke = new CPX[lwork];
            c_zheev('V','U',ndof,S_Sum_k,ndof,energies_k,worke,lworke,rwork,&iiinfo);
            if (iiinfo) return (LOGCERR, EXIT_FAILURE);
            delete[] worke;
            c_zcopy(ndofsq,S_Sum_k,1,H_Sum_k,1);
            for (int ie=0;ie<ndof;ie++) {
                double eigele=energies_k[ie];
                if (eigele<=0.0) {
                    eigele=1.0E-3;
                } else {
                    eigele=sqrt(eigele);
                }
                c_zscal(ndof,eigele,&H_Sum_k[ie*ndof],1);
            }
            c_zgemm('N','C',ndof,ndof,ndof,z_one,H_Sum_k,ndof,H_Sum_k,ndof,z_zer,S_Sum_k,ndof);
            c_zscal(ndofsq,z_zer,H_Sum_k,1);
            for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
                c_zaxpy(ndofsq,pow(kval,ibandw),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_k,1);
            }
            c_zhegv(1,'V','U',ndof,H_Sum_k,ndof,S_Sum_k,ndof,energies_k,work,lwork,rwork,&iiinfo);
            if (iiinfo) return (LOGCERR, EXIT_FAILURE);
        } else if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] work;
        delete[] rwork;
    }

    return 0;
}
