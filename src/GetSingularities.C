#include "GetSingularities.H"
#undef min
#undef max
#include <valarray>
#include <limits>
#include "Utilities.H"
#include "p_eig.H"

Singularities::Singularities(transport_parameters *parameter_sab,contact_type *pcontactvec)
{
    contactvec=pcontactvec;

    eps_singularities=parameter_sab->eps_singularity_curvatures;
    n_k=parameter_sab->n_kpoint;
    Temp=parameter_sab->temperature;
    n_mu=parameter_sab->num_contacts;

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam); 

    energy_gs=(numeric_limits<double>::max)();
    energies_vb.resize(n_mu,-(numeric_limits<double>::max)());
    energies_cb.resize(n_mu, (numeric_limits<double>::max)());

    energies_extremum.resize(n_mu);
    curvatures_extremum.resize(n_mu);
    kval_extremum.resize(n_mu);

    energies_matrix.resize(n_mu);
    derivatives_matrix.resize(n_mu);
    curvatures_matrix.resize(n_mu);

    size_bs_comm = max(1,nprocs/n_k);
    int color = iam / size_bs_comm;
    if (iam>=size_bs_comm*n_k) color = MPI_UNDEFINED;
    MPI_Comm_split(MPI_COMM_WORLD,color,iam,&bs_comm);
    MPI_Comm_rank(bs_comm,&rank_bs_comm);
    color = rank_bs_comm;
    if (iam>=size_bs_comm*n_k) {
        color = MPI_UNDEFINED;
        rank_bs_comm = -1;
    }
    MPI_Comm_split(MPI_COMM_WORLD,color,iam,&equal_bs_rank_comm);
    MPI_Comm_rank(equal_bs_rank_comm,&k_rank);
    master_ranks.resize(max(1,nprocs/size_bs_comm));
    if (!rank_bs_comm) {
        MPI_Gather(&iam,1,MPI_INT,&master_ranks[0],1,MPI_INT,0,equal_bs_rank_comm);
    }
    MPI_Bcast(&master_ranks[0],master_ranks.size(),MPI_INT,0,MPI_COMM_WORLD);
}

Singularities::~Singularities()
{
    MPI_Comm_free(&bs_comm);
    MPI_Comm_free(&equal_bs_rank_comm);
}

int Singularities::Execute(TCSR<double> *KohnSham,TCSR<double> *Overlap)
/**  \brief Initialize array energies and fill it with n_energies energy points at which there are singularities in the DOS, in addition get integration range
 *
 *   \param KohnSham      H matrix in CSR format
 *   \param Overlap       S matrix in CSR format
 */
{
    std::vector<double> k(n_k);
    for (int i=1;i<n_k;i++) k[i]=i*M_PI/(n_k-1);
    int seq_per_cpu=int(ceil(double(n_k)/master_ranks.size()));
    int kpos;
    CPX *H = NULL;
    CPX *S = NULL;
    for (int i_mu=0;i_mu<n_mu;i_mu++) {
        int inj_sign=contactvec[i_mu].inj_sign;
//        int contact_start=contactvec[i_mu].start;
        int bandwidth=contactvec[i_mu].bandwidth;
        int ndof=contactvec[i_mu].ndof;
        int noccunitcell=contactvec[i_mu].n_occ;
        int ndofsq=ndof*ndof;
        int ndofsqbw=ndofsq*(2*bandwidth+1);
        if (contactvec[i_mu].inj_sign==+1) {
            TCSR<double> *Hcut = new TCSR<double>(KohnSham,0,ndof,0,(bandwidth+1)*ndof);
            TCSR<double> *Hsp = new TCSR<double>(Hcut,&master_ranks[0],master_ranks.size(),MPI_COMM_WORLD);
            delete Hcut;
            if (!rank_bs_comm) {
                Hsp->shift_resize(0,ndof,0,(bandwidth+1)*ndof);
                H = new CPX[ndofsqbw];
                Hsp->sparse_to_cmp_full(&H[bandwidth*ndofsq],ndof,(bandwidth+1)*ndof);
                for (int ibw=1;ibw<=bandwidth;ibw++)
                    full_transpose(ndof,ndof,&H[(bandwidth+inj_sign*ibw)*ndofsq],&H[(bandwidth-inj_sign*ibw)*ndofsq]);
            }
            delete Hsp;
            TCSR<double> *Scut = new TCSR<double>(Overlap,0,ndof,0,(bandwidth+1)*ndof);
            TCSR<double> *Ssp = new TCSR<double>(Scut,&master_ranks[0],master_ranks.size(),MPI_COMM_WORLD);
            delete Scut;
            if (!rank_bs_comm) {
                Ssp->shift_resize(0,ndof,0,(bandwidth+1)*ndof);
                S = new CPX[ndofsqbw];
                Ssp->sparse_to_cmp_full(&S[bandwidth*ndofsq],ndof,(bandwidth+1)*ndof);
                for (int ibw=1;ibw<=bandwidth;ibw++)
                    full_transpose(ndof,ndof,&S[(bandwidth+inj_sign*ibw)*ndofsq],&S[(bandwidth-inj_sign*ibw)*ndofsq]);
            }
            delete Ssp;
        } else if (contactvec[i_mu].inj_sign==-1) {
            TCSR<double> *Hcut = new TCSR<double>(KohnSham,KohnSham->size_tot-ndof,ndof,KohnSham->size_tot-(bandwidth+1)*ndof,(bandwidth+1)*ndof);
            TCSR<double> *Hsp = new TCSR<double>(Hcut,&master_ranks[0],master_ranks.size(),MPI_COMM_WORLD);
            delete Hcut;
            if (!rank_bs_comm) {
                Hsp->shift_resize(KohnSham->size_tot-ndof,ndof,KohnSham->size_tot-(bandwidth+1)*ndof,(bandwidth+1)*ndof);
                H = new CPX[ndofsqbw];
                Hsp->sparse_to_cmp_full(H,ndof,(bandwidth+1)*ndof);
                for (int ibw=1;ibw<=bandwidth;ibw++)
                    full_transpose(ndof,ndof,&H[(bandwidth+inj_sign*ibw)*ndofsq],&H[(bandwidth-inj_sign*ibw)*ndofsq]);
            }
            delete Hsp;
            TCSR<double> *Scut = new TCSR<double>(Overlap,Overlap->size_tot-ndof,ndof,Overlap->size_tot-(bandwidth+1)*ndof,(bandwidth+1)*ndof);
            TCSR<double> *Ssp = new TCSR<double>(Scut,&master_ranks[0],master_ranks.size(),MPI_COMM_WORLD);
            delete Scut;
            if (!rank_bs_comm) {
                Ssp->shift_resize(Overlap->size_tot-ndof,ndof,Overlap->size_tot-(bandwidth+1)*ndof,(bandwidth+1)*ndof);
                S = new CPX[ndofsqbw];
                Ssp->sparse_to_cmp_full(S,ndof,(bandwidth+1)*ndof);
                for (int ibw=1;ibw<=bandwidth;ibw++)
                    full_transpose(ndof,ndof,&S[(bandwidth+inj_sign*ibw)*ndofsq],&S[(bandwidth-inj_sign*ibw)*ndofsq]);
            }
            delete Ssp;
        }
        double *energies_local = new double[ndof*seq_per_cpu]();
        double *derivatives_local = new double[ndof*seq_per_cpu]();
        double *curvatures_local = new double[ndof*seq_per_cpu]();
        for (int iseq=0;iseq<seq_per_cpu;iseq++)
            if ( (kpos=iseq+k_rank*seq_per_cpu)<n_k && rank_bs_comm>=0 ) {
                if (determine_velocities(H,S,k[kpos],&energies_local[iseq*ndof],&derivatives_local[iseq*ndof],&curvatures_local[iseq*ndof],ndof,bandwidth))
                    return (LOGCERR, EXIT_FAILURE);
}
        if (!rank_bs_comm) {
            delete[] H;
            delete[] S;
        }
        if (!iam) {
            energies_matrix[i_mu].resize(ndof*int(master_ranks.size())*seq_per_cpu);
            derivatives_matrix[i_mu].resize(ndof*int(master_ranks.size())*seq_per_cpu);
            curvatures_matrix[i_mu].resize(ndof*int(master_ranks.size())*seq_per_cpu);
        }
        if (!rank_bs_comm) {
            MPI_Gather(energies_local,ndof*seq_per_cpu,MPI_DOUBLE,&energies_matrix[i_mu][0],ndof*seq_per_cpu,MPI_DOUBLE,0,equal_bs_rank_comm);
            MPI_Gather(derivatives_local,ndof*seq_per_cpu,MPI_DOUBLE,&derivatives_matrix[i_mu][0],ndof*seq_per_cpu,MPI_DOUBLE,0,equal_bs_rank_comm);
            MPI_Gather(curvatures_local,ndof*seq_per_cpu,MPI_DOUBLE,&curvatures_matrix[i_mu][0],ndof*seq_per_cpu,MPI_DOUBLE,0,equal_bs_rank_comm);
        }
        delete[] energies_local;
        delete[] derivatives_local;
        delete[] curvatures_local;
        if (!iam) {
            follow_band(i_mu);
            for (int i=0;i<ndof;i++) {
                for (int j=0;j<n_k-1;j++) {
                    double xval=-derivatives_matrix[i_mu][i+j*ndof]/curvatures_matrix[i_mu][i+j*ndof];
                    if (xval<M_PI/(n_k-1) && xval>=0.0) {
                        energies_extremum[i_mu].push_back(energies_matrix[i_mu][i+j*ndof]+derivatives_matrix[i_mu][i+j*ndof]*xval/2.0);
                        double frval=xval/M_PI*(n_k-1);
                        curvatures_extremum[i_mu].push_back(frval*curvatures_matrix[i_mu][i+j*ndof]+(1-frval)*curvatures_matrix[i_mu][i+(j+1)*ndof]);
                        kval_extremum[i_mu].push_back(k[j]+xval);
                    }
                }
            }
            for (int j=0;j<n_k;j++) if (energies_matrix[i_mu][j*ndof]<energy_gs) energy_gs=energies_matrix[i_mu][j*ndof];
            for (int j=0;j<n_k;j++) if (energies_matrix[i_mu][j*ndof+noccunitcell-1]>energies_vb[i_mu]) energies_vb[i_mu]=energies_matrix[i_mu][j*ndof+noccunitcell-1];
            for (int j=0;j<n_k;j++) if (energies_matrix[i_mu][j*ndof+noccunitcell  ]<energies_cb[i_mu]) energies_cb[i_mu]=energies_matrix[i_mu][j*ndof+noccunitcell  ];
            cout << "Contact " << i_mu << " Valence band edge " << energies_vb[i_mu] << " Conduction band edge " << energies_cb[i_mu] << endl;
        }
    }
    MPI_Bcast(&energy_gs,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&energies_vb[0],n_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&energies_cb[0],n_mu,MPI_DOUBLE,0,MPI_COMM_WORLD);

    return 0;
}

std::vector< std::vector< std::vector<double> > > Singularities::get_propagating(const std::vector<CPX> &evec)
{
    if (!iam) {
        std::vector< std::vector< std::vector<double> > > prop(n_mu);
        for (int i_mu=0;i_mu<n_mu;i_mu++) {
            prop[i_mu].resize(evec.size());
            int ndof=contactvec[i_mu].ndof;
// #pragma omp parallel for num_threads(num_threads_bc)
            for (uint ie=0;ie<evec.size();ie++) {
                if (!imag(evec[ie])) {
                    for (int iband=0;iband<ndof;iband++) {
                        for (int i_k=0;i_k<n_k-1;i_k++) {
                            double val=-derivatives_matrix[i_mu][iband+i_k*ndof]/curvatures_matrix[i_mu][iband+i_k*ndof];
                            double discri=val*val+2.0/curvatures_matrix[i_mu][iband+i_k*ndof]*(real(evec[ie])-energies_matrix[i_mu][iband+i_k*ndof]);
                            if (discri>=0.0) {
                                double x1=(val+sqrt(discri))*(n_k-1)/M_PI;
                                if (x1<1.0 && x1>=0.0) {
                                    prop[i_mu][ie].push_back((i_k+x1)*M_PI/(n_k-1));
                                }
                                double x2=(val-sqrt(discri))*(n_k-1)/M_PI;
                                if (x2<1.0 && x2>=0.0) {
                                    prop[i_mu][ie].push_back((i_k+x2)*M_PI/(n_k-1));
                                }
                            }
                        }
                    }
                }
            }
        }
        return prop;
    } else {
        return std::vector< std::vector< std::vector<double> > > (); 
    }
}

double Singularities::determine_fermi(double doping,int i_mu) //slightly differs from OMEN
{
    double mu;
    if (!iam) {
//doping=2.0;
        double nocctol=max(1.0E-2/n_k*abs(doping),1.0E-4/n_k);//WHAT IS THE MAX PRECISION I CAN GET DEPENDING ON n_k?
        cout << "Fermi Level / Number of Electrons with precision " << nocctol << endl;
        int noccunitcell=contactvec[i_mu].n_occ;
        int ndof=contactvec[i_mu].ndof;
        mu=(energies_matrix[i_mu][noccunitcell-1]+energies_matrix[i_mu][noccunitcell])/2;
        double nocciter = 0.0;
        for (int j=0;j<n_k;j++) {
            for (int i=0;i<ndof;i++) {
                nocciter+=2.0/n_k*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
            }
        }
        cout << mu << " / " << nocciter << endl;
        double mu_a = mu;
        double nocc_a = nocciter;
//Find Interval
        while ((2.0*noccunitcell+doping-nocc_a)*(2.0*noccunitcell+doping-nocciter)>0 && abs(2.0*noccunitcell+doping-nocciter)>nocctol) {
            mu_a=mu;
            nocc_a=nocciter;
            mu+=(2.0*noccunitcell+doping-nocciter);
            nocciter=0.0;
            for (int j=0;j<n_k;j++) {
                for (int i=0;i<ndof;i++) {
                    nocciter+=2.0/n_k*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                }
            }
cout << "[" << mu_a << "," << mu << "]" << " / " << nocciter << " FIND" << endl;
        }
        double mu_b = mu;
        if (mu_a>mu_b) {
            swap(mu_a,mu_b);
            swap(nocc_a,nocciter);
        }
//Bisection
        while (abs(mu_a-mu_b)>K_BOLTZMANN*Temp && abs(2.0*noccunitcell+doping-nocciter)>nocctol) {
            mu=(mu_a+mu_b)/2.0;
            nocciter=0.0;
            for (int j=0;j<n_k;j++) {
                for (int i=0;i<ndof;i++) {
                    nocciter+=2.0/n_k*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                }
            }
            if ((2.0*noccunitcell+doping-nocc_a)*(2.0*noccunitcell+doping-nocciter)>0) {
                mu_a=mu;
                nocc_a=nocciter;
            } else {
                mu_b=mu;
            }
cout << "[" << mu_a << "," << mu_b << "]" << " / " << nocciter << " BISECT" << endl;
        }
        mu=(mu_a+mu_b)/2.0;
//Newton
        while (abs(2.0*noccunitcell+doping-nocciter)>nocctol) {
            double dfermi=0.0;
            for (int j=0;j<n_k;j++) {
                for (int i=0;i<ndof;i++) {
                    dfermi+=2.0/n_k*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,2);
                }
            }
            mu+=(2.0*noccunitcell+doping-nocciter)/dfermi;
            nocciter=0.0;
            for (int j=0;j<n_k;j++) {
                for (int i=0;i<ndof;i++) {
                    nocciter+=2.0/n_k*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                }
            }
cout << mu << " / " << nocciter << " NEWTON" << endl;
        }
        cout << mu << " / " << nocciter << endl;
    }
    MPI_Bcast(&mu,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    return mu;
}

void Singularities::follow_band(int i_mu)
{
    int ndof=contactvec[i_mu].ndof;
    valarray<double> searchvec(ndof);
    for (int iband=0;iband<ndof;iband++) {
        for (int i_k=1;i_k<n_k;i_k++) {
            c_dcopy(ndof,&energies_matrix[i_mu][i_k*ndof],1,&searchvec[0],1);
            double E=energies_matrix[i_mu][iband+(i_k-1)*ndof];
            double D=derivatives_matrix[i_mu][iband+(i_k-1)*ndof];
            double C=curvatures_matrix[i_mu][iband+(i_k-1)*ndof];
            searchvec=abs(searchvec-E-D*M_PI/(n_k-1)-C/2.0*M_PI/(n_k-1)*M_PI/(n_k-1));
            int iband2=distance(&searchvec[0],min_element(&searchvec[0]+iband,&searchvec[0]+ndof));
            swap(energies_matrix[i_mu][iband+i_k*ndof],energies_matrix[i_mu][iband2+i_k*ndof]);
            swap(derivatives_matrix[i_mu][iband+i_k*ndof],derivatives_matrix[i_mu][iband2+i_k*ndof]);
            swap(curvatures_matrix[i_mu][iband+i_k*ndof],curvatures_matrix[i_mu][iband2+i_k*ndof]);
        }
    }
}

void Singularities::write_bandstructure(int i_mu)
{
    if (!iam) {
        int ndof=contactvec[i_mu].ndof;
        ofstream myfile;
        stringstream mysstream;
        mysstream << "EnergiesWRTk" << i_mu;
        myfile.open(mysstream.str().c_str());
        myfile.precision(15);
        for (int iband=0;iband<ndof;iband++) {
            for (int i_k=0;i_k<n_k;i_k++) {
                myfile << " " << energies_matrix[i_mu][iband+ndof*i_k];
            }
            myfile << endl;
        }
        myfile.close();
        mysstream.str("");
        mysstream.clear();
        mysstream << "DerivativesWRTk" << i_mu;
        myfile.open(mysstream.str().c_str());
        myfile.precision(15);
        for (int iband=0;iband<ndof;iband++) {
            for (int i_k=0;i_k<n_k;i_k++) {
                myfile << " " << derivatives_matrix[i_mu][iband+ndof*i_k];
            }
            myfile << endl;
        }
        myfile.close();
        mysstream.str("");
        mysstream.clear();
        mysstream << "CurvaturesWRTk" << i_mu;
        myfile.open(mysstream.str().c_str());
        myfile.precision(15);
        for (int iband=0;iband<ndof;iband++) {
            for (int i_k=0;i_k<n_k;i_k++) {
                myfile << " " << curvatures_matrix[i_mu][iband+ndof*i_k];
            }
            myfile << endl;
        }
        myfile.close();
        mysstream.str("");
        mysstream.clear();
        mysstream << "ExtremeValues" << i_mu;
        myfile.open(mysstream.str().c_str());
        myfile.precision(15);
        for (uint ival=0;ival<energies_extremum[i_mu].size();ival++) {
            myfile << kval_extremum[i_mu][ival] << " " << energies_extremum[i_mu][ival] << " " << curvatures_extremum[i_mu][ival] << endl;
        }
        myfile.close();
    }
}

int Singularities::determine_imaginary_bandstructure(CPX *H,CPX *S,double kimag_in,CPX *lambda,int i_mu)
{
    double kval=exp(-kimag_in);
    int N=contactvec[i_mu].ndof;
    int bandwidth=contactvec[i_mu].bandwidth;
    double* areal = new double[N*N]();
    double* breal = new double[N*N]();

    for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
        c_daxpy(N*N,pow(kval,ibandw),(double*)&H[(bandwidth+ibandw)*N*N],2,areal,1);
        c_daxpy(N*N,pow(kval,ibandw),(double*)&S[(bandwidth+ibandw)*N*N],2,breal,1);
    }

    double* lambda_up_real = new double[N];
    double* lambda_up_imag = new double[N];
    double* lambda_down_real = new double[N];
    double* vl_real = new double[1];
    double* vr_real = new double[1];

    int iinfo;
    double work_test;
    c_dggev('N','N',N,areal,N,breal,N,lambda_up_real,lambda_up_imag,lambda_down_real,vl_real,1,vr_real,1,&work_test,-1,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    int lwork_real=int(work_test);
    double* work_real = new double[lwork_real];
    c_dggev('N','N',N,areal,N,breal,N,lambda_up_real,lambda_up_imag,lambda_down_real,vl_real,1,vr_real,1,work_real,lwork_real,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    delete[] work_real;
    delete[] areal;
    delete[] breal;
    for (int IN=0;IN<N;IN++) {
        if (abs(lambda_down_real[IN])>1.0E-15) {
            lambda[IN] = CPX(lambda_up_real[IN],lambda_up_imag[IN])/lambda_down_real[IN];
        } else {
            lambda[IN] = CPX(INF,INF);
        }
    }
    delete[] lambda_up_real;
    delete[] lambda_up_imag;
    delete[] lambda_down_real;
    delete[] vl_real;
    delete[] vr_real;

    return 0;
}

int Singularities::determine_velocities(CPX *H,CPX *S,double k_in,double *energies_k,double *derivatives_k,double *curvatures_k,int ndof,int bandwidth)
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

    CPX *eigvec = new CPX[ndofsq]();
    CPX *ovlmat = NULL;
    if (!rank_bs_comm) {
        ovlmat = new CPX[ndofsq]();
        for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
            c_zaxpy(ndofsq,pow(kval,ibandw),&H[(bandwidth+ibandw)*ndofsq],1,eigvec,1);
            c_zaxpy(ndofsq,pow(kval,ibandw),&S[(bandwidth+ibandw)*ndofsq],1,ovlmat,1);
        }
    }

    if (p_eig(eigvec,ovlmat,energies_k,ndof,bs_comm)) return (LOGCERR, EXIT_FAILURE);
    MPI_Bcast(eigvec,ndofsq,MPI_DOUBLE_COMPLEX,0,bs_comm);

    CPX *H_Sum_dk = new CPX[ndofsq]();
    CPX *S_Sum_dk = new CPX[ndofsq]();
    if (!rank_bs_comm) {
        for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
            c_zaxpy(ndofsq,z_img*CPX(ibandw,d_zer)*pow(kval,ibandw),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_dk,1);
            c_zaxpy(ndofsq,z_img*CPX(ibandw,d_zer)*pow(kval,ibandw),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_dk,1);
        }
    }
    MPI_Bcast(H_Sum_dk,ndofsq,MPI_DOUBLE_COMPLEX,0,bs_comm);
    MPI_Bcast(S_Sum_dk,ndofsq,MPI_DOUBLE_COMPLEX,0,bs_comm);

    CPX *vector = new CPX[ndof];
    double *tmp = new double[ndof]();
    int seq_per_cpu=int(ceil(double(ndof)/size_bs_comm));
    int ivec;
    for (int iseq=0;iseq<seq_per_cpu;iseq++) {
        if ( (ivec=iseq+rank_bs_comm*seq_per_cpu)<ndof ) {
            c_zgemv('N',ndof,ndof,z_one,H_Sum_dk,ndof,&eigvec[ivec*ndof],1,z_zer,vector,1);
            c_zgemv('N',ndof,ndof,CPX(-energies_k[ivec],d_zer),S_Sum_dk,ndof,&eigvec[ivec*ndof],1,z_one,vector,1);
            tmp[ivec]=real(c_zdotc(ndof,&eigvec[ivec*ndof],1,vector,1));
        }
    }
    MPI_Allreduce(tmp,derivatives_k,ndof,MPI_DOUBLE,MPI_SUM,bs_comm);
 
    c_dscal(ndof,d_zer,tmp,1);

    if (size_bs_comm==1) {

        c_zhemm('L','U',ndof,ndof,z_one,H_Sum_dk,ndof,eigvec,ndof,z_zer,ovlmat,ndof);
        c_zgemm('C','N',ndof,ndof,ndof,z_one,eigvec,ndof,ovlmat,ndof,z_zer,H_Sum_dk,ndof);
        c_zhemm('L','U',ndof,ndof,z_one,S_Sum_dk,ndof,eigvec,ndof,z_zer,ovlmat,ndof);
        c_zgemm('C','N',ndof,ndof,ndof,z_one,eigvec,ndof,ovlmat,ndof,z_zer,S_Sum_dk,ndof);

    } else {

        int icontxt=MPI_Comm_c2f(bs_comm);
        int nprow = size_bs_comm;
        int npcol = 1;
        int myrow, mycol;
        int nbl = 64;
        char gridr[1] = {'R'};
        Cblacs_gridinit(&icontxt,gridr,nprow,npcol);
        Cblacs_gridinfo(icontxt,&nprow,&npcol,&myrow,&mycol);
        int rloc      = max(1,c_numroc(ndof,nbl,myrow,0,nprow));
        int cloc      = c_numroc(ndof,nbl,mycol,0,npcol);

        int iinfo;
        int descH_Sum_dk[9],descS_Sum_dk[9],desceigvec[9];
        c_descinit(descH_Sum_dk,ndof,ndof,ndof,ndof,0,0,icontxt,ndof,&iinfo);
        c_descinit(descS_Sum_dk,ndof,ndof,ndof,ndof,0,0,icontxt,ndof,&iinfo);
        c_descinit(desceigvec,ndof,ndof,ndof,ndof,0,0,icontxt,ndof,&iinfo);
        int descSumLoc[9],descTmpLoc[9],descEigLoc[9];
        c_descinit(descSumLoc,ndof,ndof,nbl,nbl,0,0,icontxt,rloc,&iinfo);
        c_descinit(descTmpLoc,ndof,ndof,nbl,nbl,0,0,icontxt,rloc,&iinfo);
        c_descinit(descEigLoc,ndof,ndof,nbl,nbl,0,0,icontxt,rloc,&iinfo);

        CPX *SumLoc = new CPX[rloc*cloc];
        CPX *TmpLoc = new CPX[rloc*cloc];
        CPX *EigLoc = new CPX[rloc*cloc];

        c_pzgeadd('N',ndof,ndof,z_one,eigvec,1,1,desceigvec,z_zer,EigLoc,1,1,descEigLoc);

        c_pzgeadd('N',ndof,ndof,z_one,H_Sum_dk,1,1,descH_Sum_dk,z_zer,SumLoc,1,1,descSumLoc);
        c_pzgemm('N','N',ndof,ndof,ndof,z_one,SumLoc,1,1,descSumLoc,EigLoc,1,1,descEigLoc,z_zer,TmpLoc,1,1,descTmpLoc);
        c_pzgemm('C','N',ndof,ndof,ndof,z_one,EigLoc,1,1,descEigLoc,TmpLoc,1,1,descTmpLoc,z_zer,SumLoc,1,1,descSumLoc);
        c_pzgeadd('N',ndof,ndof,z_one,SumLoc,1,1,descSumLoc,z_zer,H_Sum_dk,1,1,descH_Sum_dk);

        c_pzgeadd('N',ndof,ndof,z_one,S_Sum_dk,1,1,descS_Sum_dk,z_zer,SumLoc,1,1,descSumLoc);
        c_pzgemm('N','N',ndof,ndof,ndof,z_one,SumLoc,1,1,descSumLoc,EigLoc,1,1,descEigLoc,z_zer,TmpLoc,1,1,descTmpLoc);
        c_pzgemm('C','N',ndof,ndof,ndof,z_one,EigLoc,1,1,descEigLoc,TmpLoc,1,1,descTmpLoc,z_zer,SumLoc,1,1,descSumLoc);
        c_pzgeadd('N',ndof,ndof,z_one,SumLoc,1,1,descSumLoc,z_zer,S_Sum_dk,1,1,descS_Sum_dk);

        delete[] SumLoc;
        delete[] TmpLoc;
        delete[] EigLoc;

        Cblacs_gridexit(icontxt);

    }

    if (!rank_bs_comm) {

        for (int ivec=0;ivec<ndof;ivec++) {
            tmp[ivec]=-2.0*derivatives_k[ivec]*real(S_Sum_dk[ivec+ndof*ivec]);
        }

        c_zcopy(ndofsq,S_Sum_dk,1,ovlmat,1);
        for (int ivec=0;ivec<ndof;ivec++) {
            c_zscal(ndof,CPX(-energies_k[ivec],d_zer),&S_Sum_dk[ivec*ndof],1);
        }
        c_zaxpy(ndofsq,z_one,S_Sum_dk,1,H_Sum_dk,1);

        for (int ii=0;ii<ndof;ii++) {
            for (int jj=0;jj<ndof;jj++) {
                double delta_e=energies_k[jj]-energies_k[ii];
                double delta_o=abs(-conj(H_Sum_dk[jj+ndof*ii])/delta_e+H_Sum_dk[ii+ndof*jj]/delta_e+ovlmat[ii+ndof*jj]);
                if (delta_o<1.0E-12 && abs(delta_e)>1.0E-6) {
                    tmp[jj]+=2.0*norm(H_Sum_dk[ii+ndof*jj])/delta_e;
                }
            }
        }

        delete[] ovlmat;

    }

    c_zscal(ndofsq,z_zer,H_Sum_dk,1);
    c_zscal(ndofsq,z_zer,S_Sum_dk,1);
    if (!rank_bs_comm) {
        for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
            c_zaxpy(ndofsq,-CPX(ibandw,d_zer)*CPX(ibandw,d_zer)*pow(kval,ibandw),&H[(bandwidth+ibandw)*ndofsq],1,H_Sum_dk,1);
            c_zaxpy(ndofsq,-CPX(ibandw,d_zer)*CPX(ibandw,d_zer)*pow(kval,ibandw),&S[(bandwidth+ibandw)*ndofsq],1,S_Sum_dk,1);
        }
    }
    MPI_Bcast(H_Sum_dk,ndofsq,MPI_DOUBLE_COMPLEX,0,bs_comm);
    MPI_Bcast(S_Sum_dk,ndofsq,MPI_DOUBLE_COMPLEX,0,bs_comm);

    for (int iseq=0;iseq<seq_per_cpu;iseq++) {
        if ( (ivec=iseq+rank_bs_comm*seq_per_cpu)<ndof ) {
            c_zgemv('N',ndof,ndof,z_one,H_Sum_dk,ndof,&eigvec[ivec*ndof],1,z_zer,vector,1);
            c_zgemv('N',ndof,ndof,CPX(-energies_k[ivec],d_zer),S_Sum_dk,ndof,&eigvec[ivec*ndof],1,z_one,vector,1);
            tmp[ivec]+=real(c_zdotc(ndof,&eigvec[ivec*ndof],1,vector,1));
        }
    }
    MPI_Allreduce(tmp,curvatures_k,ndof,MPI_DOUBLE,MPI_SUM,bs_comm);
 
    delete[] vector;
    delete[] tmp;
    delete[] eigvec;
    delete[] H_Sum_dk;
    delete[] S_Sum_dk;

    return 0;
}
