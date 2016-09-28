#include "GetSingularities.H"
#undef min
#undef max
#include <valarray>
#include <limits>
#include "Utilities.H"
#include "ParallelEig.H"

Singularities::Singularities(transport_parameters transport_params,std::vector<contact_type> pcontactvec)
{
    contactvec=pcontactvec;

    dothederivs=(transport_params.real_int_method==real_int_methods::GAUSSCHEBYSHEV);
    eps_singularities=transport_params.eps_singularity_curvatures;
    eps_mu=transport_params.eps_mu;
    n_k=transport_params.n_kpoint;
    Temp=transport_params.temperature;
    n_mu=contactvec.size();
    evfac=transport_params.evoltfactor;

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
    int k_size;
    if (!iam) MPI_Comm_size(equal_bs_rank_comm,&k_size);
    MPI_Bcast(&k_size,1,MPI_INT,0,MPI_COMM_WORLD);
    master_ranks.resize(k_size);
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

int Singularities::Execute(cp2k_csr_interop_type KohnSham,cp2k_csr_interop_type Overlap)
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
    for (int i_mu=0;i_mu<n_mu;i_mu++) {
        int inj_sign=contactvec[i_mu].inj_sign;
        int start=contactvec[i_mu].start_bs;
        int bandwidth=contactvec[i_mu].bandwidth;
        int ndof=contactvec[i_mu].ndof;
        int noccunitcell=floor(contactvec[i_mu].n_ele/2.0);
        TCSR<double> **H = new TCSR<double>*[2*bandwidth+1];
        TCSR<double> **S = new TCSR<double>*[2*bandwidth+1];
        for (int ibw=0;ibw<=bandwidth;ibw++) {
            TCSR<double> *Hcut = new TCSR<double>(KohnSham,start,ndof,start+inj_sign*ibw*ndof,ndof);
            TCSR<double> *Scut = new TCSR<double>(Overlap ,start,ndof,start+inj_sign*ibw*ndof,ndof);
            H[bandwidth+inj_sign*ibw] = new TCSR<double>(Hcut,&master_ranks[0],master_ranks.size(),MPI_COMM_WORLD);
            S[bandwidth+inj_sign*ibw] = new TCSR<double>(Scut,&master_ranks[0],master_ranks.size(),MPI_COMM_WORLD);
            c_dscal(H[bandwidth+inj_sign*ibw]->n_nonzeros,evfac,H[bandwidth+inj_sign*ibw]->nnz,1);
            delete Hcut;
            delete Scut;
        }
        for (int ibw=1;ibw<=bandwidth;ibw++) {
            H[bandwidth-inj_sign*ibw] = new TCSR<double>(H[bandwidth+inj_sign*ibw]);
            H[bandwidth-inj_sign*ibw]->sparse_transpose(H[bandwidth+inj_sign*ibw]);
            S[bandwidth-inj_sign*ibw] = new TCSR<double>(S[bandwidth+inj_sign*ibw]);
            S[bandwidth-inj_sign*ibw]->sparse_transpose(S[bandwidth+inj_sign*ibw]);
        }
        double *energies_local = new double[ndof*seq_per_cpu]();
        double *derivatives_local = new double[ndof*seq_per_cpu]();
        double *curvatures_local = new double[ndof*seq_per_cpu]();
        for (int iseq=0;iseq<seq_per_cpu;iseq++)
            if ( (kpos=iseq+k_rank*seq_per_cpu)<n_k && rank_bs_comm>=0 )
                if (determine_velocities(H,S,k[kpos],&energies_local[iseq*ndof],&derivatives_local[iseq*ndof],&curvatures_local[iseq*ndof],ndof,bandwidth))
                    return (LOGCERR, EXIT_FAILURE);
        for (int ibw=0;ibw<2*bandwidth+1;ibw++) {
            delete H[ibw];
            delete S[ibw];
        }
        delete[] H;
        delete[] S;
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
            if (dothederivs) {
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

double Singularities::determine_fermi(double n_ele,int i_mu)
{
    double mu;
    if (!iam) {
        cout << "Fermi Level / Number of Electrons" << endl;
        int ndof=contactvec[i_mu].ndof;
        mu=(energies_matrix[i_mu][floor(n_ele/2.0)-1]+energies_matrix[i_mu][floor(n_ele/2.0)])/2.0;
        double n_ele_iter = 0.0;
        for (int j=0;j<n_k;j++) {
            for (int i=0;i<ndof;i++) {
                n_ele_iter+=2.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                if (j==0 || j==n_k-1) n_ele_iter-=1.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
            }
        }
        cout << mu << " / " << n_ele_iter << endl;
        double mu_a = mu;
        double n_ele_a = n_ele_iter;
//Find Interval
        while ((n_ele-n_ele_a)*(n_ele-n_ele_iter)>0 && abs(n_ele-n_ele_iter)>eps_mu) {
            mu_a=mu;
            n_ele_a=n_ele_iter;
            mu+=(n_ele-n_ele_iter);
            n_ele_iter=0.0;
            for (int j=0;j<n_k;j++) {
                for (int i=0;i<ndof;i++) {
                    n_ele_iter+=2.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                    if (j==0 || j==n_k-1) n_ele_iter-=1.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                }
            }
cout << "[" << mu_a << "," << mu << "]" << " / " << n_ele_iter << " FIND" << endl;
        }
        double mu_b = mu;
        if (mu_a>mu_b) {
            swap(mu_a,mu_b);
            swap(n_ele_a,n_ele_iter);
        }
//Bisection
        while (abs(mu_a-mu_b)>Temp && abs(n_ele-n_ele_iter)>eps_mu) {
            mu=(mu_a+mu_b)/2.0;
            n_ele_iter=0.0;
            for (int j=0;j<n_k;j++) {
                for (int i=0;i<ndof;i++) {
                    n_ele_iter+=2.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                    if (j==0 || j==n_k-1) n_ele_iter-=1.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                }
            }
            if ((n_ele-n_ele_a)*(n_ele-n_ele_iter)>0) {
                mu_a=mu;
                n_ele_a=n_ele_iter;
            } else {
                mu_b=mu;
            }
cout << "[" << mu_a << "," << mu_b << "]" << " / " << n_ele_iter << " BISECT" << endl;
        }
        mu=(mu_a+mu_b)/2.0;
//Newton
        while (abs(n_ele-n_ele_iter)>eps_mu) {
            double dfermi=0.0;
            for (int j=0;j<n_k;j++) {
                for (int i=0;i<ndof;i++) {
                    dfermi+=2.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,2);
                    if (j==0 || j==n_k-1) dfermi-=1.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,2);
                }
            }
            mu+=(n_ele-n_ele_iter)/dfermi;
            n_ele_iter=0.0;
            for (int j=0;j<n_k;j++) {
                for (int i=0;i<ndof;i++) {
                    n_ele_iter+=2.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                    if (j==0 || j==n_k-1) n_ele_iter-=1.0/(n_k-1)*fermi(energies_matrix[i_mu][i+j*ndof],mu,Temp,0);
                }
            }
cout << mu << " / " << n_ele_iter << " NEWTON" << endl;
        }
        cout << mu << " / " << n_ele_iter << " FERMI" << endl;
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

void Singularities::write_bandstructure(int i_mu,int iter_num,int do_follow_band,int do_debug)
{
    if (!iam) {
        if (do_follow_band) follow_band(i_mu);
        int ndof=contactvec[i_mu].ndof;
        ofstream myfile;
        stringstream mysstream;
        mysstream << "EnergiesWRTk" << i_mu << "_" << iter_num;
        myfile.open(mysstream.str().c_str());
        myfile.precision(15);
        for (int iband=0;iband<ndof;iband++) {
            for (int i_k=0;i_k<n_k;i_k++) {
                myfile << " " << energies_matrix[i_mu][iband+ndof*i_k];
            }
            myfile << endl;
        }
        myfile.close();
        if (do_debug) {
            mysstream.str("");
            mysstream.clear();
            mysstream << "DerivativesWRTk" << i_mu << "_" << iter_num;
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
            mysstream << "CurvaturesWRTk" << i_mu << "_" << iter_num;
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
            mysstream << "ExtremeValues" << i_mu << "_" << iter_num;
            myfile.open(mysstream.str().c_str());
            myfile.precision(15);
            for (uint ival=0;ival<energies_extremum[i_mu].size();ival++) {
                myfile << kval_extremum[i_mu][ival] << " " << energies_extremum[i_mu][ival] << " " << curvatures_extremum[i_mu][ival] << endl;
            }
            myfile.close();
        }
    }
}

int Singularities::determine_imaginary_bandstructure(TCSR<double> **H,TCSR<double> **S,double kimag_in,CPX *lambda,int i_mu)
{
    double kval=exp(-kimag_in);
    int N=contactvec[i_mu].ndof;
    int bandwidth=contactvec[i_mu].bandwidth;
    double* areal = new double[N*N]();
    double* breal = new double[N*N]();

    for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
        H[bandwidth+ibandw]->add_sparse_to_full(areal,N,N,pow(kval,ibandw));
        S[bandwidth+ibandw]->add_sparse_to_full(breal,N,N,pow(kval,ibandw));
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

int Singularities::determine_velocities(TCSR<double> **H,TCSR<double> **S,double k_in,double *energies_k,double *derivatives_k,double *curvatures_k,int ndof,int bandwidth)
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

    CPX *eigvec = NULL;
    CPX *ovlmat = NULL;
    if (!rank_bs_comm) {
        eigvec = new CPX[ndofsq]();
        ovlmat = new CPX[ndofsq]();
        for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
            H[bandwidth+ibandw]->add_sparse_to_cmp_full(eigvec,ndof,ndof,pow(kval,ibandw));
            S[bandwidth+ibandw]->add_sparse_to_cmp_full(ovlmat,ndof,ndof,pow(kval,ibandw));
        }
    }

    if (p_eig(eigvec,ovlmat,energies_k,ndof,bs_comm)) return (LOGCERR, EXIT_FAILURE);

    if (dothederivs) {

        if (rank_bs_comm) eigvec = new CPX[ndofsq]();
        MPI_Bcast(eigvec,ndofsq,MPI_DOUBLE_COMPLEX,0,bs_comm);
 
        CPX *H_Sum_dk = new CPX[ndofsq]();
        CPX *S_Sum_dk = new CPX[ndofsq]();
        if (!rank_bs_comm) {
            for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
                H[bandwidth+ibandw]->add_sparse_to_cmp_full(H_Sum_dk,ndof,ndof,z_img*CPX(ibandw,d_zer)*pow(kval,ibandw));
                S[bandwidth+ibandw]->add_sparse_to_cmp_full(S_Sum_dk,ndof,ndof,z_img*CPX(ibandw,d_zer)*pow(kval,ibandw));
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
            int nprowcol[2]={0,0};
            MPI_Dims_create(size_bs_comm,2,nprowcol);
            int nprow = nprowcol[0];
            int npcol = nprowcol[1];
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
 
 
        }
 
        c_zscal(ndofsq,z_zer,H_Sum_dk,1);
        c_zscal(ndofsq,z_zer,S_Sum_dk,1);
        if (!rank_bs_comm) {
            for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
                H[bandwidth+ibandw]->add_sparse_to_cmp_full(H_Sum_dk,ndof,ndof,-CPX(ibandw,d_zer)*CPX(ibandw,d_zer)*pow(kval,ibandw));
                S[bandwidth+ibandw]->add_sparse_to_cmp_full(S_Sum_dk,ndof,ndof,-CPX(ibandw,d_zer)*CPX(ibandw,d_zer)*pow(kval,ibandw));
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
        delete[] H_Sum_dk;
        delete[] S_Sum_dk;
        if (rank_bs_comm) delete[] eigvec;

    }

    if (!rank_bs_comm) {
        delete[] eigvec;
        delete[] ovlmat;
    }

    return 0;
}

void Singularities::copy_full_to_cp2k_csr(cp2k_csr_interop_type& cp2kCSRmat,double* Pf,int start_i_to,int start_j_to,int length_i,int length_j)
{
    for(int r=0;r<cp2kCSRmat.nrows_local;r++){
        int i=r+cp2kCSRmat.first_row-start_i_to;
        if(i>=0 && i<length_i){
            for(int e=cp2kCSRmat.rowptr_local[r]-1;e<cp2kCSRmat.rowptr_local[r+1]-1;e++){
                int j=cp2kCSRmat.colind_local[e]-1-start_j_to;
                if(j>=0 && j<length_j){
                    cp2kCSRmat.nzvals_local[e]=Pf[i+j*length_i];
                }
            }
        }
    }
}

int Singularities::DensityFromBS(cp2k_csr_interop_type KohnSham,cp2k_csr_interop_type Overlap,cp2k_csr_interop_type* Density,std::vector<double> muvec)
{
    std::vector<double> k(n_k);
    for (int i=1;i<n_k;i++) k[i]=i*M_PI/(n_k-1);
    int seq_per_cpu=int(ceil(double(n_k)/master_ranks.size()));
    int kpos;
    for (int i_mu=0;i_mu<n_mu;i_mu++) {
        int inj_sign=contactvec[i_mu].inj_sign;
        int start=contactvec[i_mu].start_bs;
        int bandwidth=contactvec[i_mu].bandwidth;
        int ndof=contactvec[i_mu].ndof;
        TCSR<double> **H = new TCSR<double>*[2*bandwidth+1];
        TCSR<double> **S = new TCSR<double>*[2*bandwidth+1];
        for (int ibw=0;ibw<=bandwidth;ibw++) {
            TCSR<double> *Hcut = new TCSR<double>(KohnSham,start,ndof,start+inj_sign*ibw*ndof,ndof);
            TCSR<double> *Scut = new TCSR<double>(Overlap ,start,ndof,start+inj_sign*ibw*ndof,ndof);
            H[bandwidth+inj_sign*ibw] = new TCSR<double>(Hcut,&master_ranks[0],master_ranks.size(),MPI_COMM_WORLD);
            S[bandwidth+inj_sign*ibw] = new TCSR<double>(Scut,&master_ranks[0],master_ranks.size(),MPI_COMM_WORLD);
            c_dscal(H[bandwidth+inj_sign*ibw]->n_nonzeros,evfac,H[bandwidth+inj_sign*ibw]->nnz,1);
            delete Hcut;
            delete Scut;
        }
        for (int ibw=1;ibw<=bandwidth;ibw++) {
            H[bandwidth-inj_sign*ibw] = new TCSR<double>(H[bandwidth+inj_sign*ibw]);
            H[bandwidth-inj_sign*ibw]->sparse_transpose(H[bandwidth+inj_sign*ibw]);
            S[bandwidth-inj_sign*ibw] = new TCSR<double>(S[bandwidth+inj_sign*ibw]);
            S[bandwidth-inj_sign*ibw]->sparse_transpose(S[bandwidth+inj_sign*ibw]);
        }
        double **density_from_bs = new double*[2*bandwidth+1];
        for (int ibw=0;ibw<2*bandwidth+1;ibw++) {
            density_from_bs[ibw] = new double[ndof*ndof]();
        }
        for (int iseq=0;iseq<seq_per_cpu;iseq++)
            if ( (kpos=iseq+k_rank*seq_per_cpu)<n_k && rank_bs_comm>=0 )
                if (determine_density_from_bs(H,S,density_from_bs,muvec[i_mu],k[kpos],ndof,bandwidth))
                    return (LOGCERR, EXIT_FAILURE);
        for (int ibw=0;ibw<2*bandwidth+1;ibw++) {
            if (!rank_bs_comm) {
                MPI_Allreduce(MPI_IN_PLACE,density_from_bs[ibw],ndof*ndof,MPI_DOUBLE,MPI_SUM,equal_bs_rank_comm);
            }
            MPI_Bcast(density_from_bs[ibw],ndof*ndof,MPI_DOUBLE,0,MPI_COMM_WORLD);
            for (int i=0;i<bandwidth;i++) {
                int i_bw_pos=start+i*inj_sign*ndof;
                int j_bw_pos=i_bw_pos+(ibw-bandwidth)*ndof;
                int start_sigma=start;
                if (inj_sign==-1) {
                    start_sigma+=-(bandwidth-1)*ndof;
                }
                if (j_bw_pos>=start_sigma && j_bw_pos<start_sigma+bandwidth*ndof) {
                    copy_full_to_cp2k_csr(*Density,density_from_bs[ibw],i_bw_pos,j_bw_pos,ndof,ndof);
                }
            }
            delete[] density_from_bs[ibw];
            delete H[ibw];
            delete S[ibw];
        }
        delete[] density_from_bs;
        delete[] H;
        delete[] S;
    }

    return 0;
}

int Singularities::determine_density_from_bs(TCSR<double> **H,TCSR<double> **S,double **density,double mu,double k_in,int ndof,int bandwidth)
{
    double d_zer=0.0;
    CPX z_zer=CPX(d_zer,d_zer);
    CPX z_one=CPX(1.0,d_zer);
    int ndofsq=ndof*ndof;

    CPX kval=exp(CPX(d_zer,k_in));

    double *energies_k = new double[ndof]();
    CPX *eigvec = NULL;
    CPX *ovlmat = NULL;
    if (!rank_bs_comm) {
        eigvec = new CPX[ndofsq]();
        ovlmat = new CPX[ndofsq]();
        for (int ibandw=-bandwidth;ibandw<=bandwidth;ibandw++) {
            H[bandwidth+ibandw]->add_sparse_to_cmp_full(eigvec,ndof,ndof,pow(kval,ibandw));
            S[bandwidth+ibandw]->add_sparse_to_cmp_full(ovlmat,ndof,ndof,pow(kval,ibandw));
        }
    }

    if (p_eig(eigvec,ovlmat,energies_k,ndof,bs_comm)) return (LOGCERR, EXIT_FAILURE);

    if (!rank_bs_comm) {
        c_zscal(ndofsq,z_zer,ovlmat,1);
        for (int i=0;i<ndof;i++) {
            c_zgemm('N','C',ndof,ndof,1,0.5/(n_k-1)*pow(kval,bandwidth)*fermi(energies_k[i],mu,Temp,0),&eigvec[i*ndof],ndof,&eigvec[i*ndof],ndof,z_one,ovlmat,ndof);
        }
        for (int i=0;i<2*bandwidth+1;i++) {
            double prefac=1.0;
            if (k_in==0 || k_in==M_PI) prefac*=0.5;
            c_daxpy(ndofsq,prefac,(double*)ovlmat,2,density[i],1);
            c_zscal(ndofsq,1.0/kval,ovlmat,1);
        }
        delete[] eigvec;
        delete[] ovlmat;
    }
    delete[] energies_k;

    return 0;
}
