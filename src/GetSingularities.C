#include "GetSingularities.H"

Singularities::Singularities(transport_parameters *parameter_sab,int pn_mu)
{
    eps_singularities=parameter_sab->eps_singularity_curvatures;
    n_k=parameter_sab->n_kpoint;
    n_cells=parameter_sab->n_cells;
    bandwidth=parameter_sab->bandwidth;
    noccunitcell=parameter_sab->n_occ/parameter_sab->n_cells; // THIS IS AN INTEGER DIVISION
    Temp=parameter_sab->temperature;

    n_mu=pn_mu;

    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam); 

    energy_gs= (numeric_limits<double>::max)();
    energies_vb.resize(n_mu,-(numeric_limits<double>::max)());
    energies_cb.resize(n_mu, (numeric_limits<double>::max)());

    energies_extremum.resize(n_mu);
    curvatures_extremum.resize(n_mu);
    kval_extremum.resize(n_mu);

    energies_matrix.resize(n_mu);
    derivatives_matrix.resize(n_mu);
    curvatures_matrix.resize(n_mu);
}

int Singularities::Execute(TCSR<double> *KohnSham,TCSR<double> *Overlap,int *contactvec)
/**  \brief Initialize array energies and fill it with n_energies energy points at which there are singularities in the DOS, in addition get integration range
 *
 *   \param KohnSham      H matrix in CSR format
 *   \param Overlap       S matrix in CSR format
 */
{
    std::vector<double> k(n_k);
    for (int i=1;i<n_k;i++) k[i]=i*M_PI/(n_k-1);
    int seq_per_cpu=int(ceil(double(n_k)/nprocs));
    ndof=Overlap->size_tot/n_cells;
    int ndofsq=ndof*ndof;
    int ndofsqbw=ndofsq*(2*bandwidth+1);
    int kpos;
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
        double *energies_local = new double[ndof*seq_per_cpu]();
        double *derivatives_local = new double[ndof*seq_per_cpu]();
        double *curvatures_local = new double[ndof*seq_per_cpu]();
        for (int iseq=0;iseq<seq_per_cpu;iseq++)
            if ( (kpos=iseq+iam*seq_per_cpu)<n_k )
                if (determine_velocities(H,S,k[kpos],&energies_local[iseq*ndof],&derivatives_local[iseq*ndof],&curvatures_local[iseq*ndof]))
                    return (LOGCERR, EXIT_FAILURE);
        delete[] H;
        delete[] S;
        if (!iam) {
            energies_matrix[i_mu].resize(ndof*nprocs*seq_per_cpu);
            derivatives_matrix[i_mu].resize(ndof*nprocs*seq_per_cpu);
            curvatures_matrix[i_mu].resize(ndof*nprocs*seq_per_cpu);
        }
        MPI_Gather(energies_local,ndof*seq_per_cpu,MPI_DOUBLE,&energies_matrix[i_mu][0],ndof*seq_per_cpu,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(derivatives_local,ndof*seq_per_cpu,MPI_DOUBLE,&derivatives_matrix[i_mu][0],ndof*seq_per_cpu,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(curvatures_local,ndof*seq_per_cpu,MPI_DOUBLE,&curvatures_matrix[i_mu][0],ndof*seq_per_cpu,MPI_DOUBLE,0,MPI_COMM_WORLD);
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

Singularities::~Singularities()
{
}

std::vector< std::vector< std::vector<double> > > Singularities::get_propagating(const std::vector<CPX> &evec)
{
    if (!iam) {
        std::vector< std::vector< std::vector<double> > > prop(n_mu);
        for (int i_mu=0;i_mu<n_mu;i_mu++) {
            prop[i_mu].resize(evec.size());
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
        while ((2.0*noccunitcell+doping-nocc_a)*(2.0*noccunitcell+doping-nocciter)>0) {
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

int Singularities::determine_imaginary_bandstructure(CPX *H,CPX *S,double kimag_in,CPX *lambda)
{
    double kval=exp(-kimag_in);
    int N=ndof;
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
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
if (!worldrank) cout << "TIME FOR SINGULARITY DIAGONALIZATION " << get_time(sabtime) << endl;

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

    int do_extended=0;
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
