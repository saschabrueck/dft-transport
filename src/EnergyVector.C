#include "Utilities.H"
#include "Density.H"
#include "GetSingularities.H"
#include "Quadrature.H"
#include "EnergyVector.H"
#include "pole.hpp"
#include <iterator>
#include <limits>

Energyvector::Energyvector()
{
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);
}

Energyvector::~Energyvector()
{
}

int Energyvector::Execute(cp2k_csr_interop_type Overlap,cp2k_csr_interop_type KohnSham,cp2k_csr_interop_type *P,cp2k_csr_interop_type *PImag,std::vector<double> muvec,std::vector<contact_type> contactvec,std::vector<int> Bsizes,std::vector<int> orb_per_at,double *Vatom,double *rho_atom,transport_parameters transport_params)
{
    std::vector<CPX> energyvector;
    std::vector<CPX> energyvector_real;
    std::vector<CPX> stepvector;
    std::vector<CPX> stepvector_real;
    std::vector< std::vector<int> > propagating_sizes;
    if (determine_energyvector(energyvector,stepvector,energyvector_real,stepvector_real,propagating_sizes,KohnSham,Overlap,muvec,contactvec,transport_params)) return (LOGCERR, EXIT_FAILURE);
    if (!iam) cout << "Size of Energyvectors " << energyvector.size() << " " << energyvector_real.size() << endl;

    distribution_methods::distribution_method_type distribution_method_cc = distribution_methods::NO_DISTRIBUTION;
    if (energyvector.size()) {
        if (transport_params.inv_solver_method==inv_solver_methods::FULL)         distribution_method_cc = distribution_methods::CEILING_DISTRIBUTION;
        else if (transport_params.inv_solver_method==inv_solver_methods::PEXSI)   distribution_method_cc = distribution_methods::FLOOR_DISTRIBUTION;
        else if (transport_params.inv_solver_method==inv_solver_methods::PARDISO) distribution_method_cc = distribution_methods::MASTER_DISTRIBUTION;
        else if (transport_params.inv_solver_method==inv_solver_methods::RGF)     distribution_method_cc = distribution_methods::MASTER_DISTRIBUTION;
        else return (LOGCERR, EXIT_FAILURE);
    }

    distribution_methods::distribution_method_type distribution_method = distribution_methods::NO_DISTRIBUTION;
    if (energyvector_real.size()) {
        distribution_method = distribution_methods::CEILING_DISTRIBUTION;
        if (transport_params.lin_solver_method==lin_solver_methods::SPLITSOLVE) distribution_method = distribution_methods::SPLITSOLVE_DISTRIBUTION;
        if (transport_params.lin_solver_method==lin_solver_methods::PARDISO)    distribution_method = distribution_methods::MASTER_DISTRIBUTION;
        if (transport_params.lin_solver_method==lin_solver_methods::UMFPACK)    distribution_method = distribution_methods::MASTER_DISTRIBUTION;
        if (transport_params.lin_solver_method==lin_solver_methods::SUPERLU && distribution_method_cc==distribution_methods::FLOOR_DISTRIBUTION) distribution_method = distribution_methods::FLOOR_DISTRIBUTION;
        if (transport_params.lin_solver_method==lin_solver_methods::MUMPS   && distribution_method_cc==distribution_methods::FLOOR_DISTRIBUTION) distribution_method = distribution_methods::FLOOR_DISTRIBUTION;
    }

    if (transport_params.tasks_per_point==transport_params.tasks_per_point_cc && distribution_method==distribution_method_cc) {
        if (distribute_and_execute(energyvector,stepvector,energyvector_real,stepvector_real,propagating_sizes,distribution_method,transport_params.tasks_per_point,Overlap,KohnSham,P,PImag,muvec,contactvec,Bsizes,orb_per_at,Vatom,rho_atom,transport_params)) return (LOGCERR, EXIT_FAILURE);
    } else {
        if (energyvector.size()) if (distribute_and_execute(energyvector,stepvector,std::vector<CPX>(),std::vector<CPX>(),propagating_sizes,distribution_method_cc,transport_params.tasks_per_point_cc,Overlap,KohnSham,P,PImag,muvec,contactvec,Bsizes,orb_per_at,Vatom,rho_atom,transport_params)) return (LOGCERR, EXIT_FAILURE);
        double *Ptmp = NULL;
        if (energyvector.size() && energyvector_real.size()) {
            Ptmp = new double[P->nze_local];
            c_dcopy(P->nze_local,P->nzvals_local,1,Ptmp,1);
            c_dscal(P->nze_local,0.0,P->nzvals_local,1);
        }
        if (energyvector_real.size()) if (distribute_and_execute(std::vector<CPX>(),std::vector<CPX>(),energyvector_real,stepvector_real,propagating_sizes,distribution_method,transport_params.tasks_per_point,Overlap,KohnSham,P,PImag,muvec,contactvec,Bsizes,orb_per_at,Vatom,rho_atom,transport_params)) return (LOGCERR, EXIT_FAILURE);
        if (energyvector.size() && energyvector_real.size()) {
            c_daxpy(P->nze_local,1.0,Ptmp,1,P->nzvals_local,1);
            delete[] Ptmp;
            Ptmp = NULL;
        }
    }

    return 0;
}

int Energyvector::distribute_and_execute(std::vector<CPX> energyvector,std::vector<CPX> stepvector,std::vector<CPX> energyvector_real,std::vector<CPX> stepvector_real,std::vector< std::vector<int> > propagating_sizes,distribution_methods::distribution_method_type distribution_method,int tasks_per_point,cp2k_csr_interop_type Overlap,cp2k_csr_interop_type KohnSham,cp2k_csr_interop_type *P,cp2k_csr_interop_type *PImag,std::vector<double> muvec,std::vector<contact_type> contactvec,std::vector<int> Bsizes,std::vector<int> orb_per_at,double *Vatom,double *rho_atom,transport_parameters transport_params)
{
double sabtime;
    std::vector<int> Tsizes = get_tsizes(distribution_method,Overlap.nrows_total-transport_params.cutl-transport_params.cutr,Bsizes,orb_per_at,transport_params.gpus_per_point,tasks_per_point);
    if (!Tsizes.size()) return (LOGCERR, EXIT_FAILURE);
sabtime=get_time(0.0);
    MPI_Comm matrix_comm;
    TCSR<double> *OverlapCollect  = new TCSR<double>(Overlap ,MPI_COMM_WORLD,&Tsizes[0],Tsizes.size(),transport_params.cutl,transport_params.cutr,&matrix_comm);
    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD,&Tsizes[0],Tsizes.size(),transport_params.cutl,transport_params.cutr,&matrix_comm);
    if (transport_params.cutl || transport_params.cutr) {
        TCSR<double> *OverlapCollectCut  = new TCSR<double>(OverlapCollect ,0,OverlapCollect->size_tot,transport_params.cutl,OverlapCollect->size_tot);
        TCSR<double> *KohnShamCollectCut = new TCSR<double>(KohnShamCollect,0,OverlapCollect->size_tot,transport_params.cutl,OverlapCollect->size_tot);
        delete OverlapCollect;
        delete KohnShamCollect;
        OverlapCollect  = OverlapCollectCut;
        KohnShamCollect = KohnShamCollectCut;
    }
    c_dscal(KohnShamCollect->n_nonzeros,transport_params.evoltfactor,KohnShamCollect->nnz,1);
    TCSR<double> *DensReal = new TCSR<double>(OverlapCollect);
    c_dscal(DensReal->n_nonzeros,0.0,DensReal->nnz,1);
    TCSR<double> *DensImag = NULL;
    if (transport_params.cp2k_method!=cp2k_methods::LOCAL_SCF) {
        DensImag = new TCSR<double>(OverlapCollect);
        c_dscal(DensImag->n_nonzeros,0.0,DensImag->nnz,1);
    }
if (!iam) cout << "TIME FOR DISTRIBUTING MATRICES " << get_time(sabtime) << endl;

    std::vector<int> atom_of_bf;
    if (transport_params.cp2k_method==cp2k_methods::LOCAL_SCF) {
        int atom=0;
        for (int ibf=0;ibf<OverlapCollect->size_tot;ibf++) {
            if (ibf==orb_per_at[atom+1]) ++atom;
            atom_of_bf.push_back(atom);
        }
        if (++atom != transport_params.n_atoms) return (LOGCERR, EXIT_FAILURE);
        KohnShamCollect->add_pot(OverlapCollect,&atom_of_bf[0],Vatom);
    }

sabtime=get_time(0.0);
    std::vector<double> transmission(energyvector_real.size(),0.0);
    int matrix_size,matrix_rank;
    MPI_Comm_size(matrix_comm,&matrix_size);
    MPI_Comm_rank(matrix_comm,&matrix_rank);
    int matrix_id = iam/matrix_size;
    int n_mat_comm = nprocs/matrix_size;
    int transmission_warning=0;
    int propagating_warning=0;
    int degeneracy_warning=0;
    energyvector.insert(energyvector.end(),energyvector_real.begin(),energyvector_real.end());
    stepvector.insert(stepvector.end(),stepvector_real.begin(),stepvector_real.end());
    for (int iseq=0;iseq<int(ceil(double(energyvector.size())/n_mat_comm));iseq++) {
        int jpos=matrix_id+iseq*n_mat_comm;
        int propos=jpos-(energyvector.size()-energyvector_real.size());
        if (jpos<int(energyvector.size())) {
            if (abs(stepvector[jpos])>0.0) {
                std::vector<result_type> resvec(muvec.size());
                transport_methods::transport_method_type method;
                if (propos>=0) {
                    if (transport_params.negf_solver) {
                        method=transport_methods::NEGF;
                    } else {
                        method=transport_methods::WF;
                    }
                } else {
                    if (transport_params.obc) {
                        method=transport_methods::GF;
                    } else {
                        method=transport_methods::EQ;
                    }
                }
                if (density(KohnShamCollect,OverlapCollect,DensReal,DensImag,energyvector[jpos],stepvector[jpos],method,muvec,contactvec,resvec,Bsizes,orb_per_at,transport_params,matrix_comm)) return (LOGCERR, EXIT_FAILURE);
                if (!matrix_rank && propos>=0) {
                    for (uint i_mu=0;i_mu<muvec.size();i_mu++) {
                        if (resvec[i_mu].npro!=propagating_sizes[propos][i_mu] && transport_params.real_int_method==real_int_methods::GAUSSCHEBYSHEV) propagating_warning++;
                        if (resvec[i_mu].eigval_degeneracy>=0) degeneracy_warning++;
                        if (resvec[i_mu].rcond<numeric_limits<double>::epsilon()) return (LOGCERR, EXIT_FAILURE);
                    }
                    bool transmission_difference=abs(abs(resvec[0].transm)-abs(resvec[1].transm))/max(1.0,min(abs(resvec[0].transm),abs(resvec[1].transm)))<0.1;
                    bool transmission_magnitude=abs(resvec[0].transm)<*max_element(propagating_sizes[propos].begin(),propagating_sizes[propos].end())*10.0 || transport_params.real_int_method!=real_int_methods::GAUSSCHEBYSHEV;
                    if (transmission_difference && transmission_magnitude) {
                        transmission[propos]=resvec[0].transm;
                    } else {
                        transmission_warning++;
                    }
                }
            }
        }
        if (!iam) cout << "Finished " << int((iseq+1)*100.0/ceil(double(energyvector.size())/n_mat_comm)) << "%" << endl;
    }
if (!iam) cout << "TIME FOR DENSITY " << get_time(sabtime) << endl;
    MPI_Allreduce(MPI_IN_PLACE,&transmission_warning,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) if (transmission_warning) cout << "WARNING: INCORRECT TRANSMISSION FOR " << int(transmission_warning*100.0/energyvector.size()) << "%" << " OF THE ENERGY POINTS" << endl;
    MPI_Allreduce(MPI_IN_PLACE,&propagating_warning,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) if (propagating_warning) cout << "WARNING: DEVIATION IN NUMBER OF BANDS FOR " << int(propagating_warning*100.0/energyvector.size()/muvec.size()) << "%" << " OF THE ENERGY POINTS" << endl;
    MPI_Allreduce(MPI_IN_PLACE,&degeneracy_warning,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) if (degeneracy_warning) cout << "WARNING: DEGENERACY OF THE EIGENVECTORS FOR " << int(degeneracy_warning*100.0/energyvector.size()/muvec.size()) << "%" << " OF THE ENERGY POINTS" << endl;
    MPI_Allreduce(MPI_IN_PLACE,&transmission[0],energyvector_real.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (energyvector_real.size()) write_transmission_current(energyvector_real,stepvector_real,transmission,muvec,transport_params);

    if (transport_params.cutl || transport_params.cutr) {
        TCSR<double> *DensRealCollect = new TCSR<double>(*P,MPI_COMM_WORLD,&Tsizes[0],Tsizes.size(),transport_params.cutl,transport_params.cutr,&matrix_comm);
        c_dscal(DensRealCollect->n_nonzeros,double(Tsizes.size())/double(nprocs),DensRealCollect->nnz,1);
        DensReal->copy_shifted(DensRealCollect,0,OverlapCollect->size_tot,transport_params.cutl,OverlapCollect->size_tot);
        delete DensReal;
        DensReal = DensRealCollect;
#ifdef HAVE_PIMAG
        if (transport_params.cp2k_method!=cp2k_methods::LOCAL_SCF) {
            TCSR<double> *DensImagCollect = new TCSR<double>(*PImag,MPI_COMM_WORLD,&Tsizes[0],Tsizes.size(),transport_params.cutl,transport_params.cutr,&matrix_comm);
            c_dscal(DensImagCollect->n_nonzeros,double(Tsizes.size())/double(nprocs),DensImagCollect->nnz,1);
            DensImag->copy_shifted(DensImagCollect,0,OverlapCollect->size_tot,transport_params.cutl,OverlapCollect->size_tot);
            delete DensImag;
            DensImag = DensImagCollect;
        }
#endif
    }
    if (transport_params.cp2k_method!=cp2k_methods::LOCAL_SCF) {
        if (!(transport_params.cp2k_method==cp2k_methods::TRANSMISSION && transport_params.extra_scf)) {
            DensReal->distribute_back(*P,MPI_COMM_WORLD,&Tsizes[0],Tsizes.size(),transport_params.cutl,transport_params.cutr,matrix_comm);
        }
#ifdef HAVE_PIMAG
        DensImag->distribute_back(*PImag,MPI_COMM_WORLD,&Tsizes[0],Tsizes.size(),transport_params.cutl,transport_params.cutr,matrix_comm);
#endif
    } else {
        DensReal->atom_allocate(OverlapCollect,&atom_of_bf[0],rho_atom,2.0);
        MPI_Allreduce(MPI_IN_PLACE,rho_atom,transport_params.n_atoms,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
    delete KohnShamCollect;
    delete OverlapCollect;
    delete DensReal;
    if (transport_params.cp2k_method!=cp2k_methods::LOCAL_SCF) {
        delete DensImag;
    }
    MPI_Comm_free(&matrix_comm);

    return 0;
}

std::vector<int> Energyvector::get_tsizes(distribution_methods::distribution_method_type distribution_method,int nrows_total_cut,std::vector<int> Bsizes,std::vector<int> orb_per_at,int gpus_per_point,int tasks_per_point)
{
    std::vector<int> Tsizes;
    if (distribution_method==distribution_methods::SPLITSOLVE_DISTRIBUTION) {
        std::vector<int> Bmin(1,0);
        for (uint i=1;i<Bsizes.size();i++) Bmin.push_back(Bmin[i-1]+Bsizes[i-1]);
        std::vector<int> Fsizes;
        for (uint i=0;i<Bsizes.size();i++) Fsizes.push_back(orb_per_at[Bmin[i]+Bsizes[i]]-orb_per_at[Bmin[i]]);
        Tsizes.assign(tasks_per_point,0);
        double load_max=pow(double((nrows_total_cut)/gpus_per_point),3);
        int stride=tasks_per_point/gpus_per_point;
        int iblock=0;
        for (int igpu=0;igpu<gpus_per_point;igpu++) {
            int Gsize=0;
            Gsize+=Fsizes[iblock++];
            Gsize+=Fsizes[iblock++];
            while (pow(double(Gsize+Fsizes[iblock]),3)<=load_max && int(Fsizes.size())-iblock-1>=2*(gpus_per_point-igpu-1) && iblock<int(Fsizes.size())) {
                Gsize+=Fsizes[iblock++];
            }
            if (igpu==gpus_per_point-1) {
                while (iblock<int(Fsizes.size())) {
                    Gsize+=Fsizes[iblock++];
                }
            }
            int loc=stride*(2*(igpu/2)+1)-1;
            if (igpu%2) loc++;
            Tsizes[loc]=Gsize;
        }
    } else if (distribution_method==distribution_methods::FLOOR_DISTRIBUTION) {
        int loc_size=int(floor(double(nrows_total_cut)/double(tasks_per_point)));
        Tsizes.assign(tasks_per_point,loc_size);
        Tsizes[tasks_per_point-1]=nrows_total_cut-(tasks_per_point-1)*loc_size;
    } else if (distribution_method==distribution_methods::CEILING_DISTRIBUTION) {
        int loc_size=int(ceil(double(nrows_total_cut)/double(tasks_per_point)));
        Tsizes.assign(tasks_per_point,loc_size);
        Tsizes[tasks_per_point-1]=nrows_total_cut-(tasks_per_point-1)*loc_size;
    } else if (distribution_method==distribution_methods::MASTER_DISTRIBUTION) {
        Tsizes.assign(tasks_per_point,0);
        Tsizes[0]=nrows_total_cut;
    }
    return Tsizes;
}

int Energyvector::write_transmission_current(std::vector<CPX> energyvector,std::vector<CPX> stepvector,std::vector<double> transmission,std::vector<double> muvec,transport_parameters transport_params)
{
    if (!iam) {
        stringstream mysstream;
        mysstream << "Transmission_" << transport_params.cp2k_scf_iter;
        ofstream myfile(mysstream.str().c_str());
        myfile.precision(15);
        double current = 0.0;
        for (uint iele=0;iele<energyvector.size();iele++) {
            if (!imag(energyvector[iele])) {
                myfile << real(energyvector[iele]) << " " << real(stepvector[iele]) << " " << transmission[iele] << endl;
                double diff_fermi=fermi(real(energyvector[iele]),muvec[0],transport_params.temperature,0)-fermi(real(energyvector[iele]),muvec[1],transport_params.temperature,0);
                current += transport_params.conduct_quant*diff_fermi*real(stepvector[iele])*(-transmission[iele]);
            }
        }
        myfile.close();
        cout << "CURRENT IS " << current << endl;
    }
    if (!iam) {
        stringstream mysstream;
        mysstream << "CurrentFromTransmission_" << transport_params.cp2k_scf_iter;
        ofstream myfile(mysstream.str().c_str());
        myfile.precision(15);
        for (int ibias=-600;ibias<=600;ibias++) {
            double dbias=ibias*0.001;
            double current = 0.0;
            for (uint iele=0;iele<energyvector.size();iele++) {
                if (!imag(energyvector[iele])) {
                    double diff_fermi = fermi(real(energyvector[iele]),muvec[0]+dbias,transport_params.temperature,0)-fermi(real(energyvector[iele]),muvec[0],transport_params.temperature,0);
                    current += transport_params.conduct_quant*diff_fermi*real(stepvector[iele])*(-transmission[iele]);
                }
            }
            myfile << dbias << " " << current << endl;
        }
        myfile.close();
    }
    return 0;
}

int Energyvector::determine_energyvector(std::vector<CPX> &energyvector_cc,std::vector<CPX> &stepvector_cc,std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,std::vector< std::vector<int> > &propagating_sizes,cp2k_csr_interop_type KohnSham,cp2k_csr_interop_type Overlap,std::vector<double> &muvec,std::vector<contact_type> contactvec,transport_parameters transport_params)
{
    Singularities singularities(transport_params,contactvec);
    int determine_singularities = !(transport_params.real_int_method!=real_int_methods::GAUSSCHEBYSHEV && transport_params.cp2k_method==cp2k_methods::LOCAL_SCF);
    int propagating_from_bs = (transport_params.real_int_method==real_int_methods::GAUSSCHEBYSHEV);
    double bands_start=-1.0E6;
    if (determine_singularities) {
double sabtime=get_time(0.0);
        if ( singularities.Execute(KohnSham,Overlap) ) return (LOGCERR, EXIT_FAILURE);
        if (transport_params.cp2k_method!=cp2k_methods::LOCAL_SCF) for (uint i_mu=0;i_mu<muvec.size();i_mu++) muvec[i_mu]=singularities.determine_fermi(contactvec[i_mu].n_ele,i_mu);
        bands_start=singularities.energy_gs;
if (!iam) cout << "TIME FOR SINGULARITIES " << get_time(sabtime) << endl;
        int follow_bands = (transport_params.real_int_method==real_int_methods::GAUSSCHEBYSHEV);
        int debugout = 0;
        for (uint i_mu=0;i_mu<contactvec.size();i_mu++) singularities.write_bandstructure(i_mu,transport_params.cp2k_scf_iter,follow_bands,debugout);
    }
 
    double Temp=transport_params.temperature;
    double delta_eps_fermi=-log(transport_params.eps_fermi)*Temp;
    double muvec_min=*min_element(muvec.begin(),muvec.end());
    double muvec_max=*max_element(muvec.begin(),muvec.end());
    double muvec_avg=accumulate(muvec.begin(),muvec.end(),0.0)/muvec.size();
    double nonequi_start=muvec_min-delta_eps_fermi;
    double nonequi_end=muvec_max+delta_eps_fermi;

if (!iam) cout << "Fermi level difference " << muvec_max-muvec_min << endl;

// all localized states with lowest fermi level corresponding to occupation of localized states in bandgap
    if (transport_params.cp2k_method==cp2k_methods::LOCAL_SCF) {
//      double energy_vb=*max_element(singularities.energies_vb.begin(),singularities.energies_vb.end());
        double energy_cb=*min_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
        if (assign_real_axis_energies(energy_cb,nonequi_end,energyvector,stepvector,singularities.energies_extremum,muvec.size(),transport_params)) return (LOGCERR, EXIT_FAILURE);
    } else if (transport_params.cp2k_method==cp2k_methods::TRANSPORT) {
        if (assign_cmpx_cont_energies(bands_start,muvec_min,energyvector_cc,stepvector_cc,transport_params.temperature,transport_params.n_abscissae)) return (LOGCERR, EXIT_FAILURE);
        if (assign_real_axis_energies(nonequi_start,nonequi_end,energyvector,stepvector,singularities.energies_extremum,muvec.size(),transport_params)) return (LOGCERR, EXIT_FAILURE);
    } else if (transport_params.cp2k_method==cp2k_methods::TRANSMISSION) {
        if (!transport_params.extra_scf) {
            if (assign_cmpx_cont_energies(bands_start,muvec_avg,energyvector_cc,stepvector_cc,transport_params.temperature,transport_params.n_abscissae)) return (LOGCERR, EXIT_FAILURE);
        } else {
            if (assign_real_axis_energies(nonequi_start,nonequi_end,energyvector,stepvector,singularities.energies_extremum,muvec.size(),transport_params)) return (LOGCERR, EXIT_FAILURE);
        }
    } else return (LOGCERR, EXIT_FAILURE);
    if (!iam) {
        ofstream myfile("E_dat");
        myfile.precision(15);
        myfile << energyvector.size() << endl;
        for (uint iele=0;iele<energyvector.size();iele++)
            myfile << real(energyvector[iele]) << endl;
        myfile.close();
    }
// get propagating modes from bandstructure
    propagating_sizes.resize(energyvector.size());
    for (uint ie=0;ie<energyvector.size();ie++) propagating_sizes[ie].resize(contactvec.size(),-1);
    if (propagating_from_bs) {
double sabtime=get_time(0.0);
        if (!iam) {
            std::vector< std::vector< std::vector<double> > > propagating = singularities.get_propagating(energyvector);
            for (uint ie=0;ie<energyvector.size();ie++) {
                for (uint i_mu=0;i_mu<contactvec.size();i_mu++) {
                    propagating_sizes[ie][i_mu]=propagating[i_mu][ie].size();
                }
            }
        }
if (!iam) cout << "TIME FOR PROPAGATING MODES " << get_time(sabtime) << endl;
        for (uint ie=0;ie<energyvector.size();ie++) {
            MPI_Bcast(&propagating_sizes[ie][0],contactvec.size(),MPI_INT,0,MPI_COMM_WORLD);
        }
        if (!iam && contactvec.size()==muvec.size()+1) {
            stringstream mysstream;
            mysstream << "CurrentFromBandstructure_" << transport_params.cp2k_scf_iter;
            ofstream myfile(mysstream.str().c_str());
            myfile.precision(15);
            double cb_max = *max_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
            for (int ibias=-600;ibias<=600;ibias++) {
                double dbias=ibias*0.001;
                double current = 0.0;
                for (uint iele=0;iele<energyvector.size();iele++) {
                    if (!imag(energyvector[iele]) && real(energyvector[iele])>cb_max) {
                        double diff_fermi = fermi(real(energyvector[iele]),muvec[0]+dbias,transport_params.temperature,0)-fermi(real(energyvector[iele]),muvec[0],transport_params.temperature,0);
                        current += transport_params.conduct_quant*diff_fermi*real(stepvector[iele])*propagating_sizes[iele][contactvec.size()-1];
                    }
                }
                myfile << dbias << " " << current << endl;
            }
            myfile.close();
        }
        if (!iam && contactvec.size()==muvec.size()+1) {
            ofstream myfile("TransmissionFromBandstructure");
            myfile.precision(15);
            double cb_max = *max_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
            for (uint iele=0;iele<energyvector.size();iele++)
                if (!imag(energyvector[iele]))
                    myfile << real(energyvector[iele]) << " " << real(stepvector[iele]) << " " << (real(energyvector[iele])>cb_max)*propagating_sizes[iele][contactvec.size()-1] << endl;
            myfile.close();
        }
    }
    return 0;
}

int Energyvector::assign_real_axis_energies(double nonequi_start,double nonequi_end,std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,const std::vector< std::vector<double> > &energies_extremum,int muvec_size,transport_parameters transport_params)
{
    energyvector.clear();
    stepvector.clear();
    if (transport_params.real_int_method==real_int_methods::READFROMFILE) {
        ifstream evecfile("E.dat");
        if (evecfile.fail()) return (LOGCERR, EXIT_FAILURE);
        istream_iterator<double> start_evec(evecfile), end_evec;
        energyvector.assign(start_evec,end_evec);
        if (energyvector.size()==1) {
            stepvector.push_back(1.0);
        } else {
            stepvector.assign(1,(energyvector[1]-energyvector[0])/2.0);
            for (uint istep=1;istep<energyvector.size()-1;istep++) {
                stepvector.push_back((energyvector[istep+1]-energyvector[istep-1])/2.0);
            }
            stepvector.push_back((energyvector[energyvector.size()-1]-energyvector[energyvector.size()-2])/2.0);
        }
        evecfile.close();
    } else if (transport_params.real_int_method==real_int_methods::TRAPEZOIDAL) {
        for (int istep=0;istep<int(abs(nonequi_end-nonequi_start)/transport_params.energy_interval)+1;istep++) {
            energyvector.push_back(nonequi_start+istep*transport_params.energy_interval);
        }
        if (energyvector.size()==1) {
            stepvector.push_back(1.0);
        } else {
            stepvector.push_back((energyvector[1]-energyvector[0])/2.0);
            for (uint istep=1;istep<energyvector.size()-1;istep++) {
                stepvector.push_back((energyvector[istep+1]-energyvector[istep-1])/2.0);
            }
            stepvector.push_back((energyvector[energyvector.size()-1]-energyvector[energyvector.size()-2])/2.0);
        }
    } else if (transport_params.real_int_method==real_int_methods::GAUSSCHEBYSHEV) {
        std::vector<double> energylist;
        int n_energies;
        if (!iam) {
            energylist.push_back(nonequi_start);
            for (int i_mu=0;i_mu<muvec_size;i_mu++)
                for (uint i_energies=0;i_energies<energies_extremum[i_mu].size();i_energies++)
                    if (energies_extremum[i_mu][i_energies]>nonequi_start && energies_extremum[i_mu][i_energies]<nonequi_end)
                        energylist.push_back(energies_extremum[i_mu][i_energies]);
            energylist.push_back(nonequi_end);
            std::sort(energylist.begin(),energylist.end());
            n_energies=energylist.size();
        }
        MPI_Bcast(&n_energies,1,MPI_INT,0,MPI_COMM_WORLD);
        energylist.resize(n_energies);
        MPI_Bcast(&energylist[0],n_energies,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if (int(abs(nonequi_end-nonequi_start)/transport_params.energy_interval)+1<n_energies*transport_params.num_interval) return (LOGCERR, EXIT_FAILURE);
        double smallest_energy_distance=transport_params.min_interval;
        if (!iam) cout<<"Smallest enery distance "<<smallest_energy_distance<<endl;
        if (!iam) cout<<"Max number of points per small interval "<<transport_params.num_interval<<endl;
        if (!iam) cout<<"Average distance for big intervals "<<transport_params.energy_interval<<endl;
        if (!iam) cout<<"Singularities in range "<< n_energies-2 << endl;
        for (uint i_energies=1;i_energies<energylist.size();i_energies++) {
            int num_points_per_interval=max(transport_params.num_interval,int(ceil(abs(energylist[i_energies]-energylist[i_energies-1])/transport_params.energy_interval)));
            while ((energylist[i_energies]-energylist[i_energies-1])/2.0*(1.0-cos(M_PI/(2.0*num_points_per_interval)))<smallest_energy_distance && num_points_per_interval>1)
//          while ((energylist[i_energies]-energylist[i_energies-1])/2.0*(1.0-tanh(M_PI/2.0*sinh(3.0)))<smallest_energy_distance && num_points_per_interval>1)
                --num_points_per_interval;
            if (num_points_per_interval>1) {
                Quadrature quadrature(quadrature_types::GC,energylist[i_energies-1],energylist[i_energies],num_points_per_interval);
                energyvector.insert(energyvector.end(),quadrature.abscissae.begin(),quadrature.abscissae.end());
                stepvector.insert(stepvector.end(),quadrature.weights.begin(),quadrature.weights.end());
            }
        }
    } else return (LOGCERR, EXIT_FAILURE);
    return 0;
}

int Energyvector::assign_cmpx_cont_energies(double start,double mu,std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,double Temp,int num_points_on_contour)
{
    energyvector.clear();
    stepvector.clear();
    enum cc_int_method {do_pexsi,do_pole_summation,do_contour,do_line};
    cc_int_method method=do_pexsi;
    if (method==do_pexsi) {
        energyvector.resize(num_points_on_contour);
        stepvector.resize(num_points_on_contour);
        double EM=abs(mu-start); // Max|E-mu| for all Eigenvalues E of 
        if (PEXSI::GetPoleDensity(&energyvector[0],&stepvector[0],num_points_on_contour,Temp,0.0,EM,mu)) return (LOGCERR, EXIT_FAILURE);
        c_zscal(num_points_on_contour,CPX(M_PI/2.0,0.0),&stepvector[0],1);
    } else if (method==do_pole_summation) {
        double Temp_r=1.0*Temp;
        double Temp_i=1.0*Temp;
        double eps=-log((numeric_limits<double>::epsilon)())*Temp;
        double eps_r=-log((numeric_limits<double>::epsilon)())*Temp_r;
        double mu_r=start-eps_r;
        double eps_i=-log((numeric_limits<double>::epsilon)())*Temp_i;
        CPX    mu_i=CPX(mu_r-eps_r,eps_i);

        int nu;
        CPX zval;
        nu=0;
        while (imag(zval=CPX(mu,(2*nu+++1)*M_PI*Temp))<2.0*eps_i) {
            energyvector.push_back(zval);
            stepvector.push_back(-CPX(0.0,2.0*M_PI*Temp)*fermi(CPX(0.0,-1.0)*zval,CPX(0.0,-1.0)*mu_i,Temp_i,0));
        }

        nu=0;
        while (imag(zval=CPX(mu_r,(2*nu+++1)*M_PI*Temp_r))<2.0*eps_i) {
            energyvector.push_back(zval);
            stepvector.push_back(+CPX(0.0,2.0*M_PI*Temp_r)*fermi(CPX(0.0,-1.0)*zval,CPX(0.0,-1.0)*mu_i,Temp_i,0));
        }

        nu=0;
        while (real(zval=mu_i+CPX((2*nu+++1)*M_PI*Temp_i,0.0))<mu+eps) {
            energyvector.push_back(zval);
            stepvector.push_back(CPX(2.0*M_PI*Temp_i,0.0)*(fermi(zval,CPX(1.0,0.0)*mu,Temp,0)-fermi(zval,CPX(1.0,0.0)*mu_r,Temp_r,0)));
        }
    } else if (method==do_contour) {
        Quadrature quadrature(quadrature_types::CCGL,start,mu,num_points_on_contour);//mu is end here
        energyvector.insert(energyvector.end(),quadrature.abscissae.begin(),quadrature.abscissae.end());
        stepvector.insert(stepvector.end(),quadrature.weights.begin(),quadrature.weights.end());
    } else if (method==do_line) {
        Quadrature quadrature(quadrature_types::MR,start,mu,num_points_on_contour);//mu is end here
        energyvector.insert(energyvector.end(),quadrature.abscissae.begin(),quadrature.abscissae.end());
        stepvector.insert(stepvector.end(),quadrature.weights.begin(),quadrature.weights.end());
    }
    return 0;
}
