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

int Energyvector::Execute(TCSR<double> *Overlap,TCSR<double> *KohnSham,std::vector<double> muvec,std::vector<contact_type> contactvec,transport_parameters *transport_params)
{
    double sabtime;
// allocate matrices to gather on every node
sabtime=get_time(0.0);
    MPI_Comm matrix_comm;
    MPI_Comm eq_rank_matrix_comm;
    MPI_Comm_split(MPI_COMM_WORLD,iam/transport_params->tasks_per_point,iam,&matrix_comm);
    MPI_Comm_split(MPI_COMM_WORLD,iam%transport_params->tasks_per_point,iam,&eq_rank_matrix_comm);
    TCSR<double> *KohnShamCollect = new TCSR<double>(eq_rank_matrix_comm,KohnSham);
    TCSR<double> *OverlapCollect  = new TCSR<double>(eq_rank_matrix_comm,Overlap);
    TCSR<double> *Ps = new TCSR<double>(OverlapCollect);
    Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
if (!iam) cout << "TIME FOR DISTRIBUTING MATRICES " << get_time(sabtime) << endl;
    std::vector<CPX> energyvector;
    std::vector<CPX> stepvector;
    std::vector<transport_methods::transport_method> methodvector;
    std::vector< std::vector<int> > propagating_sizes;
    if (determine_energyvector(energyvector,stepvector,methodvector,propagating_sizes,KohnSham,Overlap,muvec,contactvec,transport_params)) return (LOGCERR, EXIT_FAILURE);
    std::vector<double> currentvector(energyvector.size(),0.0);
    std::vector<double> transmission(energyvector.size(),0.0);
    std::vector<double> dos_profile(energyvector.size(),0.0);
    int matrix_id, n_mat_comm;
    MPI_Comm_size(eq_rank_matrix_comm,&n_mat_comm);
    MPI_Comm_rank(eq_rank_matrix_comm,&matrix_id);
sabtime=get_time(0.0);
    unsigned int jpos;
    for (int iseq=0;iseq<int(ceil(double(energyvector.size())/n_mat_comm));iseq++)
        if ( (jpos=matrix_id+iseq*n_mat_comm)<energyvector.size())
            if (abs(stepvector[jpos])>0.0)
                if (density(KohnShamCollect,OverlapCollect,Ps,energyvector[jpos],stepvector[jpos],methodvector[jpos],muvec,contactvec,currentvector[jpos],transmission[jpos],dos_profile[jpos],propagating_sizes[jpos],transport_params,jpos,matrix_comm))
                    return (LOGCERR, EXIT_FAILURE);
if (!iam) cout << "TIME FOR DENSITY " << get_time(sabtime) << endl;
    delete KohnShamCollect;
    delete OverlapCollect;
// trPS per energy point
    std::vector<double> currentvector2(energyvector.size(),0.0);
    std::vector<double> transmission2(energyvector.size(),0.0);
    std::vector<double> dos_profile2(energyvector.size(),0.0);
    MPI_Allreduce(&currentvector[0],&currentvector2[0],energyvector.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&transmission[0],&transmission2[0],energyvector.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&dos_profile[0],&dos_profile2[0],energyvector.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) {
        ofstream myfile("DOS_Profile");
        myfile.precision(15);
        for (uint iele=0;iele<energyvector.size();iele++)
            myfile << real(energyvector[iele]) << " " << imag(energyvector[iele]) << " " << real(stepvector[iele]) << " " << imag(stepvector[iele]) << " " << dos_profile2[iele] << endl;
        myfile.close();
    }
    if (!iam) {
        ofstream myfile("Transmission");
        myfile.precision(15);
        for (uint iele=0;iele<energyvector.size();iele++)
            if (!imag(energyvector[iele]))
                myfile << real(energyvector[iele]) << " " << real(stepvector[iele]) << " " << transmission2[iele] << endl;
        myfile.close();
    }
    if (!iam) {
        ofstream myfile("CurrentFromTransmission");
        myfile.precision(15);
        for (int ibias=1;ibias<=600;ibias++) {
            double current = 0.0;
            for (uint iele=0;iele<energyvector.size();iele++) {
                if (!imag(energyvector[iele])) {
                    double diff_fermi = fermi(real(energyvector[iele]),muvec[0]+ibias*0.001,transport_params->temperature,0)-fermi(real(energyvector[iele]),muvec[0],transport_params->temperature,0);
                    current += 2.0*E_ELECTRON*E_ELECTRON/(2.0*M_PI*H_BAR)*diff_fermi*real(stepvector[iele])*(-transmission2[iele]);
                }
            }
            myfile << ibias*0.001 << " " << current << endl;
        }
        myfile.close();
    }
    if (!iam) cout << "CURRENT IS " << c_ddot(energyvector.size(),&currentvector2[0],1,(double*)&stepvector[0],2) << endl;

    MPI_Reduce(Ps->nnz,Overlap->nnz,Ps->n_nonzeros,MPI_DOUBLE,MPI_SUM,0,eq_rank_matrix_comm);
    delete Ps;
    MPI_Comm_free(&matrix_comm);
    MPI_Comm_free(&eq_rank_matrix_comm);

    return 0;
}

int Energyvector::determine_energyvector(std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,std::vector<transport_methods::transport_method> &methodvector,std::vector< std::vector<int> > &propagating_sizes,TCSR<double> *KohnSham,TCSR<double> *Overlap,std::vector<double> &muvec,std::vector<contact_type> contactvec,transport_parameters *transport_params)
{
    Singularities singularities(transport_params,contactvec);
    int determine_singularities = !(transport_params->rlaxis_integration_method!=31 && transport_params->method==2);
    int propagating_from_bs = (determine_singularities);
    double bands_start;
    if (determine_singularities) {
double sabtime=get_time(0.0);
        if ( singularities.Execute(KohnSham,Overlap) ) return (LOGCERR, EXIT_FAILURE);
        if (transport_params->method!=2) for (uint i_mu=0;i_mu<muvec.size();i_mu++) muvec[i_mu]=singularities.determine_fermi(contactvec[i_mu].n_ele,i_mu);
        bands_start=singularities.energy_gs;
if (!iam) cout << "TIME FOR SINGULARITIES " << get_time(sabtime) << endl;
        int follow_bands = (transport_params->rlaxis_integration_method==31);
        for (uint i_mu=0;i_mu<contactvec.size();i_mu++) singularities.write_bandstructure(i_mu,follow_bands);
    }
 
    double Temp=transport_params->temperature;
    double delta_eps_fermi=-log(transport_params->eps_fermi)*K_BOLTZMANN*Temp;
    double muvec_min=*min_element(muvec.begin(),muvec.end());
    double muvec_max=*max_element(muvec.begin(),muvec.end());
    double muvec_avg=accumulate(muvec.begin(),muvec.end(),0.0)/muvec.size();
    double nonequi_start=muvec_min-delta_eps_fermi;
    double nonequi_end=muvec_max+delta_eps_fermi;

// all localized states with lowest fermi level corresponding to occupation of localized states in bandgap
    if (transport_params->method==2) {
        transport_params->n_abscissae=0;
//      double energy_vb=*max_element(singularities.energies_vb.begin(),singularities.energies_vb.end());
        double energy_cb=*min_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
        if (add_real_axis_energies(energy_cb,nonequi_end,energyvector,stepvector,methodvector,singularities.energies_extremum,transport_params)) return (LOGCERR, EXIT_FAILURE);
    } else if (transport_params->method==3) {
        if (add_cmpx_cont_energies(bands_start,muvec_min,energyvector,stepvector,methodvector,transport_params)) return (LOGCERR, EXIT_FAILURE);
        if (add_real_axis_energies(nonequi_start,nonequi_end,energyvector,stepvector,methodvector,singularities.energies_extremum,transport_params)) return (LOGCERR, EXIT_FAILURE);
    } else if (transport_params->method==4) {
        if (!transport_params->extra_scf) {
            if (add_cmpx_cont_energies(bands_start,muvec_avg,energyvector,stepvector,methodvector,transport_params)) return (LOGCERR, EXIT_FAILURE);
        } else {
            if (add_real_axis_energies(nonequi_start,nonequi_end,energyvector,stepvector,methodvector,singularities.energies_extremum,transport_params)) return (LOGCERR, EXIT_FAILURE);
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
        if (!iam && contactvec.size()>muvec.size()) {
            ofstream myfile("CurrentFromBandstructure");
            myfile.precision(15);
            double cb_max = *max_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
            for (int ibias=1;ibias<=12;ibias++) {
                double current = 0.0;
                for (uint iele=0;iele<energyvector.size();iele++) {
                    if (!imag(energyvector[iele]) && real(energyvector[iele])>cb_max) {
                        double diff_fermi = fermi(real(energyvector[iele]),muvec[0]+ibias*0.05,transport_params->temperature,0)-fermi(real(energyvector[iele]),muvec[0],transport_params->temperature,0);
                        current += 2.0*E_ELECTRON*E_ELECTRON/(2.0*M_PI*H_BAR)*diff_fermi*real(stepvector[iele])*propagating_sizes[iele][contactvec.size()-1];
                    }
                }
                myfile << ibias*0.05 << " " << current << endl;
            }
            myfile.close();
        }
        if (!iam && contactvec.size()>muvec.size()) {
            ofstream myfile("TransmissionFromBandstructure");
            myfile.precision(15);
            double cb_max = *max_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
            for (uint iele=0;iele<energyvector.size();iele++)
                if (!imag(energyvector[iele]))
                    myfile << real(energyvector[iele]) << " " << real(stepvector[iele]) << " " << (real(energyvector[iele])>cb_max)*propagating_sizes[iele][contactvec.size()-1] << endl;
            myfile.close();
        }
    }
    if (!iam) cout << "Size of Energyvector " << energyvector.size() << endl;
    return 0;
}

int Energyvector::read_real_axis_energies(std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,std::vector<transport_methods::transport_method> &methodvector)
{
    ifstream evecfile("OMEN_E");
    if (evecfile.fail()) return (LOGCERR, EXIT_FAILURE);
    energyvector.clear();
    stepvector.clear();
    methodvector.clear();
    istream_iterator<double> start_evec(evecfile), end_evec;
    energyvector.assign(start_evec,end_evec);
    methodvector.assign(energyvector.size(),transport_methods::WF);
    if (energyvector.size()==1) {
        stepvector.assign(1,CPX(1.0,0.0));
    } else {
        stepvector.assign(1,(energyvector[1]-energyvector[0])/2.0);
        for (uint istep=1;istep<energyvector.size()-1;istep++) {
            stepvector.push_back((energyvector[istep+1]-energyvector[istep-1])/2.0);
        }
        stepvector.push_back((energyvector[energyvector.size()-1]-energyvector[energyvector.size()-2])/2.0);
    }
    evecfile.close();
    return 0;
}

int Energyvector::add_real_axis_energies(double nonequi_start,double nonequi_end,std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,std::vector<transport_methods::transport_method> &methodvector,const std::vector< std::vector<double> > &energies_extremum,transport_parameters *transport_params)
{
    if (transport_params->rlaxis_integration_method==33) {
        if(read_real_axis_energies(energyvector,stepvector,methodvector)) return (LOGCERR, EXIT_FAILURE);
    } else {
        int num_trapez=int(abs(nonequi_end-nonequi_start)/transport_params->energy_interval)+1;
        if (transport_params->rlaxis_integration_method==32) {
            energyvector.clear();
            stepvector.clear();
            methodvector.clear();
            for (int istep=0;istep<num_trapez;istep++) {
                energyvector.push_back(nonequi_start+istep*transport_params->energy_interval);
            }
            methodvector.assign(energyvector.size(),transport_methods::WF);
            if (energyvector.size()==1) {
                stepvector.assign(1,CPX(1.0,0.0));
            } else {
                stepvector.assign(1,(energyvector[1]-energyvector[0])/2.0);
                for (uint istep=1;istep<energyvector.size()-1;istep++) {
                    stepvector.push_back((energyvector[istep+1]-energyvector[istep-1])/2.0);
                }
                stepvector.push_back((energyvector[energyvector.size()-1]-energyvector[energyvector.size()-2])/2.0);
            }
        } else {
            std::vector<double> energylist;
            int n_energies;
            if (!iam) {
                energylist.push_back(nonequi_start);
                for (int i_mu=0;i_mu<transport_params->num_contacts;i_mu++)
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
            if (num_trapez<n_energies*transport_params->num_interval) return (LOGCERR, EXIT_FAILURE);
            double smallest_energy_distance=transport_params->min_interval;
            if (!iam) cout<<"Smallest enery distance "<<smallest_energy_distance<<endl;
            if (!iam) cout<<"Max number of points per small interval "<<transport_params->num_interval<<endl;
            if (!iam) cout<<"Average distance for big intervals "<<transport_params->energy_interval<<endl;
            if (!iam) cout<<"Singularities in range "<< n_energies-2 << endl;
            for (uint i_energies=1;i_energies<energylist.size();i_energies++) {
                int num_points_per_interval=max(transport_params->num_interval,int(ceil(abs(energylist[i_energies]-energylist[i_energies-1])/transport_params->energy_interval)));
                while ((energylist[i_energies]-energylist[i_energies-1])/2.0*(1.0-cos(M_PI/(2.0*num_points_per_interval)))<smallest_energy_distance && num_points_per_interval>1)
//              while ((energylist[i_energies]-energylist[i_energies-1])/2.0*(1.0-tanh(M_PI/2.0*sinh(3.0)))<smallest_energy_distance && num_points_per_interval>1)
                    --num_points_per_interval;
                if (num_points_per_interval>1) {
                    Quadrature quadrature(quadrature_types::GC,energylist[i_energies-1],energylist[i_energies],num_points_per_interval);
                    energyvector.insert(energyvector.end(),quadrature.abscissae.begin(),quadrature.abscissae.end());
                    stepvector.insert(stepvector.end(),quadrature.weights.begin(),quadrature.weights.end());
                    methodvector.resize(energyvector.size(),transport_methods::WF);
                }
            }
        }
    }
    return 0;
}

int Energyvector::add_cmpx_cont_energies(double start,double mu,std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,std::vector<transport_methods::transport_method> &methodvector,transport_parameters *transport_params)
{
    double Temp=transport_params->temperature;
    int num_points_on_contour=transport_params->n_abscissae;
    enum choose_method {do_pexsi,do_pole_summation,do_contour,do_line};
    choose_method method=do_pexsi;
    if (method==do_pexsi) {
        energyvector.resize(num_points_on_contour);
        stepvector.resize(num_points_on_contour);
        double EM=abs(mu-start); // Max|E-mu| for all Eigenvalues E of 
        if (PEXSI::GetPoleDensity(&energyvector[0],&stepvector[0],num_points_on_contour,K_BOLTZMANN*Temp,0.0,EM,mu)) return (LOGCERR, EXIT_FAILURE);
        c_zscal(num_points_on_contour,CPX(M_PI/2.0,0.0),&stepvector[0],1);
    } else if (method==do_pole_summation) {
        double Temp_r=1.0*Temp;
        double Temp_i=1.0*Temp;
        double eps=-log((numeric_limits<double>::epsilon)())*K_BOLTZMANN*Temp;
        double eps_r=-log((numeric_limits<double>::epsilon)())*K_BOLTZMANN*Temp_r;
        double mu_r=start-eps_r;
        double eps_i=-log((numeric_limits<double>::epsilon)())*K_BOLTZMANN*Temp_i;
        CPX    mu_i=CPX(mu_r-eps_r,eps_i);

        int nu;
        CPX zval;
        nu=0;
        while (imag(zval=CPX(mu,(2*nu+++1)*M_PI*K_BOLTZMANN*Temp))<2.0*eps_i) {
            energyvector.push_back(zval);
            stepvector.push_back(-CPX(0.0,2.0*M_PI*K_BOLTZMANN*Temp)*fermi(CPX(0.0,-1.0)*zval,CPX(0.0,-1.0)*mu_i,Temp_i,0));
        }

        nu=0;
        while (imag(zval=CPX(mu_r,(2*nu+++1)*M_PI*K_BOLTZMANN*Temp_r))<2.0*eps_i) {
            energyvector.push_back(zval);
            stepvector.push_back(+CPX(0.0,2.0*M_PI*K_BOLTZMANN*Temp_r)*fermi(CPX(0.0,-1.0)*zval,CPX(0.0,-1.0)*mu_i,Temp_i,0));
        }

        nu=0;
        while (real(zval=mu_i+CPX((2*nu+++1)*M_PI*K_BOLTZMANN*Temp_i,0.0))<mu+eps) {
            energyvector.push_back(zval);
            stepvector.push_back(CPX(2.0*M_PI*K_BOLTZMANN*Temp_i,0.0)*(fermi(zval,CPX(1.0,0.0)*mu,Temp,0)-fermi(zval,CPX(1.0,0.0)*mu_r,Temp_r,0)));
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
    if (transport_params->method==4 && !transport_params->cutout) {
        methodvector.resize(energyvector.size(),transport_methods::EQ);
    } else {
        methodvector.resize(energyvector.size(),transport_methods::GF);
    }
    return 0;
}
