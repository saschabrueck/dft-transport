#include "GetSingularities.H"
#include "EnergyVector.H"
#include <numeric>
#include <iterator>
#include <limits>
#include "SemiSelfConsistent.H"
 
int semiselfconsistent(cp2k_csr_interop_type S,cp2k_csr_interop_type KS,cp2k_csr_interop_type *P,cp2k_csr_interop_type *PImag,std::vector<double> muvec,std::vector<contact_type> contactvec,std::vector<int> Bsizes,std::vector<int> orb_per_atom,transport_parameters transport_params)
{
    int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    int procs;MPI_Comm_size(MPI_COMM_WORLD,&procs);

    int n_mu=muvec.size();
    if (n_mu != 2) return (LOGCERR, EXIT_FAILURE);

    if (FEM->NAtom != orb_per_atom.size()-1) return (LOGCERR, EXIT_FAILURE);

    double *doping_atom = new double[FEM->NAtom];
    for (int ia=0;ia<FEM->NAtom;ia++) {
        doping_atom[ia]=FEM->doping[FEM->real_at_index[ia]];
    }
    double dopingvec[2];
    for (int i_mu=0;i_mu<n_mu;i_mu++) {
        dopingvec[i_mu]=0;
        for (int ia=contactvec[i_mu].atomstart;ia<contactvec[i_mu].atomstart+contactvec[i_mu].natoms;ia++) {
            dopingvec[i_mu]+=doping_atom[ia];
        }
    }
    delete[] doping_atom;

    double *Vnew = new double[FEM->NGrid]();
    double *Vold = new double[FEM->NGrid];
    double *Vatom = new double[FEM->NAtom];
    double *rho_atom = new double[2*FEM->NAtom]();//ZERO IN THE SECOND COMPONENT
    double *rho_atom_previous = new double[2*FEM->NAtom];
    double *drho_atom_dV = new double[2*FEM->NAtom]();
    double *drho_atom_dV_previous = new double[2*FEM->NAtom];

    double Vs=voltage->Vsmin;
    double Vd=voltage->Vdmin;
    double Vg=-voltage->Vgmin[0];
    double Vm;
    double dV;
    ifstream fermilevelfile("OMEN_F");
    if (fermilevelfile) {
        for (int i_mu=0;i_mu<n_mu;i_mu++) fermilevelfile >> muvec[i_mu];
        Vm=std::accumulate(muvec.begin(),muvec.end(),0.0)/n_mu;
        for (int i_mu=0;i_mu<n_mu;i_mu++) fermilevelfile >> muvec[i_mu];
        double mu_avg=std::accumulate(muvec.begin(),muvec.end(),0.0)/n_mu;
        dV=(muvec[0]-muvec[1])/2.0;
        double min_cond_band_edge;
        fermilevelfile >> min_cond_band_edge;
        Vg+=mu_avg-min_cond_band_edge;
        muvec[0]=mu_avg-Vs;
        muvec[1]=mu_avg-Vd;
    } else {
        Singularities singularities(transport_params,contactvec);
        if ( singularities.Execute(KS,S) ) return (LOGCERR, EXIT_FAILURE);
        for (int i_mu=0;i_mu<n_mu;i_mu++) muvec[i_mu]=singularities.determine_fermi(contactvec[i_mu].n_ele,i_mu);
        Vm=std::accumulate(muvec.begin(),muvec.end(),0.0)/n_mu;
        muvec[0]=singularities.determine_fermi(contactvec[0].n_ele+dopingvec[0],0);
        muvec[1]=singularities.determine_fermi(contactvec[1].n_ele+dopingvec[1],1);
        double mu_avg=std::accumulate(muvec.begin(),muvec.end(),0.0)/n_mu;
        dV=(muvec[0]-muvec[1])/2.0;
        Vg+=mu_avg-*min_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
        muvec[0]=mu_avg-Vs;
        muvec[1]=mu_avg-Vd;
    }
    fermilevelfile.close();

if (!iam) cout << "FERMI LEVEL LEFT " << muvec[0] << " RIGHT " << muvec[1] << endl;
if (!iam) cout << "GATE POTENTIAL " << Vg << endl;

    double Ls=nanowire->Ls;
    double Lc=nanowire->Lc;

    enum choose_initial_guess {ramp,ferm,file};
    choose_initial_guess initial_guess=ramp;
    ifstream potgridinfile;
    if (initial_guess==file) potgridinfile.open("potgrid.input");
    for (int IX=0;IX<FEM->NGrid;IX++) {
        double x=FEM->grid[3*IX]-FEM->grid[0];
        if (initial_guess==ramp) {
            double dx_ramp = 1;
            if(x<(Ls-dx_ramp)){
                Vnew[IX] = -Vs-dV;
            }
            if((x>=(Ls-dx_ramp))&&(x<Ls)){
                Vnew[IX] = -Vs-dV+Vg*(x-(Ls-dx_ramp))/dx_ramp;
                //V[IX] = -Vs-dV+(Vg+Vs+dV)*(x-(Ls-dx_ramp))/dx_ramp;
            }
            if((x>=Ls)&&(x<(Ls+Lc))){
                Vnew[IX] = (Vs-Vd+2*dV)*(x-Ls)/Lc+Vg-Vs-dV;
                //V[IX] = Vg;
            }
            if((x>=(Ls+Lc))&&(x<(Ls+Lc+dx_ramp))){
                Vnew[IX] = -Vd+dV+Vg*(Ls+Lc+dx_ramp-x)/dx_ramp;
                //V[IX] = -Vd+dV+(Vg-dV+Vd)*(Ls+Lc+dx_ramp-x)/dx_ramp;
            }
            if(x>=(Ls+Lc+dx_ramp)){
                Vnew[IX] = -Vd+dV;
            }
        } else if (initial_guess==ferm) {
            double bfac = 0.5;
            Vnew[IX] = (Vs+Vg+Vd)*fermi(x-(Ls+Lc),0.0,bfac,0)*fermi(-x+Ls,0.0,bfac,0)-Vs*fermi(x-(Ls+Lc),0.0,bfac,0)-Vd*fermi(-x+Ls,0.0,bfac,0);
        } else if (initial_guess==file) {
            potgridinfile >> Vnew[IX];
        } else {
            return (LOGCERR, EXIT_FAILURE);
        }
    }
    if (initial_guess==file) potgridinfile.close();

if(!iam){
stringstream mysstream;
mysstream << "potgridfile" << 0;
ofstream potfile(mysstream.str().c_str());
for (int ig=0;ig<FEM->NGrid;ig++) potfile<<Vnew[ig]<<endl;
potfile.close();
}
if(!iam){
stringstream mysstream;
mysstream << "potfile" << 0;
ofstream rhofile(mysstream.str().c_str());
for (int ig=0;ig<FEM->NAtom;ig++) rhofile<<Vnew[FEM->real_at_index[ig]]<<endl;
rhofile.close();
}

    double residual=(numeric_limits<double>::max)();
    double density_old=(numeric_limits<double>::max)();
    double density_new;
    double density_criterion=parameter->poisson_criterion;
    int max_iter=parameter->poisson_iteration;
    for (int i_iter=1;i_iter<=max_iter;i_iter++) {

        for (int ia=0;ia<FEM->NAtom;ia++) {
            Vatom[ia]=Vnew[FEM->real_at_index[ia]];//add Vm here?
        }
        double sabtime=get_time(0.0);
        Energyvector energyvector;
        c_dscal(2*FEM->NAtom,0.0,rho_atom,1);
        if (energyvector.Execute(S,KS,P,PImag,muvec,contactvec,Bsizes,orb_per_atom,Vatom,rho_atom,transport_params)) return (LOGCERR, EXIT_FAILURE);
        if (!iam) cout << "TIME FOR SCHROEDINGER " << get_time(sabtime) << endl;

if(!iam){
stringstream mysstream;
mysstream << "rhofile" << i_iter;
ofstream rhofile(mysstream.str().c_str());
for (int ig=0;ig<FEM->NAtom;ig++) rhofile<<rho_atom[ig]<<endl;
rhofile.close();
}

        density_new=c_dasum(FEM->NAtom,rho_atom,1);
        if (abs(density_new-density_old)<density_criterion && residual<parameter->poisson_inner_criterion) {
            if (!iam) cout << "DENSITY CONVERGED" << endl;
            break;
        }
        density_old=density_new;

//        double mixing_parameter = 0.8;
        double mixing_parameter = 1.0;
        if (i_iter==1) {
            c_dcopy(2*FEM->NAtom,rho_atom,1,rho_atom_previous,1);
            c_dcopy(2*FEM->NAtom,drho_atom_dV,1,drho_atom_dV_previous,1);
        }
        c_dscal(2*FEM->NAtom,1.0-mixing_parameter,rho_atom_previous,1);
        c_dscal(2*FEM->NAtom,1.0-mixing_parameter,drho_atom_dV_previous,1);
        c_daxpy(2*FEM->NAtom,mixing_parameter,rho_atom,1,rho_atom_previous,1);
        c_daxpy(2*FEM->NAtom,mixing_parameter,drho_atom_dV,1,drho_atom_dV_previous,1);

        c_dcopy(FEM->NGrid,Vnew,1,Vold,1);
        double Temp=transport_params.temperature/transport_params.boltzmann_ev;
        OMEN_Poisson_Solver->solve(Vnew,Vold,rho_atom_previous,drho_atom_dV_previous,1,NULL,FEM,Wire,Temp,&Vg,Vs,Vs,Vd,&residual,parameter->poisson_inner_criterion,parameter->poisson_inner_iteration,1,1,MPI_COMM_SELF,MPI_COMM_WORLD,1,MPI_COMM_WORLD,0);
//        OMEN_Poisson_Solver->solve(Vnew,Vold,rho_atom,drho_atom_dV,1,NULL,FEM,Wire,Temp,&Vg,Vs,Vs,Vd,&residual,parameter->poisson_inner_criterion,parameter->poisson_inner_iteration,1,1,MPI_COMM_SELF,MPI_COMM_WORLD,1,MPI_COMM_WORLD,0);


if(!iam){
stringstream mysstream;
mysstream << "potgridfile" << i_iter;
ofstream potfile(mysstream.str().c_str());
for (int ig=0;ig<FEM->NGrid;ig++) potfile<<Vnew[ig]<<endl;
potfile.close();
}
if(!iam){
stringstream mysstream;
mysstream << "potfile" << i_iter;
ofstream rhofile(mysstream.str().c_str());
for (int ig=0;ig<FEM->NAtom;ig++) rhofile<<Vnew[FEM->real_at_index[ig]]<<endl;
rhofile.close();
}

    }

    delete[] Vnew;
    delete[] Vold;
    delete[] Vatom;
    delete[] rho_atom;
    delete[] rho_atom_previous;
    delete[] drho_atom_dV;
    delete[] drho_atom_dV_previous;

    return 0;
}
