#include "GetSingularities.H"
#include "EnergyVector.H"
#include <numeric>
#include <iterator>
#include <limits>
#include "SemiSelfConsistent.H"
 
int semiselfconsistent(TCSR<double> *Overlap,TCSR<double> *KohnSham,transport_parameters *transport_params)
{
    int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    int procs;MPI_Comm_size(MPI_COMM_WORLD,&procs);

//if (Overlap->additionalentries(KohnSham)) return (LOGCERR, EXIT_FAILURE);
//if (KohnSham->additionalentries(Overlap)) return (LOGCERR, EXIT_FAILURE);

    if ( Overlap->size_tot%transport_params->n_cells || transport_params->bandwidth<1 ) return (LOGCERR, EXIT_FAILURE);
    int n_mu=transport_params->num_contacts;
    double *muvec = new double[n_mu];
    contact_type *contactvec = new contact_type[n_mu];
    for (int i_mu=0;i_mu<n_mu;i_mu++) {
        contactvec[i_mu].bandwidth=transport_params->bandwidth;
        contactvec[i_mu].ndof=Overlap->size_tot/transport_params->n_cells; // ONLY IF ALL CELLS EQUAL
        contactvec[i_mu].n_occ=transport_params->n_occ/transport_params->n_cells; // THIS IS AN INTEGER DIVISION, IN GENERAL THE RESULT IS NOT CORRECT
    }
    contactvec[0].start=0;
    contactvec[0].inj_sign=+1;
    contactvec[1].start=Overlap->size_tot-contactvec[1].ndof*contactvec[1].bandwidth;
    contactvec[1].inj_sign=-1;

    if (FEM->NAtom != transport_params->n_atoms) return (LOGCERR, EXIT_FAILURE);

    double *Vbf = new double[Overlap->size_tot];
    double *rho_atom = new double[2*FEM->NAtom]();//ZERO IN THE SECOND COMPONENT
    double *rho_atom_previous = new double[2*FEM->NAtom];
    double *rho_dist = new double[FEM->NAtom];
    double *drho_atom_dV = new double[2*FEM->NAtom]();
    double *drho_atom_dV_previous = new double[2*FEM->NAtom];
    int *atom_of_bf = new int[Overlap->size_tot];
    for (int ibf=0;ibf<Overlap->size_tot;ibf++) {
        atom_of_bf[ibf]=ibf/(Overlap->size_tot/FEM->NAtom);
    }

    double Temp=transport_params->temperature;
    double *Vnew = new double[FEM->NGrid]();
    double *Vold = new double[FEM->NGrid];
    double *doping_atom = new double[FEM->NAtom];
    for (int ia=0;ia<FEM->NAtom;ia++) {
        doping_atom[ia]=FEM->doping[FEM->real_at_index[ia]];
    }
    double *doping_cell = new double[transport_params->n_cells]();
    for (int ia=0;ia<FEM->NAtom;ia++) {
        doping_cell[ia/(FEM->NAtom/transport_params->n_cells)]+=doping_atom[ia];
    }

if (!iam) cout << "DOPING " << doping_cell[0] << " " << doping_cell[transport_params->n_cells-1] << endl;

    double Vs=voltage->Vsmin;
    double Vd=voltage->Vdmin;
    double Vg=-voltage->Vgmin[0];
    double Vm;
    double dV;
    ifstream fermilevelfile("OMEN_F");
    if (fermilevelfile) {
        for (int i_mu=0;i_mu<n_mu;i_mu++) fermilevelfile >> muvec[i_mu];
        Vm=std::accumulate(muvec,muvec+n_mu,0.0)/n_mu;
        for (int i_mu=0;i_mu<n_mu;i_mu++) fermilevelfile >> muvec[i_mu];
        double mu_avg=std::accumulate(muvec,muvec+n_mu,0.0)/n_mu;
        dV=(muvec[0]-muvec[1])/2.0;
        double min_cond_band_edge;
        fermilevelfile >> min_cond_band_edge;
        Vg+=mu_avg-min_cond_band_edge;
        muvec[0]=mu_avg-Vs;
        muvec[1]=mu_avg-Vd;
    } else {
        Singularities singularities(transport_params,contactvec);
        if ( singularities.Execute(KohnSham,Overlap) ) return (LOGCERR, EXIT_FAILURE);
        for (int i_mu=0;i_mu<n_mu;i_mu++) muvec[i_mu]=singularities.determine_fermi(0.0,i_mu);
        Vm=std::accumulate(muvec,muvec+n_mu,0.0)/n_mu;
        muvec[0]=singularities.determine_fermi(doping_cell[0],0);
        muvec[1]=singularities.determine_fermi(doping_cell[transport_params->n_cells-1],1);
        double mu_avg=std::accumulate(muvec,muvec+n_mu,0.0)/n_mu;
        dV=(muvec[0]-muvec[1])/2.0;
        Vg+=mu_avg-*min_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
        muvec[0]=mu_avg-Vs;
        muvec[1]=mu_avg-Vd;
    }
    fermilevelfile.close();

if (!iam) cout << "FERMI LEVEL LEFT " << muvec[0] << " RIGHT " << muvec[1] << endl;
if (!iam) cout << "GATE POTENTIAL " << Vg << endl;

    std::vector<double> nuclearchargeperatom(FEM->NAtom,-(2.0*transport_params->n_occ)/FEM->NAtom);

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
            double bfac = 0.5/K_BOLTZMANN;
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

    double *Overlap_nnz = new double[Overlap->n_nonzeros];
    c_dcopy(Overlap->n_nonzeros,Overlap->nnz,1,Overlap_nnz,1);
    double *KohnSham_nnz = new double[KohnSham->n_nonzeros];
    c_dcopy(KohnSham->n_nonzeros,KohnSham->nnz,1,KohnSham_nnz,1);

    double residual=(numeric_limits<double>::max)();
    double density_old=(numeric_limits<double>::max)();
    double density_new;
    double density_criterion=parameter->poisson_criterion;
    int max_iter=parameter->poisson_iteration;
    for (int i_iter=1;i_iter<=max_iter;i_iter++) {

        double sabtime=get_time(0.0);
        c_dcopy(Overlap->n_nonzeros,Overlap_nnz,1,Overlap->nnz,1);
        c_dcopy(KohnSham->n_nonzeros,KohnSham_nnz,1,KohnSham->nnz,1);
        for (int ibf=0;ibf<Overlap->size_tot;ibf++) {
            Vbf[ibf]=Vnew[FEM->real_at_index[atom_of_bf[ibf]]];
        }
        KohnSham->add_pot(Overlap,Vbf);
        for (int ia=0;ia<FEM->NAtom;ia++) {
            Vbf[ia]=Vm+Vnew[FEM->real_at_index[ia]];
        }
        Energyvector energyvector;
        if (energyvector.Execute(Overlap,KohnSham,muvec,contactvec,transport_params)) return (LOGCERR, EXIT_FAILURE);

        if (!iam) cout << "TIME FOR SCHROEDINGER " << get_time(sabtime) << endl;

        for (int i=0;i<Overlap->n_nonzeros;i++) Overlap->nnz[i]*=Overlap_nnz[i];
        c_dscal(FEM->NAtom,0.0,rho_dist,1);
        Overlap->atom_allocate(atom_of_bf,rho_dist,2.0);
        MPI_Allreduce(rho_dist,rho_atom,FEM->NAtom,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if (!iam) {
            cout << "Number of electrons per unit cell";
            double e_total=0.0;
            for (int icell=0;icell<transport_params->n_cells;icell++) {
                double e_per_unit_cell=0.0;
                for (int iatom=icell*FEM->NAtom/transport_params->n_cells;iatom<(icell+1)*FEM->NAtom/transport_params->n_cells;iatom++) {
                    e_per_unit_cell+=rho_atom[iatom];
                    e_total+=rho_atom[iatom];
                }
                cout << " " << e_per_unit_cell;
            }
            cout << endl;
            cout << "Total number of electrons " << e_total << endl;
        }

stringstream dosstream;
dosstream << "DOS_Profile" << i_iter;
rename("DOS_Profile",dosstream.str().c_str());
stringstream trastream;
trastream << "Transmission" << i_iter;
rename("Transmission",trastream.str().c_str());

        if (transport_params->n_abscissae) c_daxpy(FEM->NAtom,1.0,&nuclearchargeperatom[0],1,rho_atom,1);

/*
double offsetminval=*min_element(rho_atom,rho_atom+FEM->NAtom);
for (int i=0;i<FEM->NAtom;i++) rho_atom[i]+=-offsetminval+0.1;
*/

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
            return 0;
        }
        density_old=density_new;

        MPI_Comm newcomm;
        MPI_Comm_split(MPI_COMM_WORLD,iam,iam,&newcomm);

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
        OMEN_Poisson_Solver->solve(Vnew,Vold,rho_atom_previous,drho_atom_dV_previous,1,NULL,FEM,Wire,Temp,&Vg,Vs,Vs,Vd,&residual,parameter->poisson_inner_criterion,parameter->poisson_inner_iteration,1,1,newcomm,MPI_COMM_WORLD,1,MPI_COMM_WORLD,0);
//        OMEN_Poisson_Solver->solve(Vnew,Vold,rho_atom,drho_atom_dV,1,NULL,FEM,Wire,Temp,&Vg,Vs,Vs,Vd,&residual,parameter->poisson_inner_criterion,parameter->poisson_inner_iteration,1,1,newcomm,MPI_COMM_WORLD,1,MPI_COMM_WORLD,0);


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

    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
}
