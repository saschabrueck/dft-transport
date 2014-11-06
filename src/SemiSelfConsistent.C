#include "FEMGrid.H"
#include "Poisson.H"
#include "GetSingularities.H"
#include "EnergyVector.H"
#include <numeric>
#include <iterator>

extern WireStructure *nanowire;
extern WireGenerator* Wire;
extern Poisson *OMEN_Poisson_Solver;
extern FEMGrid *FEM;
extern PARAM *parameter;
extern VOLTAGE *voltage;
 
int semiselfconsistent(TCSR<double> *Overlap,TCSR<double> *KohnSham,c_transport_type transport_params)
{
    int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    int procs;MPI_Comm_size(MPI_COMM_WORLD,&procs);

    c_dscal(KohnSham->n_nonzeros,transport_params.evoltfactor,KohnSham->nnz,1);

    if (transport_params.method==0) {
        TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
        TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
        if (!iam) {
            KohnShamCollect->write_CSR("KohnSham");
            OverlapCollect->write_CSR("Overlap");
            ofstream paramoutfile("TransportParams");
            paramoutfile << transport_params.method << endl;
            paramoutfile << transport_params.bandwidth << endl;
            paramoutfile << transport_params.n_cells << endl;
            paramoutfile << transport_params.n_occ << endl;
            paramoutfile << transport_params.n_abscissae << endl;
            paramoutfile << transport_params.n_kpoint << endl;
            paramoutfile << transport_params.extra_int_param1 << endl;
            paramoutfile << transport_params.extra_int_param2 << endl;
            paramoutfile << transport_params.extra_int_param3 << endl;
            paramoutfile << transport_params.evoltfactor << endl;
            paramoutfile << transport_params.colzero_threshold << endl;
            paramoutfile << transport_params.eps_limit << endl;
            paramoutfile << transport_params.eps_decay << endl;
            paramoutfile << transport_params.eps_singularities << endl;
            paramoutfile << transport_params.extra_param1 << endl;
            paramoutfile << transport_params.extra_param2 << endl;
            paramoutfile << transport_params.extra_param3 << endl;
            paramoutfile.close();
        }
        MPI_Barrier(MPI_COMM_WORLD);
        exit(0);
    }

    int do_semiself=0;
    if (transport_params.method==3) do_semiself=1;
    int NAtom_work=transport_params.extra_int_param3;
    if (do_semiself) NAtom_work=FEM->NAtom;

if (!iam) cout << "N_ATOMS " << NAtom_work << endl;

    double *Vbf = new double[Overlap->size_tot];
    double *rho_atom = new double[2*NAtom_work]();//ZERO IN THE SECOND COMPONENT
    double *rho_atom_previous = new double[2*NAtom_work];
    double *drho_atom_dV = new double[2*NAtom_work]();
    double *drho_atom_dV_previous = new double[2*NAtom_work];
    int *atom_of_bf = new int[Overlap->size_tot];
    for (int ibf=0;ibf<Overlap->size_tot;ibf++) {
        atom_of_bf[ibf]=ibf/(Overlap->size_tot/NAtom_work);
    }

    int n_mu=2;
    double *muvec = new double[n_mu];
    int *contactvec = new int[n_mu];
    contactvec[0]=1;
    contactvec[1]=2;

    if (!do_semiself) {
        Singularities singularities(transport_params,n_mu);
        if ( singularities.Execute(KohnSham,Overlap,contactvec) ) return (LOGCERR, EXIT_FAILURE);
        for (int i_mu=0;i_mu<n_mu;i_mu++) muvec[i_mu]=singularities.determine_fermi(0.0,i_mu);
        double mu_avg=std::accumulate(muvec,muvec+n_mu,0.0)/n_mu;
        muvec[0]=mu_avg;
        muvec[1]=mu_avg;
        Energyvector energyvector;
        if (energyvector.Execute(Overlap,KohnSham,n_mu,muvec,contactvec,rho_atom,drho_atom_dV,Vbf,NAtom_work,atom_of_bf,transport_params)) return (LOGCERR, EXIT_FAILURE);
if(!iam){
stringstream mysstream;
mysstream << "rhofile";
ofstream rhofile(mysstream.str().c_str());
for (int ig=0;ig<NAtom_work;ig++) rhofile<<rho_atom[ig]<<endl;
rhofile.close();
}
        return 0;
    }

    double Temp=transport_params.extra_param3;
    double *Vnew = new double[FEM->NGrid]();
    double *Vold = new double[FEM->NGrid];
    double *doping_atom = new double[FEM->NAtom];
    for (int ia=0;ia<FEM->NAtom;ia++) {
        doping_atom[ia]=FEM->doping[FEM->real_at_index[ia]];
    }
    double *doping_cell = new double[transport_params.n_cells]();
    for (int ia=0;ia<FEM->NAtom;ia++) {
        doping_cell[ia/(FEM->NAtom/transport_params.n_cells)]+=doping_atom[ia];
    }
    double *dopingvec = new double[n_mu];
    dopingvec[0]=doping_cell[0];
    dopingvec[1]=doping_cell[transport_params.n_cells-1];

if (!iam) cout << "DOPING " << dopingvec[0] << " " << dopingvec[1] << endl;

    double Vs=voltage->Vsmin;
    double Vd=voltage->Vdmin;
    double Vg=-voltage->Vgmin[0];
    double Vm;
    double dV;
    {
        Singularities singularities(transport_params,n_mu);
        if ( singularities.Execute(KohnSham,Overlap,contactvec) ) return (LOGCERR, EXIT_FAILURE);
        for (int i_mu=0;i_mu<n_mu;i_mu++) muvec[i_mu]=singularities.determine_fermi(0.0,i_mu);
        Vm=std::accumulate(muvec,muvec+n_mu,0.0)/n_mu;
        for (int i_mu=0;i_mu<n_mu;i_mu++) muvec[i_mu]=singularities.determine_fermi(dopingvec[i_mu],i_mu);
        double mu_avg=std::accumulate(muvec,muvec+n_mu,0.0)/n_mu;
        dV=(muvec[0]-muvec[1])/2.0;
        Vg+=mu_avg-*min_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
        muvec[0]=mu_avg-Vs;
        muvec[1]=mu_avg-Vd;
    }

if (!iam) cout << "FERMI LEVEL LEFT " << muvec[0] << " RIGHT " << muvec[1] << endl;
if (!iam) cout << "GATE POTENTIAL " << Vg << endl;

    std::vector<double> nuclearchargeperatom(FEM->NAtom,-(2.0*transport_params.n_occ)/FEM->NAtom);

    double Ls=nanowire->Ls;
    double Lc=nanowire->Lc;

    for (int IX=0;IX<FEM->NGrid;IX++) {
        double x=FEM->grid[3*IX]-FEM->grid[0];
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
        double bfac = 0.5/K_BOLTZMANN;
//        Vnew[IX] = (Vs+Vg+Vd)*fermi(x-(Ls+Lc),0.0,bfac,0)*fermi(-x+Ls,0.0,bfac,0)-Vs*fermi(x-(Ls+Lc),0.0,bfac,0)-Vd*fermi(-x+Ls,0.0,bfac,0);
    }
    c_dcopy(FEM->NGrid,Vnew,1,Vold,1);

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

    double residual=(numeric_limits<double>::max)();
    double density_old=(numeric_limits<double>::max)();
    double density_new;
    double density_criterion=parameter->poisson_criterion;
    int max_iter=parameter->poisson_iteration;
    for (int i_iter=1;i_iter<=max_iter;i_iter++) {

        double sabtime=get_time(0.0);
        c_dcopy(Overlap->n_nonzeros,Overlap_nnz,1,Overlap->nnz,1);
        for (int ibf=0;ibf<Overlap->size_tot;ibf++) {
            Vbf[ibf]=Vnew[FEM->real_at_index[atom_of_bf[ibf]]];
        }
        TCSR<double> *Pot = new TCSR<double>(Overlap,Vbf);
        TCSR<double> *KohnShamPot = new TCSR<double>(1.0,KohnSham,1.0,Pot);
        delete Pot;
        for (int ia=0;ia<FEM->NAtom;ia++) {
            Vbf[ia]=Vm+Vnew[FEM->real_at_index[ia]];
        }
        Energyvector energyvector;
        if (energyvector.Execute(Overlap,KohnShamPot,n_mu,muvec,contactvec,rho_atom,drho_atom_dV,Vbf,FEM->NAtom,atom_of_bf,transport_params)) return (LOGCERR, EXIT_FAILURE);
        delete KohnShamPot;

        if (!iam) cout << "TIME FOR SCHROEDINGER " << get_time(sabtime) << endl;

stringstream dosstream;
dosstream << "DOS_Profile" << i_iter;
rename("DOS_Profile",dosstream.str().c_str());
stringstream trastream;
trastream << "Transmission" << i_iter;
rename("Transmission",trastream.str().c_str());

///*
        sabtime=get_time(0.0);
        TCSR<double> *DensityCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
        double *rho_grid_tmp = new double[FEM->NGrid]();
        double* xyz_atoms = new double[3*FEM->NAtom];
        double* xyz_grid = new double[3*FEM->NGrid];
        c_dlacpy('A',3,FEM->NAtom,Wire->Layer_Matrix,7,xyz_atoms,3);
        c_dcopy(3*FEM->NGrid,FEM->grid,1,xyz_grid,1);
        int ngridproc=FEM->NGrid/procs+1;
        for (int ix=iam*ngridproc;ix<min((iam+1)*ngridproc,FEM->NGrid);ix++) {
            rho_grid_tmp[ix]=DensityCollect->collect_density(&xyz_grid[3*ix],xyz_atoms,atom_of_bf);
        }
        delete DensityCollect;
        double *rho_grid = new double[FEM->NGrid];
        MPI_Allreduce(rho_grid_tmp,rho_grid,FEM->NGrid,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if (!iam) cout << "TIME FOR COLLECT DENSITY " << get_time(sabtime) << endl;

if(!iam){
stringstream mysstream;
mysstream << "rhogridfile" << i_iter;
ofstream rgfile(mysstream.str().c_str());
for (int ig=0;ig<FEM->NGrid;ig++) rgfile<<rho_grid[ig]<<endl;
rgfile.close();
}
//*/

        if (transport_params.n_abscissae) c_daxpy(FEM->NAtom,1.0,&nuclearchargeperatom[0],1,rho_atom,1);

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
            c_dcopy(2*NAtom_work,rho_atom,1,rho_atom_previous,1);
            c_dcopy(2*NAtom_work,drho_atom_dV,1,drho_atom_dV_previous,1);
        }
        c_dscal(2*NAtom_work,1.0-mixing_parameter,rho_atom_previous,1);
        c_dscal(2*NAtom_work,1.0-mixing_parameter,drho_atom_dV_previous,1);
        c_daxpy(2*NAtom_work,mixing_parameter,rho_atom,1,rho_atom_previous,1);
        c_daxpy(2*NAtom_work,mixing_parameter,drho_atom_dV,1,drho_atom_dV_previous,1);

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

    return 0;
}
