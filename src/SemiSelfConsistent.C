#include "CSR.H"
#include "FEMGrid.H"
#include "Poisson.H"
#include "GetSingularities.H"
#include <vector>
#include <limits>

extern WireStructure *nanowire;
extern WireGenerator* Wire;
extern Poisson *OMEN_Poisson_Solver;
extern FEMGrid *FEM;
extern PARAM *parameter;
extern VOLTAGE *voltage;

int energyvector(TCSR<double>*,TCSR<double>*,int,double*,int*,double*,double*,double*,int,int*,c_transport_type);
        
int semiselfconsistent(TCSR<double> *Overlap,TCSR<double> *KohnSham,c_transport_type transport_params)
{
    int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    int procs;MPI_Comm_rank(MPI_COMM_WORLD,&procs);

    int do_semiself=0;
    if (transport_params.method==3) do_semiself=1;
    int NAtom_work=transport_params.extra_int_param3;
    if (do_semiself) NAtom_work=FEM->NAtom;

if (!iam) cout << "N_ATOMS " << NAtom_work << endl;

    double *Vbf = new double[Overlap->size_tot];
    double *rho_atom = new double[2*NAtom_work]();
    double *drho_atom_dV = new double[2*NAtom_work]();
    int *atom_of_bf = new int[Overlap->size_tot];
    for (int ibf=0;ibf<Overlap->size_tot;ibf++) {
        atom_of_bf[ibf]=ibf/(Overlap->size_tot/NAtom_work);
    }

    int n_mu=2;
    double *muvec = new double[n_mu];
    int *contactvec = new int[n_mu];
    contactvec[0]=1;
    contactvec[1]=2;
    double *dopingvec = new double[n_mu]();

    c_dscal(KohnSham->n_nonzeros,transport_params.evoltfactor,KohnSham->nnz,1);

    if (!do_semiself) {
        if (energyvector(Overlap,KohnSham,n_mu,muvec,contactvec,dopingvec,rho_atom,drho_atom_dV,NAtom_work,atom_of_bf,transport_params)) return (LOGCERR, EXIT_FAILURE);
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
    dopingvec[0]=doping_cell[0];
    dopingvec[1]=doping_cell[transport_params.n_cells-1];

if (!iam) cout << "DOPING " << dopingvec[0] << " " << dopingvec[1] << endl;

    double Vs=voltage->Vsmin;
    double Vd=voltage->Vdmin;
    double Vg=-voltage->Vgmin[0];
    {
        Singularities singularities(transport_params);
        if ( singularities.Execute(KohnSham,Overlap,n_mu,muvec,dopingvec,contactvec) ) return (LOGCERR, EXIT_FAILURE);
        Vg+=std::accumulate(muvec,muvec+n_mu,0.0)/n_mu-singularities.energy_cb;
    }

if (!iam) cout << "GATE POTENTIAL " << Vg << endl;

    std::vector<double> nuclearchargeperatom(FEM->NAtom,-4.0);

    double Ls=nanowire->Ls;
    double Lc=nanowire->Lc;

    for (int IX=0;IX<FEM->NGrid;IX++) {
        double x=FEM->grid[3*IX]-FEM->grid[0];
        double dx_ramp = 1.0;
        double dV= (muvec[0]-muvec[1])/2.0;

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
    }
    c_dcopy(FEM->NGrid,Vnew,1,Vold,1);

    double *Overlap_nnz = new double[Overlap->n_nonzeros];
    c_dcopy(Overlap->n_nonzeros,Overlap->nnz,1,Overlap_nnz,1);

    double residual=(numeric_limits<double>::max)();
    double density_old=(numeric_limits<double>::max)();
    double density_new;
    double density_criterion=parameter->poisson_criterion;
    int max_iter=parameter->poisson_iteration;
    for (int i_iter=1;i_iter<=max_iter;i_iter++) {

        c_dcopy(Overlap->n_nonzeros,Overlap_nnz,1,Overlap->nnz,1);
        for (int ibf=0;ibf<Overlap->size_tot;ibf++) {
            Vbf[ibf]=Vnew[FEM->real_at_index[atom_of_bf[ibf]]];
        }
        TCSR<double> *Pot = new TCSR<double>(Overlap,Vbf);
        for (int i_mu=0;i_mu<n_mu;i_mu++) {
            Pot->contactdensity(Pot->size_tot/transport_params.n_cells,transport_params.bandwidth,contactvec[i_mu],0,MPI_COMM_WORLD);
        }
        TCSR<double> *KohnShamPot = new TCSR<double>(1.0,KohnSham,1.0,Pot);
        delete Pot;

        if (energyvector(Overlap,KohnShamPot,n_mu,muvec,contactvec,dopingvec,rho_atom,drho_atom_dV,FEM->NAtom,atom_of_bf,transport_params)) return (LOGCERR, EXIT_FAILURE);
        delete KohnShamPot;

stringstream dosstream;
dosstream << "DOS_Profile" << i_iter;
rename("DOS_Profile",dosstream.str().c_str());
stringstream sinstream;
sinstream << "SingularityList" << i_iter;
rename("SingularityList",sinstream.str().c_str());

/*
        double *rho_grid_tmp = new double[FEM->NGrid]();
        double* xyz_atoms = new double[3*FEM->NAtom];
        double* xyz_grid = new double[3*FEM->NGrid];
        if (!iam) {
            c_dlacpy('A',3,FEM->NAtom,Wire->Layer_Matrix,7,xyz_atoms,3);
            c_dcopy(3*FEM->NGrid,FEM->grid,1,xyz_grid,1);
        }
        MPI_Bcast(xyz_atoms,3*FEM->NAtom,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(xyz_grid,3*FEM->NGrid,MPI_DOUBLE,0,MPI_COMM_WORLD);
        double sabtime=get_time(0.0);
        for (int ix=0;ix<FEM->NGrid;ix++) {
            Overlap->collect_density(rho_grid_tmp[ix],&xyz_grid[3*ix],xyz_atoms,atom_of_bf);
        }
        cout << "TIME FOR COLLECT DENSITY ON " << iam << " WITH " << Overlap->size << " ROWS IS " << get_time(sabtime) << endl;
        double *rho_grid = new double[FEM->NGrid];
        MPI_Allreduce(rho_grid_tmp,rho_grid,FEM->NGrid,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

if(!iam){
stringstream mysstream;
mysstream << "rhogridfile" << i_iter;
ofstream rgfile(mysstream.str().c_str());
for (int ig=0;ig<FEM->NGrid;ig++) rgfile<<rho_grid[ig]<<endl;
rgfile.close();
}
*/

//        c_daxpy(FEM->NAtom,1.0,&nuclearchargeperatom[0],1,rho_atom,1);

if(!iam){
stringstream mysstream;
mysstream << "rhofile" << i_iter;
ofstream rhofile(mysstream.str().c_str());
for (int ig=0;ig<FEM->NAtom;ig++) rhofile<<rho_atom[ig]<<endl;
rhofile.close();
}

        density_new=std::accumulate(rho_atom,rho_atom+FEM->NAtom,0.0);
        if (abs(density_new-density_old)<density_criterion && residual<parameter->poisson_inner_criterion) {
            if (!iam) cout << "DENSITY CONVERGED" << endl;
            return 0;
        }
        density_old=density_new;

        MPI_Comm newcomm;
        MPI_Comm_split(MPI_COMM_WORLD,iam,iam,&newcomm);
        OMEN_Poisson_Solver->solve(Vnew,Vold,rho_atom,drho_atom_dV,1,NULL,FEM,Wire,Temp,&Vg,Vs,Vs,Vd,&residual,parameter->poisson_inner_criterion,parameter->poisson_inner_iteration,1,1,newcomm,MPI_COMM_WORLD,1,MPI_COMM_WORLD,0);

if(!iam){
stringstream mysstream;
mysstream << "potfile" << i_iter;
ofstream potfile(mysstream.str().c_str());
for (int ig=0;ig<FEM->NGrid;ig++) potfile<<Vnew[ig]<<endl;
potfile.close();
}

    }

    return 0;
}
