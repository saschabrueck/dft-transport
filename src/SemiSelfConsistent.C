#include "CSR.H"
#include "FEMGrid.H"
#include "Poisson.H"
#include <vector>

extern WireStructure *nanowire;
extern WireGenerator* Wire;
extern Poisson *OMEN_Poisson_Solver;
extern FEMGrid *FEM;

int energyvector(TCSR<double>*,TCSR<double>*,int,double*,int*,double*,double*,double*,int,int*,c_transport_type);
        
int semiselfconsistent(TCSR<double> *Overlap,TCSR<double> *KohnSham,c_transport_type transport_params)
{
int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);

    int do_semiself=0;
    if (transport_params.method==3) do_semiself=1;
    int NAtom_work=transport_params.extra_int_param3;
    if (do_semiself) NAtom_work=FEM->NAtom;

if (!iam) cout << "N_ATOMS " << NAtom_work << endl;

    double Vs=0.0;
    double Vg=0.046;
    double Vd=0.0;
    double Temp=transport_params.extra_param3;
    double *Vbf = new double[Overlap->size_tot];
    double *rho_atom = new double[NAtom_work]();
    double *drho_atom_dV = new double[NAtom_work]();
    double *rho_atom_previous = new double[NAtom_work];
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

    if (!do_semiself) {
        int addpotential=0;
        if (addpotential) {
            dopingvec[0]=0.14;
            dopingvec[1]=0.14;
            int na=6;
            int nb=na+8;
            for (int ivvec=0;                                               ivvec<na*(Overlap->size_tot/transport_params.n_cells); ivvec++) Vbf[ivvec] = Vs;
            for (int ivvec=na*(Overlap->size_tot/transport_params.n_cells); ivvec<nb*(Overlap->size_tot/transport_params.n_cells); ivvec++) Vbf[ivvec] = Vg;
            for (int ivvec=nb*(Overlap->size_tot/transport_params.n_cells); ivvec<Overlap->size_tot; ivvec++)                               Vbf[ivvec] = Vd;
            TCSR<double> *Pot = new TCSR<double>(Overlap,Vbf);
            TCSR<double> *KohnShamPot = new TCSR<double>(1.0,KohnSham,1.0/transport_params.evoltfactor,Pot);
            if (energyvector(Overlap,KohnShamPot,n_mu,muvec,contactvec,dopingvec,rho_atom,drho_atom_dV,NAtom_work,atom_of_bf,transport_params)) return (LOGCERR, EXIT_FAILURE);
        } else {
            if (energyvector(Overlap,KohnSham,n_mu,muvec,contactvec,dopingvec,rho_atom,drho_atom_dV,NAtom_work,atom_of_bf,transport_params)) return (LOGCERR, EXIT_FAILURE);
        }
        return 0;
    }

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

    std::vector<double> nuclearchargeperatom(FEM->NAtom,-4.0);

    double Ls=nanowire->Ls;
    double Lg=nanowire->Lc;

    for (int ix=0;ix<FEM->NGrid;ix++) {
        double x=FEM->grid[3*ix]-FEM->grid[0];
//        if (x<=Ls) Vnew[ix]=Vs;
        if (x<=Ls-0.5) Vnew[ix]=Vs;
        else if (x>Ls-0.5 && x<=Ls) Vnew[ix]=Vs+Vg*(x-Ls+0.5)/0.5;
        else if (x>Ls && x<=Ls+Lg) Vnew[ix]=Vg;
        else if (x>Ls+Lg && x<=Ls+Lg+0.5) Vnew[ix]=Vg*(x-Ls-Lg-0.5)/(-0.5)+Vd;
        else Vnew[ix]=Vd;
    }

    double residual;
    double dens_tol = 1.0E-4;
    int i_iter;
    int max_iter=10;

    for (i_iter=1;i_iter<=max_iter;i_iter++) {

        for (int ibf=0;ibf<Overlap->size_tot;ibf++) {
            Vbf[ibf]=Vnew[FEM->real_at_index[atom_of_bf[ibf]]];
        }
        if (i_iter==0) {
            int na=6;
            int nb=na+8;
            for (int ivvec=0;                                               ivvec<na*(Overlap->size_tot/transport_params.n_cells); ivvec++) Vbf[ivvec] = Vs;
            for (int ivvec=na*(Overlap->size_tot/transport_params.n_cells); ivvec<nb*(Overlap->size_tot/transport_params.n_cells); ivvec++) Vbf[ivvec] = Vg;
            for (int ivvec=nb*(Overlap->size_tot/transport_params.n_cells); ivvec<Overlap->size_tot; ivvec++)                               Vbf[ivvec] = Vd;
        }
        TCSR<double> *Pot = new TCSR<double>(Overlap,Vbf);
        TCSR<double> *KohnShamPot = new TCSR<double>(1.0,KohnSham,1.0/transport_params.evoltfactor,Pot);
        delete Pot;

if (!iam) cout << "DOPING " << dopingvec[0] << " " << dopingvec[1] << endl;

        c_dcopy(FEM->NAtom,rho_atom,1,rho_atom_previous,1);
        if (energyvector(Overlap,KohnShamPot,n_mu,muvec,contactvec,dopingvec,rho_atom,drho_atom_dV,FEM->NAtom,atom_of_bf,transport_params)) return (LOGCERR, EXIT_FAILURE);
//for (int isab=0;isab<312;isab++) rho_atom[isab]=-2.004463;

        delete KohnShamPot;
        c_daxpy(FEM->NAtom,1.0,&nuclearchargeperatom[0],1,rho_atom,1);

        c_daxpy(FEM->NAtom,-1.0,rho_atom,1,rho_atom_previous,1);
        double condition = c_dasum(FEM->NAtom,rho_atom_previous,1)/FEM->NAtom;
        if (condition<dens_tol) break;

if(!iam){
ofstream rhofile("rhofile");
for (int ig=0;ig<FEM->NAtom;ig++) rhofile<<rho_atom[ig]<<endl;
rhofile.close();
}
        c_dscal(FEM->NGrid,0.0,Vold,1);
        MPI_Comm newcomm;
        int iam;
        MPI_Comm_rank(MPI_COMM_WORLD,&iam);
        MPI_Comm_split(MPI_COMM_WORLD,iam,iam,&newcomm);
        OMEN_Poisson_Solver->solve(Vnew,Vold,rho_atom,drho_atom_dV,1,NULL,FEM,Wire,Temp,&Vg,Vs,Vs,Vd,&residual,1.0E-4,10,1,1,newcomm,MPI_COMM_WORLD,1,MPI_COMM_WORLD,0);

if(!iam){
ofstream potfile("potfile");
for (int ig=0;ig<FEM->NGrid;ig++) potfile<<Vnew[ig]<<endl;
potfile.close();
}
        for (int ibf=0;ibf<Overlap->size_tot;ibf++)
            Vbf[ibf]=Vnew[FEM->real_at_index[atom_of_bf[ibf]]];
if(!iam){
ofstream vbffile("vbffile");
for (int ig=0;ig<Overlap->size_tot;ig++) vbffile<<Vbf[ig]<<endl;
vbffile.close();
}

//this is only for the condition

// !!!!!!!!!!! DONT USE THIS CONDITION AS THE POTENTIAL SHIFTS BUT THE DENSITY DOESNT SO USE THE DENSITY !!!!!!!!!!!!
if(!iam) cout<<"CONDITION "<<condition<<endl;
MPI_Barrier(MPI_COMM_WORLD);
//exit(0);
    }

    delete[] Vnew;
    delete[] Vold;
    delete[] Vbf;
    delete[] rho_atom;

    return 0;
}
