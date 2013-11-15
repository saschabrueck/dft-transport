#include "CSR.H"
#include "FEMGrid.H"
#include "Poisson.H"
#include <vector>

extern WireStructure *nanowire;
extern WireGenerator* Wire;
extern Poisson *OMEN_Poisson_Solver;
extern FEMGrid *FEM;

int energyvector(TCSR<double>*,TCSR<double>*,int,double*,int*,double*,double*,c_transport_type);
        
int semiselfconsistent(TCSR<double> *Overlap,TCSR<double> *KohnSham,c_transport_type transport_params)
{
    double Vs=0.0;
//    double Vg=1.5;
//    double Vd=0.5;
    double Vg=0.0;
    double Vd=0.0;
    double Ls=nanowire->Ls;
    double Lg=nanowire->Lc;
    double residual;
    double Temp=0.0;
    double *Vnew = new double[FEM->NGrid]();
    double *Vold = new double[FEM->NGrid];
    double *Vbf = new double[Overlap->size_tot];
    double *rho_atom = new double[FEM->NAtom];
    double *drho_atom_dV = new double[FEM->NAtom]();
    double condition_tol = 1.0E-2;
    double condition = 2*condition_tol;

    for (int ix=0;ix<FEM->NGrid;ix++) {
        double x=FEM->grid[3*ix]-FEM->grid[0];
        if (x<=Ls) Vnew[ix]=Vs;
        else if (x>Ls && x<=Ls+Lg) Vnew[ix]=Vg;
        else Vnew[ix]=Vd;
    }
/*
    int na=2;
    int nb=na+1;
    int nc=5;
    for (int ivvec=0;          ivvec<na*4*13*12; ivvec++) Vbf[ivvec] = Vs;
    for (int ivvec=na*4*13*12; ivvec<nb*4*13*12; ivvec++) Vbf[ivvec] = Vg;
    for (int ivvec=nb*4*13*12; ivvec<nc*4*13*12; ivvec++) Vbf[ivvec] = Vd;
*/
    std::vector<double> nuclearchargeperatom(FEM->NAtom,2.0);

    int n_mu=2;
    double *muvec = new double[n_mu];
    int *contactvec = new int[n_mu];
    contactvec[0]=1;
    contactvec[1]=2;
    double *dopingvec = new double[n_mu];
    double *doping_atom = new double[FEM->NAtom];
    for (int ia=0;ia<FEM->NAtom;ia++) {
        doping_atom[ia]=FEM->doping[FEM->real_at_index[ia]];
    }
    double *doping_cell = new double[6]();
    for (int ia=0;ia<FEM->NAtom;ia++) {
        doping_cell[ia/52]+=doping_atom[ia];
    }
    dopingvec[0]=doping_cell[0];
    dopingvec[1]=doping_cell[5];

    while (condition>condition_tol) {
        for (int ibf=0;ibf<Overlap->size_tot;ibf++)
            Vbf[ibf]=Vnew[FEM->real_at_index[ibf/12]];
        c_dcopy(FEM->NGrid,Vnew,1,Vold,1);//only needed for the condition
// add the potential 
        TCSR<double> *Pot = new TCSR<double>(Overlap,Vbf);
        TCSR<double> *KohnShamPot = new TCSR<double>(1.0,KohnSham,1.0/transport_params.evoltfactor,Pot);
        delete Pot;

int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);
if (!iam) cout << "DOPING " << dopingvec[0] << " " << dopingvec[1] << endl;

        if (energyvector(Overlap,KohnShamPot,n_mu,muvec,contactvec,dopingvec,rho_atom,transport_params)) return (LOGCERR, EXIT_FAILURE);
        delete KohnShamPot;
        c_daxpy(FEM->NAtom,1.0,&nuclearchargeperatom[0],1,rho_atom,1);
//zero for testing influence of bad mulliken wiggles, but gives weird results
//c_dcopy(FEM->NAtom,drho_atom_dV,1,rho_atom,1);

if(!iam){
ofstream rhofile("rhofile");
for (int ig=0;ig<FEM->NAtom;ig++) rhofile<<rho_atom[ig]<<endl;
rhofile.close();
}
        OMEN_Poisson_Solver->solve(Vnew,Vold,rho_atom,drho_atom_dV,1,NULL,FEM,Wire,Temp,Vg,Vs,Vs,Vd,&residual,1.0E-4,10,1,1,MPI_COMM_WORLD,MPI_COMM_WORLD,1,MPI_COMM_WORLD,0);

if(!iam){
ofstream potfile("potfile");
for (int ig=0;ig<FEM->NGrid;ig++) potfile<<Vnew[ig]<<endl;
potfile.close();
}
//this is only for the condition
        c_daxpy(FEM->NGrid,-1.0,Vnew,1,Vold,1);
        condition = c_dasum(FEM->NGrid,Vold,1)/c_dasum(FEM->NGrid,Vnew,1);
if(!iam) cout<<"CONDITION "<<condition<<endl;
//exit(0);
    }

    delete[] Vnew;
    delete[] Vold;
    delete[] Vbf;
    delete[] rho_atom;

    return 0;
}
