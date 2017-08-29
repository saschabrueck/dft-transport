/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "GetSingularities.H"
#include "EnergyVector.H"
#include <numeric>
#include <iterator>
#include <limits>
#include "SemiSelfConsistent.H"
 
int semiselfconsistent(cp2k_csr_interop_type S,cp2k_csr_interop_type KS,cp2k_csr_interop_type *P,cp2k_csr_interop_type *PImag,std::vector<double> muvec,std::vector<contact_type> contactvec,std::vector<int> Bsizes,std::vector<int> orb_per_atom,double mixing_parameter,transport_parameters transport_params)
{
    int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    int procs;MPI_Comm_size(MPI_COMM_WORLD,&procs);

    if (muvec.size() != 2) return (LOGCERR, EXIT_FAILURE);

    if (FEM->NAtom != int(orb_per_atom.size())-1) return (LOGCERR, EXIT_FAILURE);

    double dopingvec[2];
    for (uint i_mu=0;i_mu<muvec.size();i_mu++) {
        dopingvec[i_mu]=0;
        for (int ia=contactvec[i_mu].atomstart;ia<contactvec[i_mu].atomstart+contactvec[i_mu].natoms;ia++) {
            dopingvec[i_mu]+=FEM->doping[FEM->real_at_index[ia]];
        }
    }

    double *Vnew = new double[FEM->NGrid]();
    double *Vold = new double[FEM->NGrid];
    double *rho_atom = new double[2*FEM->NAtom]();//ZERO IN THE SECOND COMPONENT
    double *rho_atom_previous = new double[2*FEM->NAtom];

    int* atom_of_bf = new int[S.nrows_total];
    int atom=0;
    for (int i=0;i<S.nrows_total;i++) {
        if (i==orb_per_atom[atom+1]) ++atom;
        atom_of_bf[i]=atom;
    }

    double min_cond_band_edge;
    ifstream fermilevelfile("OMEN_F");
    if (fermilevelfile) {
        for (uint i_mu=0;i_mu<muvec.size();i_mu++) fermilevelfile >> muvec[i_mu];
        fermilevelfile >> min_cond_band_edge;
        fermilevelfile.close();
    } else {
        Singularities singularities(transport_params,contactvec);
        if ( singularities.Execute(KS,S) ) return (LOGCERR, EXIT_FAILURE);
        for (uint i_mu=0;i_mu<muvec.size();i_mu++) muvec[i_mu]=singularities.determine_fermi(contactvec[i_mu].n_ele+dopingvec[i_mu],i_mu);
        min_cond_band_edge=*min_element(singularities.energies_cb.begin(),singularities.energies_cb.end());
    }
    double dV=(muvec[0]-muvec[1])/2.0;
    double mu_avg=std::accumulate(muvec.begin(),muvec.end(),0.0)/muvec.size();
    double Vs=voltage->Vsmin;
    double Vd=voltage->Vdmin;
    double Vg=mu_avg-min_cond_band_edge-voltage->Vgmin[0]+nanowire->phi_m-nanowire->Xi_wire;
    muvec[0]=mu_avg-Vs;
    muvec[1]=mu_avg-Vd;

    enum choose_initial_guess {ramp,ferm,file};
    choose_initial_guess initial_guess=ramp;

    double Ls=nanowire->Ls;
    double Lc=nanowire->Lc;
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

        double sabtime=get_time(0.0);
        for(int i=0;i<KS.nrows_local;i++) {
            int atom_i=atom_of_bf[i+KS.first_row];
            for(int e=KS.rowptr_local[i]-1;e<KS.rowptr_local[i+1]-1;e++) {
                int atom_j=atom_of_bf[KS.colind_local[e]-1];
                KS.nzvals_local[e]+=(Vnew[FEM->real_at_index[atom_i]]+Vnew[FEM->real_at_index[atom_j]])/2.0/transport_params.evoltfactor*S.nzvals_local[e];
            }
        }
        c_dscal(2*FEM->NAtom,0.0,rho_atom,1);
        Energyvector energyvector;
        if (energyvector.Execute(S,KS,P,PImag,muvec,contactvec,Bsizes,orb_per_atom,rho_atom,transport_params)) return (LOGCERR, EXIT_FAILURE);
        for(int i=0;i<KS.nrows_local;i++) {
            int atom_i=atom_of_bf[i+KS.first_row];
            for(int e=KS.rowptr_local[i]-1;e<KS.rowptr_local[i+1]-1;e++) {
                int atom_j=atom_of_bf[KS.colind_local[e]-1];
                KS.nzvals_local[e]-=(Vnew[FEM->real_at_index[atom_i]]+Vnew[FEM->real_at_index[atom_j]])/2.0/transport_params.evoltfactor*S.nzvals_local[e];//add Vm here?
            }
        }
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

        if (i_iter==1) {
            c_dcopy(2*FEM->NAtom,rho_atom,1,rho_atom_previous,1);
        }
        c_dscal(2*FEM->NAtom,1.0-mixing_parameter,rho_atom_previous,1);
        c_daxpy(2*FEM->NAtom,mixing_parameter,rho_atom,1,rho_atom_previous,1);

        c_dcopy(FEM->NGrid,Vnew,1,Vold,1);
        double Temp=transport_params.temperature/transport_params.boltzmann_ev;
        OMEN_Poisson_Solver->solve(Vnew,Vold,rho_atom_previous,NULL,1,NULL,FEM,Wire,Temp,&Vg,Vs,Vs,Vd,&residual,parameter->poisson_inner_criterion,parameter->poisson_inner_iteration,procs,1,MPI_COMM_WORLD,MPI_COMM_WORLD,0,MPI_COMM_NULL,34);


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
    delete[] rho_atom;
    delete[] rho_atom_previous;

    delete[] atom_of_bf;

    return 0;
}
