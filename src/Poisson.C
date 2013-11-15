/************************************************************************************************
For 1D case, computation of quasi Fermi levels and charge density was changed:
- n = N*exp(...)       : before, worked well, except at low temperature
- n = N*log(1+exp(...)): new treatment, was existing, but not used
************************************************************************************************/

#include "Poisson.H"

Poisson::Poisson()
{
    N3D           = 3;
    Eps0          = 8.854e-12;
    e             = 1.6022e-19;
    kB            = 1.38e-23;
    findx         = 1;
    del_PMatrix   = 0;
    del_P1DMatrix = 0;

    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
}

/************************************************************************************************/

Poisson::~Poisson()
{
    if(del_PMatrix){
        //delete P1;
        delete P2;
	delete RP2;
    }
    if(del_P1DMatrix){
        delete P;
    }
}

/************************************************************************************************/

void Poisson::init(WireGenerator *Wire,WireStructure *nanowire,FEMGrid *FEM,int max_proc_poisson,\
		   int poisson_solver,MPI_Comm vg_comm)
{
    if(!mpi_rank){
        printf("Initializing Poisson Equation\n");
    }
    
    int vg_size,vg_rank;

    Vol = Wire->Vol_atom;

    MPI_Comm_size(vg_comm,&vg_size);
    MPI_Comm_rank(vg_comm,&vg_rank);
    
    if(nanowire->NDim==1){

        if((vg_rank<max_proc_poisson)&&(poisson_solver)){
	    create_FDM_matrix(FEM);
	    del_P1DMatrix = 1;
	}

    }else{

        if((vg_rank<max_proc_poisson)&&(poisson_solver)){

	    switch(nanowire->update_fitness){
	    case 0:
	        create_FEM_matrix(FEM->poisson_index,FEM->NPoisson,FEM->tetrahedron,\
				  FEM->tetraVol,FEM->NTetra,FEM->grid,FEM->NGrid,FEM->Eps);
		break;
	    case 1:
	        P2 = new TCSR<double>(nanowire->fit_file);
		break;
	    case 2:
	        P2 = new TCSR<double>(nanowire->fit_file,min(vg_size,max_proc_poisson),vg_rank);
		break;
	    }
	    reduce_P(FEM->poisson_index,FEM->NGrid);
	
	    del_PMatrix = 1;

	}

	FEM->delete_tetrahedron();

    }
    
    //P1->write("P1.dat");
    //P2->write("P2.dat");
    //RP2->write("RP2.dat");
    
}

/************************************************************************************************/

void Poisson::solve(double *Vnew,double *Vold,double *rho_atom,double *drho_atom_dV,\
		    double sign_rho,FermiStructure *FS,FEMGrid* FEM,WireGenerator *Wire,\
		    double Temp,double Vg,double Vground,double Vsource,double Vdrain,\
                    double *residual,double poisson_criterion,int poisson_iteration,\
                    int max_proc_poisson,int poisson_solver,MPI_Comm loc_comm,\
		    MPI_Comm glob_comm,int bcast_size,MPI_Comm bcast_comm,int no_task)
{

    if(poisson_solver){

        if(Wire->nwire->NDim==1){
	    solve_FDM(Vnew,Vold,rho_atom,drho_atom_dV,sign_rho,FS,FEM,Wire,Temp,Vg,Vground, \
		      residual,poisson_criterion,poisson_iteration,max_proc_poisson, \
		      loc_comm,glob_comm);
	}else{
	    c_daxpy(Wire->No_Atom,1.0,&rho_atom[Wire->No_Atom],1,rho_atom,1);
	    solve_FEM(Vnew,Vold,rho_atom,drho_atom_dV,sign_rho,FEM,Wire,Temp,Vg,Vground,Vsource,\
                      Vdrain,residual, poisson_criterion,poisson_iteration,max_proc_poisson,\
                      loc_comm,glob_comm,bcast_size,bcast_comm,no_task);
	}
    }else{

        c_dcopy(FEM->NGrid,Vold,1,Vnew,1);
	*residual = 0.0;

    }

}

/************************************************************************************************/

void Poisson::solve_FDM(double *Vnew,double *Vold,double *rho_atom,double *drho_atom_dV,\
			double sign_rho,FermiStructure *FS,FEMGrid* FEM,WireGenerator *Wire,\
			double Temp,double Vg,double Vground,double *residual,\
			double poisson_criterion,int poisson_iteration,int max_proc_poisson,\
			MPI_Comm loc_comm,MPI_Comm glob_comm)
{

    int mpi_loc_size  = 1;
    int mpi_loc_rank  = 0;
    int mpi_glob_size = 1;
    int mpi_glob_rank = 0;

    MPI_Comm_size(glob_comm,&mpi_glob_size);
    MPI_Comm_rank(glob_comm,&mpi_glob_rank);

    if(mpi_glob_rank<max_proc_poisson){

        MPI_Comm_size(loc_comm,&mpi_loc_size);
        MPI_Comm_rank(loc_comm,&mpi_loc_rank);

	int IC            = 0;
	int max_iter      = poisson_iteration;
	int start         = P->first_row;
	int stop          = P->first_row+P->size;
	double t0         = 0;
	double criterion  = poisson_criterion/1000;
	double beta       = (e*1.0e-18)/Eps0;
	double NC         = fabs(max_sign_abs_vec(FEM->doping,FEM->NGrid,1));
	double *rho       = new double[2*P->size];
	double *drho_dV   = new double[P->size];
	double *doping    = new double[P->size];
	double *res       = new double[P->size];
	double *Ef        = new double[2*P->size];
	double *rhsign    = new double[2*P->size];
	double *delta_V   = new double[FEM->NGrid];
	double condition  = INF;
	double loc_cond   = INF;

	Fermi *FF         = new Fermi();
	LinearSolver<double> *solver;

	if(!mpi_loc_rank){
	    printf("Solving Poisson Equation\n");
	}
    
	t0 = get_time(0.0);
    
	init_var(rho,P->size);
	init_var(drho_dV,P->size);
	init_var(res,P->size);
	init_var(delta_V,FEM->NGrid);

	init_V(Vnew,Vold,FEM->NGrid,NULL,0,0.0,NULL,0,0.0,NULL,0,0.0,NULL,0,0.0);
	c_dcopy(P->size,&FEM->doping[start],1,doping,1);

	find_1D_fermi(Ef,rhsign,Vnew,rho,drho_dV,FF,FS,rho_atom,drho_atom_dV,Wire->No_Atom,\
		      NC,FEM->grid,start,stop,Temp);

	while((condition>criterion)&&(IC<max_iter)){
	    
	    calc_1D_cd_cdd(rho,drho_dV,rhsign,Ef,Vnew,FF,FS,NC,Temp,start,stop);

	    P->mat_vec_mult(Vnew,res,1);
	    c_daxpy(P->size,beta,rho,1,res,1);
	    c_daxpy(P->size,-beta,doping,1,res,1);
	    P->update_loc_diag(drho_dV,beta);
        
	    solver = new Umfpack<double>(P,loc_comm);
	    solver->prepare();
	    solver->solve_equation(delta_V,res,1);

	    delete solver;

	    iterate_V(Vnew,delta_V,NULL,FEM->NGrid);

	    loc_cond = fabs(max_sign_abs_vec(res,P->size,1)/(beta*NC));
	    if(mpi_loc_size==1){
	        condition = loc_cond;
	    }else{
	        MPI_Allreduce(&loc_cond,&condition,1,MPI_DOUBLE,MPI_MAX,loc_comm);
	    }
	    
	    if(!mpi_loc_rank){
	        printf("Poisson Inner Loop %i: %e\n",IC+1,condition);
	    }

	    P->update_loc_diag(drho_dV,-beta);
	    
	    //if(!IC){*residual = condition;}
	    *residual = condition;
        
	    IC++;
	}

	t0 = get_time(t0);

	if(!mpi_loc_rank){
	    printf("Solving Poisson Equation took %e [s]\n",t0);
	}

	delete[] rho;
	delete[] drho_dV;
	delete[] doping;
	delete[] delta_V;
	delete[] res;
	delete[] Ef;
	delete[] rhsign;
	delete FF;
    }
    
    if(mpi_glob_size>max_proc_poisson){
        MPI_Bcast(residual,1,MPI_DOUBLE,0,glob_comm);
        MPI_Bcast(Vnew,FEM->NGrid,MPI_DOUBLE,0,glob_comm);
    }

}

/************************************************************************************************/

void Poisson::solve_FEM(double *Vnew,double *Vold,double *rho_atom,double *drho_atom_dV,\
			double sign_rho,FEMGrid* FEM,WireGenerator *Wire,double Temp,double Vg,\
			double Vground,double Vsource,double Vdrain,double *residual,\
                        double poisson_criterion,int poisson_iteration,int max_proc_poisson,\
                        MPI_Comm loc_comm,MPI_Comm glob_comm,int bcast_size,MPI_Comm bcast_comm,\
                        int no_task)
{

    MPI_Status status;
    int mpi_loc_size  = 1;
    int mpi_loc_rank  = 0;
    int mpi_glob_size = 1;
    int mpi_glob_rank = 0;

    MPI_Comm_size(glob_comm,&mpi_glob_size);
    MPI_Comm_rank(glob_comm,&mpi_glob_rank);

    if(mpi_glob_rank<max_proc_poisson){

        MPI_Comm_size(loc_comm,&mpi_loc_size);
        MPI_Comm_rank(loc_comm,&mpi_loc_rank);

	int IC            = 0;
	int max_iter      = poisson_iteration;
	int start         = RP2->first_row;
	int stop          = RP2->first_row+RP2->size;
	double t0         = 0;
	double criterion  = poisson_criterion/1000;
	double beta1      = (e*1.0e9)/Eps0;
	double beta2      = (e*1.0e-18)/Eps0;
	double NC         = fabs(max_sign_abs_vec(FEM->doping,FEM->NGrid,1));
	double *rho       = new double[RP2->size];
	double *drho_dV   = new double[RP2->size];
	double *doping    = new double[RP2->size];
	double *res       = new double[RP2->size];
	double *delta_V   = new double[FEM->NPoisson];
	double *Ef        = new double[Wire->No_Atom];
	double *rhsign    = new double[Wire->No_Atom];
	double condition  = INF;
	double loc_cond   = INF;
	LinearSolver<double> *solver;

	if(!mpi_loc_rank){
	    printf("Solving Poisson Equation\n");
	}
    
	t0 = get_time(0.0);
    
	init_var(rho,RP2->size);
	init_var(drho_dV,RP2->size);
	init_var(res,RP2->size);
	init_var(delta_V,FEM->NPoisson);


//i do not touch the old potential at all, i dont need it
	init_V(Vnew,Vold,FEM->NGrid,FEM->gate_index,FEM->NGate,Vg,FEM->ground_index,\
	       FEM->NGround,Vground,FEM->source_index,FEM->NSource,Vsource,\
               FEM->drain_index,FEM->NDrain,Vdrain);

	init_doping(doping,FEM->doping,FEM->poisson_index,FEM->NGrid,start,stop);

	        copy_cd_cdd(rho,drho_dV,rho_atom,drho_atom_dV,FEM->ratom_index,\
			    Wire->ch_conv,FEM->poisson_index,FEM->atom_index,\
			    FEM->NAtom,start,stop);

	init_var(res,RP2->size);
//also i dont understand the minus sign for the doping
	    c_daxpy(RP2->size,-beta1,rho,1,res,1);
	    c_daxpy(RP2->size,beta1,doping,1,res,1);
      
		solver = new Umfpack<double>(RP2,loc_comm);
		solver->prepare();
	    solver->solve_equation(delta_V,res,1);

	    delete solver;
       
//why is the delta V subtracted? oh i change it to copy, but still dont understand what it used to be like 
            iterate_V(Vnew,delta_V,FEM->poisson_index,FEM->NGrid);
 
	t0 = get_time(t0);

	if(!mpi_loc_rank){
	    printf("Solving Poisson Equation took %e [s]\n",t0);
	}

	delete[] rho;
	delete[] drho_dV;
	delete[] doping;
	delete[] delta_V;
	delete[] res;
	delete[] Ef;
	delete[] rhsign;

    }
    
    if(mpi_glob_size>max_proc_poisson){
        if(bcast_size==1){
	    MPI_Bcast(residual,1,MPI_DOUBLE,0,glob_comm);
	    MPI_Bcast(Vnew,FEM->NGrid,MPI_DOUBLE,0,glob_comm);
	}else{
	    if(!mpi_glob_rank){
	        for(int i_proc=1;i_proc<bcast_size;i_proc++){
		    MPI_Send(residual,1,MPI_DOUBLE,i_proc*mpi_glob_size/bcast_size,\
			     0,glob_comm);
		    MPI_Send(Vnew,FEM->NGrid,MPI_DOUBLE,i_proc*mpi_glob_size/bcast_size,\
			     1,glob_comm);
		}
	    }else{
	        if(mpi_glob_rank%(mpi_glob_size/bcast_size)==0){
		    MPI_Recv(residual,1,MPI_DOUBLE,0,0,glob_comm,&status);
		    MPI_Recv(Vnew,FEM->NGrid,MPI_DOUBLE,0,1,glob_comm,&status);
		}
	    }
	    MPI_Bcast(residual,1,MPI_DOUBLE,0,bcast_comm);
	    MPI_Bcast(Vnew,FEM->NGrid,MPI_DOUBLE,0,bcast_comm);
	}
    }
 
}

/************************************************************************************************/

void Poisson::init_V(double *Vnew,double *Vold,int NGrid,int *gate_index,int NGate,double Vg,\
		     int *ground_index,int NGround,double Vground,int *source_index,int NSource,\
                     double Vsource,int *drain_index,int NDrain,double Vdrain)
{
    int IG,IS,ID;

//    c_dcopy(NGrid,Vold,1,Vnew,1);

    for(IG=0;IG<NGate;IG++){
        Vnew[gate_index[IG]] = Vg;
    }

    for(IG=0;IG<NGround;IG++){
        Vnew[ground_index[IG]] = Vground;
    }

    for(IS=0;IS<NSource;IS++){
        Vnew[source_index[IS]] = Vsource;
    }

    for(ID=0;ID<NDrain;ID++){
        Vnew[drain_index[ID]] = Vdrain;
    }
}

/************************************************************************************************/

void Poisson::iterate_V(double *Vnew,double *delta_V,int *poisson_index,int NGrid)
{
    int IA,NPoisson = 0;
    
    if(poisson_index!=NULL){
        for(IA=0;IA<NGrid;IA++){
	    if(poisson_index[IA]){
	        Vnew[IA] = delta_V[NPoisson];
		NPoisson++;
	    }
	}
    }else{
        for(IA=0;IA<NGrid;IA++){
	  Vnew[IA] = delta_V[IA];
	}
    }
}

/************************************************************************************************/

void Poisson::find_fermi(double *Ef,double *rhsign,double *Vnew,double *rho_atom,double NC,\
			 int *real_at_index,int NAtom,double Temp)
{
    int IA;
    double drho = 1.0e-8;
    
    for(IA=0;IA<NAtom;IA++){
        rhsign[IA] = rho_atom[IA]/fabs(rho_atom[IA]);
        Ef[IA]     = Vnew[real_at_index[IA]]+rhsign[IA]*kB*Temp/e*\
	  log(exp((rho_atom[IA]+rhsign[IA]*drho)/(NC*rhsign[IA]))-1);
    }
}

/************************************************************************************************/

void Poisson::find_1D_fermi(double *Ef,double *rhsign,double *Vnew,double *rho,double *drho_dV,\
			    Fermi *FF,FermiStructure *FS,double *rho_atom,double *drho_atom_dV,\
			    int NAtom,double NC,double *grid,int start,int stop,double Temp)
{
    int IA,IG;
    int NG           = stop-start;
    double drho      = 1.0e-8;
    double *rho_act  = new double[2*NG];

    for(IG=start;IG<stop;IG++){

        /*
        //rho[IG-start]       = 0.0;
        rho_act[IG-start]    = 0.0;
        rho_act[NG+IG-start] = 0.0;
        drho_dV[IG-start]    = 0.0;
	*/

        rho[IG-start]        = 0.0;
        rho[NG+IG-start]     = 0.0;
        drho_dV[IG-start]    = 0.0;

        for(IA=Round(grid[N3D*IG+1]);IA<Round(grid[N3D*IG+1]+grid[N3D*IG+2]);IA++){
          /*
            //rho[IG-start]      = rho[IG-start]+rho_atom[IA]+rho_atom[NAtom+IA];
            rho_act[IG-start]    = rho_act[IG-start]+rho_atom[IA];
            rho_act[NG+IG-start] = rho_act[NG+IG-start]+rho_atom[NAtom+IA];
            drho_dV[IG-start]    = drho_dV[IG-start]+drho_atom_dV[IA];
          */
            rho[IG-start]        = rho[IG-start]+rho_atom[IA];
            rho[NG+IG-start]     = rho[NG+IG-start]+rho_atom[NAtom+IA];
            drho_dV[IG-start]    = drho_dV[IG-start]+drho_atom_dV[IA];
        }

    }
       
    for(IG=start;IG<stop;IG++){

        /* 
        if(IG==start){
	    rho[IG-start]     = (rho_act[IG-start]+rho_act[IG+1-start])/2.0;
	    rho[NG+IG-start]  = (rho_act[NG+IG-start]+rho_act[NG+IG+1-start])/2.0; 
	}
	if(IG==stop-1){
	    rho[IG-start]     = (rho_act[IG-start]+rho_act[IG-1-start])/2.0; 
	    rho[NG+IG-start]  = (rho_act[NG+IG-start]+rho_act[NG+IG-1-start])/2.0; 
	}
	if((IG>start)&&(IG<stop-1)){
	    rho[IG-start]     = rho_act[IG-start]/2.0+\
	      (rho_act[IG-1-start]+rho_act[IG+1-start])/4.0;
	    rho[NG+IG-start]  = rho_act[NG+IG-start]/2.0+\
	      (rho_act[NG+IG-1-start]+rho_act[NG+IG+1-start])/4.0;
	}
	*/

	rho[IG-start]       = rho[IG-start]/(grid[N3D*IG+2]*Vol);
	rho[NG+IG-start]    = rho[NG+IG-start]/(grid[N3D*IG+2]*Vol);
	drho_dV[IG-start]   = drho_dV[IG-start]/(grid[N3D*IG+2]*Vol);

	if(!FS->derivative){

	    rho[IG-start]    = rho[IG-start]+rho[NG+IG-start];
	    rhsign[IG-start] = rho[IG-start]/fabs(rho[IG-start]);

	    rho[IG-start]    = rho[IG-start]+rhsign[IG-start]*drho/(grid[N3D*IG+2]*Vol);
            
	    Ef[IG-start]      = Vnew[IG]+rhsign[IG-start]*kB*Temp/e*\
	      log(exp(rho[IG-start]/(NC*rhsign[IG-start]))-1);

            /*
	    Ef[IG-start]      = Vnew[IG]+rhsign[IG-start]*kB*Temp/e*\
	      log(rho[IG-start]/(NC*rhsign[IG-start]));
            */
	}else{

	    if(rho[IG-start]!=0){
	        rhsign[IG-start]    = rho[IG-start]/fabs(rho[IG-start]);
	    }else{
	        rhsign[IG-start]    = 1.0;
	    }
	    if(rho[NG+IG-start]!=0){
	        rhsign[NG+IG-start] = rho[NG+IG-start]/fabs(rho[NG+IG-start]);
	    }else{
	        rhsign[NG+IG-start] = 1.0;
	    }

	    rho[IG-start]       = rho[IG-start]+rhsign[IG-start]*drho/(grid[N3D*IG+2]*Vol);
	    rho[NG+IG-start]    = rho[NG+IG-start]+rhsign[NG+IG-start]*drho/(grid[N3D*IG+2]*Vol);
	    
	    Ef[IG-start]        = FF->find_fermi(fabs(rho[IG-start]),rhsign[IG-start],\
						 FS->Ekl,FS->dkx,FS->Nkx,FS->dky,\
						 abs(FS->Nky),FS->dkz,abs(FS->Nkz),\
						 FS->n_of_modes,FS->spin_factor,\
						 FS->Temp,FS->cell_area);
	    Ef[IG-start]        = Vnew[IG]+Ef[IG-start];
	    
	    Ef[NG+IG-start]     = FF->find_fermi(fabs(rho[NG+IG-start]),rhsign[NG+IG-start],\
						 FS->Ekr,FS->dkx,FS->Nkx,FS->dky,\
						 abs(FS->Nky),FS->dkz,abs(FS->Nkz),\
						 FS->n_of_modes,FS->spin_factor,\
						 FS->Temp,FS->cell_area);
	    Ef[NG+IG-start]     = Vnew[IG]+Ef[NG+IG-start];
	
	    rho[IG-start]       = rho[IG-start]+rho[NG+IG-start];

	}
    }

    delete[] rho_act;

}

/************************************************************************************************/

void Poisson::init_doping(double *doping,double *full_doping,int *poisson_index,int NGrid,\
			  int start,int stop)
{
    int NPoisson,IA,index;

    index    = 0;
    NPoisson = 0;
    for(IA=0;IA<NGrid;IA++){
        if(poisson_index[IA]){
	    if((index>=start)&&(index<stop)){
	        doping[NPoisson] = full_doping[IA];
		NPoisson++;
	    }
	    index++;
        }
    }
}

/************************************************************************************************/

void Poisson::copy_cd_cdd(double *rho,double *drho_dV,double *rho_atom,double *drho_atom_dV,\
			  int *ratom_index,int *ch_conv,int *poisson_index,int *atom_index,\
			  int NAtom,int start,int stop)
{
    int IA;
    
    for(IA=0;IA<NAtom;IA++){
        if(poisson_index[atom_index[IA]]){
	    if((ratom_index[IA]>=start)&&(ratom_index[IA]<stop)){
	        rho[ratom_index[IA]-start]     = rho_atom[ch_conv[IA]];
	        drho_dV[ratom_index[IA]-start] = drho_atom_dV[ch_conv[IA]];
	    }
	}
    }

}

/************************************************************************************************/

void Poisson::calc_cd_cdd(double *rho,double *drho_dV,double *rhsign,double *Ef,double *Vnew,\
                          double NC,int *real_at_index,int *ratom_index,int *ch_conv,\
			  int *poisson_index,int *atom_index,int NAtom,double Temp,\
			  int start,int stop)
{
    int IA,cv;

    for(IA=0;IA<NAtom;IA++){
        if(poisson_index[atom_index[IA]]){
	    if((ratom_index[IA]>=start)&&(ratom_index[IA]<stop)){

	        cv                             = ch_conv[IA];
 
		rho[ratom_index[IA]-start]     = rhsign[cv]*NC*\
		  log(1+exp(rhsign[cv]*e*(Ef[cv]-Vnew[real_at_index[cv]])/(kB*Temp)));
		drho_dV[ratom_index[IA]-start] = -rhsign[cv]*e*rhsign[cv]*NC/(kB*Temp)/\
		  (exp(rhsign[cv]*e*(Vnew[real_at_index[cv]]-Ef[cv])/(kB*Temp))+1);
	    }
	}
    }
}

/************************************************************************************************/

void Poisson::calc_1D_cd_cdd(double *rho,double *drho_dV,double *rhsign,double *Ef,double *Vnew,\
			     Fermi *FF,FermiStructure *FS,double NC,double Temp,int start,\
			     int stop)
{
    int IA,cv;
    int NG = stop-start;
    double l_drho,r_drho;

    for(IA=start;IA<stop;IA++){

        cv          = IA-start;

	if(!FS->derivative){
                
	    rho[cv]     = rhsign[cv]*NC*\
	      log(1+exp(rhsign[cv]*e*(Ef[cv]-Vnew[IA])/(kB*Temp)));
	    drho_dV[cv] = -rhsign[cv]*e*rhsign[cv]*NC/(kB*Temp)/\
	      (exp(rhsign[cv]*e*(Vnew[IA]-Ef[cv])/(kB*Temp))+1);

            /*
	    rho[cv]     = rhsign[cv]*NC*\
	      exp(rhsign[cv]*e*(Ef[cv]-Vnew[IA])/(kB*Temp));
       
	    drho_dV[cv] = -rhsign[cv]*e*rhsign[cv]*NC/(kB*Temp)*\
	      exp(rhsign[cv]*e*(Ef[cv]-Vnew[IA])/(kB*Temp));
            */
	}else{

	    FF->density(&rho[cv],rhsign[cv],FS->Ekl,FS->dkx,FS->Nkx,FS->dky,abs(FS->Nky),\
			FS->dkz,abs(FS->Nkz),FS->n_of_modes,(Ef[cv]-Vnew[IA]),\
			(Ef[cv]-Vnew[IA]),1,FS->spin_factor,FS->Temp,\
			FS->cell_area);

	    FF->density(&rho[NG+cv],rhsign[NG+cv],FS->Ekr,FS->dkx,FS->Nkx,FS->dky,abs(FS->Nky),\
			FS->dkz,abs(FS->Nkz),FS->n_of_modes,(Ef[NG+cv]-Vnew[IA]),\
			(Ef[NG+cv]-Vnew[IA]),1,FS->spin_factor,FS->Temp,\
			FS->cell_area);

	    FF->derivate(&l_drho,rhsign[cv],FS->Ekl,FS->dkx,FS->Nkx,FS->dky,abs(FS->Nky),\
			 FS->dkz,abs(FS->Nkz),FS->n_of_modes,(Ef[cv]-Vnew[IA]),\
			 (Ef[cv]-Vnew[IA]),1,FS->spin_factor,FS->Temp,\
			 FS->cell_area);

	    FF->derivate(&r_drho,rhsign[NG+cv],FS->Ekr,FS->dkx,FS->Nkx,FS->dky,abs(FS->Nky),\
			 FS->dkz,abs(FS->Nkz),FS->n_of_modes,(Ef[NG+cv]-Vnew[IA]),\
			 (Ef[NG+cv]-Vnew[IA]),1,FS->spin_factor,FS->Temp,\
			 FS->cell_area);

	    rho[cv]     = rhsign[cv]*rho[cv]+rhsign[NG+cv]*rho[NG+cv];
	    drho_dV[cv] = -rhsign[cv]*l_drho-rhsign[NG+cv]*r_drho;

	}

    }
}

/************************************************************************************************/

void Poisson::create_FEM_matrix(int *poisson_index,int NPoisson,int *tetrahedron,\
                                double *tetraVol,int NTetra,double *grid,int NGrid,double *Eps)
{
    int IT,IF,IP,no_element;
    double Volume,Eps_act,A[4][3],h[4],p[4][3];//,c0ii,c0ij;
    //double *mp1 = new double[16*NTetra];
    double *mp2 = new double[16*NTetra];
    IJPOS  *ij  = new IJPOS[16*NTetra];
    double *mp1 = NULL;

    no_element = 0;
    
    for(IT=0;IT<NTetra;IT++){    
        Volume = tetraVol[IT];
        for(IF=0;IF<N3D+1;IF++){
            get_points(&p[0][0],&tetrahedron[(N3D+1)*IT],grid,IF);
            geometry(A[IF],&h[IF],p[0],p[1],p[2],p[3]);
        }
        Eps_act = get_permitivity(Eps,&tetrahedron[(N3D+1)*IT]);
        //c0ii    = Volume/10.0;
        //c0ij    = Volume/20.0;
        for(IF=0;IF<N3D+1;IF++){
            for(IP=IF;IP<N3D+1;IP++){
                ij[no_element].i         = tetrahedron[(N3D+1)*IT+IF];
                ij[no_element].j         = tetrahedron[(N3D+1)*IT+IP];
                ij[no_element].index     = no_element;
                if (IP == IF){
		    //mp1[no_element]      = c0ii;
                    mp2[no_element]      = Eps_act*c_dnrm2(N3D,A[IF],1)/(3.0*h[IF]);
                }else{
		    //mp1[no_element]      = c0ij;
                    mp2[no_element]      = Eps_act*c_ddot(N3D,A[IF],1,A[IP],1)/(9.0*Volume);
                    no_element++;
                    ij[no_element].i     = tetrahedron[(N3D+1)*IT+IP];
                    ij[no_element].j     = tetrahedron[(N3D+1)*IT+IF];
                    ij[no_element].index = no_element;
                    //mp1[no_element]      = mp1[no_element-1];
                    mp2[no_element]      = mp2[no_element-1];
                }
                no_element++;
            }
        }
    }
    assemble_matrix(no_element,mp1,mp2,ij,poisson_index,NPoisson,NGrid);
    
    //delete[] mp1;
    delete[] mp2;
    delete ij;

}

/************************************************************************************************/

void Poisson::assemble_matrix(int no_element,double *mp1,double *mp2,IJPOS *ij,\
                              int *poisson_index,int NPoisson,int NGrid)
{
    int I,J,II=0,N=0;
    int *index   = new int[no_element];
    int *IMIN    = new int[NGrid];
    int *IMAX    = new int[NGrid];
    int *index_i = new int[NPoisson];
    int *index_j = new int[no_element];
    //double *MP1  = new double[no_element];
    double *MP2  = new double[no_element];

    init_var(index_i,NPoisson);
    
    sort_ij(ij,index,IMIN,IMAX,no_element);

    for(I=0;I<NGrid;I++){
        if(poisson_index[I]){
            for(J=IMIN[I];J<=IMAX[I];J++){
                if(J>0){
                    if((ij[J].i!=ij[J-1].i)||(ij[J].j!=ij[J-1].j)){
		        //MP1[N]      = mp1[index[J]];
                        MP2[N]      = mp2[index[J]];
                        index_j[N]  = ij[J].j+findx;
                        index_i[II] = index_i[II]+1;
                        N++;
                    }else{
		        //MP1[N-1]    = MP1[N-1]+mp1[index[J]];
                        MP2[N-1]    = MP2[N-1]+mp2[index[J]];
                    }
                }else{
		    //MP1[N]          = mp1[index[J]];
                    MP2[N]          = mp2[index[J]];
                    index_j[N]      = ij[J].j+findx;
                    index_i[II]     = index_i[II]+1;
                    N++;
                }
            }
            II++;
        }
    }

    //P1 = new TCSR<double>(NPoisson,N,findx);
    P2 = new TCSR<double>(NPoisson,N,findx);

    //c_dcopy(N,MP1,1,P1->nnz,1);
    c_dcopy(N,MP2,1,P2->nnz,1);

    //icopy(NPoisson,index_i,P1->index_i);
    icopy(NPoisson,index_i,P2->index_i);

    //icopy(N,index_j,P1->index_j);
    icopy(N,index_j,P2->index_j);

    //P1->get_row_edge();
    P2->get_row_edge();
    
    delete[] index;
    delete[] IMIN;
    delete[] IMAX;
    delete[] index_i;
    delete[] index_j;
    //delete[] MP1;
    delete[] MP2;

}

/************************************************************************************************/

void Poisson::create_FDM_matrix(FEMGrid *FEM)
{

    double dp,dm,pEps,mEps,Eps,connection[3];
    int IB,IP,p_nnz,ind_j;
    P = new TCSR<double>(FEM->NGrid,3*FEM->NGrid,findx);

    init_var(P->index_i,FEM->NGrid);
    
    p_nnz = 0;
    for(IP=0;IP<FEM->NGrid;IP++){
        if(IP==0){
	    dp            = FEM->grid[N3D*(IP+1)]-FEM->grid[N3D*IP];
	    dm            = dp;
	    Eps           = FEM->Eps[IP];
	    pEps          = Eps;
	    mEps          = 0.0;
	    connection[1] = -(pEps+Eps)/(dp*(dp+dm));
	    connection[2] = -connection[1];
	}else{
	    if(IP==FEM->NGrid-1){
	        dm            = FEM->grid[N3D*IP]-FEM->grid[N3D*(IP-1)];
		dp            = dm;
		Eps           = FEM->Eps[IP];
		mEps          = Eps;
		connection[1] = -(mEps+Eps)/(dm*(dp+dm));
		connection[0] = -connection[1];
	    }else{
	        dp            = FEM->grid[N3D*(IP+1)]-FEM->grid[N3D*IP];
		dm            = FEM->grid[N3D*IP]-FEM->grid[N3D*(IP-1)];
		Eps           = FEM->Eps[IP];
		pEps          = FEM->Eps[IP+1];
		mEps          = FEM->Eps[IP-1];
		connection[0] = (mEps+Eps)/(dm*(dp+dm));
		connection[2] = (pEps+Eps)/(dp*(dp+dm));
		connection[1] = -connection[0]-connection[2];
	    }
	}

	for(IB=0;IB<3;IB++){
	    ind_j = IP+IB-1;
	    if((ind_j>=0)&&(ind_j<FEM->NGrid)){
	        P->index_i[IP]      = P->index_i[IP]+1;
		P->index_j[p_nnz]   = ind_j+findx;
	        P->nnz[p_nnz]       = connection[IB];
		if(ind_j==IP){
		    P->diag_pos[IP] = p_nnz;
		}
		p_nnz++;
	    }
	}
    }
    P->get_row_edge();
    P->n_nonzeros = p_nnz;

}

/************************************************************************************************/

void Poisson::reduce_P(int *poisson_index,int NGrid)
{
    int I,J,n_nnz=0,*j_shift;
    
    RP2     = new TCSR<double>(P2->size,P2->n_nonzeros,P2->findx);
    j_shift = new int[NGrid];

    RP2->size_tot  = P2->size_tot;
    RP2->first_row = P2->first_row;

    init_var(RP2->index_i,RP2->size);
    init_j_shift(j_shift,poisson_index,NGrid);

    for(I=0;I<P2->size;I++){
        for(J=(P2->edge_i[I]-P2->findx);J<(P2->edge_i[I+1]-P2->findx);J++){
            if(poisson_index[P2->index_j[J]-P2->findx]){
                RP2->nnz[n_nnz]      = P2->nnz[J];
                RP2->index_j[n_nnz]  = P2->index_j[J]+j_shift[P2->index_j[J]-P2->findx];
                RP2->index_i[I]      = RP2->index_i[I]+1;
                if((RP2->index_j[n_nnz]-RP2->findx)==(I+P2->first_row)){
                    RP2->diag_pos[I] = n_nnz;
                }
                n_nnz++;
            }
        }
    }
    RP2->get_row_edge();
    RP2->n_nonzeros = n_nnz;

    delete[] j_shift;

}

/************************************************************************************************/

void Poisson::init_j_shift(int *j_shift,int *poisson_index,int NGrid)
{
    int IP;
    
    if(!poisson_index[0]){
        j_shift[0] = -1;
    }else{
        j_shift[0] = 0;
    }

    for(IP=1;IP<NGrid;IP++){
        if(!poisson_index[IP]){
            j_shift[IP] = j_shift[IP-1]-1;
        }else{
            j_shift[IP] = j_shift[IP-1];
        }
    }
    
}

/************************************************************************************************/

double Poisson::get_permitivity(double *Eps,int *tetrahedron)
{
    int i;
    double Eps_act = 0.0;

    for(i=0;i<N3D+1;i++){
        Eps_act = Eps_act + Eps[tetrahedron[i]];
    }

    return (Eps_act)/(N3D+1.0);
}

/************************************************************************************************/

void Poisson::get_points(double *p,int *tetrahedron,double *grid,int i_pos)
{
    int i,index=0;

    for(i=0;i<i_pos;i++){
        c_dcopy(N3D,&grid[N3D*tetrahedron[i]],1,&p[N3D*index],1);
        index++;
    }
    
    for(i=i_pos+1;i<N3D+1;i++){
        c_dcopy(N3D,&grid[N3D*tetrahedron[i]],1,&p[N3D*index],1);
        index++;
    }

    c_dcopy(N3D,&grid[N3D*tetrahedron[i_pos]],1,&p[N3D*index],1);
}

/************************************************************************************************/

void Poisson::geometry(double *A,double *h,double *p1,double *p2,double *p3,double *p4)
{
    double v[3],v1[3],v2[3],v3[3];
    
    c_dcopy(N3D,p2,1,v1,1);
    c_daxpy(N3D,-1.0,p1,1,v1,1);
    c_dcopy(N3D,p3,1,v2,1);
    c_daxpy(N3D,-1.0,p1,1,v2,1);

    A[0] = (v1[1]*v2[2]-v1[2]*v2[1])/2.0;
    A[1] = (-v1[0]*v2[2]+v1[2]*v2[0])/2.0;
    A[2] = (v1[0]*v2[1]-v1[1]*v2[0])/2.0;

    c_dcopy(N3D,p4,1,v3,1);
    c_daxpy(N3D,-1.0,p1,1,v3,1);
    c_dcopy(N3D,A,1,v,1);
    c_dscal(N3D,1.0/c_dnrm2(N3D,v,1),v,1);
    *h = c_ddot(N3D,v3,1,v,1);

    if(*h<0){
        *h=-*h;
        c_dscal(N3D,-1.0,A,1);
    }
}

/************************************************************************************************/

void Poisson::sort_ij(IJPOS* ij, int *index, int *IMIN, int *IMAX, int no_element)
{
    
    int jpos0=0,jpos1=1,I,NI=0;
    
    sort(ij,ij+no_element,sorti);
    
    while(jpos1<no_element){

        while(ij[jpos1].i<=(ij[jpos0].i)){
            jpos1=jpos1+1;
            if (jpos1>=no_element){
                break;
            }
        }
        sort(&ij[jpos0],&ij[min(jpos1,no_element)],sortj);
        IMIN[NI]=jpos0;
        IMAX[NI]=jpos1-1;
 	NI++;

        if(jpos1<no_element){
            jpos0 = jpos1;
            jpos1 = jpos0+1;
        }
    }

    for(I=0;I<no_element;I++){
        index[I]=ij[I].index;
    }
    
}

/************************************************************************************************/

void Poisson::generate_output(const char *filename, double *data, int NR, int NC)
{
    int i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    for(i=0;i<NR;i++){
        for(j=0;j<NC;j++){
            myfile<<data[i+j*NR]<<" ";
        }
        myfile<<"\n";
    }
    myfile.close();
}

/************************************************************************************************/

void Poisson::generate_output(const char *filename, int *data, int NR, int NC)
{
    int i,j;
    
    ofstream myfile;
    myfile.open (filename);
    for(i=0;i<NR;i++){
        for(j=0;j<NC;j++){
            myfile<<data[i+j*NR]<<" ";
        }
        myfile<<"\n";
    }
    myfile.close();
}

/************************************************************************************************/


