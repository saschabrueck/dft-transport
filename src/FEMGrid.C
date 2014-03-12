#include "FEMGrid.H"
extern "C" {
#include <qhull.h>
#include <qset.h>
}

FEMGrid::FEMGrid()
{
    N3D         = 3;

    del_all_var = 1;

    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
}

/************************************************************************************************/

FEMGrid::~FEMGrid()
{

    int ndg;

    delete[] grid;
    delete[] doping;
    delete[] Eps;
    delete[] real_at_index;
    delete[] NGate;
    if(del_all_var){
        delete[] atom_index;
	delete[] ratom_index;
	delete[] channel_index;
	delete[] poisson_index;
	delete[] ground_index;
        delete[] source_index;
        delete[] drain_index;

	for(ndg=0;ndg<no_diff_gate;ndg++){
	    delete[] gate_index[ndg];
	}
	delete[] gate_index;
    }

}

/************************************************************************************************/

void FEMGrid::delete_tetrahedron()
{
    delete[] tetrahedron;
    delete[] tetraVol;
}

/************************************************************************************************/

void FEMGrid::execute_task(WireGenerator *Wire,WireStructure *nanowire,int max_proc_poisson,\
			   int poisson_solver,MPI_Comm po_comm,MPI_Comm vg_comm)
{
    if(!mpi_rank){
        printf("Generating the %iD Grid\n",nanowire->NDim);
    }

    NB           = Wire->NB;
    no_diff_gate = nanowire->no_diff_gate;
    NGate        = new int[no_diff_gate];

    if(nanowire->NDim>1){
        make_grid(Wire,nanowire->grid_accuracy);
	delaunay(nanowire->update_tetra,nanowire->tetra_file,nanowire->update_fitness,\
		 max_proc_poisson,poisson_solver,po_comm,vg_comm);
	check_dimension(Wire,nanowire);   
	find_doping(nanowire,Wire->Vol_atom,NAtom);
	find_permittivity(nanowire);
	find_contact(nanowire);
	reduce_atom_index();
    }else{
        make_1D_grid(Wire);
	find_1D_doping(nanowire);
	find_1D_permittivity(nanowire);
	del_all_var = 0;
	NGate[0]    = 0;
        NGround     = 0;
        NSource     = 0;
        NDrain      = 0;
    }

    //write_grid("grid.dat");
    //write_tetrahedron("tetrahedron.dat");
    //write_gate("gate_index.dat");
    //write_permittivity("permittivity.dat");
 
}

/************************************************************************************************/

void FEMGrid::execute_reduced_task(WireGenerator *Wire,WireStructure *nanowire,int max_proc_poisson,\
				   MPI_Comm po_comm,MPI_Comm vg_comm)
{
    if(!mpi_rank){
        printf("Generating the %iD Grid\n",nanowire->NDim);
    }

    NB           = Wire->NB;
    no_diff_gate = nanowire->no_diff_gate;
    NGate        = new int[no_diff_gate];

    if(nanowire->NDim>1){
        make_grid(Wire,nanowire->grid_accuracy);
	check_dimension(Wire,nanowire);
	find_doping(nanowire,Wire->Vol_atom,NAtom);
	find_permittivity(nanowire);
	find_contact(nanowire);
	reduce_atom_index();
    }else{
        make_1D_grid(Wire);
	find_1D_doping(nanowire);
	find_1D_permittivity(nanowire);
	del_all_var = 0;
	NGate[0]    = 0;
        NGround     = 0;
        NSource     = 0;
        NDrain      = 0;
    }

}

/************************************************************************************************/

void FEMGrid::make_1D_grid(WireGenerator *Wire)
{

    int IA,IG,SLM;

    grid          = new double[N3D*Wire->No_Atom];
    real_at_index = new int[Wire->No_Atom];

    SLM           = Wire->SLM;
    grid[0]       = Wire->Layer_Matrix[SLM*Wire->ch_pos[0]];
    grid[1]       = 0;
    grid[2]       = 1;
    NGrid         = 1;

    for(IA=1;IA<Wire->No_Atom;IA++){
        if(abs(Wire->Layer_Matrix[SLM*Wire->ch_pos[IA]]-grid[N3D*(NGrid-1)])>10*tollim){
	    grid[N3D*NGrid]       = Wire->Layer_Matrix[SLM*Wire->ch_pos[IA]];
	    grid[N3D*NGrid+1]     = IA;
	    grid[N3D*NGrid+2]     = 1;
	    NGrid++;
	}else{
	    grid[N3D*(NGrid-1)+2] = grid[N3D*(NGrid-1)+2]+1;
	}
    }

    for(IG=0;IG<NGrid;IG++){
        for(IA=Round(grid[N3D*IG+1]);IA<Round(grid[N3D*IG+1]+grid[N3D*IG+2]);IA++){
	    real_at_index[IA] = IG;
	}
    }

}

/************************************************************************************************/

void FEMGrid::make_grid(WireGenerator *Wire,int grid_accuracy)
{
    int IA,IB,IC,NP,SLM;
    double *gpoints,vec[3];
    int Around_Atom = Wire->Around_tot;
    NAtom           = Wire->Channel_tot;
    NP              = NAtom+Around_Atom;
    SLM             = Wire->SLM;
    
    grid_accuracy   = max(grid_accuracy,0);
    gpoints         = new double[N3D*(Around_Atom+NAtom+NB*grid_accuracy*NAtom)];

    for(IA=0;IA<NAtom;IA++){
        c_dcopy(N3D,&Wire->Layer_Matrix[SLM*IA],1,&gpoints[N3D*IA],1);
        for(IB=0;IB<N3D;IB++){
            gpoints[N3D*IA+IB] = gpoints[N3D*IA+IB]-fmod(gpoints[N3D*IA+IB],tollim);
        }
    }
    c_dcopy(N3D*Around_Atom,Wire->Around_Matrix,1,&gpoints[N3D*NAtom],1);

    for(IA=0;IA<NAtom;IA++){
        for(IB=0;IB<NB;IB++){
            if(Wire->Layer_Matrix[SLM*IA+4+IB]>(IA+1)){
                c_dcopy(N3D,&Wire->Layer_Matrix[SLM*Round(Wire->Layer_Matrix[SLM*IA+4+IB]-1)],\
                        1,vec,1);
                c_daxpy(N3D,-1.0,&Wire->Layer_Matrix[SLM*IA],1,vec,1);
                c_dscal(N3D,1.0/((double)grid_accuracy+1.0),vec,1);
                for(IC=1;IC<grid_accuracy+1;IC++){
                    c_dcopy(N3D,&Wire->Layer_Matrix[SLM*IA],1,&gpoints[N3D*NP],1);
                    c_daxpy(N3D,(double)IC,vec,1,&gpoints[N3D*NP],1);
                    NP++;
                }
            }
        }
    }

    sort_grid_point(Wire,gpoints,NP,NAtom,Around_Atom,max_vec(Wire->Layer_Matrix,NAtom,SLM));
    
    delete[] gpoints;

}

/************************************************************************************************/

void FEMGrid::reduce_atom_index()
{
    int IA,IP;
    int *pos_shift = new int[NGrid];
    ratom_index    = new int[NAtom];
    
    if(!poisson_index[0]){
        pos_shift[0] = -1;
    }else{
        pos_shift[0] = 0;
    }

    for(IP=1;IP<NGrid;IP++){
        if(!poisson_index[IP]){
            pos_shift[IP] = pos_shift[IP-1]-1;
        }else{
            pos_shift[IP] = pos_shift[IP-1];
        }
    }

    for(IA=0;IA<NAtom;IA++){
        ratom_index[IA] = atom_index[IA]+pos_shift[atom_index[IA]];
    }

    delete[] pos_shift;

}

/************************************************************************************************/

void FEMGrid::sort_grid_point(WireGenerator *Wire,double *gpoints,int NP,int No_Atom,\
                              int Around_Atom,double xmax)
{
    int IA,IL,NLayer;
    XYZPOS* Grid_xyz       = new XYZPOS[NP];
    int *xyz_index         = new int[NP];
    int *Lmin              = new int[NP];
    int *Lmax              = new int[NP];

    for(IA=0;IA<NP;IA++){
        Grid_xyz[IA].x     = gpoints[N3D*IA];
        Grid_xyz[IA].y     = gpoints[N3D*IA+1];
        Grid_xyz[IA].z     = gpoints[N3D*IA+2];
        Grid_xyz[IA].index = IA;
    }
    sort_xyz(NP,Grid_xyz,xyz_index,&NLayer,Lmin,Lmax);

    grid                   = new double[N3D*NP];
    atom_index             = new int[No_Atom];
    real_at_index          = new int[Wire->No_Atom];
    channel_index          = new int[NP-Around_Atom];
    
    NGrid    = 0;
    NChannel = 0;
    for(IL=0;IL<NLayer;IL++){
        if(gpoints[N3D*xyz_index[Lmin[IL]]]<(xmax+tollim)){
            for(IA=Lmin[IL];IA<=Lmax[IL];IA++){
                c_dcopy(N3D,&gpoints[N3D*xyz_index[IA]],1,&grid[N3D*IA],1);
                if(xyz_index[IA]<No_Atom){
                    atom_index[xyz_index[IA]] = IA;
                    if(!Wire->index_channel[xyz_index[IA]]){
                        real_at_index[Wire->ch_conv[xyz_index[IA]]] = IA;
                    }
                }
                if((xyz_index[IA]<No_Atom)||(xyz_index[IA]>=(No_Atom+Around_Atom))){
                    channel_index[NChannel]   = IA;
                    NChannel++;
                }
                NGrid++;
            }
        }
    }

    delete Grid_xyz;
    delete[] xyz_index;
    delete[] Lmin;
    delete[] Lmax;

}

/************************************************************************************************/

void FEMGrid::find_contact(WireStructure *nanowire)
{

    int IP,IG,ndg;
    int condition1,condition2;
    int condition3,condition4;
    double xmin,xmax;
    double max_dist = 1.5*nanowire->a0/4;
    
    poisson_index   = new int[NGrid];
    ground_index    = new int[NGrid];
    source_index    = new int[NGrid];
    drain_index     = new int[NGrid];

    gate_index      = new int*[no_diff_gate];

    init_var(poisson_index,NGrid);

    for(IG=0;IG<nanowire->no_gate;IG++){
        if(!strcmp(nanowire->gate[IG]->type,"square")){
            nanowire->gate[IG]->face_area = calc_quad_area(nanowire->gate[IG]->p[0]->coord,\
                                                           nanowire->gate[IG]->p[1]->coord,\
                                                           nanowire->gate[IG]->p[2]->coord,\
                                                           nanowire->gate[IG]->p[3]->coord,N3D);
        }
    }

    for(IG=0;IG<nanowire->no_ground;IG++){
        if(!strcmp(nanowire->ground[IG]->type,"square")){
            nanowire->ground[IG]->face_area = calc_quad_area(nanowire->ground[IG]->p[0]->coord,\
							     nanowire->ground[IG]->p[1]->coord,\
							     nanowire->ground[IG]->p[2]->coord,\
							     nanowire->ground[IG]->p[3]->coord,N3D);
        }
    }
    
    NPoisson = 0;
    NGround  = 0;
    NSource  = 0;
    NDrain   = 0;

    for(ndg=0;ndg<no_diff_gate;ndg++){
        NGate[ndg]      = 0;
	gate_index[ndg] = new int[NGrid];
    }

    xmin     = min_vec(grid,NGrid,N3D);
    xmax     = max_vec(grid,NGrid,N3D);

    for(IP=0;IP<NGrid;IP++){
            
        condition1 = 0;
        for(IG=0;IG<nanowire->no_gate;IG++){
            if(is_in_gate(nanowire->gate[IG],&grid[N3D*IP],max_dist)){
                condition1        = 1;
		if(no_diff_gate>1){
		    gate_index[IG][NGate[IG]] = IP;
		    NGate[IG]++;
		}else{
		    gate_index[0][NGate[0]] = IP;
		    NGate[0]++;
		}
                break;
            }
        }
        
	condition2 = 0;
        for(IG=0;IG<nanowire->no_ground;IG++){
            if(is_in_gate(nanowire->ground[IG],&grid[N3D*IP],max_dist)){
                condition2            = 1;
                ground_index[NGround] = IP;
                NGround++;
                break;
            }
        }
        
        condition3 = 0;
        if(nanowire->Schottky->active){          
            if(abs(grid[N3D*IP]-xmin)<tollim){
                condition3            = 1;    
                source_index[NSource] = IP;
                NSource++;
            }
        }

        condition4 = 0;
        if(nanowire->Schottky->active){
            if(abs(grid[N3D*IP]-xmax)<tollim){
                condition4          = 1;    
                drain_index[NDrain] = IP;
                NDrain++;
            }
        }
                
        if((!condition1)&&(!condition2)&&(!condition3)&&(!condition4)){
            poisson_index[IP] = 1;
            NPoisson++;
        }
    }

}

/************************************************************************************************/

void FEMGrid::find_permittivity(WireStructure *nanowire)
{
    int IP,IM,ind,condition;
    double dVol_g,dVol_l;

    Eps = new double[NGrid];

    init_var(Eps,NGrid);

    if(nanowire->update_perm){
        FILE *F = fopen(nanowire->perm_file,"r");
	for(IP=0;IP<NGrid;IP++){
	    fscanf(F,"%lg",&Eps[IP]);
	}
	fclose(F);
    }else{
        for(IP=0;IP<NGrid;IP++){
	    for(IM=0;IM<nanowire->no_element;IM++){
	        condition = is_in_quad3D(nanowire->mat[IM],&grid[N3D*IP]);
		if(condition){
		    if(IM<nanowire->no_ch_element){
		        Eps[IP] = nanowire->Eps_wire[nanowire->mat[IM]->id_number-1];
		    }else{
		        Eps[IP] = nanowire->Eps_ox[nanowire->mat[IM]->id_number-1];
		    }
		    break;
		}
	    }
	    if(Eps[IP]==0){
	        dVol_g = INF;
		for(IM=0;IM<nanowire->no_element;IM++){
		    dVol_l = abs(check_volume(nanowire->mat[IM],&grid[N3D*IP])-\
				 nanowire->mat[IM]->volume);
		    if(dVol_l<dVol_g){
		        dVol_g = dVol_l;
			ind    = IM;
		    }
		}

		if(ind<nanowire->no_ch_element){
		    Eps[IP] = nanowire->Eps_wire[nanowire->mat[ind]->id_number-1];
		}else{
		    Eps[IP] = nanowire->Eps_ox[nanowire->mat[ind]->id_number-1];
		}
	    }
	}                    
    }
}

/************************************************************************************************/

void FEMGrid::find_1D_permittivity(WireStructure *nanowire)
{
    int IP,ind,IM;
    double dVol_g,dVol_l;

    Eps = new double[NGrid];

    init_var(Eps,NGrid);

    if(nanowire->update_perm){
        FILE *F = fopen(nanowire->perm_file,"r");
	for(IP=0;IP<NGrid;IP++){
	    fscanf(F,"%lg",&Eps[IP]);
	}
	fclose(F);
    }else{
        for(IP=0;IP<NGrid;IP++){
	    for(IM=0;IM<nanowire->no_element;IM++){
	        if((grid[N3D*IP]>=nanowire->mat[IM]->p[0]->coord[0]-tollim)&&\
		   (grid[N3D*IP]<=nanowire->mat[IM]->p[1]->coord[0]+tollim)){
		    Eps[IP] = nanowire->Eps_wire[nanowire->mat[IM]->id_number-1];
		    break;
		}
	    }
	    if(Eps[IP]==0){
	        dVol_g = INF;
		for(IM=0;IM<nanowire->no_element;IM++){
		    dVol_l = abs(check_volume(nanowire->mat[IM],&grid[N3D*IP])-\
				 nanowire->mat[IM]->volume);
		    if(dVol_l<dVol_g){
		        dVol_g = dVol_l;
			ind    = IM;
		    }
		}

		if(ind<nanowire->no_ch_element){
		    Eps[IP] = nanowire->Eps_wire[nanowire->mat[ind]->id_number-1];
		}else{
		    Eps[IP] = nanowire->Eps_ox[nanowire->mat[ind]->id_number-1];
		}
	    }
	}                
    }
}

/************************************************************************************************/

void FEMGrid::find_doping(WireStructure *nanowire,double Vol_atom,int No_Atom)
{
    FILE *F;
    int IP,ID;
    double Vol;

    doping = new double[NGrid];

    init_var(doping,NGrid);

    switch(nanowire->update_dop){

    case -1:

        calc_doping_info(nanowire);
        
	Vol = Vol_atom;

	for(IP=0;IP<NGrid;IP++){
	    for(ID=0;ID<nanowire->no_doping;ID++){
	        if(is_in_quad3D(nanowire->doping[ID],&grid[N3D*IP])){
		    doping[IP] = \
		      (nanowire->doping[ID]->ND-nanowire->doping[ID]->NA)*Vol;
		    break;
		}
	    }
	}
	break;
    case 0:

        calc_doping_info(nanowire);
        
	Vol = Vol_atom;
	
	for(IP=0;IP<No_Atom;IP++){
	    for(ID=0;ID<nanowire->no_doping;ID++){
	        if(is_in_quad3D(nanowire->doping[ID],&grid[N3D*atom_index[IP]])){
		    doping[atom_index[IP]] = \
		      (nanowire->doping[ID]->ND-nanowire->doping[ID]->NA)*Vol;
		    break;
		}
	    }
	}

	add_doping(nanowire,Vol_atom,No_Atom);

	break;
    case 1:

        F = fopen(nanowire->dop_file,"r");
        for(IP=0;IP<No_Atom;IP++){
	    fscanf(F,"%lg",&doping[atom_index[IP]]);
        }
        fclose(F);
	break;
    case 2:

        F = fopen(nanowire->dop_file,"r");
        for(IP=0;IP<NGrid;IP++){
	    fscanf(F,"%lg",&doping[IP]);
        }
        fclose(F);
	break;
    }

}

/************************************************************************************************/

void FEMGrid::find_1D_doping(WireStructure *nanowire)
{
    FILE *F;
    int IP,ID;

    doping = new double[NGrid];

    init_var(doping,NGrid);

    if(!nanowire->update_dop){

	for(IP=0;IP<NGrid;IP++){
	    for(ID=0;ID<nanowire->no_doping;ID++){
	        if((grid[N3D*IP]>=nanowire->doping[ID]->p[0]->coord[0]-tollim)&&\
		   (grid[N3D*IP]<=nanowire->doping[ID]->p[1]->coord[0]+tollim)){
		    doping[IP] = nanowire->doping[ID]->ND-nanowire->doping[ID]->NA;
		    break;
		}
	    }
	}
    }else{

        F = fopen(nanowire->dop_file,"r");
        for(IP=0;IP<NGrid;IP++){
	    fscanf(F,"%lg",&doping[IP]);
        }
        fclose(F);
    }

}

/************************************************************************************************/

void FEMGrid::add_doping(WireStructure *nanowire,double Vol_atom,int No_Atom)
{

    int IA,ID;
    int cond1,cond2,cond3;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    double x,y,z;
    double slope[3],conc;

    for(ID=0;ID<nanowire->no_doping_domain;ID++){

        xmin     = nanowire->doping_domain[ID]->xmin;
	xmax     = nanowire->doping_domain[ID]->xmax;
	ymin     = nanowire->doping_domain[ID]->ymin;
	ymax     = nanowire->doping_domain[ID]->ymax;
	zmin     = nanowire->doping_domain[ID]->zmin;
	zmax     = nanowire->doping_domain[ID]->zmax;

	slope[0] = nanowire->doping_domain[ID]->slope[0];
	slope[1] = nanowire->doping_domain[ID]->slope[1];
	slope[2] = nanowire->doping_domain[ID]->slope[2];

	conc     = nanowire->doping_domain[ID]->conc*Vol_atom;

        for(IA=0;IA<No_Atom;IA++){

	    x     = grid[N3D*atom_index[IA]];
	    y     = grid[N3D*atom_index[IA]+1];
	    z     = grid[N3D*atom_index[IA]+2];

	    cond1 = (x>=min(xmin,xmax)-tollim)&&(x<=max(xmax,xmin)+tollim);
	    cond2 = (y>=min(ymin,ymax)-tollim)&&(y<=max(ymax,ymin)+tollim);
	    cond3 = (z>=min(zmin,zmax)-tollim)&&(z<=max(zmax,zmin)+tollim);

	    if(cond1&&cond2&&cond3){
	        doping[atom_index[IA]] = conc*exp(-pow((x-xmin)/slope[0],2.0)*log(10.0))*\
		  exp(-pow((y-ymin)/slope[1],2.0)*log(10.0))*\
		  exp(-pow((z-zmin)/slope[2],2.0)*log(10.0));
	    }
	}
    }
}

/************************************************************************************************/

int FEMGrid::is_in_gate(GATE *gate,double *p0,double max_dist)
{
    int is_in_gate = 0;

    if(!strcmp(gate->type,"square")){

        double v[3],h[3],d;
        int cond1,cond2;
    
        calc_vec(v,gate->p[0]->coord,gate->p[1]->coord,gate->p[2]->coord);

        c_dcopy(N3D,p0,1,h,1);
        c_daxpy(N3D,-1.0,gate->p[0]->coord,1,h,1);

        d          = c_ddot(N3D,v,1,h,1);

        c_dcopy(N3D,p0,1,h,1);
        c_daxpy(N3D,-d,v,1,h,1);

        cond1      = is_in_quad2D(gate,h);
        cond2      = (fabs(d)<max_dist);

        is_in_gate = cond1&&cond2;
        
    }else{

        double r2,d2,phi;

	phi = atan(gate->radius[0]/gate->radius[1]*\
		   (p0[2]-gate->p[0]->coord[2])/(p0[1]-gate->p[0]->coord[1]+tollim));
        r2  = pow(p0[1]-gate->p[0]->coord[1],2.0)+pow(p0[2]-gate->p[0]->coord[2],2.0);
	d2  = pow(gate->radius[0]*cos(phi),2.0)+pow(gate->radius[1]*sin(phi),2.0);
       
        bool condition1 = p0[0]>=(gate->p[0]->coord[0]-tollim);
        bool condition2 = p0[0]<=(gate->p[1]->coord[0]+tollim);
        bool condition3 = fabs(sqrt(d2)-sqrt(r2))<max_dist;

	if(!strcmp(gate->type,"Omega")){

	    if((p0[1]-gate->p[0]->coord[1])<-tollim){
	        phi = phi+PI;
	    }

	    bool condition4 = (phi>=(2.0*PI*gate->angle[0]/360.0));
	    bool condition5 = (phi<=(2.0*PI*(gate->angle[1]/360.0)));

	    is_in_gate = condition1&&condition2&&condition3&&condition4&&condition5;

	}else{
	    is_in_gate = condition1&&condition2&&condition3;
	}
    }
    
    return is_in_gate;
}

/************************************************************************************************/

void FEMGrid::sort_xyz(int NP, XYZPOS* xyz, int *index, int* NL, int* LMI, int* LMA)
{
    
    int ypos0=0,ypos1=1,zpos0,zpos1,zlimit,IA;
    double coord_lim = 1e-6;

    NL[0]=0;
    sort(xyz,xyz+NP,sortx);
    
    while(ypos1<NP){

        while(xyz[ypos1].x<=(xyz[ypos0].x+coord_lim)){
            ypos1 = ypos1+1;
            if (ypos1>=NP){
                break;
            }
        }
        sort(&xyz[ypos0],&xyz[min(ypos1,NP)],sorty);
        LMI[*NL] = ypos0;
        LMA[*NL] = ypos1-1;
	*NL      = *NL+1;

        zpos0  = ypos0;
        zpos1  = ypos0+1;
	zlimit = min(ypos1,NP);
        while(zpos1<zlimit){
            while(xyz[zpos1].y<=(xyz[zpos0].y+coord_lim)){
                zpos1 = zpos1+1;
                if (zpos1>=zlimit){
                    break;
                }
            }
            sort(&xyz[zpos0],&xyz[min(zpos1,zlimit)],sortz);
            if(zpos1<zlimit){
                zpos0 = zpos1;
                zpos1 = zpos0+1;
            }
        }
        
        if(ypos1<NP){
            ypos0 = ypos1;
            ypos1 = ypos0+1;
        }
	
	if((ypos1==NP)&&(LMA[*NL-1]<NP-1)){
	//if((ypos1==NP)&&(LMA[*NL-1]<NP)){
	    LMI[*NL] = LMA[*NL-1]+1;
	    LMA[*NL] = NP-1;
	    *NL      = *NL+1;
	}
	
    }

    for(IA=0;IA<NP;IA++){
        index[IA] = xyz[IA].index;
    }
    
}

/************************************************************************************************/

void FEMGrid::delaunay(int update_tetra,char* tetra_file,int update_fitness,int max_proc_poisson,\
		       int poisson_solver,MPI_Comm po_comm,MPI_Comm vg_comm)
{

    int vg_rank;

    MPI_Comm_rank(vg_comm,&vg_rank);
    
    if(!vg_rank){
      
        if((!update_tetra)&&(!update_fitness)&&(poisson_solver)){

	    char dataset_name[255];
	    sprintf(dataset_name,"qhull_%i.dat",mpi_rank);
	    boolT ismalloc = 0;
	    char flags[]   = "qhull d Qt Qbb Qc";
	    FILE *outfile  = fopen(dataset_name,"w");
	    FILE *errfile  = stderr;
	    int exitcode;
	    facetT *facet;
	    int nbf,s,curlong,totlong;
	    vertexT *vertex, **vertexp;
	    double Vol,dist;
	    double Volmin  = 1e-6;
	    double maxdist = 2;
       
	    exitcode = qh_new_qhull(N3D,NGrid,grid,ismalloc,flags,outfile,errfile);
  
	    if (!exitcode) { 
	        NTetra      = 0;
		FORALLfacets { if (!facet->upperdelaunay) NTetra++; }

		tetrahedron = new int[(N3D+1)*NTetra];
		tetraVol    = new double[NTetra];

		nbf         = 0;
		FORALLfacets {
		    if (!facet->upperdelaunay) {
		        s = 0;
			FOREACHvertex_(facet->vertices) {
			    tetrahedron[nbf*(N3D+1)+s] = qh_pointid(vertex->point);
			    s++;
			}
			Vol  = calc_tri_volume(&grid[N3D*tetrahedron[nbf*(N3D+1)]],\
					       &grid[N3D*tetrahedron[nbf*(N3D+1)+1]],\
					       &grid[N3D*tetrahedron[nbf*(N3D+1)+2]],\
					       &grid[N3D*tetrahedron[nbf*(N3D+1)+3]]);
			dist = calc_max_distance(&grid[N3D*tetrahedron[nbf*(N3D+1)]],\
						 &grid[N3D*tetrahedron[nbf*(N3D+1)+1]],\
						 &grid[N3D*tetrahedron[nbf*(N3D+1)+2]],\
						 &grid[N3D*tetrahedron[nbf*(N3D+1)+3]]);
		    
			if((Vol>Volmin)&&(dist<maxdist)){
			    tetraVol[nbf] = Vol;
			    nbf++;
			}
		    }
		}
		NTetra = nbf;
        
		qh_freeqhull(!qh_ALL);
		qh_memfreeshort (&curlong, &totlong);
	    }
	    fclose(outfile);
	}else{

	    if(update_tetra){
 
	        int nbf,s;
		FILE *F = fopen(tetra_file,"r");

		fscanf(F,"%i",&NTetra);

		tetrahedron = new int[(N3D+1)*NTetra];
		tetraVol    = new double[NTetra];

		for(nbf=0;nbf<NTetra;nbf++){
		    for(s=0;s<(N3D+1);s++){
		        fscanf(F,"%i",&tetrahedron[nbf*(N3D+1)+s]);
		    }
		    fscanf(F,"%lg",&tetraVol[nbf]);
		}
		fclose(F);
	    }

	    if((update_fitness)||(!poisson_solver)){
  	        NTetra = 1;

                tetrahedron = new int[(N3D+1)*NTetra];                     
                tetraVol    = new double[NTetra];                             

	    }
	}
    }

    if(vg_rank<max_proc_poisson){
	MPI_Bcast(&NTetra,1,MPI_INT,0,po_comm);
    }else{
        NTetra = 1;
    }

    if(vg_rank){
        tetrahedron = new int[NTetra*(N3D+1)];
	tetraVol    = new double[NTetra];
    }

    if(vg_rank<max_proc_poisson){
        MPI_Bcast(tetrahedron,NTetra*(N3D+1),MPI_INT,0,po_comm);
	MPI_Bcast(tetraVol,NTetra,MPI_DOUBLE,0,po_comm);
    }

}

/************************************************************************************************/

void FEMGrid::write_grid(const char *filename)
{
    int i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    for(i=0;i<NGrid;i++){
        for(j=0;j<N3D;j++){
            myfile<<grid[N3D*i+j]<<" ";
        }
        myfile<<"\n";
    }
    myfile.close();
}

/************************************************************************************************/

void FEMGrid::write_tetrahedron(const char *filename)
{
    int i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    myfile<<NTetra<<"\n";
    for(i=0;i<NTetra;i++){
        for(j=0;j<N3D+1;j++){
            myfile<<tetrahedron[(N3D+1)*i+j]<<" ";
        }
        myfile<<tetraVol[i]<<"\n";
    }
    myfile.close();
}

/************************************************************************************************/

void FEMGrid::write_gate(const char *filename)
{
    int i,ndg;
    
    ofstream myfile;
    myfile.open (filename);
    for(ndg=0;ndg<no_diff_gate;ndg++){
        for(i=0;i<NGate[ndg];i++){
	    myfile<<gate_index[ndg][i]<<" ";
	}
	myfile<<endl;
    }
    myfile.close();
}

/************************************************************************************************/

void FEMGrid::write_permittivity(const char *filename)
{
    int i;
    
    ofstream myfile;
    myfile.open (filename);
    for(i=0;i<NGrid;i++){
        myfile<<Eps[i]<<"\n";
    }
    myfile.close();
}

/************************************************************************************************/

double FEMGrid::calc_tri_area(double* p1, double* p2, double* p3, int ND)
{
    double A,b,k;
    double* v1 = new double[ND];
    double* v2 = new double[ND];
    
   c_dcopy(ND,p2,1,v1,1);
   c_daxpy(ND,-1.0,p1,1,v1,1);
   c_dcopy(ND,p3,1,v2,1);
   c_daxpy(ND,-1.0,p1,1,v2,1);
    
    b = c_dnrm2(ND,v1,1);
    k = c_ddot(ND,v1,1,v2,1)/pow(b,2.0);
    c_daxpy(ND,-k,v1,1,v2,1);
    A = b*c_dnrm2(ND,v2,1)/2;
    
    if (v1!=NULL){
        delete[] v1;
    }
    if (v2!=NULL){
        delete[] v2;
    }
    return A;
}

/************************************************************************************************/

void FEMGrid::calc_vec(double* v, double* p1, double* p2, double* p3)
{
    double v1[3];
    double v2[3];

    c_dcopy(N3D,p2,1,v1,1);
    c_daxpy(N3D,-1.0,p1,1,v1,1);
    c_dcopy(N3D,p3,1,v2,1);
    c_daxpy(N3D,-1.0,p1,1,v2,1);

    v[0] = v1[1]*v2[2]-v1[2]*v2[1];
    v[1] = -v1[0]*v2[2]+v1[2]*v2[0];
    v[2] = v1[0]*v2[1]-v1[1]*v2[0];

    c_dscal(N3D,1/c_dnrm2(N3D,v,1),v,1);
}

/************************************************************************************************/

double FEMGrid::calc_tri_volume(double *p1, double *p2, double *p3, double *p4)
{
    double p4_p1[3];
    double v123[3];
    double A123;
    
    calc_vec(v123,p1,p2,p3);
    
    A123 = calc_tri_area(p1,p2,p3,N3D);
    
    c_dcopy(N3D,p4,1,p4_p1,1);
    c_daxpy(N3D,-1.0,p1,1,p4_p1,1);

    return fabs(c_ddot(N3D,v123,1,p4_p1,1)*A123/3.0);
}

/************************************************************************************************/

double FEMGrid::calc_max_distance(double *p1, double *p2, double *p3, double *p4)
{
    int i,j;
    double d,dmax;
    double vec[3];
    double p[4][3];

    c_dcopy(N3D,p1,1,p[0],1);
    c_dcopy(N3D,p2,1,p[1],1);
    c_dcopy(N3D,p3,1,p[2],1);
    c_dcopy(N3D,p4,1,p[3],1);

    dmax = 0.0;
    for(i=0;i<4;i++){
       for(j=i+1;j<4;j++){
	   c_dcopy(N3D,p[i],1,vec,1);
	   c_daxpy(N3D,-1.0,p[j],1,vec,1);
	   d = c_dnrm2(N3D,vec,1);
	   if(d>dmax) dmax = d;
       }
    }

    return dmax;
}

/************************************************************************************************/

double FEMGrid::calc_tri_volume(double A, double* p0, double* p1, double* v)
{
    double p1_p0[3];

    c_dcopy(N3D,p1,1,p1_p0,1);
    c_daxpy(N3D,-1.0,p0,1,p1_p0,1);

    return fabs(c_ddot(N3D,v,1,p1_p0,1)*A/3.0);
}

/************************************************************************************************/

void FEMGrid::calc_doping_info(WireStructure* nanowire)
{
    int no_doping=nanowire->no_doping,i;

    for(i=0;i<no_doping;i++){
        if(!strcmp(nanowire->doping[i]->type,"square")){
            nanowire->doping[i]->face_area[0]=calc_quad_area(nanowire->doping[i]->p[0]->coord,\
                                                             nanowire->doping[i]->p[1]->coord,\
                                                             nanowire->doping[i]->p[2]->coord,\
                                                             nanowire->doping[i]->p[3]->coord,N3D);
            c_dcopy(N3D,nanowire->doping[i]->p[0]->coord,1,nanowire->doping[i]->p_ref[0],1);
            calc_vec(nanowire->doping[i]->vec_dir[0],nanowire->doping[i]->p[0]->coord,\
                     nanowire->doping[i]->p[1]->coord,nanowire->doping[i]->p[2]->coord);
        
            nanowire->doping[i]->face_area[1]=calc_quad_area(nanowire->doping[i]->p[4]->coord,\
                                                             nanowire->doping[i]->p[5]->coord,\
                                                             nanowire->doping[i]->p[6]->coord,\
                                                             nanowire->doping[i]->p[7]->coord,N3D);
            c_dcopy(N3D,nanowire->doping[i]->p[4]->coord,1,nanowire->doping[i]->p_ref[1],1);
            calc_vec(nanowire->doping[i]->vec_dir[1],nanowire->doping[i]->p[4]->coord,\
                     nanowire->doping[i]->p[5]->coord,nanowire->doping[i]->p[6]->coord);

            nanowire->doping[i]->face_area[2]=calc_quad_area(nanowire->doping[i]->p[0]->coord,\
                                                             nanowire->doping[i]->p[4]->coord,\
                                                             nanowire->doping[i]->p[7]->coord,\
                                                             nanowire->doping[i]->p[3]->coord,N3D);
            c_dcopy(N3D,nanowire->doping[i]->p[0]->coord,1,nanowire->doping[i]->p_ref[2],1);
            calc_vec(nanowire->doping[i]->vec_dir[2],nanowire->doping[i]->p[0]->coord,\
                     nanowire->doping[i]->p[4]->coord,nanowire->doping[i]->p[7]->coord);
        
            nanowire->doping[i]->face_area[3]=calc_quad_area(nanowire->doping[i]->p[0]->coord,\
                                                             nanowire->doping[i]->p[1]->coord,\
                                                             nanowire->doping[i]->p[5]->coord,\
                                                             nanowire->doping[i]->p[4]->coord,N3D);
            c_dcopy(N3D,nanowire->doping[i]->p[0]->coord,1,nanowire->doping[i]->p_ref[3],1);
            calc_vec(nanowire->doping[i]->vec_dir[3],nanowire->doping[i]->p[0]->coord,\
                     nanowire->doping[i]->p[1]->coord,nanowire->doping[i]->p[5]->coord);
        
            nanowire->doping[i]->face_area[4]=calc_quad_area(nanowire->doping[i]->p[1]->coord,\
                                                             nanowire->doping[i]->p[2]->coord,\
                                                             nanowire->doping[i]->p[6]->coord,\
                                                             nanowire->doping[i]->p[5]->coord,N3D);
            c_dcopy(N3D,nanowire->doping[i]->p[1]->coord,1,nanowire->doping[i]->p_ref[4],1);
            calc_vec(nanowire->doping[i]->vec_dir[4],nanowire->doping[i]->p[1]->coord,\
                     nanowire->doping[i]->p[2]->coord,nanowire->doping[i]->p[6]->coord);
        
            nanowire->doping[i]->face_area[5]=calc_quad_area(nanowire->doping[i]->p[2]->coord,\
                                                             nanowire->doping[i]->p[3]->coord,\
                                                             nanowire->doping[i]->p[7]->coord,\
                                                             nanowire->doping[i]->p[6]->coord,N3D);
            c_dcopy(N3D,nanowire->doping[i]->p[2]->coord,1,nanowire->doping[i]->p_ref[5],1);
            calc_vec(nanowire->doping[i]->vec_dir[5],nanowire->doping[i]->p[2]->coord,\
                     nanowire->doping[i]->p[3]->coord,nanowire->doping[i]->p[7]->coord);

            nanowire->doping[i]->volume=calc_quad_volume(nanowire->doping[i]);

        }
        
    }
    
}

/************************************************************************************************/

double FEMGrid::calc_quad_area(double* p1, double* p2, double* p3, double* p4, int ND)
{
    return calc_tri_area(p1,p2,p4,ND)+calc_tri_area(p3,p2,p4,ND);
}

/************************************************************************************************/

double FEMGrid::calc_quad_volume(DOPING* doping)
{
    double V1,V2,V3;

    V1 = calc_tri_volume(doping->face_area[1],doping->p[0]->coord,doping->p_ref[1],\
                         doping->vec_dir[1]);
    V2 = calc_tri_volume(doping->face_area[4],doping->p[0]->coord,doping->p_ref[4],\
                         doping->vec_dir[4]);
    V3 = calc_tri_volume(doping->face_area[5],doping->p[0]->coord,doping->p_ref[5],\
                         doping->vec_dir[5]);

    return V1+V2+V3;
}

/************************************************************************************************/

int FEMGrid::is_in_quad3D(DOPING* doping, double* p)
{
    int is_in = 0;

    if(!strcmp(doping->type,"square")){
    
        double V = 0.0;
        int i    = 0;
    
        for(i=0;i<6;i++){
            V = V+calc_tri_volume(doping->face_area[i],p,doping->p_ref[i],doping->vec_dir[i]);
        }

        if((fabs(V-doping->volume)/V)<tollim) is_in = 1;
        
    }else{

        bool condition1 = p[0]>=doping->p[0]->coord[0];
        bool condition2 = p[0]<=doping->p[1]->coord[0];
        bool condition3 = (pow(p[1]-doping->p[0]->coord[1],2.0)/pow(doping->radius[0],2.0)+\
			   pow(p[2]-doping->p[0]->coord[2],2.0)/pow(doping->radius[1],2.0))<=1.0;

        is_in = condition1&&condition2&&condition3;
        
    }

    return is_in;
}

/************************************************************************************************/

int FEMGrid::is_in_quad3D(MAT* mat, double* p)
{
    int is_in = 0;

    if(!strcmp(mat->type,"square")){
    
        double V = 0.0;
        int i    = 0;
    
        for(i=0;i<6;i++){
            V = V+calc_tri_volume(mat->face_area[i],p,mat->p_ref[i],mat->vec_dir[i]);
        }

        if((fabs(V-mat->volume)/V)<tollim) is_in = 1;
        
    }else{

        bool condition1 = p[0]>=mat->p[0]->coord[0];
        bool condition2 = p[0]<=mat->p[1]->coord[0];
        bool condition3 = (pow(p[1]-mat->p[0]->coord[1],2.0)/pow(mat->radius[0],2.0)+\
			   pow(p[2]-mat->p[0]->coord[2],2.0)/pow(mat->radius[1],2.0))<=1.0;

        is_in = condition1&&condition2&&condition3;
        
    }

    return is_in;
}

/************************************************************************************************/

int FEMGrid::is_in_quad2D(GATE* gate, double* p)
{
    int i,is_in=0;
    double A=0;

    for(i=0;i<3;i++){
        A=A+calc_tri_area(gate->p[i]->coord,gate->p[i+1]->coord,p,N3D);
    }
    A=A+calc_tri_area(gate->p[3]->coord,gate->p[0]->coord,p,N3D);

    if((fabs(A-gate->face_area)/A)<tollim) is_in=1;

    return is_in;
}

/************************************************************************************************/

double FEMGrid::check_volume(MAT* mat, double* p)
{
    double V = 0.0;

    if(!strcmp(mat->type,"square")){
    
        int i;
    
        for(i=0;i<6;i++){
            V = V+calc_tri_volume(mat->face_area[i],p,mat->p_ref[i],mat->vec_dir[i]);
        }
        
    }else{

        V = (mat->p[1]->coord[0]-mat->p[0]->coord[0])*PI*\
	  sqrt(pow(p[1]-mat->p[0]->coord[1],2.0)+pow(p[0]-mat->p[0]->coord[0],2.0));
        
    }

    return V;
}

/************************************************************************************************/

void FEMGrid::check_dimension(WireGenerator *Wire,WireStructure *nanowire)
{

    if(nanowire->NDim==2){

        int no_doping = nanowire->no_doping;
        int no_gate   = nanowire->no_gate;
	int no_ground = nanowire->no_ground;
        int IE,ID,IG;
        
        for(IE=0;IE<no_doping;IE++){

            for(ID=0;ID<N3D+1;ID++){

                c_dcopy(N3D,nanowire->doping[IE]->p[ID]->coord,1,\
                        nanowire->doping[IE]->p[ID+N3D+1]->coord,1);
                nanowire->doping[IE]->p[ID]->coord[2]       = -Wire->Lz-10*tollim;
                nanowire->doping[IE]->p[ID+N3D+1]->coord[2] = 2*Wire->Lz-10*tollim;
                
            }
        }

        for(IG=0;IG<no_gate;IG++){

            for(ID=0;ID<nanowire->NDim;ID++){
                
                c_dcopy(nanowire->NDim,nanowire->gate[IG]->p[ID]->coord,1,\
                        nanowire->gate[IG]->p[N3D-ID]->coord,1);
                nanowire->gate[IG]->p[ID]->coord[2]     = -Wire->Lz-10*tollim;
                nanowire->gate[IG]->p[N3D-ID]->coord[2] = 2*Wire->Lz-10*tollim;
                
            }
        }

	for(IG=0;IG<no_ground;IG++){

            for(ID=0;ID<nanowire->NDim;ID++){
                
                c_dcopy(nanowire->NDim,nanowire->ground[IG]->p[ID]->coord,1,\
                        nanowire->ground[IG]->p[N3D-ID]->coord,1);
                nanowire->ground[IG]->p[ID]->coord[2]     = -Wire->Lz-10*tollim;
                nanowire->ground[IG]->p[N3D-ID]->coord[2] = 2*Wire->Lz-10*tollim;
                
            }
        }

        
    }
}

/************************************************************************************************/
