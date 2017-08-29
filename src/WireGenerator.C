/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/****************************************************************************************
All the matrices used in BLAS and Lapack are stored in the fortran style:
                         if A(N,M), A[i][j] stored as A[i+j*N]
*****************************************************************************************/

#include "WireGenerator.H"

WireGenerator::WireGenerator(int plattice_type,int ptransport_type)
{
    int IA,IB,ID,IP;
    int BZNP;
    
    tb             = 10;
    N3D            = 3;
    N2D            = 2;
    N1D            = 1;

    lattice_type   = plattice_type;
    transport_type = ptransport_type;

    switch(lattice_type){
    case 1: //zincblende
        NA            = 2;
	NB            = 4;
	SP            = 4;
        BZNP          = 5;
	max_neighbors = 40;
	break;
    case 3: //wurtzite
        NA            = 4;
	NB            = 4;
	SP            = 4;
        BZNP          = 8;
	max_neighbors = 40;
	break;
    case 4: //graphene
        NA            = 2;
        NB            = 3;
	SP            = 4;
        BZNP          = 4;
	max_neighbors = 18;
	break;
    case 5: //rhombohedral
        NA            = 5;
	NB            = 6;
	SP            = 4;
        BZNP          = 6;
	max_neighbors = 0;
        break;
    case 6: //cnt
        NA            = 2;
        NB            = 3;
	SP            = 4;
        BZNP          = 4;
	max_neighbors = 18;
	break;
    case 7: //bilayer graphene
        NA            = 4;
	NB            = 4;
	SP            = 4;
        BZNP          = 4;
	max_neighbors = 18;
	break;
    case 8: //multilayer graphene
        NA            = 4;
	NB            = 5;
	SP            = 4;
        BZNP          = 4;
	max_neighbors = 18;
	break;
    case 9: //rocksalt
        NA            = 2;
	NB            = 6;
	SP            = 4;
        BZNP          = 6;
	max_neighbors = 0;
	break;
    case 10: //dichalcogenide
        NA            = 3;
	NB            = 6;
	SP            = 4;
        BZNP          = 4;
	max_neighbors = 40;
	break;
    case 11: //multilayer dichalcogenide
        NA            = 6;
	NB            = 6;
	SP            = 4;
        BZNP          = 4;
	max_neighbors = 40;
	break;
    default:
        NA            = 2;
	NB            = 4;
	SP            = 4;
        BZNP          = 5;
	max_neighbors = 40;
    }

    SLM = 4+NB;

    unit_cell                = new Unit_Cell();
    unit_cell->atom          = new Atom*[NA];
    unit_cell->type          = new at_type*[NA];
    unit_cell->type_original = new at_type*[NA];
    unit_cell->axis          = new Axis*[N3D];
    unit_cell->BZPoint       = new POINT3D*[BZNP];

    unit_cell->BZNP          = BZNP;
    
    for(IA=0;IA<NA;IA++){
        unit_cell->atom[IA]                = new Atom();
        unit_cell->type[IA]                = new at_type();
        unit_cell->type_original[IA]       = new at_type();
        unit_cell->type[IA]->bond          = new Bond*[NB];
        unit_cell->type_original[IA]->bond = new Bond*[NB];

        for(IB=0;IB<NB;IB++){
            unit_cell->type[IA]->bond[IB]          = new Bond();
            unit_cell->type_original[IA]->bond[IB] = new Bond();
        }
    }
    
    for(ID=0;ID<N3D;ID++){
        unit_cell->axis[ID] = new Axis();
    }

    for(IP=0;IP<BZNP;IP++){
        unit_cell->BZPoint[IP] = new POINT3D();    
    }

    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
    
}

/************************************************************************************************/

WireGenerator::~WireGenerator()
{

    int IA,IB,ID,IP;

    if(nwire->NDim){
            
        delete_variables();

        for(IA=0;IA<NA;IA++){
            for(IB=0;IB<NB;IB++){
                delete unit_cell->type[IA]->bond[IB];
                delete unit_cell->type_original[IA]->bond[IB];
            }
            delete unit_cell->type[IA]->bond;
            delete unit_cell->type_original[IA]->bond;
            delete unit_cell->atom[IA];
            delete unit_cell->type[IA];
            delete unit_cell->type_original[IA];
        }

        for(ID=0;ID<N3D;ID++){
            delete unit_cell->axis[ID];    
        }

        for(IP=0;IP<unit_cell->BZNP;IP++){
            delete unit_cell->BZPoint[IP];        
        }

        delete unit_cell->atom;
        delete unit_cell->type;
        delete unit_cell->type_original;
        delete unit_cell->axis;
        delete unit_cell->BZPoint;
        delete[] atomic_mass;
        delete[] EMidGap;
        delete[] EGap;
        delete[] Lmin;
        delete[] Lmax;
        delete[] Around_Matrix;
    }

    delete[] Layer_Matrix;
    delete[] orb_per_at; 

}

/************************************************************************************************/

void WireGenerator::delete_variables()
{

    int i;

    delete[] Neigh_Matrix;
    delete[] Smin;
    delete[] Smax;
    delete[] Smin_tot;
    delete[] Smax_tot;
    delete[] neighbor_layer;
    delete[] index_channel;
    delete[] ch_pos;
    delete[] ch_conv;
    delete[] inv_ch_pos;
    delete[] index_arch;
    delete[] arch_pos;
    delete[] arch_conv;
    delete[] b_shift;

    for(i=0;i<2;i++){
        delete[] index_boundary[i];
        delete[] bound_pos[i];
        delete[] bound_conv[i];
        delete[] Boundary[i]->UNN;
        delete[] Boundary[i]->UNNp;
        delete[] Boundary[i]->UNNm;
        delete[] Boundary[i]->UNNneigh;
        delete Boundary[i];
    }
    delete[] index_boundary;
    delete[] bound_pos;
    delete[] bound_conv;
    delete Boundary;

    for(i=0;i<NLayer;i++){
        delete[] Enn[i]->UNN;
        delete[] Enn[i]->UNNp;
        delete[] Enn[i]->UNNneigh;
        delete Enn[i];
    }
    delete Enn;

}

/************************************************************************************************/

void WireGenerator::execute_simple(WireStructure* nanowire,MPI_Comm int_comm)
{
    if(!mpi_rank){
        printf("Generating the atomic structure\n");
    }

    wg_comm = int_comm;
    MPI_Comm_size(wg_comm,&mpi_size);
    MPI_Comm_rank(wg_comm,&mpi_rank);

    if(nanowire->hydrogen_passivation){
        nanowire->dsp3 = 0.0;
    }

    nwire   = nanowire;
    make_unit_cell(nanowire->first_atom,nanowire->a0,nanowire->c0,nanowire->u0);
    change_orientation(nanowire->x,nanowire->y,nanowire->z,nanowire->dsp3);

    if(!nanowire->NDim){
        Layer_Matrix = NULL;
        orb_per_at = NULL;
    }else{
        check_dimension(nanowire);
        calc_info(nanowire);
        make_wire(nanowire);
        make_connections(nanowire);
        hydrogen_passivation(nanowire);    
        generate_roughness(nanowire,0);
        replicate_unit_cell(nanowire->replicate_unit_cell,nanowire);
        read_atom_position(nanowire->read_atom_pos,nanowire);
	check_periodicity(nanowire);
        init_variables(nanowire);
	cut_boundary_slab(nanowire);
        update_atom_pos(nanowire,NULL);
	generate_alloy_disorder(nanowire);
        cut_layer(nanowire);
        separate_dimension(nanowire);
        convert_position(nanowire);
	update_periodicity(nanowire);
        get_qm_region(nanowire);
        orb_per_at = NULL;
        EMidGap = NULL;
        EGap = NULL;
        add_strain(nanowire->strain,nanowire->no_strain_domain,\
                   nanowire->strain_domain,nanowire->update_at);
        atomic_mass = NULL;
    }    
}

/************************************************************************************************/

void WireGenerator::execute_task(WireStructure* nanowire,Material *material,int sample_id,\
				 MPI_Comm int_comm)
{
    if(!mpi_rank){
        printf("Generating the atomic structure\n");
    }

    wg_comm = int_comm;
    MPI_Comm_size(wg_comm,&mpi_size);
    MPI_Comm_rank(wg_comm,&mpi_rank);

    if(nanowire->hydrogen_passivation){
        nanowire->dsp3 = 0.0;
    }

    nwire   = nanowire;
    make_unit_cell(nanowire->first_atom,nanowire->a0,nanowire->c0,nanowire->u0);
    change_orientation(nanowire->x,nanowire->y,nanowire->z,nanowire->dsp3);

    if(!nanowire->NDim){

        if(nanowire->strain->on){
            make_strained_unit_cell(nanowire->strain);    
        }
        make_bulk_layer_matrix(nanowire,material->no_orb);
            
    }else{
        check_dimension(nanowire);
        calc_info(nanowire);
        make_wire(nanowire);
        make_connections(nanowire);
        hydrogen_passivation(nanowire);    
        generate_roughness(nanowire,sample_id);
        replicate_unit_cell(nanowire->replicate_unit_cell,nanowire);
        read_atom_position(nanowire->read_atom_pos,nanowire);
	check_periodicity(nanowire);
        init_variables(nanowire);
	cut_boundary_slab(nanowire);
        update_atom_pos(nanowire,material);
	generate_alloy_disorder(nanowire);
        cut_layer(nanowire);
        separate_dimension(nanowire);
        convert_position(nanowire);
	update_periodicity(nanowire);
        get_qm_region(nanowire);
        get_orbital_per_atom(material->no_orb);
        get_mid_gap_energy(material);
        add_strain(nanowire->strain,nanowire->no_strain_domain,\
                   nanowire->strain_domain,nanowire->update_at);
        adapt_parameters(material);
    }    
}

/************************************************************************************************/

void WireGenerator::update_atom_pos(WireStructure* nanowire,Material *material)
{
    int IA,IP,IZ;
    int IPmax[2] = {N3D+1,SLM};

    if(nanowire->update_at){
        FILE *F = fopen(nanowire->at_file,"r");
        switch(nanowire->NDim){
	case 1:
	    for(IA=0;IA<No_Atom;IA++){
	        for(IP=0;IP<IPmax[nanowire->update_at-1];IP++){
		    fscanf(F,"%lg",&Layer_Matrix[SLM*IA+IP]);
		}
	    }
	    fclose(F);
	    break;
	case 2:
	    for(IA=0;IA<No_Atom;IA=IA+3*NzFold){
	        for(IZ=0;IZ<NzFold;IZ++){
		    for(IP=0;IP<IPmax[nanowire->update_at-1];IP++){
		        fscanf(F,"%lg",&Layer_Matrix[SLM*(IA+IZ)+IP]);
		    }
		    c_dcopy(IPmax[nanowire->update_at-1],&Layer_Matrix[SLM*(IA+IZ)],1,\
			    &Layer_Matrix[SLM*(IA+IZ+NzFold)],1);
		    c_dcopy(IPmax[nanowire->update_at-1],&Layer_Matrix[SLM*(IA+IZ)],1,\
			    &Layer_Matrix[SLM*(IA+IZ+2*NzFold)],1);
		    Layer_Matrix[SLM*(IA+IZ)+2]          = Layer_Matrix[SLM*(IA+IZ)+2]-Lz;
		    Layer_Matrix[SLM*(IA+IZ+2*NzFold)+2] = Layer_Matrix[SLM*(IA+IZ+2*NzFold)+2]+Lz;
		}
	    }
	    fclose(F);
	    break;
	case 3:
	    for(IA=0;IA<No_Atom;IA++){
	        for(IP=0;IP<IPmax[nanowire->update_at-1];IP++){
		    fscanf(F,"%lg",&Layer_Matrix[SLM*IA+IP]);
		}
	    }
	    fclose(F);
	    break;
	}
    }else{
        if (nanowire->relax_atoms){
	  /*
	    MacoptStrain *strain = new MacoptStrain(this, material);

	    strain->relax_atoms();
	    delete strain;
	  */
	}
    }

    if((nanowire->update_at)||(nanowire->relax_atoms)){
        get_cell_width();
        update_boundary_slab(nanowire);
    }
}

/************************************************************************************************/
    
void WireGenerator::make_connections(WireStructure* nanowire)
{
    MPI_Status status;
    int IA,IB,i,j,IPIV[3],INFO,type,next_type,i_proc,no_line,max_rank;
    int At_local=ceil((double)No_Atom/mpi_size),No_Atom_local,Nx,Ny,Nz,entry;
    double V[9],invV[9]={0},next_pos[3],vec_trans[3];
    double pos_inc[3] = {hxmin,hymin,hzmin};

    Nx = hxmax-hxmin+1;
    Ny = hymax-hymin+1;
    Nz = hzmax-hzmin+1;
    
    for(i=0;i<N3D;i++){
        invV[i*(N3D+1)]=1;
        for(j=0;j<N3D;j++){
            V[i+j*N3D] = unit_cell->axis[j]->vec[i];
        }
    }
    c_dgetrf(N3D,N3D,V,N3D,IPIV,&INFO);
    c_dgetrs('N',N3D,N3D,V,N3D,IPIV,invV,N3D,&INFO);
    
    No_Atom_local = min(No_Atom,(mpi_rank+1)*At_local)-mpi_rank*At_local;
    max_rank      = ceil((double)No_Atom/At_local);

    if(mpi_rank<max_rank){
      
        for(IA=mpi_rank*At_local;IA<min(No_Atom,(mpi_rank+1)*At_local);IA++){

            type = (Round(abs(Layer_Matrix[SLM*IA+3]))-1)%NA;
	    
            for(IB=0;IB<NB;IB++){

	        next_type = unit_cell->type[type]->bond[IB]->neigh-1;

		if(next_type>=0){

		    c_dcopy(N3D,&Layer_Matrix[SLM*IA],1,next_pos,1);
		    c_daxpy(N3D,-1.0,unit_cell->atom[next_type]->coord,1,next_pos,1);
		    c_daxpy(N3D,1.0,unit_cell->type[type]->bond[IB]->vec,1,next_pos,1);
		    c_dcopy(N3D,next_pos,1,vec_trans,1);
            
		    c_dgemv('N',N3D,N3D,1.0,invV,N3D,vec_trans,1,0.0,next_pos,1);
		    c_daxpy(N3D,-1,pos_inc,1,next_pos,1);

		    entry = Round(next_type*Nx*Ny*Nz+next_pos[0]*Ny*Nz+next_pos[1]*Nz+next_pos[2]);

		    if((entry>=0)&&(entry<NA*Nx*Ny*Nz)){
		        Layer_Matrix[SLM*IA+4+IB] = HXYZ[entry];
		    }
		}
            }
        }
	
        if(mpi_rank){
            MPI_Send(&No_Atom_local,1,MPI_INT,0,0,wg_comm);
            MPI_Send(Layer_Matrix+(SLM*mpi_rank*At_local),SLM*No_Atom_local,MPI_DOUBLE,0,1,\
                     wg_comm);
        }else{
            for(i_proc=1;i_proc<max_rank;i_proc++){
                MPI_Recv(&no_line,1,MPI_INT,i_proc,0,wg_comm,&status);
                MPI_Recv(Layer_Matrix+(SLM*i_proc*At_local),SLM*no_line,MPI_DOUBLE,i_proc,\
                         1,wg_comm,&status);
            }
        }
    }

    MPI_Bcast(Layer_Matrix,SLM*No_Atom,MPI_DOUBLE,0,wg_comm);

    delete[] HXYZ;
    delete[] XYZH;

}

/************************************************************************************************/

void WireGenerator::make_wire(WireStructure* nanowire)
{
    
    int Nx,Ny,Nz,Nx_local,IX,IY,IZ,IM,IR,IA,max_rank;
    int no_element,no_ch_element,no_ox_element,APNA,APANA;
    int condition,ro_condition;
    int At_local,Around_local,Rough_local;
    double shiftx[3]={0},shiftxy[3],shiftxyz[3],at_pos[4];
    double *LM_local,*XYZH_local,*AM_local;

    no_element    = nanowire->no_element;
    no_ch_element = nanowire->no_ch_element;
    no_ox_element = nanowire->no_ox_element;
    
    hxmin = INF;
    hymin = INF;
    hzmin = INF;
    hzmin = INF;
    hxmax = -INF;
    hymax = -INF;
    hzmax = -INF;
    hzmax = -INF;
    find_3D_boundary(&hxmin,&hxmax,&hymin,&hymax,&hzmin,&hzmax,nanowire);

    Nx           = hxmax-hxmin+1;
    Ny           = hymax-hymin+1;
    Nz           = hzmax-hzmin+1;
    Nx_local     = ceil((double)Nx/mpi_size);
    max_rank     = ceil((double)Nx/Nx_local);

    APNA         = Round(3.0*volume/(Vol_atom*1e27));
    APANA        = Round(3.0*(volume_tot-volume)/(Vol_atom*1e27));

    LM_local     = new double[APNA*SLM];
    XYZH_local   = new double[APNA*N3D];
    AM_local     = new double[max(APANA,1)*N3D];
    
    At_local     = 0;
    Around_local = 0;
    Rough_local  = 0;
    
    init_var(LM_local,APNA*SLM);

    if(mpi_rank<max_rank){

        c_daxpy(N3D,(hxmin+mpi_rank*Nx_local),unit_cell->axis[0]->vec,1,shiftx,1);

        for(IX=0;IX<Nx_local;IX++){

            c_dcopy(N3D,shiftx,1,shiftxy,1);
            c_daxpy(N3D,hymin,unit_cell->axis[1]->vec,1,shiftxy,1);

            for(IY=0;IY<Ny;IY++){

                c_dcopy(N3D,shiftxy,1,shiftxyz,1);
                c_daxpy(N3D,hzmin,unit_cell->axis[2]->vec,1,shiftxyz,1);
            
                for(IZ=0;IZ<Nz;IZ++){

                    c_dcopy(N3D,shiftxyz,1,at_pos,1);
		    
                    for(IA=0;IA<NA;IA++){

                        c_daxpy(N3D,1.0,unit_cell->atom[IA]->coord,1,at_pos,1);

                        for(IM=0;IM<no_element;IM++){
			  
                            condition = is_in_quad3D(nanowire->mat[IM],at_pos);
			  
                            if(condition){

			        at_pos[N3D] = NA*(nanowire->mat[IM]->id_number-1)+IA+1;

                                if(IM<no_ch_element){
                                    At_local++;
                                    c_dcopy(N3D+1,at_pos,1,&LM_local[SLM*(At_local-1)],1);
                                    XYZH_local[N3D*(At_local-1)]   = IX+mpi_rank*Nx_local;
                                    XYZH_local[N3D*(At_local-1)+1] = IY;
                                    XYZH_local[N3D*(At_local-1)+2] = IZ;
                                }else{
				    ro_condition = 0;
				    for(IR=no_ch_element+no_ox_element;IR<no_element;IR++){
				        ro_condition = is_in_quad3D(nanowire->mat[IR],at_pos);
					if(ro_condition){
					    At_local++;
					    Rough_local++;
                                            c_dcopy(N3D+1,at_pos,1,&LM_local[SLM*(At_local-1)],1);
					    LM_local[SLM*(At_local-1)+3]   = -LM_local[SLM*(At_local-1)+3];
                                            XYZH_local[3*(At_local-1)]     = IX+mpi_rank*Nx_local;
                                            XYZH_local[3*(At_local-1)+1]   = IY;
                                            XYZH_local[3*(At_local-1)+2]   = IZ;
					    break;
					}
				    }
				    if(!ro_condition){
				        Around_local++;
                                        c_dcopy(N3D,at_pos,1,&AM_local[3*(Around_local-1)],1);
				    }
                                }
                                break;
                            }
                        }
                        c_dcopy(N3D,shiftxyz,1,at_pos,1);
                    }
                    c_daxpy(N3D,1.0,unit_cell->axis[2]->vec,1,shiftxyz,1);
                }
                c_daxpy(N3D,1.0,unit_cell->axis[1]->vec,1,shiftxy,1);
            }
            c_daxpy(N3D,1.0,unit_cell->axis[0]->vec,1,shiftx,1);
        }
        exchange_wire_info(At_local,LM_local,XYZH_local,max_rank);
        exchange_oxide_info(Around_local,Rough_local,AM_local,max_rank);
    }

    //Copy the information to all processors
    bcast_wire_info(Nx,Ny,Nz);
    bcast_oxide_info();
    
    delete[] LM_local;
    delete[] XYZH_local;
    delete[] AM_local;

}

/************************************************************************************************/

void WireGenerator::bcast_wire_info(int Nx,int Ny,int Nz)
{
    int IA,pos_ind;
    
    MPI_Bcast(&No_Atom,1,MPI_INT,0,wg_comm);
    if(mpi_rank){
        Layer_Matrix = new double[SLM*No_Atom];
        XYZH         = new double[N3D*No_Atom];
    }
    MPI_Bcast(Layer_Matrix,SLM*No_Atom,MPI_DOUBLE,0,wg_comm);
    MPI_Bcast(XYZH,N3D*No_Atom,MPI_DOUBLE,0,wg_comm);

    MPI_Bcast(&NLayer,1,MPI_INT,0,wg_comm);
    if(mpi_rank){
        Lmin = new int[NLayer];
        Lmax = new int[NLayer];
    }
    MPI_Bcast(Lmin,NLayer,MPI_INT,0,wg_comm);
    MPI_Bcast(Lmax,NLayer,MPI_INT,0,wg_comm);

    HXYZ = new double[NA*Nx*Ny*Nz];
    init_var(HXYZ,NA*Nx*Ny*Nz);
    
    for(IA=0;IA<No_Atom;IA++){
        pos_ind = Round(abs(Layer_Matrix[SLM*IA+3])-1)%NA;
        HXYZ[Round(pos_ind*Nx*Ny*Nz+XYZH[N3D*IA]*Ny*Nz+XYZH[N3D*IA+1]*Nz+XYZH[N3D*IA+2])] = IA+1;
    }
}

/************************************************************************************************/

void WireGenerator::bcast_oxide_info()
{

    MPI_Bcast(&Around_Atom,1,MPI_INT,0,wg_comm);
    MPI_Bcast(&Rough_Atom,1,MPI_INT,0,wg_comm);
    if(mpi_rank){
        Around_Matrix = new double[max(N3D*(Around_Atom+2*Rough_Atom),1)];
    }
    if(Around_Atom){
        MPI_Bcast(Around_Matrix,N3D*Around_Atom,MPI_DOUBLE,0,wg_comm);
    }
}

/************************************************************************************************/

void WireGenerator::exchange_wire_info(int At_local,double *LM_local,double *XYZH_local,\
                                       int max_rank)
{
    MPI_Status status;
    int IA,i_proc,LM_size,XYZH_size,*At_vec;

    At_vec = new int[mpi_size];

    if(mpi_rank){

        MPI_Send(&At_local,1,MPI_INT,0,0,wg_comm);
        if(At_local){
            MPI_Send(LM_local,SLM*At_local,MPI_DOUBLE,0,1,wg_comm);
            MPI_Send(XYZH_local,N3D*At_local,MPI_DOUBLE,0,2,wg_comm);
        }
    }else{

        At_vec[0]    = At_local;
        No_Atom      = At_local;

        for(i_proc=1;i_proc<max_rank;i_proc++){
            MPI_Recv(&At_vec[i_proc],1,MPI_INT,i_proc,0,wg_comm,&status);
            No_Atom  = No_Atom+At_vec[i_proc];
        }

        Layer_Matrix = new double[SLM*No_Atom];
        XYZH         = new double[N3D*No_Atom];

        XYZH_size    = N3D*At_vec[0];
        LM_size      = SLM*At_vec[0];
        c_dcopy(SLM*At_vec[0],LM_local,1,Layer_Matrix,1);
        c_dcopy(N3D*At_vec[0],XYZH_local,1,XYZH,1);

        for(i_proc=1;i_proc<max_rank;i_proc++){

            if(At_vec[i_proc]){
                MPI_Recv(&Layer_Matrix[LM_size],SLM*At_vec[i_proc],MPI_DOUBLE,i_proc,1,\
                         wg_comm,&status);
                MPI_Recv(&XYZH[XYZH_size],N3D*At_vec[i_proc],MPI_DOUBLE,i_proc,2,\
                         wg_comm,&status);
                LM_size   = LM_size+SLM*At_vec[i_proc];
                XYZH_size = XYZH_size+N3D*At_vec[i_proc];
            }
        }

        XYZPOS* LM_xyz  = new XYZPOS[No_Atom];
        int *xyz_index  = new int[No_Atom];
        int *Lmin_local = new int[No_Atom];
        int *Lmax_local = new int[No_Atom];

        for(IA=0;IA<No_Atom;IA++){
            LM_xyz[IA].x     = Layer_Matrix[SLM*IA];
            LM_xyz[IA].y     = Layer_Matrix[SLM*IA+1];
            LM_xyz[IA].z     = Layer_Matrix[SLM*IA+2];
            LM_xyz[IA].index = IA;
        }
        sort_xyz(LM_xyz,xyz_index,&NLayer,Lmin_local,Lmax_local);
        dreshape(No_Atom,SLM,Layer_Matrix,xyz_index);
        dreshape(No_Atom,N3D,XYZH,xyz_index);

        Lmin = new int[NLayer];
        Lmax = new int[NLayer];

        icopy(NLayer,Lmin_local,Lmin);
        icopy(NLayer,Lmax_local,Lmax);
        
        delete LM_xyz;
        delete[] xyz_index;
        delete[] Lmin_local;
        delete[] Lmax_local;
    }

    delete[] At_vec;
}

/************************************************************************************************/

void WireGenerator::exchange_oxide_info(int Around_local,int Rough_local,double *AM_local,\
					int max_rank)
{
    MPI_Status status;
    int AM_size,i_proc,*Around_vec;

    Around_vec = new int[mpi_size];

    if(mpi_rank){

        MPI_Send(&Around_local,1,MPI_INT,0,0,wg_comm);
	MPI_Send(&Rough_local,1,MPI_INT,0,1,wg_comm);
        if(Around_local){
            MPI_Send(AM_local,N3D*Around_local,MPI_DOUBLE,0,2,wg_comm);
        }
    }else{

        Around_vec[0]    = Around_local;
        Around_Atom      = Around_local;
	Rough_Atom       = Rough_local;

        for(i_proc=1;i_proc<max_rank;i_proc++){
            MPI_Recv(&Around_vec[i_proc],1,MPI_INT,i_proc,0,wg_comm,&status);
	    MPI_Recv(&Rough_local,1,MPI_INT,i_proc,1,wg_comm,&status);
            Around_Atom  = Around_Atom+Around_vec[i_proc];
	    Rough_Atom   = Rough_Atom+Rough_local;
        }

        Around_Matrix = new double[max(N3D*(Around_Atom+2*Rough_Atom),1)];

        AM_size       = N3D*Around_vec[0];
        c_dcopy(AM_size,AM_local,1,Around_Matrix,1);

        for(i_proc=1;i_proc<max_rank;i_proc++){

            if(Around_vec[i_proc]){
                MPI_Recv(&Around_Matrix[AM_size],N3D*Around_vec[i_proc],MPI_DOUBLE,i_proc,2,\
                         wg_comm,&status);
                AM_size   = AM_size+N3D*Around_vec[i_proc];
            }
        }
        
        delete[] Around_vec;

    }
}

/************************************************************************************************/

void WireGenerator::check_periodicity(WireStructure *nanowire)
{
  
    if(nanowire->periodic_system){

        int IA,IB,IS,IN;
	int NA;
	int index;
	int atom_shift_pos;
	double Lperiodic;
	double *LM;
     
        init_variables(nanowire);

	//determine size of the structure and number of added atoms

	Lperiodic = 0.0;
	for(IS=0;IS<NSlab-1;IS++){
	    Lperiodic = Lperiodic+(Layer_Matrix[SLM*Smin[IS+1]]-Layer_Matrix[SLM*Smin[IS]]);
	}

	Lperiodic      = Lperiodic/(NSlab-1.0)*NSlab;
	
	NA             = Smax[NSlab-1]-Smin[0]+Smax[0]-Smin[0]+Smax[NSlab-1]-Smin[NSlab-1]+3;
	
	atom_shift_pos = Smax[NSlab-1]-Smin[NSlab-1]+1;
 
	for(IA=0;IA<No_Atom;IA++){
	    for(IB=0;IB<NB;IB++){
	        if(abs(Layer_Matrix[SLM*IA+4+IB])>tollim){
		    Layer_Matrix[SLM*IA+4+IB] = Layer_Matrix[SLM*IA+4+IB]+atom_shift_pos;
		}
	    }
	}

	//copy last unit cell to the beginning and first unit cell to the end

        LM = new double[NA*SLM];

	init_var(LM,NA*SLM);

	c_dcopy(No_Atom*SLM,Layer_Matrix,1,&LM[SLM*(Smax[NSlab-1]-Smin[NSlab-1]+1)],1);

	for(IA=Smin[0];IA<=Smax[0];IA++){
	    index = Smax[NSlab-1]-Smin[NSlab-1]+1+Smax[NSlab-1]-Smin[0]+1+IA-Smin[0];
	    c_dcopy(SLM-NB,&Layer_Matrix[SLM*(IA-Smin[0])],1,&LM[SLM*index],1);
	    LM[SLM*index] = LM[SLM*index]+Lperiodic;
	}

	for(IA=Smin[NSlab-1];IA<=Smax[NSlab-1];IA++){
	    index = IA-Smin[NSlab-1];
	    c_dcopy(SLM-NB,&Layer_Matrix[SLM*(IA-Smin[0])],1,&LM[SLM*index],1);
	    LM[SLM*index] = LM[SLM*index]-Lperiodic;
	}

	//make the missing connections in Layer_Matrix

	No_Atom = NA;
	
	delete[] Layer_Matrix;
	Layer_Matrix = new double[No_Atom*SLM];

	c_dcopy(SLM*No_Atom,LM,1,Layer_Matrix,1);

	int neighbor,no_neighbor;
	int ind1               = Smax[NSlab-1]-Smin[NSlab-1]+1;
	int ind2               = ind1+Smax[0]-Smin[0]+1;
	int ind3               = ind1+Smin[NSlab-1];
	int ind4               = ind3+Smax[NSlab-1]-Smin[NSlab-1]+1;
	int at_boundary_min[2] = {ind1,ind3};
	int at_boundary_max[2] = {ind2,ind4};
	int se_boundary_min[2] = {0,ind4};
	int se_boundary_max[2] = {ind1,No_Atom};
	double deformation;
	double at_pos[3],neigh_pos[3];

	deformation  = nanowire->max_bond_deformation*\
	  c_dnrm2(N3D,unit_cell->type_original[0]->bond[0]->vec,1);

	for(IS=0;IS<2;IS++){
	
	    for(IA=at_boundary_min[IS];IA<at_boundary_max[IS];IA++){

	        for(IB=0;IB<NB;IB++){

		    if(abs(Layer_Matrix[IA*SLM+N3D+1+IB])<tollim){

		        c_dcopy(N3D,&Layer_Matrix[IA*SLM],1,at_pos,1);
			c_daxpy(N3D,1.0,unit_cell->type[Round(Layer_Matrix[SLM*IA+N3D]-1)%NA]-> \
				bond[IB]->vec,1,at_pos,1);
			no_neighbor = 0;
                
			for(IN=se_boundary_min[IS];IN<se_boundary_max[IS];IN++){

			    c_dcopy(N3D,at_pos,1,neigh_pos,1);
			    c_daxpy(N3D,-1.0,&Layer_Matrix[SLM*IN],1,neigh_pos,1);
                        
			    if(c_dnrm2(N3D,neigh_pos,1)<=deformation){
			        no_neighbor = no_neighbor+1;
				if(no_neighbor>1){
				    printf("Please reduce max_bond_deformation\n");
				    exit(0);
				}else{
				    neighbor = IN;
				}    
			    }
			}
                
			if(no_neighbor==1){
			    Layer_Matrix[SLM*IA+N3D+1+IB] = neighbor+1.0;
			}
		    }
		}
	    }
	}

	delete_variables();

	//update the layer indices Lmin and Lmax

	XYZPOS* LM_xyz  = new XYZPOS[No_Atom];
        int *xyz_index  = new int[No_Atom];
        int *Lmin_local = new int[No_Atom];
        int *Lmax_local = new int[No_Atom];
   
        for(IA=0;IA<No_Atom;IA++){
	    LM_xyz[IA].x     = Layer_Matrix[SLM*IA];
            LM_xyz[IA].y     = Layer_Matrix[SLM*IA+1];
            LM_xyz[IA].z     = Layer_Matrix[SLM*IA+2];
            LM_xyz[IA].index = IA;
        }
        sort_xyz(LM_xyz,xyz_index,&NLayer,Lmin_local,Lmax_local);

	delete[] Lmin;
	delete[] Lmax;

	Lmin = new int[NLayer];
	Lmax = new int[NLayer];

	c_icopy(NLayer,Lmin_local,1,Lmin,1);
	c_icopy(NLayer,Lmax_local,1,Lmax,1);
	
	delete[] LM;
	delete[] xyz_index;
	delete[] Lmin_local;
	delete[] Lmax_local;
	delete LM_xyz;
    }
}

/************************************************************************************************/

void WireGenerator::update_periodicity(WireStructure *nanowire)
{

    if(nanowire->periodic_system){

        int IA,IL,IS;
	int NCH;
	int atom_shift_pos;
	double Lperiodic;
	double *LM;

	//determine size of the structure and number of added atoms

	Lperiodic = 0.0;
	for(IS=0;IS<NSlab-3;IS++){
	    Lperiodic = Lperiodic+(Layer_Matrix[SLM*Smin[IS+1]]-Layer_Matrix[SLM*Smin[IS]]);
	}

	Lperiodic      = Lperiodic/(NSlab-3.0)*(NSlab-2.0);

	atom_shift_pos = Smax[NSlab-1]-Smin[NSlab-1]+1;

	//shift atom coordinates

        LM = new double[No_Atom*SLM];

	c_dcopy(No_Atom*SLM,Layer_Matrix,1,LM,1);

	if(nanowire->update_at){

	    for(IA=0;IA<(Smax[NSlab-2]-Smin[1]+1);IA++){
	        c_dcopy(N3D+1,&LM[SLM*IA],1,&Layer_Matrix[SLM*(IA+atom_shift_pos)],1);
	    }

	    for(IA=Smin[0];IA<=Smax[0];IA++){
	        c_dcopy(N3D+1,&LM[SLM*IA],1,&Layer_Matrix[SLM*(Smin[NSlab-1]+IA-Smin[0])],1);
		Layer_Matrix[SLM*(Smin[NSlab-1]+IA-Smin[0])] = \
		  Layer_Matrix[SLM*(Smin[NSlab-1]+IA-Smin[0])]+Lperiodic;
	    }

	    for(IA=Smin[NSlab-3];IA<=Smax[NSlab-3];IA++){
	        c_dcopy(SLM-NB,&LM[SLM*IA],1,&Layer_Matrix[SLM*(IA-Smin[NSlab-3])],1);
		Layer_Matrix[SLM*(IA-Smin[NSlab-3])] = Layer_Matrix[SLM*(IA-Smin[NSlab-3])]-Lperiodic;
	    }
	}

	atom_shift_pos = (Smax[NSlab-2]-Smin[1]+1)-(Smax[0]-Smin[0]+1);
	for(IA=Smin[0];IA<=Smax[0];IA++){
	    index_channel[IA] = -1;
	    ch_conv[IA]       = IA+atom_shift_pos;
	}

	NCH = 0;
	for(IA=Smin[1];IA<=Smax[NSlab-2];IA++){
	    ch_pos[NCH]    = IA;
	    inv_ch_pos[IA] = NCH;
	    ch_conv[IA]    = NCH;
	    NCH++;
	}

	atom_shift_pos = (Smax[NSlab-2]-Smin[1]+1)+(Smax[0]-Smin[0]+1);
	for(IA=Smin[NSlab-1];IA<=Smax[NSlab-1];IA++){
	    index_channel[IA] = 1;
	    ch_conv[IA]       = IA-atom_shift_pos;
	}

	Channel_tot    = No_Atom;
        No_Atom        = NCH;

	NSlab          = NSlab-2;
	NLayer         = NLayer-2*layer_per_slab;

	atom_shift_pos = Smax[0]-Smin[0]+1;

	for(IS=0;IS<NSlab;IS++){
	    Smin[IS] = Smin[IS+1]-atom_shift_pos;
	    Smax[IS] = Smax[IS+1]-atom_shift_pos;
	}

	for(IL=0;IL<NLayer;IL++){
	    Lmin[IL] = Lmin[IL+layer_per_slab]-atom_shift_pos;
	    Lmax[IL] = Lmax[IL+layer_per_slab]-atom_shift_pos;
	}

	delete[] LM;
    }
}

/************************************************************************************************/

void WireGenerator::hydrogen_passivation(WireStructure *nanowire)
{

    //to test the model: TED 52, 1097 (2005) x=[100] y=[011] z=[0-11]

    if(nanowire->hydrogen_passivation){

        int No_hydro,NH;
	int IA,IB,IL;
	int atype;
	int ext_nlayer;
	int *atom_shift_pos;
	double *LM;
	double neigh_pos[3];
     
        init_variables(nanowire);
	cut_boundary_slab(nanowire);
	delete_variables();

	No_hydro = 0;

	for(IA=0;IA<No_Atom;IA++){
	    for(IB=0;IB<NB;IB++){
	        if(abs(Layer_Matrix[IA*SLM+4+IB])<tollim){
		    No_hydro++;
		}
	    }
	}

        LM             = new double[(No_Atom+No_hydro)*SLM];
	atom_shift_pos = new int[No_Atom+1];

	init_var(LM,(No_Atom+No_hydro)*SLM);
	init_var(atom_shift_pos,No_Atom+1);

	//The H atoms are placed at a distance bond_length/4.0 from the semiconducting atoms

	NH = 0;
	for(IA=0;IA<No_Atom;IA++){
	  
	    c_dcopy(SLM,&Layer_Matrix[IA*SLM],1,&LM[IA*SLM],1);

	    atype                = Round(Layer_Matrix[IA*SLM+3]-1)%NA;

	    atom_shift_pos[IA+1] = atom_shift_pos[IA];

	    for(IB=0;IB<NB;IB++){

	        if(abs(Layer_Matrix[IA*SLM+4+IB])<tollim){
		  
		    c_dcopy(N3D,&Layer_Matrix[IA*SLM],1,neigh_pos,1);
		    c_daxpy(N3D,0.25,unit_cell->type[atype]->bond[IB]->vec,1,neigh_pos,1);

		    c_dcopy(N3D,neigh_pos,1,&LM[(No_Atom+NH)*SLM],1);

		    LM[(No_Atom+NH)*SLM+3]    = 4-(atype)%2;
		    LM[(No_Atom+NH)*SLM+4+IB] = IA+1;
		    LM[IA*SLM+4+IB]           = No_Atom+NH+1;
		    atom_shift_pos[IA+1]      = atom_shift_pos[IA+1]+1;

		    NH++;

		}
	    }

	}

	No_Atom = No_Atom+No_hydro;

        XYZPOS* LM_xyz  = new XYZPOS[No_Atom];
        int *xyz_index  = new int[No_Atom];
	int *at_index   = new int[No_Atom];
        int *Lmin_local = new int[No_Atom];
        int *Lmax_local = new int[No_Atom];
   
        for(IA=0;IA<No_Atom;IA++){
	    LM_xyz[IA].x     = LM[SLM*IA];
            LM_xyz[IA].y     = LM[SLM*IA+1];
            LM_xyz[IA].z     = LM[SLM*IA+2];
            LM_xyz[IA].index = IA;
	    at_index[IA]     = IA;
        }
        sort_xyz(LM_xyz,xyz_index,&ext_nlayer,Lmin_local,Lmax_local);
        dreshape(No_Atom,SLM,LM,xyz_index);
       
	for(IA=0;IA<No_Atom;IA++){
	    at_index[xyz_index[IA]] = IA;
	}

	for(IA=0;IA<No_Atom;IA++){
	    for(IB=0;IB<NB;IB++){
	        if(LM[IA*SLM+4+IB]>tollim){
		    LM[IA*SLM+4+IB] = at_index[Round(LM[IA*SLM+4+IB])-1]+1;
		}
		if(abs(LM[IA*SLM+4+IB])<tollim){
		    LM[IA*SLM+4+IB] = -2.0; 
		}
	    }
	}

	for(IL=0;IL<NLayer;IL++){
	    Lmin[IL] = Lmin[IL]+atom_shift_pos[Lmin[IL]];
	    Lmax[IL] = Lmax[IL]+atom_shift_pos[Lmax[IL]+1];
	}

	delete[] Layer_Matrix;
	Layer_Matrix = new double[No_Atom*SLM];

	c_dcopy(SLM*No_Atom,LM,1,Layer_Matrix,1);

        delete LM_xyz;
	delete[] atom_shift_pos;
        delete[] xyz_index;
        delete[] Lmin_local;
        delete[] Lmax_local;
	delete[] at_index;
	delete[] LM;
    }
  
}

/************************************************************************************************/

void WireGenerator::roll_cnt()
{

    int IA,IB,IN,neigh,IA_shift;
    int update_boundary,IAs,INs,indb;
    double angle;
    double at_pos[3],neigh_pos[3];
    double yref = Layer_Matrix[1];
    double R0   = cnt_radius;
    double Lx   = get_length(0);
    
    for(IA=0;IA<No_Atom;IA++){
        for(IB=0;IB<NB;IB++){
	    if(Layer_Matrix[IA*SLM+4+IB]>0){
	        neigh = Layer_Matrix[IA*SLM+4+IB]-1;
		if(abs(Layer_Matrix[neigh*SLM+1]-yref)>tollim){
		    angle = (Layer_Matrix[neigh*SLM+1]-yref)/cnt_radius;
		    R0    = cnt_radius*angle/(2.0*sin(angle/2.0));
		    IB    = NB+1;
		    IA    = No_Atom+1;
		}
	    }
	}
    }

    IA_shift = 0;
    for(IA=0;IA<No_Atom;IA++){
        if(abs(Layer_Matrix[IA*SLM]-Layer_Matrix[0])<Lx){
	    IA_shift++;
	}else{
	    break;
	}
    }
   
    for(IA=0;IA<No_Atom;IA++){
        for(IB=0;IB<NB;IB++){
	    if(Layer_Matrix[IA*SLM+4+IB]==0){
	        c_dcopy(N3D,&Layer_Matrix[IA*SLM],1,at_pos,1);
		at_pos[1] = at_pos[1]+2.0*PI*cnt_radius;
		c_daxpy(N3D,1.0,unit_cell->type[Round(Layer_Matrix[SLM*IA+3]-1)%NA]->\
			bond[IB]->vec,1,at_pos,1);
		for(IN=max(IA-IA_shift,0);IN<min(IA+IA_shift,No_Atom);IN++){
		    c_dcopy(N3D,at_pos,1,neigh_pos,1);
		    c_daxpy(N3D,-1.0,&Layer_Matrix[SLM*IN],1,neigh_pos,1);
		    if(c_dnrm2(N3D,neigh_pos,1)<=tollim){
		        Layer_Matrix[IA*SLM+4+IB] = IN+1;
			Layer_Matrix[IN*SLM+4+IB] = IA+1;
			update_boundary           = 0;
			if(IA<Smin[1]){
			    indb = 0;
			    if(IN<Smin[1]){
			        update_boundary = 1;
				IAs             = IA;
				INs             = IN;
			    }else{
			        update_boundary = 2;
				IAs             = IA;
				INs             = IN-Smin[1];
			    }
			}
			if((IA>=Smin[1])&&(IA<Smin[2])){
			    indb = 0;
			    if(IN<Smin[1]){
			        update_boundary = 2;
				IAs             = IN; 
				INs             = IA-Smin[1];
			    }
			}
			if((IA>=Smin[NSlab-2])&&(IA<Smin[NSlab-1])){
			    indb = 1;
			    if(IN>=Smin[NSlab-1]){
			        update_boundary = 2;
				IAs             = IA-Smin[NSlab-2];
				INs             = IN-Smin[NSlab-1];
			    }
			}
			if(IA>=Smin[NSlab-1]){
			    indb = 1;
			    if(IN<Smin[NSlab-1]){
			        update_boundary = 2;
				IAs             = IN-Smin[NSlab-2];
				INs             = IA-Smin[NSlab-1]; 
			    }else{
			        update_boundary = 1;
				IAs             = IA-Smin[NSlab-1];
				INs             = IN-Smin[NSlab-1];
			    }
			}
			if(update_boundary==1){
			    Boundary[indb]->UNN[IAs*SLM+4+IB]         = INs+1;
			    Boundary[indb]->UNNneigh[IAs*(NB+1)+1+IB] = \
			      Layer_Matrix[IN*SLM+N3D];
			    Boundary[indb]->UNN[INs*SLM+4+IB]         = IAs+1;
			    Boundary[indb]->UNNneigh[INs*(NB+1)+1+IB] = \
				  Layer_Matrix[IA*SLM+N3D];
			}
			if(update_boundary==2){
			    Boundary[indb]->UNNp[IAs*SLM+4+IB]        = INs+1;
			    Boundary[indb]->UNN[IAs*SLM+4+IB]         = -1;
			    Boundary[indb]->UNNneigh[IAs*(NB+1)+1+IB] = \
			      Layer_Matrix[IN*SLM+N3D];
			    Boundary[indb]->UNNm[INs*SLM+4+IB]        = IAs+1;
			    Boundary[indb]->UNN[INs*SLM+4+IB]         = -1;
			    Boundary[indb]->UNNneigh[INs*(NB+1)+1+IB] = \
				  Layer_Matrix[IA*SLM+N3D];
			}
			break;
		    }
		}
	    }
	}
    }
    
    for(IA=0;IA<No_Atom;IA++){
        angle = (Layer_Matrix[IA*SLM+1]-yref)/cnt_radius;
	Layer_Matrix[IA*SLM+1] = R0*cos(angle);
	Layer_Matrix[IA*SLM+2] = R0*sin(angle);
	if(IA<Smin[1]){
	   Boundary[0]->UNN[IA*SLM+1]  = R0*cos(angle);
	   Boundary[0]->UNN[IA*SLM+2]  = R0*sin(angle);
	   Boundary[0]->UNNm[IA*SLM+1] = R0*cos(angle);
	   Boundary[0]->UNNm[IA*SLM+2] = R0*sin(angle);
	   Boundary[0]->UNNp[IA*SLM+1] = R0*cos(angle);
	   Boundary[0]->UNNp[IA*SLM+2] = R0*sin(angle);
	}
	if(IA>=Smin[NSlab-1]){
	   Boundary[1]->UNN[(IA-Smin[NSlab-1])*SLM+1]  = R0*cos(angle);
	   Boundary[1]->UNN[(IA-Smin[NSlab-1])*SLM+2]  = R0*sin(angle);
	   Boundary[1]->UNNm[(IA-Smin[NSlab-1])*SLM+1] = R0*cos(angle);
	   Boundary[1]->UNNm[(IA-Smin[NSlab-1])*SLM+2] = R0*sin(angle);
	   Boundary[1]->UNNp[(IA-Smin[NSlab-1])*SLM+1] = R0*cos(angle);
	   Boundary[1]->UNNp[(IA-Smin[NSlab-1])*SLM+2] = R0*sin(angle);
	}
    }
    
}

/************************************************************************************************/

void WireGenerator::make_bulk_layer_matrix(WireStructure *nanowire,int *no_orb)
{

    int IA,IB;
        
    Layer_Matrix = new double[NA*SLM];
    orb_per_at   = new int[NA+1];

    Channel_tot  = NA;

    init_var(Layer_Matrix,NA*SLM);
    init_var(orb_per_at,NA+1);

    orb_per_at[0] = 0;

    for(IA=0;IA<NA;IA++){

        c_dcopy(N3D,unit_cell->atom[IA]->coord,1,&Layer_Matrix[SLM*IA],1);
        
        Layer_Matrix[SLM*IA+3] = NA*(nanowire->bulk_mat_id-1)+IA+1;
        orb_per_at[IA+1]       = orb_per_at[IA]+no_orb[IA];
            
        for(IB=0;IB<NB;IB++){
                
            Layer_Matrix[SLM*IA+4+IB] = unit_cell->type[IA]->bond[IB]->neigh;        
                
        }    
    }
        
}

/************************************************************************************************/

void WireGenerator::replicate_unit_cell(int replicate_unit_cell,WireStructure *nanowire)
{

    if(replicate_unit_cell){

        int IA,IB,IC,IS,IL,IN;
        int naps,lps,ns,na_uc,*p;
	int no_neighbor,neighbor;
	int Lshift;
	double x0,cw,*uc_matrix;
	double at_pos[3],neigh_pos[3];
	double deformation;

	deformation  = nanowire->max_bond_deformation*\
	  c_dnrm2(N3D,unit_cell->type_original[0]->bond[0]->vec,1);

        naps         = make_slab(nanowire);
	cw           = Layer_Matrix[SLM*naps]-Layer_Matrix[0];
	p            = find(Lmax,Lmax+NLayer,naps-1);
	lps          = (p-Lmax)+1;
	ns           = (int)(NLayer/lps);
	NLayer       = ns*lps;

	na_uc        = (Lmax[lps*replicate_unit_cell-1]-Lmin[lps*(replicate_unit_cell-1)]+1);

	uc_matrix    = new double[SLM*na_uc];

	init_var(uc_matrix,SLM*na_uc);
	
	for(IC=0;IC<(N3D+1);IC++){
	    c_dcopy(na_uc,&Layer_Matrix[SLM*Lmin[lps*(replicate_unit_cell-1)]+IC],SLM,\
		    &uc_matrix[IC],SLM);
	}

	x0           = uc_matrix[0];
      
	for(IA=0;IA<na_uc;IA++){
	    uc_matrix[SLM*IA] = uc_matrix[SLM*IA]-x0;
	}

	delete[] Layer_Matrix;

	No_Atom      = ns*na_uc;
	Layer_Matrix = new double[No_Atom*SLM];

	init_var(Layer_Matrix,No_Atom*SLM);

	Lshift       = Lmin[lps*(replicate_unit_cell-1)];

	for(IS=0;IS<ns;IS++){
	    c_dcopy(SLM*na_uc,uc_matrix,1,&Layer_Matrix[IS*na_uc*SLM],1);
	    for(IA=0;IA<na_uc;IA++){
	        Layer_Matrix[(IS*na_uc+IA)*SLM] = Layer_Matrix[(IS*na_uc+IA)*SLM]+IS*cw; 
	    }
	    for(IL=0;IL<lps;IL++){
	        Lmin[IS*lps+IL] = Lmin[lps*(replicate_unit_cell-1)+IL]-Lshift+IS*na_uc;
		Lmax[IS*lps+IL] = Lmax[lps*(replicate_unit_cell-1)+IL]-Lshift+IS*na_uc;
	    }
	    if(IS==(replicate_unit_cell-1)){
	        Lshift = Lmin[lps*(replicate_unit_cell-1)];
	    }
	}

	for(IL=0;IL<NLayer;IL++){

	    for(IA=Lmin[IL];IA<=Lmax[IL];IA++){

	        for(IB=0;IB<NB;IB++){

		    if(abs(Layer_Matrix[IA*SLM+N3D+1+IB])<tollim){

		        c_dcopy(N3D,&Layer_Matrix[IA*SLM],1,at_pos,1);
			c_daxpy(N3D,1.0,unit_cell->type[Round(Layer_Matrix[SLM*IA+N3D]-1)%NA]->	\
				bond[IB]->vec,1,at_pos,1);
			no_neighbor = 0;
                
			for(IN=IA+1;IN<=Lmax[min(IL+lps,NLayer-1)];IN++){

			    c_dcopy(N3D,at_pos,1,neigh_pos,1);
			    c_daxpy(N3D,-1.0,&Layer_Matrix[SLM*IN],1,neigh_pos,1);
                        
			    if(c_dnrm2(N3D,neigh_pos,1)<=deformation){
			        no_neighbor = no_neighbor+1;
				if(no_neighbor>1){
				    printf("Please reduce max_bond_deformation\n");
				    exit(0);
				}else{
				    neighbor = IN;
				}
			    }    
			}
                
			if(no_neighbor==1){
			    Layer_Matrix[SLM*IA+N3D+1+IB]       = neighbor+1.0;
			    Layer_Matrix[SLM*neighbor+N3D+1+IB] = IA+1.0;
			}
		    }
		}
	    }
	}

	delete[] uc_matrix;
    }

}

/************************************************************************************************/

void WireGenerator::read_atom_position(int read_atom_pos,WireStructure *nanowire)
{

    if(read_atom_pos){

        int IA,IB,IP,IL,IN;
	int naps,lps,*p;
	int no_neighbor,neighbor;
	double at_pos[3],neigh_pos[3];

	double deformation  = nanowire->max_bond_deformation*\
	  c_dnrm2(N3D,unit_cell->type_original[0]->bond[0]->vec,1);

	FILE *F = fopen(nanowire->at_file,"r");

	fscanf(F,"%i",&No_Atom);

	delete[] Layer_Matrix;

	Layer_Matrix = new double[No_Atom*SLM];

	init_var(Layer_Matrix,No_Atom*SLM);
	for(IA=0;IA<No_Atom;IA++){
	    for(IP=0;IP<(N3D+1);IP++){
	        fscanf(F,"%lg",&Layer_Matrix[SLM*IA+IP]);
	    }
	}
	fclose(F);

	delete[] Lmin;
	delete[] Lmax;

	XYZPOS* LM_xyz  = new XYZPOS[No_Atom];
        int *xyz_index  = new int[No_Atom];
        int *Lmin_local = new int[No_Atom];
        int *Lmax_local = new int[No_Atom];

        for(IA=0;IA<No_Atom;IA++){
            LM_xyz[IA].x     = Layer_Matrix[SLM*IA];
            LM_xyz[IA].y     = Layer_Matrix[SLM*IA+1];
            LM_xyz[IA].z     = Layer_Matrix[SLM*IA+2];
            LM_xyz[IA].index = IA;
        }

        sort_xyz(LM_xyz,xyz_index,&NLayer,Lmin_local,Lmax_local);
        dreshape(No_Atom,SLM,Layer_Matrix,xyz_index);

        Lmin = new int[NLayer];
        Lmax = new int[NLayer];

        icopy(NLayer,Lmin_local,Lmin);
        icopy(NLayer,Lmax_local,Lmax);

        delete LM_xyz;
        delete[] xyz_index;
        delete[] Lmin_local;
        delete[] Lmax_local;

	naps = make_slab(nanowire);
	p    = find(Lmax,Lmax+NLayer,naps-1);
	lps  = (p-Lmax)+1;

	for(IL=0;IL<NLayer;IL++){

	    for(IA=Lmin[IL];IA<=Lmax[IL];IA++){

	        for(IB=0;IB<NB;IB++){

		    if(abs(Layer_Matrix[IA*SLM+N3D+1+IB])<tollim){

		        c_dcopy(N3D,&Layer_Matrix[IA*SLM],1,at_pos,1);
			c_daxpy(N3D,1.0,unit_cell->type[Round(Layer_Matrix[SLM*IA+N3D]-1)%NA]->	\
				bond[IB]->vec,1,at_pos,1);
			no_neighbor = 0;
                
			for(IN=IA+1;IN<=Lmax[min(IL+lps,NLayer-1)];IN++){

			    c_dcopy(N3D,at_pos,1,neigh_pos,1);
			    c_daxpy(N3D,-1.0,&Layer_Matrix[SLM*IN],1,neigh_pos,1);
                        
			    if(c_dnrm2(N3D,neigh_pos,1)<=deformation){
			        no_neighbor = no_neighbor+1;
				if(no_neighbor>1){
				    printf("Please reduce max_bond_deformation\n");
				    exit(0);
				}else{
				    neighbor = IN;
				}
			    }    
			}
                
			if(no_neighbor==1){
			    Layer_Matrix[SLM*IA+N3D+1+IB]       = neighbor+1.0;
			    Layer_Matrix[SLM*neighbor+N3D+1+IB] = IA+1.0;
			}
		    }
		}
	    }
	}
    }
}

/************************************************************************************************/

void WireGenerator::get_orbital_per_atom(int *no_orb)
{

    int IA;

    orb_per_at = new int[No_Atom+1];

    orb_per_at[0] = 0;
    for(IA=0;IA<No_Atom;IA++){
        orb_per_at[IA+1] = orb_per_at[IA]+no_orb[Round(Layer_Matrix[SLM*ch_pos[IA]+3])-1];
    }

}

/************************************************************************************************/

void WireGenerator::get_mid_gap_energy(Material *material)
{
    int IA,IB,neigh_pos;
    int typeA,typeB;
    double no_neigh;

    EMidGap = new double[No_Atom];
    EGap    = new double[No_Atom];

    init_var(EMidGap,No_Atom);
    init_var(EGap,No_Atom);
    
    EVmaxL = 0.0;
    EVmaxR = 0.0;
    ECminL = 0.0;
    ECminR = 0.0;
    
    EMGL   = 0.0;
    EMGR   = 0.0;

    for(IA=0;IA<No_Atom;IA++){
        typeA    = Round(Layer_Matrix[SLM*ch_pos[IA]+3])-1;
        no_neigh = 0.0;
        for(IB=0;IB<NB;IB++){
	    if(Layer_Matrix[SLM*ch_pos[IA]+4+IB]>0){
	        neigh_pos   = Round(Layer_Matrix[SLM*ch_pos[IA]+4+IB])-1;
	        typeB       = Round(Layer_Matrix[SLM*neigh_pos+3])-1;
		EMidGap[IA] = EMidGap[IA]+material->mid_gap_energy[typeA][typeB];
		EGap[IA]    = EGap[IA]+material->band_gap_table[typeA][typeB];
		no_neigh    = no_neigh+1.0;
	    }
	}
	if(no_neigh>0){
	    EMidGap[IA] = EMidGap[IA]/no_neigh;
	    EGap[IA]    = EGap[IA]/no_neigh;
	}else{
	    typeB       = (typeA+1)%2;
	    EMidGap[IA] = material->mid_gap_energy[typeA][typeB];
	    EGap[IA]    = material->band_gap_table[typeA][typeB];
	}
	if((IA>=Smin[0])&&(IA<=Smax[0])){
	    EVmaxL = EVmaxL+(EMidGap[IA]-EGap[IA]/2.0);
	    ECminL = ECminL+(EMidGap[IA]+EGap[IA]/2.0);
	    EMGL   = EMGL+EMidGap[IA];
	}
	if((IA>=Smin[NSlab-1])&&(IA<=Smax[NSlab-1])){
	    EVmaxR = EVmaxR+(EMidGap[IA]-EGap[IA]/2.0);
	    ECminR = ECminR+(EMidGap[IA]+EGap[IA]/2.0);
	    EMGR   = EMGR+EMidGap[IA];
	}
    }

    EVmaxL = EVmaxL/(Smax[0]-Smin[0]+1.0);
    ECminL = ECminL/(Smax[0]-Smin[0]+1.0);
    EMGL   = EMGL/(Smax[0]-Smin[0]+1.0);

    EVmaxR = EVmaxR/(Smax[NSlab-1]-Smin[NSlab-1]+1.0);
    ECminR = ECminR/(Smax[NSlab-1]-Smin[NSlab-1]+1.0);
    EMGR   = EMGR/(Smax[NSlab-1]-Smin[NSlab-1]+1.0);
}

/************************************************************************************************/

void WireGenerator::generate_roughness(WireStructure *nanowire,int sample_id)
{

    if((nanowire->rough->on)&&(!mpi_rank)){
      
        int IA,IB,IF,JF;
	int IL,IN;
        int no_surface_atom,no_mat_face,no_face;
        int neigh,no_neighbor,atom_index;
	int *face,*face_info,*surface_info;
	int is_surface_atom;
	int seed;
	int shift;
	int iteration;
	int max_iteration = 10;
	double *roughness,*roughness_factor;
	double dist,p0[3],p1[3];

	if(!strcmp(nanowire->rough->type,"square")){
	    no_mat_face = 6;
	}else{
	    no_mat_face = 1;
	}

        no_surface_atom = 0;
        for(IA=0;IA<No_Atom;IA++){
	    if(Layer_Matrix[SLM*IA+3]>0){
	        for(IB=0;IB<NB;IB++){
	            if(Layer_Matrix[SLM*IA+4+IB]>0){
	                if(Layer_Matrix[SLM*(Round(Layer_Matrix[SLM*IA+4+IB])-1)+3]<0){
	                    no_surface_atom++;
	                    break;
	                }
	            }
	        }
            }
	}

        roughness        = new double[no_surface_atom];
	roughness_factor = new double[No_Atom];
	surface_info     = new int[2*no_surface_atom];
	face             = new int[2*no_mat_face*nanowire->no_ch_element];
	face_info        = new int[no_mat_face*nanowire->no_ch_element];

	init_var(face,2*no_mat_face*nanowire->no_ch_element);

	atom_index = 0;
        no_face    = 0;
        for(IA=0;IA<No_Atom;IA++){
	    if(Layer_Matrix[SLM*IA+3]>0){
	        for(IB=0;IB<NB;IB++){
	            if((Layer_Matrix[SLM*IA+4+IB]>0)){
	                if(Layer_Matrix[SLM*(Round(Layer_Matrix[SLM*IA+4+IB])-1)+3]<0){
		            get_face_info(&no_face,face,&surface_info[atom_index],IA,\
				          &Layer_Matrix[SLM*IA],no_surface_atom,no_mat_face,\
				          nanowire);
	                    atom_index++;
	                    break;
	                }
	            }
	        }
            }
	}

	for(IF=0;IF<no_face;IF++){
	    for(JF=0;JF<nanowire->no_ch_element*no_mat_face;JF++){
	        if((face[JF]-1)==IF){
	            face_info[IF] = face[nanowire->no_ch_element*no_mat_face+JF];
	            break;
	        }
	    }
	}	
	
	if(nanowire->rough->seed<0){
	    seed = (sample_id/100.0+1.0)*time(0);
	}else{
	    seed = nanowire->rough->seed;
	}

	srand(seed);
	generate_random_surface(roughness,surface_info,no_surface_atom,no_face,face_info,\
				nanowire);
	printf("The seed number for the rough surface is %i\n",seed);

	init_var(roughness_factor,No_Atom);

	atom_index = 0;
	for(IL=0;IL<NLayer;IL++){

	    for(IA=Lmin[IL];IA<=Lmax[IL];IA++){

	        is_surface_atom = 0;

	        if(Layer_Matrix[SLM*IA+3]>0){

		    for(IB=0;IB<NB;IB++){

		        if(Layer_Matrix[SLM*IA+4+IB]>0){

			    if(Layer_Matrix[SLM*(Round(Layer_Matrix[SLM*IA+4+IB])-1)+3]<0){
			        is_surface_atom = 1;
				break;
			    }
			}
		    }
		}
		
		if(is_surface_atom){

		    c_dcopy(N3D,&Layer_Matrix[SLM*IA],1,p0,1);

		    for(IN=Lmin[max(IL-1,0)];IN<=Lmax[min(IL+1,NLayer-1)];IN++){

		        c_dcopy(N3D,&Layer_Matrix[SLM*IN],1,p1,1);

			dist = sqrt(pow(p1[0]-p0[0],2.0)+pow(p1[1]-p0[1],2.0)+\
				    pow(p1[2]-p0[2],2.0));

			if(dist<=abs(roughness[atom_index])){

			    if(IN!=IA){
			        roughness_factor[IN] = roughness_factor[IN]+roughness[atom_index]/dist;
			    }else{
			        roughness_factor[IN] = roughness_factor[IN]+10.0*roughness[atom_index];
			    }
			}
		    }
		    atom_index++;
		}
	    }
	}

	shift = 2*(Lmax[0]-Lmin[0]+1);

	for(IA=max(surface_info[0]-shift,0);IA<min(surface_info[no_surface_atom-1]+shift,No_Atom);IA++){

	    if(Layer_Matrix[SLM*IA+3]<0){

	        if(roughness_factor[IA]>0){
		    Layer_Matrix[SLM*IA+3] = -Layer_Matrix[SLM*IA+3];
		}
	    }else{

	        if(roughness_factor[IA]<0){
		    Layer_Matrix[SLM*IA+3] = -Layer_Matrix[SLM*IA+3];
		}
	    }
	}

	//to avoid atoms with less than 1 connection;
	iteration = 0;
	while(iteration<max_iteration){
	    for(IA=max(surface_info[0]-shift,0);IA<min(surface_info[no_surface_atom-1]+shift,No_Atom);IA++){
	        no_neighbor = 0;
		if(Layer_Matrix[SLM*IA+3]>0){
		    for(IB=0;IB<NB;IB++){
		        if((Layer_Matrix[SLM*IA+4+IB]>0)){
			    neigh = Round(Layer_Matrix[SLM*IA+4+IB])-1;
			    if(Layer_Matrix[SLM*neigh+3]>0){
			        no_neighbor++;
			    }
			}
		    }
		    if(no_neighbor<=1){
		        Layer_Matrix[SLM*IA+3] = -Layer_Matrix[SLM*IA+3];
		    }
		}
	    }
	    iteration++;
	}

	delete[] roughness;
	delete[] roughness_factor;
	delete[] surface_info;
	delete[] face;
	delete[] face_info;
    }

    if(!mpi_rank){

        int IL,IA,IB,atom_index,no_rm_atom,start,stop,*position;

	position = new int[No_Atom];

        atom_index = 0;
	no_rm_atom = 0;
        for(IL=0;IL<NLayer;IL++){

	    start    = Lmin[IL];
	    stop     = Lmax[IL]+1;
	    Lmin[IL] = Lmin[IL]-no_rm_atom;

	    for(IA=start;IA<stop;IA++){

	        if(Layer_Matrix[SLM*IA+3]>0){

		    c_dcopy(SLM,&Layer_Matrix[SLM*IA],1,&Layer_Matrix[SLM*atom_index],1);
		    position[IA] = atom_index;
		    atom_index++;

	        }else{

		    //c_dcopy(N3D,&Layer_Matrix[SLM*IA],1,&Around_Matrix[N3D*Around_Atom],1);
		    position[IA] = -1;
		    no_rm_atom++;
		    No_Atom--;
		    //Around_Atom++;
		}
	    }

	    Lmax[IL] = Lmax[IL]-no_rm_atom;
        }

	for(IA=0;IA<No_Atom;IA++){

	    for(IB=0;IB<NB;IB++){

	        if(Layer_Matrix[SLM*IA+4+IB]>0){

	            atom_index = Round(Layer_Matrix[SLM*IA+4+IB])-1;

	            if(position[atom_index]!=-1){

		        Layer_Matrix[SLM*IA+4+IB] = (double)(position[atom_index]+1);

		    }else{

		        Layer_Matrix[SLM*IA+4+IB] = 0.0;

		    }
		}
	    }
	}

	delete[] position;
    }

    MPI_Bcast(&No_Atom,1,MPI_INT,0,wg_comm);
    MPI_Bcast(&Around_Atom,1,MPI_INT,0,wg_comm);
    MPI_Bcast(Layer_Matrix,SLM*No_Atom,MPI_DOUBLE,0,wg_comm);
    MPI_Bcast(Lmin,NLayer,MPI_INT,0,wg_comm);
    MPI_Bcast(Lmax,NLayer,MPI_INT,0,wg_comm);
    if(Around_Atom){
        MPI_Bcast(Around_Matrix,N3D*Around_Atom,MPI_DOUBLE,0,wg_comm);
    }

}

/************************************************************************************************/

void WireGenerator::get_face_info(int *no_face,int *face,int *surface_info,int atom_index,\
				  double *p_coord,int no_surface_atom,int no_mat_face,\
				  WireStructure *nanowire)
{
    int IM,IF;
    double max_dist = nanowire->a0*sqrt(3.0)/4.0;

    for(IM=0;IM<nanowire->no_ch_element;IM++){
        for(IF=0;IF<no_mat_face;IF++){
	    if(is_in_face(p_coord,nanowire->mat[IM],IF,nanowire->rough->type,max_dist)){
	        if(!face[IM*no_mat_face+IF]){
		    face[IM*no_mat_face+IF] = (*no_face)+1;
		    (*no_face)              = (*no_face)+1;
		}
		face[(IM+nanowire->no_ch_element)*no_mat_face+IF]++;
		surface_info[0]               = atom_index;
		surface_info[no_surface_atom] = face[IM*no_mat_face+IF]-1;
	        break;
	    }
        }
    }

}

/************************************************************************************************/

int WireGenerator::is_in_face(double *p0,MAT *mat,int index,char *rtype,double max_dist)
{
    int is_in_face = 0;

    if(!strcmp(rtype,"square")){

        double h[3],d;

        c_dcopy(N3D,p0,1,h,1);
        c_daxpy(N3D,-1.0,mat->p_ref[index],1,h,1);

        d          = c_ddot(N3D,mat->vec_dir[index],1,h,1);

        is_in_face = (fabs(d)<max_dist);
        
    }else{

        bool condition1 = p0[0]>=mat->p[0]->coord[0];
        bool condition2 = p0[0]<=mat->p[1]->coord[0];

        is_in_face = condition1&&condition2;
        
    }
    
    return is_in_face;
}

/************************************************************************************************/

void WireGenerator::generate_random_surface(double *roughness,int *surface_info,int no_surface_atom,\
					    int no_face,int *face_info,WireStructure *nanowire)
{
    int IF,I,J,index_i,index_j;
    double *ACVMatrix,*Sigma,*SR;
    double p0[3],p1[3],face_center[3];

    for(IF=0;IF<no_face;IF++){

        ACVMatrix = new double[face_info[IF]*face_info[IF]];
	Sigma     = new double[face_info[IF]];
	SR        = new double[face_info[IF]];

	get_face_center(face_center,IF,face_info[IF],no_surface_atom,surface_info);

	index_i = 0;
	for(I=0;I<no_surface_atom;I++){

	    if(surface_info[no_surface_atom+I]==IF){

	        index_j = index_i;
	        for(J=I;J<no_surface_atom;J++){

	            if(surface_info[no_surface_atom+J]==IF){

		        c_dcopy(N3D,&Layer_Matrix[SLM*surface_info[I]],1,&p0[0],1);
		        c_dcopy(N3D,&Layer_Matrix[SLM*surface_info[J]],1,&p1[0],1);

			get_acv_value(&ACVMatrix[index_i+index_j*face_info[IF]],p0,p1,\
				      face_center,nanowire->rough);
			
			index_j++;

	            }
	        }
		index_i++;
	    }
	}

	eig(ACVMatrix,Sigma,face_info[IF]);

	for(I=0;I<face_info[IF];I++){
	    Sigma[I] = sqrt(abs(Sigma[I]))*randn();
	} 

	c_dgemm('N','N',face_info[IF],1,face_info[IF],1.0,ACVMatrix,face_info[IF],\
		Sigma,face_info[IF],0.0,SR,face_info[IF]);

	index_i = 0;
	for(I=0;I<no_surface_atom;I++){

	    if(surface_info[no_surface_atom+I]==IF){

	        roughness[I] = SR[index_i];
	        index_i++;

	    }
	}

	delete[] ACVMatrix;
	delete[] Sigma;
	delete[] SR;
    }

}

/************************************************************************************************/

void WireGenerator::get_face_center(double *face_center,int face_index,int face_info,\
				    int no_surface_atom,int *surface_info)
{
    int IA;

    init_var(face_center,N3D);

    for(IA=0;IA<no_surface_atom;IA++){
        if(surface_info[no_surface_atom+IA]==face_index){
            face_center[0] = face_center[0]+Layer_Matrix[SLM*surface_info[IA]];
            face_center[1] = face_center[1]+Layer_Matrix[SLM*surface_info[IA]+1];
            face_center[2] = face_center[2]+Layer_Matrix[SLM*surface_info[IA]+2];
        }
    }

    c_dscal(N3D,1.0/face_info,face_center,1);

}

/************************************************************************************************/

void WireGenerator::get_acv_value(double *ACVM,double *p0,double *p1,double *face_center,\
				  Roughness *rough)
{

    double x;

    if(!strcmp(rough->type,"square")){

        x = sqrt(pow(p1[0]-p0[0],2.0)+pow(p1[1]-p0[1],2.0)+pow(p1[2]-p0[2],2.0));
        
    }else{

      double alpha,cos_alpha,v0[2],v1[2];

      v0[0] = p0[1]-face_center[1];
      v0[1] = p0[2]-face_center[2];
      v1[0] = p1[1]-face_center[1];
      v1[1] = p1[2]-face_center[2];
      
      cos_alpha = c_ddot(2,v0,1,v1,1)/(c_dnrm2(2,v0,1)*c_dnrm2(2,v1,1));
      cos_alpha = cos_alpha-fmod(cos_alpha,tollim);
      alpha     = acos(cos_alpha);

      x     = sqrt(pow(p1[0]-p0[0],2.0)+pow((c_dnrm2(2,v0,1)+c_dnrm2(2,v1,1))*alpha/2.0,2.0));
      
    }

    *ACVM = pow(rough->rms,2.0)*exp(-x*sqrt(2.0)/rough->Lms);

}

/************************************************************************************************/

void WireGenerator::eig(double *A,double *lambda,int N)
{
    int info;
    char jobz    = 'V';
    char uplo    = 'U';
    int lwork    = 20*N;
    double *work = new double[lwork];

    c_dsyev(jobz,uplo,N,A,N,lambda,work,lwork,&info);

    delete[] work;
}

/************************************************************************************************/

void WireGenerator::generate_alloy_disorder(WireStructure *nanowire)
{

    if((nanowire->alloy->on)&&(!mpi_rank)){

        int seed;
	int IA;
        double r_number;

	if(nanowire->alloy->seed<0){
	    seed = time(0);
	}else{
	    seed = nanowire->alloy->seed;
	}

	printf("The seed number for alloy disorder is %i\n",seed);

	srand(seed);

	for(IA=0;IA<No_Atom;IA++){

	    r_number = (double)rand()/RAND_MAX;

	    if(r_number>=nanowire->alloy->composition){
	        Layer_Matrix[IA*SLM+3] = Layer_Matrix[IA*SLM+3]+2;
	    }
	}

	for(IA=Smin[1];IA<=Smax[1];IA++){
	    Layer_Matrix[IA*SLM+3] = Layer_Matrix[(IA-Smin[1]+Smin[0])*SLM+3];
	}

	for(IA=Smin[NSlab-2];IA<=Smax[NSlab-2];IA++){
	    Layer_Matrix[IA*SLM+3] = Layer_Matrix[(IA-Smin[NSlab-2]+Smin[NSlab-1])*SLM+3];
	}
    }

    if(nanowire->alloy->on){
        MPI_Bcast(Layer_Matrix,SLM*No_Atom,MPI_DOUBLE,0,wg_comm);

        get_cell_width();
        update_boundary_slab(nanowire);
    }
}

/************************************************************************************************/

void WireGenerator::init_variables(WireStructure *nanowire)
{
    int IL,IA,IB,IS,*p;
    
    Neigh_Matrix   = new int[(NB+1)*No_Atom];
    
    NA_per_slice   = make_slab(nanowire);
    cell_width     = Layer_Matrix[SLM*NA_per_slice]-Layer_Matrix[0];
    p              = find(Lmax,Lmax+NLayer,NA_per_slice-1);
    layer_per_slab = (p-Lmax)+1;
    l_per_slab_ref = layer_per_slab;

    if(layer_per_slab == 2){
        NA_per_slice   = 2*NA_per_slice;
        cell_width     = 2.0*cell_width;
        layer_per_slab = 2*layer_per_slab;
	l_per_slab_ref = 2*l_per_slab_ref;
    }

    //NA_per_slice   = nanowire->NxFold*NA_per_slice;
    //cell_width     = nanowire->NxFold*cell_width;
    //layer_per_slab = nanowire->NxFold*layer_per_slab;

    if(nanowire->open_system||nanowire->periodic_system){
        NSlab  = (int)(NLayer/layer_per_slab);
	NLayer = NSlab*layer_per_slab;
    }else{
        NSlab  = (ceil)((double)NLayer/layer_per_slab);
    }
    No_Atom  = Lmax[NLayer-1]+1;

    Smin     = new int[NSlab];
    Smax     = new int[NSlab];
    Smin_tot = new int[NSlab];
    Smax_tot = new int[NSlab];

    for(IS=0;IS<NSlab;IS++){
        Smin[IS] = Lmin[IS*layer_per_slab];
	Smax[IS] = Lmax[min((IS+1)*layer_per_slab,NLayer)-1];
    }

    index_channel  = new int[No_Atom];
    ch_pos         = new int[No_Atom];
    ch_conv        = new int[No_Atom];
    inv_ch_pos     = new int[No_Atom];
    index_arch     = new int[max(Around_Atom,1)];
    arch_pos       = new int[max(Around_Atom,1)];
    arch_conv      = new int[max(Around_Atom,1)];
    neighbor_layer = new int[NLayer];
    Enn            = new BOUNDARY_ENN*[NLayer];

    for(IL=0;IL<NLayer;IL++){

        Enn[IL]           = new BOUNDARY_ENN();
        Enn[IL]->NA       = Lmax[IL]-Lmin[IL]+1;
        Enn[IL]->UNN      = new double[SLM*Enn[IL]->NA];
        Enn[IL]->UNNp     = new double[SLM*Enn[IL]->NA];
        Enn[IL]->UNNneigh = new int[(NB+1)*Enn[IL]->NA];

        neighbor_layer[IL] = 0;
        for(IA=Lmin[IL];IA<=Lmax[IL];IA++){
            for(IB=0;IB<NB;IB++){
                if(Layer_Matrix[SLM*IA+4+IB]>0){
                    p = find_if(Lmax,Lmax+NLayer,bind2nd(greater_equal<int>(),\
							 Layer_Matrix[SLM*IA+4+IB]-1));
                    if((p-Lmax-IL)>neighbor_layer[IL]){
                        neighbor_layer[IL] = (p-Lmax-IL);
                    }
                }
            }
        }
    }
    get_boundary_size(nanowire);
 
    int index[2]               = {0,NSlab-1};
    
    Boundary                   = new BOUNDARY_ENN*[2];
    index_boundary             = new int*[2];
    bound_pos                  = new int*[2];
    bound_conv                 = new int*[2];
    b_shift                    = new int[2];
    
    for(IB=0;IB<2;IB++){

        Boundary[IB]           = new BOUNDARY_ENN();
        Boundary[IB]->NA       = Smax[index[IB]]-Smin[index[IB]]+1;
	Boundary[IB]->NA_tot   = Boundary[IB]->NA;

        Boundary[IB]->UNN      = new double[SLM*(Smax[index[IB]]-Smin[index[IB]]+1)];
        Boundary[IB]->UNNp     = new double[SLM*(Smax[index[IB]]-Smin[index[IB]]+1)];
        Boundary[IB]->UNNm     = new double[SLM*(Smax[index[IB]]-Smin[index[IB]]+1)];
        Boundary[IB]->UNNneigh = new int[(NB+1)*(Smax[index[IB]]-Smin[index[IB]]+1)];
        index_boundary[IB]     = new int[Boundary[IB]->NA];
        bound_pos[IB]          = new int[Boundary[IB]->NA];
        bound_conv[IB]         = new int[Boundary[IB]->NA];

        b_shift[IB]            = Smin[index[IB]];
    }
 
    get_cell_width();
    
}

/************************************************************************************************/

void WireGenerator::get_cell_width()
{

    int IB,IS;
    int index[2]               = {0,NSlab-2};

    /*
    Boundary[0]->cell_width    = min_vec(&Layer_Matrix[SLM*Smin[1]],Smax[1]-Smin[1]+1,SLM)-\
                                 min_vec(&Layer_Matrix[SLM*Smin[0]],Smax[0]-Smin[0]+1,SLM);
    Boundary[1]->cell_width    = max_vec(&Layer_Matrix[SLM*Smin[NSlab-1]],Smax[NSlab-1]-Smin[NSlab-1]+1,SLM)-\
                                 max_vec(&Layer_Matrix[SLM*Smin[NSlab-2]],Smax[NSlab-2]-Smin[NSlab-2]+1,SLM);
    */

    if(NSlab>1){
        for(IB=0;IB<2;IB++){
            Boundary[IB]->cell_width = 0.0;
            for(IS=Smin[index[IB]+1];IS<=Smax[index[IB]+1];IS++){
                Boundary[IB]->cell_width = Boundary[IB]->cell_width+Layer_Matrix[SLM*IS]/\
                    (Smax[index[IB]+1]-Smin[index[IB]+1]+1);
            }
            for(IS=Smin[index[IB]];IS<=Smax[index[IB]];IS++){
                Boundary[IB]->cell_width = Boundary[IB]->cell_width-Layer_Matrix[SLM*IS]/\
                    (Smax[index[IB]]-Smin[index[IB]]+1);
            }
        }
    }else{
        for(IB=0;IB<2;IB++){
            Boundary[IB]->cell_width = cell_width;
        }
    }

    for(IB=0;IB<2;IB++){
        Boundary[IB]->cell_area = Boundary[IB]->NA*Vol_atom/(Boundary[IB]->cell_width*1e-27);
    }

}

/************************************************************************************************/

int WireGenerator::make_slab(WireStructure* nanowire)
{

    double Lx    = nanowire->NxFold*get_length(0);
    double x0    = min_vec(Layer_Matrix,No_Atom,SLM);
    double xN    = x0+Lx;
    int At_local = 0;
    int IA;
    
    for(IA=0;IA<No_Atom;IA++){
        if(Layer_Matrix[SLM*IA]<xN-tollim){
	    At_local++;
	}else{
	    break;
	}
    }

    /*
    int Ny,Nz,At_local=0,IY,IZ,IS,IA,condition,uc_axis[2];
    int no_element,no_ch_element,hymi, hyma, hzmi, hzma;
    double shifty[2]={0},shiftyz[2],at_pos[2];

    no_element    = nanowire->no_element;
    no_ch_element = nanowire->no_ch_element;

    hymi = INF;
    hzmi = INF;
    hyma = -INF;
    hzma = -INF;
    find_2D_boundary(&hymi,&hyma,&hzmi,&hzma,nanowire);
    Ny = hyma-hymi+1;
    Nz = hzma-hzmi+1;

    check_uc_axis(uc_axis);
    
    c_daxpy(N2D,hymi,&unit_cell->axis[uc_axis[0]]->vec[1],1,shifty,1);
    
    for(IY=0;IY<Ny;IY++){
        
        c_dcopy(N2D,shifty,1,shiftyz,1);
        c_daxpy(N2D,hzmi,&unit_cell->axis[uc_axis[1]]->vec[1],1,shiftyz,1);

        for(IZ=0;IZ<Nz;IZ++){

            c_dcopy(N2D,shiftyz,1,at_pos,1);

            for(IA=0;IA<NA;IA++){

                c_daxpy(N2D,1.0,&unit_cell->atom[IA]->coord[1],1,at_pos,1);

                for(IS=0;IS<no_element;IS++){

                    if(!strcmp(nanowire->surf[IS]->cross_section,"yes")){
                        
                        condition = is_in_quad2D(nanowire->surf[IS],at_pos);

                        if(condition){
                            if(IS<no_ch_element){At_local=At_local+1;}
                            break;
                        }
                        
                    }
                    
                }
                c_dcopy(N2D,shiftyz,1,at_pos,1);
            }
            c_daxpy(N2D,1.0,&unit_cell->axis[uc_axis[1]]->vec[1],1,shiftyz,1);
        }
        c_daxpy(N2D,1.0,&unit_cell->axis[uc_axis[0]]->vec[1],1,shifty,1);
    }
    */

    return At_local;
}

/************************************************************************************************/

void WireGenerator::check_dimension(WireStructure* nanowire)
{

    y_width = get_length(1);
    Ly      = 1e9;
    Lz      = 1e9;
     
    if(nanowire->NDim==1){

        int no_element = nanowire->no_element;
	int IE,ID;

	NyFold = nanowire->NyFold;
	Ly     = NyFold*get_length(1);
	
        if((Ly<4.0*nanowire->a0/5.0)&&(transport_type==1)){
	    Ly     = 2.0*Ly;
	    NyFold = 2*NyFold;
	}
        
        for(IE=0;IE<no_element;IE++){
	    c_dcopy(N3D,nanowire->mat[IE]->p[0]->coord,1,\
		    nanowire->mat[IE]->p[3]->coord,1);
	    c_dcopy(N3D,nanowire->mat[IE]->p[1]->coord,1,\
		    nanowire->mat[IE]->p[2]->coord,1);
	    for(ID=0;ID<N3D+1;ID++){
	        if(ID<2){
		    nanowire->mat[IE]->p[ID]->coord[1] = -Ly-100*tollim;
		}else{
		    nanowire->mat[IE]->p[ID]->coord[1] = 2*Ly-100*tollim;
		}
		nanowire->surf[IE]->p[ID]->coord[0] = nanowire->mat[IE]->p[ID]->coord[1];   
            }
        }
    }
   
    if(nanowire->NDim<=2){

        int no_element = nanowire->no_element;
	int IE,ID;

	NzFold = nanowire->NzFold;
        Lz     = NzFold*get_length(2);
        if((Lz<4.0*nanowire->a0/5.0)&&(transport_type==1)){
	    Lz     = 2.0*Lz;
	    NzFold = 2*NzFold;
	}

        for(IE=0;IE<no_element;IE++){
            for(ID=0;ID<N3D+1;ID++){
                c_dcopy(N3D,nanowire->mat[IE]->p[ID]->coord,1,\
                        nanowire->mat[IE]->p[ID+N3D+1]->coord,1);
                nanowire->mat[IE]->p[ID]->coord[2]       = -Lz-100*tollim;
                nanowire->mat[IE]->p[ID+N3D+1]->coord[2] = 2*Lz-100*tollim;
                nanowire->surf[IE]->p[ID]->coord[0] = nanowire->mat[IE]->p[ID]->coord[1];
                if((ID == 0)||(ID == 3)){
                    nanowire->surf[IE]->p[ID]->coord[1] = -Lz-100*tollim;
                }else{
                    nanowire->surf[IE]->p[ID]->coord[1] = 2*Lz-100*tollim;
                }
            }
        }

    }
    
    if(nanowire->NDim==3){

        if(lattice_type==6){

	    int ID;
	    double radius   = nanowire->mat[0]->radius[0];
	    double width    = 2.0*PI*radius;
	    double uc_width = get_length(1);
	    double Ltot     = nanowire->mat[0]->p[1]->coord[0];

	    width           = Round(width/uc_width)*uc_width-uc_width/10.0;
	    cnt_radius      = Round(width/uc_width)*uc_width/(2.0*PI);
	    
	    init_var(nanowire->mat[0]->p[0]->coord,N3D);
	    init_var(nanowire->mat[0]->p[1]->coord,N3D);
	    nanowire->mat[0]->p[1]->coord[0] = Ltot;
	    c_dcopy(N3D,nanowire->mat[0]->p[1]->coord,1,nanowire->mat[0]->p[2]->coord,1);
	    nanowire->mat[0]->p[2]->coord[1] = width;
	    c_dcopy(N3D,nanowire->mat[0]->p[2]->coord,1,nanowire->mat[0]->p[3]->coord,1);
	    nanowire->mat[0]->p[3]->coord[0] = 0.0;

	    for(ID=0;ID<N3D+1;ID++){
	        nanowire->mat[0]->p[ID]->coord[2] = ceil(10.0/unit_cell->axis[2]->vec[2])*\
		  unit_cell->axis[2]->vec[2];
	        c_dcopy(N3D,nanowire->mat[0]->p[ID]->coord,1,\
			nanowire->mat[0]->p[N3D+1+ID]->coord,1);
		nanowire->mat[0]->p[N3D+1+ID]->coord[2] = \
		  nanowire->mat[0]->p[N3D+1+ID]->coord[2]+unit_cell->axis[2]->vec[2]/2.0;
	    }

	    strcpy(nanowire->mat[0]->type,"square");

	}
    }
}

/************************************************************************************************/

double WireGenerator::get_length(int axis_ind)
{  
     
   int ID,INFO,IPIV[3];
    double Length,vec_trans[3],inv_R[9],factor;

    for(ID=0;ID<N3D;ID++){
        c_dcopy(N3D,unit_cell->axis[ID]->vec,1,&inv_R[N3D*ID],1);
    }
    
    init_var(vec_trans,N3D);
    vec_trans[axis_ind] = 1.0;

    c_dgetrf(N3D,N3D,inv_R,N3D,IPIV,&INFO);
    c_dgetrs('N',N3D,1,inv_R,N3D,IPIV,vec_trans,N3D,&INFO);

    factor = INF;
    for(ID=0;ID<N3D;ID++){
        if((abs(vec_trans[ID])<abs(factor))&&(abs(vec_trans[ID])>tollim)){
	    factor = vec_trans[ID];
	}
    }
 
    Length = 0.0;
    for(ID=0;ID<N3D;ID++){
        Length = Length+vec_trans[ID]/factor*unit_cell->axis[ID]->vec[axis_ind];
    }

    return abs(Length); 

}

/************************************************************************************************/

void WireGenerator::calc_info(WireStructure* nanowire)
{
    int no_element = nanowire->no_element,i;

    cell_area  = 0.0;
    volume     = 0.0;
    volume_tot = 0.0;
    
    for(i=0;i<no_element;i++){
  
        if(!strcmp(nanowire->mat[i]->type,"square")){

            nanowire->mat[i]->face_area[0]=calc_quad_area(nanowire->mat[i]->p[0]->coord,\
                                                          nanowire->mat[i]->p[1]->coord,\
                                                          nanowire->mat[i]->p[2]->coord,\
                                                          nanowire->mat[i]->p[3]->coord,N3D);
            c_dcopy(N3D,nanowire->mat[i]->p[0]->coord,1,nanowire->mat[i]->p_ref[0],1);
            calc_vec(nanowire->mat[i]->vec_dir[0],nanowire->mat[i]->p[0]->coord,\
                     nanowire->mat[i]->p[1]->coord,nanowire->mat[i]->p[2]->coord);
        
            nanowire->mat[i]->face_area[1]=calc_quad_area(nanowire->mat[i]->p[4]->coord,\
                                                          nanowire->mat[i]->p[5]->coord,\
                                                          nanowire->mat[i]->p[6]->coord,\
                                                          nanowire->mat[i]->p[7]->coord,N3D);
            c_dcopy(N3D,nanowire->mat[i]->p[4]->coord,1,nanowire->mat[i]->p_ref[1],1);
            calc_vec(nanowire->mat[i]->vec_dir[1],nanowire->mat[i]->p[4]->coord,\
                     nanowire->mat[i]->p[5]->coord,nanowire->mat[i]->p[6]->coord);

            nanowire->mat[i]->face_area[2]=calc_quad_area(nanowire->mat[i]->p[0]->coord,\
                                                          nanowire->mat[i]->p[4]->coord,\
                                                          nanowire->mat[i]->p[7]->coord,\
                                                          nanowire->mat[i]->p[3]->coord,N3D);
            c_dcopy(N3D,nanowire->mat[i]->p[0]->coord,1,nanowire->mat[i]->p_ref[2],1);
            calc_vec(nanowire->mat[i]->vec_dir[2],nanowire->mat[i]->p[0]->coord,\
                     nanowire->mat[i]->p[4]->coord,nanowire->mat[i]->p[7]->coord);
        
            nanowire->mat[i]->face_area[3]=calc_quad_area(nanowire->mat[i]->p[0]->coord,\
                                                          nanowire->mat[i]->p[1]->coord,\
                                                          nanowire->mat[i]->p[5]->coord,\
                                                          nanowire->mat[i]->p[4]->coord,N3D);
            c_dcopy(N3D,nanowire->mat[i]->p[0]->coord,1,nanowire->mat[i]->p_ref[3],1);
            calc_vec(nanowire->mat[i]->vec_dir[3],nanowire->mat[i]->p[0]->coord,\
                     nanowire->mat[i]->p[1]->coord,nanowire->mat[i]->p[5]->coord);
        
            nanowire->mat[i]->face_area[4]=calc_quad_area(nanowire->mat[i]->p[1]->coord,\
                                                          nanowire->mat[i]->p[2]->coord,\
                                                          nanowire->mat[i]->p[6]->coord,\
                                                          nanowire->mat[i]->p[5]->coord,N3D);
            c_dcopy(N3D,nanowire->mat[i]->p[1]->coord,1,nanowire->mat[i]->p_ref[4],1);
            calc_vec(nanowire->mat[i]->vec_dir[4],nanowire->mat[i]->p[1]->coord,\
                     nanowire->mat[i]->p[2]->coord,nanowire->mat[i]->p[6]->coord);
        
            nanowire->mat[i]->face_area[5]=calc_quad_area(nanowire->mat[i]->p[2]->coord,\
                                                          nanowire->mat[i]->p[3]->coord,\
                                                          nanowire->mat[i]->p[7]->coord,\
                                                          nanowire->mat[i]->p[6]->coord,N3D);
            c_dcopy(N3D,nanowire->mat[i]->p[2]->coord,1,nanowire->mat[i]->p_ref[5],1);
            calc_vec(nanowire->mat[i]->vec_dir[5],nanowire->mat[i]->p[2]->coord,\
                     nanowire->mat[i]->p[3]->coord,nanowire->mat[i]->p[7]->coord);

            nanowire->mat[i]->volume = calc_quad_volume(nanowire->mat[i]);
	    
	    if(!strcmp(nanowire->surf[i]->type,"square")){
	        nanowire->surf[i]->face_area = calc_quad_area(nanowire->surf[i]->p[0]->coord,\
							      nanowire->surf[i]->p[1]->coord,\
							      nanowire->surf[i]->p[2]->coord,\
							      nanowire->surf[i]->p[3]->coord,N2D);
	    }else{
	        nanowire->surf[i]->face_area = PI*nanowire->surf[i]->radius[0]*\
		  nanowire->surf[i]->radius[1];
	    }
	    
            if((i<nanowire->no_ch_element)&&(!strcmp(nanowire->surf[i]->cross_section,"yes"))){
                cell_area = cell_area+nanowire->surf[i]->face_area;
            }
	    
	    volume_tot = volume_tot+nanowire->mat[i]->volume;
 
            if((i<nanowire->no_ch_element)||(i>=(nanowire->no_ch_element+nanowire->no_ox_element))){
		volume = volume+nanowire->mat[i]->volume;
	    }
            
        }

	if(!strcmp(nanowire->mat[i]->type,"circle")){
	
	    nanowire->mat[i]->volume = PI*nanowire->surf[i]->radius[0]*nanowire->surf[i]->radius[1]*
	      (nanowire->mat[i]->p[1]->coord[0]-nanowire->mat[i]->p[0]->coord[0]);

            if((i<nanowire->no_ch_element)&&(!strcmp(nanowire->surf[i]->cross_section,"yes"))){
                cell_area = cell_area+PI*nanowire->surf[i]->radius[0]*\
		  nanowire->surf[i]->radius[1];
            }

            if((i<nanowire->no_ch_element)||(i>=(nanowire->no_ch_element+nanowire->no_ox_element))){
		volume = volume+nanowire->mat[i]->volume;
	    }

	    volume_tot = volume_tot+PI*nanowire->surf[i]->radius[0]*nanowire->surf[i]->radius[1]*
	      (nanowire->mat[i]->p[1]->coord[0]-nanowire->mat[i]->p[0]->coord[0]);
            
        }

	if(!strcmp(nanowire->mat[i]->type,"sphere")){
	
	    nanowire->mat[i]->volume = 4.0/3.0*PI*pow(nanowire->surf[i]->radius[0],3.0);

            if((i<nanowire->no_ch_element)&&(!strcmp(nanowire->surf[i]->cross_section,"yes"))){
                cell_area = cell_area+PI*nanowire->surf[i]->radius[0]*\
		  nanowire->surf[i]->radius[1];
            }

            if((i<nanowire->no_ch_element)||(i>=(nanowire->no_ch_element+nanowire->no_ox_element))){
		volume = volume+nanowire->mat[i]->volume;
	    }

	    volume_tot = volume_tot+4.0/3.0*PI*pow(nanowire->surf[i]->radius[0],3.0);
        }
 
    }

}

/************************************************************************************************/

void WireGenerator::find_3D_boundary(int* hxmin, int* hxmax, int* hymin, int* hymax, \
                                     int* hzmin, int* hzmax, WireStructure* nanowire)
{
    int no_element=nanowire->no_element,i,j,IPIV[3],INFO;
    double xmin=INF,xmax=-INF,ymin=INF,ymax=-INF,zmin=INF,zmax=-INF;

    for(i=0;i<no_element;i++){
        
        if(!strcmp(nanowire->mat[i]->type,"square")){
            
            for(j=0;j<8;j++){
                if(nanowire->mat[i]->p[j]->coord[0]<xmin){
                    xmin = nanowire->mat[i]->p[j]->coord[0];
                }
                if(nanowire->mat[i]->p[j]->coord[0]>xmax){
                    xmax = nanowire->mat[i]->p[j]->coord[0];
                }
                if(nanowire->mat[i]->p[j]->coord[1]<ymin){
                    ymin = nanowire->mat[i]->p[j]->coord[1];
                }
                if(nanowire->mat[i]->p[j]->coord[1]>ymax){
                    ymax = nanowire->mat[i]->p[j]->coord[1];
                }
                if(nanowire->mat[i]->p[j]->coord[2]<zmin){
                    zmin = nanowire->mat[i]->p[j]->coord[2];
                }
                if(nanowire->mat[i]->p[j]->coord[2]>zmax){
                    zmax = nanowire->mat[i]->p[j]->coord[2];
                }
            }
        }
       
        if(!strcmp(nanowire->mat[i]->type,"circle")){

            if(nanowire->mat[i]->p[0]->coord[0]<xmin){
                xmin = nanowire->surf[i]->p[0]->coord[0];
            }
            if(nanowire->mat[i]->p[1]->coord[0]>xmax){
                xmax = nanowire->mat[i]->p[1]->coord[0];
            }
            if((nanowire->mat[i]->p[0]->coord[1]-nanowire->mat[i]->radius[0])<ymin){
                ymin = nanowire->mat[i]->p[0]->coord[1]-nanowire->mat[i]->radius[0];
            }
            if((nanowire->mat[i]->p[0]->coord[1]+nanowire->mat[i]->radius[0])>ymax){
                ymax = nanowire->mat[i]->p[0]->coord[1]+nanowire->mat[i]->radius[0];
            }
            if((nanowire->mat[i]->p[0]->coord[2]-nanowire->mat[i]->radius[1])<zmin){
                zmin = nanowire->mat[i]->p[0]->coord[2]-nanowire->mat[i]->radius[1];
            }
            if((nanowire->mat[i]->p[0]->coord[2]+nanowire->mat[i]->radius[1])>zmax){
                zmax = nanowire->mat[i]->p[0]->coord[2]+nanowire->mat[i]->radius[1];
            }   
        }

        if(!strcmp(nanowire->mat[i]->type,"sphere")){

            if(nanowire->mat[i]->p[0]->coord[0]-nanowire->mat[i]->radius[0]<xmin){
                xmin = nanowire->mat[i]->p[0]->coord[0]-nanowire->mat[i]->radius[0];
            }
            if(nanowire->mat[i]->p[1]->coord[0]+nanowire->mat[i]->radius[0]>xmax){
                xmax = nanowire->mat[i]->p[0]->coord[0]+nanowire->mat[i]->radius[0];
            }
            if((nanowire->mat[i]->p[0]->coord[1]-nanowire->mat[i]->radius[0])<ymin){
                ymin = nanowire->mat[i]->p[0]->coord[1]-nanowire->mat[i]->radius[0];
            }
            if((nanowire->mat[i]->p[0]->coord[1]+nanowire->mat[i]->radius[0])>ymax){
                ymax = nanowire->mat[i]->p[0]->coord[1]+nanowire->mat[i]->radius[0];
            }
            if((nanowire->mat[i]->p[0]->coord[2]-nanowire->mat[i]->radius[0])<zmin){
                zmin = nanowire->mat[i]->p[0]->coord[2]-nanowire->mat[i]->radius[0];
            }
            if((nanowire->mat[i]->p[0]->coord[2]+nanowire->mat[i]->radius[0])>zmax){
                zmax = nanowire->mat[i]->p[0]->coord[2]+nanowire->mat[i]->radius[0];
            }
            
        }
            
    }

    lenx = xmax - xmin;
    leny = ymax - ymin;
    lenz = zmax - zmin;

    double P[24]={xmin,ymin,zmin,xmax,ymin,zmin,xmax,ymax,zmin,\
                  xmin,ymax,zmin,xmin,ymin,zmax,xmax,ymin,zmax,\
                  xmax,ymax,zmax,xmin,ymax,zmax};
    double V[9];

    for(i=0;i<N3D;i++){
        for(j=0;j<N3D;j++){
            V[i+j*N3D] = unit_cell->axis[j]->vec[i];
        }
    }
    c_dgetrf(N3D,N3D,V,N3D,IPIV,&INFO);
    c_dgetrs('N',N3D,8,V,N3D,IPIV,P,N3D,&INFO);

    for(i=0;i<8;i++){
        if(Round(P[i*N3D]-1)<*hxmin){
            *hxmin = Round(P[i*N3D]-1);
        }
        if(Round(P[i*N3D]+1)>*hxmax){
            *hxmax = Round(P[i*N3D]+1);
        }
        if(Round(P[1+i*N3D]-1)<*hymin){
            *hymin = Round(P[1+i*N3D]-1);
        }
        if(Round(P[1+i*N3D]+1)>*hymax){
            *hymax = Round(P[1+i*N3D]+1);
        }
        if(Round(P[2+i*N3D]-1)<*hzmin){
            *hzmin = Round(P[2+i*N3D]-1);
        }
        if(Round(P[2+i*N3D]+1)>*hzmax){
            *hzmax = Round(P[2+i*N3D]+1);
        }
    }
    
}

/************************************************************************************************/

void WireGenerator::find_2D_boundary(int* hymin, int* hymax, int* hzmin, int* hzmax, \
                                     WireStructure* nanowire)
{
    int no_element=nanowire->no_element,i,j,IPIV[2],INFO,uc_axis[2];
    double ymin=INF,ymax=-INF,zmin=INF,zmax=-INF;

    for(i=0;i<no_element;i++){
        
        if(!strcmp(nanowire->surf[i]->type,"square")){

            for(j=0;j<4;j++){
                if(nanowire->surf[i]->p[j]->coord[0]<ymin){
                    ymin = nanowire->surf[i]->p[j]->coord[0];
                }
                if(nanowire->surf[i]->p[j]->coord[0]>ymax){
                    ymax = nanowire->surf[i]->p[j]->coord[0];
                }
                if(nanowire->surf[i]->p[j]->coord[1]<zmin){
                    zmin = nanowire->surf[i]->p[j]->coord[1];
                }
                if(nanowire->surf[i]->p[j]->coord[1]>zmax){
                    zmax = nanowire->surf[i]->p[j]->coord[1];
                }
            }
            
        }else{
            
            if((nanowire->surf[i]->p[0]->coord[0]-nanowire->surf[i]->radius[0])<ymin){
                ymin = nanowire->surf[i]->p[0]->coord[0]-nanowire->surf[i]->radius[0];
            }
            if((nanowire->surf[i]->p[0]->coord[0]+nanowire->surf[i]->radius[0])>ymax){
                ymax = nanowire->surf[i]->p[0]->coord[0]+nanowire->surf[i]->radius[0];
            }
            if((nanowire->surf[i]->p[0]->coord[1]-nanowire->surf[i]->radius[0])<zmin){
                zmin = nanowire->surf[i]->p[0]->coord[1]-nanowire->surf[i]->radius[0];
            }
            if((nanowire->surf[i]->p[0]->coord[1]+nanowire->surf[i]->radius[0])>zmax){
                zmax = nanowire->surf[i]->p[0]->coord[1]+nanowire->surf[i]->radius[0];
            }
            
        }
    }

    double P[8]={ymin,zmin,ymax,zmin,ymax,zmax,ymin,zmax};
    double V[4];

    check_uc_axis(uc_axis);
    
    for(i=0;i<N2D;i++){
        for(j=0;j<N2D;j++){
            V[i+j*N2D]=unit_cell->axis[uc_axis[j]]->vec[i+1];
        }
    }
    c_dgetrf(N2D,N2D,V,N2D,IPIV,&INFO);
    c_dgetrs('N',N2D,4,V,N2D,IPIV,P,N2D,&INFO);

    for(i=0;i<4;i++){
        if(Round(P[i*N2D]-1)<*hymin){
            *hymin = Round(P[i*N2D]-1);
        }
        if(Round(P[i*N2D]+1)>*hymax){
            *hymax = Round(P[i*N2D]+1);
        }
        if(Round(P[1+i*N2D]-1)<*hzmin){
            *hzmin = Round(P[1+i*N2D]-1);
        }
        if(Round(P[1+i*N2D]+1)>*hzmax){
            *hzmax = Round(P[1+i*N2D]+1);
        }
    }
}

/************************************************************************************************/

void WireGenerator::check_uc_axis(int *uc_axis)
{
    int a1,a2;
    double d;

    for(a1=0;a1<2;a1++){
        for(a2=a1+1;a2<3;a2++){
            d=unit_cell->axis[a1]->vec[1]*unit_cell->axis[a2]->vec[2]-\
                unit_cell->axis[a1]->vec[2]*unit_cell->axis[a2]->vec[1];
            if(abs(d)>tollim){
                uc_axis[0] = a1;
                uc_axis[1] = a2;
            }
        }
    }
}

/************************************************************************************************/

void WireGenerator::make_unit_cell(const char* first_atom,double a0,double c0,double u0)
{

    int IA,IB,IP,ID,JD;
    int neigh[6][6];
    double v[6][3];
    double a[3][3],b[6][3],s[6];
    double delta_c;
    double a_factor[3]={a0,a0,c0};
    double bshift[6][6];
    double bond_pass;
    double BZPoint[10][3];
   
    switch(lattice_type){
    case 1:

	v[0][0] = 0.0;  v[0][1] = 0.0;  v[0][2] = 0.0;
	v[1][0] = 0.25; v[1][1] = 0.25; v[1][2] = 0.25;

	a[0][0] = 0.5; a[0][1] = 0.5; a[0][2] = 0.0;
	a[1][0] = 0.5; a[1][1] = 0.0; a[1][2] = 0.5;
	a[2][0] = 0.0; a[2][1] = 0.5; a[2][2] = 0.5;

	b[0][0] = 0.25;  b[0][1] = 0.25;  b[0][2] = 0.25;
	b[1][0] = 0.25;  b[1][1] = -0.25; b[1][2] = -0.25;
	b[2][0] = -0.25; b[2][1] = 0.25;  b[2][2] = -0.25;
	b[3][0] = -0.25; b[3][1] = -0.25; b[3][2] = 0.25;

	s[0] = 1.0; s[1] = -1.0;

	neigh[0][0] = 2; neigh[0][1] = 2; neigh[0][2] = 2; neigh[0][3] = 2;
	neigh[1][0] = 1; neigh[1][1] = 1; neigh[1][2] = 1; neigh[1][3] = 1;

	bshift[0][0] = 0.0; bshift[0][1] = 0.0; bshift[0][2] = 0.0; bshift[0][3] = 0.0;
	bshift[1][0] = 0.0; bshift[1][1] = 0.0; bshift[1][2] = 0.0; bshift[1][3] = 0.0;

        BZPoint[0][0] = PI;     BZPoint[0][1] = PI;     BZPoint[0][2] = PI;
        BZPoint[1][0] = 0.0;    BZPoint[1][1] = 0.0;    BZPoint[1][2] = 0.0;
        BZPoint[2][0] = 2.0*PI; BZPoint[2][1] = 0.0;    BZPoint[2][2] = 0.0;
        BZPoint[3][0] = 2.0*PI; BZPoint[3][1] = PI/2.0; BZPoint[3][2] = PI/2.0;
        BZPoint[4][0] = 0.0;    BZPoint[4][1] = 0.0;    BZPoint[4][2] = 0.0;
        
	bond_pass = 1.0;

	Vol_atom = a0*a0*a0/8.0*1e-27;

	break;

    case 3:

        delta_c = u0;

	v[0][0] = 0.0; v[0][1] = 1.0/sqrt(3.0); v[0][2] = (1.0/2.0+delta_c);
	v[1][0] = 0.0; v[1][1] = 0.0;           v[1][2] = 0.0;
	v[2][0] = 0.0; v[2][1] = 0.0;           v[2][2] = delta_c;
	v[3][0] = 0.0; v[3][1] = 1.0/sqrt(3.0); v[3][2] = 1.0/2.0;

	a[0][0] = 1.0;      a[0][1] = 0.0;           a[0][2] = 0.0;
	a[1][0] = -1.0/2.0; a[1][1] = sqrt(3.0)/2.0; a[1][2] = 0.0;
	a[2][0] = 0.0;      a[2][1] = 0.0;           a[2][2] = 1.0;

	b[0][0] = 0.0;      b[0][1] = -1.0/sqrt(3.0); b[0][2] = 1.0/2.0-delta_c;
	b[1][0] = -1.0/2.0; b[1][1] = sqrt(3.0)/6.0;  b[1][2] = 1.0/2.0-delta_c;
	b[2][0] = 1.0/2.0;  b[2][1] = sqrt(3.0)/6.0;  b[2][2] = 1.0/2.0-delta_c;
	b[3][0] = 0.0;      b[3][1] = 0.0;            b[3][2] = -delta_c;

	s[0] = 1.0; s[1] = -1.0; s[2] = -1.0; s[3] = 1.0;

	neigh[0][0] = 2; neigh[0][1] = 2; neigh[0][2] = 2; neigh[0][3] = 4;
	neigh[1][0] = 1; neigh[1][1] = 1; neigh[1][2] = 1; neigh[1][3] = 3;
	neigh[2][0] = 4; neigh[2][1] = 4; neigh[2][2] = 4; neigh[2][3] = 2;
	neigh[3][0] = 3; neigh[3][1] = 3; neigh[3][2] = 3; neigh[3][3] = 1;

	bshift[0][0] = 0.0;              bshift[0][1] = 0.0;              bshift[0][2] = 0.0; 
	    bshift[0][3] = 0.0;
	bshift[1][0] = 0.0;              bshift[1][1] = 0.0;              bshift[1][2] = 0.0; 
	    bshift[1][3] = 0.0;
	bshift[2][0] = -1.0+2.0*delta_c; bshift[2][1] = -1.0+2.0*delta_c; bshift[2][2] = -1.0+2.0*delta_c;
	    bshift[2][3] = 2.0*delta_c;
	bshift[3][0] = -1.0+2.0*delta_c; bshift[3][1] = -1.0+2.0*delta_c; bshift[3][2] = -1.0+2.0*delta_c;
	    bshift[3][3] = 2.0*delta_c;

        BZPoint[0][0] = 0.0;         BZPoint[0][1] = 0.0;              BZPoint[0][2] = 0.0;
        BZPoint[1][0] = -2.0*PI/3.0; BZPoint[1][1] = 2.0*PI/sqrt(3.0); BZPoint[1][2] = 0.0;
        BZPoint[2][0] = 0.0;         BZPoint[2][1] = 2.0*PI/sqrt(3.0); BZPoint[2][2] = 0.0;
        BZPoint[3][0] = 0.0;         BZPoint[3][1] = 0.0;              BZPoint[3][2] = 0.0;
        BZPoint[4][0] = 0.0;         BZPoint[4][1] = 0.0;              BZPoint[4][2] = PI;
        BZPoint[5][0] = -2.0*PI/3.0; BZPoint[5][1] = 2.0*PI/sqrt(3.0); BZPoint[5][2] = PI;
        BZPoint[6][0] = 0.0;         BZPoint[6][1] = 2.0*PI/sqrt(3.0); BZPoint[6][2] = PI;
        BZPoint[7][0] = 0.0;         BZPoint[7][1] = 0.0;              BZPoint[7][2] = PI;

	bond_pass = 1.0;

	Vol_atom = sqrt(3.0)*a0*a0*c0/8.0*1e-27;
 
        break;

    case 4:

	v[0][0] = 0.0; v[0][1] = 0.0; v[0][2] = 0.0;
	v[1][0] = 0.0; v[1][1] = 1.0; v[1][2] = 0.0;

	a[0][0] = sqrt(3.0)/2.0;  a[0][1] = 3.0/2.0; a[0][2] = 0.0;
	a[1][0] = -sqrt(3.0)/2.0; a[1][1] = 3.0/2.0; a[1][2] = 0.0;
	a[2][0] = 0.0;            a[2][1] = 3.0/2.0; a[2][2] = 3.0/2.0;

	b[0][0] = 0.0;            b[0][1] = 1.0;      b[0][2] = 0.0;
	b[1][0] = sqrt(3.0)/2.0;  b[1][1] = -1.0/2.0; b[1][2] = 0.0;
	b[2][0] = -sqrt(3.0)/2.0; b[2][1] = -1.0/2.0; b[2][2] = 0.0;

	s[0] = 1.0; s[1] = -1.0;

	neigh[0][0] = 2; neigh[0][1] = 2; neigh[0][2] = 2; neigh[0][3] = 2;
	neigh[1][0] = 1; neigh[1][1] = 1; neigh[1][2] = 1; neigh[1][3] = 1;

	bshift[0][0] = 0.0; bshift[0][1] = 0.0; bshift[0][2] = 0.0;
	bshift[1][0] = 0.0; bshift[1][1] = 0.0; bshift[1][2] = 0.0;

        BZPoint[0][0] = PI/sqrt(3.0);           BZPoint[0][1] = PI/3.0; BZPoint[0][2] = 0.0;
        BZPoint[1][0] = 0.0;                    BZPoint[1][1] = 0.0;    BZPoint[1][2] = 0.0;
        BZPoint[2][0] = 4.0*PI/(3.0*sqrt(3.0)); BZPoint[2][1] = 0.0;    BZPoint[2][2] = 0.0;
        BZPoint[3][0] = PI/sqrt(3.0);           BZPoint[3][1] = PI/3.0; BZPoint[3][2] = 0.0;

	bond_pass = 0.0;

	Vol_atom = 9.0*sqrt(3.0)/8.0*a0*a0*a0*1e-27;

	break;

    case 5:

        delta_c = 0.0094;

        v[0][0] = 0.0; v[0][1] = 0.0; v[0][2] = 0.0;
	v[1][0] = 0.0; v[1][1] = 0.0; v[1][2] = 0.2+delta_c;
	v[2][0] = 0.0; v[2][1] = 0.0; v[2][2] = 0.4;
	v[3][0] = 0.0; v[3][1] = 0.0; v[3][2] = 0.6;
	v[4][0] = 0.0; v[4][1] = 0.0; v[4][2] = 0.8-delta_c;

        a[0][0] = 1.0/sqrt(3.0);        a[0][1] = 0.0;  a[0][2] = 1.0/3.0;
	a[1][0] = -1.0/(2.0*sqrt(3.0)); a[1][1] = 0.5;  a[1][2] = 1.0/3.0;
	a[2][0] = -1.0/(2.0*sqrt(3.0)); a[2][1] = -0.5; a[2][2] = 1.0/3.0;

	b[0][0] = 1.0/sqrt(3.0);        b[0][1] = 0.0;  b[0][2] = -1.0/15.0;
	b[1][0] = -1.0/(2.0*sqrt(3.0)); b[1][1] = 0.5;  b[1][2] = -1.0/15.0;
	b[2][0] = -1.0/(2.0*sqrt(3.0)); b[2][1] = -0.5; b[2][2] = -1.0/15.0;
	b[3][0] = -1.0/sqrt(3.0);       b[3][1] = 0.0;  b[3][2] = 1.0/15.0;
	b[4][0] = 1.0/(2.0*sqrt(3.0));  b[4][1] = -0.5; b[4][2] = 1.0/15.0;
	b[5][0] = 1.0/(2.0*sqrt(3.0));  b[5][1] = 0.5;  b[5][2] = 1.0/15.0;

	s[0] = 1.0; s[1] = 1.0; s[2] = 1.0; s[3] = 1.0; s[4] = 1.0;

	neigh[0][0] = 4; neigh[0][1] = 4; neigh[0][2] = 4; neigh[0][3] = 3; neigh[0][4] = 3; neigh[0][5] = 3;
	neigh[1][0] = 5; neigh[1][1] = 5; neigh[1][2] = 5; neigh[1][3] = 4; neigh[1][4] = 4; neigh[1][5] = 4;
	neigh[2][0] = 1; neigh[2][1] = 1; neigh[2][2] = 1; neigh[2][3] = 5; neigh[2][4] = 5; neigh[2][5] = 5;
	neigh[3][0] = 2; neigh[3][1] = 2; neigh[3][2] = 2; neigh[3][3] = 1; neigh[3][4] = 1; neigh[3][5] = 1;
	neigh[4][0] = 3; neigh[4][1] = 3; neigh[4][2] = 3; neigh[4][3] = 2; neigh[4][4] = 2; neigh[4][5] = 2;

	bshift[0][0] = 0.0;        bshift[0][1] = 0.0;        bshift[0][2] = 0.0;        bshift[0][3] = 0.0;
	      bshift[0][4] = 0.0;       bshift[0][5] = 0.0;
	bshift[1][0] = -2*delta_c; bshift[1][1] = -2*delta_c; bshift[1][2] = -2*delta_c; bshift[1][3] = -delta_c;
	      bshift[1][4] = -delta_c;  bshift[1][5] = -delta_c;
	bshift[2][0] = 0.0;        bshift[2][1] = 0.0;        bshift[2][2] = 0.0;        bshift[2][3] = -delta_c;  
              bshift[2][4] = -delta_c;  bshift[2][5] = -delta_c;
	bshift[3][0] = delta_c;    bshift[3][1] = delta_c;    bshift[3][2] = delta_c;    bshift[3][3] = 0.0;       
              bshift[3][4] = 0.0;       bshift[3][5] = 0.0;
	bshift[4][0] = delta_c;    bshift[4][1] = delta_c;    bshift[4][2] = delta_c;    bshift[4][3] = 2*delta_c; 
              bshift[4][4] = 2*delta_c; bshift[4][5] = 2*delta_c;

        BZPoint[0][0] = 1.1201*PI;        BZPoint[0][1] = 0.0; BZPoint[0][2] = 3.0*PI;
        BZPoint[1][0] = 0.0;              BZPoint[1][1] = 0.0; BZPoint[1][2] = 0.0;
        BZPoint[2][0] = 0.0;              BZPoint[2][1] = 0.0; BZPoint[2][2] = 3.0*PI;
        BZPoint[3][0] = 2.0*PI/sqrt(3.0); BZPoint[3][1] = 0.0; BZPoint[3][2] = -2.0*PI;
        BZPoint[4][0] = 0.0;              BZPoint[4][1] = 0.0; BZPoint[4][2] = 0.0;
        BZPoint[5][0] = 2.0*PI/sqrt(3.0); BZPoint[5][1] = 0.0; BZPoint[5][2] = PI;      

	bond_pass = 1.0;

	Vol_atom = 1.0/(10.0*sqrt(3.0))*a0*a0*c0*1e-27;

        break;

    case 6:

	v[0][0] = 0.0; v[0][1] = 0.0; v[0][2] = 0.0;
	v[1][0] = 0.0; v[1][1] = 1.0; v[1][2] = 0.0;

	a[0][0] = sqrt(3.0)/2.0;  a[0][1] = 3.0/2.0; a[0][2] = 0.0;
	a[1][0] = -sqrt(3.0)/2.0; a[1][1] = 3.0/2.0; a[1][2] = 0.0;
	a[2][0] = 0.0;            a[2][1] = 3.0/2.0; a[2][2] = 1.0;

	b[0][0] = 0.0;            b[0][1] = 1.0;  b[0][2] = 0.0;
	b[1][0] = sqrt(3.0)/2.0;  b[1][1] = -0.5; b[1][2] = 0.0;
	b[2][0] = -sqrt(3.0)/2.0; b[2][1] = -0.5; b[2][2] = 0.0;

	s[0] = 1.0; s[1] = -1.0;

	neigh[0][0] = 2; neigh[0][1] = 2; neigh[0][2] = 2; neigh[0][3] = 2;
	neigh[1][0] = 1; neigh[1][1] = 1; neigh[1][2] = 1; neigh[1][3] = 1;

	bshift[0][0] = 0.0; bshift[0][1] = 0.0; bshift[0][2] = 0.0;
	bshift[1][0] = 0.0; bshift[1][1] = 0.0; bshift[1][2] = 0.0;

        BZPoint[0][0] = PI/sqrt(3.0);           BZPoint[0][1] = PI/3.0; BZPoint[0][2] = 0.0;
        BZPoint[1][0] = 0.0;                    BZPoint[1][1] = 0.0;    BZPoint[1][2] = 0.0;
        BZPoint[2][0] = 4.0*PI/(3.0*sqrt(3.0)); BZPoint[2][1] = 0.0;    BZPoint[2][2] = 0.0;
        BZPoint[3][0] = PI/sqrt(3.0);           BZPoint[3][1] = PI/3.0; BZPoint[3][2] = 0.0;
        
	bond_pass = 0.0;

	Vol_atom = 3.0*sqrt(3.0)/4.0*a0*a0*a0*1e-27;
	
	break;

    case 7:

	v[0][0] = 0.0; v[0][1] = 0.0;  v[0][2] = 0.0;
	v[1][0] = 0.0; v[1][1] = 1.0;  v[1][2] = 0.0;
	v[2][0] = 0.0; v[2][1] = -1.0; v[2][2] = 1.0;
	v[3][0] = 0.0; v[3][1] = 0.0;  v[3][2] = 1.0;

	a[0][0] = sqrt(3.0)/2.0;  a[0][1] = 3.0/2.0; a[0][2] = 0.0;
	a[1][0] = -sqrt(3.0)/2.0; a[1][1] = 3.0/2.0; a[1][2] = 0.0;
	a[2][0] = 0.0;            a[2][1] = 3.0/2.0; a[2][2] = 3.0/2.0;

	b[0][0] = 0.0;            b[0][1] = 1.0;      b[0][2] = 0.0;
	b[1][0] = sqrt(3.0)/2.0;  b[1][1] = -1.0/2.0; b[1][2] = 0.0;
	b[2][0] = -sqrt(3.0)/2.0; b[2][1] = -1.0/2.0; b[2][2] = 0.0;
	b[3][0] = 0.0;            b[3][1] = 0.0;      b[3][2] = 1.0;

	s[0] = 1.0; s[1] = -1.0; s[2] = 1.0; s[3] = -1.0;

	neigh[0][0] = 2; neigh[0][1] = 2; neigh[0][2] = 2; neigh[0][3] = 4;
	neigh[1][0] = 1; neigh[1][1] = 1; neigh[1][2] = 1; neigh[1][3] = 0;
	neigh[2][0] = 4; neigh[2][1] = 4; neigh[2][2] = 4; neigh[2][3] = 0;
	neigh[3][0] = 3; neigh[3][1] = 3; neigh[3][2] = 3; neigh[3][3] = 1;

	bshift[0][0] = 0.0; bshift[0][1] = 0.0; bshift[0][2] = 0.0; bshift[0][3] = 0.0;
	bshift[1][0] = 0.0; bshift[1][1] = 0.0; bshift[1][2] = 0.0; bshift[1][3] = 0.0;
	bshift[2][0] = 0.0; bshift[2][1] = 0.0; bshift[2][2] = 0.0; bshift[2][3] = 0.0;
	bshift[3][0] = 0.0; bshift[3][1] = 0.0; bshift[3][2] = 0.0; bshift[3][3] = 0.0;

        BZPoint[0][0] = PI/sqrt(3.0);           BZPoint[0][1] = PI/3.0; BZPoint[0][2] = 0.0;
        BZPoint[1][0] = 0.0;                    BZPoint[1][1] = 0.0;    BZPoint[1][2] = 0.0;
        BZPoint[2][0] = 4.0*PI/(3.0*sqrt(3.0)); BZPoint[2][1] = 0.0;    BZPoint[2][2] = 0.0;
        BZPoint[3][0] = PI/sqrt(3.0);           BZPoint[3][1] = PI/3.0; BZPoint[3][2] = 0.0;
        
	bond_pass = 0.0;

	Vol_atom = 9.0*sqrt(3.0)/8.0*a0*c0*a0*1e-27;

	break;

    case 8:

	v[0][0] = 0.0; v[0][1] = 0.0;  v[0][2] = 0.0;
	v[1][0] = 0.0; v[1][1] = 1.0;  v[1][2] = 0.0;
	v[2][0] = 0.0; v[2][1] = -1.0; v[2][2] = 1.0;
	v[3][0] = 0.0; v[3][1] = 0.0;  v[3][2] = 1.0;

	a[0][0] = sqrt(3.0)/2.0;  a[0][1] = 3.0/2.0; a[0][2] = 0.0;
	a[1][0] = -sqrt(3.0)/2.0; a[1][1] = 3.0/2.0; a[1][2] = 0.0;
	a[2][0] = 0.0;            a[2][1] = 0.0;     a[2][2] = 2.0;

	b[0][0] = 0.0;            b[0][1] = 1.0;      b[0][2] = 0.0;
	b[1][0] = sqrt(3.0)/2.0;  b[1][1] = -1.0/2.0; b[1][2] = 0.0;
	b[2][0] = -sqrt(3.0)/2.0; b[2][1] = -1.0/2.0; b[2][2] = 0.0;
	b[3][0] = 0.0;            b[3][1] = 0.0;      b[3][2] = 1.0;
        b[4][0] = 0.0;            b[4][1] = 0.0;      b[4][2] = -1.0;

	s[0] = 1.0; s[1] = -1.0; s[2] = 1.0; s[3] = -1.0;

	neigh[0][0] = 2; neigh[0][1] = 2; neigh[0][2] = 2; neigh[0][3] = 4; neigh[0][4] = 4;
	neigh[1][0] = 1; neigh[1][1] = 1; neigh[1][2] = 1; neigh[1][3] = 0; neigh[1][4] = 0;
	neigh[2][0] = 4; neigh[2][1] = 4; neigh[2][2] = 4; neigh[2][3] = 0; neigh[2][4] = 0;
	neigh[3][0] = 3; neigh[3][1] = 3; neigh[3][2] = 3; neigh[3][3] = 1; neigh[3][4] = 1;

	bshift[0][0] = 0.0; bshift[0][1] = 0.0; bshift[0][2] = 0.0; bshift[0][3] = 0.0; bshift[0][4] = 0.0;
	bshift[1][0] = 0.0; bshift[1][1] = 0.0; bshift[1][2] = 0.0; bshift[1][3] = 0.0; bshift[1][4] = 0.0;
	bshift[2][0] = 0.0; bshift[2][1] = 0.0; bshift[2][2] = 0.0; bshift[2][3] = 0.0; bshift[2][4] = 0.0;
	bshift[3][0] = 0.0; bshift[3][1] = 0.0; bshift[3][2] = 0.0; bshift[3][3] = 0.0; bshift[3][4] = 0.0;

        BZPoint[0][0] = PI/sqrt(3.0);           BZPoint[0][1] = PI/3.0; BZPoint[0][2] = 0.0;
        BZPoint[1][0] = 0.0;                    BZPoint[1][1] = 0.0;    BZPoint[1][2] = 0.0;
        BZPoint[2][0] = 4.0*PI/(3.0*sqrt(3.0)); BZPoint[2][1] = 0.0;    BZPoint[2][2] = 0.0;
        BZPoint[3][0] = PI/sqrt(3.0);           BZPoint[3][1] = PI/3.0; BZPoint[3][2] = 0.0;
        
	bond_pass = 0.0;

	Vol_atom = 9.0*sqrt(3.0)/8.0*a0*a0*c0*1e-27;

	break;

    case 9:

	v[0][0] = 0.0;  v[0][1] = 0.0;  v[0][2] = 0.0;
	v[1][0] = 0.5;  v[1][1] = 0.0;  v[1][2] = 0.0;

	a[0][0] = 0.5; a[0][1] = 0.5; a[0][2] = 0.0;
	a[1][0] = 0.5; a[1][1] = 0.0; a[1][2] = 0.5;
	a[2][0] = 0.0; a[2][1] = 0.5; a[2][2] = 0.5;

	b[0][0] = 0.5;  b[0][1] = 0.0;  b[0][2] = 0.0;
	b[1][0] = 0.0;  b[1][1] = 0.5;  b[1][2] = 0.0;
	b[2][0] = 0.0;  b[2][1] = 0.0;  b[2][2] = 0.5;
	b[3][0] = -0.5; b[3][1] = 0.0;  b[3][2] = 0.0;
	b[4][0] = 0.0;  b[4][1] = -0.5; b[4][2] = 0.0;
	b[5][0] = 0.0;  b[5][1] = 0.0;  b[5][2] = -0.5;

	s[0] = 1.0; s[1] = -1.0;

	neigh[0][0] = 2; neigh[0][1] = 2; neigh[0][2] = 2; neigh[0][3] = 2; neigh[0][4] = 2; neigh[0][5] = 2;
	neigh[1][0] = 1; neigh[1][1] = 1; neigh[1][2] = 1; neigh[1][3] = 1; neigh[1][4] = 1; neigh[1][5] = 1;

	bshift[0][0] = 0.0; bshift[0][1] = 0.0; bshift[0][2] = 0.0; bshift[0][3] = 0.0; bshift[0][4] = 0.0; bshift[0][5] = 0.0;
	bshift[1][0] = 0.0; bshift[1][1] = 0.0; bshift[1][2] = 0.0; bshift[1][3] = 0.0; bshift[1][4] = 0.0; bshift[1][5] = 0.0;

        BZPoint[0][0] = 2.0*PI; BZPoint[0][1] = PI;     BZPoint[0][2] = 0.0;
        BZPoint[1][0] = PI;     BZPoint[1][1] = PI;     BZPoint[1][2] = PI;
        BZPoint[2][0] = 0.0;    BZPoint[2][1] = 0.0;    BZPoint[2][2] = 0.0;
	BZPoint[3][0] = 2.0*PI; BZPoint[3][1] = 0.0;    BZPoint[3][2] = 0.0;
        BZPoint[4][0] = 2.0*PI; BZPoint[4][1] = PI/2.0; BZPoint[4][2] = PI/2.0;
        BZPoint[5][0] = 0.0;    BZPoint[5][1] = 0.0;    BZPoint[5][2] = 0.0;
        
	bond_pass = 0.0;

	Vol_atom = a0*a0*a0/8.0*1e-27;

	break;

    case 10:

        v[0][0] = 0.0; v[0][1] = 1.0/sqrt(3.0); v[0][2] = 0.0;
	v[1][0] = 0.0; v[1][1] = 0.0; 		v[1][2] = 0.5;
	v[2][0] = 0.0; v[2][1] = 1.0/sqrt(3.0); v[2][2] = 1.0;

	a[0][0] = 1.0/2.0;  a[0][1] = sqrt(3.0)/2.0; a[0][2] = 0.0;
	a[1][0] = -1.0/2.0; a[1][1] = sqrt(3.0)/2.0; a[1][2] = 0.0;
	a[2][0] = 0.0;      a[2][1] = 0.0;           a[2][2] = 1.5;

	b[0][0] = 0.0;  b[0][1] = 1.0/sqrt(3.0);      b[0][2] = 0.5;
	b[1][0] = 0.5;  b[1][1] = -1.0/(2*sqrt(3.0)); b[1][2] = 0.5;
	b[2][0] = -0.5; b[2][1] = -1.0/(2*sqrt(3.0)); b[2][2] = 0.5;
	b[3][0] = 0.0;  b[3][1] = 1.0/sqrt(3.0);      b[3][2] = -0.5;
	b[4][0] = 0.5;  b[4][1] = -1.0/(2*sqrt(3.0)); b[4][2] = -0.5;
	b[5][0] = -0.5; b[5][1] = -1.0/(2*sqrt(3.0)); b[5][2] = -0.5;

	s[0] = -1; s[1] = 1; s[2] = -1;

	neigh[0][0] = 0; neigh[0][1] = 0; neigh[0][2] = 0; neigh[0][3] = 2; neigh[0][4] = 2; neigh[0][5] = 2;
	neigh[1][0] = 3; neigh[1][1] = 3; neigh[1][2] = 3; neigh[1][3] = 1; neigh[1][4] = 1; neigh[1][5] = 1;
	neigh[2][0] = 2; neigh[2][1] = 2; neigh[2][2] = 2; neigh[2][3] = 0; neigh[2][4] = 0; neigh[2][5] = 0;

	bshift[0][0] = 0.0; bshift[0][1] = 0.0; bshift[0][2] = 0.0; bshift[0][3] = 0.0; bshift[0][4] = 0.0; bshift[0][5] = 0.0;
	bshift[1][0] = 0.0; bshift[1][1] = 0.0; bshift[1][2] = 0.0; bshift[1][3] = 0.0; bshift[1][4] = 0.0; bshift[1][5] = 0.0;
	bshift[2][0] = 0.0; bshift[2][1] = 0.0; bshift[2][2] = 0.0; bshift[2][3] = 0.0; bshift[2][4] = 0.0; bshift[2][5] = 0.0;

        BZPoint[0][0] = 2*PI/3.0; BZPoint[0][1] = 2*PI/sqrt(3.0); BZPoint[0][2] = 0.0;
        BZPoint[1][0] = 0.0;      BZPoint[1][1] = 0.0;    	  BZPoint[1][2] = 0.0;
        BZPoint[2][0] = PI;	  BZPoint[2][1] = PI/(sqrt(3.0)); BZPoint[2][2] = 0.0;
        BZPoint[3][0] = 2*PI/3.0; BZPoint[3][1] = 2*PI/sqrt(3.0); BZPoint[3][2] = 0.0;

	bond_pass = 0.0;

	Vol_atom = sqrt(3.0)/6*c0*a0*a0*1e-27;

        break;

    case 11:
      
        v[0][0] = 0.0; v[0][1] = 1.0/sqrt(3.0); v[0][2] = 0.0;
	v[1][0] = 0.0; v[1][1] = 0.0; 		v[1][2] = 0.5;
	v[2][0] = 0.0; v[2][1] = 1.0/sqrt(3.0); v[2][2] = 1.0;
        v[3][0] = 0.0; v[3][1] = 0.0;		v[3][2] = u0/2.0;
	v[4][0] = 0.0; v[4][1] = 1.0/sqrt(3.0);	v[4][2] = u0/2.0+0.5;
	v[5][0] = 0.0; v[5][1] = 0.0;		v[5][2] = u0/2.0+1.0;

	a[0][0] = 1.0/2.0;  a[0][1] = sqrt(3.0)/2.0; a[0][2] = 0.0;
	a[1][0] = -1.0/2.0; a[1][1] = sqrt(3.0)/2.0; a[1][2] = 0.0;
	a[2][0] = 0.0;      a[2][1] = 0.0;           a[2][2] = u0;

	b[0][0] = 0.0;  b[0][1] = 1.0/sqrt(3.0);      b[0][2] = 0.5;
	b[1][0] = 0.5;  b[1][1] = -1.0/(2*sqrt(3.0)); b[1][2] = 0.5;
	b[2][0] = -0.5; b[2][1] = -1.0/(2*sqrt(3.0)); b[2][2] = 0.5;
	b[3][0] = 0.0;  b[3][1] = 1.0/sqrt(3.0);      b[3][2] = -0.5;
	b[4][0] = 0.5;  b[4][1] = -1.0/(2*sqrt(3.0)); b[4][2] = -0.5;
	b[5][0] = -0.5; b[5][1] = -1.0/(2*sqrt(3.0)); b[5][2] = -0.5;

	s[0] = -1; s[1] = 1; s[2] = -1; s[3] = 1; s[4] = -1; s[5] = 1;

	neigh[0][0] = 0; neigh[0][1] = 0; neigh[0][2] = 0; neigh[0][3] = 2; neigh[0][4] = 2; neigh[0][5] = 2;
	neigh[1][0] = 3; neigh[1][1] = 3; neigh[1][2] = 3; neigh[1][3] = 1; neigh[1][4] = 1; neigh[1][5] = 1;
	neigh[2][0] = 2; neigh[2][1] = 2; neigh[2][2] = 2; neigh[2][3] = 0; neigh[2][4] = 0; neigh[2][5] = 0;
	neigh[3][0] = 5; neigh[3][1] = 5; neigh[3][2] = 5; neigh[3][3] = 0; neigh[3][4] = 0; neigh[3][5] = 0;
	neigh[4][0] = 4; neigh[4][1] = 4; neigh[4][2] = 4; neigh[4][3] = 6; neigh[4][4] = 6; neigh[4][5] = 6;
	neigh[5][0] = 0; neigh[5][1] = 0; neigh[5][2] = 0; neigh[5][3] = 5; neigh[5][4] = 5; neigh[5][5] = 5;

	bshift[0][0] = 0.0; bshift[0][1] = 0.0; bshift[0][2] = 0.0; bshift[0][3] = 0.0; bshift[0][4] = 0.0; bshift[0][5] = 0.0;
	bshift[1][0] = 0.0; bshift[1][1] = 0.0; bshift[1][2] = 0.0; bshift[1][3] = 0.0; bshift[1][4] = 0.0; bshift[1][5] = 0.0;
	bshift[2][0] = 0.0; bshift[2][1] = 0.0; bshift[2][2] = 0.0; bshift[2][3] = 0.0; bshift[2][4] = 0.0; bshift[2][5] = 0.0;
	bshift[3][0] = 0.0; bshift[3][1] = 0.0; bshift[3][2] = 0.0; bshift[3][3] = 0.0; bshift[3][4] = 0.0; bshift[3][5] = 0.0;
	bshift[4][0] = 0.0; bshift[4][1] = 0.0; bshift[4][2] = 0.0; bshift[4][3] = 0.0; bshift[4][4] = 0.0; bshift[4][5] = 0.0;
	bshift[5][0] = 0.0; bshift[5][1] = 0.0; bshift[5][2] = 0.0; bshift[5][3] = 0.0; bshift[5][4] = 0.0; bshift[5][5] = 0.0;

        BZPoint[0][0] = 2*PI/3.0; BZPoint[0][1] = 2*PI/sqrt(3.0); BZPoint[0][2] = 0.0;
        BZPoint[1][0] = 0.0;      BZPoint[1][1] = 0.0;    	  BZPoint[1][2] = 0.0;
        BZPoint[2][0] = PI;	  BZPoint[2][1] = PI/(sqrt(3.0)); BZPoint[2][2] = 0.0;
        BZPoint[3][0] = 2*PI/3.0; BZPoint[3][1] = 2*PI/sqrt(3.0); BZPoint[3][2] = 0.0;

	bond_pass = 0.0;

	Vol_atom = sqrt(3.0)/12*c0*a0*a0*1e-27;

        break;
    }
        
    if (strcmp(first_atom,"cation")==0){
        for(IA=0;IA<NA;IA++){
	    for(ID=0;ID<N3D;ID++){
	        unit_cell->atom[IA]->coord[ID] = a_factor[ID]*(v[IA][ID]-v[1][ID]);
            }
        }   
    }else{
        for(IA=0;IA<NA;IA++){
	    for(ID=0;ID<N3D;ID++){
	        unit_cell->atom[IA]->coord[ID] = a_factor[ID]*v[IA][ID];
            }
        }
    }

    for(IA=0;IA<NA;IA++){
        unit_cell->atom[IA]->type = IA+1;
	init_var(unit_cell->atom[IA]->int_disp,N3D);
    }
    
    for(ID=0;ID<N3D;ID++){
        for(JD=0;JD<N3D;JD++){
            unit_cell->axis[ID]->vec[JD] = a_factor[JD]*a[ID][JD];
        }
    }
    
    for(IA=0;IA<NA;IA++){
        for (IB=0;IB<NB;IB++){
            for(ID=0;ID<N3D;ID++){
                unit_cell->type[IA]->bond[IB]->vec[ID] = a_factor[ID]*s[IA]*b[IB][ID];
            }
	    unit_cell->type[IA]->bond[IB]->vec[2] = unit_cell->type[IA]->bond[IB]->vec[2]+\
	      a_factor[2]*s[IA]*bshift[IA][IB];
	    unit_cell->type[IA]->bond[IB]->neigh = neigh[IA][IB];
        }
    }

    init_var(unit_cell->strain_matrix,N3D*N3D);
    for(ID=0;ID<N3D;ID++){
        unit_cell->strain_matrix[ID*(N3D+1)] = 1.0;
    }

    for(IP=0;IP<unit_cell->BZNP;IP++){
        for(ID=0;ID<N3D;ID++){
            unit_cell->BZPoint[IP]->coord[ID] = BZPoint[IP][ID]/a_factor[ID];        
        }    
    }

    unit_cell->bond_pass      = bond_pass;
    unit_cell->passivation[0] = 1.0/2.0;
    unit_cell->passivation[1] = sqrt(3.0)/2.0;
    unit_cell->passivation[2] = sqrt(3.0)/2.0;
    unit_cell->passivation[3] = sqrt(3.0)/2.0;

}

/************************************************************************************************/

void WireGenerator::change_orientation(double* x,double* y,double* z,double dsp3)
{
    double x_norm=0.0,y_norm=0.0,z_norm=0.0;
    double R[9]={0.0},inv_R[9];
    double slmn[4];
    double vec_trans[3];
    double factor;
    int IPIV[3],INFO;
    int IA,IB,IB1,IB2,ID,IP;
    int index;

    x_norm = c_dnrm2(N3D,x,1);
    c_dscal(N3D,1/x_norm,x,1);
    y_norm = c_dnrm2(N3D,y,1);
    c_dscal(N3D,1/y_norm,y,1);
    z_norm = c_dnrm2(N3D,z,1);
    c_dscal(N3D,1/z_norm,z,1);
    
    for(ID=0;ID<N3D;ID++){
        inv_R[ID]       = x[ID];
        inv_R[ID+N3D]   = y[ID];
        inv_R[ID+2*N3D] = z[ID];
        R[ID*(N3D+1)] = 1.0;
    }
    c_dgetrf(N3D,N3D,inv_R,N3D,IPIV,&INFO);
    c_dgetrs('N',N3D,N3D,inv_R,N3D,IPIV,R,N3D,&INFO);
    
    for(IA=0;IA<NA;IA++){

        for(IB=0;IB<NB;IB++){
            
            c_dcopy(N3D,unit_cell->type[IA]->bond[IB]->vec,1,vec_trans,1);
            c_dgemv('N',N3D,N3D,1.0,R,N3D,vec_trans,1,0.0,\
                  unit_cell->type[IA]->bond[IB]->vec,1);
            c_dcopy(N3D,unit_cell->type[IA]->bond[IB]->vec,1,\
                  unit_cell->type_original[IA]->bond[IB]->vec,1);
            unit_cell->type_original[IA]->bond[IB]->d0 =\
                c_dnrm2(N3D,unit_cell->type_original[IA]->bond[IB]->vec,1);
            
	    slmn[0] = 1;
	    slmn[1] = unit_cell->type[IA]->bond[IB]->vec[0]/\
	      c_dnrm2(N3D,unit_cell->type[IA]->bond[IB]->vec,1);
	    slmn[2] = unit_cell->type[IA]->bond[IB]->vec[1]/\
	      c_dnrm2(N3D,unit_cell->type[IA]->bond[IB]->vec,1);
	    slmn[3] = unit_cell->type[IA]->bond[IB]->vec[2]/\
	      c_dnrm2(N3D,unit_cell->type[IA]->bond[IB]->vec,1);
	    
	    for(IB1=0;IB1<SP;IB1++){
                for(IB2=0;IB2<SP;IB2++){
                    unit_cell->type[IA]->bond[IB]->SP3[IB1*SP+IB2] =\
		      unit_cell->bond_pass*(dsp3*unit_cell->passivation[IB1]*slmn[IB1]*\
					    unit_cell->passivation[IB2]*slmn[IB2]+5*tollim);
		}
	    }
        }
        
        c_dcopy(N3D,unit_cell->atom[IA]->coord,1,vec_trans,1);
        c_dgemv('N',N3D,N3D,1.0,R,N3D,vec_trans,1,0.0,\
              unit_cell->atom[IA]->coord,1);
    }
    
    for(ID=0;ID<N3D;ID++){
        c_dcopy(N3D,unit_cell->axis[ID]->vec,1,vec_trans,1);
        c_dgemv('N',N3D,N3D,1.0,R,N3D,vec_trans,1,0.0,\
              unit_cell->axis[ID]->vec,1);
    }

    for(IP=0;IP<unit_cell->BZNP;IP++){
        c_dcopy(N3D,unit_cell->BZPoint[IP]->coord,1,vec_trans,1);        
        c_dgemv('N',N3D,N3D,1.0,R,N3D,vec_trans,1,0.0,\
                unit_cell->BZPoint[IP]->coord,1);
    }


    for(ID=0;ID<N3D;ID++){
        c_dcopy(N3D,unit_cell->axis[ID]->vec,1,&inv_R[N3D*ID],1);
    }
    
    init_var(vec_trans,N3D);
    vec_trans[0] = 1.0;
    
    c_dgetrf(N3D,N3D,inv_R,N3D,IPIV,&INFO);
    c_dgetrs('N',N3D,1,inv_R,N3D,IPIV,vec_trans,N3D,&INFO);

    factor = INF;
    for(ID=0;ID<N3D;ID++){
        if((abs(vec_trans[ID])<factor)&&(abs(vec_trans[ID])>tollim)){
            factor = vec_trans[ID];
            index  = ID;
        }
    }
    
    for(ID=0;ID<N3D;ID++){
        if(ID!=index){
            c_daxpy(N3D,vec_trans[ID]/factor,unit_cell->axis[ID]->vec,1,\
                    unit_cell->axis[index]->vec,1);
        }
    }

}

/************************************************************************************************/

void WireGenerator::get_d_matrix(double *d_matrix,double l,double m,double n)
{

    int ID1      = 0;
    int ID2      = 1;
    int ID3      = 2;
    int ID4      = 3;
    int ID5      = 4;
    int NDorb    = 5;
    double l2    = l*l;
    double m2    = m*m;
    double n2    = n*n;

    if(abs(n)<1-tollim){

        d_matrix[ID1+NDorb*ID1] = n*(l2-m2)/(1.0-n2);
	d_matrix[ID2+NDorb*ID1] = -l;
	d_matrix[ID3+NDorb*ID1] = m;
	d_matrix[ID4+NDorb*ID1] = -2.0*l*m*n/(1.0-n2);
	d_matrix[ID5+NDorb*ID1] = 0.0;
	d_matrix[ID1+NDorb*ID2] = (m2-l2)/sqrt(1.0-n2);
	d_matrix[ID2+NDorb*ID2] = -n*l/sqrt(1.0-n2);
	d_matrix[ID3+NDorb*ID2] = n*m/sqrt(1.0-n2);
	d_matrix[ID4+NDorb*ID2] = 2.0*l*m/sqrt(1.0-n2);
	d_matrix[ID5+NDorb*ID2] = 0.0;
	d_matrix[ID1+NDorb*ID3] = -2.0*l*m*n/sqrt(1.0-n2);
	d_matrix[ID2+NDorb*ID3] = m*(-2.0*n2+1.0)/sqrt(1.0-n2);
	d_matrix[ID3+NDorb*ID3] = l*(1-2.0*n2)/sqrt(1.0-n2);
	d_matrix[ID4+NDorb*ID3] = -n*(l2-m2)/sqrt(1.0-n2);
	d_matrix[ID5+NDorb*ID3] = sqrt(3.0)*n*sqrt(1.0-n2);
	d_matrix[ID1+NDorb*ID4] = l*m*(1.0+n2)/(1.0-n2);
	d_matrix[ID2+NDorb*ID4] = -m*n;
	d_matrix[ID3+NDorb*ID4] = -l*n;
	d_matrix[ID4+NDorb*ID4] = (1.0+n2)*(l2-m2)/(2.0*(1.0-n2));
	d_matrix[ID5+NDorb*ID4] = (1.0-n2)/2.0*sqrt(3.0);
	d_matrix[ID1+NDorb*ID5] = m*l*sqrt(3.0);
	d_matrix[ID2+NDorb*ID5] = m*n*sqrt(3.0);
	d_matrix[ID3+NDorb*ID5] = l*n*sqrt(3.0);
	d_matrix[ID4+NDorb*ID5] = sqrt(3.0)*(l2-m2)/2.0;
	d_matrix[ID5+NDorb*ID5] = (3.0*n2-1.0)/2.0;

    }else{

        init_var(d_matrix,NDorb*NDorb);

	d_matrix[ID1+NDorb*ID1] = 1.0;
	d_matrix[ID2+NDorb*ID2] = 1.0;
	d_matrix[ID3+NDorb*ID3] = 1.0;
	d_matrix[ID4+NDorb*ID4] = 1.0;
	d_matrix[ID5+NDorb*ID5] = 1.0;

    }
    
}

/************************************************************************************************/

void WireGenerator::make_strained_unit_cell(Strain *strain)
{
    int ID,IP,IPIV[3],INFO;
    double vec_trans[3];
    int IA,IB,i,j,k,l;
    double strain_matrix[9]={1.0+strain->Eps_xx,strain->Eps_xy,strain->Eps_xz, \
			     strain->Eps_xy,1.0+strain->Eps_yy,strain->Eps_yz, \
			     strain->Eps_xz,strain->Eps_yz,1.0+strain->Eps_zz};
    double inv_strain_matrix[9],WORK[9];
    double A[9],int_disp[3];
    double bfactor[2]={1.0,-1.0};
    double afactor[2]={0.0,1.0};
    
    c_dcopy(N3D*N3D,strain_matrix,1,unit_cell->strain_matrix,1);

    init_var(int_disp,N3D);

    /*Internal displacement*/

    if(lattice_type==1){

        for(i=0;i<N3D;i++){

	    for(j=0;j<N3D;j++){

	        A[i+j*N3D] = 0.0;

		for(k=0;k<N3D;k++){

		    A[i+j*N3D]  = A[i+j*N3D]+2*(unit_cell->type[0]->bond[0]->vec[k]- \
						unit_cell->type[0]->bond[1+i]->vec[k])*strain_matrix[j+k*N3D];
		
		    for(l=0;l<N3D;l++){

		        int_disp[i] = int_disp[i]-(unit_cell->type[0]->bond[0]->vec[j]*strain_matrix[j+k*N3D]*strain_matrix[k+l*N3D]*unit_cell->type[0]->bond[0]->vec[l]- \
						   unit_cell->type[0]->bond[1+i]->vec[j]*strain_matrix[j+k*N3D]*strain_matrix[k+l*N3D]*unit_cell->type[0]->bond[1+i]->vec[l])*strain->zeta;
		    }
		}
	    }
	}

	c_dgetrf(N3D,N3D,A,N3D,IPIV,&INFO);
	c_dgetrs('N',N3D,1,A,N3D,IPIV,int_disp,N3D,&INFO);
    }
    /*End*/

    for(IA=0;IA<NA;IA++){

        for(IB=0;IB<NB;IB++){
	    c_dcopy(N3D,unit_cell->type[IA]->bond[IB]->vec,1,vec_trans,1);
	    c_dgemv('N',N3D,N3D,1.0,strain_matrix,N3D,vec_trans,1,0.0,\
		    unit_cell->type[IA]->bond[IB]->vec,1);
	    c_daxpy(N3D,bfactor[IA],int_disp,1,unit_cell->type[IA]->bond[IB]->vec,1);
	}

	c_dcopy(N3D,unit_cell->atom[IA]->coord,1,vec_trans,1);
	c_dgemv('N',N3D,N3D,1.0,strain_matrix,N3D,vec_trans,1,0.0,\
		unit_cell->atom[IA]->coord,1);
	c_daxpy(N3D,afactor[IA],int_disp,1,unit_cell->atom[IA]->coord,1);
	c_daxpy(N3D,afactor[IA],int_disp,1,unit_cell->atom[IA]->int_disp,1);
    }

    for(ID=0;ID<N3D;ID++){
        c_dcopy(N3D,unit_cell->axis[ID]->vec,1,vec_trans,1);
	c_dgemv('N',N3D,N3D,1.0,strain_matrix,N3D,vec_trans,1,0.0,\
		unit_cell->axis[ID]->vec,1);
    }

    c_dcopy(N3D*N3D,strain_matrix,1,inv_strain_matrix,1);
    c_dgetrf(N3D,N3D,inv_strain_matrix,N3D,IPIV,&INFO);
    c_dgetri(N3D,inv_strain_matrix,N3D,IPIV,WORK,N3D*N3D,&INFO);
    
    for(IP=0;IP<unit_cell->BZNP;IP++){
        c_dcopy(N3D,unit_cell->BZPoint[IP]->coord,1,vec_trans,1);        
        c_dgemv('N',N3D,N3D,1.0,inv_strain_matrix,N3D,vec_trans,1,0.0,\
                unit_cell->BZPoint[IP]->coord,1);
    }
}

/************************************************************************************************/

double WireGenerator::calc_tri_area(double* p1, double* p2, double* p3, int ND)
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
        delete [] v1;
    }
    if (v2!=NULL){
         delete [] v2;
    }
    return A;
}

/************************************************************************************************/

double WireGenerator::calc_quad_area(double* p1, double* p2, double* p3, double* p4, int ND)
{
    return calc_tri_area(p1,p2,p4,ND)+calc_tri_area(p3,p2,p4,ND);
}

/************************************************************************************************/

void WireGenerator::calc_vec(double* v, double* p1, double* p2, double* p3)
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

double WireGenerator::calc_tri_volume(double A, double* p0, double* p1, double* v)
{
    double p1_p0[3];

    c_dcopy(N3D,p1,1,p1_p0,1);
    c_daxpy(N3D,-1.0,p0,1,p1_p0,1);

    return fabs(c_ddot(N3D,v,1,p1_p0,1)*A/3.0);
}

/************************************************************************************************/

double WireGenerator::calc_quad_volume(MAT* mat)
{
    double V1,V2,V3;

    V1 = calc_tri_volume(mat->face_area[1],mat->p[0]->coord,mat->p_ref[1],mat->vec_dir[1]);
    V2 = calc_tri_volume(mat->face_area[4],mat->p[0]->coord,mat->p_ref[4],mat->vec_dir[4]);
    V3 = calc_tri_volume(mat->face_area[5],mat->p[0]->coord,mat->p_ref[5],mat->vec_dir[5]);

    return V1+V2+V3;
}

/************************************************************************************************/

int WireGenerator::is_in_quad3D(MAT* mat, double* p)
{
    int is_in = 0;

    if(!strcmp(mat->type,"square")){

        double V = 0.0;
        int i    = 0;
        
        for(i=0;i<6;i++){
            V = V+calc_tri_volume(mat->face_area[i],p,mat->p_ref[i],mat->vec_dir[i]);
        }

        if((fabs(V-mat->volume)/V)<tollim) is_in = 1;
        
    }

    if(!strcmp(mat->type,"circle")){
        
        bool condition1 = (p[0]+tollim)>=mat->p[0]->coord[0];
        bool condition2 = (p[0]-tollim)<=mat->p[1]->coord[0];
        bool condition3 = (pow(p[1]-mat->p[0]->coord[1],2.0)/pow(mat->radius[0],2.0)+\
			   pow(p[2]-mat->p[0]->coord[2],2.0)/pow(mat->radius[1],2.0))<=1.0;

        is_in = condition1&&condition2&&condition3;
        
    }

    if(!strcmp(mat->type,"sphere")){
        
        is_in = (pow(p[0]-mat->p[0]->coord[0],2.0)/pow(mat->radius[0],2.0)+\
		 pow(p[1]-mat->p[0]->coord[1],2.0)/pow(mat->radius[0],2.0)+\
                 pow(p[2]-mat->p[0]->coord[2],2.0)/pow(mat->radius[0],2.0))<=1.0;
        
    }

    return is_in;
}

/************************************************************************************************/

int WireGenerator::is_in_quad2D(SURFACE* surf, double* p)
{
    int is_in = 0;

    if(!strcmp(surf->type,"square")){

        double A = 0.0;
        int i    = 0;

        for(i=0;i<3;i++){
            A=A+calc_tri_area(surf->p[i]->coord,surf->p[i+1]->coord,p,N2D);
        }
        A=A+calc_tri_area(surf->p[3]->coord,surf->p[0]->coord,p,N2D);

        if((fabs(A-surf->face_area)/A)<tollim) is_in = 1;
        
    }else{

        double r2,d2,phi;

	phi   = atan(surf->radius[0]/surf->radius[1]*\
		     (p[2]-surf->p[0]->coord[1])/(p[1]-surf->p[0]->coord[0]+tollim));
        r2    = pow(p[1]-surf->p[0]->coord[0],2.0)+pow(p[2]-surf->p[0]->coord[1],2.0);
	d2    = pow(surf->radius[0]*cos(phi),2.0)+pow(surf->radius[1]*sin(phi),2.0);

	is_in = r2<=d2;

	/*
        is_in = (pow(p[0]-surf->p[0]->coord[0],2.0)+pow(p[1]-surf->p[0]->coord[1],\
                                                        2.0))<=pow(surf->radius,2.0);
	*/
    }

    return is_in;
}

/************************************************************************************************/

void WireGenerator::sort_xyz(XYZPOS* xyz, int *index, int* NL, int* LMI, int* LMA)
{
    
    int ypos0 = 0;
    int ypos1 = 1;
    int zpos0,zpos1,zlimit,IA;

    NL[0] = 0;
    sort(xyz,xyz+No_Atom,sortx);
    
    while(ypos1<No_Atom){
        while(xyz[ypos1].x<=(xyz[ypos0].x+tollim)){
            ypos1 = ypos1+1;
            if (ypos1>=No_Atom){
                break;
            }
        }
        sort(&xyz[ypos0],&xyz[min(ypos1,No_Atom)],sorty);
        LMI[*NL] = ypos0;
        LMA[*NL] = ypos1-1;
        *NL      = *NL+1;
        
        zpos0    = ypos0;
        zpos1    = ypos0+1;
        zlimit   = min(ypos1,No_Atom);
       
        while(zpos1<zlimit){
            while(xyz[zpos1].y<=(xyz[zpos0].y+tollim)){
                zpos1=zpos1+1;
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
        
        if(ypos1<No_Atom){
            ypos0 = ypos1;
            ypos1 = ypos0+1;
        }
    }

    for(IA=0;IA<No_Atom;IA++){
        index[IA] = xyz[IA].index;
    }
    
}

/************************************************************************************************/

void WireGenerator::cut_boundary_slab(WireStructure* nanowire)
{
    int i,ib,IA,IB,IN,index[2]={0,NSlab-1},Sshift,no_neighbor,neighbor,neigh_ib;
    double at_pos0[3],at_pos[3],neigh_pos[3],bond_vec[3];
    double max_bond_deformation = nanowire->max_bond_deformation,deformation;

    deformation = max_bond_deformation*\
        c_dnrm2(N3D,unit_cell->type_original[0]->bond[0]->vec,1);
    
    for(i=0;i<2;i++){

        Sshift          = Smin[index[i]];
        Boundary[i]->NC = 0;
        
        for(IA=Smin[index[i]];IA<=Smax[index[i]];IA++){    
       
            c_dcopy(N3D+1,&Layer_Matrix[SLM*IA],1,&Boundary[i]->UNN[SLM*(IA-Sshift)],1);
            c_dcopy(N3D+1,&Layer_Matrix[SLM*IA],1,&Boundary[i]->UNNp[SLM*(IA-Sshift)],1);
            c_dcopy(N3D+1,&Layer_Matrix[SLM*IA],1,&Boundary[i]->UNNm[SLM*(IA-Sshift)],1);
            c_dcopy(N3D,&Layer_Matrix[SLM*IA],1,at_pos0,1);
            
            Boundary[i]->UNNp[SLM*(IA-Sshift)]        = Boundary[i]->UNNp[SLM*(IA-Sshift)]+\
                Boundary[i]->cell_width;
            Boundary[i]->UNNm[SLM*(IA-Sshift)]        = Boundary[i]->UNNm[SLM*(IA-Sshift)]-\
                Boundary[i]->cell_width;
            Boundary[i]->UNNneigh[(NB+1)*(IA-Sshift)] = Round(Layer_Matrix[SLM*IA+3]);
            at_pos0[0]                                = at_pos0[0]-Boundary[i]->cell_width;
	    
            for(IB=0;IB<NB;IB++){
                
                Boundary[i]->UNN[SLM*(IA-Sshift)+4+IB]         = 0;
                Boundary[i]->UNNp[SLM*(IA-Sshift)+4+IB]        = 0;
                Boundary[i]->UNNm[SLM*(IA-Sshift)+4+IB]        = 0;
                Boundary[i]->UNNneigh[(NB+1)*(IA-Sshift)+1+IB] = 0;
                
                if((Layer_Matrix[SLM*IA+4+IB]-1)>=Smin[index[i]]&\
                   (Layer_Matrix[SLM*IA+4+IB]-1)<=Smax[index[i]]){
                    
                    Boundary[i]->UNN[SLM*(IA-Sshift)+4+IB] = Layer_Matrix[SLM*IA+4+IB]-Sshift;
                    Boundary[i]->UNNneigh[(NB+1)*(IA-Sshift)+1+IB]=\
                            Round(Layer_Matrix[SLM*Round(Layer_Matrix[SLM*IA+4+IB]-1)+3]);
                    Boundary[i]->NC = Boundary[i]->NC+1;
                }else{
                    c_dcopy(N3D,at_pos0,1,at_pos,1);
                    c_daxpy(N3D,1.0,unit_cell->type[Round(Layer_Matrix[SLM*IA+3]-1)%NA]->\
                          bond[IB]->vec,1,at_pos,1);
                    no_neighbor = 0;
                
                    for(IN=Smin[index[i]];IN<=Smax[index[i]];IN++){

                        c_dcopy(N3D,at_pos,1,neigh_pos,1);
                        c_daxpy(N3D,-1.0,&Layer_Matrix[SLM*IN],1,neigh_pos,1);
                        
                        if(c_dnrm2(N3D,neigh_pos,1)<=deformation){
                            no_neighbor = no_neighbor+1;
                            if(no_neighbor>1){
                                printf("Please reduce max_bond_deformation\n");
                                exit(0);
                            }else{
                                neighbor = IN-Sshift;
                            }
                        }    
                    }
                
                    if(no_neighbor==1){
                        Boundary[i]->UNNp[SLM*(IA-Sshift)+4+IB] = neighbor+1.0;
                        Boundary[i]->UNNneigh[(NB+1)*(IA-Sshift)+1+IB] =\
                            Round(Layer_Matrix[SLM*(neighbor+Sshift)+3]);
                        Boundary[i]->NC = Boundary[i]->NC+1;
                        Boundary[i]->UNN[SLM*(IA-Sshift)+4+IB] = -1;
                        if(IA>=Smin[NSlab-1]){
			    if(nanowire->open_system){
			        Layer_Matrix[SLM*IA+4+IB] = -1;
			    }else{
			        Layer_Matrix[SLM*IA+4+IB] = 0;
			    }
                        }

			for(ib=0;ib<NB;ib++){
			    c_dcopy(N3D,unit_cell->type[Round(Layer_Matrix[SLM*IA+3]-1)%NA]->bond[IB]->vec,1,\
				    bond_vec,1);
			    c_daxpy(N3D,1.0,unit_cell->type[Round(Layer_Matrix[SLM*neighbor+3]-1)%NA]->bond[ib]->vec,\
				    1,bond_vec,1);
			    if(c_dnrm2(N3D,bond_vec,1)<=tollim){
			        neigh_ib = ib;
				break;
			    }
			}
                        
                        Boundary[i]->UNNm[SLM*neighbor+4+neigh_ib] = IA-Sshift+1;
                        Boundary[i]->UNNneigh[(NB+1)*neighbor+1+neigh_ib] =\
                            Round(Layer_Matrix[SLM*IA+3]);
                        Boundary[i]->NC = Boundary[i]->NC+1;
                        Boundary[i]->UNN[SLM*neighbor+4+neigh_ib] = -1;
                        if(IA<=Smax[0]){
			    if(nanowire->open_system){
			        Layer_Matrix[SLM*(neighbor+Sshift)+4+neigh_ib] = -1;
			    }else{
			        Layer_Matrix[SLM*(neighbor+Sshift)+4+neigh_ib] = 0;
			    }
                        }
                    }
                }
            }
        }
    }
  
    if(lattice_type==6){
        roll_cnt();
    }

 }

/************************************************************************************************/

void WireGenerator::update_boundary_slab(WireStructure* nanowire)
{
    int i,IA,IB,index[2]={0,NSlab-1},Sshift;
   
    for(i=0;i<2;i++){
        
        Sshift = Smin[index[i]];
        
        for(IA=Smin[index[i]];IA<=Smax[index[i]];IA++){    
       
            c_dcopy(N3D+1,&Layer_Matrix[SLM*IA],1,&Boundary[i]->UNN[SLM*(IA-Sshift)],1);
            c_dcopy(N3D+1,&Layer_Matrix[SLM*IA],1,&Boundary[i]->UNNp[SLM*(IA-Sshift)],1);
            c_dcopy(N3D+1,&Layer_Matrix[SLM*IA],1,&Boundary[i]->UNNm[SLM*(IA-Sshift)],1);
            
            Boundary[i]->UNNp[SLM*(IA-Sshift)]        = Boundary[i]->UNNp[SLM*(IA-Sshift)]+\
                Boundary[i]->cell_width;
            Boundary[i]->UNNm[SLM*(IA-Sshift)]        = Boundary[i]->UNNm[SLM*(IA-Sshift)]-\
                Boundary[i]->cell_width;
            Boundary[i]->UNNneigh[(NB+1)*(IA-Sshift)] = Round(Layer_Matrix[SLM*IA+3]);

	}

	for(IA=Smin[index[i]];IA<=Smax[index[i]];IA++){
	   
            for(IB=0;IB<NB;IB++){

	        if(Layer_Matrix[SLM*IA+4+IB]>0){
		    Boundary[i]->UNNneigh[(NB+1)*(IA-Sshift)+1+IB]=\
		      Round(Layer_Matrix[SLM*Round(Layer_Matrix[SLM*IA+4+IB]-1)+3]);
		}

		if(Layer_Matrix[SLM*IA+4+IB]==-1){
		    if(Boundary[i]->UNNp[SLM*(IA-Sshift)+4+IB]>0){
		        Boundary[i]->UNNneigh[(NB+1)*(IA-Sshift)+1+IB] = \
			  Round(Boundary[i]->UNN[SLM*Round(Boundary[i]->UNNp[SLM*(IA-Sshift)+4+IB]-1)+3]);
		    }else{
		        Boundary[i]->UNNneigh[(NB+1)*(IA-Sshift)+1+IB] = \
			  Round(Boundary[i]->UNN[SLM*Round(Boundary[i]->UNNm[SLM*(IA-Sshift)+4+IB]-1)+3]);
		    }
		}
	    }
	}
    }
      
 }

/************************************************************************************************/

void WireGenerator::cut_layer(WireStructure* nanowire)
{
    int IL,IA,IB,Lshift;

    for(IL=0;IL<NLayer;IL++){
        
        Lshift      = Lmin[IL];
        Enn[IL]->NC = 0;
        
        for(IA=Lmin[IL];IA<=Lmax[IL];IA++){
	  
            c_dcopy(N3D+1,&Layer_Matrix[SLM*IA],1,&Enn[IL]->UNN[SLM*(IA-Lshift)],1);
            c_dcopy(N3D+1,&Layer_Matrix[SLM*IA],1,&Enn[IL]->UNNp[SLM*(IA-Lshift)],1);
            Enn[IL]->UNNneigh[(NB+1)*(IA-Lshift)] = Round(Layer_Matrix[SLM*IA+3]);

            Neigh_Matrix[(NB+1)*IA] = Round(Layer_Matrix[SLM*IA+3]);

            for(IB=0;IB<NB;IB++){
	      
                Enn[IL]->UNN[SLM*(IA-Lshift)+4+IB]      = 0;
                Enn[IL]->UNNp[SLM*(IA-Lshift)+4+IB]     = 0;
                Enn[IL]->UNNneigh[(NB+1)*(IA-Lshift)+1+IB] = 0;
              
                Neigh_Matrix[(NB+1)*IA+1+IB]            = 0;

                if((Layer_Matrix[SLM*IA+4+IB]-1)>=Lmin[IL]&\
                   (Layer_Matrix[SLM*IA+4+IB]-1)<=Lmax[IL]){
                    Enn[IL]->UNN[SLM*(IA-Lshift)+4+IB]=Layer_Matrix[SLM*IA+4+IB]-Lshift;
                    Enn[IL]->UNNneigh[(NB+1)*(IA-Lshift)+1+IB]=\
                            Round(Layer_Matrix[SLM*Round(Layer_Matrix[SLM*IA+4+IB]-1)+3]);
                    Enn[IL]->NC = Enn[IL]->NC+1;
                    
                    Neigh_Matrix[(NB+1)*IA+1+IB]=\
                        Round(Layer_Matrix[SLM*Round(Layer_Matrix[SLM*IA+4+IB]-1)+3]);
                }else{
                    if((Layer_Matrix[SLM*IA+4+IB]-1)>Lmax[IL]){
                        Enn[IL]->UNNp[SLM*(IA-Lshift)+4+IB]        = \
			  Layer_Matrix[SLM*IA+4+IB]-Lmax[IL]-1;
                        Enn[IL]->UNN[SLM*(IA-Lshift)+4+IB]         = -1;
                        Enn[IL]->UNNneigh[(NB+1)*(IA-Lshift)+1+IB] = \
                            Round(Layer_Matrix[SLM*Round(Layer_Matrix[SLM*IA+4+IB]-1)+3]);
                        Enn[IL]->NC = Enn[IL]->NC+1;

                        Neigh_Matrix[(NB+1)*IA+1+IB] = \
                            Round(Layer_Matrix[SLM*Round(Layer_Matrix[SLM*IA+4+IB]-1)+3]);
                    }else{
                        if((Layer_Matrix[SLM*IA+4+IB]-1)<Lmin[IL]&\
                            (Layer_Matrix[SLM*IA+4+IB]-1)>=0){
                            Enn[IL]->UNN[SLM*(IA-Lshift)+4+IB]         = -1;
                            Enn[IL]->UNNneigh[(NB+1)*(IA-Lshift)+1+IB] = \
                                Round(Layer_Matrix[SLM*Round(Layer_Matrix[SLM*IA+4+IB]-1)+3]);
                            Enn[IL]->NC = Enn[IL]->NC+1;

                            Neigh_Matrix[(NB+1)*IA+1+IB] = \
                                Round(Layer_Matrix[SLM*Round(Layer_Matrix[SLM*IA+4+IB]-1)+3]);
                        }
                    }
                }
                if(IA<=Smax[0]){
                    if(Boundary[0]->UNNm[SLM*(IA-Smin[0])+4+IB]>0){
                        Enn[IL]->UNN[SLM*(IA-Lshift)+4+IB]         = -1;
                        Enn[IL]->UNNneigh[(NB+1)*(IA-Lshift)+1+IB] = \
                                Round(Boundary[0]->UNN[SLM*Round(Boundary[0]->\
                                                         UNNm[SLM*(IA-Smin[0])+4+IB]-1)+3]);
                        Enn[IL]->NC = Enn[IL]->NC+1;

			if(nanowire->open_system){
			    Neigh_Matrix[(NB+1)*IA+1+IB] =\
			      Round(Boundary[0]->UNN[SLM*Round(Boundary[0]-> \
							       UNNm[SLM*(IA-Smin[0])+4+IB]-1)+3]);
			}
                    }
                }
                if(IA>=Smin[NSlab-1]){
                    if(Boundary[1]->UNNp[SLM*(IA-Smin[NSlab-1])+4+IB]>0){
                        Enn[IL]->UNN[SLM*(IA-Lshift)+4+IB]         = -1;
                        Enn[IL]->UNNneigh[(NB+1)*(IA-Lshift)+1+IB] = \
                                Round(Boundary[1]->UNN[SLM*Round(Boundary[1]->\
                                                         UNNp[SLM*(IA-Smin[NSlab-1])+\
                                                             4+IB]-1)+3]);
                        Enn[IL]->NC = Enn[IL]->NC+1;

			if(nanowire->open_system){
			    Neigh_Matrix[(NB+1)*IA+1+IB] =\
			      Round(Boundary[1]->UNN[SLM*Round(Boundary[1]-> \
							       UNNp[SLM*(IA-Smin[NSlab-1])+ \
								    4+IB]-1)+3]);
			}
                    }
                }
            }
        }
    }
    
}

/************************************************************************************************/

void WireGenerator::separate_dimension(WireStructure *nanowire)
{

    int IB,IA,factor;
    int NCH,NACH,NBOUND;

    init_var(index_channel,No_Atom);
    c_icopy(NSlab,Smin,1,Smin_tot,1);
    c_icopy(NSlab,Smax,1,Smax_tot,1);

    if(nanowire->NDim==1){
        
        for(IA=Smin[1];IA<Smin[NSlab-1];IA++){
	    if((Layer_Matrix[SLM*IA+1]>-100*tollim)&&(Layer_Matrix[SLM*IA+1]<(Ly-100*tollim))){
	        index_channel[IA] = 0;
            }
            else{
                if(Layer_Matrix[SLM*IA+1]<-100*tollim){
                    factor = -3;
                }
                if(Layer_Matrix[SLM*IA+1]>(Ly-100*tollim)){
                    factor = 3;
                }
                index_channel[IA] = factor;
           }
        }

    }
    
    if(nanowire->NDim<=2){

        inv_bound_pos    = new int*[2];
        inv_bound_pos[0] = new int[Boundary[0]->NA];

        NCH              = 0;
        NACH             = 0;
        NBOUND           = 0;

        boundary_dimension(&NBOUND,&NCH,nanowire->NDim,0);
        
        for(IA=Smin[1];IA<Smin[NSlab-1];IA++){
            if((Layer_Matrix[SLM*IA+2]>-100*tollim)&&(Layer_Matrix[SLM*IA+2]<(Lz-100*tollim))&&\
	       (!index_channel[IA])){
                index_channel[IA] = 0;
                ch_pos[NCH]       = IA;
                inv_ch_pos[IA]    = NCH;
                NCH++;
            }
            else{
	        factor = 0;
                if(Layer_Matrix[SLM*IA+2]<-100*tollim){
                    factor = -1;
                }
                if(Layer_Matrix[SLM*IA+2]>(Lz-100*tollim)){
                    factor = 1;
                }
                index_channel[IA] = index_channel[IA]+factor;
            }
        }

        NBOUND           = 0;
        inv_bound_pos[1] = new int[Boundary[1]->NA];

        boundary_dimension(&NBOUND,&NCH,nanowire->NDim,1);

        for(IA=0;IA<Around_Atom;IA++){
            if((Around_Matrix[3*IA+2]>-100*tollim)&&(Around_Matrix[3*IA+2]<(Lz-100*tollim))){
                index_arch[IA]   = 0;
                arch_pos[NACH]   = IA;
                arch_conv[IA]    = NACH;
                NACH++;
            }
            else{
                index_arch[IA]   = 1;
            }
        }

        Channel_tot = No_Atom;
        Around_tot  = Around_Atom;
        No_Atom     = NCH;
        Around_Atom = NACH;

	if(nanowire->NDim==2){
	    cell_area  = cell_area/3.0;
	}else{
	    cell_area  = cell_area/9.0;
	}
        
    }else{
        for(IB=0;IB<2;IB++){
            init_var(index_boundary[IB],Boundary[IB]->NA);
            for(IA=0;IA<Boundary[IB]->NA;IA++){
                bound_pos[IB][IA]  = IA;
                bound_conv[IB][IA] = IA;
            }
        }
        init_var(index_channel,No_Atom);
        init_var(index_arch,Around_Atom);
        for(IA=0;IA<No_Atom;IA++){
            ch_pos[IA]     = IA;
            ch_conv[IA]    = IA;
            inv_ch_pos[IA] = IA;
        }
        for(IA=0;IA<Around_Atom;IA++){
            arch_pos[IA]   = IA;
            arch_conv[IA]  = IA;
        }
        Channel_tot = No_Atom;
        Around_tot  = Around_Atom;
    }
    
}

/************************************************************************************************/

void WireGenerator::boundary_dimension(int *NBOUND,int *NCH,int NDim,int IB)
{

    int IA,factor,ind[2]={0,NSlab-1};
    
    if(NDim==1){
      
        for(IA=0;IA<Boundary[IB]->NA;IA++){
        
	    if((Boundary[IB]->UNN[SLM*IA+1]>-100*tollim)&&\
	       (Boundary[IB]->UNN[SLM*IA+1]<(Ly-100*tollim))){
		index_channel[Smin[ind[IB]]+IA] = 0;
	    }
	    else{
	        if(Boundary[IB]->UNN[SLM*IA+1]<-100*tollim){
		    factor = -3;
		}
		if(Boundary[IB]->UNN[SLM*IA+1]>(Ly-100*tollim)){
		    factor = 3;
		}
		index_channel[Smin[ind[IB]]+IA] = factor;
	    }
	}
    }
   
    for(IA=0;IA<Boundary[IB]->NA;IA++){
        if((Boundary[IB]->UNN[SLM*IA+2]>-100*tollim)&&	\
           (Boundary[IB]->UNN[SLM*IA+2]<(Lz-100*tollim))&&\
	   (!index_channel[Smin[ind[IB]]+IA])){
            index_boundary[IB][IA]          = 0;
            index_channel[Smin[ind[IB]]+IA] = 0;
            bound_pos[IB][*NBOUND]          = IA;
            ch_pos[*NCH]                     = Smin[ind[IB]]+IA;
            inv_bound_pos[IB][IA]           = *NBOUND;
            inv_ch_pos[Smin[ind[IB]]+IA]    = *NCH;
            *NBOUND                         = *NBOUND+1;
            *NCH                            = *NCH+1;
        }
        else{
	    factor = 0;
            if(Boundary[IB]->UNN[SLM*IA+2]<-100*tollim){
                factor = -1;
            }
            if(Boundary[IB]->UNN[SLM*IA+2]>(Lz-100*tollim)){
                factor = 1;
            }
            index_channel[Smin[ind[IB]]+IA] = index_channel[Smin[ind[IB]]+IA]+factor;
	    index_boundary[IB][IA]          = index_channel[Smin[ind[IB]]+IA];
         }
    }
    Boundary[IB]->NA        = *NBOUND;
    
    Boundary[IB]->cell_area = Boundary[IB]->NA*Vol_atom/(Boundary[IB]->cell_width*1e-27);

}

/************************************************************************************************/

void WireGenerator::convert_position(WireStructure *nanowire)
{
    int IS,IL,IA,IN;
    int n_element;
    double atom_pos[3],conv_pos[3];
    double shift[18]={-Ly,-Lz,-Ly,0.0,-Ly,Lz,0.0,-Lz,0.0,0.0,0.0,Lz,Ly,-Lz,Ly,0.0,Ly,Lz};
    
    if(nanowire->NDim<=2){

        for(IS=0;IS<NSlab;IS++){

            n_element = 0;
            
            for(IA=Smin[IS];IA<=Smax[IS];IA++){

	        c_dcopy(N3D,&Layer_Matrix[SLM*IA],1,atom_pos,1);
		c_daxpy(N2D,-1.0,&shift[2*(index_channel[IA]+4)],1,&atom_pos[1],1);
		
		for(IN=max(IA-NA_per_slice,Smin[IS]);IN<=min(IA+NA_per_slice,Smax[IS]);IN++){

		    c_dcopy(N3D,atom_pos,1,conv_pos,1);
		    c_daxpy(N3D,-1.0,&Layer_Matrix[SLM*IN],1,conv_pos,1);

		    if(c_dnrm2(N3D,conv_pos,1)<=1e-6){		      		    
		        ch_conv[IA] = inv_ch_pos[IN];
			if(IS==0){
			    bound_conv[0][IA-Smin[IS]] = inv_bound_pos[0][IN-Smin[IS]];
			}
			if(IS==NSlab-1){
			    bound_conv[1][IA-Smin[IS]] = inv_bound_pos[1][IN-Smin[IS]];
			}
			break;
		    }
		}

                if(!abs(index_channel[IA])) n_element++;
                
            }
            
            if(IS>0) Smin[IS] = Smax[IS-1]+1;
            Smax[IS] = Smin[IS]+n_element-1;
            
        }
	
        for(IL=0;IL<NLayer;IL++){

            n_element = 0;
            
            for(IA=Lmin[IL];IA<=Lmax[IL];IA++){

                if(!abs(index_channel[IA])) n_element++;
                
            }

            if(IL>0) Lmin[IL] = Lmax[IL-1]+1;
            Lmax[IL] = Lmin[IL]+n_element-1;
            
        }

        get_boundary_size(nanowire);
        
        for(IS=0;IS<2;IS++){
            delete[] inv_bound_pos[IS];
        }
        delete[] inv_bound_pos;
    }
}

/************************************************************************************************/

void WireGenerator::add_strain(Strain *strain,int no_strain_domain,StrainDomain **strain_domain,\
			       int update_at)
{

    if(strain->on||no_strain_domain){
      
        int IB,ID,IL;
	double strain_matrix[9];
	double xmin = 0.0;
	double xmax = INF;
	double ymin = 0.0;
	double ymax = INF;
	double zmin = 0.0;
	double zmax = INF;

	if(strain->on){
	    make_strained_unit_cell(strain);
	    no_strain_domain = 1;
	    c_dcopy(N3D*N3D,unit_cell->strain_matrix,1,strain_matrix,1);
	}

	for(ID=0;ID<no_strain_domain;ID++){

	    if(!update_at){

	        if(!strain->on){
		    get_strain_matrix(strain_matrix,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax,\
				      strain_domain[ID]);
		}

		strain_coordinate(Layer_Matrix,Channel_tot,SLM,strain_matrix,xmin,xmax,\
				  ymin,ymax,zmin,zmax,strain->on);

		for(IB=0;IB<2;IB++){
		    strain_coordinate(Boundary[IB]->UNN,Boundary[IB]->NA_tot,SLM,strain_matrix,\
				      xmin,xmax,ymin,ymax,zmin,zmax,strain->on);
		    strain_coordinate(Boundary[IB]->UNNp,Boundary[IB]->NA_tot,SLM,strain_matrix, \
				      xmin,xmax,ymin,ymax,zmin,zmax,strain->on);
		    strain_coordinate(Boundary[IB]->UNNm,Boundary[IB]->NA_tot,SLM,strain_matrix, \
				      xmin,xmax,ymin,ymax,zmin,zmax,strain->on);
		    Boundary[IB]->cell_width = unit_cell->strain_matrix[0]*Boundary[IB]->cell_width;
		    Boundary[IB]->cell_area  = unit_cell->strain_matrix[N3D+1]* \
		      unit_cell->strain_matrix[2*(N3D+1)]*Boundary[IB]->cell_area;
		}

		for(IL=0;IL<NLayer;IL++){
		    strain_coordinate(Enn[IL]->UNN,Enn[IL]->NA,SLM,strain_matrix,\
				      xmin,xmax,ymin,ymax,zmin,zmax,strain->on);
		    strain_coordinate(Enn[IL]->UNNp,Enn[IL]->NA,SLM,strain_matrix,\
				      xmin,xmax,ymin,ymax,zmin,zmax,strain->on);
		}

		cell_width = unit_cell->strain_matrix[0]*cell_width;
		cell_area  = unit_cell->strain_matrix[N3D+1]*unit_cell->strain_matrix[2*(N3D+1)]*cell_area;
		y_width    = unit_cell->strain_matrix[N3D+1]*y_width;
		Ly         = unit_cell->strain_matrix[N3D+1]*Ly;
		Lz         = unit_cell->strain_matrix[2*(N3D+1)]*Lz;
	    }

	    strain_coordinate(Around_Matrix,Around_tot,N3D,strain_matrix,xmin,xmax,ymin,\
			      ymax,zmin,zmax,strain->on);
	}
    }
}

/************************************************************************************************/

void WireGenerator::strain_coordinate(double *coord_matrix,int NR,int NC,double *strain_matrix,\
				      double xmin,double xmax,double ymin,double ymax,\
				      double zmin,double zmax,int strain_on)
{

    int IR;
    int cond1,cond2,cond3;
    int done;
    double coord[3];
    double shift[3];

    for(IR=0;IR<NR;IR++){

        shift[0] = xmin;
	shift[1] = ymin;
	shift[2] = zmin;
	done     = 0;

        cond1 = (coord_matrix[IR*NC]>=xmin-tollim)&&(coord_matrix[IR*NC]<=xmax+tollim);
	cond2 = (coord_matrix[IR*NC+1]>=ymin-tollim)&&(coord_matrix[IR*NC+1]<=ymax+tollim);
	cond3 = (coord_matrix[IR*NC+2]>=zmin-tollim)&&(coord_matrix[IR*NC+2]<=zmax+tollim);

	if((cond1&&cond2&&cond3&&(!done))||strain_on){

	    c_dcopy(N3D,&coord_matrix[IR*NC],1,coord,1);
	    c_daxpy(N3D,-1.0,shift,1,coord,1);

	    c_dgemv('N',N3D,N3D,1.0,strain_matrix,N3D,coord,1,0.0, \
		    &coord_matrix[IR*NC],1);

	    c_daxpy(N3D,1.0,shift,1,&coord_matrix[IR*NC],1);

	    if(NC>N3D){
	        c_daxpy(N3D,1.0,unit_cell->atom[Round(coord_matrix[IR*NC+N3D]-1)%NA]->int_disp,1,\
			&coord_matrix[IR*NC],1);
	    }

	    done = 1;
	}

	cond1 = (coord_matrix[IR*NC]>xmax);
	cond2 = (coord_matrix[IR*NC+1]>=ymin-tollim)&&(coord_matrix[IR*NC+1]<=ymax+tollim);
	cond3 = (coord_matrix[IR*NC+2]>=zmin-tollim)&&(coord_matrix[IR*NC+2]<=zmax+tollim);

	if(cond1&&cond2&&cond3&&(!done)&&(!strain_on)){

	    c_dcopy(N3D,&coord_matrix[IR*NC],1,coord,1);
	    c_daxpy(N3D,-1.0,shift,1,coord,1);

	    coord[0] = (xmax-xmin);
	    shift[0] = coord_matrix[IR*NC]-(xmax-xmin);

	    c_dgemv('N',N3D,N3D,1.0,strain_matrix,N3D,coord,1,0.0, \
		    &coord_matrix[IR*NC],1);

	    c_daxpy(N3D,1.0,shift,1,&coord_matrix[IR*NC],1);

	    done = 1;
	}

	cond1 = (coord_matrix[IR*NC]>=xmin-tollim)&&(coord_matrix[IR*NC]<=xmax+tollim);
	cond2 = (coord_matrix[IR*NC+1]>ymax);
	cond3 = (coord_matrix[IR*NC+2]>=zmin-tollim)&&(coord_matrix[IR*NC+2]<=zmax+tollim);

	if(cond1&&cond2&&cond3&&(!done)&&(!strain_on)){

	    c_dcopy(N3D,&coord_matrix[IR*NC],1,coord,1);
	    c_daxpy(N3D,-1.0,shift,1,coord,1);

	    coord[1] = (ymax-ymin);
	    shift[1] = coord_matrix[IR*NC+1]-(ymax-ymin);

	    c_dgemv('N',N3D,N3D,1.0,strain_matrix,N3D,coord,1,0.0, \
		    &coord_matrix[IR*NC],1);

	    c_daxpy(N3D,1.0,shift,1,&coord_matrix[IR*NC],1);

	    done = 1;
	}

	cond1 = (coord_matrix[IR*NC]>=xmin-tollim)&&(coord_matrix[IR*NC]<=xmax+tollim);
	cond2 = (coord_matrix[IR*NC+1]>=ymin-tollim)&&(coord_matrix[IR*NC+1]<=ymax+tollim);
	cond3 = (coord_matrix[IR*NC+2]>zmax+tollim);

	if(cond1&&cond2&&cond3&&(!done)&&(!strain_on)){

	    c_dcopy(N3D,&coord_matrix[IR*NC],1,coord,1);
	    c_daxpy(N3D,-1.0,shift,1,coord,1);

	    coord[2] = (zmax-zmin);
	    shift[2] = coord_matrix[IR*NC+2]-(zmax-zmin);

	    c_dgemv('N',N3D,N3D,1.0,strain_matrix,N3D,coord,1,0.0, \
		    &coord_matrix[IR*NC],1);

	    c_daxpy(N3D,1.0,shift,1,&coord_matrix[IR*NC],1);

	    done = 1;
	}
    }
}

/************************************************************************************************/

void WireGenerator::get_strain_matrix(double *strain_matrix,double *xmin,double *xmax,\
				      double *ymin,double *ymax,double *zmin,double *zmax,\
				      StrainDomain *strain_domain)
{

    init_var(strain_matrix,N3D*N3D);

    strain_matrix[0] = 1.0+strain_domain->Eps_vec[0];
    strain_matrix[1] = strain_domain->Eps_vec[3];
    strain_matrix[2] = strain_domain->Eps_vec[4];
    strain_matrix[3] = strain_domain->Eps_vec[3];
    strain_matrix[4] = 1.0+strain_domain->Eps_vec[1];
    strain_matrix[5] = strain_domain->Eps_vec[5];
    strain_matrix[6] = strain_domain->Eps_vec[4];
    strain_matrix[7] = strain_domain->Eps_vec[5];
    strain_matrix[8] = 1.0+strain_domain->Eps_vec[2];

    *xmin            = strain_domain->xmin;
    *xmax            = strain_domain->xmax;
    *ymin            = strain_domain->ymin;
    *ymax            = strain_domain->ymax;
    *zmin            = strain_domain->zmin;
    *zmax            = strain_domain->zmax;

}

/************************************************************************************************/

void WireGenerator::get_boundary_size(WireStructure *nanowire)
{
    int left_max=-INF,right_min,IL;
    int l_p_s = layer_per_slab/nanowire->NxFold;

    for(IL=0;IL<l_p_s;IL++){
        if(IL+neighbor_layer[IL]-l_p_s>left_max){
            left_max = IL+neighbor_layer[IL]-l_p_s;
        }
    }
    LB_size = Lmax[left_max]-Lmin[0]+1;
    
    for(IL=NLayer-2*l_p_s;IL<NLayer-l_p_s;IL++){
        if(IL+neighbor_layer[IL]>=NLayer-l_p_s){
            right_min = IL+l_p_s;
            break;
        }
    }
    RB_size = Lmax[NLayer-1]-Lmin[right_min]+1;

    if(nanowire->robust_numerics){
        LB_size = Lmax[layer_per_slab-1]-Lmin[0]+1;
	RB_size = Lmax[NLayer-1]-Lmin[NLayer-layer_per_slab]+1;
    }

}

/************************************************************************************************/

void WireGenerator::get_qm_region(WireStructure* nanowire)
{

    int IL;

    switch(nanowire->QMregion){
    case 0:
        QMfirst = 0;
	QMlast  = NLayer-1;
        break;
    case 2:
        QMfirst = 0;
	QMlast  = NLayer-1;
        for(IL=0;IL<NLayer;IL++){
	    if(Layer_Matrix[SLM*ch_pos[Lmin[IL]]]>=nanowire->QMstart){
	        QMfirst = IL;
		break;
	    }
        }
	for(IL=QMfirst;IL<NLayer;IL++){
	    if(Layer_Matrix[SLM*ch_pos[Lmin[IL]]]>=nanowire->QMstop){
	        QMlast = IL;
		break;
	    }
	}
        break;
    default:
      if(!mpi_rank){
	  printf("There is a problem with the definition of the Quantum Mechanical Region\n");
	  abort();
      }
    }

}

/************************************************************************************************/

void WireGenerator::adapt_parameters(Material *material)
{

    int IL,IA;

    int vec_size = max((transport_type==1)*No_Atom,1);
    atomic_mass  = new double[vec_size];

    if(transport_type==1){

        for(IL=0;IL<NLayer;IL++){
	    neighbor_layer[IL] = layer_per_slab;
	}

	orb_per_at[0] = 0;
	for(IA=0;IA<No_Atom;IA++){
	    orb_per_at[IA+1] = orb_per_at[IA]+N3D;
	    atomic_mass[IA]  = material->atomic_mass[Round(Layer_Matrix[SLM*ch_pos[IA]+3])-1];
	}

    }

}

/************************************************************************************************/

void WireGenerator::write_LM(const char *filename)
{
    int i,j;
    int local_rank;

    MPI_Comm_rank(MPI_COMM_WORLD,&local_rank);
    
    if(!local_rank){
        ofstream myfile;
	myfile.open(filename);
	myfile.precision(8);
	for(i=0;i<No_Atom;i++){
	    for(j=0;j<SLM;j++){
	        if(j<4){
		    myfile<<Layer_Matrix[SLM*ch_pos[i]+j]<<" ";
		}else{
		    if(Layer_Matrix[SLM*ch_pos[i]+j]>0){
		        myfile<<ch_conv[Round(Layer_Matrix[SLM*ch_pos[i]+j])-1]+1<<" ";
		    }else{
		        myfile<<Layer_Matrix[SLM*ch_pos[i]+j]<<" ";
		    }
		}
	    }
	    myfile<<"\n";
	}
	myfile.close();
    }

}

/************************************************************************************************/

void WireGenerator::write_TLM(const char *filename)
{
    int i,j;
    
    ofstream myfile;
    myfile.open(filename);
    myfile.precision(8);
    for(i=0;i<Channel_tot;i++){
        for(j=0;j<SLM;j++){
            myfile<<Layer_Matrix[SLM*i+j]<<" ";
        }
        myfile<<"\n";
    }
    myfile.close();

}

/************************************************************************************************/

void WireGenerator::write_AM(const char *filename)
{
    int i,j;
    
    ofstream myfile;
    myfile.open (filename);
    for(i=0;i<Around_Atom;i++){
        for(j=0;j<3;j++){
            myfile<<Around_Matrix[3*arch_pos[i]+j]<<" ";
        }
        myfile<<"\n";
    }
    myfile.close();
}

/************************************************************************************************/
                    
void WireGenerator::write_positions(const char *filename)
{                   
    int i, j;   
                   
    ofstream myfile;
    myfile.open(filename);
    myfile.precision(8);
    for(i=0;i<No_Atom;i++){
        for(j=0;j<4;j++){
            myfile<<Layer_Matrix[SLM*ch_pos[i]+j]<<" ";
        }
        myfile<<"\n";
    }
    myfile.close();
}

/************************************************************************************************/
