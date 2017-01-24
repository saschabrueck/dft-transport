/***************************************************************************** 
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This application is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)

Written by:  Chris Bowen
             Gerhard Klimeck
             Fabiano Oyafuso
             Seungwon Lee
             Olga Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/base/MPI_Timing.h,v 1.13 2007/08/10 15:03:32 baeh Exp $ 
*****************************************************************************/

#ifndef MPI_TIMING_H
#define MPI_TIMING_H

#ifdef MPI_TIMING
#include <mpi.h>

int mpi_omen_id;

class MPI_Timing{
    public:
	    // TODO: Add any timing variables to store
        double total_0, total;
        double setup_0, setup;
        double work_0, work;
        double fileio_0, fileio;

        double matmul_spds;
        double Hv_tot;
        double h3ddiag_small;
        double h3doffd_small;

        // MPI profile for strain caculation
		// sendrecv, allreduce, barrier, bcast, send, and recv
		double strain_mpi[6];

        // MPI profile for electronic structure
        double electr_mpi[6];

		// MPI profile for rest of the code
		double other_mpi[6];

		// MPI profile location
		double *section;

        double start_time;
        double start_h3ddiag;
        double start_h3doffd;

		// Lanczos /////////////
		double lanczos_single_iter;
		double lanczos_eigenvalue;
		double lanczos_eigenvector;
		double lanczos_total;
		////////////////////////

        double matmul_spds_small;
        double krylov;
        double h3ddiag;
        double h3doffd;
        double strain;

        double start_time5;
        double test5;
        double start_time6;
        double test6;
        double start_time7;
        double test7;
        double start_time8;
        double test8;
        double start_time9;
        double test9;

        double parpack_function;
        double rayleigh_function;
  
        // This is for matrix matrix multiply sendrecv routine
        double matmult_sendrecv;

        int out_data_length;

        void Initialize() {
            total_0 = total = 0.0;
            setup_0 = setup = 0.0;
            work_0 = work = 0.0;
            
            matmul_spds=0.0;
            h3ddiag_small=0.0;
            h3doffd_small=0.0;

            start_time=0.0;
            start_h3ddiag=0.0;
            start_h3doffd=0.0;

            fileio_0 = fileio =0.0;

		    lanczos_single_iter=0.0;
		    lanczos_eigenvalue=0.0;
		    lanczos_eigenvector=0.0;
		    lanczos_total=0.0;

            matmul_spds_small=0.0;
            Hv_tot=0.0;
            krylov=0.0;
            h3ddiag=0.0;
            h3doffd=0.0;
            strain =0.0;
      
            start_time5=0.0;
            test5=0.0;
            start_time6=0.0;
            test6=0.0;
            start_time7=0.0;
            test7=0.0;
            start_time8=0.0;
            test8=0.0;
            start_time9=0.0;
            test9=0.0;
   
            parpack_function=0.0;
            rayleigh_function=0.0;
    
            matmult_sendrecv = 0.0;

            for (int i=0; i<6; i++) {
              strain_mpi[i] = 0.0;
              electr_mpi[i] = 0.0;
              other_mpi[i] = 0.0;
            }
            section = other_mpi;
            out_data_length=25;
        }

        MPI_Timing() { Initialize(); }

		#if 1 // TODO: Construct output format on your own.
        void print() {
            printf("------------------  CPU %2d ------------------\n", mpi_omen_id);
			printf("cpu(%3d) LANCZOS: total=%5.1f iter=%5.1f eigenvalue=%5.1f eigenvector=%5.1f\n",mpi_omen_id,lanczos_total,lanczos_single_iter,lanczos_eigenvalue,lanczos_eigenvector);
        }
		#else
        void print() {
            printf("------------------  CPU %2d ------------------\n", mpi_omen_id);
            printf("cpu: %3d TOTAL:        total=%5.1f\t setup=%5.1f\t strain=%5.1f electronics=%5.1f fileio=%5.1f\n", 
                           mpi_omen_id,total, setup, strain, work, fileio);
            printf("cpu: %3d MATRIX:       h3ddia=%5.1f\t h3doff=%5.1f\t multiply_nocom=%5.1f\t  multiply_tot=%5.1f\n",
                           mpi_omen_id,h3ddiag, h3doffd, matmul_spds, Hv_tot);
            printf("cpu: %3d COMM1:         sendrec=%5.1f\t allred=%5.1f\t barrier=%5.1f\t bcast=%5.1f\t send=%5.1f\t receive=%5.1f\n",
                           mpi_omen_id, 
              other_mpi[0]+strain_mpi[0]+electr_mpi[0],
              other_mpi[1]+strain_mpi[1]+electr_mpi[1],
              other_mpi[2]+strain_mpi[2]+electr_mpi[2],
              other_mpi[3]+strain_mpi[3]+electr_mpi[3],
              other_mpi[4]+strain_mpi[4]+electr_mpi[4],
              other_mpi[5]+strain_mpi[5]+electr_mpi[5]);
            printf("cpu: %3d COMM2:      matmult_sendrecv = %5.1f\n", mpi_omen_id,matmult_sendrecv);
            printf("cpu: %3d STRAIN:     %5.1f\n", mpi_omen_id,strain);
            printf("cpu: %3d LANCZ:      %5.1f\n", mpi_omen_id,lanczos_function); 
            printf("cpu: %3d RAYLEIGH:   %5.1f\n", mpi_omen_id,rayleigh_function); 
            printf("cpu: %3d PARPACK:    %5.1f\n", mpi_omen_id,parpack_function); 
        }
		#endif

        //replace this after complete porting done from NEMO
        //rvectr GetData() {
        double *GetData() {
            //rvectr out_data = Rvectr(out_data_length);
            double *out_data = new double [out_data_length];
            out_data[0]  = (double) mpi_omen_id;
            out_data[1]  = total;
            out_data[2]  = setup;
            out_data[3]  = strain;
            out_data[4]  = work;
            out_data[5]  = fileio;
            out_data[12] = matmult_sendrecv;
            for(int i=0; i<6; i++) {
                out_data[6+i]=other_mpi[i];
                out_data[13+i]=strain_mpi[i];
                out_data[19+i]=electr_mpi[i];
            }	
            return out_data;
        }
}; /* class MPI_Timing */

// Global variable
// Just include MPI_Timing.H to each file that you want to use TIC/TOC
MPI_Timing mpiTiming;
///////////////////////////////////////////////////////////////////////

// Wrapper macros for mpi profiling
#define MPI_Sendrecv(sa,sb,sc,sd,se,ra,rb,rc,rd,re,f,g) { \
	double t = MPI_Wtime(); \
	MPI_Sendrecv(sa,sb,sc,sd,se,ra,rb,rc,rd,re,f,g); \
	mpiTiming.section[0] += MPI_Wtime()-t;}

#define MPI_Allreduce(a,b,c,d,e,f) { \
	double t = MPI_Wtime(); \
	MPI_Allreduce(a,b,c,d,e,f); \
	mpiTiming.section[1] += MPI_Wtime()-t;}

#define MPI_Barrier(a) { \
	double t = MPI_Wtime(); \
	MPI_Barrier(a); \
	mpiTiming.section[2] += MPI_Wtime()-t;} 

#define MPI_Bcast(a,b,c,d,e) { \
	double t = MPI_Wtime(); \
	MPI_Bcast(a,b,c,d,e); \
	mpiTiming.section[3] += MPI_Wtime()-t;}

#define MPI_Send(a,b,c,d,e,f) { \
	double t = MPI_Wtime(); \
	MPI_Send(a,b,c,d,e,f); \
	mpiTiming.section[4] += MPI_Wtime()-t;}

#define MPI_Recv(a,b,c,d,e,f,g) { \
	double t = MPI_Wtime(); \
	MPI_Recv(a,b,c,d,e,f,g); \
	mpiTiming.section[5] += MPI_Wtime()-t;}

// NOTE: TIC/TOC usage example
// double mpi_time;
// MPI_TIC(mpi_time);
// ....
// MPI_TOC("member variable from class MPI_Timing",mpi_time);
#define MPI_TIME_INIT(T) double T=0.0
#define MPI_TIC(T)       T = MPI_Wtime()
#define MPI_TOC(T,T0)    T += (MPI_Wtime() - T0)

#endif /* MPI_TIMING */
#endif /* MPI_TIMING_H */
