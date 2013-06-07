#ifndef cDBCSR_H
#define cDBCSR_H

extern "C"
{
   typedef struct{ 
           int nblks, n_nze, dblksize, data_type,
               nblkrows_local, nblkrows_total, nblkcols_total, 
               fullmatrix_nrows, fullmatrix_ncols, myid;
           int *row_blk_size, *col_blk_size, *row_blk_offset, *col_blk_offset,
               *row_p, *col_i, *blk_p, *local_rows, *row_dist, *col_dist;
           double *dblock; 
           } c_DBCSR;
}

#endif
