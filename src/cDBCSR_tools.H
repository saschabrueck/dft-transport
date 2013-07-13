#ifndef cDBCSR_tools_H
#define cDBCSR_tools_H

#include <numeric>
#include "cDBCSR.H"
#include "vectorTypes.H"

//template <class T>
void c_DBCSR_init (c_DBCSR& c_dbcsr_mat, int nblks, int n_nze, int data_area_size,
                   int data_type, int nblkrows_local, int nblkrows_total, int nblkcols_total,
                   int fullmatrix_nrows, int fullmatrix_ncols, int myid);

//template <class T>
void c_DBCSR_create (c_DBCSR& c_dbcsr_mat, int myid, int nblkrows_local, int data_type,
                     Vector1D<int>& rbs, Vector1D<int>& cbs, Vector1D<int>& rbo, Vector1D<int>& cbo,
                     Vector1D<int>& row_p, Vector1D<int>& col_i, Vector1D<int>& blk_p,
                     Vector1D<int>& local_rows, Vector1D<double>& dblock,
                     Vector1D<int>& row_dist, Vector1D<int>& col_dist);

void c_DBCSR_delete (c_DBCSR& c_dbcsr_mat);

#endif
