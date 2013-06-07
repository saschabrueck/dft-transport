#include "cDBCSR_tools.H"

//template <class T>
void c_DBCSR_init (c_DBCSR& c_dbcsr_mat, int nblks, int n_nze, int data_area_size,
                   int data_type, int nblkrows_local, int nblkrows_total, int nblkcols_total, 
                   int fullmatrix_nrows, int fullmatrix_ncols, int myid)
{
   c_dbcsr_mat.row_blk_size = new int[nblkrows_total];
   c_dbcsr_mat.col_blk_size = new int[nblkcols_total];
   c_dbcsr_mat.row_blk_offset = new int[nblkrows_total+1];
   c_dbcsr_mat.col_blk_offset = new int[nblkcols_total+1];
   c_dbcsr_mat.row_p = new int[nblkrows_total+1];
   c_dbcsr_mat.col_i = new int[nblks];
   c_dbcsr_mat.blk_p = new int[nblks];
   c_dbcsr_mat.local_rows = new int[nblkrows_local];
   c_dbcsr_mat.dblock = new double[data_area_size];
   c_dbcsr_mat.row_dist = new int[nblkrows_total];
   c_dbcsr_mat.col_dist = new int[nblkcols_total];

   c_dbcsr_mat.nblks = nblks;
   c_dbcsr_mat.n_nze = n_nze;
   c_dbcsr_mat.dblksize = data_area_size;
   c_dbcsr_mat.data_type = data_type;
   c_dbcsr_mat.nblkrows_local = nblkrows_local;
   c_dbcsr_mat.nblkrows_total = nblkrows_total;
   c_dbcsr_mat.nblkcols_total = nblkcols_total;
   c_dbcsr_mat.fullmatrix_nrows = fullmatrix_nrows;
   c_dbcsr_mat.fullmatrix_ncols = fullmatrix_ncols;
   c_dbcsr_mat.myid = myid;
}

//template <class T>
void c_DBCSR_create (c_DBCSR& c_dbcsr_mat, int myid, int nblkrows_local, int data_type, 
                     Vector1D<int>& rbs, Vector1D<int>& cbs, Vector1D<int>& rbo, Vector1D<int>& cbo,
                     Vector1D<int>& row_p, Vector1D<int>& col_i, Vector1D<int>& blk_p,
                     Vector1D<int>& local_rows, Vector1D<double>& dblock, 
                     Vector1D<int>& row_dist, Vector1D<int>& col_dist)
{
   int nblks, n_nze, data_area_size, nblkrows_total, nblkcols_total, 
       fullmatrix_nrows, fullmatrix_ncols;

   nblks = col_i.size(); 
   n_nze = dblock.size();
   data_area_size = dblock.size();
   nblkrows_total = rbs.size();
   nblkcols_total = cbs.size();
   fullmatrix_nrows = std::accumulate(rbs.begin(),rbs.end(),0);
   fullmatrix_ncols = std::accumulate(cbs.begin(),cbs.end(),0);

   c_dbcsr_mat.row_blk_size = new int[nblkrows_total];
   c_dbcsr_mat.col_blk_size = new int[nblkcols_total];
   c_dbcsr_mat.row_blk_offset = new int[nblkrows_total+1];
   c_dbcsr_mat.col_blk_offset = new int[nblkcols_total+1];
   c_dbcsr_mat.row_p = new int[nblkrows_total+1];
   c_dbcsr_mat.col_i = new int[nblks];
   c_dbcsr_mat.blk_p = new int[nblks];
   c_dbcsr_mat.local_rows = new int[nblkrows_local];
   c_dbcsr_mat.dblock = new double[data_area_size];
   c_dbcsr_mat.row_dist = new int[nblkrows_total];
   c_dbcsr_mat.col_dist = new int[nblkcols_total];

   std::copy(rbs.begin(),rbs.end(),c_dbcsr_mat.row_blk_size);
   std::copy(cbs.begin(),cbs.end(),c_dbcsr_mat.col_blk_size);
   std::copy(rbo.begin(),rbo.end(),c_dbcsr_mat.row_blk_offset);
   std::copy(cbo.begin(),cbo.end(),c_dbcsr_mat.col_blk_offset);
   std::copy(row_p.begin(),row_p.end(),c_dbcsr_mat.row_p);
   std::copy(col_i.begin(),col_i.end(),c_dbcsr_mat.col_i);
   std::copy(blk_p.begin(),blk_p.end(),c_dbcsr_mat.blk_p);
   std::copy(local_rows.begin(),local_rows.end(),c_dbcsr_mat.local_rows);
   std::copy(dblock.begin(),dblock.end(),c_dbcsr_mat.dblock);
   std::copy(row_dist.begin(),row_dist.end(),c_dbcsr_mat.row_dist);
   std::copy(col_dist.begin(),col_dist.end(),c_dbcsr_mat.col_dist);

   c_dbcsr_mat.nblks = nblks;
   c_dbcsr_mat.n_nze = n_nze;
   c_dbcsr_mat.dblksize = data_area_size;
   c_dbcsr_mat.data_type = data_type;
   c_dbcsr_mat.nblkrows_local = nblkrows_local;
   c_dbcsr_mat.nblkrows_total = nblkrows_total;
   c_dbcsr_mat.nblkcols_total = nblkcols_total;
   c_dbcsr_mat.fullmatrix_nrows = fullmatrix_nrows;
   c_dbcsr_mat.fullmatrix_nrows = fullmatrix_ncols;
   c_dbcsr_mat.myid = myid;
}

void c_DBCSR_delete (c_DBCSR& c_dbcsr_mat)
{
   delete[] c_dbcsr_mat.row_blk_size;
   delete[] c_dbcsr_mat.col_blk_size;
   delete[] c_dbcsr_mat.row_blk_offset;
   delete[] c_dbcsr_mat.col_blk_offset;
   delete[] c_dbcsr_mat.row_p;
   delete[] c_dbcsr_mat.col_i;
   delete[] c_dbcsr_mat.blk_p;
   delete[] c_dbcsr_mat.local_rows;
   delete[] c_dbcsr_mat.dblock;
   delete[] c_dbcsr_mat.row_dist;
   delete[] c_dbcsr_mat.col_dist;

   c_dbcsr_mat.nblks = 0;
   c_dbcsr_mat.n_nze = 0;
   c_dbcsr_mat.dblksize = 0;
   c_dbcsr_mat.data_type = 0;
   c_dbcsr_mat.nblkrows_local = 0;
   c_dbcsr_mat.nblkrows_total = 0;
   c_dbcsr_mat.nblkcols_total = 0;
   c_dbcsr_mat.fullmatrix_nrows = 0;
   c_dbcsr_mat.fullmatrix_ncols = 0;
   c_dbcsr_mat.myid = 0;
}

