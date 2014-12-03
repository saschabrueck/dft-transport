#ifndef conversion_H
#define conversion_H

#include <set>

#include "CSR.H"
#include "cDBCSR.H"
#include "cDBCSR_tools.H"
#include "vector_tools.H"

#include "array_tools.H"

#define unused(x) ((void)(x))

/*===============================================================
edge_i_block: row pointer vector for every block
edge_i_blkr_tmp: sum over edge_i_block for every block row
=================================================================*/
template <class T>
void cDBCSR_to_CSR(c_DBCSR& dbcsr_mat, TCSR<T>* csr_mat)
{
   int m, n, mn, bc, block_p, fm_i, fm_j,
       nze_counter, nze_blk_row,
       old_row, new_row,
       N, n_nnz;

   T mat_element;
   T *block, *nnz_arr;
   int *edge_i_arr, *index_i_arr, *index_j_arr; 

   unused(fm_i);

   Vector1D<int> edge_i_blkr_tmp, edge_i_block, edge_i, index_j, index_i_vec;
   Vector1D<T> data_block, nnz;
   Vector2D<int> index_j_blkr;
   Vector2D<T> nnz_blkr;

   data_block.assign(dbcsr_mat.dblock, dbcsr_mat.dblock + dbcsr_mat.dblksize);
   edge_i.push_back(0);
   for(int br = 0; br < dbcsr_mat.nblkrows_total; br++){
      m = dbcsr_mat.row_blk_size[br];
      initialize_Vector2D (index_j_blkr, m);
      initialize_Vector2D (nnz_blkr, m);
      for(int blk = dbcsr_mat.row_p[br]+1; blk <= dbcsr_mat.row_p[br+1]; blk++){
         bc = dbcsr_mat.col_i[blk-1];
         n = dbcsr_mat.col_blk_size[bc-1];
         mn = m*n;
         block_p = abs(dbcsr_mat.blk_p[blk-1]);
         if(block_p != 0){
            block = new T[mn];
            for(int i=0; i<=mn-1; i++)
               block[i] = data_block[block_p+i-1];
            nze_counter = -1;
            old_row = 0;
            for(int r = 0; r < m; r++){
               for(int c = 0; c < n; c++){
//                  if(fabs(block[r+c*m]) >= 1e-60){
                  nze_counter = nze_counter + 1;
                  fm_i = r+dbcsr_mat.row_blk_offset[br]-1;
                  fm_j = c+dbcsr_mat.col_blk_offset[bc-1]-1;
                  mat_element = block[r+c*m];
                  new_row = fm_i+1;

                  index_j_blkr[r].push_back(fm_j);
                  nnz_blkr[r].push_back(mat_element);
                  if (new_row != old_row) edge_i_block.push_back(nze_counter); 
                  old_row = new_row;
//                  }; //if 
               }; //r row
            }; //c column
            edge_i_block.push_back(nze_counter+1);
            edge_i_blkr_tmp.resize(edge_i_block.size()); 
            std::transform (edge_i_blkr_tmp.begin(), edge_i_blkr_tmp.end(), edge_i_block.begin(), edge_i_blkr_tmp.begin(), std::plus<int>());
            edge_i_block.clear();
            delete [] block;
         }; //if
      }; //blk blocks
      nze_blk_row = edge_i.back();
      std::transform(edge_i_blkr_tmp.begin(), edge_i_blkr_tmp.end(), edge_i_blkr_tmp.begin(), bind2nd(std::plus<int>(), nze_blk_row));

      if (!edge_i_blkr_tmp.empty()) 
         edge_i.insert(edge_i.end(), edge_i_blkr_tmp.begin()+1, edge_i_blkr_tmp.end());

      convert_Vector2D_to_Vector1D(index_j_blkr, index_j);
      convert_Vector2D_to_Vector1D(nnz_blkr, nnz);
      index_j_blkr.clear();
      nnz_blkr.clear();
      edge_i_blkr_tmp.clear();
   }; //br block row

   N = edge_i.size()-1;
   n_nnz = nnz.size();
   edge_i_arr = new int[N+1];
   index_i_arr = new int[N];
   index_j_arr = new int[n_nnz];
   nnz_arr = new T[n_nnz];

   index_i_vec.resize(edge_i.size());
   adjacent_difference(edge_i.begin(), edge_i.end(), index_i_vec.begin());
   index_i_vec.erase(index_i_vec.begin());
   Vector1D_abs(index_i_vec);

   std::copy(edge_i.begin(),edge_i.end(),edge_i_arr);
   std::copy(index_i_vec.begin(),index_i_vec.end(),index_i_arr);
   std::copy(index_j.begin(),index_j.end(),index_j_arr);
   std::copy(nnz.begin(),nnz.end(),nnz_arr);

   csr_mat->set_nnz(n_nnz,nnz_arr);
   csr_mat->set_index_i(N,index_i_arr);
   csr_mat->set_index_j(n_nnz,index_j_arr);
   csr_mat->set_edge_i(N,edge_i_arr);

   delete [] edge_i_arr;
   delete [] index_i_arr;
   delete [] index_j_arr;
   delete [] nnz_arr;
};

/*===================================================================================================
outer for loop over the local block rows , 

col_ind (index_j) =0 4 5 1 4 5 2 4 5 3 4 5|4 5|5 ...
row_ptr (index_i) =0 3 6 9|12|14|15 18 21 24|27 ... 
rbs = cbs = 4 1 1 4 1 1 ...
******* NOTE: cbs can be different from rbs ********** 

cblkind : column block indices corresponding to global column indices
col_block_index : column block indices sorted according to the csr col_index   
cblkind_curr_rblk : column block indices for the current row block

index_i : 0 3 3 3 3 2 1 3 3 3 3 2 1 ...
rblk_p : pointer to the first element of each row block repeated by the number of rows of every row block
         0 0 0 0|4|5|6 6 6 6|10|11| ...
  
from index_i we get temp_vec1:
temp_vec1:          0 0 0 1 1 1 2 2 2 3 3 3|4 4|5|6 6 6 7 7 7 8 8 8 9 9 9|10 10|11|12 12 12 13 13 13 14 14 14 15 15 15
temp_vec2 (before): 0 0 0 0 0 0 0 0 0 0 0 0|4 4|5|6 6 6 6 6 6 6 6 6 6 6 6|10 10|11|12 12 12 12 12 12 12 12 12 12 12 12
temp_vec2 (after):  0 0 0 1 1 1 2 2 2 3 3 3|0 0|0|0 0 0 1 1 1 2 2 2 3 3 3|0  0 |0 |0  0  0  1  1  1  2  2  2  3  3  3              
from temp_vec2 and using order :
lrow_index:         0 1 2 3 0 1 2 3 0 1 2 3
                    0 0
                    0 
                    0 1 2 3 0 1 2 3 0 1 2 3 
                    0 0 
                    0

gcol_index  0 1 2 3 4 4 4 4 5 5 5 5 | 4 5 | 5 | 6 7 8 9 10 10 10 10 11 11 11 11 
col_diff    0 0 0 0 4 4 4 4 5 5 5 5 | 4 5 | 5 | 6 6 6 6 10 10 10 10 11 11 11 11
lcol_index  0 1 2 3 0 0 0 0 0 0 0 0 | 0 0 | 0 | 0 1 2 3 0  0  0  0  0  0  0  0 

local_nze : number of non-zero elements in every block blonging to the current row block
=====================================================================================================*/
template <class T>
void CSR_to_cDBCSR(TCSR<T>* csr_mat, c_DBCSR& dbcsr_mat,
                   Vector1D<int>& rbs, Vector1D<int>& cbs, 
                   Vector1D<int>& row_dist, Vector1D<int>& col_dist, 
                   Vector1D<int>& local_rows, Vector1D<int>& nblkrows_local_all)
{
   int            row_block_index, counter_old, counter_new, difference, ctr_old, ctr_new, sum,
                  count_blocks, column_difference, blocksize, blknze, rowblock_size,
                  row_p_head_size, row_p_tail_size, first_row, last_row, nblkrows_local,
                  rank;//, num_nodes, sendbuf;
   Vector1D<int>  edge_i, index_j, row_p_temp, row_p_head, row_p_tail,
                  row_p, col_i, blk_p, order,
                  gcol_index, col_i_temp, blk_p_temp,
                  lcol_index, col_diff, lrow_index, local_index, local_nze, 
                  cblkind, col_block_index, cblkind_curr_rblk, cblkind_curr_rblk_uniq, 
                  index_i, rblk_p, rbs_part_sum,
                  temp_vec1, temp_vec2,
                  rbo, cbo;//, nblkrows_local_all; 
   Vector1D<T>    dblock_temp, dblock, datablock;
   std::set<int>  temp_set ;
   Vector1D<int>::iterator iter, it1, it2, iter1, iter2;
   MPI::Intercomm Comm;

   Comm = MPI::COMM_WORLD;  
   rank = Comm.Get_rank();
   nblkrows_local = nblkrows_local_all[rank];
   row_p_temp.push_back(0);

   index_j.resize(csr_mat->n_nonzeros);
   std::copy(csr_mat->index_j , csr_mat->index_j+csr_mat->n_nonzeros, index_j.begin());
   index_i.resize(csr_mat->size);
   std::copy(csr_mat->index_i, csr_mat->index_i+csr_mat->size, index_i.begin());

   for (Vector1D<int>::iterator i = cbs.begin(); i != cbs.end(); i++) 
      cblkind.insert(cblkind.end(),*i,std::distance(cbs.begin(),i)+1);
   for (Vector1D<int>::iterator i = index_j.begin(); i != index_j.end(); i++) 
      col_block_index.push_back(cblkind[*i]);

   edge_i.resize(csr_mat->size+1);
   std::copy(csr_mat->edge_i, csr_mat->edge_i+csr_mat->size+1, edge_i.begin());

   sum = 0;

   if (local_rows.size() != 0){
      first_row = local_rows[0]-1;
      last_row = local_rows[nblkrows_local-1];
   }
   else{
      first_row = 0;
      last_row = 0;
   }

   for (Vector1D<int>::iterator i = rbs.begin()+first_row; i != rbs.begin()+last_row; i++){ 
      rblk_p.insert(rblk_p.end(),*i, sum);
      sum = sum + *i;
   }
   for (Vector1D<int>::iterator i = index_i.begin(); i != index_i.end(); i++){
      temp_vec1.insert(temp_vec1.end(),*i,std::distance(index_i.begin(),i));
      temp_vec2.insert(temp_vec2.end(),*i,rblk_p[std::distance(index_i.begin(),i)]);
   }
   std::transform(temp_vec1.begin(),temp_vec1.end(),temp_vec2.begin(),temp_vec2.begin(),std::minus<int>());

   counter_old = 0;
   ctr_old = 0;
   for (int index = 0 ; index < nblkrows_local; index ++ ) {
      row_block_index = local_rows[index]-1;
      rowblock_size = rbs[row_block_index];

      counter_new = counter_old + rowblock_size;
      difference = csr_mat->edge_i[counter_new] - csr_mat->edge_i[counter_old]; 
      ctr_new = ctr_old + difference; 

      cblkind_curr_rblk.insert(cblkind_curr_rblk.begin(),col_block_index.begin() + ctr_old,col_block_index.begin() + ctr_new);
      /* REMOVE DUPLICATE ELEMENTS OF cblkind_curr_rblk AND STORE THE UNIQUE ELEMENTS IN cblk_curr_rblk_uniq */
      cblkind_curr_rblk_uniq.resize(cblkind_curr_rblk.size());
      std::copy(cblkind_curr_rblk.begin(),cblkind_curr_rblk.end(),cblkind_curr_rblk_uniq.begin());
      for(iter1 = cblkind_curr_rblk_uniq.begin(), iter2 = cblkind_curr_rblk_uniq.begin(); iter1 != cblkind_curr_rblk_uniq.end(); ++iter1)
         if(temp_set.insert(*iter1).second) *iter2++ = *iter1;
      cblkind_curr_rblk_uniq.erase(iter2, cblkind_curr_rblk_uniq.end());
      /* FIND THE INDICES WHERE THE UNIQUE ELEMENTS ARE LOCATED IN cblkind_curr_blk, THE RESULT PLUS
         THE FIRST ELEMENT OF THE row_p FOR THE CURRENT ROW BLOCK WOULD BE THE VECTOR ORDER */ 
      for(it2 = cblkind_curr_rblk_uniq.begin(); it2 != cblkind_curr_rblk_uniq.end(); ++it2){
         blknze = 0;
         for(it1 = cblkind_curr_rblk.begin(); it1 != cblkind_curr_rblk.end(); ++it1){
            if(*it1 == *it2) {
               order.push_back(csr_mat->edge_i[counter_old] + std::distance(cblkind_curr_rblk.begin(),it1)); 
               ++ blknze;
            }
        }
        local_nze.push_back(blknze); 
      }  

      count_blocks = cblkind_curr_rblk_uniq.size();
      row_p_temp.push_back(row_p_temp.back()+count_blocks);

      for (Vector1D<int>::iterator j = order.begin() ; j != order.end(); j ++) {
         gcol_index.push_back(index_j[*j]);
         lrow_index.push_back(temp_vec2[*j]);
         dblock_temp.push_back(csr_mat->nnz[*j]);
         col_i_temp.push_back(cblkind[gcol_index.back()]);
         column_difference = std::accumulate(cbs.begin(),cbs.begin()+col_i_temp.back()-1,0);
         col_diff.push_back(column_difference);
      }

      lcol_index.resize(col_diff.size());
      std::transform(gcol_index.begin(),gcol_index.end(),col_diff.begin(),lcol_index.begin(),std::minus<int>());

      for (unsigned i = 0; i < lrow_index.size(); ++i)
         local_index.push_back(lrow_index[i]*cbs[col_i_temp[i]-1]+lcol_index[i]);

      iter = std::unique(col_i_temp.begin(), col_i_temp.end());
      col_i_temp.resize(std::distance(col_i_temp.begin(),iter));
      col_i.insert(col_i.end(),col_i_temp.begin(),col_i_temp.end());

      for (unsigned i = 0; i < local_nze.size(); ++i) {
         blocksize = rbs[row_block_index]*cbs[col_i_temp[i]-1];
         datablock.assign(blocksize,0);
         for (int ind = 0; ind < local_nze[i]; ind++)
            datablock[local_index[ind]] = dblock_temp[ind];
         
         dblock_temp.erase(dblock_temp.begin(),dblock_temp.begin()+local_nze[i]);
         local_index.erase(local_index.begin(),local_index.begin()+local_nze[i]);
         dblock.insert(dblock.end(),datablock.begin(),datablock.end());
      }

      for (Vector1D<int>::iterator k = col_i_temp.begin(); k != col_i_temp.end(); k++) 
         blk_p_temp.push_back(rbs[*k-1]);
      /* NOTE: row_block_index is used for cbs */ 
      std::transform(blk_p_temp.begin(), blk_p_temp.end(), blk_p_temp.begin(), bind2nd(std::multiplies<int>(), cbs[row_block_index])); 
      blk_p.insert(blk_p.end(),blk_p_temp.begin(),blk_p_temp.end());

      counter_old = counter_new; 
      ctr_old = ctr_new;

      cblkind_curr_rblk.clear();
      cblkind_curr_rblk_uniq.clear();
      temp_set.clear();
      order.clear();
      lrow_index.clear();
      local_nze.clear();
      gcol_index.clear();
      col_diff.clear();
      local_index.clear();
      col_i_temp.clear();
      lcol_index.clear();
      blk_p_temp.clear();
      dblock_temp.clear();
      datablock.clear();

   } // for row_block_index

   blk_p.insert(blk_p.begin(),1);
   partial_sum(blk_p.begin(), blk_p.end(), blk_p.begin());
   blk_p.pop_back();

   row_p_tail_size = std::accumulate(nblkrows_local_all.begin()+rank+1,nblkrows_local_all.end(),0);
   row_p_head_size = std::accumulate(nblkrows_local_all.begin(),nblkrows_local_all.begin()+rank,0);

   row_p_head.assign(row_p_head_size,row_p_temp.front());
   row_p_tail.assign(row_p_tail_size,row_p_temp.back());
   row_p.insert(row_p.end(),row_p_head.begin(),row_p_head.end());
   row_p.insert(row_p.end(),row_p_temp.begin(),row_p_temp.end());
   row_p.insert(row_p.end(),row_p_tail.begin(),row_p_tail.end());

   rbo = rbs;
   cbo = cbs;
   rbo.insert(rbo.begin(),1);
   cbo.insert(cbo.begin(),1);
   std::partial_sum(rbo.begin(), rbo.end(), rbo.begin(), plus<int>());
   std::partial_sum(cbo.begin(), cbo.end(), cbo.begin(), plus<int>());

   int data_type = 3; /* dbcsr_types : dbcsr_type_real_8 = 3 */
   c_DBCSR_create (dbcsr_mat, rank, nblkrows_local, data_type, rbs, cbs, rbo, cbo, 
                   row_p, col_i, blk_p, local_rows, dblock, row_dist, col_dist);
}; 

#endif
