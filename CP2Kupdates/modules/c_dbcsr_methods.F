!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief 
!>        
!> \par History
!>       2013.05 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! *****************************************************************************
MODULE c_dbcsr_methods

  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_new
  USE c_dbcsr_types,                   ONLY: c_DBCSR
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_create,&
                                             cp_dbcsr_get_info,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_release
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE dbcsr_data_methods,              ONLY: dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_get_data
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: default_string_length
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_col_dist,&
                                             dbcsr_distribution_init,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_row_dist,&
                                             dbcsr_get_nze,&
                                             dbcsr_mp_new
  USE dbcsr_transformations,           ONLY: dbcsr_fill,&
                                             dbcsr_redistribute
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_mp_obj
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

  USE, INTRINSIC :: ISO_C_BINDING

#include "cp_common_uses.h"

  IMPLICIT NONE

  PUBLIC :: cDBCSR_new, cDBCSR_to_dbcsr, cDBCSR_matrix_get_info

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'c_dbcsr_methods'

CONTAINS

SUBROUTINE cDBCSR_new (dbcsr_mat, c_DBCSR_mat, error)

  USE, INTRINSIC :: ISO_C_BINDING,        ONLY: C_ASSOCIATED, C_LOC, C_PTR, C_DOUBLE, C_INT
  IMPLICIT NONE

  TYPE(c_DBCSR), INTENT(INOUT)            :: c_DBCSR_mat
  TYPE(cp_dbcsr_type), INTENT(IN)         :: dbcsr_mat
  TYPE(cp_error_type), INTENT(INOUT)      :: error
  
  INTEGER                                 :: nblkrows_local, nblkrows_total, nblkcols_total,&
                                             nze, data_type, nblks, handle
  TYPE(cp_logger_type), POINTER           :: logger
  TYPE(array_i1d_obj)                     :: rbs_array, cbs_array, rbo_array, cbo_array, r_dist_array, c_dist_array
  CHARACTER                               :: matrix_type
  CHARACTER(LEN=default_string_length)    :: matrix_name
  TYPE(dbcsr_distribution_obj)            :: distribution
  REAL(C_DOUBLE), DIMENSION(:), POINTER   :: c_blkdata

  INTEGER(C_INT), DIMENSION(:), POINTER   :: col_blk_offset, col_blk_size,&
                                             row_blk_offset, row_blk_size,&
                                             row_p, col_i, blk_p, &
                                             local_rows, row_dist, col_dist
  REAL(dp), DIMENSION(:), POINTER         :: r_dp
  TYPE(dbcsr_data_obj)                    :: data_area
  TYPE(C_PTR)                             :: ptr_to_rowp, ptr_to_coli, ptr_to_blkp,&
                                             ptr_to_rowblksize, ptr_to_colblksize,&
                                             ptr_to_rowblkoffset, ptr_to_colblkoffset,&
                                             ptr_to_cblkdata, ptr_to_localrows,&
                                             ptr_to_row_dist, ptr_to_col_dist

  CHARACTER(LEN=*), PARAMETER :: routineN = 'cDBCSR_new', &
      routineP = moduleN//':'//routineN
!INTEGER :: n_nze
!INTEGER :: br, m , blk, bc, n, mn, block_p, c, r 
!REAL(dp), DIMENSION(:), POINTER      :: matrix

  CALL timeset(routineN,handle)
  logger => cp_error_get_logger(error)

  CALL cp_dbcsr_get_info(dbcsr_mat, nblkrows_total=nblkrows_total, nblkcols_total=nblkcols_total,&
                         nblkrows_local=nblkrows_local, local_rows=local_rows,& 
                         row_blk_size=rbs_array, col_blk_size=cbs_array,&
                         row_blk_offset=rbo_array, col_blk_offset=cbo_array, distribution = distribution, & 
                         name=matrix_name, data_area=data_area, matrix_type=matrix_type, data_type=data_type)
  nze = dbcsr_get_nze(dbcsr_mat%matrix)
  nblks = dbcsr_mat%matrix%m%nblks
  c_DBCSR_mat%n_nze = nze
  c_DBCSR_mat%nblks = nblks
  c_DBCSR_mat%data_type = data_type
  c_DBCSR_mat%nblkrows_local = nblkrows_local
  c_DBCSR_mat%nblkrows_total = nblkrows_total
  c_DBCSR_mat%nblkcols_total = nblkcols_total
  c_DBCSR_mat%myid = logger%para_env%mepos

  r_dist_array = dbcsr_distribution_row_dist(distribution)
  c_dist_array = dbcsr_distribution_col_dist(distribution)
  row_dist => array_data (r_dist_array)
  col_dist => array_data (c_dist_array)

  CALL dbcsr_get_data (data_area, r_dp)
  c_blkdata => r_dp
  row_blk_size => array_data (rbs_array)
  col_blk_size => array_data (cbs_array)
  row_blk_offset => array_data (rbo_array)
  col_blk_offset => array_data (cbo_array)
  row_p => dbcsr_mat%matrix%m%row_p
  col_i => dbcsr_mat%matrix%m%col_i
  blk_p => dbcsr_mat%matrix%m%blk_p
  c_DBCSR_mat%dblksize = SIZE(r_dp)
  c_DBCSR_mat%fullmatrix_nrows = SUM(row_blk_size)
  c_DBCSR_mat%fullmatrix_ncols = SUM(col_blk_size)

  IF (.NOT. C_ASSOCIATED(ptr_to_rowp)) ptr_to_rowp = C_LOC(row_p(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_coli)) ptr_to_coli = C_LOC(col_i(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_blkp)) ptr_to_blkp = C_LOC(blk_p(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_localrows)) ptr_to_localrows = C_LOC(local_rows(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_rowblksize)) ptr_to_rowblksize = C_LOC(row_blk_size(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_colblksize)) ptr_to_colblksize = C_LOC(col_blk_size(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_rowblkoffset)) ptr_to_rowblkoffset = C_LOC(row_blk_offset(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_colblkoffset)) ptr_to_colblkoffset = C_LOC(col_blk_offset(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_cblkdata)) ptr_to_cblkdata = C_LOC(c_blkdata(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_row_dist)) ptr_to_row_dist = C_LOC(row_dist(1))
  IF (.NOT. C_ASSOCIATED(ptr_to_col_dist)) ptr_to_col_dist = C_LOC(col_dist(1))

  c_DBCSR_mat%row_p = ptr_to_rowp
  c_DBCSR_mat%col_i = ptr_to_coli
  c_DBCSR_mat%blk_p = ptr_to_blkp
  c_DBCSR_mat%local_rows = ptr_to_localrows
  c_DBCSR_mat%row_blk_size = ptr_to_rowblksize
  c_DBCSR_mat%col_blk_size = ptr_to_colblksize
  c_DBCSR_mat%row_blk_offset = ptr_to_rowblkoffset
  c_DBCSR_mat%col_blk_offset = ptr_to_colblkoffset
  c_DBCSR_mat%dblock = ptr_to_cblkdata
  c_DBCSR_mat%row_dist = ptr_to_row_dist
  c_DBCSR_mat%col_dist = ptr_to_col_dist

!n_nze = 0
!DO br = 1, nblkrows_total
!   m = row_blk_size(br)
!   DO blk = row_p(br)+1, row_p(br+1)
!      bc = col_i(blk)
!      n = col_blk_size(bc)
!      mn = m*n
!      block_p = ABS(blk_p(blk))
!      block_exists: IF (block_p .ne. 0) THEN
!         matrix => r_dp(block_p:block_p+mn-1)
!         DO c = 1,n
!            DO r = 1,m
!               IF(ABS(matrix(r+(c-1)*m)) .gt. 1e-18_real_8) THEN 
!!                 write(*,*) r+row_blk_offset(br)-1, c+col_blk_offset(bc)-1 , matrix(r+(c-1)*m)
!                 n_nze = n_nze + 1
!               END IF 
!            END DO
!         END DO
!      END IF block_exists
!  END DO
!END DO
!write(*,*) "here ... " , n_nze

  CALL timestop(handle)
END SUBROUTINE cDBCSR_new 

SUBROUTINE cDBCSR_matrix_get_info (c_DBCSR_mat, nblks, n_nze, data_block_size, data_type,&
                                   nblkrows_local, nblkrows_total, nblkcols_total,&
                                   fullmatrix_nrows, fullmatrix_ncols, mpi_id,&
                                   row_blk_size, col_blk_size, row_blk_offset, col_blk_offset,& 
                                   row_p, col_i, blk_p, local_rows, row_dist, col_dist, r_dp, error)

  USE, INTRINSIC :: ISO_C_BINDING,        ONLY: C_F_POINTER
  IMPLICIT NONE

  TYPE(c_DBCSR), INTENT(IN)                    :: c_DBCSR_mat
  TYPE(cp_error_type), INTENT(INOUT)           :: error
  INTEGER, INTENT(OUT)                         :: n_nze, nblks, data_block_size, nblkrows_local, nblkrows_total, nblkcols_total,&
                                                  fullmatrix_nrows, fullmatrix_ncols, mpi_id, data_type
  INTEGER, DIMENSION(:), POINTER, INTENT(OUT)  :: row_p, col_i, blk_p
  TYPE(array_i1d_obj), INTENT(OUT)             :: row_blk_size, col_blk_size, row_blk_offset, col_blk_offset, &
                                                  local_rows, row_dist, col_dist
  REAL(dp), DIMENSION(:), POINTER, INTENT(OUT) :: r_dp

  INTEGER, DIMENSION(:), POINTER     :: rbs, cbs, rbo, cbo, l_rows, r_dist, c_dist
  TYPE(cp_logger_type), POINTER      :: logger
  INTEGER                            :: handle

  CHARACTER(LEN=*), PARAMETER :: routineN = 'cDBCSR_matrix_get_info', &
      routineP = moduleN//':'//routineN
  CALL timeset(routineN,handle)

  logger => cp_error_get_logger(error)

  nblks = c_DBCSR_mat%nblks 
  n_nze = c_DBCSR_mat%n_nze 
  data_block_size = c_DBCSR_mat%dblksize 
  data_type = c_DBCSR_mat%data_type 
  nblkrows_local = c_DBCSR_mat%nblkrows_local 
  nblkrows_total = c_DBCSR_mat%nblkrows_total 
  nblkcols_total = c_DBCSR_mat%nblkcols_total 
  fullmatrix_nrows = c_DBCSR_mat%fullmatrix_nrows 
  fullmatrix_ncols = c_DBCSR_mat%fullmatrix_ncols 
  mpi_id = c_DBCSR_mat%myid 

  CALL C_F_POINTER (c_DBCSR_mat%row_blk_size, rbs, [nblkrows_total])
  CALL C_F_POINTER (c_DBCSR_mat%col_blk_size, cbs, [nblkcols_total])
  CALL C_F_POINTER (c_DBCSR_mat%row_blk_offset, rbo, [nblkrows_total+1])
  CALL C_F_POINTER (c_DBCSR_mat%col_blk_offset, cbo, [nblkcols_total+1])
  CALL C_F_POINTER (c_DBCSR_mat%row_p, row_p, [nblkrows_total+1])
  CALL C_F_POINTER (c_DBCSR_mat%col_i, col_i, [nblks])
  CALL C_F_POINTER (c_DBCSR_mat%blk_p, blk_p, [nblks])
  CALL C_F_POINTER (c_DBCSR_mat%local_rows, l_rows, [nblkrows_local])
  CALL C_F_POINTER (c_DBCSR_mat%row_dist, r_dist, [nblkrows_total])
  CALL C_F_POINTER (c_DBCSR_mat%col_dist, c_dist, [nblkcols_total])
  CALL C_F_POINTER (c_DBCSR_mat%dblock, r_dp, [data_block_size])

  CALL array_new (row_blk_size, rbs, lb=1)
  CALL array_new (col_blk_size, cbs, lb=1) 
  CALL array_new (row_blk_offset, rbo, lb=1)
  CALL array_new (col_blk_offset, cbo, lb=1)
  CALL array_new (local_rows, l_rows, lb=1)
  CALL array_new (row_dist, r_dist, lb=1)
  CALL array_new (col_dist, c_dist, lb=1)

  CALL timestop(handle)
END SUBROUTINE cDBCSR_matrix_get_info 

SUBROUTINE cDBCSR_to_dbcsr (c_DBCSR_mat, dbcsr_mat,  matrix_name, matrix_type, distribution, error)

  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  TYPE(c_DBCSR), INTENT(IN)                 :: c_DBCSR_mat
  TYPE(cp_dbcsr_type), INTENT(INOUT)        :: dbcsr_mat
  CHARACTER(LEN=*), INTENT(IN)              :: matrix_name
  CHARACTER, INTENT(IN)                     :: matrix_type
  TYPE(dbcsr_distribution_obj), INTENT(IN)  :: distribution
  TYPE(cp_error_type), INTENT(INOUT)        :: error

  TYPE(cp_logger_type), POINTER         :: logger
  TYPE(cp_dbcsr_type)                   :: dbcsr_mat_brd
  INTEGER                               :: n_nze, nblks, data_block_size, &
                                           nblkrows_local, nblkrows_total, nblkcols_total,&
                                           fullmatrix_nrows, fullmatrix_ncols, mpi_id, &
                                           i, mp_group, mynode, num_p, data_type, handle
  INTEGER, DIMENSION(:), POINTER        :: row_p, col_i, blk_p
  REAL(dp), DIMENSION(:), POINTER       :: r_dp
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: pgrid
  TYPE(array_i1d_obj)                   :: row_blk_size, col_blk_size,&
                                           row_blk_offset, col_blk_offset, &
                                           local_rows, row_dist, col_dist
  TYPE(dbcsr_data_obj)                  :: data_area
  TYPE(dbcsr_distribution_obj)          :: brd_dist
  TYPE(dbcsr_mp_obj)                    :: mp_env_new
  TYPE(dbcsr_error_type)                :: err

  CHARACTER(LEN=*), PARAMETER :: routineN = 'cDBCSR_to_dbcsr', &
      routineP = moduleN//':'//routineN

!CHARACTER :: matrix_type
!matrix_type = dbcsr_type_no_symmetry

REAL(dp) :: diff_norm

  CALL timeset(routineN,handle)
  logger => cp_error_get_logger(error)
  num_p = logger%para_env%num_pe
  mynode = logger%para_env%mepos

  CALL cDBCSR_matrix_get_info (c_DBCSR_mat, nblks, n_nze, data_block_size, data_type, &
                               nblkrows_local, nblkrows_total, nblkcols_total,&
                               fullmatrix_nrows, fullmatrix_ncols, mpi_id, &
                               row_blk_size, col_blk_size, row_blk_offset, col_blk_offset,&
                               row_p, col_i, blk_p, local_rows, row_dist, col_dist, r_dp, error)

  ALLOCATE (pgrid(num_p,1))
  pgrid = RESHAPE((/(i, i = 0, num_p-1),0/),(/num_p,1/))
  mp_group = logger%para_env%group
  CALL dbcsr_mp_new(mp_env_new, pgrid, mp_group, mynode, numnodes=num_p)
  CALL dbcsr_distribution_init(brd_dist)
  CALL dbcsr_distribution_new(brd_dist,mp_env_new,row_dist,col_dist)

  CALL dbcsr_data_init (data_area)
  CALL dbcsr_data_new (data_area, data_type)
  data_area%d%r_dp => r_dp

  CALL cp_dbcsr_init (dbcsr_mat_brd, error=error)
  CALL cp_dbcsr_create (dbcsr_mat_brd, "OMEN dbcsr_brd", brd_dist, matrix_type, &
                        row_blk_size, col_blk_size, error=error)
  CALL dbcsr_fill (data_area, row_p, col_i, blk_p, nblks, data_block_size, dbcsr_mat_brd%matrix, err)

  CALL cp_dbcsr_create(dbcsr_mat, TRIM(matrix_name)//" OMEN", distribution, matrix_type, &
                       row_blk_size, col_blk_size, error=error)
  CALL dbcsr_redistribute(dbcsr_mat_brd%matrix, dbcsr_mat%matrix, error=err)

  CALL cp_dbcsr_release(dbcsr_mat_brd,error=error)
  DEALLOCATE(pgrid)

  CALL timestop(handle)
END SUBROUTINE cDBCSR_to_dbcsr 

END MODULE c_dbcsr_methods

