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
MODULE c_dbcsr_types

  USE timings,                         ONLY: timeset,&
                                             timestop

  USE, INTRINSIC :: ISO_C_BINDING

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'c_dbcsr_types'

  PUBLIC :: c_DBCSR

  TYPE, BIND(C) :: c_DBCSR
     INTEGER(C_INT) :: nblks, n_nze, dblksize, data_type,&
                       nblkrows_local, nblkrows_total, nblkcols_total,&
                       fullmatrix_nrows, fullmatrix_ncols, myid
     TYPE(C_PTR)    :: row_blk_size, col_blk_size,&
                       row_blk_offset, col_blk_offset,&
                       row_p, col_i, blk_p, local_rows,&
                       row_dist, col_dist
     TYPE(C_PTR)    :: dblock
  END TYPE c_DBCSR

  PUBLIC :: cDBCSR_nullify
  CONTAINS

  SUBROUTINE cDBCSR_nullify(c_DBCSR_mat)

    USE, INTRINSIC :: ISO_C_BINDING,      ONLY: C_NULL_PTR
    IMPLICIT NONE

    TYPE(c_DBCSR), INTENT(INOUT)          :: c_DBCSR_mat

    INTEGER                               :: handle
    CHARACTER(LEN=*), PARAMETER :: routineN = 'cDBCSR_nullify', &
        routineP = moduleN//':'//routineN
    CALL timeset(routineN,handle)

    c_DBCSR_mat%nblks = 0
    c_DBCSR_mat%n_nze = 0
    c_DBCSR_mat%dblksize = 0 
    c_DBCSR_mat%data_type = 0 
    c_DBCSR_mat%nblkrows_local = 0
    c_DBCSR_mat%nblkrows_total = 0
    c_DBCSR_mat%nblkcols_total = 0
    c_DBCSR_mat%fullmatrix_nrows = 0
    c_DBCSR_mat%fullmatrix_ncols = 0
    c_DBCSR_mat%myid = 0
    c_DBCSR_mat%row_p = C_NULL_PTR
    c_DBCSR_mat%col_i = C_NULL_PTR
    c_DBCSR_mat%blk_p = C_NULL_PTR
    c_DBCSR_mat%local_rows = C_NULL_PTR
    c_DBCSR_mat%row_blk_size = C_NULL_PTR
    c_DBCSR_mat%col_blk_size = C_NULL_PTR
    c_DBCSR_mat%row_blk_offset = C_NULL_PTR
    c_DBCSR_mat%col_blk_offset = C_NULL_PTR
    c_DBCSR_mat%dblock = C_NULL_PTR
    c_DBCSR_mat%row_dist = C_NULL_PTR
    c_DBCSR_mat%col_dist = C_NULL_PTR

    CALL timestop(handle)

  END SUBROUTINE cDBCSR_nullify 

END MODULE c_dbcsr_types
