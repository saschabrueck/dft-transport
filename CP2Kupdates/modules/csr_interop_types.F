!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief C-interoperable CSR matrices 
!> \par History
!>       05.2013 created [Hossein Bani-Hashemian]
!>       11.2014 revised into CSR matrices [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE csr_interop_types

  USE timings,                         ONLY: timeset,&
                                             timestop

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_INT, C_NULL_PTR, C_PTR

#include "./common/cp_common_uses.f90"
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csr_interop_types'

  PUBLIC :: cp2k_csr_interop_type

!> DO NOT change the ORDERING or the NAMES in the following data type
  TYPE, BIND(C) :: cp2k_csr_interop_type
    INTEGER(C_INT) :: nrows_total
    INTEGER(C_INT) :: ncols_total
    INTEGER(C_INT) :: nze_total
    INTEGER(C_INT) :: nze_local
    INTEGER(C_INT) :: nrows_local
    INTEGER(C_INT) :: data_type
    INTEGER(C_INT) :: first_row
    TYPE(C_PTR)    :: rowptr_local 
    TYPE(C_PTR)    :: colind_local
    TYPE(C_PTR)    :: nzerow_local
    TYPE(C_PTR)    :: nzvals_local
  ENDTYPE cp2k_csr_interop_type

  PUBLIC :: csr_interop_nullify

  CONTAINS

! *****************************************************************************
!> \brief nullifies (and zeroizes) a C-interoperable CSR matrix
!> \param[inout] csr_interop_mat the matrix to be nullified 
!> \param[inout] error           CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE csr_interop_nullify(csr_interop_mat, error)

    TYPE(cp2k_csr_interop_type), &
      INTENT(INOUT)                          :: csr_interop_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csr_interop_nullify', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    csr_interop_mat%nrows_total = 0
    csr_interop_mat%ncols_total = 0
    csr_interop_mat%nze_total   = 0
    csr_interop_mat%nze_local   = 0
    csr_interop_mat%nrows_local = 0
    csr_interop_mat%data_type   = 0
    csr_interop_mat%first_row   = 0
    csr_interop_mat%rowptr_local = C_NULL_PTR
    csr_interop_mat%colind_local = C_NULL_PTR
    csr_interop_mat%nzerow_local = C_NULL_PTR
    csr_interop_mat%nzvals_local = C_NULL_PTR

    CALL timestop(handle)

  END SUBROUTINE csr_interop_nullify

END MODULE csr_interop_types
