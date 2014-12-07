!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief methods to work with C-interoperable CSR matrices 
!> \par History
!>       05.2013 created [Hossein Bani-Hashemian]
!>       11.2014 revised into CSR matrices [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE csr_interop_methods

  USE cp_dbcsr_interface,              ONLY: cp_convert_csr_to_dbcsr,&
                                             cp_convert_dbcsr_to_csr,&
                                             cp_dbcsr_type,&
                                             csr_type
  USE csr_interop_types,               ONLY: cp2k_csr_interop_type,&
                                             csr_interop_matrix_get_info
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_allgather
  USE timings,                         ONLY: timeset,&
                                             timestop

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_ASSOCIATED, C_LOC, C_PTR, &
                                             C_DOUBLE, C_INT, C_NULL_PTR

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PUBLIC :: convert_dbcsr_to_csr_interop, convert_csr_interop_to_dbcsr

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csr_interop_methods'

CONTAINS

! *****************************************************************************
!> \brief converts a DBCSR matrix to a C-interoperable CSR matrix
!> \param[in]    dbcsr_mat  DBCSR matrix to be converted 
!> \param[inout] csr_mat    auxiliary CSR matrix
!> \param[inout] csr_interop_mat C-interoperable CSR matrix
!> \param[inout] error      CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE convert_dbcsr_to_csr_interop(dbcsr_mat, csr_mat, csr_interop_mat, error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: dbcsr_mat
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(cp2k_csr_interop_type), INTENT(INOUT)      :: csr_interop_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_dbcsr_to_csr_interop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, mepos, num_pe, mp_group
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nrows_local_all, first_row_all
    INTEGER(C_INT), DIMENSION(:), POINTER    :: colind_local, rowptr_local, nzerow_local
    REAL(C_DOUBLE), DIMENSION(:), POINTER    :: nzvals_local
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    mp_group = logger%para_env%group
    mepos = logger%para_env%mepos
    num_pe = logger%para_env%num_pe

!> dbcsr to csr
    CALL cp_convert_dbcsr_to_csr(dbcsr_mat, csr_mat, error=error)

!> csr to csr_interop
    rowptr_local => csr_mat%rowptr_local
    colind_local => csr_mat%colind_local
    nzerow_local => csr_mat%nzerow_local
    nzvals_local => csr_mat%nzval_local%r_dp !> support real double percision for now

    IF (.NOT. C_ASSOCIATED(csr_interop_mat%rowptr_local)) THEN 
       IF (SIZE(rowptr_local) .EQ. 0) THEN 
          csr_interop_mat%rowptr_local = C_NULL_PTR 
       ELSE
          csr_interop_mat%rowptr_local = C_LOC(rowptr_local(1))  
       END IF 
    END IF

    IF (.NOT. C_ASSOCIATED(csr_interop_mat%colind_local)) THEN 
       IF (SIZE(colind_local) .EQ. 0) THEN 
          csr_interop_mat%colind_local = C_NULL_PTR 
       ELSE
          csr_interop_mat%colind_local = C_LOC(colind_local(1))  
       END IF 
    END IF

    IF (.NOT. C_ASSOCIATED(csr_interop_mat%nzerow_local)) THEN 
       IF (SIZE(nzerow_local) .EQ. 0) THEN 
          csr_interop_mat%nzerow_local = C_NULL_PTR 
       ELSE
          csr_interop_mat%nzerow_local = C_LOC(nzerow_local(1))  
       END IF 
    END IF

    IF (.NOT. C_ASSOCIATED(csr_interop_mat%nzvals_local)) THEN 
       IF (SIZE(nzvals_local) .EQ. 0) THEN 
          csr_interop_mat%nzvals_local = C_NULL_PTR 
       ELSE
          csr_interop_mat%nzvals_local = C_LOC(nzvals_local(1))  
       END IF 
    END IF

    ALLOCATE(nrows_local_all(0:num_pe-1), first_row_all(0:num_pe-1))
    CALL mp_allgather(csr_mat%nrows_local, nrows_local_all, mp_group) 
    CALL cumsum_i(nrows_local_all, first_row_all)   

    IF (mepos .EQ. 0) THEN
       csr_interop_mat%first_row = 0
    ELSE
       csr_interop_mat%first_row = first_row_all(mepos - 1)
    END IF
    csr_interop_mat%nrows_total = csr_mat%nrows_total
    csr_interop_mat%ncols_total = csr_mat%ncols_total
    csr_interop_mat%nze_local = csr_mat%nze_local
    csr_interop_mat%nze_total = csr_mat%nze_total
    csr_interop_mat%nrows_local = csr_mat%nrows_local
    csr_interop_mat%data_type = csr_mat%nzval_local%data_type

    CALL timestop(handle)

    CONTAINS
! *****************************************************************************
!> \brief cumulative sum of a 1d array of integers
!> \param[in]  arr    input array
!> \param[out] cumsum cumulative sum of the input array
! *****************************************************************************
    SUBROUTINE cumsum_i(arr,cumsum)
    INTEGER, DIMENSION(:), INTENT(IN)        :: arr
    INTEGER, DIMENSION(SIZE(arr)), &
      INTENT(OUT)                            :: cumsum

    INTEGER                                  :: i

       cumsum(1)=arr(1)
       DO i=2,SIZE(arr)
          cumsum(i)=cumsum(i-1)+arr(i)
       END DO
    END SUBROUTINE cumsum_i

  END SUBROUTINE convert_dbcsr_to_csr_interop

! *****************************************************************************
!> \brief converts a C-interoperable CSR matrix to a DBCSR matrix
!> \param[in] csr_interop_mat C-interoperable CSR matrix
!> \param[inout] csr_mat         auxiliary CSR matrix
!> \param[inout] dbcsr_mat       DBCSR matrix
!> \param[inout] error           CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE convert_csr_interop_to_dbcsr(csr_interop_mat, csr_mat, dbcsr_mat, error)

    TYPE(cp2k_csr_interop_type), INTENT(IN)  :: csr_interop_mat
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: dbcsr_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_csr_interop_to_dbcsr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_type, handle, &
                                                ncols_total, nrows_local, &
                                                nrows_total, nze_local, &
                                                nze_total
    INTEGER, DIMENSION(:), POINTER           :: colind_local, nzerow_local, &
                                                rowptr_local
    REAL(dp), DIMENSION(:), POINTER          :: nzvals_local

    CALL timeset(routineN,handle)

!> csr_interop to csr
    CALL csr_interop_matrix_get_info(csr_interop_mat, &
                   nrows_total = nrows_total, ncols_total = ncols_total, nze_local = nze_local, &
                   nze_total = nze_total, nrows_local = nrows_local, data_type = data_type,&
                   rowptr_local = rowptr_local, colind_local = colind_local, &
                   nzerow_local = nzerow_local, nzvals_local = nzvals_local, &
                   error=error)

    csr_mat%nrows_total = nrows_total
    csr_mat%ncols_total = ncols_total
    csr_mat%nze_local = nze_local
    csr_mat%nze_total = nze_total
    csr_mat%nrows_local = nrows_local
    csr_mat%nzval_local%data_type = data_type

    csr_mat%rowptr_local = rowptr_local
    csr_mat%colind_local = colind_local
    csr_mat%nzerow_local = nzerow_local
    csr_mat%nzval_local%r_dp = nzvals_local

!> csr to dbcsr
    CALL cp_convert_csr_to_dbcsr(dbcsr_mat, csr_mat, error=error)

!> release memory allocated for csr_interop_mat
    DEALLOCATE(rowptr_local, colind_local, nzerow_local, nzvals_local)

    CALL timestop(handle)

  END SUBROUTINE convert_csr_interop_to_dbcsr 

END MODULE csr_interop_methods

