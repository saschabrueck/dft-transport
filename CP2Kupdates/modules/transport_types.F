!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief CP2K transport environment and its C-interoperable version
!> \par History
!>       07.2013 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! *****************************************************************************
MODULE transport_types
  USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                             C_INT,&
                                             c_funptr
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'transport_types'

  PUBLIC :: transport_type, c_transport_type

  TYPE transport_type
     TYPE(c_funptr)        :: ext_method_ptr
     INTEGER               :: method, bandwidth, n_cells,&
                              n_abscissae, n_kpoint,&
                              extra_int_param1, extra_int_param2, extra_int_param3
     REAL(dp)              :: colzero_threshold,&
                              eps_limit, eps_decay, eps_singularities,&
                              extra_param1, extra_param2, extra_param3
  END TYPE transport_type

  TYPE, BIND(C) :: c_transport_type
     INTEGER(C_INT)    :: method, bandwidth, n_cells, n_occ,&
                          n_abscissae, n_kpoint,&
                          extra_int_param1, extra_int_param2, extra_int_param3
     REAL(C_DOUBLE)    :: evoltfactor, colzero_threshold,&
                          eps_limit, eps_decay, eps_singularities, &
                          extra_param1, extra_param2, extra_param3
  END TYPE c_transport_type

  PUBLIC :: transport_env_release

CONTAINS

  SUBROUTINE transport_env_release(transport_env,error)
    TYPE(transport_type), POINTER            :: transport_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transport_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CPPostcondition(ASSOCIATED(transport_env),cp_failure_level,routineP,error,failure)

    DEALLOCATE(transport_env,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE transport_env_release

END MODULE transport_types


