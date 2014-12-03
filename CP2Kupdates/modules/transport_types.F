!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief CP2K transport environment and its C-interoperable version
!> \par History
!>       07.2013 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE transport_types
  USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                             C_INT,&
                                             C_BOOL,&
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
     TYPE(c_funptr)    :: ext_c_method_ptr
     INTEGER           :: method, bandwidth, n_cells,&
                          n_abscissae, n_kpoint,&
                          num_interval, num_contacts, ndof, tasks_per_point
     LOGICAL           :: real_axis 
     REAL(dp)          :: colzero_threshold,&
                          eps_limit, eps_decay, eps_singularity_curvatures, &
                          eps_mu, eps_eigval_degen, energy_interval, min_interval,&
                          temperature
  END TYPE transport_type

  TYPE, BIND(C) :: c_transport_type
     INTEGER(C_INT)    :: n_occ, n_atoms
     REAL(C_DOUBLE)    :: evoltfactor
     INTEGER(C_INT)    :: method, bandwidth, n_cells,&
                          n_abscissae, n_kpoint,&
                          num_interval, num_contacts, ndof, tasks_per_point
     LOGICAL(C_BOOL)   :: real_axis 
     REAL(C_DOUBLE)    :: colzero_threshold,&
                          eps_limit, eps_decay, eps_singularity_curvatures, &
                          eps_mu, eps_eigval_degen, energy_interval, min_interval,&
                          temperature
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


