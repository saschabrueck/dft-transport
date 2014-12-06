!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief CP2K transport environment and its C-interoperable version
!> \par History
!>       07.2013 created [Hossein Bani-Hashemian]
!>       11.2014 revised [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE transport_env_types

  USE cp_dbcsr_interface,              ONLY: cp_csr_destroy,&
                                             cp_dbcsr_release,&
                                             cp_dbcsr_type,&
                                             csr_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_DOUBLE, C_INT, C_NULL_FUNPTR, &
                                             C_FUNPTR
                                             
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'transport_env_types'

  PUBLIC :: transport_env_type, cp2k_transport_parameters

!> DO NOT change the ORDERING or the NAMES in the following data type
  TYPE, BIND(C) :: cp2k_transport_parameters
     INTEGER(C_INT)  :: n_occ
     INTEGER(C_INT)  :: n_atoms
     REAL(C_DOUBLE)  :: evoltfactor
     INTEGER(C_INT)  :: method
     INTEGER(C_INT)  :: bandwidth
     INTEGER(C_INT)  :: n_cells
     INTEGER(C_INT)  :: n_abscissae
     INTEGER(C_INT)  :: n_kpoint
     INTEGER(C_INT)  :: num_interval
     INTEGER(C_INT)  :: num_contacts
     INTEGER(C_INT)  :: ndof
     INTEGER(C_INT)  :: tasks_per_point
     INTEGER(C_INT)  :: cores_per_node
     REAL(C_DOUBLE)  :: colzero_threshold
     REAL(C_DOUBLE)  :: eps_limit
     REAL(C_DOUBLE)  :: eps_decay
     REAL(C_DOUBLE)  :: eps_singularity_curvatures
     REAL(C_DOUBLE)  :: eps_mu
     REAL(C_DOUBLE)  :: eps_eigval_degen
     REAL(C_DOUBLE)  :: energy_interval
     REAL(C_DOUBLE)  :: min_interval
     REAL(C_DOUBLE)  :: temperature
  END TYPE cp2k_transport_parameters

  TYPE transport_env_type
     TYPE(C_FUNPTR)                   :: ext_c_method_ptr
     TYPE(cp2k_transport_parameters)  :: params 
     TYPE(cp_dbcsr_type)              :: template_matrix_sym
     TYPE(cp_dbcsr_type)              :: template_matrix_nosym
     TYPE(csr_type)                   :: s_matrix 
     TYPE(csr_type)                   :: ks_matrix 
     TYPE(csr_type)                   :: p_matrix
  END TYPE transport_env_type

  PUBLIC :: transport_env_release

CONTAINS

! *****************************************************************************
!> \brief releases the transport_env
!> \param[inout] transport_env the transport_env to be released
!> \param[inout] error         CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE transport_env_release(transport_env,error)
    TYPE(transport_env_type), POINTER        :: transport_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transport_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CPPostcondition(ASSOCIATED(transport_env),cp_failure_level,routineP,error,failure)

    CALL cp_csr_destroy(transport_env%p_matrix, error)
    CALL cp_csr_destroy(transport_env%ks_matrix, error)
    CALL cp_csr_destroy(transport_env%s_matrix, error)
    CALL cp_dbcsr_release(transport_env%template_matrix_sym, error)
    CALL cp_dbcsr_release(transport_env%template_matrix_nosym, error)

    transport_env%ext_c_method_ptr = C_NULL_FUNPTR

    DEALLOCATE(transport_env, stat=istat)
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE transport_env_release

END MODULE transport_env_types


