!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines for the transport environment and its C-interoperable version  
!> \par History
!>       07.2013 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE transport

  USE input_constants,                 ONLY: do_transport,&
                                             exper_code,&
                                             kpoint_integration,&
                                             misc_method,&
                                             scalapack_diagonalization
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE transport_types,                 ONLY: c_transport_type,&
                                             transport_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'transport'

  PUBLIC :: transport_env_create, c_transport_type_init

CONTAINS

  SUBROUTINE transport_env_create(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'transport_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: input
    TYPE(transport_type), POINTER            :: transport_env

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env, transport_env=transport_env, input=input, error=error)

    CPPostcondition(.NOT.ASSOCIATED(transport_env),cp_failure_level,routineP,error,failure)

    ALLOCATE(transport_env,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL transport_init_read_input(input,transport_env,error)
    CALL set_qs_env(qs_env,transport_env=transport_env,error=error)

    CALL timestop(handle)

  END SUBROUTINE transport_env_create

  SUBROUTINE transport_init_read_input(input,transport_env,error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(transport_type), INTENT(INOUT)      :: transport_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transport_init_read_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: transport_section

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    failure=.FALSE.

    transport_section => section_vals_get_subs_vals(input,"DFT%TRANSPORT",error=error)
    
    ! from input
    CALL section_vals_val_get(transport_section,"DENSITY_MATRIX_CONSTRUCTION_METHOD",i_val=transport_env%method,error=error)
    CALL section_vals_val_get(transport_section,"BANDWIDTH",i_val=transport_env%bandwidth,error=error)
    CALL section_vals_val_get(transport_section,"N_CELLS",i_val=transport_env%n_cells,error=error)
    CALL section_vals_val_get(transport_section,"N_ABSCISSAE",i_val=transport_env%n_abscissae,error=error)
    CALL section_vals_val_get(transport_section,"N_KPOINTS",i_val=transport_env%n_kpoint,error=error)
    CALL section_vals_val_get(transport_section,"NUM_INTERVAL",i_val=transport_env%num_interval,error=error)
    CALL section_vals_val_get(transport_section,"NUM_CONTACTS",i_val=transport_env%num_contacts,error=error)
    CALL section_vals_val_get(transport_section,"N_DOF",i_val=transport_env%ndof,error=error)
    CALL section_vals_val_get(transport_section,"TASKS_PER_POINT",i_val=transport_env%tasks_per_point,error=error)
    CALL section_vals_val_get(transport_section,"REAL_AXIS",l_val=transport_env%real_axis,error=error)
    CALL section_vals_val_get(transport_section,"COLZERO_THRESHOLD",r_val=transport_env%colzero_threshold,error=error)
    CALL section_vals_val_get(transport_section,"EPS_LIMIT",r_val=transport_env%eps_limit,error=error)
    CALL section_vals_val_get(transport_section,"EPS_DECAY",r_val=transport_env%eps_decay,error=error)
    CALL section_vals_val_get(transport_section,"EPS_SINGULARITY_CURVATURES",&
                                                r_val=transport_env%eps_singularity_curvatures,error=error)
    CALL section_vals_val_get(transport_section,"EPS_MU",r_val=transport_env%eps_mu,error=error)
    CALL section_vals_val_get(transport_section,"EPS_EIGVAL_DEGEN",r_val=transport_env%eps_eigval_degen,error=error)
    CALL section_vals_val_get(transport_section,"ENERGY_INTERVAL",r_val=transport_env%energy_interval,error=error)
    CALL section_vals_val_get(transport_section,"MIN_INTERVAL",r_val=transport_env%min_interval,error=error)
    CALL section_vals_val_get(transport_section,"TEMPERATURE",r_val=transport_env%temperature,error=error)

    SELECT CASE(transport_env%method)
    CASE(misc_method)
    CASE(scalapack_diagonalization)
    CASE(kpoint_integration)
    CASE(exper_code)
    CASE(do_transport)
    CASE DEFAULT
        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE transport_init_read_input

  SUBROUTINE c_transport_type_init(transport_env, n_occupied, natoms, c_transport_env, error)
    TYPE(transport_type), POINTER            :: transport_env
    INTEGER, INTENT(IN)                      :: n_occupied, natoms
    TYPE(c_transport_type), INTENT(INOUT)    :: c_transport_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'c_transport_type_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)   

    c_transport_env%method=transport_env%method
    c_transport_env%bandwidth=transport_env%bandwidth
    c_transport_env%n_cells=transport_env%n_cells
    c_transport_env%n_abscissae=transport_env%n_abscissae
    c_transport_env%n_kpoint=transport_env%n_kpoint
    c_transport_env%num_interval=transport_env%num_interval
    c_transport_env%num_contacts=transport_env%num_contacts
    c_transport_env%ndof=transport_env%ndof
    c_transport_env%tasks_per_point=transport_env%tasks_per_point
    c_transport_env%real_axis=transport_env%real_axis
    c_transport_env%colzero_threshold=transport_env%colzero_threshold
    c_transport_env%eps_limit=transport_env%eps_limit
    c_transport_env%eps_decay=transport_env%eps_decay
    c_transport_env%eps_singularity_curvatures=transport_env%eps_singularity_curvatures
    c_transport_env%eps_mu=transport_env%eps_mu
    c_transport_env%eps_eigval_degen=transport_env%eps_eigval_degen
    c_transport_env%energy_interval=transport_env%energy_interval
    c_transport_env%min_interval=transport_env%min_interval
    c_transport_env%temperature=transport_env%temperature
 
    c_transport_env%n_occ = n_occupied
    c_transport_env%n_atoms = natoms
    c_transport_env%evoltfactor = evolt

    CALL timestop(handle)
 
  END SUBROUTINE c_transport_type_init

END MODULE transport


