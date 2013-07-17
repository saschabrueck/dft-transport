!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  
!> \par History
!>       2013.07 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! *****************************************************************************
MODULE transport_init

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
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'transport_init'

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

    ! get basic quantities from the qs_env
    CALL get_qs_env(qs_env, transport_env=transport_env, input=input, error=error)

    CPPostcondition(.NOT.ASSOCIATED(transport_env),cp_failure_level,routineP,error,failure)

    ALLOCATE(transport_env,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! copy some basic stuff

    ! parse the transport section and set derived quantities
    CALL transport_init_read_input(input,transport_env,error)

    ! put the ls_scf_env in qs_env
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

    ! should come from input
    CALL section_vals_val_get(transport_section,"DENSITY_MATRIX_CONSTRUCTION_METHOD",i_val=transport_env%method,error=error)
    CALL section_vals_val_get(transport_section,"BANDWIDTH",i_val=transport_env%bandwidth,error=error)
    CALL section_vals_val_get(transport_section,"N_CELLS",i_val=transport_env%n_cells,error=error)
    CALL section_vals_val_get(transport_section,"N_ABSCISSAE",i_val=transport_env%n_abscissae,error=error)
    CALL section_vals_val_get(transport_section,"N_KPOINTS",i_val=transport_env%n_kpoint,error=error)
    CALL section_vals_val_get(transport_section,"EXTRA_INT_PARAM1",i_val=transport_env%extra_int_param1,error=error)
    CALL section_vals_val_get(transport_section,"EXTRA_INT_PARAM2",i_val=transport_env%extra_int_param2,error=error)
    CALL section_vals_val_get(transport_section,"EXTRA_INT_PARAM3",i_val=transport_env%extra_int_param3,error=error)
    CALL section_vals_val_get(transport_section,"COLZERO_THRESHOLD",r_val=transport_env%colzero_threshold,error=error)
    CALL section_vals_val_get(transport_section,"EPS_LIMIT",r_val=transport_env%eps_limit,error=error)
    CALL section_vals_val_get(transport_section,"EPS_DECAY",r_val=transport_env%eps_decay,error=error)
    CALL section_vals_val_get(transport_section,"EPS_SINGULARITIES",r_val=transport_env%eps_singularities,error=error)
    CALL section_vals_val_get(transport_section,"EXTRA_PARAM1",r_val=transport_env%extra_param1,error=error)
    CALL section_vals_val_get(transport_section,"EXTRA_PARAM2",r_val=transport_env%extra_param2,error=error)
    CALL section_vals_val_get(transport_section,"EXTRA_PARAM3",r_val=transport_env%extra_param3,error=error)

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

  SUBROUTINE c_transport_type_init(transport_env, n_occupied, c_transport_env, error)
    TYPE(transport_type), POINTER            :: transport_env
    INTEGER                                  :: n_occupied
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
    c_transport_env%extra_int_param1=transport_env%extra_int_param1
    c_transport_env%extra_int_param2=transport_env%extra_int_param2
    c_transport_env%extra_int_param3=transport_env%extra_int_param3
    c_transport_env%colzero_threshold=transport_env%colzero_threshold
    c_transport_env%eps_limit=transport_env%eps_limit
    c_transport_env%eps_decay=transport_env%eps_decay
    c_transport_env%eps_singularities=transport_env%eps_singularities
    c_transport_env%extra_param1=transport_env%extra_param1
    c_transport_env%extra_param2=transport_env%extra_param2
    c_transport_env%extra_param3=transport_env%extra_param3
 
    c_transport_env%n_occ = n_occupied
    c_transport_env%evoltfactor = evolt

    CALL timestop(handle)
 
  END SUBROUTINE c_transport_type_init

END MODULE transport_init


