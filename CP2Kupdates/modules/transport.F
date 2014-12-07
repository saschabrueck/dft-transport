!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines for DFT+NEGF calculations (coupling with the quantum transport code OMEN)
!> \par History
!>       07.2013 created [Hossein Bani-Hashemian]
!>       11.2014 revised [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE transport

  USE cp_dbcsr_interface,              ONLY: &
       cp_convert_dbcsr_to_csr, cp_csr_create, cp_csr_create_from_dbcsr, &
       cp_dbcsr_copy, cp_dbcsr_copy_into_existing, cp_dbcsr_desymmetrize, &
       cp_dbcsr_has_symmetry, cp_dbcsr_init, cp_dbcsr_type
  USE csr_interop_methods,             ONLY: convert_csr_interop_to_dbcsr,&
                                             convert_dbcsr_to_csr_interop
  USE csr_interop_types,               ONLY: cp2k_csr_interop_type,&
                                             csr_interop_nullify
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE transport_env_types,             ONLY: cp2k_transport_parameters,&
                                             transport_env_type

  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_ASSOCIATED, C_F_PROCPOINTER, &
                                             C_PTR, C_INT

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'transport'

  PUBLIC :: transport_env_create, transport_initialize, external_scf_method 

!> interface between C/C++ and FORTRAN
  INTERFACE c_func_interface
! *****************************************************************************
!> \brief C routine that takes the S and H matrices as input and outputs a P matrix
!> \param cp2k_transport_params transport parameters read form a CP2K input file
!> \param s_mat  C-interoperable overlap matrix 
!> \param ks_mat C-interoperable Kohn-Sham matrix
!> \param p_mat  C-interoperable density matrix
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
    SUBROUTINE c_scf_routine(cp2k_transport_params, s_mat, ks_mat, p_mat) BIND(C)
      IMPORT :: C_INT, C_PTR, cp2k_csr_interop_type, cp2k_transport_parameters
      IMPLICIT NONE
      TYPE(cp2k_transport_parameters), VALUE, INTENT(IN) :: cp2k_transport_params
      TYPE(cp2k_csr_interop_type), VALUE, INTENT(IN)     :: s_mat
      TYPE(cp2k_csr_interop_type), VALUE, INTENT(IN)     :: ks_mat
      TYPE(cp2k_csr_interop_type), INTENT(INOUT)         :: p_mat
    END SUBROUTINE c_scf_routine
  END INTERFACE c_func_interface

CONTAINS

! *****************************************************************************
!> \brief creates the transport environment
!> \param[inout] qs_env the qs_env containing the transport_env
!> \param[inout] error  CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE transport_env_create(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'transport_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: input
    TYPE(transport_env_type), POINTER        :: transport_env

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env, transport_env=transport_env, input=input, error=error)

    CPPostcondition(.NOT.ASSOCIATED(transport_env),cp_failure_level,routineP,error,failure)

    ALLOCATE(transport_env,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL transport_init_read_input(input,transport_env,error)
    CALL cp_dbcsr_init(transport_env%template_matrix_sym, error)
    CALL cp_dbcsr_init(transport_env%template_matrix_nosym, error)
    CALL set_qs_env(qs_env,transport_env=transport_env,error=error)

    CALL timestop(handle)

  END SUBROUTINE transport_env_create

! *****************************************************************************
!> \brief intitializes all fields of transport_env using the parameters read from 
!>        the corresponding input section 
!> \param[inout] input         the input file 
!> \param[inout] transport_env the transport_env to be initialized
!> \param[inout] error         CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE transport_init_read_input(input,transport_env,error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(transport_env_type), INTENT(INOUT)  :: transport_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transport_init_read_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(section_vals_type), POINTER         :: transport_section

    CALL timeset(routineN,handle)

    transport_section => section_vals_get_subs_vals(input,"DFT%TRANSPORT",error=error)
    
!> read from input
    CALL section_vals_val_get(transport_section,"DENSITY_MATRIX_CONSTRUCTION_METHOD",i_val=transport_env%params%method,error=error)
    CALL section_vals_val_get(transport_section,"BANDWIDTH",i_val=transport_env%params%bandwidth,error=error)
    CALL section_vals_val_get(transport_section,"N_CELLS",i_val=transport_env%params%n_cells,error=error)
    CALL section_vals_val_get(transport_section,"N_ABSCISSAE",i_val=transport_env%params%n_abscissae,error=error)
    CALL section_vals_val_get(transport_section,"N_KPOINTS",i_val=transport_env%params%n_kpoint,error=error)
    CALL section_vals_val_get(transport_section,"NUM_INTERVAL",i_val=transport_env%params%num_interval,error=error)
    CALL section_vals_val_get(transport_section,"NUM_CONTACTS",i_val=transport_env%params%num_contacts,error=error)
    CALL section_vals_val_get(transport_section,"N_DOF",i_val=transport_env%params%ndof,error=error)
    CALL section_vals_val_get(transport_section,"TASKS_PER_POINT",i_val=transport_env%params%tasks_per_point,error=error)
    CALL section_vals_val_get(transport_section,"CORES_PER_NODE",i_val=transport_env%params%cores_per_node,error=error)
    CALL section_vals_val_get(transport_section,"COLZERO_THRESHOLD",r_val=transport_env%params%colzero_threshold,error=error)
    CALL section_vals_val_get(transport_section,"EPS_LIMIT",r_val=transport_env%params%eps_limit,error=error)
    CALL section_vals_val_get(transport_section,"EPS_DECAY",r_val=transport_env%params%eps_decay,error=error)
    CALL section_vals_val_get(transport_section,"EPS_SINGULARITY_CURVATURES",&
                                                r_val=transport_env%params%eps_singularity_curvatures,error=error)
    CALL section_vals_val_get(transport_section,"EPS_MU",r_val=transport_env%params%eps_mu,error=error)
    CALL section_vals_val_get(transport_section,"EPS_EIGVAL_DEGEN",r_val=transport_env%params%eps_eigval_degen,error=error)
    CALL section_vals_val_get(transport_section,"ENERGY_INTERVAL",r_val=transport_env%params%energy_interval,error=error)
    CALL section_vals_val_get(transport_section,"MIN_INTERVAL",r_val=transport_env%params%min_interval,error=error)
    CALL section_vals_val_get(transport_section,"TEMPERATURE",r_val=transport_env%params%temperature,error=error)

    CALL timestop(handle)

  END SUBROUTINE transport_init_read_input

! *****************************************************************************
!> \brief initializes the transport environment
!> \param[inout] transport_env the transport env to be initialized 
!> \param[in]    template_matrix   template matrix to keep the sparsity of matrices fixed
!> \param[inout] error             CP2K error
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE transport_initialize(transport_env, template_matrix, error)
    TYPE(transport_env_type), INTENT(INOUT)  :: transport_env
    TYPE(cp_dbcsr_type), INTENT(IN)          :: template_matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'transport_initialize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)  

    IF (cp_dbcsr_has_symmetry(template_matrix)) THEN
       CALL cp_dbcsr_copy(transport_env%template_matrix_sym, template_matrix, error=error)
       CALL cp_dbcsr_desymmetrize(transport_env%template_matrix_sym, transport_env%template_matrix_nosym, error)
    ELSE 
       CALL cp_dbcsr_copy(transport_env%template_matrix_nosym, template_matrix, error=error)
       CALL cp_dbcsr_copy(transport_env%template_matrix_sym, template_matrix, error=error)
    END IF

    CALL cp_csr_create_from_dbcsr(transport_env%template_matrix_nosym, transport_env%s_matrix, equal_dist=.FALSE., error=error)
    CALL cp_convert_dbcsr_to_csr(transport_env%template_matrix_nosym, transport_env%s_matrix, error=error)

    CALL cp_csr_create(transport_env%ks_matrix, transport_env%s_matrix, error)
    CALL cp_csr_create(transport_env%p_matrix , transport_env%s_matrix, error)

    CALL timestop(handle)
 
  END SUBROUTINE transport_initialize

! *****************************************************************************
!> \brief SCF calcualtion with an externally evaluated density matrix
!> \param[inout] transport_env  transport environment 
!> \param[in]    matrix_s       DBCSR overlap matrix
!> \param[in]    matrix_ks      DBCSR Kohn-Sham matrix
!> \param[inout] matrix_p       DBCSR density matrix
!> \param[in]    nelectron_spin number of electrons 
!> \param[in]    natoms         number of atoms 
!> \param[inout] error          CP2K error
!> \par History
!>       12.2012 created [Hossein Bani-Hashemian]
!>       12.2014 revised [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
   SUBROUTINE external_scf_method(transport_env, matrix_s, matrix_ks, matrix_p, nelectron_spin, natoms, error)

    TYPE(transport_env_type), INTENT(INOUT)  :: transport_env
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s, matrix_ks
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_p
    INTEGER, INTENT(IN)                      :: nelectron_spin, natoms
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'external_scf_method', &
      routineP = moduleN//':'//routineN

    TYPE(cp2k_csr_interop_type)              :: ks_mat, p_mat, s_mat

    PROCEDURE(c_scf_routine), POINTER        :: c_method
    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    CALL C_F_PROCPOINTER(transport_env%ext_c_method_ptr, c_method)
    CALL cp_assert(C_ASSOCIATED(transport_env%ext_c_method_ptr),cp_fatal_level,cp_assertion_failed,routineP,&
                  "MISSING C/C++ ROUTINE: The TRANSPORT section is meant to be used together with an external "//&
                  "program, e.g. the quantum transport code OMEN, that provides CP2K with a density matrix.",&
                  error, failure)

    transport_env%params%n_occ = nelectron_spin
    transport_env%params%n_atoms = natoms
    transport_env%params%evoltfactor = evolt

    CALL csr_interop_nullify(s_mat, error)
    CALL csr_interop_nullify(ks_mat, error)
    CALL csr_interop_nullify(p_mat, error)

    CALL cp_dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_s, error) 
    CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%s_matrix, s_mat, error=error)

    CALL cp_dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_ks, error)
    CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%ks_matrix, ks_mat, error=error)

    CALL c_method(transport_env%params, s_mat, ks_mat, p_mat)

    CALL convert_csr_interop_to_dbcsr(p_mat, transport_env%p_matrix, transport_env%template_matrix_nosym, error)
    CALL cp_dbcsr_copy(matrix_p, transport_env%template_matrix_nosym, error=error)

    CALL timestop(handle)

  END SUBROUTINE external_scf_method

END MODULE transport


