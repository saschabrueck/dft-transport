!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input section for transport calculations
!>
!> \par History
!>       07.2013 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! *****************************************************************************
MODULE input_cp2k_transport

  USE f77_blas
  USE input_constants,                 ONLY: do_transport,&
                                             exper_code,&
                                             kpoint_integration,&
                                             misc_method,&
                                             scalapack_diagonalization
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_transport'

  PUBLIC :: create_transport_section

CONTAINS

  SUBROUTINE create_transport_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_transport_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,"TRANSPORT",&
            description="Specifies the parameters for transport",&
            n_keywords=15, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY (keyword)

       CALL keyword_create(keyword, name="DENSITY_MATRIX_CONSTRUCTION_METHOD",&
            description="Method used for constructing the density matrix",&
            usage="DENSITY_MATRIX_CONSTRUCTION_METHOD DIAG",&
            default_i_val=scalapack_diagonalization,&
            enum_c_vals=s2a( "DIAG", "KPOINT", "TRANSPORT", "EXPER", "MISC"),&
            enum_desc=s2a("diagonalization using ScaLapack routines",&
                          "k-point integration method",&
                          "transport code",&
                          "experimental code",&
                          "miscellaneous method"),&
            enum_i_vals=(/scalapack_diagonalization, kpoint_integration, do_transport, exper_code, misc_method/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BANDWIDTH",&
            description="The number of neighboring unit cells that one unit cell interacts with",&
            usage="BANDWIDTH 2", default_i_val=2,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="N_CELLS",&
            description="The number of unit cells",&
            usage="N_CELLS 5", default_i_val=5,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="N_ABSCISSAE",&
            description="The number of abscissae per integration interval on the real axis",&
            usage="N_ABSCISSAE 10", default_i_val=10,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="N_KPOINTS",&
            description="The number of k points for determination of the singularities",&
            usage="N_KPOINTS 10", default_i_val=10,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRA_INT_PARAM1",&
            description="extra integer parameter just in case ",&
            usage="EXTRA_INT_PARAM1 0", default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRA_INT_PARAM2",&
            description="extra integer parameter just in case ",&
            usage="EXTRA_INT_PARAM2 0", default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRA_INT_PARAM3",&
            description="extra integer parameter just in case ",&
            usage="EXTRA_INT_PARAM3 0", default_i_val=0,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COLZERO_THRESHOLD",&
            description="The smallest number that is not zero in the full diagonalization part",&
            usage="COLZERO_THRESHOLD 1.0E-12", default_r_val=1.0E-12_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_LIMIT",&
            description="The smallest eigenvalue that is kept",&
            usage="EPS_LANCZOS 1.0E-6", default_r_val=1.0E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_DECAY",&
            description="The smallest imaginary part that a decaying eigenvalue may have not to be considered as propagating",&
            usage="EPS_DECAY 1.0E-6", default_r_val=1.0E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_SINGULARITIES",&
            description="Threshold for ...",&
            usage="EPS_SINGULARITIES 1.0E-3", default_r_val=1.0E-3_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRA_PARAM1",&
            description="extra double parameter just in case ",&
            usage="EXTRA_PARAM1 0.0", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRA_PARAM2",&
            description="extra double parameter just in case ",&
            usage="EXTRA_PARAM2 0.0", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRA_PARAM3",&
            description="extra double parameter just in case ",&
            usage="EXTRA_PARAM3 0.0", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_transport_section

END MODULE input_cp2k_transport
