!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input for the transport
!> \author Hossein Bani-Hashemian
! *****************************************************************************
MODULE input_cp2k_transport

  USE f77_blas
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
  USE kinds,                           ONLY: dp
  USE input_constants,                 ONLY: misc_method,&
                                             scalapack_diagonalization,&
                                             kpoint_integration,&
                                             exper_code,&
                                             do_transport
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_transport'

  PUBLIC :: create_transport_section

CONTAINS
! *****************************************************************************
!> \brief creates the transport section
!> \author Hossein Bani-Hashemian [2013-07]
! *****************************************************************************
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
            n_keywords=10, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
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
            description="bandwidth ... ",&
            usage="BANDWIDTH 3", default_i_val=3,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NCELLS",&
            description="ncells ...",&
            usage="NCELLS 5", default_i_val=5,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COLZERO_THRESHOLD",&
            description="Threshold used for ...",&
            usage="COLZERO_THRESHOLD 1.0E-12", default_r_val=1.0E-12_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_LIMIT",&
            description="Threshold used for ...",&
            usage="EPS_LANCZOS 1.0E-6", default_r_val=1.0E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_DECAY",&
            description="Threshold used for ....",&
            usage="EPS_DECAY 1.0E-6", default_r_val=1.0E-6_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_SINGULARITIES",&
            description="threshold for ...",&
            usage="EPS_SINGULARITIES 1.0E-3", default_r_val=1.0E-3_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EXTRA_PARAM1",&
            description="extra double parameter just in case ",&
            usage="EXTRA_PARAM1 0.0", default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_transport_section

END MODULE input_cp2k_transport
