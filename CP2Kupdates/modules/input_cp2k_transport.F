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
            n_keywords=19, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
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
            description="The number of neighboring unit cells that one unit cell interacts with.",&
            usage="BANDWIDTH <integer>", default_i_val=2, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="N_CELLS",&
            description="The number of unit cells.",&
            usage="N_CELLS <integer>", default_i_val=5, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="N_ABSCISSAE",&
            description="The number of abscissae per integration interval on the real axis.",&
            usage="N_ABSCISSAE <integer>", default_i_val=10, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="N_KPOINTS",&
            description="The number of k points for determination of the singularities.",&
            usage="N_KPOINTS <integer>", default_i_val=10, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="NUM_INTERVAL",&
            description="Max number of energy points per small interval.",&
            usage="NUM_INTERVAL <integer>", default_i_val=0, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="NUM_CONTACTS",&
            description="The number of contacts.",&
            usage="NUM_CONTACTS <integer>", default_i_val=2, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="N_DOF",&
            description="Number of degrees of freedom for the contact unit cell.",&
            usage="N_DOF <integer>", default_i_val=0, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="TASKS_PER_POINT",&
            description="Number of tasks per energy point.",&
            usage="TASKS_PER_POINT <integer>", default_i_val=0, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="REAL_AXIS",&
            description="Integration over real axis only.",&
            usage="REAL_AXIS .TRUE.", default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="COLZERO_THRESHOLD",&
            description="The smallest number that is not zero in the full diagonalization part.",&
            usage="COLZERO_THRESHOLD <real>", default_r_val=1.0E-12_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="EPS_LIMIT",&
            description="The smallest eigenvalue that is kept.",&
            usage="EPS_LIMIT <real>", default_r_val=1.0E-6_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="EPS_DECAY",&
            description="The smallest imaginary part that a decaying eigenvalue may have not to be considered as propagating.",&
            usage="EPS_DECAY <real>", default_r_val=1.0E-6_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EPS_SINGULARITY_CURVATURES",&
            description="Filter for degenerate bands in the bandstructure.",&
            usage="EPS_SINGULARITY_CURVATURES <real>", default_r_val=1.0E-3_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="EPS_MU",&
            description="Accuracy to which the Fermi level should be determined.",&
            usage="EPS_MU <real>", default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="EPS_EIGVAL_DEGEN",&
            description="Filter for degenerate bands in the injection vector.",&
            usage="EPS_EIGVAL_DEGEN <real>", default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="ENERGY_INTERVAL",&
            description="Average distance for big intervals in energy vector.",&
            usage="ENERGY_INTERVAL <real>", default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="MIN_INTERVAL",&
            description="Smallest enery distance in energy vector.",&
            usage="MIN_INTERVAL <real>", default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

       CALL keyword_create(keyword, name="TEMPERATURE",&
            description="Temperature.",&
            usage="TEMPERATURE [K] 300.0", default_r_val=300.0_dp, error=error)
       CALL section_add_keyword(section, keyword, error=error)
       CALL keyword_release(keyword, error=error)

    END IF

  END SUBROUTINE create_transport_section

END MODULE input_cp2k_transport
