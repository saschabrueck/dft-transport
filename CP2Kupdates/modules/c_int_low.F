!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief interface (based on f77_int_low.F) to use cp2k as library in C/C++
!> \par History
!>      12.2012 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! *****************************************************************************
SUBROUTINE cp_init_cp2k(init_mpi,ierr) BIND(C,name="cp_c_init_cp2k")
  USE f77_interface,               ONLY: icp => init_cp2k
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int)  :: init_mpi,ierr
  LOGICAL         :: i_mpi
  i_mpi=init_mpi/=0
  CALL icp(i_mpi,ierr)
END SUBROUTINE cp_init_cp2k

! *****************************************************************************
SUBROUTINE cp_finalize_cp2k(finalize_mpi,ierr) BIND(C,name="cp_c_finalize_cp2k")
  USE f77_interface,               ONLY: kcp => finalize_cp2k
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int) :: finalize_mpi, ierr
  LOGICAL        :: f_mpi
  f_mpi=finalize_mpi/=0
  CALL kcp(f_mpi,ierr)
END SUBROUTINE cp_finalize_cp2k

! *****************************************************************************
SUBROUTINE cp_create_fenv(new_env_id,input_file_path,output_file_path,&
     ierr) BIND(C,name="cp_c_create_fenv")
  USE f77_interface,               ONLY: cfe => create_force_env
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_char, c_null_char
  USE kinds,                       ONLY: default_path_length
  IMPLICIT NONE
  CHARACTER(len=1,kind=c_char) :: input_file_path(*), output_file_path(*)
  INTEGER(c_int)               :: new_env_id, ierr
  CHARACTER(len=default_path_length) :: ifp
  CHARACTER(len=default_path_length) :: ofp
  INTEGER  :: i,j

  ifp = " " ; ofp = " "
  i = 1 ; j = 1
  DO WHILE (input_file_path(i) .NE. c_null_char)
     ifp(i:i) = input_file_path(i)
     i = i+1
  END DO
  DO WHILE (output_file_path(j) .NE. c_null_char)
     ofp(j:j) = output_file_path(j)
     j = j+1
  END DO

  CALL cfe(new_env_id,ifp,ofp,ierr=ierr)
END SUBROUTINE cp_create_fenv

! *****************************************************************************
SUBROUTINE cp_create_fenv_comm(new_env_id,input_file_path,output_file_path,&
     mpi_comm,ierr) BIND(C,name="cp_c_create_fenv_comm")
  USE f77_interface,               ONLY: cfe => create_force_env
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_char, c_null_char
  USE kinds,                       ONLY: default_path_length
  IMPLICIT NONE
  CHARACTER(len=1,kind=c_char) :: input_file_path(*), output_file_path(*)
  INTEGER(c_int)               :: new_env_id, mpi_comm, ierr
  CHARACTER(len=default_path_length) :: ifp
  CHARACTER(len=default_path_length) :: ofp
  INTEGER  :: i,j

  ifp = " " ; ofp = " "
  i = 1 ; j = 1
  DO WHILE (input_file_path(i) .NE. c_null_char)
     ifp(i:i) = input_file_path(i)
     i = i+1
  END DO
  DO WHILE (output_file_path(j) .NE. c_null_char)
     ofp(j:j) = output_file_path(j)
     j = j+1
  END DO

  CALL cfe(new_env_id,ifp,ofp,mpi_comm,ierr=ierr)
END SUBROUTINE cp_create_fenv_comm

! *****************************************************************************
RECURSIVE SUBROUTINE cp_destroy_fenv(env_id, ierr) BIND(C,name="cp_c_destroy_fenv") 
  USE f77_interface,               ONLY: dfe => destroy_force_env
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int) :: env_id, ierr

  CALL dfe(env_id,ierr)
END SUBROUTINE cp_destroy_fenv

! *****************************************************************************
SUBROUTINE cp_set_pos(env_id, new_pos, n_el, ierr) BIND(C,name="cp_c_set_pos")
  USE kinds,                       ONLY: dp
  USE f77_interface,               ONLY: sp => set_pos
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_double
  IMPLICIT NONE
  INTEGER(c_int) :: n_el, env_id, ierr
  REAL(c_double), DIMENSION(0:n_el-1) :: new_pos    
  REAL(dp), DIMENSION(1:n_el)  :: fortran_new_pos 

  fortran_new_pos = new_pos 
  CALL sp(env_id,fortran_new_pos,n_el,ierr)
END SUBROUTINE cp_set_pos

! *****************************************************************************
SUBROUTINE cp_set_vel(env_id, new_vel, n_el, ierr) BIND(C,name="cp_c_set_vel")
  USE kinds,                       ONLY: dp
  USE f77_interface,               ONLY: sv => set_vel
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_double
  IMPLICIT NONE
  INTEGER(c_int) :: n_el, env_id, ierr
  REAL(c_double), DIMENSION(0:n_el-1) :: new_vel
  REAL(dp), DIMENSION(1:n_el)  :: fortran_new_vel

  fortran_new_vel = new_vel
  CALL sv(env_id,fortran_new_vel,n_el,ierr)
END SUBROUTINE cp_set_vel

! *****************************************************************************
SUBROUTINE cp_get_result_r1(env_id, description, N, RESULT, ierr) BIND(C,name="cp_c_get_result_r1")
  USE kinds,                       ONLY: dp, default_string_length
  USE f77_interface,               ONLY: get_result_r1
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_char, c_double, c_null_char
  IMPLICIT NONE
  CHARACTER(len=1,kind=c_char) :: description(*)
  INTEGER(c_int) :: N
  REAL(c_double), DIMENSION(0:N-1) :: RESULT
  REAL(dp), DIMENSION(1:N)  :: fortran_RESULT
  INTEGER(c_int) :: env_id, ierr
!  CHARACTER(LEN=default_string_length) :: description_low
  CHARACTER(len=default_string_length) :: desc_low
  INTEGER  :: i

!  description_low=description
  desc_low = " "
  i = 1
  DO WHILE (description(i) .NE. c_null_char)
     desc_low(i:i) = description(i)
     i = i+1
  END DO

  CALL get_result_r1(env_id,desc_low, N, fortran_RESULT, ierr)
  RESULT = fortran_RESULT
END SUBROUTINE cp_get_result_r1

! *****************************************************************************
SUBROUTINE cp_get_natom(env_id, natom, ierr) BIND(C,name="cp_c_get_natom")
  USE f77_interface,               ONLY: gna => get_natom
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int) :: natom, env_id, ierr

  CALL gna(env_id,natom,ierr)
END SUBROUTINE cp_get_natom

! *****************************************************************************
SUBROUTINE cp_get_nparticle(env_id, nparticle, ierr) BIND(C,name="cp_c_get_nparticle")
  USE f77_interface,               ONLY: gnp => get_nparticle
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int) :: nparticle, env_id, ierr

  CALL gnp(env_id,nparticle,ierr)
END SUBROUTINE cp_get_nparticle

! *****************************************************************************
SUBROUTINE cp_get_pos(env_id, pos, n_el, ierr) BIND(C,name="cp_c_get_pos")
  USE kinds,                       ONLY: dp
  USE f77_interface,               ONLY: gp => get_pos
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_double
  IMPLICIT NONE
  INTEGER(c_int) :: n_el, env_id, ierr
  REAL(c_double), DIMENSION(0:n_el-1) :: pos
  REAL(dp), DIMENSION(1:n_el) :: fortran_pos

  CALL gp(env_id,fortran_pos,n_el,ierr)
  pos = fortran_pos
END SUBROUTINE cp_get_pos

! *****************************************************************************
SUBROUTINE cp_get_force(env_id,force, n_el, ierr) BIND(C,name="cp_c_get_force")
  USE kinds,                       ONLY: dp
  USE f77_interface,               ONLY: gf => get_force
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_double
  IMPLICIT NONE
  INTEGER(c_int) :: n_el, env_id, ierr
  REAL(c_double), DIMENSION(0:n_el-1) :: force
  REAL(dp), DIMENSION(1:n_el) :: fortran_force

  CALL gf(env_id,fortran_force,n_el,ierr)
  force = fortran_force
END SUBROUTINE cp_get_force

! *****************************************************************************
RECURSIVE SUBROUTINE cp_get_energy(env_id, e_pot, ierr) BIND(C,name="cp_c_get_energy")
  USE kinds,                       ONLY: dp
  USE f77_interface,               ONLY: ge => get_energy
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_double
  IMPLICIT NONE
  INTEGER(c_int) :: env_id, ierr
  REAL(c_double) :: e_pot

  CALL ge(env_id,e_pot,ierr)
END SUBROUTINE cp_get_energy

! *****************************************************************************
RECURSIVE SUBROUTINE cp_calc_energy_force(env_id,calc_force,ierr) BIND(C,name="cp_c_calc_energy_force")
  USE f77_interface,               ONLY: cef => calc_energy_force
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_bool
  IMPLICIT NONE
  INTEGER(c_int) :: calc_force, env_id, ierr
  LOGICAL(c_bool) :: cf
  LOGICAL(kind=4) :: fortran_cf

  cf=(calc_force/=0)
  fortran_cf = LOGICAL(cf,kind=4)
  CALL cef(env_id,fortran_cf,ierr)
END SUBROUTINE cp_calc_energy_force

! *****************************************************************************
RECURSIVE SUBROUTINE cp_calc_energy(env_id,pos,n_el,e_pot,ierr) BIND(C,name="cp_c_calc_energy")
  USE kinds,                       ONLY: dp
  USE f77_interface,               ONLY: ce => calc_energy
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_double
  IMPLICIT NONE
  INTEGER(c_int) :: env_id, ierr, n_el
  REAL(c_double) :: e_pot
  REAL(c_double), DIMENSION(0:n_el-1) :: pos
  REAL(dp), DIMENSION(1:n_el) :: fortran_pos

  fortran_pos = pos
!write(*,*) fortran_pos
  CALL ce(env_id,fortran_pos,n_el,e_pot,ierr)
END SUBROUTINE cp_calc_energy

! *****************************************************************************
SUBROUTINE cp_calc_force(env_id,pos,n_el_pos,e_pot,force,n_el_force,ierr)&
  BIND(C,name="cp_c_calc_force")
  USE kinds,                       ONLY: dp
  USE f77_interface,               ONLY: cf => calc_force
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_double
  IMPLICIT NONE
  INTEGER(c_int) :: env_id, ierr, n_el_pos, n_el_force
  REAL(c_double) :: e_pot
  REAL(c_double), DIMENSION(0:n_el_pos-1) :: pos
  REAL(c_double), DIMENSION(0:n_el_force-1) :: force
  REAL(dp), DIMENSION(1:n_el_pos) :: fortran_pos
  REAL(dp), DIMENSION(1:n_el_force) :: fortran_force

  fortran_pos = pos
  fortran_force = force
  CALL cf(env_id,fortran_pos,n_el_pos,e_pot,fortran_force,n_el_force,ierr)
  force = fortran_force
END SUBROUTINE cp_calc_force

! *****************************************************************************
SUBROUTINE cp_run_input(input_file_path,output_file_path,ierr) BIND(C,name="cp_c_run_input")
  USE cp2k_runs,                   ONLY: rin => run_input
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_char, c_null_char
  USE kinds,                       ONLY: default_path_length
  IMPLICIT NONE
  CHARACTER(len=1,kind=c_char) :: input_file_path(*), output_file_path(*)
  INTEGER(c_int) :: ierr
  CHARACTER(len=default_path_length) :: ifp
  CHARACTER(len=default_path_length) :: ofp
  INTEGER  :: i,j

  ifp = " " ; ofp = " "
  i = 1 ; j = 1
  DO WHILE (input_file_path(i) .NE. c_null_char)
     ifp(i:i) = input_file_path(i)
     i = i+1
  END DO
  DO WHILE (output_file_path(j) .NE. c_null_char)
     ofp(j:j) = output_file_path(j)
     j = j+1
  END DO

  CALL rin(ifp,ofp,ierr=ierr)
END SUBROUTINE cp_run_input

! *****************************************************************************
RECURSIVE SUBROUTINE cp_run_input_comm(input_file_path,output_file_path,&
     mpi_comm,ierr) BIND(C,name="cp_c_run_input_comm")
  USE cp2k_runs,                   ONLY: rin => run_input
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_char, c_null_char
  USE kinds,                       ONLY: default_path_length
  IMPLICIT NONE
  CHARACTER(len=1,kind=c_char) :: input_file_path(*), output_file_path(*)
  INTEGER(c_int) :: mpi_comm, ierr
  CHARACTER(len=default_path_length) :: ifp
  CHARACTER(len=default_path_length) :: ofp
  INTEGER  :: i,j

  ifp = " " ; ofp = " "
  i = 1 ; j = 1
  DO WHILE (input_file_path(i) .NE. c_null_char)
     ifp(i:i) = input_file_path(i)
     i = i+1
  END DO
  DO WHILE (output_file_path(j) .NE. c_null_char)
     ofp(j:j) = output_file_path(j)
     j = j+1
  END DO

  CALL rin(ifp,ofp,ierr,mpi_comm)
END SUBROUTINE cp_run_input_comm

! *****************************************************************************
RECURSIVE SUBROUTINE cp_rep_init(rep_env_id,ierr) BIND(C,name="cp_c_rep_init")
  USE replica_methods,             ONLY: rpin => rep_env_init_low
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int) :: rep_env_id,ierr

  CALL rpin(rep_env_id,ierr)
END SUBROUTINE cp_rep_init

! *****************************************************************************
RECURSIVE SUBROUTINE cp_rep_destroy(rep_env_id,ierr) BIND(C,name="cp_c_rep_destroy")
  USE replica_methods,             ONLY: rpin => rep_env_destroy_low
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int) :: rep_env_id,ierr

  CALL rpin(rep_env_id,ierr)
END SUBROUTINE cp_rep_destroy

! *****************************************************************************
RECURSIVE SUBROUTINE cp_rep_calc_e_f(rep_env_id,calc_f,ierr) BIND(C,name="cp_c_rep_calc_e_f")
  USE replica_methods,             ONLY: rpef => rep_env_calc_e_f_low
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int) :: rep_env_id,calc_f,ierr

  CALL rpef(rep_env_id,calc_f,ierr)
END SUBROUTINE cp_rep_calc_e_f

! *****************************************************************************
RECURSIVE SUBROUTINE cp_ep_init(ep_env_id,ierr) BIND(C,name="cp_c_ep_init")
  USE ep_methods,                  ONLY: epin => ep_env_init_low
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int) :: ep_env_id, ierr

  CALL epin(ep_env_id, ierr)
END SUBROUTINE cp_ep_init

! *****************************************************************************
RECURSIVE SUBROUTINE cp_ep_calc_e_f(ep_env_id,calc_f,ierr) BIND(C,name="cp_c_ep_calc_e_f")
  USE ep_methods,                  ONLY: epef => ep_env_calc_e_f_low
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int
  IMPLICIT NONE
  INTEGER(c_int) :: ep_env_id, calc_f, ierr

  CALL epef(ep_env_id, calc_f, ierr)
END SUBROUTINE cp_ep_calc_e_f

! *****************************************************************************
SUBROUTINE cp_do_shake(f_env_id,dt,shake_tol,ierr) BIND(C,name="cp_c_do_shake")
  USE kinds,                       ONLY: dp
  USE f77_interface,               ONLY: do_shake
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_double
  IMPLICIT NONE
  INTEGER(c_int) :: f_env_id, ierr
  REAL(c_double) :: dt, shake_tol

  CALL do_shake(f_env_id,dt,shake_tol,ierr)
END SUBROUTINE cp_do_shake

! *****************************************************************************
SUBROUTINE cp_c_ext_scf_set_ptr(f_env_id,func,ierr) BIND(C,name="cp_c_ext_scf_set_ptr")
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_funptr
  USE f77_interface,               ONLY: f_env_add_defaults, f_env_rm_defaults, f_env_type 
  USE cp_error_handling,           ONLY: cp_error_type
  USE force_env_types,             ONLY: force_env_get, use_qs_force
  IMPLICIT NONE
  TYPE(c_funptr), VALUE :: func
  INTEGER(c_int) :: ierr
  INTEGER(c_int) :: f_env_id

  LOGICAL                   :: failure
  TYPE(cp_error_type)       :: error
  TYPE(f_env_type), POINTER :: f_env
  INTEGER                   :: in_use

  failure = .FALSE.
  NULLIFY (f_env)
  CALL f_env_add_defaults(f_env_id,f_env,error,failure)
  IF (.NOT. failure) THEN
     CALL force_env_get(f_env%force_env,in_use=in_use,error=error)
     IF (in_use .EQ. use_qs_force) THEN
         f_env%force_env%qs_env%ls_scf_env%ext_scf_ptr = func
     END IF
  END IF
  CALL f_env_rm_defaults(f_env,error,ierr) 
END SUBROUTINE cp_c_ext_scf_set_ptr

