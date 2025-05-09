!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE qs_tddfpt_eigensolver
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_symm,&
                                              cp_fm_trace
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pools_create_fm_vect,&
                                              fm_pools_give_back_fm_vect
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_element,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_to_string
   USE input_constants,                 ONLY: tddfpt_davidson,&
                                              tddfpt_lanczos
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_matrix_pools,                 ONLY: mpools_get
   USE qs_p_env_methods,                ONLY: p_op_l1,&
                                              p_op_l2,&
                                              p_postortho,&
                                              p_preortho
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_tddfpt_types,                 ONLY: tddfpt_env_type
   USE qs_tddfpt_utils,                 ONLY: co_initial_guess,&
                                              normalize,&
                                              reorthogonalize
#include "./base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt_eigensolver'

   PRIVATE

   PUBLIC :: eigensolver

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param p_env ...
!> \param qs_env ...
!> \param t_env ...
! **************************************************************************************************
   SUBROUTINE eigensolver(p_env, qs_env, t_env)

      TYPE(qs_p_env_type)                                :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_env_type), INTENT(INOUT)               :: t_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eigensolver'

      INTEGER                                            :: handle, n_ev, nspins, output_unit, &
                                                            restarts
      LOGICAL                                            :: do_kernel_save
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ievals
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      NULLIFY (dft_control)

      output_unit = cp_logger_get_default_io_unit()

      CALL get_qs_env(qs_env, dft_control=dft_control)
      n_ev = dft_control%tddfpt_control%n_ev
      nspins = dft_control%nspins

      ALLOCATE (ievals(n_ev))

      !---------------!
      ! initial guess !
      !---------------!
      do_kernel_save = dft_control%tddfpt_control%do_kernel
      dft_control%tddfpt_control%do_kernel = .FALSE.
      IF (output_unit > 0) THEN
         WRITE (output_unit, *) " Generating initial guess"
         WRITE (output_unit, *)
      END IF
      IF (ASSOCIATED(dft_control%tddfpt_control%lumos)) THEN
         CALL co_initial_guess(t_env%evecs, ievals, n_ev, qs_env)
      ELSE
         IF (output_unit > 0) WRITE (output_unit, *) "LUMOS are needed in TDDFPT!"
         CPABORT("")
      END IF
      DO restarts = 1, dft_control%tddfpt_control%n_restarts
         IF (iterative_solver(ievals, t_env, p_env, qs_env, ievals)) EXIT
         IF (output_unit > 0) THEN
            WRITE (output_unit, *) " Restarting"
            WRITE (output_unit, *)
         END IF
      END DO
      dft_control%tddfpt_control%do_kernel = do_kernel_save

      !-----------------!
      ! call the solver !
      !-----------------!
      IF (output_unit > 0) THEN
         WRITE (output_unit, *)
         WRITE (output_unit, *) " Doing TDDFPT calculation"
         WRITE (output_unit, *)
      END IF
      DO restarts = 1, dft_control%tddfpt_control%n_restarts
         IF (iterative_solver(ievals, t_env, p_env, qs_env, t_env%evals)) EXIT
         IF (output_unit > 0) THEN
            WRITE (output_unit, *) " Restarting"
            WRITE (output_unit, *)
         END IF
      END DO

      !---------!
      ! cleanup !
      !---------!
      DEALLOCATE (ievals)

      CALL timestop(handle)

   END SUBROUTINE eigensolver

   ! in_evals  : approximations to the eigenvalues for the preconditioner
   ! t_env     : TD-DFT environment values
   ! p_env     : perturbation environment values
   ! qs_env    : general Quickstep environment values
   ! out_evals : the resulting eigenvalues
   ! error     : used for error handling
   !
   ! res       : the function will return wheter the eigenvalues are converged or not

! **************************************************************************************************
!> \brief ...
!> \param in_evals ...
!> \param t_env ...
!> \param p_env ...
!> \param qs_env ...
!> \param out_evals ...
!> \return ...
! **************************************************************************************************
   FUNCTION iterative_solver(in_evals, &
                             t_env, p_env, qs_env, &
                             out_evals) RESULT(res)

      REAL(KIND=dp), DIMENSION(:)                        :: in_evals
      TYPE(tddfpt_env_type), INTENT(INOUT)               :: t_env
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(kind=dp), DIMENSION(:), OPTIONAL              :: out_evals
      LOGICAL                                            :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'iterative_solver', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: mode
      INTEGER                                            :: col, handle, i, iev, iter, j, k, &
                                                            max_krylovspace_dim, max_kv, n_ev, &
                                                            n_kv, nspins, output_unit, row, spin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: must_improve
      REAL(dp)                                           :: Atilde_ij, convergence, tmp, tmp2
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_difference, evals_tmp
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
      TYPE(cp_fm_struct_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: kv_fm_struct
      TYPE(cp_fm_struct_type), POINTER                   :: tilde_fm_struct
      TYPE(cp_fm_type)                                   :: Atilde, Us
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: R, X
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: Ab, b, Sb
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(tddfpt_control_type), POINTER                 :: tddfpt_control

      res = .FALSE.

      CALL timeset(routineN, handle)

      NULLIFY (ao_mo_fm_pools, tddfpt_control, &
               tilde_fm_struct, matrix_s, dft_control, &
               para_env, blacs_env)

      CALL get_qs_env(qs_env, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      tddfpt_control => dft_control%tddfpt_control
      output_unit = cp_logger_get_default_io_unit()
      n_ev = tddfpt_control%n_ev
      nspins = dft_control%nspins

      IF (dft_control%tddfpt_control%diag_method == tddfpt_lanczos) THEN
         mode = 'L'
      ELSE IF (dft_control%tddfpt_control%diag_method == tddfpt_davidson) THEN
         mode = 'D'
      END IF

      !-----------------------------------------!
      ! determine the size of the problem       !
      ! and how many krylov space vetors to use !
      !-----------------------------------------!
      max_krylovspace_dim = SUM(p_env%n_ao(1:nspins)*p_env%n_mo(1:nspins))
      max_kv = tddfpt_control%max_kv
      IF (max_krylovspace_dim <= max_kv) THEN
         max_kv = max_krylovspace_dim
         IF (output_unit > 0) THEN
            WRITE (output_unit, *) "  Setting the maximum number of krylov vectors to ", max_kv, "!"
         END IF
      END IF

      !----------------------!
      ! allocate the vectors !
      !----------------------!
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
      CALL fm_pools_create_fm_vect(ao_mo_fm_pools, X, name=routineP//":X")
      CALL fm_pools_create_fm_vect(ao_mo_fm_pools, R, name=routineP//":R")

      ALLOCATE (evals_difference(n_ev))

      ALLOCATE (must_improve(n_ev))

      ALLOCATE (evals(max_kv, 0:max_kv))
      ALLOCATE (evals_tmp(max_kv))

      ALLOCATE (b(max_kv, nspins), Ab(max_kv, nspins), &
                Sb(max_kv, nspins))

      ALLOCATE (kv_fm_struct(nspins))

      DO spin = 1, nspins
         CALL cp_fm_struct_create(kv_fm_struct(spin)%struct, para_env, blacs_env, &
                                  p_env%n_ao(spin), p_env%n_mo(spin))
      END DO

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(2X,A,T69,A)') &
            "nvec", "Convergence"
         WRITE (output_unit, '(2X,A)') &
            "-----------------------------------------------------------------------------"
      END IF

      iter = 1
      k = 0
      n_kv = n_ev
      iteration: DO

         CALL allocate_krylov_vectors(b, "b-", k + 1, n_kv, nspins, kv_fm_struct)
         CALL allocate_krylov_vectors(Ab, "Ab-", k + 1, n_kv, nspins, kv_fm_struct)
         CALL allocate_krylov_vectors(Sb, "Sb-", k + 1, n_kv, nspins, kv_fm_struct)

         DO i = 1, n_kv
            k = k + 1

            IF (k <= SIZE(t_env%evecs, 1)) THEN ! the first iteration

               ! take the initial guess
               DO spin = 1, nspins
                  CALL cp_fm_to_fm(t_env%evecs(k, spin), b(k, spin))
               END DO

            ELSE

               ! create a new vector
               IF (mode == 'L') THEN

                  DO spin = 1, nspins
                     IF (tddfpt_control%invert_S) THEN
                        CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                                        1.0_dp, t_env%invS(spin), Ab(k - 1, spin), &
                                        0.0_dp, b(k, spin))
                     ELSE
                        CALL cp_fm_to_fm(Ab(k - 1, spin), b(k, spin))
                     END IF
                  END DO

               ELSE IF (mode == 'D') THEN

                  iev = must_improve(i)
                  ! create the new davidson vector
                  DO spin = 1, nspins

                     CALL cp_fm_set_all(R(spin), 0.0_dp)
                     DO j = 1, k - i
                        CALL cp_fm_to_fm(Ab(j, spin), X(spin))
                        CALL cp_fm_scale_and_add(1.0_dp, X(spin), &
                                                 -evals(iev, iter - 1), Sb(j, spin))
                        CALL cp_fm_get_element(Us, j, iev, tmp)
                        CALL cp_fm_scale_and_add(1.0_dp, R(spin), &
                                                 tmp, X(spin))
                     END DO

                     IF (tddfpt_control%invert_S) THEN
                        CALL cp_fm_symm('L', 'U', p_env%n_ao(spin), p_env%n_mo(spin), &
                                        1.0_dp, t_env%invS(spin), R(spin), &
                                        0.0_dp, X(spin))
                     ELSE
                        CALL cp_fm_to_fm(R(spin), X(spin))
                     END IF

                     !----------------!
                     ! preconditioner !
                     !----------------!
                     IF (dft_control%tddfpt_control%precond) THEN
                        DO col = 1, p_env%n_mo(spin)
                           IF (col <= n_ev) THEN
                              tmp2 = ABS(evals(iev, iter - 1) - in_evals(col))
                           ELSE
                              tmp2 = ABS(evals(iev, iter - 1) - (in_evals(n_ev) + 10.0_dp))
                           END IF
                           ! protect against division by 0 by a introducing a cutoff.
                           tmp2 = MAX(tmp2, 100*EPSILON(1.0_dp))
                           DO row = 1, p_env%n_ao(spin)
                              CALL cp_fm_get_element(X(spin), row, col, tmp)
                              CALL cp_fm_set_element(b(k, spin), row, col, tmp/tmp2)
                           END DO
                        END DO
                     ELSE
                        CALL cp_fm_to_fm(X(spin), b(k, spin))
                     END IF

                  END DO

               ELSE
                  IF (output_unit > 0) WRITE (output_unit, *) "unknown mode"
                  CPABORT("")
               END IF

            END IF

            CALL p_preortho(p_env, qs_env, b(k, :))
            DO j = 1, tddfpt_control%n_reortho
               CALL reorthogonalize(b(k, :), b, Sb, R, k - 1) ! R is temp
            END DO
            CALL normalize(b(k, :), R, matrix_s) ! R is temp
            DO spin = 1, nspins
               CALL cp_fm_to_fm(b(k, spin), X(spin))
            END DO
            CALL apply_op(X, Ab(k, :), p_env, qs_env, &
                          dft_control%tddfpt_control%do_kernel)
            CALL p_postortho(p_env, qs_env, Ab(k, :))
            DO spin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, &
                                            b(k, spin), &
                                            Sb(k, spin), &
                                            p_env%n_mo(spin))
            END DO
         END DO

         !--------------------------------------------!
         ! deallocate memory for the reduced matrices !
         !--------------------------------------------!
         CALL cp_fm_release(Atilde)
         CALL cp_fm_release(Us)
         IF (ASSOCIATED(tilde_fm_struct)) CALL cp_fm_struct_release(tilde_fm_struct)

         !------------------------------------------!
         ! allocate memory for the reduced matrices !
         !------------------------------------------!
         CALL cp_fm_struct_create(tilde_fm_struct, para_env, blacs_env, k, k)
         CALL cp_fm_create(Atilde, &
                           tilde_fm_struct, &
                           routineP//"Atilde")
         CALL cp_fm_create(Us, &
                           tilde_fm_struct, &
                           routineP//"Us")

         !---------------------------------------!
         ! calc the matrix Atilde = transp(b)*Ab !
         !---------------------------------------!
         DO i = 1, k
            DO j = 1, k
               Atilde_ij = 0.0_dp
               DO spin = 1, nspins
                  CALL cp_fm_trace(b(i, spin), Ab(j, spin), tmp)
                  Atilde_ij = Atilde_ij + tmp
               END DO
               CALL cp_fm_set_element(Atilde, i, j, Atilde_ij)
            END DO
         END DO

         !--------------------!
         ! diagonalize Atilde !
         !--------------------!
         evals_tmp(:) = evals(:, iter)
         CALL cp_fm_syevd(Atilde, Us, evals_tmp(:))
         evals(:, iter) = evals_tmp(:)

         !-------------------!
         ! check convergence !
         !-------------------!
         evals_difference = 1.0_dp
         IF (iter /= 1) THEN

            evals_difference(:) = ABS((evals(1:n_ev, iter - 1) - evals(1:n_ev, iter)))
            ! For debugging
            IF (output_unit > 0) THEN
               WRITE (output_unit, *)
               DO i = 1, n_ev
                  WRITE (output_unit, '(2X,F10.7,T69,ES11.4)') evals(i, iter)*evolt, evals_difference(i)
               END DO
               WRITE (output_unit, *)
            END IF

            convergence = MAXVAL(evals_difference)
            IF (output_unit > 0) WRITE (output_unit, '(2X,I4,T69,ES11.4)') k, convergence

            IF (convergence < tddfpt_control%tolerance) THEN
               res = .TRUE.
               EXIT iteration
            END IF
         END IF

         IF (mode == 'L') THEN
            n_kv = 1
         ELSE
            must_improve = 0
            DO i = 1, n_ev
               IF (evals_difference(i) > tddfpt_control%tolerance) must_improve(i) = 1
            END DO
!! Set must_improve to 1 if all the vectors should
!! be updated in one iteration.
!!          must_improve = 1
            n_kv = SUM(must_improve)
            j = 1
            DO i = 1, n_ev
               IF (must_improve(i) == 1) THEN
                  must_improve(j) = i
                  j = j + 1
               END IF
            END DO
         END IF

         IF (k + n_kv > max_kv) EXIT iteration

         iter = iter + 1

      END DO iteration

      IF (PRESENT(out_evals)) THEN
         out_evals(1:n_ev) = evals(1:n_ev, iter)
      END IF

      DO spin = 1, nspins
         DO j = 1, n_ev
            CALL cp_fm_set_all(t_env%evecs(j, spin), 0.0_dp)
            DO i = 1, k
               CALL cp_fm_get_element(Us, i, j, tmp)
               CALL cp_fm_scale_and_add(1.0_dp, t_env%evecs(j, spin), &
                                        tmp, b(i, spin))
            END DO
         END DO
      END DO

      !----------!
      ! clean up !
      !----------!
      CALL cp_fm_release(Atilde)
      CALL cp_fm_release(Us)
      IF (ASSOCIATED(tilde_fm_struct)) CALL cp_fm_struct_release(tilde_fm_struct)
      CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, X)
      CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, R)
      DO spin = 1, nspins
         CALL cp_fm_struct_release(kv_fm_struct(spin)%struct)
      END DO
      CALL cp_fm_release(b)
      CALL cp_fm_release(Ab)
      CALL cp_fm_release(Sb)
      DEALLOCATE (evals, evals_tmp, evals_difference, must_improve, kv_fm_struct)

      CALL timestop(handle)

   END FUNCTION iterative_solver

   ! X        : the vector on which to apply the op
   ! R        : the result
   ! t_env    : td-dft environment (mainly control information)
   ! p_env    : perturbation environment (variables)
   !            both of these carry info for the tddfpt calculation
   ! qs_env   : info about a quickstep ground state calculation

! **************************************************************************************************
!> \brief ...
!> \param X ...
!> \param R ...
!> \param p_env ...
!> \param qs_env ...
!> \param do_kernel ...
! **************************************************************************************************
   SUBROUTINE apply_op(X, R, p_env, qs_env, do_kernel)

      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: X, R
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: do_kernel

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apply_op'

      INTEGER                                            :: handle, nspins, spin
      INTEGER, SAVE                                      :: counter = 0
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (dft_control)

      CALL timeset(routineN, handle)

      counter = counter + 1
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      !------------!
      ! R = HX-SXL !
      !------------!
      CALL p_op_l1(p_env, qs_env, X, R) ! acts on both spins, result in R

      !-----------------!
      ! calc P1 and     !
      ! R = R + K(P1)*C !
      !-----------------!
      IF (do_kernel) THEN
         DO spin = 1, nspins
            CALL dbcsr_set(p_env%p1(spin)%matrix, 0.0_dp) ! optimize?
            CALL cp_dbcsr_plus_fm_fm_t(p_env%p1(spin)%matrix, &
                                       matrix_v=p_env%psi0d(spin), &
                                       matrix_g=X(spin), &
                                       ncol=p_env%n_mo(spin), &
                                       symmetry_mode=1)
         END DO
         DO spin = 1, nspins
            CALL cp_fm_set_all(X(spin), 0.0_dp)
         END DO
         CALL p_op_l2(p_env, qs_env, p_env%p1, X, &
                      alpha=1.0_dp, beta=0.0_dp) ! X = beta*X + alpha*K(P1)*C
         DO spin = 1, nspins
            CALL cp_fm_scale_and_add(1.0_dp, R(spin), &
                                     1.0_dp, X(spin)) ! add X to R
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_op

! **************************************************************************************************
!> \brief ...
!> \param vectors ...
!> \param vectors_name ...
!> \param startv ...
!> \param n_v ...
!> \param nspins ...
!> \param fm_struct ...
! **************************************************************************************************
   SUBROUTINE allocate_krylov_vectors(vectors, vectors_name, &
                                      startv, n_v, nspins, fm_struct)

      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: vectors
      CHARACTER(LEN=*), INTENT(IN)                       :: vectors_name
      INTEGER, INTENT(IN)                                :: startv, n_v, nspins
      TYPE(cp_fm_struct_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: fm_struct

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_krylov_vectors', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: mat_name
      INTEGER                                            :: index, spin

      DO spin = 1, nspins
         DO index = startv, startv + n_v - 1
            mat_name = routineP//vectors_name//TRIM(cp_to_string(index)) &
                       //","//TRIM(cp_to_string(spin))
            CALL cp_fm_create(vectors(index, spin), &
                              fm_struct(spin)%struct, mat_name)
            IF (.NOT. ASSOCIATED(vectors(index, spin)%matrix_struct)) &
               CPABORT("Could not allocate "//TRIM(mat_name)//".")
         END DO
      END DO

   END SUBROUTINE allocate_krylov_vectors

END MODULE qs_tddfpt_eigensolver
