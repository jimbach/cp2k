!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Update a QM/MM calculations with ADAPTIVE force mixing
!> \par History
!>      5.2004 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE adqmmmx_update
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_release,&
                                              cp_subsys_type
   USE particle_types,                  ONLY: particle_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_restart_force_eval,        ONLY: update_force_eval
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_val_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals3,&
                                              section_vals_release,&
                                              section_vals_type
   USE qmmm_create,                     ONLY: qmmm_env_create
   USE qmmm_types,                      ONLY: qmmm_env_get
   USE adqmmmx_types,                   ONLY: adqmmmx_env_release,&
                                              adqmmmx_env_get,&
                                              adqmmmx_env_type
   USE adqmmmx_util,                    ONLY: setup_adqmmmx_sections,&
                                              update_adqmmmx_labels,&
                                              update_adqmmmx_values,&
!                                             print_adqmmmx_output,&
                                              copy_adqmmmx_env,&
                                              g_of_lambda
                                              !old_update_adqmmmx_labels,&
   USE kinds,                           ONLY: dp
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'adqmmmx_update'

   PUBLIC :: adqmmmx_update_force_env, adqmmmx_vv_first, adqmmmx_vv_second

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \note Updated names (qmmmx->adqmmmx)
!> \param root_section ...
!> \par History
!>      12.2021 created
! **************************************************************************************************
   SUBROUTINE adqmmmx_update_force_env(force_env, root_section, dt)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(cp_subsys_type), POINTER                      :: subsys,subsys_new!, subsys_oldsys,&
                                                            !subsys_newsys
!     TYPE(particle_type), DIMENSION(:), POINTER         :: particles_oldsys, particles_newsys
      TYPE(distribution_1d_type), POINTER                :: local_particles, new_local_particles
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds, new_atomic_kinds
      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_update_force_env', routineP = moduleN//':'//routineN
      LOGICAL                                            :: file_exists,labels_changed, adqmmmx_active,&
                                                            lam_int_disc
      REAL(dp)                                           :: dt, switch_t, Lam_Error, swprecision
      REAL(dp), POINTER                                  :: lambda, vlambda
      REAL(KIND=dp), PARAMETER                           :: pi = 3.14159265358979323846264338_dp
!      TYPE(adqmmmx_env_type), POINTER                    :: new_adqmmmx_env
! JB MODS 28 03
      TYPE(adqmmmx_env_type)                             :: new_adqmmmx_env
      TYPE(section_vals_type), POINTER                   :: qmmm_oldsys_section, &
                                                            qmmm_newsys_section, &
                                                            adqmmmx_section, qmmm_section, &
                                                            subsys_section
      INTEGER                                            :: max_n_sw

      IF (.NOT. ASSOCIATED(force_env)) RETURN
      IF (ASSOCIATED(force_env%adqmmmx_env)) THEN

      
      swprecision=0.01d0
      WRITE(*,*) "this is a test"

      INQUIRE(FILE = "update_log.out", EXIST=file_exists)
      IF (file_exists) THEN
         OPEN (UNIT=1337, FILE = "update_log.out", STATUS = "old", POSITION = "append", ACTION = "write")
      ELSE
         OPEN (UNIT=1337, FILE = "update_log.out", STATUS = "new", ACTION = "write")
      END IF
     
!     IF ( force_env%adqmmmx_env%oldsys .EQ. force_env%adqmmmx_env%newsys) THEN
!        WRITE (1337,*) "OLDSYS=NEWSYS ==> SKIP update routine"
!        force_env%adqmmmx_env%switch_active=.FALSE.
!        RETURN
      IF ( force_env%adqmmmx_env%dpot == 0.0) THEN
         WRITE (1337,*) "DPOT = 0 ==> SKIP update routine and set switch_active to False"
         force_env%adqmmmx_env%switch_active=.FALSE.
         force_env%adqmmmx_env%lambda_adqmmmx=1.d0
         force_env%adqmmmx_env%lambda_velocity=0.d0
!        RETURN
      END IF 
      ELSE IF (.NOT. ASSOCIATED(force_env%adqmmmx_env)) THEN
              RETURN
      END IF
      NULLIFY(lambda)
      lambda => force_env%adqmmmx_env%lambda_adqmmmx
      NULLIFY(vlambda)
      vlambda => force_env%adqmmmx_env%lambda_velocity

!     !MN Print out Output one time when t=0
!     IF (force_env%adqmmmx_env%prop_time==0) CALL print_adqmmmx_output(force_env%adqmmmx_env)
      !MN Propagate lambda timer
      force_env%adqmmmx_env%prop_time = force_env%adqmmmx_env%prop_time + dt
      !MN Check if discrete or analytical calculation of Lambda is needed
      qmmm_section => section_vals_get_subs_vals(force_env%force_env_section, "QMMM", can_return_null=.TRUE.)
      adqmmmx_section => section_vals_get_subs_vals(qmmm_section, "ADAPTIVE_QMMM", can_return_null=.TRUE.) 
      CALL section_vals_val_get(adqmmmx_section, "DISCRETE_INTEGRATION", l_val=lam_int_disc)


      !MN Get adaptive section input values
      ! check everything for not null, because sometimes (e.g. metadynamics in parallel) it happens
!     IF (.NOT. ASSOCIATED(force_env%force_env_section)) RETURN

      !MN Check if MAX_N_SWITCH has been reached
      CALL section_vals_val_get(adqmmmx_section, "MAX_N_SWITCH",i_val=max_n_sw)
      IF ((force_env%adqmmmx_env%c_n_sw>max_n_sw) .AND. (max_n_sw/=-1)) THEN
         WRITE(1337,*) "NO SWITCHING BECAUSE MAXIMUM WAS REACHED", max_n_sw,&
                        force_env%adqmmmx_env%c_n_sw
         RETURN
      END IF



      !JB Insert Lambda check
      IF ((0 <= lambda) .AND. ((lambda-1.d0) >= swprecision) &
             .AND. (force_env%adqmmmx_env%switch_active)) THEN
         WRITE (1337,*) "UPDATE_FORCE_ENV: Switching still active..."
         IF (force_env%adqmmmx_env%lambda_mass == 0.d0) THEN
            CALL section_vals_val_get(adqmmmx_section, "SWITCH_T", r_val=switch_t)
            !MN This definition only covers a quarter of a transition, multiply by 4!
            WRITE (1337,*) "UPDATE_FORCE_ENV: First and only calculation of lambda_mass"
            IF (lam_int_disc) THEN
               switch_t = 4.0*(switch_t) - dt
               force_env%adqmmmx_env%lambda_mass = (2*DABS(force_env%adqmmmx_env%dpot)*switch_t**2)/&
                                                (4*pi**2)
               ELSE
               switch_t = 4.0*(switch_t)
               force_env%adqmmmx_env%lambda_mass = (2*DABS(force_env%adqmmmx_env%dpot)*switch_t**2)/&
                                                (4*pi**2)
            END IF
            WRITE (1337,*) "MASS: ", force_env%adqmmmx_env%lambda_mass
         END IF
      !JB modify check
      ELSE IF ((lambda-1.d0) <= swprecision .OR. (.NOT. force_env%adqmmmx_env%switch_active)) THEN
         WRITE (1337,*) "UPDATE_FORCE_ENV: Switching complete, start new update cycle..."
         WRITE (1337,*) "Switch Active:", force_env%adqmmmx_env%switch_active

         IF (force_env%adqmmmx_env%energy_add) THEN
            !JB Calculate Error function for total lambda energy only once!
            Lam_Error = force_env%adqmmmx_env%dpot*(lambda - 1.d0)**2
            WRITE (1337,*) "g of dlambda", g_of_lambda(lambda),(lambda-1.d0)**2
            WRITE (1337,*) "Total E_l is updated. OLD: ", force_env%adqmmmx_env%Tot_Elambda,&
                           "Added E: ", force_env%adqmmmx_env%Elambda,&
                           "Added dE:", Lam_Error
            force_env%adqmmmx_env%Tot_Elambda     = force_env%adqmmmx_env%Tot_Elambda+&
                                          force_env%adqmmmx_env%Elambda + Lam_Error+&
                                          force_env%adqmmmx_env%Elamcorr
            force_env%adqmmmx_env%Elamcorr        = 0.d0
            force_env%adqmmmx_env%Elambda         = 0.d0
            !Added JB 16.04.25
            force_env%adqmmmx_env%lambda_force    = 0.d0           
            lambda                                = 1.d0
            force_env%adqmmmx_env%lambda_velocity = 0.d0
            force_env%adqmmmx_env%energy_add      =.FALSE.

            force_env%adqmmmx_env%switch_active   = .FALSE.

         !JB add a check for lambda > 1
         ELSE IF ((lambda > 1)) THEN
            CPABORT("Lambda > 1 overshot, choose a longer switching time")
         END IF
         ELSE
            CPABORT("LAMBDA IS LOWER THAN 0 OR NaN, SOMETHING WENT WRONG!")
         END IF
      
      IF (.NOT. ASSOCIATED(adqmmmx_section)) RETURN
      CALL section_vals_get(adqmmmx_section, explicit=adqmmmx_active)
      WRITE (1337,*) "adqmmmx_active= ", adqmmmx_active
      IF (.NOT. ASSOCIATED(force_env%adqmmmx_env)) CPABORT("force_env%adqmmmx_env not associated")

      IF (.NOT. ASSOCIATED(qmmm_section) .OR. .NOT. ASSOCIATED(adqmmmx_section)) THEN
         WRITE (1337,*) "UPDATE_FORCE_ENV: PROBLEM WITH QMMM section", ASSOCIATED(qmmm_section), ASSOCIATED(adqmmmx_section)
         CPABORT("Problem with QMMM section")
      END IF


      IF ( .NOT. force_env%adqmmmx_env%switch_active ) THEN 
         WRITE (1337,*) "Starting label update"
         labels_changed = .TRUE. !MN Test if logical variable is given to routine correctly
         CALL force_env_get(force_env, subsys=subsys)
         
         !MN Create new_adqmmmx_env
         !JB changed 28.4.24
         !ALLOCATE (new_adqmmmx_env)
         !MN Create sections for oldsys using the old updated indices/labels
         CALL setup_adqmmmx_sections(subsys, qmmm_section, qmmm_oldsys_section)

         subsys_section => section_vals_get_subs_vals(force_env%force_env_section, "SUBSYS")
         ![ADAPT] no sure about use_motion_section
         !MN Create qmmm_env with oldsys_section
         !JB ADD 4 02 2025
         !DEALLOCATE (new_adqmmmx_env%oldsys)
         ALLOCATE (new_adqmmmx_env%oldsys)
         CALL qmmm_env_create(new_adqmmmx_env%oldsys, &
                              force_env%root_section, force_env%para_env, force_env%globenv, &
                              force_env%force_env_section, qmmm_oldsys_section, subsys_section, use_motion_section=.TRUE., &
                              prev_subsys=subsys, ignore_outside_box=.TRUE.)

         ! [NB] need to copy wiener process data, since it's not recreated when
         ! fist subsys is recreated by qmmm_env_create
         CALL qmmm_env_get(force_env%adqmmmx_env%newsys, subsys=subsys)
         CALL qmmm_env_get(new_adqmmmx_env%oldsys, subsys=subsys)
         CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, local_particles=local_particles)
         !MN Copy wiener process for oldsys
         CALL qmmm_env_get(new_adqmmmx_env%oldsys, subsys=subsys_new)
         CALL cp_subsys_get(subsys_new, atomic_kinds=new_atomic_kinds, local_particles=new_local_particles)
         IF (ASSOCIATED(local_particles%local_particle_set)) THEN
            CALL copy_wiener_process(atomic_kinds, local_particles, new_atomic_kinds, new_local_particles)
         ENDIF

         !MN Call update routine
         CALL update_adqmmmx_labels(subsys,qmmm_section,labels_changed=labels_changed,&
               include_nc=force_env%adqmmmx_env%include_nc,c_n_sw=force_env%adqmmmx_env%c_n_sw)
         WRITE (1337,*) "Labels changed?", labels_changed

         IF (.NOT. labels_changed) THEN
               WRITE (1337,*) "Labels have not changed or MAX_N_SWITCH reached"
               force_env%adqmmmx_env%switch_active=.FALSE.
               ! JB MOD 4.4.25
               CALL adqmmmx_env_release(new_adqmmmx_env)
               RETURN               
         END IF

         !MN Check if Indices and labels from new_*** are correctly handled
         CPWARN("Adaptive QMMM labels changed, rebuilding QM/MM calculations! ")
         CALL update_force_eval(force_env, root_section, .FALSE.) !MN F=No binary restart file

         WRITE(1337,*) "Update_force_eval called. Now create new adqmmmx environment"

         !MN Create sections for newsys using the new updated indices/labels

         CALL setup_adqmmmx_sections(subsys, qmmm_section, qmmm_newsys_section)
         !JB ADD 4 02 2025
         !DEALLOCATE (new_adqmmmx_env%newsys)
         ALLOCATE (new_adqmmmx_env%newsys)
         CALL qmmm_env_create(new_adqmmmx_env%newsys, &
                              force_env%root_section, force_env%para_env, force_env%globenv, &
                              force_env%force_env_section, qmmm_newsys_section, subsys_section, use_motion_section=.TRUE., &
                              prev_subsys=subsys, ignore_outside_box=.TRUE.)
         WRITE (1337,*) "NEWSYS ASSOCIATED", ASSOCIATED(new_adqmmmx_env%newsys), new_adqmmmx_env%newsys%ref_count

         

         ! [NB] need to copy wiener process data, since it's not recreated when
         ! fist subsys is recreated by qmmm_env_create
         CALL qmmm_env_get(force_env%adqmmmx_env%newsys, subsys=subsys_new)
         CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, local_particles=local_particles)
         CALL qmmm_env_get(new_adqmmmx_env%newsys, subsys=subsys_new)
         CALL cp_subsys_get(subsys_new, atomic_kinds=new_atomic_kinds, local_particles=new_local_particles)
         IF (ASSOCIATED(local_particles%local_particle_set)) THEN
            CALL copy_wiener_process(atomic_kinds, local_particles, new_atomic_kinds, new_local_particles)
         ENDIF

         CALL section_vals_release(qmmm_oldsys_section)
         CALL section_vals_release(qmmm_newsys_section)

         !MN Copy adaptive switching parameters and release old environemnt
         
         CALL copy_adqmmmx_env(new_adqmmmx_env, force_env%adqmmmx_env)
         CALL adqmmmx_env_release(force_env%adqmmmx_env)
         !force_env%adqmmmx_env => new_adqmmmx_env
         ! JB MOD 28 03
         force_env%adqmmmx_env = new_adqmmmx_env

      END IF

      !JB Write lambda propagation output
      IF ( force_env%adqmmmx_env%switch_active ) THEN
         !MN adqmmmx_section empty?
         CALL section_vals_val_get(adqmmmx_section, "SWITCH_T", r_val=switch_t)
         !MN This definition only covers a quarter of a transition, multiply by 4!
         switch_t = 4.0*switch_t
         WRITE(1337,*) "ADQMMMX_SECTION ASS:", ASSOCIATED(adqmmmx_section), "SWITCH_T:", switch_t/4
         force_env%adqmmmx_env%do_discrete = lam_int_disc
         WRITE(1337,*) "LAM_INT_DISC: ", lam_int_disc,force_env%adqmmmx_env%do_discrete
         WRITE (1337,*) "Elambda = ", force_env%adqmmmx_env%Elambda, "TOT_Elambda = ", force_env%adqmmmx_env%Tot_Elambda

      ELSE
         WRITE (1337,*) "First Lambda-Step 'skipped'"
         force_env%adqmmmx_env%switch_active = .TRUE.
      END IF

      !MN adqmmmx values in restart file
      CALL update_adqmmmx_values(adqmmmx_section, force_env%adqmmmx_env, 1337)

      CLOSE (1337)


   END SUBROUTINE adqmmmx_update_force_env


! **************************************************************************************************
!> \brief ...
!> \param from_local_particle_kinds ...
!> \param from_local_particles ...
!> \param to_local_particle_kinds ...
!> \param to_local_particles ...
! **************************************************************************************************
   SUBROUTINE lambda_propagation(adqmmmx_env, switch_t, dt, output_unit)
      !TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      TYPE(adqmmmx_env_type), INTENT(INOUT)                 :: adqmmmx_env
      CHARACTER(len=*), PARAMETER :: routineN = 'lambda_propagation', routineP = moduleN//':'//routineN
      REAL(dp)                                           :: lambda, prop_time, dpot, vlambda, mlambda
      REAL(dp)                                           :: dt, omega, omega1, switch_t, A, B
      INTEGER                                            :: output_unit
      REAL(KIND=dp), PARAMETER                           :: pi = 3.14159265358979323846264338_dp


      lambda = adqmmmx_env%lambda_adqmmmx
      prop_time = adqmmmx_env%prop_time
      vlambda = adqmmmx_env%lambda_velocity
      dpot = adqmmmx_env%dpot
      mlambda = adqmmmx_env%lambda_mass
      !MN Omega is variable
      omega  = (2*DABS(dpot))/(mlambda) 
      WRITE (output_unit,*) "OMEGA**2", omega
      WRITE (output_unit,*) "OLD LAMBDA", lambda
      omega  = DSQRT(omega) 
      WRITE (output_unit,*) "OMEGA", omega
      omega1 = (2*pi)/(switch_t) !MN NaN Hack if omega is too small; Talk with Nikos about this
      adqmmmx_env%omega = omega
      
      IF (omega/=0.0) THEN
         A = vlambda / omega
!     ELSE
!        A = 0
      END IF

      B = lambda - 1
      
!     print *, "Old Lambda", lambda, "Old VLAMBDA", vlambda
      
      !MN At first we used prop_time here, but dt is correct when considering
      !Nikos program!
      !MN Propagate only if dpot/= 0
      IF ( mlambda/=0.0) THEN
          lambda  = 1 + A*DSIN(omega*dt) + B*DCOS(omega*dt)
          vlambda = omega*A*DCOS(omega*dt) - omega*B*DSIN(omega*dt)
      ELSE
         print *,"mlambda is 0, skip propagation"
      END IF

      adqmmmx_env%lambda_adqmmmx   = lambda
      adqmmmx_env%lambda_velocity  = vlambda

      WRITE (output_unit,*) "DPOT", dpot, "SWITCH_T", switch_t, "DT", dt
      WRITE (output_unit,*) "MASS", mlambda, "Omega", omega
      WRITE (output_unit,*) "A", A, "B", B, "New Lambda", lambda, "New VLAMBDA", vlambda
      WRITE (output_unit,*) "l= ", lambda, "v= ", vlambda, "t= ", prop_time, "w= ", omega
      WRITE (output_unit,*) "OMEGA CHECK"
      WRITE (output_unit,*) "mass calc:", omega, "switch_t calc:", omega1
   END SUBROUTINE
! **************************************************************************************************
!> \brief Integrate lambda full and lambda velocity half a step
!> \param from_local_particle_kinds ...
!> \param from_local_particles ...
!> \param to_local_particle_kinds ...
!> \param to_local_particles ...
! **************************************************************************************************
   SUBROUTINE lambda_propagation_discrete(adqmmmx_env, dt, output_unit)
      !TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      TYPE(adqmmmx_env_type), INTENT(INOUT)                 :: adqmmmx_env
      CHARACTER(len=*), PARAMETER :: routineN = 'lambda_propagation_discrete', routineP = moduleN//':'//routineN
      REAL(dp)                                           :: dt, test_F
      INTEGER                                            :: output_unit
     
      WRITE (output_unit,*) "Discrete Lambda Propagation (first step)" 
      test_F = -2*adqmmmx_env%dpot*(adqmmmx_env%lambda_adqmmmx - 1)
      adqmmmx_env%lambda_velocity = adqmmmx_env%lambda_velocity + &
                           0.5*dt*(adqmmmx_env%lambda_force/adqmmmx_env%lambda_mass) 
      adqmmmx_env%lambda_adqmmmx = adqmmmx_env%lambda_adqmmmx + dt*adqmmmx_env%lambda_velocity
      WRITE (output_unit,*) "Half V_Lam: ", adqmmmx_env%lambda_velocity
      WRITE (output_unit,*) "Acc.: ", adqmmmx_env%lambda_force/adqmmmx_env%lambda_mass


   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param from_local_particle_kinds ...
!> \param from_local_particles ...
!> \param to_local_particle_kinds ...
!> \param to_local_particles ...
! **************************************************************************************************
   SUBROUTINE adqmmmx_vv_first(force_env, dt, output_unit)
      TYPE(force_env_type), POINTER                    :: force_env
      REAL(dp)                                         :: dt
      INTEGER                                          :: output_unit

      IF (force_env%adqmmmx_env%lambda_adqmmmx==0.0) THEN
         WRITE(output_unit,*) "Skipped first vv_second"
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(force_env%adqmmmx_env) ) RETURN
      WRITE (output_unit,*) "Discrete Lambda Propagation (second step)",&
                             force_env%adqmmmx_env%do_discrete,&
                             force_env%adqmmmx_env%energy_add

      
      force_env%adqmmmx_env%lambda_velocity = force_env%adqmmmx_env%lambda_velocity + &
                           0.5d0*dt*(force_env%adqmmmx_env%lambda_force/force_env%adqmmmx_env%lambda_mass) 
      force_env%adqmmmx_env%lambda_adqmmmx = force_env%adqmmmx_env%lambda_adqmmmx + dt*force_env%adqmmmx_env%lambda_velocity
      
      WRITE(output_unit,*) "Full V_Lam: ", force_env%adqmmmx_env%lambda_velocity
      WRITE (output_unit,*) "New F_Lam: ", force_env%adqmmmx_env%lambda_force
      
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param from_local_particle_kinds ...
!> \param from_local_particles ...
!> \param to_local_particle_kinds ...
!> \param to_local_particles ...
! **************************************************************************************************
   SUBROUTINE adqmmmx_vv_second(force_env, dt, output_unit)
      TYPE(force_env_type), POINTER                    :: force_env
      REAL(dp)                                         :: dt
      INTEGER                                          :: output_unit

      IF (force_env%adqmmmx_env%lambda_adqmmmx==0.0) THEN
         WRITE(output_unit,*) "Skipped first vv_second"
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(force_env%adqmmmx_env) ) RETURN
      WRITE (output_unit,*) "Discrete Lambda Propagation (second step)",&
                             force_env%adqmmmx_env%do_discrete,&
                             force_env%adqmmmx_env%energy_add


      force_env%adqmmmx_env%lambda_velocity = force_env%adqmmmx_env%lambda_velocity +&
      0.5*dt*(force_env%adqmmmx_env%lambda_force/force_env%adqmmmx_env%lambda_mass) 
      WRITE(output_unit,*) "Full V_Lam: ", force_env%adqmmmx_env%lambda_velocity
      WRITE (output_unit,*) "New F_Lam: ", force_env%adqmmmx_env%lambda_force
      
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param from_local_particle_kinds ...
!> \param from_local_particles ...
!> \param to_local_particle_kinds ...
!> \param to_local_particles ...
! **************************************************************************************************
   SUBROUTINE copy_wiener_process(from_local_particle_kinds, from_local_particles, &
                                  to_local_particle_kinds, to_local_particles)
      TYPE(atomic_kind_list_type), POINTER               :: from_local_particle_kinds
      TYPE(distribution_1d_type), POINTER                :: from_local_particles
      TYPE(atomic_kind_list_type), POINTER               :: to_local_particle_kinds
      TYPE(distribution_1d_type), POINTER                :: to_local_particles

      CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_wiener_process', &
         routineP = moduleN//':'//routineN

      INTEGER :: from_iparticle_kind, from_iparticle_local(1), from_nparticle_kind, &
         from_nparticle_local, handle, to_iparticle_global, to_iparticle_kind, to_iparticle_local, &
         to_nparticle_kind, to_nparticle_local, tot_from_nparticle_local, tot_to_nparticle_local
      LOGICAL                                            :: found_it
      
      !MN Is this necessary?
      
      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(from_local_particles))
      CPASSERT(ASSOCIATED(to_local_particles))

      IF (.NOT. ASSOCIATED(from_local_particles%local_particle_set)) RETURN
      CPASSERT(.NOT. ASSOCIATED(to_local_particles%local_particle_set))

      from_nparticle_kind = from_local_particle_kinds%n_els
      to_nparticle_kind = to_local_particle_kinds%n_els

      ! make sure total number of particles hasn't changed, even if particle kinds have
      tot_from_nparticle_local = 0
      DO from_iparticle_kind = 1, from_nparticle_kind
         tot_from_nparticle_local = tot_from_nparticle_local + from_local_particles%n_el(from_iparticle_kind)
      END DO
      tot_to_nparticle_local = 0
      DO to_iparticle_kind = 1, to_nparticle_kind
         tot_to_nparticle_local = tot_to_nparticle_local + to_local_particles%n_el(to_iparticle_kind)
      END DO
      CPASSERT(tot_from_nparticle_local == tot_to_nparticle_local)

      ALLOCATE (to_local_particles%local_particle_set(to_nparticle_kind))
      DO to_iparticle_kind = 1, to_nparticle_kind

         to_nparticle_local = to_local_particles%n_el(to_iparticle_kind)
         ALLOCATE (to_local_particles%local_particle_set(to_iparticle_kind)%rng(to_nparticle_local))

         DO to_iparticle_local = 1, to_nparticle_local
            to_iparticle_global = to_local_particles%list(to_iparticle_kind)%array(to_iparticle_local)
            ALLOCATE (to_local_particles%local_particle_set(to_iparticle_kind)%rng(to_iparticle_local)%stream)

            found_it = .FALSE.
            ! find the matching kind/index where this particle was before
            DO from_iparticle_kind = 1, from_nparticle_kind
               from_nparticle_local = from_local_particles%n_el(from_iparticle_kind)
               IF (MINVAL(ABS(from_local_particles%list(from_iparticle_kind)%array(1:from_nparticle_local) - &
                              to_iparticle_global)) == 0) THEN
                  from_iparticle_local = &
                     MINLOC(ABS(from_local_particles%list(from_iparticle_kind)%array(1:from_nparticle_local) - &
                                to_iparticle_global))
                  to_local_particles%local_particle_set(to_iparticle_kind)%rng(to_iparticle_local)%stream = &
                     from_local_particles%local_particle_set(from_iparticle_kind)%rng(from_iparticle_local(1))%stream
                  found_it = .TRUE.
                  EXIT
               ENDIF
            END DO
            CPASSERT(found_it)

         END DO ! to_iparticle_local

      END DO ! to_iparticle_kind
      CALL timestop(handle)

   END SUBROUTINE copy_wiener_process

END MODULE adqmmmx_update
