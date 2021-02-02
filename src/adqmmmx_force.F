!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates QM/MM energy and forces with Force-Mixing
!> \par History
!>      2015 Factored out of force_env_methods.F
!> \author Ole Schuett, modified by Jim Bachmann and Marvin Nyenhuis
! **************************************************************************************************
MODULE adqmmmx_force
   USE cell_types,                      ONLY: cell_type
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE fist_environment_types,          ONLY: fist_env_get
   USE input_constants,                 ONLY: do_fm_mom_conserv_QM,&
                                              do_fm_mom_conserv_buffer,&
                                              do_fm_mom_conserv_core,&
                                              do_fm_mom_conserv_equal_a,&
                                              do_fm_mom_conserv_equal_f,&
                                              do_fm_mom_conserv_none
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_types,                  ONLY: particle_type
   USE qmmm_force,                      ONLY: qmmm_calc_energy_force
   USE qmmm_types,                      ONLY: qmmm_env_get,&
                                              qmmm_env_type
   USE qmmm_types_low,                  ONLY: adqmmmx_label_QM_core,&
                                              adqmmmx_label_QM_dynamics,&
                                              adqmmmx_label_buffer
   USE qmmm_util,                       ONLY: apply_qmmm_unwrap,&
                                              apply_qmmm_wrap
   USE adqmmmx_types,                   ONLY: adqmmmx_env_type
   USE adqmmmx_util,                    ONLY: apply_adqmmmx_translate
   USE qs_environment_types,            ONLY: get_qs_env
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'adqmmmx_force'

   PUBLIC :: adqmmmx_calc_energy_force

CONTAINS

! **************************************************************************************************
!> \brief calculates the qm/mm energy and forces
!> \param adqmmmx_env ...
!> \param calc_force if also the forces should be calculated
!> \par History
!>      05.2004 created [fawzi]
!>      08.2020 modified by JB and MN
!> \author Fawzi Mohamed, modified by Jim Bachmann and Marvin Nyenhuis
! **************************************************************************************************
   SUBROUTINE adqmmmx_calc_energy_force(adqmmmx_env, calc_force, consistent_energies, linres, &
                                      require_consistent_energy_force)
      TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      LOGICAL, INTENT(IN)                                :: calc_force, consistent_energies, linres
      LOGICAL, INTENT(IN), OPTIONAL :: require_consistent_energy_force

      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_calc_energy_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ip, mom_conserv_min_label, &
                                                            mom_conserv_n, mom_conserv_region, &
                                                            mom_conserv_type
      INTEGER, POINTER                                   :: cur_indices(:), cur_labels(:)
      ! JB added lambda_adqmmmx epot1 epot2
      REAL(dp)                                           :: delta_a(3), delta_f(3), &
                                                            mom_conserv_mass, total_f(3), lambda_adqmmmx &
                                                            epot1, epot2
      TYPE(cp_subsys_type), POINTER                      :: subsys_primary, subsys_qmmm_core, &
                                                            subsys_qmmm_extended
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_primary, particles_qmmm_core, &
                                                            particles_qmmm_extended
      TYPE(section_vals_type), POINTER                   :: force_env_section

      IF (PRESENT(require_consistent_energy_force)) THEN
         IF (require_consistent_energy_force) &
            CALL cp_abort(__LOCATION__, &
                          "qmmmx_energy_and_forces got require_consistent_energy_force but force mixing is active. ")
      ENDIF

      ! Possibly translate the system
      CALL apply_adqmmmx_translate(adqmmmx_env)

      ! actual energy force calculation
      CALL adqmmmx_calc_energy_force_low(adqmmmx_env%ext, calc_force, consistent_energies, linres, "ext")
      CALL adqmmmx_calc_energy_force_low(adqmmmx_env%core, calc_force, consistent_energies, linres, "core")

      ! get forces from subsys of each sub force env
      ! JB added epot1 epot2
      CALL qmmm_env_get(adqmmmx_env%core, subsys=subsys_qmmm_core, potential_energy=epot1)
      CALL qmmm_env_get(adqmmmx_env%ext, subsys=subsys_qmmm_extended, potential_energy=epot2)

      CALL get_qs_env(adqmmmx_env%ext%qs_env, input=force_env_section)
      CALL section_vals_val_get(force_env_section, "QMMM%ADAPTIVE_QMMM%RESTART_INFO%INDICES", i_vals=cur_indices)
      CALL section_vals_val_get(force_env_section, "QMMM%ADAPTIVE_QMMM%RESTART_INFO%LABELS", i_vals=cur_labels)

      particles_qmmm_extended => subsys_qmmm_extended%particles%els
      particles_qmmm_core => subsys_qmmm_core%particles%els
      DO ip = 1, SIZE(cur_indices)
         IF (cur_labels(ip) >= adqmmmx_label_QM_dynamics) THEN ! this is a QM atom
            ! copy (QM) force from extended calculation
            particles_qmmm_core(cur_indices(ip))%f = particles_qmmm_extended(cur_indices(ip))%f
         ENDIF
      END DO

      ! JB zero momentum not needed for adqmmx
      ! zero momentum
      !CALL section_vals_val_get(force_env_section, "QMMM%ADAPTIVE_QMMM%MOMENTUM_CONSERVATION_TYPE", &
      !                          i_val=mom_conserv_type)
      !IF (mom_conserv_type /= do_fm_mom_conserv_none) THEN
      !   CALL section_vals_val_get(force_env_section, "QMMM%ADAPTIVE_QMMM%MOMENTUM_CONSERVATION_REGION", &
      !                             i_val=mom_conserv_region)
      !
      !   IF (mom_conserv_region == do_fm_mom_conserv_core) THEN
      !      mom_conserv_min_label = adqmmmx_label_QM_core
      !   ELSEIF (mom_conserv_region == do_fm_mom_conserv_QM) THEN
      !      mom_conserv_min_label = adqmmmx_label_QM_dynamics
      !   ELSEIF (mom_conserv_region == do_fm_mom_conserv_buffer) THEN
      !      mom_conserv_min_label = adqmmmx_label_buffer
      !   ELSE
      !      CPABORT("Got unknown MOMENTUM_CONSERVATION_REGION (not CORE, QM, or BUFFER) !")
      !   ENDIF

      !   total_f = 0.0_dp
      !   DO ip = 1, SIZE(particles_qmmm_core)
      !      total_f(1:3) = total_f(1:3) + particles_qmmm_core(ip)%f(1:3)
      !   END DO
      !   IF (mom_conserv_type == do_fm_mom_conserv_equal_f) THEN
      !      mom_conserv_n = COUNT(cur_labels >= mom_conserv_min_label)
      !      delta_f = total_f/mom_conserv_n
      !      DO ip = 1, SIZE(cur_indices)
      !         IF (cur_labels(ip) >= mom_conserv_min_label) THEN
      !            particles_qmmm_core(cur_indices(ip))%f = particles_qmmm_core(cur_indices(ip))%f - delta_f
      !         ENDIF
      !      END DO
      !   ELSE IF (mom_conserv_type == do_fm_mom_conserv_equal_a) THEN
      !      mom_conserv_mass = 0.0_dp
      !      DO ip = 1, SIZE(cur_indices)
      !         IF (cur_labels(ip) >= mom_conserv_min_label) &
      !            mom_conserv_mass = mom_conserv_mass + particles_qmmm_core(cur_indices(ip))%atomic_kind%mass
      !      END DO
      !      delta_a = total_f/mom_conserv_mass
      !      DO ip = 1, SIZE(cur_indices)
      !         IF (cur_labels(ip) >= mom_conserv_min_label) THEN
      !            particles_qmmm_core(cur_indices(ip))%f = particles_qmmm_core(cur_indices(ip))%f - &
      !                                                     particles_qmmm_core(cur_indices(ip))%atomic_kind%mass*delta_a
      !         ENDIF
      !      END DO
      !   ENDIF
      !ENDIF

      !JB mod get the adaptive lambda
      lambda_adqmmmx = adqmmx_env%lambda_adqmmmx

      !JB here we need to calc the force on lambda
      adqmmx_env%lambda_force = -dg(lambda_adqmmmx) * DABS(epot1 - epot2)

      !JB mixing the forces
      CALL qmmm_env_get(adqmmmx_env%ext, subsys=subsys_primary)
      particles_primary => subsys_primary%particles%els
      DO ip = 1, SIZE(particles_qmmm_core)
         particles_primary(ip)%f = g(lambda_adqmmmx)*particles_qmmm_core(ip)%f &
                                   + (1.d0 - g(lambda_adqmmmx))*particles_qmmm_extended(ip)%f
      END DO

   END SUBROUTINE adqmmmx_calc_energy_force

! **************************************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \param calc_force ...
!> \param consistent_energies ...
!> \param linres ...
!> \param label ...
! **************************************************************************************************
   SUBROUTINE adqmmmx_calc_energy_force_low(qmmm_env, calc_force, consistent_energies, linres, label)
      TYPE(qmmm_env_type), POINTER                       :: qmmm_env
      LOGICAL, INTENT(IN)                                :: calc_force, consistent_energies, linres
      CHARACTER(*)                                       :: label

      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_calc_energy_force_low', &
         routineP = moduleN//':'//routineN

      CHARACTER(default_string_length)                   :: new_restart_fn, new_restart_hist_fn, &
                                                            old_restart_fn, old_restart_hist_fn
      INTEGER, DIMENSION(:), POINTER                     :: qm_atom_index
      LOGICAL                                            :: saved_do_translate
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: saved_pos
      TYPE(cell_type), POINTER                           :: mm_cell
      TYPE(cp_subsys_type), POINTER                      :: subsys_mm, subsys_qm
      TYPE(section_vals_type), POINTER                   :: force_env_section

      NULLIFY (mm_cell, subsys_qm, subsys_mm, qm_atom_index)

      CALL get_qs_env(qmmm_env%qs_env, input=force_env_section)

      ! rewrite RESTART%FILENAME
      CALL section_vals_val_get(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME", &
                                c_val=old_restart_fn)
      new_restart_fn = TRIM(old_restart_fn)//"-"//TRIM(label)
      CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME", &
                                c_val=new_restart_fn)

      ! rewrite RESTART_HISTORY%FILENAME
      CALL section_vals_val_get(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME", &
                                c_val=old_restart_hist_fn)
      new_restart_hist_fn = TRIM(old_restart_hist_fn)//"-"//TRIM(label)
      CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME", &
                                c_val=new_restart_hist_fn)

      ! wrap positions before QM/MM calculation.
      ! Required if diffusion causes atoms outside of periodic box get added to QM
      CALL fist_env_get(qmmm_env%fist_env, cell=mm_cell, subsys=subsys_mm)
      CALL get_qs_env(qmmm_env%qs_env, cp_subsys=subsys_qm)
      qm_atom_index => qmmm_env%qm%qm_atom_index
      CALL apply_qmmm_wrap(subsys_mm, mm_cell, subsys_qm, qm_atom_index, saved_pos)

      ! Turn off box translation, it was already performed by apply_adqmmmx_translate(),
      ! the particles coordinates will still be copied from MM to QM.
      saved_do_translate = qmmm_env%qm%do_translate
      qmmm_env%qm%do_translate = .FALSE.

      ! actual energy force calculation
      CALL qmmm_calc_energy_force(qmmm_env, calc_force, consistent_energies, linres)

      ! restore do_translate
      qmmm_env%qm%do_translate = saved_do_translate

      ! restore unwrapped positions
      CALL apply_qmmm_unwrap(subsys_mm, subsys_qm, qm_atom_index, saved_pos)

      ! restore RESTART filenames
      CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME", &
                                c_val=old_restart_fn)
      CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME", &
                                c_val=old_restart_hist_fn)

   END SUBROUTINE adqmmmx_calc_energy_force_low

END MODULE adqmmmx_force
