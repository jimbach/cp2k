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
   USE physcon,                         ONLY: femtoseconds
   USE mathconstants,                   ONLY: pi
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_subsys_types,                 ONLY: cp_subsys_type,&
   !MN Remove after connectivity is checked   
                                              cp_subsys_get
   USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                              get_molecule_kind, &
                                              bond_type, bend_type
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_type
   USE molecule_types,                  ONLY: molecule_type, &
                                              get_molecule
   USE molecule_list_types,             ONLY: molecule_list_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   !------------------------------------------------------------
   USE fist_environment_types,          ONLY: fist_env_get,&
                                              fist_environment_type
   !MN Added for output
   USE fist_energy_types,               ONLY: fist_energy_type
   USE simpar_types,                    ONLY: simpar_type
   USE force_env_types,                 ONLY: force_env_type
   USE input_constants,                 ONLY: do_fm_mom_conserv_QM,&
                                              do_fm_mom_conserv_buffer,&
                                              do_fm_mom_conserv_core,&
                                              do_fm_mom_conserv_equal_a,&
                                              do_fm_mom_conserv_equal_f,&
                                              do_fm_mom_conserv_none
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_get_subs_vals,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_types,                  ONLY: particle_type
   USE qmmm_force,                      ONLY: qmmm_calc_energy_force
   USE qmmm_types,                      ONLY: qmmm_env_get,&
                                              qmmm_env_type
   USE qmmm_types_low,                  ONLY: adqmmmx_label_QM_core,&
                                              adqmmmx_label_QM_dynamics
   USE qmmm_util,                       ONLY: apply_qmmm_unwrap,&
                                              apply_qmmm_wrap
   USE adqmmmx_types,                   ONLY: adqmmmx_env_type
   USE adqmmmx_util,                    ONLY: apply_adqmmmx_translate,&
                                              g_of_lambda, dg_of_lambda
   USE qs_environment_types,            ONLY: get_qs_env
   !MN Temporarily added output_unit
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr
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
!>      08.2021 modified by MN
!> \author Fawzi Mohamed, modified by Jim Bachmann and Marvin Nyenhuis
! **************************************************************************************************
   SUBROUTINE adqmmmx_calc_energy_force(force_env, calc_force, consistent_energies, linres)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      TYPE(fist_energy_type), POINTER                    :: thermo
      LOGICAL, INTENT(IN)                                :: calc_force, consistent_energies, linres
      LOGICAL                                            :: file_exists
      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_calc_energy_force', &
         routineP = moduleN//':'//routineN

      !MN Added output_unit
      INTEGER                                            :: ip, output_unit, iatoms !,mom_conserv_min_label, &
                                                            !MN mom_conserv_n, mom_conserv_region, &
                                                            !MN mom_conserv_type !unused variables
      !INTEGER, POINTER                                   :: cur_indices(:), cur_labels(:)
      REAL(dp)                                           :: epot1, epot2, glam, totpot,&
                                                            kin_e, mass,&
                                                            md_kin_e!,&
                                                            !oldf_l, newf_l, alpha, dist
      REAL(dp), DIMENSION(3)                             :: old_part, new_part,vel!, distv
      TYPE(cp_subsys_type), POINTER                      :: subsys_primary, subsys_qmmm_oldsys, &
                                                            subsys_qmmm_newsys
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_primary, particles_qmmm_oldsys, &
                                                            particles_qmmm_newsys, particles_fist
      TYPE(section_vals_type), POINTER                   :: force_env_section

      output_unit = cp_logger_get_default_io_unit()
       
      adqmmmx_env => force_env%adqmmmx_env

      !MN Inserted output file entry to control force calculation
      IF (ASSOCIATED(force_env%adqmmmx_env)) THEN
      INQUIRE(FILE = "force_log.out", EXIST=file_exists)
      IF (file_exists) THEN
         OPEN (UNIT=1335, FILE = "force_log.out", STATUS = "old", POSITION = "append", ACTION = "write")
      ELSE
         OPEN (UNIT=1335, FILE = "force_log.out", STATUS = "new", ACTION = "write")
      END IF
      END IF

      WRITE(1335, *) "SWITCHING Active?", adqmmmx_env%switch_active
      ! Possibly translate the system
      CALL apply_adqmmmx_translate(adqmmmx_env, output_unit)
      !MN Calculate energies for old and new systems
      IF (adqmmmx_env%switch_active) THEN
          CALL adqmmmx_calc_energy_force_low(adqmmmx_env%oldsys, calc_force, consistent_energies, linres, "oldsys")
          CALL adqmmmx_calc_energy_force_low(adqmmmx_env%newsys, calc_force, consistent_energies, linres, "newsys")

          ! MN Obtain potential energies from old and new systems
          CALL qmmm_env_get(adqmmmx_env%oldsys, subsys=subsys_qmmm_oldsys, potential_energy=epot1)
          CALL qmmm_env_get(adqmmmx_env%newsys, subsys=subsys_qmmm_newsys, potential_energy=epot2)

          CALL get_qs_env(adqmmmx_env%newsys%qs_env, input=force_env_section)

          !MN Get forces of new and old systems
          particles_qmmm_newsys => subsys_qmmm_newsys%particles%els
          particles_qmmm_oldsys => subsys_qmmm_oldsys%particles%els
     
          IF (epot2 /= epot1) THEN
             adqmmmx_env%dpot         = epot1 - epot2
          ELSE
              adqmmmx_env%dpot = SIGN(0.d0, epot1 - epot2)
              WRITE (1335,*) "DPOT WAS SET TO 0"
          END IF

          !MN Get sign of first DeltaV
          IF (adqmmmx_env%dpot_sgn==0.0) THEN
              adqmmmx_env%dpot_sgn = SIGN(1.d0,epot1 - epot2)
              WRITE (1335,*) "Creating dpot_sgn: ", adqmmmx_env%dpot_sgn
          END IF

          !MN here we need to calc the force on lambda
          adqmmmx_env%lambda_force=-dg_of_lambda(adqmmmx_env%lambda_adqmmmx)*DABS(epot1 - epot2)
          WRITE(1335,*) "Lambda Force: ", adqmmmx_env%lambda_force

          !MN Get particle environment
          CALL qmmm_env_get(adqmmmx_env%newsys, subsys=subsys_primary)
          particles_primary => subsys_primary%particles%els


          glam = g_of_lambda(adqmmmx_env%lambda_adqmmmx)
          totpot = glam*epot1  + (1-glam)*epot2
          adqmmmx_env%totpot = totpot
     
          !MN Calculate kinetic energy
          NULLIFY (thermo)
    
          CALL fist_env_get(adqmmmx_env%newsys%fist_env, thermo=thermo,&
                            particle_set=particles_fist)
          kin_e = thermo%kin_shell
!         e_bonded = thermo%e_bonded
!         pot = thermo%pot
!         kin = thermo%kin

          !MN Get connectivity
!         CALL get_connectivity(adqmmmx_env%oldsys%fist_env, "oldsys", 1335)
!         CALL get_connectivity(adqmmmx_env%newsys%fist_env, "newsys", 1335)
          DO iatoms = 1,SIZE(particles_fist)
             mass  = particles_fist(iatoms)%atomic_kind%mass
             vel   = particles_fist(iatoms)%v
             kin_e = kin_e + 0.5*mass*DOT_PRODUCT(vel,vel)

          END DO
          !MN Ugly hack to have consistent kinetic energies
          IF (kin_e .EQ. 0.0) kin_e = md_kin_e
          adqmmmx_env%kin_e = kin_e
          adqmmmx_env%md_pot_e = epot2

          !MN Mix forces for particles
          DO ip = 1, SIZE(particles_qmmm_oldsys)
             old_part = glam*particles_qmmm_oldsys(ip)%f
             new_part = (1.0-glam)*particles_qmmm_newsys(ip)%f
             particles_primary(ip)%f = old_part + new_part
!            oldf_l = DSQRT(DOT_PRODUCT(particles_qmmm_oldsys(ip)%f,particles_qmmm_oldsys(ip)%f))
!            newf_l = DSQRT(DOT_PRODUCT(particles_qmmm_newsys(ip)%f,particles_qmmm_newsys(ip)%f))
!            alpha  = DOT_PRODUCT(particles_qmmm_oldsys(ip)%f,particles_qmmm_newsys(ip)%f)
!            alpha  = alpha/oldf_l
!            alpha  = alpha/newf_l
!            alpha  = DACOS(alpha)*180.0/pi
!            distv  = particles_qmmm_oldsys(ip)%r - particles_qmmm_newsys(ip)%r
!            dist   = DSQRT(DOT_PRODUCT(distv, distv))
!            WRITE(1335,*) "Distance: ", dist
!            IF (oldf_l/newf_l<0.9 .OR. oldf_l/newf_l>1.1) THEN
!                   WRITE(1335,*) "Force ratio of ", ip, ": ",&
!                          oldf_l/newf_l, alpha
!            END IF
         
!         WRITE(1335,*) "Thermo stats: ", e_bonded, pot, kin
          !MN Track epot1, epot2 and totpot
!         WRITE(1335,*) "T; ", adqmmmx_env%prop_time*femtoseconds, "L: ",&
!                       adqmmmx_env%lambda_adqmmmx,"Pot1: ", epot1,&
!                       "Pot2: ", epot2, "TPot: ", totpot
          END DO

      ELSE
          WRITE (1335,*) "ONLY CALCULATE NEWSYS FORCES"
          CALL adqmmmx_calc_energy_force_low(adqmmmx_env%newsys, calc_force, consistent_energies, linres, "newsys")
          CALL qmmm_env_get(adqmmmx_env%newsys, subsys=subsys_qmmm_newsys, potential_energy=epot2)

          CALL get_qs_env(adqmmmx_env%newsys%qs_env, input=force_env_section)
          !MN Get forces of new system
          particles_qmmm_newsys => subsys_qmmm_newsys%particles%els
          !MN Get particle environment
          CALL qmmm_env_get(adqmmmx_env%newsys, subsys=subsys_primary)
          particles_primary => subsys_primary%particles%els

          totpot = epot2
          adqmmmx_env%totpot = totpot
     
          !MN Calculate kinetic energy
          NULLIFY (thermo)
    
          CALL fist_env_get(adqmmmx_env%newsys%fist_env, thermo=thermo,&
                            particle_set=particles_fist)
          kin_e = thermo%kin_shell

          DO iatoms = 1,SIZE(particles_fist)
             mass  = particles_fist(iatoms)%atomic_kind%mass
             vel   = particles_fist(iatoms)%v
             kin_e = kin_e + 0.5*mass*DOT_PRODUCT(vel,vel)

          END DO
          !MN Ugly hack to have consistent kinetic energies
          IF (kin_e .EQ. 0.0) kin_e = md_kin_e
          adqmmmx_env%kin_e = kin_e
          adqmmmx_env%md_pot_e = epot2

          !MN Mix forces for particles
          DO ip = 1, SIZE(particles_qmmm_newsys)
             particles_primary(ip)%f = particles_qmmm_newsys(ip)%f
          END DO
         
      END IF

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
      INTEGER                                            :: output_unit
      LOGICAL                                            :: saved_do_translate
!     REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: saved_pos
      TYPE(cell_type), POINTER                           :: mm_cell
      TYPE(cp_subsys_type), POINTER                      :: subsys_mm, subsys_qm
      TYPE(section_vals_type), POINTER                   :: force_env_section

      NULLIFY (mm_cell, subsys_qm, subsys_mm, qm_atom_index)
      

      output_unit = cp_logger_get_default_io_unit()
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
      !MN Try without wrap
      !CALL apply_qmmm_wrap(subsys_mm, mm_cell, subsys_qm, qm_atom_index, saved_pos)

      ! Turn off box translation, it was already performed by apply_adqmmmx_translate(),
      ! the particles coordinates will still be copied from MM to QM.
      saved_do_translate = qmmm_env%qm%do_translate
      qmmm_env%qm%do_translate = .FALSE.

      ! actual energy force calculation
      CALL qmmm_calc_energy_force(qmmm_env, calc_force, consistent_energies, linres)

      ! restore do_translate
      qmmm_env%qm%do_translate = saved_do_translate

      ! restore unwrapped positions
      !CALL apply_qmmm_unwrap(subsys_mm, subsys_qm, qm_atom_index, saved_pos)
      WRITE(output_unit,*) label, " QM:", SIZE(qmmm_env%qm%qm_atom_index),&
                           "MM:", SIZE(qmmm_env%qm%mm_atom_index)
      ! restore RESTART filenames
      CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME", &
                                c_val=old_restart_fn)
      CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME", &
                                c_val=old_restart_hist_fn)

   END SUBROUTINE adqmmmx_calc_energy_force_low

! **************************************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \param calc_force ...
!> \param consistent_energies ...
!> \param linres ...
!> \param label ...
! **************************************************************************************************
   SUBROUTINE get_connectivity(fist_env, label, output_unit)
      TYPE(fist_environment_type), POINTER              :: fist_env
      TYPE(molecule_kind_type), POINTER                 :: molecule_kind_set(:)
      TYPE(molecule_kind_type), POINTER                 :: molecule_kind
      TYPE(molecule_type), POINTER                      :: molecule
      TYPE(distribution_1d_type), POINTER               :: local_molecules
      TYPE(molecule_type), POINTER                      :: molecule_set(:)
      TYPE(particle_type), POINTER                      :: particle_set(:)
      TYPE(cell_type), POINTER                          :: cell
      TYPE(bond_type), POINTER                          :: bond_list(:)
      TYPE(bend_type), POINTER                          :: bend_list(:)
      INTEGER                                           :: nkind, ikind, nmol, imol, output_unit, &
                                                           a, b, c, ibond, nbonds, ibend, nbends, &
                                                           first_atom, i
      REAL(dp), DIMENSION(3)                            :: rij
      CHARACTER(*)                                      :: label

      WRITE(1335,*) "Connectivity for", label
      CALL fist_env_get(fist_env, molecule_kind_set=molecule_kind_set,&
                        molecule_set=molecule_set,local_molecules=local_molecules,&
                        particle_set=particle_set, cell=cell)

      nkind = SIZE(molecule_kind_set)
      DO ikind = 1, nkind
         nmol = local_molecules%n_el(ikind)

         DO imol = 1, nmol
            i = local_molecules%list(ikind)%array(imol)
            molecule => molecule_set(i)
            molecule_kind => molecule%molecule_kind
            CALL get_molecule_kind(molecule_kind, nbend=nbends, nbond=nbonds,&
                                   bond_list=bond_list, bend_list=bend_list)
            CALL get_molecule(molecule, first_atom=first_atom)
            !BONDS
            DO ibond = 1, nbonds
               a = bond_list(ibond)%a + first_atom - 1
               b = bond_list(ibond)%b + first_atom - 1
               rij = particle_set(a)%r - particle_set(b)%r
               rij = pbc(rij, cell)
               WRITE(output_unit,*) "B: ", a, b, DSQRT(DOT_PRODUCT(rij,rij))
            END DO
           !BENDS
            DO ibend = 1, nbends
               a = bend_list(ibend)%a + first_atom - 1
               b = bend_list(ibend)%b + first_atom - 1
               c = bend_list(ibend)%c + first_atom - 1
               WRITE(output_unit,*) "A: ", a, b, c
            END DO
         END DO         

      END DO

   END SUBROUTINE get_connectivity

END MODULE adqmmmx_force
