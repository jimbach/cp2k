!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines used for force-mixing QM/MM calculations
!> \par History
!>      2.2012 created [noam]
!> \author Noam Bernstein
! **************************************************************************************************
MODULE adqmmmx_util
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              set_atomic_kind
   USE cell_types,                      ONLY: cell_copy,&
                                              cell_type,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr
!   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_deallocate,&
                                              fist_neighbor_type,&
                                              neighbor_kind_pairs_type
   USE fist_neighbor_lists,             ONLY: build_fist_neighbor_lists
   USE input_section_types,             ONLY: &
        section_vals_add_values, section_vals_duplicate, section_vals_get, &
        section_vals_get_subs_vals, section_vals_get_subs_vals3, section_vals_release, &
        section_vals_remove_values, section_vals_set_subs_vals, section_vals_type, &
        section_vals_val_get, section_vals_val_set, section_vals_write
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE physcon,                         ONLY: femtoseconds, kcalmol
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_list_types,             ONLY: molecule_list_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
   USE qmmm_types,                      ONLY: qmmm_env_get
   USE qmmm_types_low,                  ONLY: adqmmmx_label_QM_core,&
                                              adqmmmx_label_QM_core_list,&
                                              adqmmmx_label_QM_dynamics,&
                                              adqmmmx_label_QM_dynamics_list,&
                                              adqmmmx_label_QM_starting,&
                                              adqmmmx_label_QM_starting_list,&
                                              adqmmmx_label_none,&
                                              adqmmmx_label_termination
   USE qmmm_util,                       ONLY: apply_qmmm_translate
   USE adqmmmx_types,                   ONLY: adqmmmx_env_type
   USE force_env_types,                 ONLY: force_env_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE string_utilities,                ONLY: compress
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'adqmmmx_util'

   PUBLIC :: setup_adqmmmx_sections, update_adqmmmx_labels, &
             apply_adqmmmx_translate, dg_of_lambda, g_of_lambda, &
             get_adqmmmx_values, update_adqmmmx_values, &
             print_adqmmmx_output, copy_adqmmmx_env
             !&,old_update_adqmmmx_labels

CONTAINS

! **************************************************************************************************
!> \brief Apply translation to the full system in order to center the QM
!>      system into the QM box
!> \param adqmmmx_env ...
!> \par History
!>      08.2007 created [tlaino] - Zurich University
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE apply_adqmmmx_translate(adqmmmx_env, output_unit)
      TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_adqmmmx_translate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ip, output_unit
      TYPE(cell_type), POINTER                           :: cell_oldsys, cell_newsys
      TYPE(cp_subsys_type), POINTER                      :: subsys_oldsys, subsys_newsys
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles_oldsys, particles_newsys

      NULLIFY (cell_oldsys, cell_newsys)
      NULLIFY (subsys_oldsys, subsys_newsys)
      NULLIFY (particles_oldsys, particles_newsys)

      ! want to center newsys and make oldsys consistent with that
      CALL apply_qmmm_translate(adqmmmx_env%newsys)


!     !MN QMMM Wrapping is problematic, lets try to align old and new coordinates and wrap both systems
!     !according to QM labeling

      WRITE (output_unit,*) "SUBSYS CHECK", ASSOCIATED(subsys_oldsys),ASSOCIATED(subsys_newsys)
      WRITE (output_unit,*) "CELL CHECK", ASSOCIATED(cell_oldsys),ASSOCIATED(cell_newsys)
      WRITE (output_unit,*) "ENV CHECK ", ASSOCIATED(adqmmmx_env)


      ! translate core fist particles
      CALL qmmm_env_get(adqmmmx_env%newsys, subsys=subsys_newsys)
      CALL cp_subsys_get(subsys_newsys, cell=cell_newsys)

!     WRITE (output_unit,*) "NEW SUBSYS CHECK", ASSOCIATED(subsys_newsys)
!     WRITE (output_unit,*) "NEW CELL CHECK", ASSOCIATED(cell_newsys)

      CALL qmmm_env_get(adqmmmx_env%oldsys, subsys=subsys_oldsys)
      CALL cp_subsys_get(subsys_oldsys, cell=cell_oldsys)
!     WRITE (output_unit,*) "OLD SUBSYS CHECK", ASSOCIATED(subsys_oldsys)
!     WRITE (output_unit,*) "OLD CELL CHECK", ASSOCIATED(cell_oldsys)

      particles_newsys => subsys_newsys%particles%els
      particles_oldsys => subsys_oldsys%particles%els

      DO ip = 1, SIZE(particles_newsys)
      !MN Overwrite oldsys coords with newsys coords
         particles_oldsys(ip)%r = particles_newsys(ip)%r
         particles_oldsys(ip)%v = particles_newsys(ip)%v
      END DO

      CALL cell_copy(cell_newsys, cell_oldsys)

!     !Maybe one wrap of newsys is sufficient
!     CALL apply_qmmm_translate(adqmmmx_env%newsys)
!     CALL apply_qmmm_translate(adqmmmx_env%oldsys)

      ! The core QM particles will be updated the regular call
      ! to apply_qmmm_translate() from within qmmm_calc_energy_force()

   END SUBROUTINE apply_adqmmmx_translate


! **************************************************************************************************
!> \brief ...
!> \param subsys ...
!> \par History
!>      12.2021 created
!> \author Marvin Nyenhuis
! **************************************************************************************************
   SUBROUTINE update_adqmmmx_labels(subsys,qmmm_section,labels_changed,include_nc, c_n_sw)
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: qmmm_section, adqmmmx_section,&
                                                            restart_section,qm_kind_section
      REAL(dp), POINTER                                  :: r_qm(:)
      INTEGER, POINTER                                   :: cur_indices(:), cur_labels(:),&
                                                            mm_index_entry(:),&
                                                            new_indices(:),&
                                                            new_labels(:), qm_core(:),&
                                                            starting_qm_indices(:),&
                                                            starting_qm_labels(:)
      LOGICAL                                            :: file_exists, file_exists2, mol_in!,&
                                                            ! explicit
      LOGICAL, OPTIONAL                                  :: labels_changed
      LOGICAL                                            :: include_nc
      INTEGER                                            :: max_n_qm,&
                                                            natoms,&
                                                            ip,i_rep_section,i_rep_val,&
                                                            n_rep_section,n_rep_val,&
                                                            i,  i_start,i_end,&
                                                            ind_counter, &
                                                            starting_qm_counter,&
                                                            c_n_sw, max_n_sw
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(fist_neighbor_type), POINTER                  :: nlist

      !MN Create temporary output files
      INQUIRE(FILE = "util_log.out", EXIST=file_exists)
      IF (file_exists) THEN
         OPEN (UNIT=1338, FILE = "util_log.out", STATUS = "old", POSITION = "append", ACTION = "write")
      ELSE
         OPEN (UNIT=1338, FILE = "util_log.out", STATUS = "new", ACTION = "write")
      END IF

      INQUIRE(FILE = "natrium_distance.out", EXIST=file_exists2)
      IF (file_exists2) THEN
         OPEN (UNIT=1339, FILE = "natrium_distance.out", STATUS = "old", POSITION = "append", ACTION = "write")
      ELSE
         OPEN (UNIT=1339, FILE = "natrium_distance.out", STATUS = "new", ACTION = "write")
      END IF

      !MN Get current indices and labels of QM region
      CALL get_adqmmmx_indices(qmmm_section, cur_indices,cur_labels)

      !MN Get atom information
      NULLIFY (particles, molecules, cell)
      CALL cp_subsys_get(subsys=subsys, particles=particles, molecules=molecules, cell=cell)
      particle_set => particles%els
      molecule_set => molecules%els
      natoms = SIZE(particle_set)

      !MN Get input parameter r_qm for adaptive QMMMX
      adqmmmx_section => section_vals_get_subs_vals3(qmmm_section,"ADAPTIVE_QMMM")
      NULLIFY(r_qm)
      CALL section_vals_val_get(adqmmmx_section, "R_QM", r_vals=r_qm)
      CALL section_vals_val_get(adqmmmx_section, "MAX_N_QM", i_val=max_n_qm)
      WRITE(1338,*) "R_QM is ", r_qm, "MAX_N_QM ", max_n_qm

      !MN Create tmp_indices with length of natoms

      NULLIFY (new_indices, new_labels, starting_qm_indices, starting_qm_labels)
      CALL reallocate (new_indices, 1, natoms)
      CALL reallocate (new_labels, 1, natoms)
      CALL reallocate (starting_qm_indices, 1, natoms)
      CALL reallocate (starting_qm_labels, 1, natoms)
      new_indices = 0
      new_labels  = adqmmmx_label_none
      starting_qm_indices = 0
      starting_qm_labels  = adqmmmx_label_none

      ind_counter = 1
      starting_qm_counter = 1

      WRITE(*,*) "I AM HERE IREP"

      !MN Try to get Core section directly from input
      !MN Get New indices, sort it
      NULLIFY (mm_index_entry, qm_core)
      CALL section_vals_val_get(adqmmmx_section, "QM_CORE", i_vals=qm_core)
      qm_kind_section => section_vals_get_subs_vals3(qmmm_section, "QM_KIND")
      CALL section_vals_get(qm_kind_section, n_repetition=n_rep_section)
      DO i_rep_section = 1, n_rep_section
         CALL section_vals_val_get(qm_kind_section, "MM_INDEX", i_rep_section=i_rep_section, n_rep_val=n_rep_val)
         DO i_rep_val = 1, n_rep_val
            CALL section_vals_val_get(qm_kind_section, "MM_INDEX", i_rep_section=i_rep_section, i_rep_val=i_rep_val, &
                                      i_vals=mm_index_entry)
            DO ip = 1, SIZE(mm_index_entry)
               IF (ANY(qm_core==mm_index_entry(ip))) THEN
                  WRITE (1338,*) "Core atom found ",mm_index_entry(ip) 
                  new_indices(ind_counter) = mm_index_entry(ip)
                  new_labels(ind_counter)  = adqmmmx_label_QM_core
                  ind_counter = ind_counter + 1
               ELSE IF ( include_nc ) THEN
                  WRITE (1338,*) "NC atom found ",mm_index_entry(ip) 
                  starting_qm_indices(starting_qm_counter) = mm_index_entry(ip)
                  starting_qm_labels(starting_qm_counter) = adqmmmx_label_QM_starting
                  starting_qm_counter = starting_qm_counter + 1
               END IF
            END DO ! ip
         END DO ! i_rep_val
      END DO ! i_rep_section

      starting_qm_counter = starting_qm_counter - 1
      CALL reallocate (starting_qm_indices, 1, starting_qm_counter)
      CALL reallocate (starting_qm_labels, 1, starting_qm_counter)
 
      WRITE(1338,*) "Starting QM Indices"
      WRITE(1338,*) starting_qm_indices
      !MN Get neighbor list
      NULLIFY (nlist)
      CALL make_neighbor_list(adqmmmx_section, subsys, cell, r_qm(2), nlist)

      !MN Create new index list
      CALL add_new_index(nlist, qm_core, particle_set, cell, 1339,&
                                "ADQMMMX NEIGHBOR LIST", "angstrom", new_indices,&
                                new_labels, ind_counter)

      !MN Add all atoms of a molecule containing QM center to new_indices
      DO i=1,SIZE(molecule_set) 
         mol_in = .FALSE.
         i_start = molecule_set(i)%first_atom
         i_end   = molecule_set(i)%last_atom

         IF ( (max_n_qm/=0) .AND. (ind_counter + i_end - i_start) > max_n_qm) THEN
            CPWARN("TOO MANY QM PARTICLES IN BOX. STOPPED ADDING NEW MOLECULES")
            EXIT
         END IF

         DO ip=i_start,i_end
            IF (ANY(new_indices==ip)) mol_in = .TRUE. 
         END DO    
         IF (mol_in) THEN
         WRITE(1338,*) "Molecule in QM region: ", i
         DO ip=i_start,i_end
            IF (.NOT. ANY(new_indices==ip)) THEN
               IF (ANY(starting_qm_indices==ip)) THEN
               new_indices(ind_counter) = ip
               new_labels(ind_counter) = adqmmmx_label_QM_starting
               WRITE (1338,*) "Atom of QM molecule added to new_indices ",& 
                              ip," With label ", adqmmmx_label_QM_starting
               ELSE
                  new_indices(ind_counter) = ip
                  new_labels(ind_counter) = adqmmmx_label_QM_dynamics
               WRITE (1338,*) "Atom of QM molecule added to new_indices ",& 
                              ip," With label ", adqmmmx_label_QM_dynamics
               END IF
               ind_counter = ind_counter + 1
            END IF 
         END DO    
         END IF
      END DO

      !MN I think i have to retract the first +1 of new_counter
      ind_counter = ind_counter - 1

      CALL reallocate (new_indices, 1, ind_counter)
      CALL reallocate (new_labels, 1, ind_counter)
      !MN Sorting algorithm 
      CALL sort_indlab(cur_indices, cur_labels)
      CALL sort_indlab(new_indices, new_labels)

      !MN Catches the situations of true updates
      IF ( PRESENT(labels_changed) ) THEN
         CALL index_vec_compare(cur_indices,new_indices,labels_changed)
         CALL section_vals_val_get(adqmmmx_section, "MAX_N_SWITCH", i_val=max_n_sw)
         IF ((c_n_sw>=max_n_sw) .AND. (max_n_sw/=-1)) THEN
            WRITE (1338,*) "MAX_N_SWITCH REACHED", c_n_sw, max_n_sw
            labels_changed = .FALSE.
         END IF
         IF ( labels_changed ) THEN
            WRITE (1338,*) "Labels changed at time "
            WRITE (1338,*) "Cur: ", cur_indices
            WRITE (1338,*) "New: ", new_indices
            WRITE (1338,*) "SIZES Indic:", SIZE(cur_indices), SIZE(new_indices)
            WRITE (1338,*) "SIZES Labels:", SIZE(cur_labels), SIZE(new_labels)
            !MN Overwrite oldsys labels and indices with current indices/labels
            restart_section => section_vals_get_subs_vals(qmmm_section, "ADAPTIVE_QMMM%RESTART_INFO")
            CALL section_vals_val_set(restart_section, "OLD_INDICES", i_vals_ptr=cur_indices)
            CALL section_vals_val_set(restart_section, "OLD_LABELS", i_vals_ptr=cur_labels)
            CALL section_vals_val_set(restart_section, "INDICES", i_vals_ptr=new_indices)
            CALL section_vals_val_set(restart_section, "LABELS", i_vals_ptr=new_labels)
            CALL section_vals_write(qmmm_section,unit_nr=1338)

         ELSE
            WRITE (1338,*) "Labels didnt change"
            WRITE (1338,*) "Cur: ", cur_indices
            WRITE (1338,*) "New: ", new_indices
            WRITE (1338,*) "SIZES:", SIZE(cur_indices), SIZE(new_indices)
         END IF
      !MN Catches the case of adqmmmx_create
      ELSE IF ( include_nc .AND. (.NOT. PRESENT(labels_changed)) ) THEN
            WRITE (1338,*) "First creation of old and new lists"
            CALL get_adqmmmx_indices(qmmm_section, cur_indices, cur_labels)
            CALL sort_indlab(cur_indices, cur_labels)
            CALL sort_indlab(new_indices, new_labels)
            WRITE (1338,*) "Cur: ", cur_indices
            WRITE (1338,*) "New: ", new_indices
            WRITE (1338,*) "SIZES Indic:", SIZE(cur_indices), SIZE(new_indices)
            WRITE (1338,*) "SIZES Labels:", SIZE(cur_labels), SIZE(new_labels)
            !MN Overwrite oldsys labels and indices with current indices/labels
            restart_section => section_vals_get_subs_vals(qmmm_section, "ADAPTIVE_QMMM%RESTART_INFO")
            CALL section_vals_val_set(restart_section, "OLD_INDICES", i_vals_ptr=cur_indices)
            CALL section_vals_val_set(restart_section, "OLD_LABELS", i_vals_ptr=cur_labels)
            CALL section_vals_val_set(restart_section, "INDICES", i_vals_ptr=new_indices)
            CALL section_vals_val_set(restart_section, "LABELS", i_vals_ptr=new_labels)
            CALL section_vals_write(qmmm_section,unit_nr=1338)

      END IF

      WRITE(1338,*) "SUMMARY OF NEW QM ATOMS"

      
      ! save info in input structure
!     restart_section => section_vals_get_subs_vals(qmmm_section, "ADAPTIVE_QMMM%RESTART_INFO")
!     CALL section_vals_get(restart_section, explicit=explicit) !MN 02.08. changed
!     IF (explicit) CALL section_vals_remove_values(restart_section) 
      !MN Save all QM centers in restart information

      IF (1338 > 0) THEN
         WRITE (unit=1338, fmt='(A,A,I6,A,I6)') &
            "ADQMMMX final count (not including links): ", &
            " N_QM core     ", COUNT(new_labels == adqmmmx_label_QM_core), &
            " N_QM starting ", COUNT(new_labels == adqmmmx_label_QM_starting), &
            " N_QM dynamic  ", COUNT(new_labels == adqmmmx_label_QM_dynamics)
      ENDIF

!     CLOSE(1339)
!     CLOSE(1338)

!     DEALLOCATE (cur_indices, cur_labels)
      CALL fist_neighbor_deallocate(nlist)


   END SUBROUTINE update_adqmmmx_labels

! **************************************************************************************************
!> \brief ...
!> \param subsys ...
!> \par History
!>      12.2021 created
!> \author Marvin Nyenhuis
! **************************************************************************************************
  SUBROUTINE sort_indlab(ind, lab)
      INTEGER, POINTER                                   :: ind(:), lab(:)
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: temp_ind, temp_lab
      INTEGER                                            :: i, loc

      IF (SIZE(ind)/=SIZE(lab)) CPABORT("sort_indlab SIZES DONT MATCH")
      ALLOCATE(temp_ind(SIZE(ind)), temp_lab(SIZE(lab)))

      DO i=1,SIZE(ind)
         temp_ind(i) = ind(i)
         temp_lab(i) = lab(i)
      END DO

      !MN Sorting starts here
      DO i=1, SIZE(ind)
         loc = MINLOC(temp_ind, 1)
         ind(i) = temp_ind(loc)
         lab(i) = temp_lab(loc)
         temp_ind(loc) = HUGE(0)
      END DO

  END SUBROUTINE sort_indlab

! **************************************************************************************************
!> \brief ...
!> \param subsys ...
!> \par History
!>      12.2021 created
!> \author Marvin Nyenhuis
! **************************************************************************************************
  SUBROUTINE get_adqmmmx_values(qmmm_section, adqmmmx_env, output_unit)
     !TYPE(adqmmmx_env_type), POINTER, OPTIONAL          :: adqmmmx_env
     ! JB Mod 7.02
     !TYPE(adqmmmx_env_type), INTENT(INOUT), OPTIONAL    :: adqmmmx_env
     ! JB Mod 10.02
     TYPE(adqmmmx_env_type), INTENT(INOUT), OPTIONAL      :: adqmmmx_env
     TYPE(section_vals_type), POINTER                   :: qmmm_section,&
                                                           adqmmmx_section,&
                                                           restart_section
     LOGICAL                                            :: explicit
     INTEGER                                            :: output_unit, max_n_sw
     REAL (dp)                                          :: dummy_lambda


     WRITE(output_unit,*) "GET_ADQMMMX_VALUES"
     restart_section => section_vals_get_subs_vals(qmmm_section, "ADAPTIVE_QMMM%RESTART_INFO")
     adqmmmx_section => section_vals_get_subs_vals(qmmm_section, "ADAPTIVE_QMMM", can_return_null=.TRUE.)
     CALL section_vals_get(restart_section, explicit=explicit)
!    CALL section_vals_write(restart_section, unit_nr=output_unit)
     CALL section_vals_val_get(adqmmmx_section, "MAX_N_SWITCH", i_val=max_n_sw)
     IF (max_n_sw/=-1) THEN
       WRITE(output_unit,*) "Setting C_N_SW to 0"
       adqmmmx_env%c_n_sw=1
     END IF

     IF (explicit) THEN
        WRITE(output_unit,*) "RESTART SECTION IS PRESENT", SIZE(restart_section%subs_vals)
        CALL section_vals_val_get(restart_section, "LAMBDA", r_val=dummy_lambda)  
        IF (dummy_lambda >= 1) THEN 
           !MN Push Elambda to TOTElambda and obain propagation time           
           RETURN
        ELSE
           adqmmmx_env%lambda_adqmmmx = dummy_lambda
        END IF
        CALL section_vals_val_get(restart_section, "VLAMBDA", r_val=adqmmmx_env%lambda_velocity)  
        CALL section_vals_val_get(restart_section, "PROP_TIME", r_val=adqmmmx_env%prop_time)  
        CALL section_vals_val_get(restart_section, "LMASS", r_val=adqmmmx_env%lambda_mass)  
        CALL section_vals_val_get(restart_section, "OMEGA", r_val=adqmmmx_env%omega)  
        CALL section_vals_val_get(restart_section, "DPOT", r_val=adqmmmx_env%dpot)  
        CALL section_vals_val_get(restart_section, "LFORCE", r_val=adqmmmx_env%lambda_force)  
        CALL section_vals_val_get(restart_section, "ELAMBDA", r_val=adqmmmx_env%Elambda)  
        CALL section_vals_val_get(restart_section, "TOTELAMBDA", r_val=adqmmmx_env%Tot_Elambda)  
        !MN Set booleans
        WRITE (output_unit,*) "UPDATED VALUES", adqmmmx_env%lambda_adqmmmx, &
                              adqmmmx_env%lambda_velocity, adqmmmx_env%lambda_mass
        IF (0.0 < adqmmmx_env%lambda_adqmmmx) THEN
           adqmmmx_env%switch_active=.TRUE.
           adqmmmx_env%labels_changed=.FALSE.
           adqmmmx_env%energy_add=.TRUE.
           adqmmmx_env%output_energies=.FALSE.
        END IF 
     ELSE
         WRITE(output_unit,*) "RESTART SECTION IS NOT PRESENT"
         RETURN
     END IF

  END SUBROUTINE get_adqmmmx_values
! **************************************************************************************************
!> \brief ...
!> \param subsys ...
!> \par History
!>      12.2021 created
!> \author Marvin Nyenhuis
! **************************************************************************************************
  SUBROUTINE update_adqmmmx_values(adqmmmx_section, adqmmmx_env, output_unit)
     !TYPE(adqmmmx_env_type), POINTER, OPTIONAL          :: adqmmmx_env
     ! JB Mod 7.02
     TYPE(adqmmmx_env_type), INTENT(INOUT), OPTIONAL          :: adqmmmx_env
     TYPE(section_vals_type), POINTER                   :: adqmmmx_section,&
                                                           restart_section
     LOGICAL                                            :: explicit
     INTEGER                                            :: output_unit

     restart_section => section_vals_get_subs_vals(adqmmmx_section, "RESTART_INFO")
     CALL section_vals_get(restart_section, explicit=explicit)

     IF (explicit) THEN
        !CALL section_vals_remove_values(restart_section)
        WRITE (output_unit,*) "RESTART SECTION IS PRESENT, REPLACE VALUES" 

        CALL section_vals_val_set(restart_section, "LAMBDA", r_val=adqmmmx_env%lambda_adqmmmx)  
        CALL section_vals_val_set(restart_section, "VLAMBDA", r_val=adqmmmx_env%lambda_velocity)  
        CALL section_vals_val_set(restart_section, "PROP_TIME", r_val=adqmmmx_env%prop_time)  
        CALL section_vals_val_set(restart_section, "LMASS", r_val=adqmmmx_env%lambda_mass)  
        CALL section_vals_val_set(restart_section, "OMEGA", r_val=adqmmmx_env%omega)  
        CALL section_vals_val_set(restart_section, "DPOT", r_val=adqmmmx_env%dpot)  
        CALL section_vals_val_set(restart_section, "LFORCE", r_val=adqmmmx_env%lambda_force)  
        CALL section_vals_val_set(restart_section, "ELAMBDA", r_val=adqmmmx_env%Elambda)  
        CALL section_vals_val_set(restart_section, "ECORR", r_val=adqmmmx_env%Elamcorr)  
        CALL section_vals_val_set(restart_section, "DPOT_SGN", r_val=adqmmmx_env%dpot_sgn)  
        CALL section_vals_val_set(restart_section, "TOTELAMBDA", r_val=adqmmmx_env%Tot_Elambda) 

        WRITE (output_unit,*) "RESTART SECTION WAS UPDATED", adqmmmx_env%lambda_adqmmmx

     END IF
          
  END SUBROUTINE update_adqmmmx_values
! **************************************************************************************************
!> \brief ...
!> \param subsys ...
!> \par History
!>      12.2021 created
!> \author Marvin Nyenhuis
! **************************************************************************************************
  SUBROUTINE copy_adqmmmx_env(new_adqmmmx_env, adqmmmx_env)
     ! JB Mod 7.02
     TYPE(adqmmmx_env_type), INTENT(IN), OPTIONAL          :: adqmmmx_env
     TYPE(adqmmmx_env_type), INTENT(OUT), OPTIONAL          :: new_adqmmmx_env
      !   TYPE(adqmmmx_env_type), POINTER                    :: new_adqmmmx_env,&
      !                                                         adqmmmx_env


         new_adqmmmx_env%lambda_adqmmmx  = 0.0 
         new_adqmmmx_env%prop_time       = adqmmmx_env%prop_time
         new_adqmmmx_env%prop_time       = adqmmmx_env%prop_time
         new_adqmmmx_env%lambda_velocity = 0.0 
         new_adqmmmx_env%lambda_mass     = 0.0 
         new_adqmmmx_env%omega           = 0.0
         new_adqmmmx_env%dpot            = 0.0
         new_adqmmmx_env%totpot          = 0.0
         new_adqmmmx_env%kin_e           = adqmmmx_env%md_kin_e
         new_adqmmmx_env%md_kin_e        = 0.0
         new_adqmmmx_env%md_pot_e        = 0.0
         new_adqmmmx_env%md_corr         = 0.0
         new_adqmmmx_env%lambda_force    = 0.0
         new_adqmmmx_env%Elambda         = 0.0
         new_adqmmmx_env%Tot_Elambda     = adqmmmx_env%Tot_Elambda
         new_adqmmmx_env%Elamcorr        = adqmmmx_env%Elamcorr
         new_adqmmmx_env%dpot_sgn        = adqmmmx_env%dpot_sgn
         new_adqmmmx_env%switch_active   = .FALSE.
         new_adqmmmx_env%do_discrete    = adqmmmx_env%do_discrete
         new_adqmmmx_env%include_nc      = .FALSE.

         !MN This is only called when switching procedure started
         IF (.NOT. adqmmmx_env%c_n_sw==-1) THEN
           new_adqmmmx_env%c_n_sw = adqmmmx_env%c_n_sw + 1
           print *, "COPY_ADQMMMX| Add 1 to c_n_sw"
         END IF

  END SUBROUTINE copy_adqmmmx_env
! **************************************************************************************************
!> \brief Copied from neighbor_lists.F
!> \param nonbonded ...
!> \param particle_set ...
!> \param cell ...
!> \param para_env ...
!> \param output_unit ...
!> \param name ...
!> \param unit_str ...
!> \par History
!>      03.22 created 
!> \author MN
! **************************************************************************************************
   SUBROUTINE index_vec_compare(cur_indices,new_indices,labels_changed)
      INTEGER                                            :: i
      INTEGER, POINTER                                   :: cur_indices(:), new_indices(:)
      LOGICAL                                            :: labels_changed, equal
     
      equal = SIZE(cur_indices)==SIZE(new_indices)
      WRITE (1337,*) "IVC Step I", equal
      IF (equal) THEN
         labels_changed = .FALSE. 
         WRITE (1337,*) "IVC Step II", equal, labels_changed
         DO i=1, SIZE(cur_indices)
            equal = cur_indices(i)==new_indices(I)
            IF (.NOT. equal) THEN
               labels_changed = .TRUE.
               WRITE (1337,*) "IVC Step III", equal, labels_changed
               RETURN
            END IF
         END DO
      ELSE
         labels_changed = .TRUE.
         RETURN
      END IF
 
      WRITE (1337,*) "IVC Step IV", equal, labels_changed
      
!     IF (cur_indices==new_indices) THEN
!        labels_changed = .FALSE.
!        WRITE (1337,*) "Index_vec_compare case I"
!        RETURN
!     END IF      
!     IF (SIZE(cur_indices)/=SIZE(new_indices)) THEN
!        labels_changed = .TRUE.
!        WRITE (1337,*) "Index_vec_compare case II"
!        RETURN
!     END IF
!     
!     DO i=1, SIZE(cur_indices)
!        IF (new_indices(i)/=cur_indices(i)) THEN
!           labels_changed=.TRUE.
!           WRITE (1337,*) "Index_vec_compare case III"
!        END IF
!        WRITE (1337,*) "Index_vec_compare case IV"
!        RETURN
!     END DO                                               

!     DO i=1, SIZE(cur_indices)
!        IF (.NOT. ANY(new_indices==cur_indices(i))) THEN
!           labels_changed = .TRUE.
!           RETURN
!        ELSE
!           labels_changed = .FALSE.
!        END IF
!     END DO

!     DO i=1, SIZE(new_indices)
!        IF (.NOT. ANY(cur_indices==new_indices(i))) THEN
!           labels_changed = .TRUE.
!           RETURN
!        ELSE
!           labels_changed = .FALSE.
!        END IF
!     END DO

   END SUBROUTINE index_vec_compare

! **************************************************************************************************
!> \brief Copied from neighbor_lists.F
!> \param nonbonded ...
!> \param particle_set ...
!> \param cell ...
!> \param para_env ...
!> \param output_unit ...
!> \param name ...
!> \param unit_str ...
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE add_new_index(nonbonded, qm_core, particle_set, cell,  output_unit, &
                                   name, unit_str, new_indices, new_labels, ind_counter)

      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: output_unit
      CHARACTER(LEN=*), INTENT(IN)                       :: name, unit_str

      CHARACTER(LEN=default_string_length)               :: string
      INTEGER                                            :: atom_a, atom_b, iab, ilist, mype, &
                                                            nneighbor, ind_counter, new_ind
      INTEGER, POINTER                                   :: new_indices(:), new_labels(:),&
                                                            qm_core(:)
      LOGICAL                                            :: print_headline, file_exists
      REAL(dp)                                           :: conv, dab
      REAL(dp), DIMENSION(3)                             :: cell_v, ra, rab, rb
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair


      INQUIRE(FILE = "natrium_distance.out", EXIST=file_exists)
      IF (file_exists) THEN
         OPEN (UNIT=output_unit, FILE = "natrium_distance.out", STATUS = "old", POSITION = "append", ACTION = "write")
      ELSE
         OPEN (UNIT=output_unit, FILE = "natrium_distance.out", STATUS = "new", ACTION = "write")
      END IF

      mype = 0
      ! Print headline
      string = ""
      WRITE (UNIT=string, FMT="(A,I5,A)") &
         TRIM(name)//" IN "//TRIM(unit_str)//" (PROCESS", mype, ")"
      CALL compress(string)
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") TRIM(string)
      WRITE (output_unit,*) "QM CORE: ", qm_core
      print_headline = .TRUE.
      nneighbor = 0
      conv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_str))
      DO iab = 1, SIZE(nonbonded%neighbor_kind_pairs)
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(iab)
         cell_v = MATMUL(cell%hmat, REAL(neighbor_kind_pair%cell_vector, KIND=dp))
         DO ilist = 1, neighbor_kind_pair%npairs
            nneighbor = nneighbor + 1
            IF (output_unit > 0) THEN
               ! Print second part of the headline
               atom_a = neighbor_kind_pair%list(1, ilist)
               atom_b = neighbor_kind_pair%list(2, ilist)
!              WRITE (output_unit,*) "PROBABLE PAIR : ", atom_a, atom_b
               IF ((ANY(qm_core==atom_a)) .OR. (ANY(qm_core==atom_b))) THEN
!                 WRITE (output_unit,*) "VALID PAIR"
               ELSE
                  CYCLE
               END IF
               !MN Get only relevant interactions and add them to tmp_indices
               IF (ANY( new_indices==atom_a)) THEN
                  new_ind = atom_b
!                 WRITE (output_unit,*) "VALID PAIR: ", atom_b, atom_a
               ELSE IF (ANY( new_indices==atom_b)) THEN
                  new_ind = atom_a
!                 WRITE (output_unit,*) "VALID PAIR: ", atom_a, atom_b
               ELSE
                  CYCLE
               END IF
               !MN Remove self-interaction of core particles
               IF (ANY(new_indices==new_ind)) CYCLE


               !MN Handle tmp_indices vector
               new_indices(ind_counter) = new_ind
               new_labels(ind_counter)  = adqmmmx_label_QM_dynamics
               ind_counter = ind_counter + 1

               IF (print_headline) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,2(A6,3(5X,A,5X)),1X,A11,10X,A8)") &
                     "Atom-A", "X", "Y", "Z", "Atom-B", "X", "Y", "Z", "Cell(i,j,k)", &
                     "Distance"
                  print_headline = .FALSE.
               END IF

               ra(:) = pbc(particle_set(atom_a)%r, cell)
               rb(:) = pbc(particle_set(atom_b)%r, cell)
               rab = rb(:) - ra(:) + cell_v
               dab = SQRT(DOT_PRODUCT(rab, rab))
               IF (ilist <= neighbor_kind_pair%nscale) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,2(I6,3(1X,F10.6)),3(1X,I3),10X,F8.4,L4,F11.5,F9.5)") &
                     atom_a, ra(1:3)*conv, &
                     atom_b, rb(1:3)*conv, &
                     neighbor_kind_pair%cell_vector, &
                     dab*conv, &
                     neighbor_kind_pair%is_onfo(ilist), &
                     neighbor_kind_pair%vdw_scale(ilist), &
                     neighbor_kind_pair%ei_scale(ilist)
               ELSE
                  WRITE (UNIT=output_unit, FMT="(T3,2(I6,3(1X,F10.6)),3(1X,I3),10X,F8.4)") &
                     atom_a, ra(1:3)*conv, &
                     atom_b, rb(1:3)*conv, &
                     neighbor_kind_pair%cell_vector, &
                     dab*conv
               END IF
            END IF
         END DO ! ilist
      END DO ! iab

      string = ""
      WRITE (UNIT=string, FMT="(A,I12,A,I12)") &
         "Total number of neighbor interactions for process", mype, ":", &
         nneighbor

      CALL compress(string)
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,T2,A)") TRIM(string)

   END SUBROUTINE add_new_index

! **************************************************************************************************
!> \brief Copied from neighbor_lists.F
!> \param nonbonded ...
!> \param particle_set ...
!> \param cell ...
!> \param para_env ...
!> \param output_unit ...
!> \param name ...
!> \param unit_str ...
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE write_neighbor_lists(nonbonded, particle_set, cell, para_env, output_unit, &
                                   name, unit_str, cur_indices, new_indices)

      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
   !   TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: output_unit
      CHARACTER(LEN=*), INTENT(IN)                       :: name, unit_str

      CHARACTER(LEN=default_string_length)               :: string
      INTEGER                                            :: atom_a, atom_b, iab, ilist, mype, &
                                                            nneighbor, new_ind
      INTEGER, ALLOCATABLE                               :: new_indices(:)
      INTEGER, POINTER                                   :: cur_indices(:)
      LOGICAL                                            :: print_headline, file_exists
      REAL(dp)                                           :: conv, dab
      REAL(dp), DIMENSION(3)                             :: cell_v, ra, rab, rb
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(mp_para_env_type), POINTER                    :: para_env


      INQUIRE(FILE = "natrium_distance.out", EXIST=file_exists)
      IF (file_exists) THEN
         OPEN (UNIT=output_unit, FILE = "natrium_distance.out", STATUS = "old", POSITION = "append", ACTION = "write")
      ELSE
         OPEN (UNIT=output_unit, FILE = "natrium_distance.out", STATUS = "new", ACTION = "write")
      END IF

      mype = para_env%mepos
      ! Print headline
      string = ""
      WRITE (UNIT=string, FMT="(A,I5,A)") &
         TRIM(name)//" IN "//TRIM(unit_str)//" (PROCESS", mype, ")"
      CALL compress(string)
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") TRIM(string)

      WRITE(output_unit,*) "cur_indices associated", ASSOCIATED(cur_indices)      

      print_headline = .TRUE.
      nneighbor = 0
      conv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_str))
      DO iab = 1, SIZE(nonbonded%neighbor_kind_pairs)
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(iab)
         cell_v = MATMUL(cell%hmat, REAL(neighbor_kind_pair%cell_vector, KIND=dp))
         DO ilist = 1, neighbor_kind_pair%npairs
            nneighbor = nneighbor + 1
            IF (output_unit > 0) THEN
               ! Print second part of the headline
               atom_a = neighbor_kind_pair%list(1, ilist)
               atom_b = neighbor_kind_pair%list(2, ilist)

               !MN Get only relevant interactions and add them to new_indices
               IF (ANY( cur_indices==atom_a )) THEN
                  new_ind = atom_b
               ELSE IF (ANY( cur_indices==atom_b )) THEN
                  new_ind = atom_a
               ELSE
                  CYCLE
               END IF
               !MN Remove self-interaction of core particles
               IF (ANY(cur_indices==new_ind)) CYCLE


               !MN Handle new_indices vector
               CALL addtolist(new_indices, new_ind)

               IF (print_headline) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,2(A6,3(5X,A,5X)),1X,A11,10X,A8)") &
                     "Atom-A", "X", "Y", "Z", "Atom-B", "X", "Y", "Z", "Cell(i,j,k)", &
                     "Distance"
                  print_headline = .FALSE.
               END IF

               ra(:) = pbc(particle_set(atom_a)%r, cell)
               rb(:) = pbc(particle_set(atom_b)%r, cell)
               rab = rb(:) - ra(:) + cell_v
               dab = SQRT(DOT_PRODUCT(rab, rab))
               IF (ilist <= neighbor_kind_pair%nscale) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,2(I6,3(1X,F10.6)),3(1X,I3),10X,F8.4,L4,F11.5,F9.5)") &
                     atom_a, ra(1:3)*conv, &
                     atom_b, rb(1:3)*conv, &
                     neighbor_kind_pair%cell_vector, &
                     dab*conv, &
                     neighbor_kind_pair%is_onfo(ilist), &
                     neighbor_kind_pair%vdw_scale(ilist), &
                     neighbor_kind_pair%ei_scale(ilist)
               ELSE
                  WRITE (UNIT=output_unit, FMT="(T3,2(I6,3(1X,F10.6)),3(1X,I3),10X,F8.4)") &
                     atom_a, ra(1:3)*conv, &
                     atom_b, rb(1:3)*conv, &
                     neighbor_kind_pair%cell_vector, &
                     dab*conv
               END IF
            END IF
         END DO ! ilist
      END DO ! iab

      string = ""
      WRITE (UNIT=string, FMT="(A,I12,A,I12)") &
         "Total number of neighbor interactions for process", mype, ":", &
         nneighbor
      WRITE (UNIT=string, FMT="(A,I12)") &
         "Total amount of new QM centers:", SIZE(new_indices)

      DO iab=1,SIZE(new_indices)
         WRITE (output_unit,*) "QM center:", new_indices(iab)
      END DO
      CALL compress(string)
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,T2,A)") TRIM(string)

   END SUBROUTINE write_neighbor_lists

! **************************************************************************************************
!> \brief Small routine to add new element to list
! **************************************************************************************************
   SUBROUTINE addtolist(list, element)
      CHARACTER(len=*), PARAMETER :: routineN = 'addtolist', routineP = moduleN//':'//routineN
      INTEGER                                            :: i_size
      INTEGER, ALLOCATABLE                               :: list(:)
      INTEGER, INTENT(IN)                                :: element
      INTEGER, ALLOCATABLE                               :: t_list(:)
!     LOGICAL                                            :: file_exists    
 
!     INQUIRE(FILE = "addtolist.out", EXIST=file_exists)
!     IF (file_exists) THEN
!        OPEN (UNIT=1340, FILE = "addtolist.out", STATUS = "old", POSITION = "append", ACTION = "write")
!     ELSE
!        OPEN (UNIT=1340, FILE = "addtolist.out", STATUS = "new", ACTION = "write")
!     END IF
!     WRITE (1340,*) "Start of addtolist"
!     WRITE (1340,*) "List Allocated:", ALLOCATED(list) 

      !MN Added to prevent double counting
      IF (ANY(list==element)) RETURN

      IF (ALLOCATED(list)) THEN

         i_size = SIZE(list)
         ALLOCATE (t_list(i_size+1))
         t_list(i_size+1)   = element
         DEALLOCATE (list)
         WRITE (1340,*) "2nd check: List allocated", ALLOCATED(list), ALLOCATED(t_list)

         CALL move_alloc(t_list, list)
         WRITE (1340,*) "Added new element", element
         WRITE (1340,*) "3nd check: List allocated", ALLOCATED(list), ALLOCATED(t_list)
      ELSE
         ALLOCATE (list(1))
         list(1) = element
!        WRITE (1340,*) "List created"

      END IF

!     CLOSE(1340)

   END SUBROUTINE addtolist
! **************************************************************************************************
!> \brief ...
!> \param ip ...
!> \param label ...
!> \param n_new ...
!> \param new_indices ...
!> \param new_labels ...
!> \param new_full_labels ...
!> \param max_n_qm ...
! **************************************************************************************************
   SUBROUTINE add_new_label(ip, label, n_new, new_indices, new_labels, new_full_labels, max_n_qm)
      INTEGER                                            :: ip, label, n_new
      INTEGER, POINTER                                   :: new_indices(:), new_labels(:)
      INTEGER                                            :: new_full_labels(:), max_n_qm

      CHARACTER(len=*), PARAMETER :: routineN = 'add_new_label', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, old_index

      IF (new_full_labels(ip) > adqmmmx_label_none) THEN ! already marked, just change mark
         old_index = -1
         DO i = 1, n_new
            IF (new_indices(i) == ip) THEN
               old_index = i
               EXIT
            ENDIF
         END DO
         IF (old_index <= 0) &
            CALL cp_abort(__LOCATION__, &
                          "add_new_label found atom with a label "// &
                          "already set, but not in new_indices array")
         new_labels(old_index) = label
      ELSE
         n_new = n_new + 1
         IF (n_new > max_n_qm) &
            CALL cp_abort(__LOCATION__, &
                          "add_new_label tried to add more atoms "// &
                          "than allowed by &ADAPTIVE_QMMM&MAX_N_QM!")
         IF (n_new > SIZE(new_indices)) CALL reallocate(new_indices, 1, n_new + 9)
         IF (n_new > SIZE(new_labels)) CALL reallocate(new_labels, 1, n_new + 9)
         new_indices(n_new) = ip
         new_labels(n_new) = label
      ENDIF
      new_full_labels(ip) = label
   END SUBROUTINE add_new_label
! **************************************************************************************************
!> \brief ...
!> \param force_mixing_section ...
!> \param subsys ...
!> \param cell ...
!> \param r_max ...
!> \param nlist ...
! **************************************************************************************************
   SUBROUTINE make_neighbor_list(adqmmmx_section, subsys, cell, r_max, nlist)
      TYPE(section_vals_type), POINTER                   :: adqmmmx_section
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp)                                           :: r_max
      TYPE(fist_neighbor_type), POINTER                  :: nlist

      CHARACTER(len=*), PARAMETER :: routineN = 'make_neighbor_list', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: kind_name
      CHARACTER(LEN=default_string_length), POINTER      :: kind_name_a(:)
      INTEGER                                            :: ik
      LOGICAL                                            :: skip_kind
      REAL(dp), ALLOCATABLE                              :: r_max_a(:, :), r_minsq_a(:, :)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

      ALLOCATE (r_max_a(SIZE(subsys%atomic_kinds%els), SIZE(subsys%atomic_kinds%els)))
      ALLOCATE (r_minsq_a(SIZE(subsys%atomic_kinds%els), SIZE(subsys%atomic_kinds%els)))
      r_max_a = r_max
      r_minsq_a = EPSILON(1.0_dp)

      ! save kind names
      ALLOCATE (kind_name_a(SIZE(subsys%atomic_kinds%els)))
      DO ik = 1, SIZE(subsys%atomic_kinds%els)
         atomic_kind => subsys%atomic_kinds%els(ik)
         CALL get_atomic_kind(atomic_kind=atomic_kind, name=kind_name)
         kind_name_a(ik) = kind_name
      END DO

      ! overwrite kind names so that none are QM, and so excluding QM-QM interactions
      ! (which is not what we want) will not happen
      DO ik = 1, SIZE(subsys%atomic_kinds%els)
         atomic_kind => subsys%atomic_kinds%els(ik)
         CALL get_atomic_kind(atomic_kind=atomic_kind, name=kind_name)
         ! when atom is QM atom, kind_name is replaced with original
         ! mm kind name, and return status is logical .TRUE.
         skip_kind = qmmm_ff_precond_only_qm(kind_name)
         CALL set_atomic_kind(atomic_kind=atomic_kind, name=kind_name)
      END DO

      NULLIFY (nlist)
      CALL build_fist_neighbor_lists(subsys%atomic_kinds%els, subsys%particles%els, &
                                     cell=cell, r_max=r_max_a, r_minsq=r_minsq_a, &
                                     ei_scale14=1.0_dp, vdw_scale14=1.0_dp, nonbonded=nlist, &
                                     para_env=subsys%para_env, build_from_scratch=.TRUE., geo_check=.FALSE., &
                                     mm_section=adqmmmx_section)

      DEALLOCATE (r_max_a, r_minsq_a)

      ! restore kind names
      DO ik = 1, SIZE(subsys%atomic_kinds%els)
         CALL set_atomic_kind(atomic_kind=atomic_kind, name=kind_name_a(ik))
      END DO
      DEALLOCATE (kind_name_a)

   END SUBROUTINE make_neighbor_list

! **************************************************************************************************
!> \brief ...
!> \param subsys ...
!> \param qmmm_section ...
!> \param qmmm_core_section ...
!> \param qmmm_extended_section ...
!> \par History
!>      02.2012 created [noam]
!>      01.2022 adjusted MN
!> \author Noam Bernstein
! **************************************************************************************************
   SUBROUTINE setup_adqmmmx_sections(subsys, qmmm_section, new_qmmm_section, use_old_indlab)
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: qmmm_section, &
                                                            new_qmmm_section

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_adqmmmx_sections', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length), POINTER      :: elem_mapping(:, :), elem_mapping_entry(:)
      INTEGER                                            :: i_rep_section, i_rep_val,&
                                                            qm_kind_counter, &
                                                            ielem, ip, n_elements, output_unit
      INTEGER, POINTER                                   :: cur_indices(:), cur_labels(:)
      LOGICAL                                            :: new_element, mapped
      LOGICAL, OPTIONAL                                  :: use_old_indlab
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(section_vals_type), POINTER                   :: adqmmmx_section, &
                                                            qm_kind_section
      
      qm_kind_counter = 0 
      output_unit = cp_logger_get_default_unit_nr()
      !MN Duplicate qmmm section and remove starting QM_KINDS
      NULLIFY(new_qmmm_section)
      CALL section_vals_duplicate(qmmm_section, new_qmmm_section)
      qm_kind_section => section_vals_get_subs_vals3(new_qmmm_section, "QM_KIND")
      CALL section_vals_remove_values(qm_kind_section)

      !MN get QM_KIND_ELEMENT_MAPPING and current indices/labels
      adqmmmx_section => section_vals_get_subs_vals(qmmm_section, "ADAPTIVE_QMMM")
      CALL section_vals_val_get(adqmmmx_section, "QM_KIND_ELEMENT_MAPPING", n_rep_val=n_elements)
      ALLOCATE (elem_mapping(2, n_elements))
      DO ielem = 1, n_elements
         CALL section_vals_val_get(adqmmmx_section, "QM_KIND_ELEMENT_MAPPING", i_rep_val=ielem, c_vals=elem_mapping_entry)
         elem_mapping(1:2, ielem) = elem_mapping_entry(1:2)
      END DO

      IF (PRESENT(use_old_indlab) .AND. use_old_indlab) THEN !MN Only used for adqmmmx_create.F
         CALL get_adqmmmx_indices(qmmm_section, cur_indices, cur_labels, use_old_indlab)
      ELSE
         CALL get_adqmmmx_indices(qmmm_section, cur_indices, cur_labels)
      END IF
      IF (SIZE(cur_indices) <= 0) &
         CPABORT("cur_indices is empty, found no QM atoms")

      !MN Output current indices and labels
      WRITE (output_unit, *) "USE_OLD_INDLAB?", PRESENT(use_old_indlab)
      WRITE (output_unit, *) "cur_indices ", cur_indices
      WRITE (output_unit, *) "cur_labels ", cur_labels

      ! loop through elements and atoms, and check Mapping
      particles => subsys%particles%els

      DO ip = 1, SIZE(cur_indices)
!       IF (cur_labels(ip) > adqmmmx_label_none .AND. cur_labels(ip) /= adqmmmx_label_QM_core .AND. &
!           cur_labels(ip) /= adqmmmx_label_QM_starting .AND. cur_labels(ip) /= adqmmmx_label_termination) THEN
        IF (cur_labels(ip) > adqmmmx_label_termination) THEN
           mapped = .FALSE.
            DO ielem = 1, n_elements
              IF (TRIM(particles(cur_indices(ip))%atomic_kind%element_symbol) == TRIM(elem_mapping(1, ielem))) THEN
                 mapped = .TRUE.
                 !MN Print out mapping
                  WRITE (output_unit, *) "cur ind. ", cur_indices(ip), "Mapping ", elem_mapping(1, ielem),&
                                         elem_mapping(2,ielem)
                 EXIT
              ENDIF
            END DO
           IF (.NOT. mapped) &
              CALL cp_abort(__LOCATION__, &
                            "ADQMMMX failed to find QM_KIND mapping for atom of type "// &
                            TRIM(particles(cur_indices(ip))%atomic_kind%element_symbol)// &
                            "! ")
        ENDIF
      END DO

      !MN Create counter for sections
      i_rep_section = 0

      DO ielem = 1, n_elements !MN Cycle over all elements
         new_element = .TRUE. !MN Checker for new elements in QM_KINDs
         DO ip = 1, SIZE(cur_indices) !MN Cycle over all indices
            IF (TRIM(particles(cur_indices(ip))%atomic_kind%element_symbol) /= TRIM(elem_mapping(1, ielem))) CYCLE
            !MN Check if correct labeling has happend
            IF ((cur_labels(ip) == adqmmmx_label_QM_core) .OR. &
                (cur_labels(ip) == adqmmmx_label_QM_dynamics) .OR. &
                (cur_labels(ip) == adqmmmx_label_QM_starting)) THEN
               !MN ADD QM_KIND_SECTION for center with element ielem and index ip
                  qm_kind_section => section_vals_get_subs_vals3(new_qmmm_section, "QM_KIND")
                  qm_kind_counter = qm_kind_counter + 1
                  IF (new_element) THEN ! add new QM_KIND section for this element
                     i_rep_section = i_rep_section + 1
                     CALL section_vals_add_values(qm_kind_section)
                     CALL section_vals_val_set(qm_kind_section, "_SECTION_PARAMETERS_", i_rep_section=i_rep_section, &
                                                              c_val=elem_mapping(2, ielem))
                     
                     i_rep_val = 0
                     new_element = .FALSE.
                  ENDIF
                  i_rep_val = i_rep_val + 1
                  CALL section_vals_val_set(qm_kind_section, "MM_INDEX", i_rep_section=i_rep_section, &
                                               i_rep_val=i_rep_val, i_val=cur_indices(ip))
               
 
            END IF

         END DO ! atom index ip
      END DO ! element index ielem
      WRITE(output_unit, *) "Len cur_ind", size(cur_indices), "Len QM_KIND", qm_kind_counter
!     !MN Write old qmmm section
!     WRITE (output_unit, *) "old_qmmm_section"
!     CALL section_vals_write(qmmm_section, unit_nr=output_unit)
!     !MN Write new qmmm section
!     WRITE (output_unit, *) "new_qmmm_section"
!     CALL section_vals_write(new_qmmm_section, unit_nr=output_unit)

   END SUBROUTINE setup_adqmmmx_sections

!
! **************************************************************************************************
!> \brief ...
!> \param force_mixing_section ...
!> \param indices ...
!> \param labels ...
! **************************************************************************************************
   SUBROUTINE get_adqmmmx_indices(qmmm_section, indices, labels, use_old_indlab)
      TYPE(section_vals_type), POINTER                   :: qmmm_section,adqmmmx_section,&
                                                            qm_kind_section
      INTEGER, POINTER                                   :: indices(:), labels(:), test_ind(:),&
                                                            mm_index_entry(:)

      CHARACTER(len=*), PARAMETER :: routineN = 'get_adqmmmx_indices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i_rep_val, ip,&
                                                            n_indices, n_labels, n_reps,&
                                                            i_rep_section,n_rep_section,&
                                                            n_rep_val,ind_counter
      INTEGER, POINTER                                   :: indices_entry(:), labels_entry(:),&
                                                            qm_core(:)
      LOGICAL                                            :: explicit
      LOGICAL, OPTIONAL                                  :: use_old_indlab
      TYPE(section_vals_type), POINTER                   :: restart_section
                         

      adqmmmx_section => section_vals_get_subs_vals3(qmmm_section, "ADAPTIVE_QMMM")

      restart_section => section_vals_get_subs_vals(adqmmmx_section, "RESTART_INFO")
      CALL section_vals_get(restart_section, explicit=explicit)
      !MN first call
      NULLIFY (indices, labels)
      IF (.NOT. explicit) THEN ! no old indices, labels, hardcode lists of SIZE100
         n_indices = 0
         n_labels = 0
         GO TO 100
      END IF

      CALL section_vals_val_get(restart_section, "OLD_INDICES",explicit=explicit, i_vals=test_ind)
      WRITE(1338,*) "Get ADQMMMX INDICES EXPLICIT? ", explicit, test_ind
      if ( .NOT. explicit) use_old_indlab = .FALSE.

     ![NB] maybe switch to reallocatable array
     IF ((PRESENT(use_old_indlab) .AND. (use_old_indlab))) THEN !MN Only used in adqmmmx_create.F
        CALL section_vals_val_get(restart_section, "OLD_INDICES", n_rep_val=n_reps)
        n_indices = 0
        n_labels = 0
        DO i_rep_val = 1, n_reps
           CALL section_vals_val_get(restart_section, "OLD_INDICES", &
                                     i_rep_val=i_rep_val, i_vals=indices_entry)
           n_indices = n_indices + SIZE(indices_entry)
        END DO
        ALLOCATE (indices(n_indices))
        n_indices = 0
        n_labels = 0
        DO i_rep_val = 1, n_reps
           CALL section_vals_val_get(restart_section, "OLD_INDICES", &
                                     i_rep_val=i_rep_val, i_vals=indices_entry)
           indices(n_indices + 1:n_indices + SIZE(indices_entry)) = indices_entry
           n_indices = n_indices + SIZE(indices_entry)
        END DO
 
        CALL section_vals_val_get(restart_section, "OLD_LABELS", n_rep_val=n_reps)
        n_labels = 0
        DO i_rep_val = 1, n_reps
           CALL section_vals_val_get(restart_section, "OLD_LABELS", &
                                     i_rep_val=i_rep_val, i_vals=labels_entry)
           n_labels = n_labels + SIZE(labels_entry)
        END DO
        ALLOCATE (labels(n_labels))
        n_labels = 0
        DO i_rep_val = 1, n_reps
           CALL section_vals_val_get(restart_section, "OLD_LABELS", &
                                     i_rep_val=i_rep_val, i_vals=labels_entry)
           labels(n_labels + 1:n_labels + SIZE(labels_entry)) = labels_entry
           n_labels = n_labels + SIZE(labels_entry)
        END DO
 
     ELSE
        CALL section_vals_val_get(restart_section, "INDICES", n_rep_val=n_reps)
        n_indices = 0
        DO i_rep_val = 1, n_reps
           CALL section_vals_val_get(restart_section, "INDICES", &
                                     i_rep_val=i_rep_val, i_vals=indices_entry)
           n_indices = n_indices + SIZE(indices_entry)
        END DO
        ALLOCATE (indices(n_indices))
        n_indices = 0
        DO i_rep_val = 1, n_reps
           CALL section_vals_val_get(restart_section, "INDICES", &
                                     i_rep_val=i_rep_val, i_vals=indices_entry)
           indices(n_indices + 1:n_indices + SIZE(indices_entry)) = indices_entry
           n_indices = n_indices + SIZE(indices_entry)
        END DO
 
        CALL section_vals_val_get(restart_section, "LABELS", n_rep_val=n_reps)
        n_labels = 0
        DO i_rep_val = 1, n_reps
           CALL section_vals_val_get(restart_section, "LABELS", &
                                     i_rep_val=i_rep_val, i_vals=labels_entry)
           n_labels = n_labels + SIZE(labels_entry)
        END DO
        ALLOCATE (labels(n_labels))
        n_labels = 0
        DO i_rep_val = 1, n_reps
           CALL section_vals_val_get(restart_section, "LABELS", &
                                     i_rep_val=i_rep_val, i_vals=labels_entry)
           labels(n_labels + 1:n_labels + SIZE(labels_entry)) = labels_entry
           n_labels = n_labels + SIZE(labels_entry)
        END DO
        
     END IF

      !MN Case for empty restart section, probably not necessary!
      100 CONTINUE
      IF ((n_indices==0) .OR. (n_labels==0)) THEN
         WRITE(1338, *) "ADQMMMX | No Restart available, create cur_indices"
         n_indices = 0
         n_labels  = 0

         NULLIFY (qm_core)
         CALL section_vals_val_get(adqmmmx_section, "QM_CORE", i_vals=qm_core)
         WRITE(1338,*) "CORE", qm_core
         IF (SIZE(qm_core) > 0) THEN
            qm_kind_section => section_vals_get_subs_vals3(qmmm_section, "QM_KIND")

            !Get size of labels/indices
            CALL section_vals_get(qm_kind_section, n_repetition=n_rep_section)
            DO i_rep_section = 1, n_rep_section
               CALL section_vals_val_get(qm_kind_section, "MM_INDEX",&
                                         i_rep_section=i_rep_section, n_rep_val=n_rep_val)
               DO i_rep_val = 1, n_rep_val
                  CALL section_vals_val_get(qm_kind_section, "MM_INDEX",&
                                      i_rep_section=i_rep_section, i_rep_val=i_rep_val, &
                                      i_vals=mm_index_entry)
                  DO ip = 1, SIZE(mm_index_entry)
                     n_indices = n_indices + 1
                     n_labels = n_labels + 1
                  END DO ! ip
               END DO ! i_rep_val
            END DO ! i_rep_section
            !MN Create labels/indices for core and starting QM region
            ALLOCATE(indices(n_indices))
            ALLOCATE(labels(n_labels))
            ind_counter = 1
            DO i_rep_section = 1, n_rep_section
               CALL section_vals_val_get(qm_kind_section, "MM_INDEX",&
                                         i_rep_section=i_rep_section, n_rep_val=n_rep_val)
               DO i_rep_val = 1, n_rep_val
                  CALL section_vals_val_get(qm_kind_section, "MM_INDEX",&
                                      i_rep_section=i_rep_section, i_rep_val=i_rep_val, &
                                      i_vals=mm_index_entry)
                  WRITE(1338,*) "REP_VAL ", n_rep_val, "MMIE ",SIZE(mm_index_entry)
                  DO ip = 1, SIZE(mm_index_entry)
                     IF (ANY(qm_core==mm_index_entry(ip))) THEN
                        indices(ind_counter) = mm_index_entry(ip)
                        labels(ind_counter)  = adqmmmx_label_QM_core
                        ind_counter = ind_counter + 1
                     ELSE
                        indices(ind_counter) = mm_index_entry(ip)
                        labels(ind_counter)  = adqmmmx_label_QM_starting
                        ind_counter = ind_counter + 1
                     END IF
                  END DO ! ip
               END DO ! i_rep_val
            END DO ! i_rep_section
            WRITE(1338,*) "NIND: ", n_indices, "NLAB: ", n_labels
            WRITE(1338,*) "Indices ", indices
            WRITE(1338,*) "Labels ", labels
         ELSE
            CPABORT("NO CORE REGION DEFINED FOR ADAPTIVE QMMMX!")
         ENDIF
      END IF

      IF (.NOT. ((n_labels - n_indices)==0)) THEN
         WRITE(1338,*) "NIND (LIF): ", n_indices, "NLAB (LIF): ", n_labels
         CPWARN("got unequal numbers of adqmmmx indices and labels!")
      END IF
   END SUBROUTINE get_adqmmmx_indices

! **************************************************************************************************
!> \brief Calculation of lambda_mass
!> \param lambda ...
!> \par History
!>      03.2022
!> \author Marvin Nyenhuis
! **************************************************************************************************
   SUBROUTINE print_adqmmmx_output(adqmmmx_env)

      TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      REAL(dp), POINTER                                  :: lambda, time, vlambda, mlambda, flambda,&
                                                            dpot, totpot, md_corr,kin_e, &
                                                            elambda, elamcorr, tot_elambda, constant
      REAL(dp)                                           :: conserved
      LOGICAL                                            :: file_exists, adqmmmx_first_call=.FALSE.
      
      IF(.NOT. ASSOCIATED(adqmmmx_env)) RETURN
      !MN Output for switching energies
      INQUIRE(FILE = "energies.out", EXIST=file_exists)
      IF (file_exists) THEN
         OPEN (UNIT=1336, FILE = "energies.out", STATUS = "old", POSITION = "append", ACTION = "write")
      ELSE
         OPEN (UNIT=1336, FILE = "energies.out", STATUS = "new", ACTION = "write")
         WRITE (1336,201) "Time [fs]", "Kin.",  "Mx_Pot",& 
                          "Lam_E", "ELC",&
                          "TLE","Cons."
         adqmmmx_first_call = .TRUE.
      END IF

      IF ((adqmmmx_env%output_energies) .OR. (adqmmmx_first_call)) THEN
         time        => adqmmmx_env%prop_time
         dpot        => adqmmmx_env%dpot
         totpot      => adqmmmx_env%totpot
         md_corr     => adqmmmx_env%md_corr
         elambda     => adqmmmx_env%Elambda
         elamcorr    => adqmmmx_env%Elamcorr
         tot_elambda => adqmmmx_env%Tot_Elambda
         constant    => adqmmmx_env%constant
         IF (adqmmmx_env%kin_e >= adqmmmx_env%md_kin_e) THEN
                 kin_e => adqmmmx_env%kin_e
         ELSE
                 kin_e => adqmmmx_env%md_kin_e
         END IF

!        conserved = kin_e + totpot + elambda + tot_elambda + elamcorr
         WRITE (1336,202) time*femtoseconds, kin_e, totpot,&
                          elambda, elamcorr,&
                          tot_elambda,constant

      END IF
      CLOSE (1336)

      !MN Output for MD energies
      INQUIRE(FILE = "MD_energies.out", EXIST=file_exists)
      IF (file_exists) THEN
         OPEN (UNIT=1336, FILE = "MD_energies.out", STATUS = "old", POSITION = "append", ACTION = "write")
      ELSE
         OPEN (UNIT=1336, FILE = "MD_energies.out", STATUS = "new", ACTION = "write")
         WRITE (1336,203) "Time", "Kin.",  "Pot", "DPOT", "Cons."
         adqmmmx_first_call = .TRUE.
      END IF

      IF ((adqmmmx_env%output_energies) .OR. (adqmmmx_first_call)) THEN
         NULLIFY(dpot,kin_e, totpot, elambda, tot_elambda)
         conserved = 0.0_dp
         kin_e   => adqmmmx_env%md_kin_e
         md_corr => adqmmmx_env%md_corr
         dpot    => adqmmmx_env%dpot
         totpot  => adqmmmx_env%md_pot_e

         conserved   =  kin_e + totpot + md_corr

         WRITE (1336,204) time*femtoseconds, kin_e, totpot, dpot, conserved
      END IF
      CLOSE (1336)

      !MN Output for MD energies
      INQUIRE(FILE = "lambda.out", EXIST=file_exists)
      IF (file_exists) THEN
         OPEN (UNIT=1336, FILE = "lambda.out", STATUS = "old", POSITION = "append", ACTION = "write")
      ELSE
         OPEN (UNIT=1336, FILE = "lambda.out", STATUS = "new", ACTION = "write")
         WRITE (1336,205) "Time", "Lambda", "Vlambda", "FLambda", "G(Lambda)", "Mass"
         adqmmmx_first_call = .TRUE.
      END IF

         vlambda     => adqmmmx_env%lambda_velocity
         mlambda     => adqmmmx_env%lambda_mass
         flambda     => adqmmmx_env%lambda_force

         WRITE (1336,206) time*femtoseconds, lambda, vlambda, flambda, g_of_lambda(lambda),&
                          mlambda

!     END IF
      CLOSE (1336)

   !Format for energies (only MD)
   201 FORMAT (1X, A19, 1X, A19, 1X, A19, 1X, A19, 1X, A19, 1X, A19, 1X, A19)
   202 FORMAT (1X, F19.9, 1X, F19.9, 1X, F19.9, 1X, F19.9, 1X, F19.9,1X, F19.9, 1X, F19.9)
   203 FORMAT (1X, A19, 1X, A19, 1X, A19, 1X, A19,1X,A19)
   204 FORMAT (1X, F19.9, 1X, F19.9, 1X, F19.9, 1X, F19.9, 1X, F19.9,1X,F19.9)
   205 FORMAT (1X, A19, 1X, A19, 1X, A19, 1X, A19, 1X, A19, 1X, A19, 1X, A19)
   206 FORMAT (1X, F19.9, 1X, F19.9, 1X, F19.9, 1X, F19.9, 1X, F19.9, 1X, F19.9, 1X, F19.9)
      
   END SUBROUTINE print_adqmmmx_output
! **************************************************************************************************
!> \brief routine for the switching function g
!> \param lambda
!> \par History
!>      05.2018 created [jim]
!> \author Jim Bachmann
! **************************************************************************************************
   FUNCTION g_of_lambda(lambda) result(g_lam)
      REAL(dp), intent(in)                               :: lambda
      REAL(dp)                                           :: g_lam
      !MN Maybe use the output unit here instead of the terminal (WRITE)
      !print *, "g(lambda) is calculated right now"
      g_lam = (lambda - 1.d0)**2
   END FUNCTION

! **************************************************************************************************
!> \brief routine for the derivative of the switching function
!> \param lambda ...
!> \par History
!>      05.2018 created [jim]
!> \author Jim Bachmann
! **************************************************************************************************
   FUNCTION dg_of_lambda(lambda) result(dg)
      REAL(dp), intent(in)                               :: lambda
      REAL(dp)                                           :: dg
!     INTEGER                                            :: output_unit

!     output_unit = cp_logger_get_default_io_unit()
!     WRITE (output_unit,*) "dg(lambda) is calculated right now"
      dg = 2*(lambda - 1.0)
   END FUNCTION


END MODULE adqmmmx_util
