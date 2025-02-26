!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Initialize a QM/MM calculation with adaptive QM/MM
!> \note changed names (qmmmx->adqmmmx)
!> \author Ole Schuett
! **************************************************************************************************
MODULE adqmmmx_create
   !USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE global_types,                    ONLY: global_environment_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_type, &
                                              section_vals_write
   USE message_passing,                 ONLY: mp_para_env_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr
   USE qmmm_create,                     ONLY: qmmm_env_create
   USE qmmm_types,                      ONLY: qmmm_env_get,&
                                              qmmm_env_release,&
                                              qmmm_env_type
   USE adqmmmx_types,                   ONLY: adqmmmx_env_type
   USE adqmmmx_util,                    ONLY: setup_adqmmmx_sections,&
                                              update_adqmmmx_labels,&
                                              get_adqmmmx_values
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'adqmmmx_create'

   PUBLIC :: adqmmmx_env_create

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param adqmmmx_env ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \note adjusted names (qmmmx->adqmmmx)
!> \par History
!>      02.2012 created [noam]
!> \author Noam Bernstein
! **************************************************************************************************
   SUBROUTINE adqmmmx_env_create(adqmmmx_env, root_section, para_env, globenv, &
                               force_env_section, subsys_section, use_motion_section)
      TYPE(adqmmmx_env_type), INTENT(OUT)                :: adqmmmx_env
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: force_env_section, subsys_section
      LOGICAL, INTENT(IN)                                :: use_motion_section
 
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(qmmm_env_type), POINTER                       :: dummy_qmmm_env
      TYPE(section_vals_type), POINTER                   :: new_qmmm_section, &
                                                            qmmm_section
      INTEGER                                            :: output_unit
 
      NULLIFY (dummy_qmmm_env)
      
      !JB Mod 27.2
      CPASSERT(ASSOCIATED(adqmmmx_env%oldsys))
      CPASSERT(ASSOCIATED(adqmmmx_env%newsys))
      
      output_unit = cp_logger_get_default_io_unit()
      qmmm_section => section_vals_get_subs_vals(force_env_section, "QMMM")
 
      !MN Get restart information and if Lambda>0, set booleans accordingly
      !ALLOCATE (adqmmmx_env)
      CALL get_adqmmmx_values(qmmm_section, adqmmmx_env, output_unit)

       WRITE (output_unit,*) "ADQMMMX Lambda starting"

      ALLOCATE (dummy_qmmm_env)
      IF (adqmmmx_env%lambda_adqmmmx > 0.0) THEN
         WRITE (output_unit,*) "ADQMMMX Lambda", adqmmmx_env%lambda_adqmmmx
         WRITE (output_unit,*) "CREATE ADQMMMX_ENV FROM RESTART FILE"
         CALL qmmm_env_create(dummy_qmmm_env, root_section, para_env, globenv, &
                           force_env_section, qmmm_section, subsys_section, use_motion_section, &
                           ignore_outside_box=.true.)
         CALL qmmm_env_get(dummy_qmmm_env, subsys=subsys)
          !MN Get oldsys using old indices and labels
         CALL setup_adqmmmx_sections(subsys, qmmm_section, new_qmmm_section, use_old_indlab=.TRUE.)
         
         
         !JB Add 7.02.2025
         WRITE(*,*) "I am here 1"
         !JB ADD 4 02 2025
         ALLOCATE(adqmmmx_env%oldsys)
         CALL qmmm_env_create(adqmmmx_env%oldsys, root_section, para_env, globenv, &
                              force_env_section, new_qmmm_section, subsys_section, use_motion_section, &
                              ignore_outside_box=.TRUE.)
         !MN Get newsys using new indices and labels
         CALL setup_adqmmmx_sections(subsys, qmmm_section, new_qmmm_section)
         !JB ADD 4 02 2025
         ALLOCATE(adqmmmx_env%newsys)
         CALL qmmm_env_create(adqmmmx_env%newsys, root_section, para_env, globenv, &
                              force_env_section, new_qmmm_section, subsys_section, use_motion_section, &
                              ignore_outside_box=.TRUE.)
      ELSE
 
         WRITE (output_unit,*) "ADQMMMX Lambda", adqmmmx_env%lambda_adqmmmx
         WRITE (output_unit,*) "CREATING A NEW ADQMMMX_ENV"
         CALL qmmm_env_create(dummy_qmmm_env, root_section, para_env, globenv, &
                           force_env_section, qmmm_section, subsys_section, use_motion_section, &
                           ignore_outside_box=.true.)
         CALL qmmm_env_get(dummy_qmmm_env, subsys=subsys)
         
         !MN Here indices are updated, maybe grep old updates and store them in adqmmmx_env_type?
         adqmmmx_env%include_nc = .TRUE.
         CALL update_adqmmmx_labels(subsys, qmmm_section, include_nc=.TRUE., c_n_sw=adqmmmx_env%c_n_sw)
   
         !JB Add 7.02.2025
         WRITE(*,*) "I am here 2"
   
         ! using CUR_INDICES and CUR_LABELS, create appropriate QM_KIND sections for two QM/MM calculations
         !MN And here I can create a whole new environment for the old subsys or push it on core?
         CALL setup_adqmmmx_sections(subsys, qmmm_section, new_qmmm_section)      
         WRITE (output_unit,*) "ADQMMMX | old_sys"
         CALL section_vals_write(qmmm_section, unit_nr=output_unit)
         ALLOCATE (adqmmmx_env%oldsys)
         CALL qmmm_env_create(adqmmmx_env%oldsys, root_section, para_env, globenv, &
                              force_env_section, qmmm_section, subsys_section, use_motion_section, &
                              ignore_outside_box=.TRUE.)
   
         WRITE (output_unit,*) "ADQMMMX | new_sys"
         CALL section_vals_write(new_qmmm_section, unit_nr=output_unit)
         ALLOCATE (adqmmmx_env%newsys)
         CALL qmmm_env_create(adqmmmx_env%newsys, root_section, para_env, globenv, &
                              force_env_section, new_qmmm_section, subsys_section, use_motion_section, &
                              ignore_outside_box=.TRUE.)
         
         !JB Add 7.02.2025
         WRITE(*,*) "I am here 3"
 
         CALL section_vals_release(new_qmmm_section)
         WRITE (output_unit,*) "ADQMMMX | adqmmmx_env_create"
      END IF
      
      CALL qmmm_env_release(dummy_qmmm_env)
      DEALLOCATE (dummy_qmmm_env)
 
   END SUBROUTINE adqmmmx_env_create

END MODULE adqmmmx_create 
