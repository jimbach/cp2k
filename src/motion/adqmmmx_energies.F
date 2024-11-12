!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   !
!   Copyright (C) 2000 - 2019  CP2K developers group
!   !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines used for getting ADQMMMMX energies
!> \par History
!>      2.2012 created [noam]
!> \author Noam Bernstein
! **************************************************************************************************
MODULE adqmmmx_energies
   USE input_section_types,             ONLY: &
        section_vals_add_values, section_vals_duplicate, section_vals_get, &
        section_vals_get_subs_vals, section_vals_get_subs_vals3,section_vals_release, &
        section_vals_remove_values, section_vals_set_subs_vals,section_vals_type, &
        section_vals_val_get, section_vals_val_set, section_vals_write
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE adqmmmx_types,                   ONLY: adqmmmx_env_type
   USE force_env_types,                 ONLY: force_env_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE md_ener_types,                   ONLY: md_ener_type
   USE md_environment_types,            ONLY: get_md_env,&
                                              md_environment_type
   USE simpar_types,                    ONLY: simpar_type
   USE input_constants,                 ONLY: nve_ensemble, nvt_ensemble
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE thermostat_utils,                ONLY: get_thermostat_energies
   USE thermostat_types,                ONLY: thermostat_type
#include "../base/base_uses.f90"


   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'adqmmmx_util'

   PUBLIC :: adqmmmx_calc_constant
CONTAINS

! **************************************************************************************************
!> \brief Calculation of kinetic, potential and constant energy of NVE
!ensembles
!> \param lambda ...
!> \par History
!>      03.2022
!> \author Marvin Nyenhuis
! **************************************************************************************************
   SUBROUTINE adqmmmx_calc_constant(md_env, md_ener, simpar, para_env)
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(md_ener_type), POINTER                        :: md_ener
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(cp_para_env_type), POINTER                    :: para_env
!     TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      CHARACTER(LEN=*), PARAMETER :: routineN = 'adqmmmx_calc_constant', &
      routineP = moduleN//':'//routineN

      SELECT CASE (simpar%ensemble)    
      CASE DEFAULT
         CPABORT('Unknown ensemble')
      CASE (nve_ensemble)
         CALL adqmmmx_calc_nve_constant(md_env, md_ener, para_env) 
      CASE (nvt_ensemble)
         CALL adqmmmx_calc_nvt_constant(md_env, md_ener, para_env) 
      END SELECT

   END SUBROUTINE adqmmmx_calc_constant

! **************************************************************************************************
!> \brief Calculation of kinetic, potential and constant energy of NVE
!ensembles
!> \param lambda ...
!> \par History
!>      03.2022
!> \author Marvin Nyenhuis
! **************************************************************************************************
   SUBROUTINE adqmmmx_calc_nve_constant(md_env, md_ener, para_env)
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(md_ener_type), INTENT(inout)                  :: md_ener
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'adqmmmx_calc_nve_constant',&
                                     routineP = moduleN//':'//routineN
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      TYPE(thermostat_type), POINTER                     :: thermostat_coeff,&
                                                            thermostat_shell
      REAL(dp)                                           :: Elam=0.d0

      NULLIFY (force_env, thermostat_coeff, thermostat_shell)

      CALL get_md_env(md_env, force_env=force_env, thermostat_coeff=thermostat_coeff, &
                      thermostat_shell=thermostat_shell)
      adqmmmx_env => force_env%adqmmmx_env
      adqmmmx_env%md_kin_e = md_ener%ekin
      adqmmmx_env%md_pot_e = md_ener%epot
      adqmmmx_env%md_corr  = md_ener%ekin_shell
      CALL get_thermostat_energies(thermostat_shell,md_ener%thermostat_shell_pot, &
                                   md_ener%thermostat_shell_kin, para_env)
      adqmmmx_env%md_corr  = adqmmmx_env%md_corr + md_ener%thermostat_shell_kin +&
                            md_ener%thermostat_shell_pot
      !MN Check if calculated kin and pot are 0
      IF (adqmmmx_env%kin_e==0.0) THEN
         print *, "kinetic energy is 0"
         adqmmmx_env%kin_e = md_ener%ekin
      END IF
      IF (adqmmmx_env%totpot==0.0) THEN
         print *, "totpot energy is 0"
         adqmmmx_env%totpot = md_ener%epot
      END IF
      !MN Get Lambda energy
      CALL lambda_energy_calc(adqmmmx_env)      
      Elam = adqmmmx_env%Elambda  +adqmmmx_env%Elamcorr + adqmmmx_env%Tot_Elambda 
      !MN Calculate constant quantity
      adqmmmx_env%constant = adqmmmx_env%kin_e + adqmmmx_env%totpot +& 
                             adqmmmx_env%md_corr + Elam

      IF (adqmmmx_env%constant==0.0) adqmmmx_env%constant = md_ener%constant 

   END SUBROUTINE adqmmmx_calc_nve_constant

! **************************************************************************************************
!> \brief Calculation of kinetic, potential and constant energy of NVT
!ensembles
!> \param lambda ...
!> \par History
!>      03.2022
!> \author Marvin Nyenhuis
! **************************************************************************************************
   SUBROUTINE adqmmmx_calc_nvt_constant(md_env, md_ener, para_env)
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(md_ener_type), INTENT(inout)                  :: md_ener
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'adqmmmx_calc_nve_constant',&
                                     routineP = moduleN//':'//routineN
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      TYPE(thermostat_type), POINTER                     :: thermostat_coeff,&
                                                            thermostat_shell,&
                                                            thermostat_part
      REAL(dp)                                           :: Elam=0.d0

      NULLIFY (force_env, thermostat_coeff, thermostat_shell, thermostat_part)

      CALL get_md_env(md_env, force_env=force_env, thermostat_part=thermostat_part, &
                thermostat_coeff=thermostat_coeff, thermostat_shell=thermostat_shell)
      CALL get_thermostat_energies(thermostat_part, md_ener%thermostat_part_pot,&
                                   md_ener%thermostat_part_kin, para_env)

      adqmmmx_env => force_env%adqmmmx_env
      adqmmmx_env%md_kin_e = md_ener%ekin
      adqmmmx_env%md_pot_e = md_ener%epot
      adqmmmx_env%md_corr  = md_ener%ekin_shell +&
                             md_ener%thermostat_part_kin +&
                             md_ener%thermostat_part_pot

      CALL get_thermostat_energies(thermostat_shell, md_ener%thermostat_shell_pot, &
                                   md_ener%thermostat_shell_kin, para_env)

      adqmmmx_env%md_corr  = adqmmmx_env%md_corr + md_ener%thermostat_shell_kin +&
                            md_ener%thermostat_shell_pot
      !MN Check if calculated kin and pot are 0
      IF (adqmmmx_env%kin_e==0.0) THEN
         print *, "kinetic energy is 0"
         adqmmmx_env%kin_e = md_ener%ekin
      END IF
      IF (adqmmmx_env%totpot==0.0) THEN
         print *, "totpot energy is 0"
         adqmmmx_env%totpot = md_ener%epot
      END IF
      !MN Get Lambda energy
      CALL lambda_energy_calc(adqmmmx_env)      
      Elam = adqmmmx_env%Elambda  +adqmmmx_env%Elamcorr + adqmmmx_env%Tot_Elambda 
      !MN Calculate constant quantity
      adqmmmx_env%constant = adqmmmx_env%kin_e + adqmmmx_env%totpot +& 
                             adqmmmx_env%md_corr + Elam
 
      IF (adqmmmx_env%constant==0.0) adqmmmx_env%constant = md_ener%constant 

   END SUBROUTINE adqmmmx_calc_nvt_constant
! **************************************************************************************************
!> \brief Calculation of lambda_mass
!> \param lambda ...
!> \par History
!>      03.2022
!> \author Marvin Nyenhuis
! **************************************************************************************************
   SUBROUTINE lambda_energy_calc(adqmmmx_env)
      TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      REAL(dp)                                           :: addElam

      addElam = 0.5*SIGN(1.d0,adqmmmx_env%dpot)*(adqmmmx_env%lambda_mass)*(adqmmmx_env%lambda_velocity**2)
      adqmmmx_env%Elambda = addElam
      IF (adqmmmx_env%dpot_sgn/=SIGN(1.d0, adqmmmx_env%dpot)) THEN 
         adqmmmx_env%Elamcorr = adqmmmx_env%Elamcorr - 2*addElam
         adqmmmx_env%dpot_sgn = SIGN(1.d0, adqmmmx_env%dpot)
         CPWARN("SIGN CHANGED")
      END IF

   END SUBROUTINE lambda_energy_calc


END MODULE adqmmmx_energies
