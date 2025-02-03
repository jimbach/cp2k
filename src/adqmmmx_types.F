!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Basic container type for adaptive QM/MM mixing.
!> \note adjusted names (qmmmx->adqmmmx), added lambda values.
!> \author Ole Schuett, modified by Jim Bachmann and Marvin Nyenhuis
! **************************************************************************************************
MODULE adqmmmx_types
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE kinds,                           ONLY: dp
   USE qmmm_types,                      ONLY: qmmm_env_get,&
                                              qmmm_env_release,&
                                              qmmm_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'adqmmmx_types'

   PUBLIC :: adqmmmx_env_type, adqmmmx_env_get, adqmmmx_env_release
! JB Mod , adqmmmx_env_retain deleted


   TYPE adqmmmx_env_type
      INTEGER                                                 :: ref_count = 1,&
                                                                 c_n_sw=-1

      !MN Core -> Old system, Ext -> New system
      TYPE(qmmm_env_type), POINTER                            :: oldsys => Null()
      TYPE(qmmm_env_type), POINTER                            :: newsys => Null()
!     TYPE(adqmmmx_energy_type), POINTER                      :: 
      REAL(KIND=dp)                                           :: lambda_adqmmmx = 0.0      
      REAL(KIND=dp)                                           :: prop_time = 0.0
      REAL(KIND=dp)                                           :: lambda_velocity = 0.0
      REAL(KIND=dp)                                           :: lambda_mass = 0.0
      REAL(KIND=dp)                                           :: omega = 0.0
      REAL(KIND=dp)                                           :: dpot
      REAL(KIND=dp)                                           :: constant=0.0
      REAL(KIND=dp)                                           :: totpot=0.0
      REAL(KIND=dp)                                           :: kin_e=0.0
      REAL(KIND=dp)                                           :: md_kin_e
      REAL(KIND=dp)                                           :: md_pot_e
      REAL(KIND=dp)                                           :: md_corr=0.0
      REAL(KIND=dp)                                           :: lambda_force=0.0
      REAL(KIND=dp)                                           :: Elambda=0.0      !MN Energy of current lambda
      REAL(KIND=dp)                                           :: Elamcorr=0.0      !MN Sign correction
      REAL(KIND=dp)                                           :: dpot_sgn=0.0      !MN Original sign
      REAL(KIND=dp)                                           :: Tot_Elambda=0.0  !MN Energy of all lambdas
      LOGICAL                                                 :: switch_active=.TRUE.,&
                                                                 labels_changed=.TRUE.,&
                                                                 energy_add=.TRUE.,&
                                                                 output_energies=.FALSE.,&
                                                                 do_discrete,&
                                                                 include_nc=.FALSE. 
   END TYPE adqmmmx_env_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \note Adjustes names (qmmmx->adqmmmx)
!> \param adqmmmx_env ...
!> \param subsys ...
!> \param potential_energy ...
!> \param kinetic_energy ...
!> \author MN
! **************************************************************************************************
   SUBROUTINE adqmmmx_env_get(adqmmmx_env, subsys, potential_energy, kinetic_energy)
      TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: potential_energy, kinetic_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_env_get', routineP = moduleN//':'//routineN

      ! JB Mod 3 02 2025
      !CPASSERT(ASSOCIATED(adqmmmx_env))
      !CPASSERT(adqmmmx_env%ref_count > 0)
      !MN Get Newsys and oldsys
      CALL qmmm_env_get(adqmmmx_env%oldsys, &
                        kinetic_energy=kinetic_energy, &
                        potential_energy=potential_energy, &
                        subsys=subsys)
      CALL qmmm_env_get(adqmmmx_env%newsys, &
                        kinetic_energy=kinetic_energy, &
                        potential_energy=potential_energy, &
                        subsys=subsys)

   END SUBROUTINE adqmmmx_env_get

! JB Mod 3 02 2025 this is not needed any more
! **************************************************************************************************
!> \brief ...
!> \note Adjustes names (qmmmx->adqmmmx)
!> \param qmmmx_env ...
!> \author MN
! **************************************************************************************************
!   SUBROUTINE adqmmmx_env_retain(adqmmmx_env)
!      TYPE(adqmmmx_env_type), POINTER                      :: adqmmmx_env
!
!      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_env_retain', &
!         routineP = moduleN//':'//routineN

!      CPASSERT(ASSOCIATED(adqmmmx_env))
!      CPASSERT(adqmmmx_env%ref_count > 0)
!      adqmmmx_env%ref_count = adqmmmx_env%ref_count + 1
!   END SUBROUTINE adqmmmx_env_retain

! **************************************************************************************************
!> \brief releases the given adqmmmx_env
!> \note Adjustes names (qmmmx->adqmmmx)
!> \param adqmmmx_env the object to release
!> \author MN
! **************************************************************************************************
   SUBROUTINE adqmmmx_env_release(adqmmmx_env)
      TYPE(adqmmmx_env_type), POINTER                      :: adqmmmx_env
!      TYPE(adqmmmx_env_type), INTENT(INOUT)                 :: adqmmmx_env

      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_env_release', &
         routineP = moduleN//':'//routineN

!      CALL qmmm_env_release(adqmmmx_env%oldsys)
!      DEALLOCATE (adqmmmx_env%oldsys)
!      CALL qmmm_env_release(adqmmmx_env%newsys)
!      DEALLOCATE (adqmmmx_env%newsys)
       IF (ASSOCIATED(adqmmmx_env)) THEN
            CPASSERT(adqmmmx_env%ref_count > 0)
            adqmmmx_env%ref_count = adqmmmx_env%ref_count - 1
            IF (adqmmmx_env%ref_count == 0) THEN
              CALL qmmm_env_release(adqmmmx_env%oldsys)
              CALL qmmm_env_release(adqmmmx_env%newsys)
              DEALLOCATE (adqmmmx_env%newsys)
              DEALLOCATE (adqmmmx_env%oldsys)
            END IF
       END IF
       NULLIFY (adqmmmx_env)

  ! this was working
  !    IF (ASSOCIATED(adqmmmx_env)) THEN
  !       CPASSERT(adqmmmx_env%ref_count > 0)
  !       adqmmmx_env%ref_count = adqmmmx_env%ref_count - 1
  !       IF (adqmmmx_env%ref_count == 0) THEN
  !          CALL qmmm_env_release(adqmmmx_env%oldsys)
  !          CALL qmmm_env_release(adqmmmx_env%newsys)
  !          DEALLOCATE (adqmmmx_env)
  !       END IF
  !    END IF
  !    NULLIFY (adqmmmx_env)


   END SUBROUTINE adqmmmx_env_release

END MODULE adqmmmx_types
