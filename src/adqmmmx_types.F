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

   PUBLIC :: adqmmmx_env_type, adqmmmx_env_get, adqmmmx_env_retain, adqmmmx_env_release

   TYPE adqmmmx_env_type
      INTEGER                                                 :: ref_count = 1
      !JB Core -> Old system, Ext -> New system
      TYPE(qmmm_env_type), POINTER                            :: core => Null()
      TYPE(qmmm_env_type), POINTER                            :: ext => Null()
      REAL(KIND=dp)                                           :: lambda_adqmmmx = 0.0   !JB ADQMMMX mods
                                                              :: lambda_velocity = 0.0
                                                              :: lambda_mass = 0.0
                                                              :: conserved_quantity = 0.0
                                                              :: switch_t = 0.0			!JB switching time
																															:: elambda = 0.d0
																															:: elambda_not_one = 0.d0
																															:: elambda_sum = 0.d0
			INTEGER																									:: switch_active = 0
   END TYPE adqmmmx_env_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \note Adjustes names (qmmmx->adqmmmx)
!> \param adqmmmx_env ...
!> \param subsys ...
!> \param potential_energy ...
!> \param kinetic_energy ...
!> \author JB
! **************************************************************************************************
   SUBROUTINE adqmmmx_env_get(adqmmmx_env, subsys, potential_energy, kinetic_energy)
      TYPE(adqmmmx_env_type), POINTER                    :: adqmmmx_env
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: potential_energy, kinetic_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_env_get', routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(adqmmmx_env))
      CPASSERT(adqmmmx_env%ref_count > 0)

      ! get the underlying energies from primary subsys.  This is the only subsys
      ! for conventional QM/MM, and force-mixing knows to put relevant energy there.
      CALL qmmm_env_get(adqmmmx_env%ext, &
                        kinetic_energy=kinetic_energy, &
                        potential_energy=potential_energy, &
                        subsys=subsys)

   END SUBROUTINE adqmmmx_env_get

! **************************************************************************************************
!> \brief ...
!> \note Adjustes names (qmmmx->adqmmmx)
!> \param qmmmx_env ...
!> \author JB
! **************************************************************************************************
   SUBROUTINE adqmmmx_env_retain(adqmmmx_env)
      TYPE(adqmmmx_env_type), POINTER                      :: adqmmmx_env

      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_env_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(adqmmmx_env))
      CPASSERT(adqmmmx_env%ref_count > 0)
      adqmmmx_env%ref_count = adqmmmx_env%ref_count + 1
   END SUBROUTINE adqmmmx_env_retain

! **************************************************************************************************
!> \brief releases the given adqmmmx_env
!> \note Adjustes names (qmmmx->adqmmmx)
!> \param adqmmmx_env the object to release
!> \author JB
! **************************************************************************************************
   SUBROUTINE adqmmmx_env_release(adqmmmx_env)
      TYPE(adqmmmx_env_type), POINTER                      :: adqmmmx_env

      CHARACTER(len=*), PARAMETER :: routineN = 'adqmmmx_env_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(adqmmmx_env)) THEN
         CPASSERT(adqmmmx_env%ref_count > 0)
         adqmmmx_env%ref_count = adqmmmx_env%ref_count - 1
         IF (adqmmmx_env%ref_count == 0) THEN
            CALL qmmm_env_release(adqmmmx_env%core)
            CALL qmmm_env_release(adqmmmx_env%ext)
            DEALLOCATE (adqmmmx_env)
         END IF
      END IF
      NULLIFY (adqmmmx_env)
   END SUBROUTINE adqmmmx_env_release

END MODULE adqmmmx_types
