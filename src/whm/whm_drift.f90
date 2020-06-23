!**********************************************************************************************************************************
!
!  Unit Name   : whm_drift
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Loop through planets and call Danby drift routine
!
!  Input
!    Arguments : npl      : number of planets
!                whm_pl1P : pointer to head of WHM planet structure linked-list
!                dt       : time step
!                c2       : inverse speed of light squared
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL whm_drift(npl, whm_pl1P, dt, c2)
!
!  Notes       : Adapted from Hal Levison's Swift routine drift.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_drift(npl, whm_pl1P, dt, c2)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_drift
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: dt, c2
     TYPE(whm_pl), POINTER    :: whm_pl1P

! Internals
     INTEGER(I4B)              :: i, iflag
     REAL(DP)                  :: etajm1, etaj, mu, msun
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(whm_pl), POINTER     :: whm_plP
     REAL(DP)                  :: dtp, energy, vmag2, rmag

! Executable code
     whm_plP => whm_pl1P
     swifter_plP => whm_plP%swifter
     msun = swifter_plP%mass
     etajm1 = msun
     DO i = 2, npl
          whm_plP => whm_plP%nextP
          swifter_plP => whm_plP%swifter
          etaj = etajm1 + swifter_plP%mass
          mu = msun*etaj/etajm1
          rmag = SQRT(DOT_PRODUCT(whm_plP%xj(:), whm_plP%xj(:)))
          vmag2 = DOT_PRODUCT(whm_plP%vj(:), whm_plp%vj(:))
          energy = 0.5_DP*vmag2 - mu/rmag
          dtp = dt * (1.0_DP + 3 * c2 * energy)
          CALL drift_one(mu, whm_plP%xj(:), whm_plP%vj(:), dtp, iflag)
          IF (iflag /= 0) THEN
               WRITE(*, *) " Planet ", swifter_plP%id, " is lost!!!!!!!!!!"
               WRITE(*, *) mu, dt
               WRITE(*, *) whm_plP%xj(:)
               WRITE(*, *) whm_plP%vj(:)
               WRITE(*, *) " STOPPING "
               CALL util_exit(FAILURE)
          END IF
          etajm1 = etaj
     END DO

     RETURN

END SUBROUTINE whm_drift
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
