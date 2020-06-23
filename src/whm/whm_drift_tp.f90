!**********************************************************************************************************************************
!
!  Unit Name   : whm_drift_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Loop through test particles and call Danby drift routine
!
!  Input
!    Arguments : ntp      : number of active test particles
!                whm_tp1P : pointer to head of active WHM test particle structure linked-list
!                mu       : mass of the Sun
!                dt       : time step
!                c2       : inverse speed of light squared
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_tp1P : pointer to head of active WHM test particle structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL whm_drift_tp(ntp, whm_tp1P, mu, dt, c2)
!
!  Notes       : Adapted from Hal Levison's Swift routine drift_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_drift_tp(ntp, whm_tp1P, mu, dt, c2)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_drift_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: mu, dt, c2
     TYPE(whm_tp), POINTER    :: whm_tp1P

! Internals
     INTEGER(I4B)              :: i, iflag
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(whm_tp), POINTER     :: whm_tpP
     REAL(DP)                  :: dtp, energy, vmag2, rmag

! Executable code
     whm_tpP => whm_tp1P
     DO i = 1, ntp
          swifter_tpP => whm_tpP%swifter
          IF (swifter_tpP%status == ACTIVE) THEN
               rmag = SQRT(DOT_PRODUCT(swifter_tpP%xh(:), swifter_tpP%xh(:)))
               vmag2 = DOT_PRODUCT(swifter_tpP%vh(:), swifter_tpP%vh(:))
               energy = 0.5_DP*vmag2 - mu/rmag
               dtp = dt * (1.0_DP + 3 * c2 * energy)
               CALL drift_one(mu, swifter_tpP%xh(:), swifter_tpP%vh(:), dt, iflag)
               IF (iflag /= 0) THEN
                    swifter_tpP%status = DISCARDED_DRIFTERR
                    WRITE(*, *) "Particle ", swifter_tpP%id, " lost due to error in Danby drift"
               END IF
          END IF
          whm_tpP => whm_tpP%nextP
     END DO

     RETURN

END SUBROUTINE whm_drift_tp
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
