!**********************************************************************************************************************************
!
!  Unit Name   : discard_sun
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : discard
!  Language    : Fortran 90/95
!
!  Description : Check to see if test particles should be discarded based on their positions relative to the Sun
!                or because they are unbound from the system
!
!  Input
!    Arguments : t            : time
!                ntp          : number of active test particles
!                msys         : total system mass
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!                rmin         : minimum heliocentric radius for test particle
!                rmax         : maximum heliocentric radius for test particle
!                rmaxu        : maximum unbound heliocentric radius for test particle
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : status messages
!    File      : none
!
!  Invocation  : CALL discard_sun(t, ntp, msys, swifter_tp1P, rmin, rmax, rmaxu)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_sun.f
!
!**********************************************************************************************************************************
SUBROUTINE discard_sun(t, ntp, msys, swifter_tp1P, rmin, rmax, rmaxu)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => discard_sun
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: ntp
     REAL(DP), INTENT(IN)      :: t, msys, rmin, rmax, rmaxu
     TYPE(swifter_tp), POINTER :: swifter_tp1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     rmin2 = rmin*rmin
     rmax2 = rmax*rmax
     rmaxu2 = rmaxu*rmaxu
     swifter_tpP => swifter_tp1P
     DO i = 1, ntp
          IF (swifter_tpP%status == ACTIVE) THEN
               rh2 = DOT_PRODUCT(swifter_tpP%xh(:), swifter_tpP%xh(:))
               IF ((rmax >= 0.0_DP) .AND. (rh2 > rmax2)) THEN
                    swifter_tpP%status = DISCARDED_RMAX
                    WRITE(*, *) "Particle ", swifter_tpP%id, " too far from Sun at t = ", t
               ELSE IF ((rmin >= 0.0_DP) .AND. (rh2 < rmin2)) THEN
                    swifter_tpP%status = DISCARDED_RMIN
                    WRITE(*, *) "Particle ", swifter_tpP%id, " too close to Sun at t = ", t
               ELSE IF (rmaxu >= 0.0_DP) THEN
                    rb2 = DOT_PRODUCT(swifter_tpP%xb(:), swifter_tpP%xb(:))
                    vb2 = DOT_PRODUCT(swifter_tpP%vb(:), swifter_tpP%vb(:))
                    energy = 0.5_DP*vb2 - msys/SQRT(rb2)
                    IF ((energy > 0.0_DP) .AND. (rb2 > rmaxu2)) THEN
                         swifter_tpP%status = DISCARDED_RMAXU
                         WRITE(*, *) "Particle ", swifter_tpP%id, " is unbound and too far from barycenter at t = ", t
                    END IF
               END IF
          END IF
          swifter_tpP => swifter_tpP%nextP
     END DO

     RETURN

END SUBROUTINE discard_sun
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
