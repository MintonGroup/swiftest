!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_sun_pl
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Check to see if planets should be discarded based on their positions relative to the Sun
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                msys         : total system mass
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                rmin         : minimum allowed heliocentric radius
!                rmax         : maximum allowed heliocentric radius
!                rmaxu        : maximum allowed heliocentric radius for unbound planets
!                ldiscards    : logical flag indicating whether any planets are discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                ldiscards    : logical flag indicating whether any planets are discarded
!    Terminal  : status message
!    File      : none
!
!  Invocation  : CALL symba_discard_sun_pl(t, npl, msys, swifter_pl1P, rmin, rmax, rmaxu, ldiscards)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_massive5.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_sun_pl(t, npl, msys, swifter_pl1P, rmin, rmax, rmaxu, ldiscards)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_sun_pl
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(INOUT) :: ldiscards
     INTEGER(I4B), INTENT(IN)    :: npl
     REAL(DP), INTENT(IN)        :: t, msys, rmin, rmax, rmaxu
     TYPE(swifter_pl), POINTER   :: swifter_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2
     TYPE(swifter_pl), POINTER :: swifter_plP

! Executable code
     rmin2 = rmin*rmin
     rmax2 = rmax*rmax
     rmaxu2 = rmaxu*rmaxu
     swifter_plP => swifter_pl1P
     DO i = 2, npl
          swifter_plP => swifter_plP%nextP
          IF (swifter_plP%status == ACTIVE) THEN
               rh2 = DOT_PRODUCT(swifter_plP%xh(:), swifter_plP%xh(:))
               IF ((rmax >= 0.0_DP) .AND. (rh2 > rmax2)) THEN
                    ldiscards = .TRUE.
                    swifter_plP%status = DISCARDED_RMAX
                    WRITE(*, *) "Particle ", swifter_plP%id, " too far from Sun at t = ", t
               ELSE IF ((rmin >= 0.0_DP) .AND. (rh2 < rmin2)) THEN
                    ldiscards = .TRUE.
                    swifter_plP%status = DISCARDED_RMIN
                    WRITE(*, *) "Particle ", swifter_plP%id, " too close to Sun at t = ", t
               ELSE IF (rmaxu >= 0.0_DP) THEN
                    rb2 = DOT_PRODUCT(swifter_plP%xb(:), swifter_plP%xb(:))
                    vb2 = DOT_PRODUCT(swifter_plP%vb(:), swifter_plP%vb(:))
                    energy = 0.5_DP*vb2 - msys/SQRT(rb2)
                    IF ((energy > 0.0_DP) .AND. (rb2 > rmaxu2)) THEN
                         ldiscards = .TRUE.
                         swifter_plP%status = DISCARDED_RMAXU
                         WRITE(*, *) "Particle ", swifter_plP%id, " is unbound and too far from barycenter at t = ", t
                    END IF
               END IF
          END IF
     END DO

     RETURN

END SUBROUTINE symba_discard_sun_pl
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
