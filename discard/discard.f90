!**********************************************************************************************************************************
!
!  Unit Name   : discard
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : discard
!  Language    : Fortran 90/95
!
!  Description : Check to see if test particles should be discarded based on their positions or because they are unbound from
!                the system
!
!  Input
!    Arguments : t              : time
!                dt             : time step
!                npl            : number of planets
!                ntp            : number of test particles
!                swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P   : pointer to head of active Swifter test particle structure linked-list
!                rmin           : minimum heliocentric radius for test particle
!                rmax           : maximum heliocentric radius for test particle
!                rmaxu          : maximum unbound heliocentric radius for test particle
!                qmin           : minimum pericenter distance for test particle
!                qmin_alo       : minimum semimajor axis for qmin
!                qmin_ahi       : maximum semimajor axis for qmin
!                qmin_coord     : coordinate frame to use for qmin
!                lclose         : logical flag indicating whether to check for planet-test particle encounters
!                lrhill_present : logical flag indicating whether Hill sphere radii for planets are present
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P   : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL discard(t, dt, npl, ntp, swifter_pl1P, swifter_tp1P, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi,
!                             qmin_coord, lclose, lrhill_present)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard.f
!
!**********************************************************************************************************************************
SUBROUTINE discard(t, dt, npl, ntp, swifter_pl1P, swifter_tp1P, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, qmin_coord, lclose,  &
     lrhill_present)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => discard
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)  :: lclose, lrhill_present
     INTEGER(I4B), INTENT(IN)  :: npl, ntp
     REAL(DP), INTENT(IN)      :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
     CHARACTER(*), INTENT(IN)  :: qmin_coord
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tp1P

! Internals
     REAL(DP) :: msys

! Executable code
     IF ((rmin >= 0.0_DP) .OR. (rmax >= 0.0_DP) .OR. (rmaxu >= 0.0_DP) .OR. ((qmin >= 0.0_DP) .AND. (qmin_coord == "BARY"))) THEN
          CALL coord_h2b(npl, swifter_pl1P, msys)
          CALL coord_h2b_tp(ntp, swifter_tp1P, swifter_pl1P)
     END IF
     IF ((rmin >= 0.0_DP) .OR. (rmax >= 0.0_DP) .OR. (rmaxu >= 0.0_DP)) CALL discard_sun(t, ntp, msys, swifter_tp1P, rmin, rmax,  &
          rmaxu)
     IF (qmin >= 0.0_DP) CALL discard_peri(t, npl, ntp, swifter_pl1P, swifter_tp1P, msys, qmin, qmin_alo, qmin_ahi, qmin_coord,   &
          lrhill_present)
     IF (lclose) CALL discard_pl(t, dt, npl, ntp, swifter_pl1P, swifter_tp1P)

     RETURN

END SUBROUTINE discard
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
