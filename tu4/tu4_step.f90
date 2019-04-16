!**********************************************************************************************************************************
!
!  Unit Name   : tu4_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : tu4
!  Language    : Fortran 90/95
!
!  Description : Step planets and active test particles ahead in heliocentric coordinates
!
!  Input
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                tu4_pl1P     : pointer to head of TU4 planet structure linked-list
!                tu4_tp1P     : pointer to head of active TU4 test particle structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                tu4_pl1P     : pointer to head of TU4 planet structure linked-list
!                tu4_tp1P     : pointer to head of active TU4 test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL tu4_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, tu4_pl1P, tu4_tp1P, j2rp2, j4rp4, dt)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine tu4_step.f
!
!**********************************************************************************************************************************
SUBROUTINE tu4_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, tu4_pl1P, tu4_tp1P, j2rp2, j4rp4, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_tu4
     USE module_interfaces, EXCEPT_THIS_ONE => tu4_step
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lextra_force
     LOGICAL(LGT), INTENT(INOUT) :: lfirst
     INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
     TYPE(tu4_pl), POINTER       :: tu4_pl1P
     TYPE(tu4_tp), POINTER       :: tu4_tp1P

! Internals
     INTEGER(I4B)              :: j
     REAL(DP)                  :: msys, adt, bdt, tm
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tp1P
     TYPE(tu4_pl), POINTER     :: tu4_plP

! Executable code
     IF (lfirst) THEN
          w1 = 1.0_DP/(2.0_DP - (2.0_DP**(1.0_DP/3.0_DP)))
          w0 = 1.0_DP - 2.0_DP*w1
          asymp(1) = 0.5_DP*w1
          asymp(2) = 0.5_DP*(w0 + w1)
          asymp(3) = asymp(2)
          asymp(4) = asymp(1)
          bsymp(1) = 0.0_DP
          bsymp(2) = w1
          bsymp(3) = w0
          bsymp(4) = w1
          lfirst = .FALSE.
     END IF
     swifter_pl1P => tu4_pl1P%swifter
     CALL coord_h2b(npl, swifter_pl1P, msys)
     IF (ntp > 0) THEN
          swifter_tp1P => tu4_tp1P%swifter
          CALL coord_h2b_tp(ntp, swifter_tp1P, swifter_pl1P)
     END IF
     tm = t
     adt = asymp(1)*dt
     CALL tu4_ldrift(npl, ntp, swifter_pl1P, swifter_tp1P, adt)
     tm = tm + adt
     DO j = 2, 4
          CALL tu4_getaccb(lextra_force, tm, npl, nplmax, tu4_pl1P, j2rp2, j4rp4)
          IF (ntp > 0) CALL tu4_getaccb_tp(lextra_force, tm, npl, ntp, ntpmax, tu4_pl1P, tu4_tp1P, j2rp2, j4rp4)
          bdt = bsymp(j)*dt
          CALL tu4_kickvb(npl, ntp, tu4_pl1P, tu4_tp1P, bdt)
          adt = asymp(j)*dt
          CALL tu4_ldrift(npl, ntp, swifter_pl1P, swifter_tp1P, adt)
          tm = tm + adt
     END DO
     CALL coord_b2h(npl, swifter_pl1P)
     IF (ntp > 0) CALL coord_b2h_tp(ntp, swifter_tp1P, swifter_pl1P)

     RETURN

END SUBROUTINE tu4_step
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
