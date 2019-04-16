!**********************************************************************************************************************************
!
!  Unit Name   : helio_step_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Step active test particles ahead using Democratic Heliocentric method
!
!  Input
!    Arguments : lfirsttp     : logical flag indicating whether current invocation is the first
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                helio_pl1P   : pointer to head of helio planet structure linked-list
!                helio_tp1P   : pointer to head of active helio test particle structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!                xbeg         : heliocentric planet positions at beginning of time step
!                xend         : heliocentric planet positions at end of time step
!                ptb          : negative barycentric velocity of the Sun at beginning of time step
!                pte          : negative barycentric velocity of the Sun at end of time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirsttp     : logical flag indicating whether current invocation is the first
!                helio_tp1P   : pointer to head of active helio test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, j2rp2, j4rp4, dt,
!                                   xbeg, xend, ptb, pte)
!
!  Notes       : Adapted from Hal Levison's Swift routine helio_step_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, j2rp2, j4rp4, dt, xbeg,     &
     xend, ptb, pte)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_step_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                   :: lextra_force
     LOGICAL(LGT), INTENT(INOUT)                :: lfirsttp
     INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4, dt
     REAL(DP), DIMENSION(NDIM), INTENT(IN)      :: ptb, pte
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xbeg, xend
     TYPE(helio_pl), POINTER                    :: helio_pl1P
     TYPE(helio_tp), POINTER                    :: helio_tp1P

! Internals
     LOGICAL(LGT)              :: lflag
     REAL(DP)                  :: dth, mu
     TYPE(swifter_tp), POINTER :: swifter_tp1P

! Executable code
     dth = 0.5_DP*dt
     lflag = lfirsttp
     mu = helio_pl1P%swifter%mass
     swifter_tp1P => helio_tp1P%swifter
     IF (lfirsttp) THEN
          CALL coord_vh2vb_tp(ntp, swifter_tp1P, -ptb)
          lfirsttp = .FALSE.
     END IF
     CALL helio_lindrift_tp(ntp, swifter_tp1P, dth, ptb)
     CALL helio_getacch_tp(lflag, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, xbeg, j2rp2, j4rp4)
     lflag = .TRUE.
     CALL helio_kickvb_tp(ntp, helio_tp1P, dth)
     CALL helio_drift_tp(ntp, swifter_tp1P, mu, dt)
     CALL helio_getacch_tp(lflag, lextra_force, t+dt, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, xend, j2rp2, j4rp4)
     CALL helio_kickvb_tp(ntp, helio_tp1P, dth)
     CALL helio_lindrift_tp(ntp, swifter_tp1P, dth, pte)
     CALL coord_vb2vh_tp(ntp, swifter_tp1P, -pte)

     RETURN

END SUBROUTINE helio_step_tp
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
