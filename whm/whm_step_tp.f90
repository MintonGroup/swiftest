!**********************************************************************************************************************************
!
!  Unit Name   : whm_step_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Step active test particles ahead using kick-drift-kick algorithm
!
!  Input
!    Arguments : lfirsttp     : logical flag indicating whether current invocation is the first
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                whm_pl1P     : pointer to head of WHM planet structure linked-list
!                whm_tp1P     : pointer to head of active WHM test particle structure linked-list
!                xbeg         : heliocentric planet positions at beginning of time step
!                xend         : heliocentric planet positions at end of time step
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_tp1P     : pointer to head of active WHM test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xbeg, xend, j2rp2,
!                                 j4rp4, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine step_kdk_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xbeg, xend, j2rp2, j4rp4, dt)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_step_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                   :: lfirsttp, lextra_force
     INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4, dt
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xbeg, xend
     TYPE(whm_pl), POINTER                      :: whm_pl1P
     TYPE(whm_tp), POINTER                      :: whm_tp1P

! Internals
     REAL(DP) :: dth

! Executable code
     dth = 0.5_DP*dt
     IF (lfirsttp) CALL whm_getacch_tp(lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xbeg, j2rp2, j4rp4)
     CALL whm_kickvh_tp(ntp, whm_tp1P, dth)
     CALL whm_drift_tp(ntp, whm_tp1P, whm_pl1P%swifter%mass, dt)
     CALL whm_getacch_tp(lextra_force, t+dt, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xend, j2rp2, j4rp4)
     CALL whm_kickvh_tp(ntp, whm_tp1P, dth)

     RETURN

END SUBROUTINE whm_step_tp
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
