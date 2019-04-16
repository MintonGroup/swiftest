!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_step_in_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Step all active test particles closely encountering a given planet ahead one substep in inner integration region
!
!  Input
!    Arguments : index        : inner substep number within current set
!                lfirst       : logical flag indicating whether current invocation is the first
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                nenc         : number of test particles encountering current planet
!                ntpmax       : maximum allowed number of test particles
!                rmvs_pl1P    : pointer to head of RMVS planet structure linked-list
!                rmvs_pleP    : pointer to RMVS planet structure of planet being closely encountered
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                rmvs_pleP    : pointer to RMVS planet structure of planet being closely encountered
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_step_in_tp(index, lfirst, lextra_force, t, npl, nplmax, nenc, ntpmax, rmvs_pl1P, rmvs_pleP, j2rp2,
!                                     j4rp4, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine rmvs_step_in_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_step_in_tp(index, lfirst, lextra_force, t, npl, nplmax, nenc, ntpmax, rmvs_pl1P, rmvs_pleP, j2rp2, j4rp4, dt)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_step_in_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lextra_force
     LOGICAL(LGT), INTENT(INOUT) :: lfirst
     INTEGER(I4B), INTENT(IN)    :: index, npl, nplmax, nenc, ntpmax
     REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
     TYPE(rmvs_pl), POINTER      :: rmvs_pl1P, rmvs_pleP

! Internals
     INTEGER(I4B)           :: i
     REAL(DP)               :: dth, mu
     TYPE(rmvs_tp), POINTER :: rmvs_tpenc1P, rmvs_tpP

! Executable code
     dth = 0.5_DP*dt
     mu = rmvs_pleP%whm%swifter%mass
     rmvs_tpenc1P => rmvs_pleP%tpenc1P
     IF (lfirst) THEN
          CALL rmvs_getaccp_tp(index-1, lextra_force, t, npl, nplmax, nenc, ntpmax, rmvs_pl1P, rmvs_pleP, j2rp2, j4rp4)
          lfirst = .FALSE.
     END IF
     CALL rmvs_kickvp_tp(nenc, rmvs_tpenc1P, dth)
     CALL rmvs_drift_tp(nenc, rmvs_tpenc1P, mu, dt)
     rmvs_tpP => rmvs_tpenc1P
     DO i = 1, nenc
          rmvs_tpP%whm%swifter%xh(:) = rmvs_pleP%xin(:, index) + rmvs_tpP%xpc(:)
          rmvs_tpP => rmvs_tpP%tpencP
     END DO
     CALL rmvs_getaccp_tp(index, lextra_force, t+dt, npl, nplmax, nenc, ntpmax, rmvs_pl1P, rmvs_pleP, j2rp2, j4rp4)
     CALL rmvs_kickvp_tp(nenc, rmvs_tpenc1P, dth)
     rmvs_tpP => rmvs_tpenc1P
     DO i = 1, nenc
          rmvs_tpP%whm%swifter%vh(:) = rmvs_pleP%vin(:, index) + rmvs_tpP%vpc(:)
          rmvs_tpP => rmvs_tpP%tpencP
     END DO

     RETURN

END SUBROUTINE rmvs_step_in_tp
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
