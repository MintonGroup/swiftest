!**********************************************************************************************************************************
!
!  Unit Name   : ra15_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : ra15
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
!                ra15_pl1P    : pointer to head of RA15 planet structure linked-list
!                ra15_tp1P    : pointer to head of active RA15 test particle structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : eps          : local truncation error control parameter (maximum size of the last terms in the series expansions
!                                                                         of the dependent variables)
!    File      : none
!
!  Output
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                ra15_pl1P    : pointer to head of RA15 planet structure linked-list
!                ra15_tp1P    : pointer to head of active RA15 test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL ra15_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, ra15_pl1P, ra15_tp1P, j2rp2, j4rp4, dt)
!
!  Notes       : Adapted from Edgar Everhart's RADAU15 routine ra15 (reference given in main program notes)
!
!**********************************************************************************************************************************
SUBROUTINE ra15_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, ra15_pl1P, ra15_tp1P, j2rp2, j4rp4, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_ra15
     USE module_interfaces, EXCEPT_THIS_ONE => ra15_step
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lextra_force
     LOGICAL(LGT), INTENT(INOUT) :: lfirst
     INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
     TYPE(ra15_pl), POINTER      :: ra15_pl1P
     TYPE(ra15_tp), POINTER      :: ra15_tp1P

! Internals
     INTEGER(I4B)                :: k, l, la, lb, lc, ld, le, n
     INTEGER(I4B), SAVE          :: niter
     REAL(DP)                    :: ww, x, hdid, hnext, msys
     REAL(DP), SAVE              :: hh
     TYPE(swifter_pl), POINTER   :: swifter_pl1P, swifter_plP
     TYPE(swifter_tp), POINTER   :: swifter_tp1P, swifter_tpP
     TYPE(ra15_pl), POINTER      :: ra15_plP
     TYPE(ra15_tp), POINTER      :: ra15_tpP

! Executable code
     swifter_pl1P => ra15_pl1P%swifter
     IF (ntp > 0) swifter_tp1P => ra15_tp1P%swifter
     IF (lfirst) THEN
          WRITE(*, 100, ADVANCE = "NO") "Enter the value of eps: "
 100      FORMAT(A)
          READ(*, *) eps
          WRITE(*, *) " eps = ", eps
          CALL coord_h2b(npl, swifter_pl1P, msys)
          IF (ntp > 0) CALL coord_h2b_tp(ntp, swifter_tp1P, swifter_pl1P)
          DO n = 2, 8
               ww = n + n*n
               w(n-1) = ONE/ww
               ww = n
               u(n-1) = ONE/ww
          END DO
          ra15_plP => ra15_pl1P
          DO n = 1, npl
               DO k = 1, NDIM
                    ra15_plP%bd(:, k) = (/ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP /)
                    ra15_plP%b(:, k) = (/ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP /)
               END DO
               ra15_plP => ra15_plP%nextP
          END DO
          ra15_tpP => ra15_tp1P
          DO n = 1, ntp
               DO k = 1, NDIM
                    ra15_tpP%bd(:, k) = (/ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP /)
                    ra15_tpP%b(:, k) = (/ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP /)
               END DO
               ra15_tpP => ra15_tpP%nextP
          END DO
          w1 = HALF
          c(1) = -H(2)
          d(1) = H(2)
          r(1) = ONE/(H(3) - H(2))
          la = 1
          lc = 1
          DO k = 3, 7
               lb = la
               la = lc + 1
               lc = NW(k+1)
               c(la) = -H(k)*c(lb)
               c(lc) = c(la-1) - H(k)
               d(la) = H(2)*d(lb)
               d(lc) = -c(lc)
               r(la) = ONE/(H(k+1) - H(2))
               r(lc) = ONE/(H(k+1) - H(k))
               IF (k == 3) CYCLE
               DO l = 4, k
                    ld = la + l - 3
                    le = lb + l - 4
                    c(ld) = c(le) - H(k)*c(le+1)
                    d(ld) = d(le) + H(l-1)*d(le+1)
                    r(ld) = ONE/(H(k+1) - H(l-1))
               END DO
          END DO
          hh = dt
          niter = 6
          lfirst = .FALSE.
     END IF
     x = t
     DO WHILE ((ABS(x - t - dt)/dt) > DELTARA15)
          CALL ra15_getaccb(lextra_force, t, npl, nplmax, ra15_pl1P, j2rp2, j4rp4)
          IF (ntp > 0) CALL ra15_getaccb_tp(lextra_force, t, npl, ntp, ntpmax, ra15_pl1P, ra15_tp1P, j2rp2, j4rp4)
          ra15_plP => ra15_pl1P
          DO n = 1, npl
               swifter_plP => ra15_plP%swifter
               ra15_plP%xbsav(:) = swifter_plP%xb(:)
               ra15_plP%vbsav(:) = swifter_plP%vb(:)
               ra15_plP%absav(:) = ra15_plP%ab(:)
               ra15_plP => ra15_plP%nextP
          END DO
          ra15_tpP => ra15_tp1P
          DO n = 1, ntp
               swifter_tpP => ra15_tpP%swifter
               ra15_tpP%xbsav(:) = swifter_tpP%xb(:)
               ra15_tpP%vbsav(:) = swifter_tpP%vb(:)
               ra15_tpP%absav(:) = ra15_tpP%ab(:)
               ra15_tpP => ra15_tpP%nextP
          END DO
          IF ((x + hh - t - dt)*(x + hh - t) > 0.0_DP) hh = t + dt - x
          CALL ra15_sequence(niter, npl, nplmax, ntp, ntpmax, ra15_pl1P, ra15_tp1P, x, hh, hdid, hnext, j2rp2, j4rp4, lextra_force)
          hh = hnext
          niter = 2
     END DO
     CALL coord_b2h(npl, swifter_pl1P)
     CALL coord_b2h_tp(ntp, swifter_tp1P, swifter_pl1P)

     RETURN

END SUBROUTINE ra15_step
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
