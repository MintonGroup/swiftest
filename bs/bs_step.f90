!**********************************************************************************************************************************
!
!  Unit Name   : bs_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
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
!                bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : eps          : local truncation error control parameter
!    File      : none
!
!  Output
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL bs_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine bs_step.f
!
!**********************************************************************************************************************************
SUBROUTINE bs_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_step
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lextra_force
     LOGICAL(LGT), INTENT(INOUT) :: lfirst
     INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
     TYPE(bs_pl), POINTER        :: bs_pl1P
     TYPE(bs_tp), POINTER        :: bs_tp1P

! Internals
     INTEGER(I4B)                :: i
     REAL(DP)                    :: x, h, hdid, hnext, msys
     REAL(DP), SAVE              :: eps
     TYPE(swifter_pl), POINTER   :: swifter_pl1P, swifter_plP
     TYPE(swifter_tp), POINTER   :: swifter_tp1P, swifter_tpP
     TYPE(bs_pl), POINTER        :: bs_plP
     TYPE(bs_tp), POINTER        :: bs_tpP

! Executable code
     swifter_pl1P => bs_pl1P%swifter
     swifter_tp1P => bs_tp1P%swifter
     IF (lfirst) THEN
          WRITE(*, 100, ADVANCE = "NO") "Enter the value of eps: "
 100      FORMAT(A)
          READ(*, *) eps
          WRITE(*, *) " eps = ", eps
          CALL coord_h2b(npl, swifter_pl1P, msys)
          CALL coord_h2b_tp(ntp, swifter_tp1P, swifter_pl1P)
          bs_plP => bs_pl1P
          DO i = 1, npl
               swifter_plP => bs_plP%swifter
               bs_plP%y(1:NDIM) = swifter_plP%xb(:)
               bs_plP%y(NDIM+1:NDIM2) = swifter_plP%vb(:)
               bs_plP => bs_plP%nextP
          END DO
          bs_tpP => bs_tp1P
          DO i = 1, ntp
               swifter_tpP => bs_tpP%swifter
               bs_tpP%y(1:NDIM) = swifter_tpP%xb(:)
               bs_tpP%y(NDIM+1:NDIM2) = swifter_tpP%vb(:)
               bs_tpP => bs_tpP%nextP
          END DO
          lfirst = .FALSE.
     END IF
     x = t
     h = dt
     DO WHILE ((ABS(x - t - dt)/dt) > DELTABS)
          CALL bs_derivs(lextra_force, x, npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)
          bs_plP => bs_pl1P
          DO i = 1, npl
               bs_plP%yscal(:) = ABS(bs_plP%y(:)) + ABS(h*bs_plP%dydx(:)) + TINYBS
               bs_plP => bs_plP%nextP
          END DO
          bs_tpP => bs_tp1P
          DO i = 1, ntp
               bs_tpP%yscal(:) = ABS(bs_tpP%y(:)) + ABS(h*bs_tpP%dydx(:)) + TINYBS
               bs_tpP => bs_tpP%nextP
          END DO
          IF ((x + h - t - dt)*(x + h - t) > 0.0_DP) h = t + dt - x
          CALL bs_bsstep(npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, x, h, eps, hdid, hnext, j2rp2, j4rp4, lextra_force, bs_derivs)
          h = hnext
     END DO
     bs_plP => bs_pl1P
     DO i = 1, npl
          swifter_plP => bs_plP%swifter
          swifter_plP%xb(:) = bs_plP%y(1:NDIM)
          swifter_plP%vb(:) = bs_plP%y(NDIM+1:NDIM2)
          bs_plP => bs_plP%nextP
     END DO
     bs_tpP => bs_tp1P
     DO i = 1, ntp
          swifter_tpP => bs_tpP%swifter
          swifter_tpP%xb(:) = bs_tpP%y(1:NDIM)
          swifter_tpP%vb(:) = bs_tpP%y(NDIM+1:NDIM2)
          bs_tpP => bs_tpP%nextP
     END DO
     CALL coord_b2h(npl, swifter_pl1P)
     CALL coord_b2h_tp(ntp, swifter_tp1P, swifter_pl1P)

     RETURN

END SUBROUTINE bs_step
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
