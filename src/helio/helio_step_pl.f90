!**********************************************************************************************************************************
!
!  Unit Name   : helio_step_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Step planets ahead Democratic Heliocentric method
!
!  Input
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                helio_pl1P   : pointer to head of helio planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                helio_pl1P   : pointer to head of helio planet structure linked-list
!                xbeg         : heliocentric planet positions at beginning of time step
!                xend         : heliocentric planet positions at end of time step
!                ptb          : negative barycentric velocity of the Sun at beginning of time step
!                pte          : negative barycentric velocity of the Sun at end of time step
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_step_pl(lfirst, lextra_force, t, npl, nplmax, helio_pl1P, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)
!
!  Notes       : Adapted from Hal Levison's Swift routine helio_step_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_step_pl(lfirst, lextra_force, t, npl, nplmax, helio_plA, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)

! Modules
     USE swiftest
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_step_pl
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                    :: lextra_force
     LOGICAL(LGT), INTENT(INOUT)                 :: lfirst
     INTEGER(I4B), INTENT(IN)                    :: npl, nplmax
     REAL(DP), INTENT(IN)                        :: t, j2rp2, j4rp4, dt
     REAL(DP), DIMENSION(NDIM), INTENT(OUT)      :: ptb, pte
     REAL(DP), DIMENSION(NDIM, npl), INTENT(OUT) :: xbeg, xend
     TYPE(helio_pl), INTENT(INOUT)               :: helio_plA

! Internals
     LOGICAL(LGT)              :: lflag
     INTEGER(I4B)              :: i
     REAL(DP)                  :: dth, msys

! Executable code
     dth = 0.5_DP*dt
     lflag = lfirst
     IF (lfirst) THEN
          CALL coord_vh2vb(npl, helio_plA%swiftest, msys)
          lfirst = .FALSE.
     END IF
     CALL helio_lindrift(npl, helio_plA%swiftest, dth, ptb)
     CALL helio_getacch(lflag, lextra_force, t, npl, nplmax, helio_plA, j2rp2, j4rp4)
     lflag = .TRUE.
     CALL helio_kickvb(npl, helio_plA, dth)
     DO i = 2, npl
          xbeg(:, i) = helio_plA%swiftest%xh(:,i)
     END DO
     CALL helio_drift(npl, helio_plA%swiftest, dt)
     DO i = 2, npl
          xend(:, i) = helio_plA%swiftest%xh(:,i)
     END DO
     CALL helio_getacch(lflag, lextra_force, t+dt, npl, nplmax, helio_plA, j2rp2, j4rp4)
     CALL helio_kickvb(npl, helio_plA, dth)
     CALL helio_lindrift(npl, helio_plA%swiftest, dth, pte)
     CALL coord_vb2vh(npl, helio_plA%swiftest)

     RETURN

END SUBROUTINE helio_step_pl
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
