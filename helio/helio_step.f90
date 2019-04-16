!**********************************************************************************************************************************
!
!  Unit Name   : helio_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
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
!                helio_pl1P   : pointer to head of WHM planet structure linked-list
!                helio_tp1P   : pointer to head of active WHM test particle structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                helio_pl1P   : pointer to head of WHM planet structure linked-list
!                helio_tp1P   : pointer to head of active WHM test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, j2rp2, j4rp4, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine helio_step.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, j2rp2, j4rp4, dt)

! Modules
     USE module_parameters
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_step
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lextra_force
     LOGICAL(LGT), INTENT(INOUT) :: lfirst
     INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
     TYPE(helio_pl), POINTER     :: helio_pl1P
     TYPE(helio_tp), POINTER     :: helio_tp1P

! Internals
     LOGICAL(LGT)                                 :: lfirsttp
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     REAL(DP), DIMENSION(NDIM)                    :: ptb, pte
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xbeg, xend

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(xbeg(NDIM, nplmax), xend(NDIM, nplmax))
          lmalloc = .FALSE.
     END IF
     lfirsttp = lfirst
     CALL helio_step_pl(lfirst, lextra_force, t, npl, nplmax, helio_pl1P, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)
     IF (ntp > 0) CALL helio_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, j2rp2, j4rp4,   &
          dt, xbeg, xend, ptb, pte)

     RETURN

END SUBROUTINE helio_step
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
