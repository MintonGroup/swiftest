!**********************************************************************************************************************************
!
!  Unit Name   : whm_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
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
!                whm_pl1P     : pointer to head of WHM planet structure linked-list
!                whm_tp1P     : pointer to head of active WHM test particle structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                whm_pl1P     : pointer to head of WHM planet structure linked-list
!                whm_tp1P     : pointer to head of active WHM test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, j2rp2, j4rp4, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine step_kdk.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, j2rp2, j4rp4, dt)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_step
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lextra_force
     LOGICAL(LGT), INTENT(INOUT) :: lfirst
     INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
     TYPE(whm_pl), POINTER       :: whm_pl1P
     TYPE(whm_tp), POINTER       :: whm_tp1P

! Internals
     LOGICAL(LGT)                                 :: lfirsttp
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xbeg, xend
     TYPE(whm_pl), POINTER                        :: whm_plP

! Executable code
     lfirsttp = lfirst
     IF (ntp > 0) THEN
          IF (lmalloc) THEN
               ALLOCATE(xbeg(NDIM, nplmax), xend(NDIM, nplmax))
               lmalloc = .FALSE.
          END IF
          whm_plP => whm_pl1P
          DO i = 2, npl
               whm_plP => whm_plP%nextP
               xbeg(:, i) = whm_plP%swifter%xh(:)
          END DO
     END IF
     CALL whm_step_pl(lfirst, lextra_force, t, npl, nplmax, whm_pl1P, j2rp2, j4rp4, dt)
     IF (ntp > 0) THEN
          whm_plP => whm_pl1P
          DO i = 2, npl
               whm_plP => whm_plP%nextP
               xend(:, i) = whm_plP%swifter%xh(:)
          END DO
          CALL whm_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xbeg, xend, j2rp2, j4rp4, dt)
     END IF

     RETURN

END SUBROUTINE whm_step
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
