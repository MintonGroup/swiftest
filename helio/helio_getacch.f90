!**********************************************************************************************************************************
!
!  Unit Name   : helio_getacch
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Compute heliocentric accelerations of planets
!
!  Input
!    Arguments : lflag        : logical flag indicating whether to recompute direct cross term accelerations
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                helio_pl1P   : pointer to head of helio planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_pl1P   : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_getacch(lflag, lextra_force, t, npl, nplmax, helio_pl1P, j2rp2, j4rp4)
!
!  Notes       : Adapted from Hal Levison's Swift routine helio_getacch.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_getacch(lflag, lextra_force, t, npl, nplmax, helio_pl1P, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_getacch
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lflag, lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, nplmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
     TYPE(helio_pl), POINTER  :: helio_pl1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP)                                     :: r2, fac
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P, swifter_plP
     TYPE(helio_pl), POINTER                      :: helio_plP

! Executable code
     IF (lflag) THEN
          helio_plP => helio_pl1P
          DO i = 2, npl
               helio_plP => helio_plP%nextP
               helio_plP%ahi(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          END DO
          CALL helio_getacch_int(npl, helio_pl1P)
     END IF
     IF (j2rp2 /= 0.0_DP) THEN
          swifter_pl1P => helio_pl1P%swifter
          IF (lmalloc) THEN
               ALLOCATE(xh(NDIM, nplmax), aobl(NDIM, nplmax), irh(nplmax))
               lmalloc = .FALSE.
          END IF
          swifter_plP => swifter_pl1P
          DO i = 2, npl
               swifter_plP => swifter_plP%nextP
               xh(:, i) = swifter_plP%xh(:)
               r2 = DOT_PRODUCT(xh(:, i), xh(:, i))
               irh(i) = 1.0_DP/SQRT(r2)
          END DO
          CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
          helio_plP => helio_pl1P
          DO i = 2, npl
               helio_plP => helio_plP%nextP
               helio_plP%ah(:) = helio_plP%ahi(:) + aobl(:, i) - aobl(:, 1)
          END DO
     ELSE
          helio_plP => helio_pl1P
          DO i = 2, npl
               helio_plP => helio_plP%nextP
               helio_plP%ah(:) = helio_plP%ahi(:)
          END DO
     END IF
     IF (lextra_force) CALL helio_user_getacch(t, npl, helio_pl1P)

     RETURN

END SUBROUTINE helio_getacch
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
