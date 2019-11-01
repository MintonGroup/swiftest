!**********************************************************************************************************************************
!
!  Unit Name   : symba_helio_getacch
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Compute heliocentric accelerations of planets
!
!  Input
!    Arguments : lflag        : logical flag indicating whether to recompute direct cross term heliocentric accelerations
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplm         : number of planets with mass > mtiny
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
!  Invocation  : CALL symba_helio_getacch(lflag, lextra_force, t, npl, nplm, nplmax, helio_pl1P, j2rp2, j4rp4)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_helio_getacch.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_helio_getacch(lflag, lextra_force, t, npl, nplm, nplmax, helio_plA, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_helio_getacch
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                     :: lflag, lextra_force
     INTEGER(I4B), INTENT(IN)                     :: npl, nplm, nplmax
     REAL(DP), INTENT(IN)                         :: t, j2rp2, j4rp4
     TYPE(helio_pl), INTENT(INOUT)                :: helio_plA

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP)                                     :: r2, fac
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl


! Executable code
     IF (lflag) THEN
          DO i = 2, npl
               helio_plA%ahi(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          END DO
          CALL symba_helio_getacch_int(npl, nplm, helio_plA) 
     END IF
     IF (j2rp2 /= 0.0_DP) THEN
          IF (lmalloc) THEN
               ALLOCATE(xh(NDIM, nplmax), aobl(NDIM, nplmax), irh(nplmax))
               lmalloc = .FALSE.
          END IF
          DO i = 2, npl
               xh(:, i) = helio_plA%swiftest%xh(:,i)
               r2 = DOT_PRODUCT(xh(:, i), xh(:, i))
               irh(i) = 1.0_DP/SQRT(r2)
          END DO
          CALL obl_acc(npl, helio_plA%swiftest, j2rp2, j4rp4, xh, irh, aobl) 
          DO i = 2, npl
               helio_plA%ah(:,i) = helio_plA%ahi(:,i) + aobl(:, i) - aobl(:, 1)
          END DO
     ELSE
          DO i = 2, npl
               helio_plA%ah(:,i) = helio_plA%ahi(:,i)
          END DO
     END IF
     IF (lextra_force) CALL helio_user_getacch(t, npl, helio_plA) 

     RETURN

END SUBROUTINE symba_helio_getacch
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
