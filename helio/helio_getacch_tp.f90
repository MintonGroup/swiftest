!**********************************************************************************************************************************
!
!  Unit Name   : helio_getacch_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Compute heliocentric accelerations of test particles
!
!  Input
!    Arguments : lflag        : logical flag indicating whether to recompute direct cross term accelerations
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                helio_pl1P   : pointer to head of helio planet structure linked-list
!                helio_tp1P   : pointer to head of active helio test particle structure linked-list
!                xh           : heliocentric positions of planets at time t
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_tp1P   : pointer to head of active helio test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_getacch_tp(lflag, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, xh, j2rp2, j4rp4)
!
!  Notes       : Adapted from Hal Levison's Swift routine helio_getacch_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_getacch_tp(lflag, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, xh, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_getacch_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                   :: lflag, lextra_force
     INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
     TYPE(helio_pl), INTENT(INOUT) :: helio_plA
     TYPE(helio_tp), INTENT(INOUT) :: helio_tpA

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP)                                     :: r2, fac, mu
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh, irht
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: aobl, xht, aoblt

! Executable code
     IF (lflag) THEN
          DO i = 1, ntp
               helio_tpA%ahi(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          END DO
          CALL helio_getacch_int_tp(npl, ntp, swifter_plA, helio_tpA, xh)
     END IF
     IF (j2rp2 /= 0.0_DP) THEN
          IF (lmalloc) THEN
               ALLOCATE(aobl(NDIM, nplmax), irh(nplmax), xht(NDIM, ntpmax), aoblt(NDIM, ntpmax), irht(ntpmax))
               lmalloc = .FALSE.
          END IF
          DO i = 2, npl
               r2 = DOT_PRODUCT(xh(:, i), xh(:, i))
               irh(i) = 1.0_DP/SQRT(r2)
          END DO
          CALL obl_acc(npl, swifter_plA, j2rp2, j4rp4, xh, irh, aobl)
          mu = helio_plA%swiftest%mass(1)
          DO i = 1, ntp
               xht(:, i) = helio_tpA%swiftest%xh(:,i)
               r2 = DOT_PRODUCT(xht(:, i), xht(:, i))
               irht(i) = 1.0_DP/SQRT(r2)
          END DO
          CALL obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
          DO i = 1, ntp
               helio_tpA%ah(:,i) = helio_tpA%ahi(:,i) + aoblt(:, i) - aobl(:, 1)
          END DO
     ELSE
          DO i = 1, ntp
               helio_tpA%ah(:,i) = helio_tpA%ahi(:,i)
          END DO
     END IF
     IF (lextra_force) CALL helio_user_getacch_tp(t, ntp, helio_tpA)

     RETURN

END SUBROUTINE helio_getacch_tp
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
