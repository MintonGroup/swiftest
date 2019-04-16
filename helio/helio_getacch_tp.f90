!**********************************************************************************************************************************
!
!  Unit Name   : helio_getacch_tp
!  Unit Type   : subroutine
!  Project     : Swifter
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
SUBROUTINE helio_getacch_tp(lflag, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, xh, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_getacch_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                   :: lflag, lextra_force
     INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
     TYPE(helio_pl), POINTER                    :: helio_pl1P
     TYPE(helio_tp), POINTER                    :: helio_tp1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP)                                     :: r2, fac, mu
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh, irht
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: aobl, xht, aoblt
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P
     TYPE(swifter_tp), POINTER                    :: swifter_tpP
     TYPE(helio_tp), POINTER                      :: helio_tpP

! Executable code
     swifter_pl1P => helio_pl1P%swifter
     IF (lflag) THEN
          helio_tpP => helio_tp1P
          DO i = 1, ntp
               helio_tpP%ahi(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
               helio_tpP => helio_tpP%nextP
          END DO
          CALL helio_getacch_int_tp(npl, ntp, swifter_pl1P, helio_tp1P, xh)
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
          CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
          mu = swifter_pl1P%mass
          swifter_tpP => helio_tp1P%swifter
          DO i = 1, ntp
               xht(:, i) = swifter_tpP%xh(:)
               r2 = DOT_PRODUCT(xht(:, i), xht(:, i))
               irht(i) = 1.0_DP/SQRT(r2)
               swifter_tpP => swifter_tpP%nextP
          END DO
          CALL obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
          helio_tpP => helio_tp1P
          DO i = 1, ntp
               helio_tpP%ah(:) = helio_tpP%ahi(:) + aoblt(:, i) - aobl(:, 1)
               helio_tpP => helio_tpP%nextP
          END DO
     ELSE
          helio_tpP => helio_tp1P
          DO i = 1, ntp
               helio_tpP%ah(:) = helio_tpP%ahi(:)
               helio_tpP => helio_tpP%nextP
          END DO
     END IF
     IF (lextra_force) CALL helio_user_getacch_tp(t, ntp, helio_tp1P)

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
