!**********************************************************************************************************************************
!
!  Unit Name   : discard_peri
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : discard
!  Language    : Fortran 90/95
!
!  Description : Check to see if a test particle should be discarded because its perihelion distance becomes too small
!
!  Input
!    Arguments : t              : time
!                npl            : number of planets
!                ntp            : number of active test particles
!                swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P   : pointer to head of active Swifter test particle structure linked-list
!                msys           : total system mass
!                qmin           : minimum pericenter distance for test particle
!                qmin_alo       : minimum semimajor axis for qmin
!                qmin_ahi       : maximum semimajor axis for qmin
!                qmin_coord     : coordinate frame to use for qmin
!                lrhill_present : logical flag indicating whether Hill sphere radii for planets are present
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P   : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : status message
!    File      : none
!
!  Invocation  : CALL discard_peri(t, npl, ntp, swifter_pl1P, swifter_tp1P, msys, qmin, qmin_alo, qmin_ahi, qmin_coord,
!                                  lrhill_present)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_peri.f
!
!**********************************************************************************************************************************
SUBROUTINE discard_peri(t, npl, ntp, swiftest_plA, swiftest_tpA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, lrhill_present)

! Modules
     USE swiftest, EXCEPT_THIS_ONE => discard_peri
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)  :: lrhill_present
     INTEGER(I4B), INTENT(IN)  :: npl, ntp
     REAL(DP), INTENT(IN)      :: t, msys, qmin, qmin_alo, qmin_ahi
     CHARACTER(*), INTENT(IN)  :: qmin_coord
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
     TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA

! Internals
     LOGICAL(LGT), SAVE        :: lfirst = .TRUE.
     INTEGER(I4B)              :: i, j, ih
     REAL(DP)                  :: r2
     REAL(DP), DIMENSION(NDIM) :: dx


! Executable code
     IF (lfirst) THEN
          IF (.NOT. lrhill_present) CALL util_hills(npl, swiftest_plA)
          CALL util_peri(lfirst, ntp, swiftest_tpA, swiftest_plA%mass(1), msys, qmin_coord)
          lfirst = .FALSE.
     ELSE
          CALL util_peri(lfirst, ntp, swiftest_tpA, swiftest_plA%mass(1), msys, qmin_coord)
          DO i = 1, ntp
               IF (swiftest_tpA%status(i) == ACTIVE) THEN
                    IF (swiftest_tpA%isperi(i) == 0) THEN
                         ih = 1
                         DO j = 2, npl
                              dx(:) = swiftest_tpA%xh(:,i) - swiftest_plA%xh(:,j)
                              r2 = DOT_PRODUCT(dx(:), dx(:))
                              IF (r2 <= swiftest_plA%rhill(j)*swiftest_plA%rhill(j)) ih = 0
                         END DO
                         IF (ih == 1) THEN
                              IF ((swiftest_tpA%atp(i) >= qmin_alo) .AND.      &
                                   (swiftest_tpA%atp(i) <= qmin_ahi) .AND.      &                 
                                  (swiftest_tpA%peri(i) <= qmin)) THEN
                                   swiftest_tpA%status(i) = DISCARDED_PERI
                                   WRITE(*, *) "Particle ", swiftest_tpA%name(i), " perihelion distance too small at t = ", t
                                   ldiscard_tp = .TRUE.
                              END IF
                         END IF
                    END IF
               END IF
          END DO
     END IF

     RETURN

END SUBROUTINE discard_peri
!**********************************************************************************************************************************
!
!  Author(s)   : Davname E. Kaufmann
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
