!**********************************************************************************************************************************
!
!  Unit Name   : util_peri
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Determine system pericenter passages for test particles
!
!  Input
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                ntp          : number of active test particles
!                swiftest_tp1P : pointer to head of active SWIFTER test particle structure linked-list
!                mu           : G * (m1 + m2) = mass of the Sun in this routine
!                msys         : total system mass
!                qmin_coord   : coordinate frame for qmin
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swiftest_tp1P : pointer to head of active SWIFTER test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_peri(lfirst, ntp, swiftest_tp1P, mu, msys, qmin_coord)
!
!  Notes       : Adapted from Hal Levison's Swift routine util_peri.f
!
!                If the coordinate system used is barycentric, then this routine assumes that the barycentric coordinates in the
!                test particle structures are up-to-date and are not recomputed
!
!**********************************************************************************************************************************
SUBROUTINE util_peri(lfirst, ntp, swiftest_tpA, mu, msys, qmin_coord)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => util_peri
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)  :: lfirst
     INTEGER(I4B), INTENT(IN)  :: ntp
     REAL(DP), INTENT(IN)      :: mu, msys
     CHARACTER(*), INTENT(IN)  :: qmin_coord
     TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: vdotr, e

! Executable code
     IF (lfirst) THEN
          IF (qmin_coord == "HELIO") THEN
               DO i = 1, ntp
                    IF (swiftest_tpA%status(i) == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i))
                         IF (vdotr > 0.0_DP) THEN
                              swiftest_tpA%isperi(i) = 1
                         ELSE
                              swiftest_tpA%isperi(i) = -1
                         END IF
                    END IF
               END DO
          ELSE
               DO i = 1, ntp
                    IF (swiftest_tpA%status(i) == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(swiftest_tpA%xb(:,i), swiftest_tpA%vb(:,i))
                         IF (vdotr > 0.0_DP) THEN
                              swiftest_tpA%isperi(i) = 1
                         ELSE
                              swiftest_tpA%isperi(i) = -1
                         END IF
                    END IF
               END DO
          END IF
     ELSE
          IF (qmin_coord == "HELIO") THEN
               DO i = 1, ntp
                    IF (swiftest_tpA%status(i) == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i))
                         IF (swiftest_tpA%isperi(i) == -1) THEN
                              IF (vdotr >= 0.0_DP) THEN
                                   swiftest_tpA%isperi(i) = 0
                                   CALL orbel_xv2aeq(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i), mu, & 
                                        swiftest_tpA%atp(i), e, swiftest_tpA%peri(i))
                              END IF
                         ELSE
                              IF (vdotr > 0.0_DP) THEN
                                   swiftest_tpA%isperi(i) = 1
                              ELSE
                                   swiftest_tpA%isperi(i) = -1
                              END IF
                         END IF
                    END IF
               END DO
          ELSE
               DO i = 1, ntp
                    IF (swiftest_tpA%status(i) == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(swiftest_tpA%xb(:,i), swiftest_tpA%vb(:,i))
                         IF (swiftest_tpA%isperi(i) == -1) THEN
                              IF (vdotr >= 0.0_DP) THEN
                                   swiftest_tpA%isperi(i) = 0
                                   CALL orbel_xv2aeq(swiftest_tpA%xb(:,i), swiftest_tpA%vb(:,i), msys, & 
                                    swiftest_tpA%atp(i), e, swiftest_tpA%peri(i))
                              END IF
                         ELSE
                              IF (vdotr > 0.0_DP) THEN
                                   swiftest_tpA%isperi(i) = 1
                              ELSE
                                   swiftest_tpA%isperi(i) = -1
                              END IF
                         END IF
                    END IF
               END DO
          END IF
     END IF

     RETURN

END SUBROUTINE util_peri
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
