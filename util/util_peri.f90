!**********************************************************************************************************************************
!
!  Unit Name   : util_peri
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Determine system pericenter passages for test particles
!
!  Input
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active SWIFTER test particle structure linked-list
!                mu           : G * (m1 + m2) = mass of the Sun in this routine
!                msys         : total system mass
!                qmin_coord   : coordinate frame for qmin
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active SWIFTER test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_peri(lfirst, ntp, swifter_tp1P, mu, msys, qmin_coord)
!
!  Notes       : Adapted from Hal Levison's Swift routine util_peri.f
!
!                If the coordinate system used is barycentric, then this routine assumes that the barycentric coordinates in the
!                test particle structures are up-to-date and are not recomputed
!
!**********************************************************************************************************************************
SUBROUTINE util_peri(lfirst, ntp, swifter_tp1P, mu, msys, qmin_coord)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => util_peri
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)  :: lfirst
     INTEGER(I4B), INTENT(IN)  :: ntp
     REAL(DP), INTENT(IN)      :: mu, msys
     CHARACTER(*), INTENT(IN)  :: qmin_coord
     TYPE(swifter_tp), POINTER :: swifter_tp1P

! Internals
     INTEGER(I4B)              :: i, j
     REAL(DP)                  :: vdotr, e
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     swifter_tpP => swifter_tp1P
     IF (lfirst) THEN
          IF (qmin_coord == "HELIO") THEN
               DO i = 1, ntp
                    IF (swifter_tpP%status == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(swifter_tpP%xh(:), swifter_tpP%vh(:))
                         IF (vdotr > 0.0_DP) THEN
                              swifter_tpP%isperi = 1
                         ELSE
                              swifter_tpP%isperi = -1
                         END IF
                    END IF
                    swifter_tpP => swifter_tpP%nextP
               END DO
          ELSE
               DO i = 1, ntp
                    IF (swifter_tpP%status == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(swifter_tpP%xb(:), swifter_tpP%vb(:))
                         IF (vdotr > 0.0_DP) THEN
                              swifter_tpP%isperi = 1
                         ELSE
                              swifter_tpP%isperi = -1
                         END IF
                    END IF
                    swifter_tpP => swifter_tpP%nextP
               END DO
          END IF
     ELSE
          IF (qmin_coord == "HELIO") THEN
               DO i = 1, ntp
                    IF (swifter_tpP%status == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(swifter_tpP%xh(:), swifter_tpP%vh(:))
                         IF (swifter_tpP%isperi == -1) THEN
                              IF (vdotr >= 0.0_DP) THEN
                                   swifter_tpP%isperi = 0
                                   CALL orbel_xv2aeq(swifter_tpP%xh(:), swifter_tpP%vh(:), mu, swifter_tpP%atp, e,                &
                                        swifter_tpP%peri)
                              END IF
                         ELSE
                              IF (vdotr > 0.0_DP) THEN
                                   swifter_tpP%isperi = 1
                              ELSE
                                   swifter_tpP%isperi = -1
                              END IF
                         END IF
                    END IF
                    swifter_tpP => swifter_tpP%nextP
               END DO
          ELSE
               DO i = 1, ntp
                    IF (swifter_tpP%status == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(swifter_tpP%xb(:), swifter_tpP%vb(:))
                         IF (swifter_tpP%isperi == -1) THEN
                              IF (vdotr >= 0.0_DP) THEN
                                   swifter_tpP%isperi = 0
                                   CALL orbel_xv2aeq(swifter_tpP%xb(:), swifter_tpP%vb(:), msys, swifter_tpP%atp, e,              &
                                        swifter_tpP%peri)
                              END IF
                         ELSE
                              IF (vdotr > 0.0_DP) THEN
                                   swifter_tpP%isperi = 1
                              ELSE
                                   swifter_tpP%isperi = -1
                              END IF
                         END IF
                    END IF
                    swifter_tpP => swifter_tpP%nextP
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
