!**********************************************************************************************************************************
!
!  Unit Name   : helio_drift_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Loop through test particles and call Danby drift routine
!
!  Input
!    Arguments : ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL helio_drift_tp(ntp, swifter_tp1P, mu, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine drift_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_drift_tp(ntp, swiftest_tpA, mu, dt)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => helio_drift_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)         :: ntp
     REAL(DP), INTENT(IN)             :: mu, dt
     TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA

! Internals
     INTEGER(I4B)                     :: i, iflag

! Executable code
     DO i = 1, ntp
          IF (swiftest_tpA%status(i) == ACTIVE) THEN
               CALL drift_one(mu, swiftest_tpA%xh(:,i), swiftest_tpA%vb(:,i), dt, iflag)
               IF (iflag /= 0) THEN
                    swiftest_tpA%status(i) = DISCARDED_DRIFTERR
                    WRITE(*, *) "Particle ", swiftest_tpA%id(i), " lost due to error in Danby drift"
               END IF
          END IF
     END DO

     RETURN

END SUBROUTINE helio_drift_tp
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
