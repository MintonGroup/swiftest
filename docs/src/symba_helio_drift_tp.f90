!**********************************************************************************************************************************
!
!  Unit Name   : symba_helio_drift_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Loop through test particles and call Danby drift routine
!
!  Input
!    Arguments : irec       : input recursion level
!                ntp        : number of active test particles
!                symba_tp1P : pointer to head of active SyMBA test particle structure linked-list
!                mu         : mass of the Sun
!                dt         : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_tp1P : pointer to head of active SyMBA test particle structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL symba_helio_drift_tp(irec, ntp, symba_tp1P, mu, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_helio_drift.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_helio_drift_tp(irec, ntp, symba_tpA, mu, dt)

! Modules
     USE swiftest
     USE module_swiftest
     USE helio
     USE symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_helio_drift_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)      :: irec, ntp
     REAL(DP), INTENT(IN)          :: mu, dt
     TYPE(symba_tp), INTENT(INOUT) :: symba_tpA

! Internals
     INTEGER(I4B)              :: i, iflag

! Executable code
     DO i = 1, ntp
          IF ((symba_tpA%levelg(i) == irec) .AND. (symba_tpA%helio%swiftest%status(i) == ACTIVE)) THEN
               CALL drift_one(mu, symba_tpA%helio%swiftest%xh(:,i), symba_tpA%helio%swiftest%vb(:,i), dt, iflag)
               IF (iflag /= 0) THEN
                    symba_tpA%helio%swiftest%status(i) = DISCARDED_DRIFTERR
                    WRITE(*, *) "Particle ", symba_tpA%helio%swiftest%name(i), " lost due to error in Danby drift"
               END IF
          END IF
     END DO

     RETURN

END SUBROUTINE symba_helio_drift_tp
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
