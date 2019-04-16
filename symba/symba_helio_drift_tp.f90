!**********************************************************************************************************************************
!
!  Unit Name   : symba_helio_drift_tp
!  Unit Type   : subroutine
!  Project     : Swifter
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
SUBROUTINE symba_helio_drift_tp(irec, ntp, symba_tp1P, mu, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_helio_drift_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: irec, ntp
     REAL(DP), INTENT(IN)     :: mu, dt
     TYPE(symba_tp), POINTER  :: symba_tp1P

! Internals
     INTEGER(I4B)              :: i, iflag
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(symba_tp), POINTER   :: symba_tpP

! Executable code
     !Removed by D. Minton
     !symba_tpP => symba_tp1P
     !^^^^^^^^^^^^^^^^^^^^^
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) & 
     !$OMP PRIVATE(i,swifter_tpP,symba_tpP,iflag) &
     !$OMP SHARED(ntp,symba_tp1P,mu,dt,irec) 
     DO i = 1, ntp
          symba_tpP => symba_tp1P%symba_tpPA(i)%thisP
          swifter_tpP => symba_tpP%helio%swifter
          IF ((symba_tpP%levelg == irec) .AND. (swifter_tpP%status == ACTIVE)) THEN
               CALL drift_one(mu, swifter_tpP%xh(:), swifter_tpP%vb(:), dt, iflag)
               IF (iflag /= 0) THEN
                    swifter_tpP%status = DISCARDED_DRIFTERR
                    WRITE(*, *) "Particle ", swifter_tpP%id, " lost due to error in Danby drift"
               END IF
          END IF
          !Removed by D. Minton
          !symba_tpP => symba_tpP%nextP
          !^^^^^^^^^^^^^^^^^^^^
     END DO
     !$OMP END PARALLEL DO

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
