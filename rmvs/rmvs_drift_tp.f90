!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_drift_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Loop through test particles closely encountering a planet and call Danby drift routine
!
!  Input
!    Arguments : nenc         : number of test particles encountering current planet
!                rmvs_tpenc1P : pointer to RMVS test particle structure of first test particle encountering planet
!                mu           : mass of planet being encountered
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_tpenc1P : pointer to RMVS test particle structure of first test particle encountering planet
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL rmvs_drift_tp(nenc, rmvs_tpenc1P, mu, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine drift_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_drift_tp(nenc, rmvs_tpenc1P, mu, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_drift_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: nenc
     REAL(DP), INTENT(IN)     :: mu, dt
     TYPE(rmvs_tp), POINTER   :: rmvs_tpenc1P

! Internals
     INTEGER(I4B)              :: i, iflag
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(rmvs_tp), POINTER    :: rmvs_tpP

! Executable code
     rmvs_tpP => rmvs_tpenc1P
     DO i = 1, nenc
          swifter_tpP => rmvs_tpP%whm%swifter
          IF (swifter_tpP%status == ACTIVE) THEN
               CALL drift_one(mu, rmvs_tpP%xpc(:), rmvs_tpP%vpc(:), dt, iflag)
               IF (iflag /= 0) THEN
                    swifter_tpP%status = DISCARDED_DRIFTERR
                    WRITE(*, *) "Particle ", swifter_tpP%id, " lost due to error in Danby drift"
               END IF
          END IF
          rmvs_tpP => rmvs_tpP%tpencP
     END DO

     RETURN

END SUBROUTINE rmvs_drift_tp
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
