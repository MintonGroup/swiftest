!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_kickvp_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Kick planetocentric velocities of active test particles in inner encounter
!
!  Input
!    Arguments : nenc         : number of test particles encountering current planet
!                rmvs_tpenc1P : pointer to RMVS test particle structure of first test particle encountering planet
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_tpenc1P : pointer to RMVS test particle structure of first test particle encountering planet
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_kickvp_tp(nenc, rmvs_tpenc1P, dt)
!
!  Notes       : Adapted from Martin Duncan's Swift routine kickvh_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_kickvp_tp(nenc, rmvs_tpenc1P, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_kickvp_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: nenc
     REAL(DP), INTENT(IN)     :: dt
     TYPE(rmvs_tp), POINTER   :: rmvs_tpenc1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(rmvs_tp), POINTER    :: rmvs_tpP

! Executable code
     rmvs_tpP => rmvs_tpenc1P
     DO i = 1, nenc
          swifter_tpP => rmvs_tpP%whm%swifter
          IF (swifter_tpP%status == ACTIVE) rmvs_tpP%vpc(:) = rmvs_tpP%vpc(:) + rmvs_tpP%apc(:)*dt
          rmvs_tpP => rmvs_tpP%tpencP
     END DO

     RETURN

END SUBROUTINE rmvs_kickvp_tp
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
