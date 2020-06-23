!**********************************************************************************************************************************
!
!  Unit Name   : whm_kickvh_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Kick heliocentric velocities of active test particles
!
!  Input
!    Arguments : ntp      : number of active test particles
!                whm_tp1P : pointer to head of active WHM test particle structure linked-list
!                dt       : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_tp1P : pointer to head of active WHM test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_kickvh_tp(ntp, whm_tp1P, dt)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine kickvh_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_kickvh_tp(ntp, whm_tp1P, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_kickvh_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: dt
     TYPE(whm_tp), POINTER    :: whm_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(whm_tp), POINTER     :: whm_tpP
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     whm_tpP => whm_tp1P
     DO i = 1, ntp
          swifter_tpP => whm_tpP%swifter
          IF (swifter_tpP%status == ACTIVE) swifter_tpP%vh(:) = swifter_tpP%vh(:) + whm_tpP%ah(:)*dt
          whm_tpP => whm_tpP%nextP
     END DO

     RETURN

END SUBROUTINE whm_kickvh_tp
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
