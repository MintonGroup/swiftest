!**********************************************************************************************************************************
!
!  Unit Name   : helio_kickvb_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Kick barycentric velocities of active test particles
!
!  Input
!    Arguments : ntp        : number of active test particles
!                helio_tp1P : pointer to head of active helio test particle structure linked-list
!                dt         : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_tp1P : pointer to head of active helio test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_kickvb_tp(ntp, helio_tp1P, dt)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine kickvh_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_kickvb_tp(ntp, helio_tp1P, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_kickvb_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: dt
     TYPE(helio_tp), POINTER  :: helio_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(helio_tp), POINTER   :: helio_tpP

! Executable code
     helio_tpP => helio_tp1P
     DO i = 1, ntp
          swifter_tpP => helio_tpP%swifter
          IF (swifter_tpP%status == ACTIVE) swifter_tpP%vb(:) = swifter_tpP%vb(:) + helio_tpP%ah(:)*dt
          helio_tpP => helio_tpP%nextP
     END DO

     RETURN

END SUBROUTINE helio_kickvb_tp
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
