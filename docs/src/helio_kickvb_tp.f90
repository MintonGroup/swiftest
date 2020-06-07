!**********************************************************************************************************************************
!
!  Unit Name   : helio_kickvb_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
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
SUBROUTINE helio_kickvb_tp(ntp, helio_tpA, dt)

! Modules
     USE swiftest
     USE module_swiftest
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_kickvb_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: dt
     TYPE(helio_tp), INTENT(INOUT)  :: helio_tpA

! Internals
     INTEGER(I4B)              :: i

! Executable code
     DO i = 1, ntp
          IF (helio_tpA%swiftest%status(i) == ACTIVE) helio_tpA%swiftest%vb(:,i) = helio_tpA%swiftest%vb(:,i) + helio_tpA%ah(:,i)*dt
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
