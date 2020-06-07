!**********************************************************************************************************************************
!
!  Unit Name   : helio_kickvb
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Kick barycentric velocities of planets
!
!  Input
!    Arguments : npl        : number of planets
!                helio_plA : pointer to head of helio planet structure linked-list
!                dt         : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_plA : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_kickvb(npl, helio_plA, dt)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_kickvb(npl, helio_plA, dt)

! Modules
     USE swiftest
     USE module_swiftest
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_kickvb
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: dt
     TYPE(helio_pl), INTENT(INOUT) :: helio_plA

! Internals
     INTEGER(I4B)              :: i

! Executable code
     DO i = 2, npl
          helio_plA%swiftest%vb(:,i) = helio_plA%swiftest%vb(:,i) + helio_plA%ah(:,i)*dt
     END DO

     RETURN

END SUBROUTINE helio_kickvb
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
