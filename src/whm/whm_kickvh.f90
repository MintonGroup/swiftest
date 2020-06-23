!**********************************************************************************************************************************
!
!  Unit Name   : whm_kickvh
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Kick heliocentric velocities of planets
!
!  Input
!    Arguments : npl      : number of planets
!                whm_pl1P : pointer to head of WHM planet structure linked-list
!                dt       : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_kickvh(npl, whm_pl1P, dt)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_kickvh(npl, whm_pl1P, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_kickvh
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: dt
     TYPE(whm_pl), POINTER    :: whm_pl1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(whm_pl), POINTER     :: whm_plP

! Executable code
     whm_plP => whm_pl1P
     DO i = 2, npl
          whm_plP => whm_plP%nextP
          swifter_plP => whm_plP%swifter
          swifter_plP%vh(:) = swifter_plP%vh(:) + whm_plP%ah(:)*dt
     END DO

     RETURN

END SUBROUTINE whm_kickvh
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
