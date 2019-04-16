!**********************************************************************************************************************************
!
!  Unit Name   : tu4_kickvb
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : tu4
!  Language    : Fortran 90/95
!
!  Description : Kick barycentric velocities of planets and test particles
!
!  Input
!    Arguments : npl      : number of planets
!                ntp      : number of active test particles
!                tu4_pl1P : pointer to head of TU4 planet structure linked-list
!                tu4_tp1P : pointer to head of active TU4 test particle structure linked-list
!                dt       : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : tu4_pl1P : pointer to head of TU4 planet structure linked-list
!                tu4_tp1P : pointer to head of active TU4 test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL tu4_kickvb(npl, ntp, tu4_pl1P, tu4_tp1P, dt)
!
!  Notes       : Adapted from Martin Duncan's Swift routine tu4_vkickb.f
!
!**********************************************************************************************************************************
SUBROUTINE tu4_kickvb(npl, ntp, tu4_pl1P, tu4_tp1P, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_tu4
     USE module_interfaces, EXCEPT_THIS_ONE => tu4_kickvb
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl, ntp
     REAL(DP), INTENT(IN)     :: dt
     TYPE(tu4_pl), POINTER    :: tu4_pl1P
     TYPE(tu4_tp), POINTER    :: tu4_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(tu4_pl), POINTER     :: tu4_plP
     TYPE(tu4_tp), POINTER     :: tu4_tpP

! Executable code
     tu4_plP => tu4_pl1P
     DO i = 1, npl
          swifter_plP => tu4_plP%swifter
          swifter_plP%vb(:) = swifter_plP%vb(:) + tu4_plP%ab(:)*dt
          tu4_plP => tu4_plP%nextP
     END DO
     tu4_tpP => tu4_tp1P
     DO i = 1, ntp
          swifter_tpP => tu4_tpP%swifter
          swifter_tpP%vb(:) = swifter_tpP%vb(:) + tu4_tpP%ab(:)*dt
          tu4_tpP => tu4_tpP%nextP
     END DO

     RETURN

END SUBROUTINE tu4_kickvb
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
