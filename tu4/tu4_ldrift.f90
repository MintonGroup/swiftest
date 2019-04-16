!**********************************************************************************************************************************
!
!  Unit Name   : tu4_ldrift
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : tu4
!  Language    : Fortran 90/95
!
!  Description : Perform linear drift in barycentric coordinates for planets and test particles
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!                swifter_tp1P : pointer to head of active SWIFTER test particle structure linked-list
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!                swifter_tp1P : pointer to head of active SWIFTER test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL tu4_ldrift(npl, ntp, swifter_pl1P, swifter_tp1P, dt)
!
!  Notes       : Adapted from Martin Duncan's Swift routine tu4_ldrift.f
!
!**********************************************************************************************************************************
SUBROUTINE tu4_ldrift(npl, ntp, swifter_pl1P, swifter_tp1P, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => tu4_ldrift
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl, ntp
     REAL(DP), INTENT(IN)      :: dt
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     swifter_plP => swifter_pl1P
     DO i = 1, npl
          swifter_plP%xb(:) = swifter_plP%xb(:) + swifter_plP%vb(:)*dt
          swifter_plP => swifter_plP%nextP
     END DO
     swifter_tpP => swifter_tp1P
     DO i = 1, ntp
          swifter_tpP%xb(:) = swifter_tpP%xb(:) + swifter_tpP%vb(:)*dt
          swifter_tpP => swifter_tpP%nextP
     END DO

     RETURN

END SUBROUTINE tu4_ldrift
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
