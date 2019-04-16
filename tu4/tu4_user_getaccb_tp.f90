!**********************************************************************************************************************************
!
!  Unit Name   : tu4_user_getaccb_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : tu4
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied barycentric accelerations to test particles
!
!  Input
!    Arguments : t        : time
!                ntp      : number of active test particles
!                tu4_tp1P : pointer to head of active TU4 test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : tu4_tp1P : pointer to head of active TU4 test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL tu4_user_getaccb_tp(t, ntp, tu4_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE tu4_user_getaccb_tp(t, ntp, tu4_tp1P)

! Modules
     USE module_parameters
     USE module_tu4
     USE module_interfaces, EXCEPT_THIS_ONE => tu4_user_getaccb_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: t
     TYPE(tu4_tp), POINTER    :: tu4_tp1P

! Internals

! Executable code

     RETURN

END SUBROUTINE tu4_user_getaccb_tp
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
