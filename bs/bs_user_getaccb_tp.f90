!**********************************************************************************************************************************
!
!  Unit Name   : bs_user_getaccb_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied barycentric accelerations to test particles
!
!  Input
!    Arguments : t       : time
!                ntp     : number of active test particles
!                bs_tp1P : pointer to head of active BS test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : bs_tp1P : pointer to head of active BS test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL bs_user_getaccb_tp(t, ntp, bs_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE bs_user_getaccb_tp(t, ntp, bs_tp1P)

! Modules
     USE module_parameters
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_user_getaccb_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: t
     TYPE(bs_tp), POINTER     :: bs_tp1P

! Internals

! Executable code

     RETURN

END SUBROUTINE bs_user_getaccb_tp
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
