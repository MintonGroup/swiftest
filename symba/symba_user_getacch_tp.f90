!**********************************************************************************************************************************
!
!  Unit Name   : symba_user_getacch_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied heliocentric accelerations to test particles
!
!  Input
!    Arguments : t          : time
!                ntp        : number of active test particles
!                symba_tp1P : pointer to head of active SyMBA test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : symba_tp1P : pointer to head of active SyMBA test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL symba_user_getacch_tp(t, ntp, symba_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE symba_user_getacch_tp(t, ntp, symba_tp1P)

! Modules
     USE module_parameters
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_user_getacch_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: t
     TYPE(symba_tp), POINTER  :: symba_tp1P

! Internals

! Executable code

     RETURN

END SUBROUTINE symba_user_getacch_tp
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
