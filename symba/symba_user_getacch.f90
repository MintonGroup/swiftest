!**********************************************************************************************************************************
!
!  Unit Name   : symba_user_getacch
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied heliocentric accelerations to planets
!
!  Input
!    Arguments : t          : time
!                npl        : number of planets
!                symba_pl1P : pointer to head of SyMBA planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : symba_pl1P : pointer to head of SyMBA planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL symba_user_getacch(t, npl, symba_pl1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE symba_user_getacch(t, npl, symba_pl1P)

! Modules
     USE module_parameters
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_user_getacch
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: t
     TYPE(symba_pl), POINTER  :: symba_pl1P

! Internals

! Executable code

     RETURN

END SUBROUTINE symba_user_getacch
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
