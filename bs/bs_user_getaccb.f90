!**********************************************************************************************************************************
!
!  Unit Name   : bs_user_getaccb
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied barycentric accelerations to planets
!
!  Input
!    Arguments : t       : time
!                npl     : number of planets
!                bs_pl1P : pointer to head of BS planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : bs_pl1P : pointer to head of BS planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL bs_user_getaccb(t, npl, bs_pl1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE bs_user_getaccb(t, npl, bs_pl1P)

! Modules
     USE module_parameters
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_user_getaccb
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: t
     TYPE(bs_pl), POINTER     :: bs_pl1P

! Internals

! Executable code

     RETURN

END SUBROUTINE bs_user_getaccb
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
