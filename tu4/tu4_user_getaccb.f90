!**********************************************************************************************************************************
!
!  Unit Name   : tu4_user_getaccb
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : tu4
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied barycentric accelerations to planets
!
!  Input
!    Arguments : t        : time
!                npl      : number of planets
!                tu4_pl1P : pointer to head of TU4 planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : tu4_pl1P : pointer to head of TU4 planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL tu4_user_getaccb(t, npl, tu4_pl1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE tu4_user_getaccb(t, npl, tu4_pl1P)

! Modules
     USE module_parameters
     USE module_tu4
     USE module_interfaces, EXCEPT_THIS_ONE => tu4_user_getaccb
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: t
     TYPE(tu4_pl), POINTER    :: tu4_pl1P

! Internals

! Executable code

     RETURN

END SUBROUTINE tu4_user_getaccb
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
