!**********************************************************************************************************************************
!
!  Unit Name   : helio_user_getacch
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied heliocentric accelerations to planets
!
!  Input
!    Arguments : t          : time
!                npl        : number of planets
!                helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL helio_user_getacch(t, npl, helio_pl1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE helio_user_getacch(t, npl, helio_pl1P)

! Modules
     USE module_parameters
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_user_getacch
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: t
     TYPE(helio_pl), POINTER  :: helio_pl1P

! Internals

! Executable code

     RETURN

END SUBROUTINE helio_user_getacch
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
