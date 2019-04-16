!**********************************************************************************************************************************
!
!  Unit Name   : whm_user_getacch
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied heliocentric accelerations to planets
!
!  Input
!    Arguments : t        : time
!                npl      : number of planets
!                whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL whm_user_getacch(t, npl, whm_pl1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE whm_user_getacch(t, npl, whm_pl1P)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_user_getacch
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: t
     TYPE(whm_pl), POINTER    :: whm_pl1P

! Internals

! Executable code

     RETURN

END SUBROUTINE whm_user_getacch
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
