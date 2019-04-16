!**********************************************************************************************************************************
!
!  Unit Name   : whm_user_getacch_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied heliocentric accelerations to test particles
!
!  Input
!    Arguments : t        : time
!                ntp      : number of active test particles
!                whm_tp1P : pointer to head of active WHM test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : whm_tp1P : pointer to head of active WHM test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL whm_user_getacch_tp(t, ntp, whm_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE whm_user_getacch_tp(t, ntp, whm_tp1P)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_user_getacch_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: t
     TYPE(whm_tp), POINTER    :: whm_tp1P

! Internals

! Executable code

     RETURN

END SUBROUTINE whm_user_getacch_tp
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
