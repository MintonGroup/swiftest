!**********************************************************************************************************************************
!
!  Unit Name   : ra15_user_getaccb_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : ra15
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied barycentric accelerations to test particles
!
!  Input
!    Arguments : t         : time
!                ntp       : number of active test particles
!                ra15_tp1P : pointer to head of active RA15 test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : ra15_tp1P : pointer to head of active RA15 test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL ra15_user_getaccb_tp(t, ntp, ra15_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE ra15_user_getaccb_tp(t, ntp, ra15_tp1P)

! Modules
     USE module_parameters
     USE module_ra15
     USE module_interfaces, EXCEPT_THIS_ONE => ra15_user_getaccb_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: t
     TYPE(ra15_tp), POINTER   :: ra15_tp1P

! Internals

! Executable code

     RETURN

END SUBROUTINE ra15_user_getaccb_tp
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
