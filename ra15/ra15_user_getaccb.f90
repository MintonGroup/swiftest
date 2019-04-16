!**********************************************************************************************************************************
!
!  Unit Name   : ra15_user_getaccb
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : ra15
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied barycentric accelerations to planets
!
!  Input
!    Arguments : t         : time
!                npl       : number of planets
!                ra15_pl1P : pointer to head of RA15 planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : ra15_pl1P : pointer to head of RA15 planet structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL ra15_user_getaccb(t, npl, ra15_pl1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE ra15_user_getaccb(t, npl, ra15_pl1P)

! Modules
     USE module_parameters
     USE module_ra15
     USE module_interfaces, EXCEPT_THIS_ONE => ra15_user_getaccb
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: t
     TYPE(ra15_pl), POINTER   :: ra15_pl1P

! Internals

! Executable code

     RETURN

END SUBROUTINE ra15_user_getaccb
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
