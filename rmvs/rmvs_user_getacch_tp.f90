!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_user_getacch_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Add user-supplied heliocentric accelerations to test particles encountering a planet
!
!  Input
!    Arguments : t            : time
!                nenc         : number of test particles encountering planet
!                rmvs_tpenc1P : pointer to head of encountering RMVS test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : rmvs_tpenc1P : pointer to head of encountering RMVS test particle structure linked-list
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Invocation  : CALL rmvs_user_getacch_tp(t, nenc, rmvs_tpenc1P)
!
!  Notes       : In this routine only loop over test particles interacting with the current planet
!
!                Proceed from the first test particle to the next on the interaction list by setting rmvs_tpP=>rmvs_tpenc1P to
!                begin (outside the loop), then rmvs_tpP=>rmvs_tpP%tpencP inside the loop
!
!                Use the current heliocentric position rmvs_tpP%whm%swifter%xh(:) to compute the acceleration components to add
!
!                Add the accelerations to rmvs_tpP%apc(:), the planetocentric acceleration
!
!                To add extra accelerations in RMVS to non-interacting test particles use whm_user_getacch_tp()
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_user_getacch_tp(t, nenc, rmvs_tpenc1P)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_user_getacch_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: nenc
     REAL(DP), INTENT(IN)     :: t
     TYPE(rmvs_tp), POINTER   :: rmvs_tpenc1P

! Internals

! Executable code

     RETURN

END SUBROUTINE rmvs_user_getacch_tp
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
