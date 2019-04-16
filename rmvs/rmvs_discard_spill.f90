!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_discard_spill
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Move spilled (discarded) RMVS test particle structure from active list to discard list
!
!  Input
!    Arguments : ntp        : number of active test particles
!                nsp        : number of spilled test particles
!                rmvs_tp1P  : pointer to head of active RMVS test particle structure linked-list
!                rmvs_tpd1P : pointer to head of discard RMVS test particle structure linked-list
!                rmvs_tpspP : pointer to RMVS test particle structure to be discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ntp        : number of active test particles
!                nsp        : number of spilled test particles
!                rmvs_tp1P  : pointer to head of active RMVS test particle structure linked-list
!                rmvs_tpd1P : pointer to head of discard RMVS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_discard_spill(ntp, nsp, rmvs_tp1P, rmvs_tpd1P, rmvs_tpspP)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_discard_spill(ntp, nsp, rmvs_tp1P, rmvs_tpd1P, rmvs_tpspP)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_discard_spill
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
     TYPE(rmvs_tp), POINTER      :: rmvs_tp1P, rmvs_tpd1P, rmvs_tpspP

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpspP
     TYPE(whm_tp), POINTER     :: whm_tpspP
     TYPE(rmvs_tp), POINTER    :: rmvs_tpP

! Executable code
     whm_tpspP => rmvs_tpspP%whm
     swifter_tpspP => whm_tpspP%swifter
     IF (nsp == 0) THEN
          rmvs_tpd1P => rmvs_tpspP
     ELSE
          rmvs_tpP => rmvs_tpd1P
          DO i = 1, nsp - 1
               rmvs_tpP => rmvs_tpP%nextP
          END DO
          rmvs_tpP%nextP => rmvs_tpspP
          rmvs_tpP%whm%nextP => whm_tpspP
          rmvs_tpP%whm%swifter%nextP => swifter_tpspP
     END IF
     IF (ASSOCIATED(rmvs_tpspP%prevP)) THEN
          rmvs_tpspP%prevP%nextP => rmvs_tpspP%nextP
          whm_tpspP%prevP%nextP => whm_tpspP%nextP
          swifter_tpspP%prevP%nextP => swifter_tpspP%nextP
     ELSE
          rmvs_tp1P => rmvs_tpspP%nextP
     END IF
     IF (ASSOCIATED(rmvs_tpspP%nextP)) THEN
          rmvs_tpspP%nextP%prevP => rmvs_tpspP%prevP
          whm_tpspP%nextP%prevP => whm_tpspP%prevP
          swifter_tpspP%nextP%prevP => swifter_tpspP%prevP
     END IF
     IF (nsp == 0) THEN
          NULLIFY(rmvs_tpspP%prevP)
          NULLIFY(whm_tpspP%prevP)
          NULLIFY(swifter_tpspP%prevP)
     ELSE
          rmvs_tpspP%prevP => rmvs_tpP
          whm_tpspP%prevP => rmvs_tpP%whm
          swifter_tpspP%prevP => rmvs_tpP%whm%swifter
     END IF
     NULLIFY(rmvs_tpspP%nextP)
     NULLIFY(whm_tpspP%nextP)
     NULLIFY(swifter_tpspP%nextP)
     nsp = nsp + 1
     ntp = ntp - 1

     RETURN

END SUBROUTINE rmvs_discard_spill
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
