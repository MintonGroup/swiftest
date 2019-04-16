!**********************************************************************************************************************************
!
!  Unit Name   : whm_discard_spill
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Move spilled (discarded) WHM test particle structure from active list to discard list
!
!  Input
!    Arguments : ntp       : number of active test particles
!                nsp       : number of spilled test particles
!                whm_tp1P  : pointer to head of active WHM test particle structure linked-list
!                whm_tpd1P : pointer to head of discard WHM test particle structure linked-list
!                whm_tpspP : pointer to WHM test particle structure to be discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ntp       : number of active test particles
!                nsp       : number of spilled test particles
!                whm_tp1P  : pointer to head of active WHM test particle structure linked-list
!                whm_tpd1P : pointer to head of discard WHM test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_discard_spill(ntp, nsp, whm_tp1P, whm_tpd1P, whm_tpspP)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE whm_discard_spill(ntp, nsp, whm_tp1P, whm_tpd1P, whm_tpspP)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_discard_spill
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
     TYPE(whm_tp), POINTER       :: whm_tp1P, whm_tpd1P, whm_tpspP

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpspP
     TYPE(whm_tp), POINTER     :: whm_tpP

! Executable code
     swifter_tpspP => whm_tpspP%swifter
     IF (nsp == 0) THEN
          whm_tpd1P => whm_tpspP
     ELSE
          whm_tpP => whm_tpd1P
          DO i = 1, nsp - 1
               whm_tpP => whm_tpP%nextP
          END DO
          whm_tpP%nextP => whm_tpspP
          whm_tpP%swifter%nextP => swifter_tpspP
     END IF
     IF (ASSOCIATED(whm_tpspP%prevP)) THEN
          whm_tpspP%prevP%nextP => whm_tpspP%nextP
          swifter_tpspP%prevP%nextP => swifter_tpspP%nextP
     ELSE
          whm_tp1P => whm_tpspP%nextP
     END IF
     IF (ASSOCIATED(whm_tpspP%nextP)) THEN
          whm_tpspP%nextP%prevP => whm_tpspP%prevP
          swifter_tpspP%nextP%prevP => swifter_tpspP%prevP
     END IF
     IF (nsp == 0) THEN
          NULLIFY(whm_tpspP%prevP)
          NULLIFY(swifter_tpspP%prevP)
     ELSE
          whm_tpspP%prevP => whm_tpP
          swifter_tpspP%prevP => whm_tpP%swifter
     END IF
     NULLIFY(whm_tpspP%nextP)
     NULLIFY(swifter_tpspP%nextP)
     nsp = nsp + 1
     ntp = ntp - 1

     RETURN

END SUBROUTINE whm_discard_spill
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
