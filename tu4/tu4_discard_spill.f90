!**********************************************************************************************************************************
!
!  Unit Name   : tu4_discard_spill
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : tu4
!  Language    : Fortran 90/95
!
!  Description : Move spilled (discarded) TU4 test particle structure from active list to discard list
!
!  Input
!    Arguments : ntp       : number of active test particles
!                nsp       : number of spilled test particles
!                tu4_tp1P  : pointer to head of active TU4 test particle structure linked-list
!                tu4_tpd1P : pointer to head of discard TU4 test particle structure linked-list
!                tu4_tpspP : pointer to TU4 test particle structure to be discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ntp       : number of active test particles
!                nsp       : number of spilled test particles
!                tu4_tp1P  : pointer to head of active TU4 test particle structure linked-list
!                tu4_tpd1P : pointer to head of discard TU4 test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL tu4_discard_spill(ntp, nsp, tu4_tp1P, tu4_tpd1P, tu4_tpspP)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE tu4_discard_spill(ntp, nsp, tu4_tp1P, tu4_tpd1P, tu4_tpspP)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_tu4
     USE module_interfaces, EXCEPT_THIS_ONE => tu4_discard_spill
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
     TYPE(tu4_tp), POINTER       :: tu4_tp1P, tu4_tpd1P, tu4_tpspP

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpspP
     TYPE(tu4_tp), POINTER     :: tu4_tpP

! Executable code
     swifter_tpspP => tu4_tpspP%swifter
     IF (nsp == 0) THEN
          tu4_tpd1P => tu4_tpspP
     ELSE
          tu4_tpP => tu4_tpd1P
          DO i = 1, nsp - 1
               tu4_tpP => tu4_tpP%nextP
          END DO
          tu4_tpP%nextP => tu4_tpspP
          tu4_tpP%swifter%nextP => swifter_tpspP
     END IF
     IF (ASSOCIATED(tu4_tpspP%prevP)) THEN
          tu4_tpspP%prevP%nextP => tu4_tpspP%nextP
          swifter_tpspP%prevP%nextP => swifter_tpspP%nextP
     ELSE
          tu4_tp1P => tu4_tpspP%nextP
     END IF
     IF (ASSOCIATED(tu4_tpspP%nextP)) THEN
          tu4_tpspP%nextP%prevP => tu4_tpspP%prevP
          swifter_tpspP%nextP%prevP => swifter_tpspP%prevP
     END IF
     IF (nsp == 0) THEN
          NULLIFY(tu4_tpspP%prevP)
          NULLIFY(swifter_tpspP%prevP)
     ELSE
          tu4_tpspP%prevP => tu4_tpP
          swifter_tpspP%prevP => tu4_tpP%swifter
     END IF
     NULLIFY(tu4_tpspP%nextP)
     NULLIFY(swifter_tpspP%nextP)
     nsp = nsp + 1
     ntp = ntp - 1

     RETURN

END SUBROUTINE tu4_discard_spill
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
