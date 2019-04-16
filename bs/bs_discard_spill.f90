!**********************************************************************************************************************************
!
!  Unit Name   : bs_discard_spill
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Move spilled (discarded) BS test particle structure from active list to discard list
!
!  Input
!    Arguments : ntp      : number of active test particles
!                nsp      : number of spilled test particles
!                bs_tp1P  : pointer to head of active BS test particle structure linked-list
!                bs_tpd1P : pointer to head of discard BS test particle structure linked-list
!                bs_tpspP : pointer to BS test particle structure to be discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ntp      : number of active test particles
!                nsp      : number of spilled test particles
!                bs_tp1P  : pointer to head of active BS test particle structure linked-list
!                bs_tpd1P : pointer to head of discard BS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL bs_discard_spill(ntp, nsp, bs_tp1P, bs_tpd1P, bs_tpspP)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE bs_discard_spill(ntp, nsp, bs_tp1P, bs_tpd1P, bs_tpspP)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_discard_spill
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
     TYPE(bs_tp), POINTER        :: bs_tp1P, bs_tpd1P, bs_tpspP

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpspP
     TYPE(bs_tp), POINTER      :: bs_tpP

! Executable code
     swifter_tpspP => bs_tpspP%swifter
     IF (nsp == 0) THEN
          bs_tpd1P => bs_tpspP
     ELSE
          bs_tpP => bs_tpd1P
          DO i = 1, nsp - 1
               bs_tpP => bs_tpP%nextP
          END DO
          bs_tpP%nextP => bs_tpspP
          bs_tpP%swifter%nextP => swifter_tpspP
     END IF
     IF (ASSOCIATED(bs_tpspP%prevP)) THEN
          bs_tpspP%prevP%nextP => bs_tpspP%nextP
          swifter_tpspP%prevP%nextP => swifter_tpspP%nextP
     ELSE
          bs_tp1P => bs_tpspP%nextP
     END IF
     IF (ASSOCIATED(bs_tpspP%nextP)) THEN
          bs_tpspP%nextP%prevP => bs_tpspP%prevP
          swifter_tpspP%nextP%prevP => swifter_tpspP%prevP
     END IF
     IF (nsp == 0) THEN
          NULLIFY(bs_tpspP%prevP)
          NULLIFY(swifter_tpspP%prevP)
     ELSE
          bs_tpspP%prevP => bs_tpP
          swifter_tpspP%prevP => bs_tpP%swifter
     END IF
     NULLIFY(bs_tpspP%nextP)
     NULLIFY(swifter_tpspP%nextP)
     nsp = nsp + 1
     ntp = ntp - 1

     RETURN

END SUBROUTINE bs_discard_spill
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
