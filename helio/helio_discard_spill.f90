!**********************************************************************************************************************************
!
!  Unit Name   : helio_discard_spill
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Move spilled (discarded) helio test particle structure from active list to discard list
!
!  Input
!    Arguments : ntp         : number of active test particles
!                nsp         : number of spilled test particles
!                helio_tp1P  : pointer to head of active helio test particle structure linked-list
!                helio_tpd1P : pointer to head of discard helio test particle structure linked-list
!                helio_tpspP : pointer to helio test particle structure to be discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ntp         : number of active test particles
!                nsp         : number of spilled test particles
!                helio_tp1P  : pointer to head of active helio test particle structure linked-list
!                helio_tpd1P : pointer to head of discard helio test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_discard_spill(ntp, nsp, helio_tp1P, helio_tpd1P, helio_tpspP)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE helio_discard_spill(ntp, nsp, helio_tp1P, helio_tpd1P, helio_tpspP)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_discard_spill
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
     TYPE(helio_tp), POINTER     :: helio_tp1P, helio_tpd1P, helio_tpspP

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpspP
     TYPE(helio_tp), POINTER   :: helio_tpP

! Executable code
     swifter_tpspP => helio_tpspP%swifter
     IF (nsp == 0) THEN
          helio_tpd1P => helio_tpspP
     ELSE
          helio_tpP => helio_tpd1P
          DO i = 1, nsp - 1
               helio_tpP => helio_tpP%nextP
          END DO
          helio_tpP%nextP => helio_tpspP
          helio_tpP%swifter%nextP => swifter_tpspP
     END IF
     IF (ASSOCIATED(helio_tpspP%prevP)) THEN
          helio_tpspP%prevP%nextP => helio_tpspP%nextP
          swifter_tpspP%prevP%nextP => swifter_tpspP%nextP
     ELSE
          helio_tp1P => helio_tpspP%nextP
     END IF
     IF (ASSOCIATED(helio_tpspP%nextP)) THEN
          helio_tpspP%nextP%prevP => helio_tpspP%prevP
          swifter_tpspP%nextP%prevP => swifter_tpspP%prevP
     END IF
     IF (nsp == 0) THEN
          NULLIFY(helio_tpspP%prevP)
          NULLIFY(swifter_tpspP%prevP)
     ELSE
          helio_tpspP%prevP => helio_tpP
          swifter_tpspP%prevP => helio_tpP%swifter
     END IF
     NULLIFY(helio_tpspP%nextP)
     NULLIFY(swifter_tpspP%nextP)
     nsp = nsp + 1
     ntp = ntp - 1

     RETURN

END SUBROUTINE helio_discard_spill
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
