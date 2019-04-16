!**********************************************************************************************************************************
!
!  Unit Name   : ra15_discard_spill
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : ra15
!  Language    : Fortran 90/95
!
!  Description : Move spilled (discarded) RA15 test particle structure from active list to discard list
!
!  Input
!    Arguments : ntp        : number of active test particles
!                nsp        : number of spilled test particles
!                ra15_tp1P  : pointer to head of active RA15 test particle structure linked-list
!                ra15_tpd1P : pointer to head of discard RA15 test particle structure linked-list
!                ra15_tpspP : pointer to RA15 test particle structure to be discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ntp        : number of active test particles
!                nsp        : number of spilled test particles
!                ra15_tp1P  : pointer to head of active RA15 test particle structure linked-list
!                ra15_tpd1P : pointer to head of discard RA15 test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL ra15_discard_spill(ntp, nsp, ra15_tp1P, ra15_tpd1P, ra15_tpspP)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE ra15_discard_spill(ntp, nsp, ra15_tp1P, ra15_tpd1P, ra15_tpspP)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_ra15
     USE module_interfaces, EXCEPT_THIS_ONE => ra15_discard_spill
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
     TYPE(ra15_tp), POINTER      :: ra15_tp1P, ra15_tpd1P, ra15_tpspP

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpspP
     TYPE(ra15_tp), POINTER    :: ra15_tpP

! Executable code
     swifter_tpspP => ra15_tpspP%swifter
     IF (nsp == 0) THEN
          ra15_tpd1P => ra15_tpspP
     ELSE
          ra15_tpP => ra15_tpd1P
          DO i = 1, nsp - 1
               ra15_tpP => ra15_tpP%nextP
          END DO
          ra15_tpP%nextP => ra15_tpspP
          ra15_tpP%swifter%nextP => swifter_tpspP
     END IF
     IF (ASSOCIATED(ra15_tpspP%prevP)) THEN
          ra15_tpspP%prevP%nextP => ra15_tpspP%nextP
          swifter_tpspP%prevP%nextP => swifter_tpspP%nextP
     ELSE
          ra15_tp1P => ra15_tpspP%nextP
     END IF
     IF (ASSOCIATED(ra15_tpspP%nextP)) THEN
          ra15_tpspP%nextP%prevP => ra15_tpspP%prevP
          swifter_tpspP%nextP%prevP => swifter_tpspP%prevP
     END IF
     IF (nsp == 0) THEN
          NULLIFY(ra15_tpspP%prevP)
          NULLIFY(swifter_tpspP%prevP)
     ELSE
          ra15_tpspP%prevP => ra15_tpP
          swifter_tpspP%prevP => ra15_tpP%swifter
     END IF
     NULLIFY(ra15_tpspP%nextP)
     NULLIFY(swifter_tpspP%nextP)
     nsp = nsp + 1
     ntp = ntp - 1

     RETURN

END SUBROUTINE ra15_discard_spill
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
