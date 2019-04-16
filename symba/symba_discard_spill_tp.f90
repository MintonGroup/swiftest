!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_spill_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Move spilled (discarded) SyMBA test particle structure from active list to discard list
!
!  Input
!    Arguments : ntp         : number of active test particles
!                nsp         : number of spilled test particles
!                symba_tp1P  : pointer to head of active SyMBA test particle structure linked-list
!                symba_tpd1P : pointer to head of discard SyMBA test particle structure linked-list
!                symba_tpspP : pointer to SyMBA test particle structure to be discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ntp         : number of active test particles
!                nsp         : number of spilled test particles
!                symba_tp1P  : pointer to head of active SyMBA test particle structure linked-list
!                symba_tpd1P : pointer to head of discard SyMBA test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_discard_spill_tp(ntp, nsp, symba_tp1P, symba_tpd1P, symba_tpspP)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_spill_tp(ntp, nsp, symba_tp1P, symba_tpd1P, symba_tpspP)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_spill_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
     TYPE(symba_tp), POINTER     :: symba_tp1P, symba_tpd1P, symba_tpspP

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpspP
     TYPE(helio_tp), POINTER   :: helio_tpspP
     TYPE(symba_tp), POINTER   :: symba_tpP

! Executable code
     helio_tpspP => symba_tpspP%helio
     swifter_tpspP => helio_tpspP%swifter
     IF (nsp == 0) THEN
          symba_tpd1P => symba_tpspP
     ELSE
          symba_tpP => symba_tpd1P
          DO i = 1, nsp - 1
               symba_tpP => symba_tpP%nextP
          END DO
          symba_tpP%nextP => symba_tpspP
          symba_tpP%helio%nextP => helio_tpspP
          symba_tpP%helio%swifter%nextP => swifter_tpspP
     END IF
     IF (ASSOCIATED(symba_tpspP%prevP)) THEN
          symba_tpspP%prevP%nextP => symba_tpspP%nextP
          helio_tpspP%prevP%nextP => helio_tpspP%nextP
          swifter_tpspP%prevP%nextP => swifter_tpspP%nextP
     ELSE
          symba_tp1P => symba_tpspP%nextP
     END IF
     IF (ASSOCIATED(symba_tpspP%nextP)) THEN
          symba_tpspP%nextP%prevP => symba_tpspP%prevP
          helio_tpspP%nextP%prevP => helio_tpspP%prevP
          swifter_tpspP%nextP%prevP => swifter_tpspP%prevP
     END IF
     IF (nsp == 0) THEN
          NULLIFY(symba_tpspP%prevP)
          NULLIFY(helio_tpspP%prevP)
          NULLIFY(swifter_tpspP%prevP)
     ELSE
          symba_tpspP%prevP => symba_tpP
          helio_tpspP%prevP => symba_tpP%helio
          swifter_tpspP%prevP => symba_tpP%helio%swifter
     END IF
     NULLIFY(symba_tpspP%nextP)
     NULLIFY(helio_tpspP%nextP)
     NULLIFY(swifter_tpspP%nextP)
     nsp = nsp + 1
     ntp = ntp - 1

     RETURN

END SUBROUTINE symba_discard_spill_tp
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
