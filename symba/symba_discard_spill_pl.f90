!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_spill_pl
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Move spilled (discarded) SyMBA planet structure from active list to discard list
!
!  Input
!    Arguments : npl         : number of planets
!                nsp         : number of spilled planets
!                symba_pld1P : pointer to head of discard SyMBA planet structure linked-list
!                symba_plspP : pointer to SyMBA planet structure to be discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl         : number of planets
!                nsp         : number of spilled planets
!                symba_pld1P : pointer to head of discard SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_discard_spill_pl(npl, nsp, symba_pld1P, symba_plspP)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_spill_pl(npl, nsp, symba_pld1P, symba_plspP)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_spill_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: npl, nsp
     TYPE(symba_pl), POINTER     :: symba_pld1P, symba_plspP

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plspP
     TYPE(helio_pl), POINTER   :: helio_plspP
     TYPE(symba_pl), POINTER   :: symba_plP

! Executable code
     helio_plspP => symba_plspP%helio
     swifter_plspP => helio_plspP%swifter
     IF (nsp == 0) THEN
          symba_pld1P => symba_plspP
     ELSE
          symba_plP => symba_pld1P
          DO i = 1, nsp - 1
               symba_plP => symba_plP%nextP
          END DO
          symba_plP%nextP => symba_plspP
          symba_plP%helio%nextP => helio_plspP
          symba_plP%helio%swifter%nextP => swifter_plspP
     END IF
     symba_plspP%prevP%nextP => symba_plspP%nextP
     helio_plspP%prevP%nextP => helio_plspP%nextP
     swifter_plspP%prevP%nextP => swifter_plspP%nextP
     IF (ASSOCIATED(symba_plspP%nextP)) THEN
          symba_plspP%nextP%prevP => symba_plspP%prevP
          helio_plspP%nextP%prevP => helio_plspP%prevP
          swifter_plspP%nextP%prevP => swifter_plspP%prevP
     END IF
     IF (nsp == 0) THEN
          NULLIFY(symba_plspP%prevP)
          NULLIFY(helio_plspP%prevP)
          NULLIFY(swifter_plspP%prevP)
     ELSE
          symba_plspP%prevP => symba_plP
          helio_plspP%prevP => symba_plP%helio
          swifter_plspP%prevP => symba_plP%helio%swifter
     END IF
     NULLIFY(symba_plspP%nextP)
     NULLIFY(helio_plspP%nextP)
     NULLIFY(swifter_plspP%nextP)
     nsp = nsp + 1
     npl = npl - 1

     RETURN

END SUBROUTINE symba_discard_spill_pl
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
