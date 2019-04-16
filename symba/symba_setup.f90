!**********************************************************************************************************************************
!
!  Unit Name   : symba_setup
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Set up pointers within SyMBA, helio and Swifter planet and test particle structure linked-lists
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                symba_plA    : SyMBA planet structure array
!                symba_tpA    : SyMBA test particle structure array
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_plA    : SyMBA planet structure array
!                symba_tpA    : SyMBA test particle structure array
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P   : pointer to head of active SyMBA test particle structure linked-list
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1P, symba_tp1P, swifter_pl1P, swifter_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1P, symba_tp1P, swifter_pl1P, swifter_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_setup
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                            :: npl, ntp
     TYPE(swifter_pl), POINTER                           :: swifter_pl1P
     TYPE(swifter_tp), POINTER                           :: swifter_tp1P
     TYPE(symba_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: symba_plA
     TYPE(symba_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: symba_tpA
     TYPE(symba_pl), POINTER                             :: symba_pl1P
     TYPE(symba_tp), POINTER                             :: symba_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(helio_pl), POINTER   :: helio_pl1P, helio_plP
     TYPE(helio_tp), POINTER   :: helio_tp1P, helio_tpP
     TYPE(symba_pl), POINTER   :: symba_plP
     TYPE(symba_tp), POINTER   :: symba_tpP

! Executable code
     symba_pl1P => symba_plA(1)
     helio_pl1P => symba_plA(1)%helio
     swifter_pl1P => symba_plA(1)%helio%swifter
     NULLIFY(symba_pl1P%prevP)
     NULLIFY(helio_pl1P%prevP)
     NULLIFY(swifter_pl1P%prevP)
     IF (npl == 1) THEN
          NULLIFY(symba_pl1P%nextP)
          NULLIFY(helio_pl1P%nextP)
          NULLIFY(swifter_pl1P%nextP)
     ELSE
          symba_pl1P%nextP => symba_plA(2)
          helio_pl1P%nextP => symba_plA(2)%helio
          swifter_pl1P%nextP => symba_plA(2)%helio%swifter
          DO i = 2, npl - 1
               symba_plA(i)%prevP => symba_plA(i-1)
               symba_plA(i)%nextP => symba_plA(i+1)
               helio_plP => symba_plA(i)%helio
               helio_plP%prevP => symba_plA(i-1)%helio
               helio_plP%nextP => symba_plA(i+1)%helio
               swifter_plP => symba_plA(i)%helio%swifter
               swifter_plP%prevP => symba_plA(i-1)%helio%swifter
               swifter_plP%nextP => symba_plA(i+1)%helio%swifter
          END DO
          symba_plA(npl)%prevP => symba_plA(npl-1)
          symba_plP => symba_plA(npl)
          NULLIFY(symba_plP%nextP)
          helio_plP => symba_plA(npl)%helio
          helio_plP%prevP => symba_plA(npl-1)%helio
          NULLIFY(helio_plP%nextP)
          swifter_plP => symba_plA(npl)%helio%swifter
          swifter_plP%prevP => symba_plA(npl-1)%helio%swifter
          NULLIFY(swifter_plP%nextP)
     END IF
     NULLIFY(symba_tp1P)
     NULLIFY(helio_tp1P)
     NULLIFY(swifter_tp1P)
     IF (ntp > 0) THEN
          symba_tp1P => symba_tpA(1)
          helio_tp1P => symba_tpA(1)%helio
          swifter_tp1P => symba_tpA(1)%helio%swifter
          NULLIFY(symba_tp1P%prevP)
          NULLIFY(helio_tp1P%prevP)
          NULLIFY(swifter_tp1P%prevP)
          IF (ntp == 1) THEN
               NULLIFY(symba_tp1P%nextP)
               NULLIFY(helio_tp1P%nextP)
               NULLIFY(swifter_tp1P%nextP)
          ELSE
               symba_tp1P%nextP => symba_tpA(2)
               helio_tp1P%nextP => symba_tpA(2)%helio
               swifter_tp1P%nextP => symba_tpA(2)%helio%swifter
               DO i = 2, ntp - 1
                    symba_tpA(i)%prevP => symba_tpA(i-1)
                    symba_tpA(i)%nextP => symba_tpA(i+1)
                    helio_tpP => symba_tpA(i)%helio
                    helio_tpP%prevP => symba_tpA(i-1)%helio
                    helio_tpP%nextP => symba_tpA(i+1)%helio
                    swifter_tpP => symba_tpA(i)%helio%swifter
                    swifter_tpP%prevP => symba_tpA(i-1)%helio%swifter
                    swifter_tpP%nextP => symba_tpA(i+1)%helio%swifter
               END DO
               symba_tpA(ntp)%prevP => symba_tpA(ntp-1)
               symba_tpP => symba_tpA(ntp)
               NULLIFY(symba_tpP%nextP)
               helio_tpP => symba_tpA(ntp)%helio
               helio_tpP%prevP => symba_tpA(ntp-1)%helio
               NULLIFY(helio_tpP%nextP)
               swifter_tpP => symba_tpA(ntp)%helio%swifter
               swifter_tpP%prevP => symba_tpA(ntp-1)%helio%swifter
               NULLIFY(swifter_tpP%nextP)
          END IF
     END IF

     RETURN

END SUBROUTINE symba_setup
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
