!**********************************************************************************************************************************
!
!  Unit Name   : helio_setup
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Set up pointers within helio and Swifter planet and test particle structure linked-lists
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                helio_plA    : helio planet structure array
!                helio_tpA    : helio test particle structure array
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_plA    : helio planet structure array
!                helio_tpA    : helio test particle structure array
!                helio_pl1P   : pointer to head of helio planet structure linked-list
!                helio_tp1P   : pointer to head of active helio test particle structure linked-list
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_setup(npl, ntp, helio_plA, helio_tpA, helio_pl1P, helio_tp1P, swifter_pl1P, swifter_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE helio_setup(npl, ntp, helio_plA, helio_tpA, helio_pl1P, helio_tp1P, swifter_pl1P, swifter_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_setup
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                            :: npl, ntp
     TYPE(swifter_pl), POINTER                           :: swifter_pl1P
     TYPE(swifter_tp), POINTER                           :: swifter_tp1P
     TYPE(helio_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: helio_plA
     TYPE(helio_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: helio_tpA
     TYPE(helio_pl), POINTER                             :: helio_pl1P
     TYPE(helio_tp), POINTER                             :: helio_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(helio_pl), POINTER   :: helio_plP
     TYPE(helio_tp), POINTER   :: helio_tpP

! Executable code
     helio_pl1P => helio_plA(1)
     swifter_pl1P => helio_plA(1)%swifter
     NULLIFY(helio_pl1P%prevP)
     NULLIFY(swifter_pl1P%prevP)
     IF (npl == 1) THEN
          NULLIFY(helio_pl1P%nextP)
          NULLIFY(swifter_pl1P%nextP)
     ELSE
          helio_pl1P%nextP => helio_plA(2)
          swifter_pl1P%nextP => helio_plA(2)%swifter
          DO i = 2, npl - 1
               helio_plA(i)%prevP => helio_plA(i-1)
               helio_plA(i)%nextP => helio_plA(i+1)
               swifter_plP => helio_plA(i)%swifter
               swifter_plP%prevP => helio_plA(i-1)%swifter
               swifter_plP%nextP => helio_plA(i+1)%swifter
          END DO
          helio_plA(npl)%prevP => helio_plA(npl-1)
          helio_plP => helio_plA(npl)
          NULLIFY(helio_plP%nextP)
          swifter_plP => helio_plA(npl)%swifter
          swifter_plP%prevP => helio_plA(npl-1)%swifter
          NULLIFY(swifter_plP%nextP)
     END IF
     NULLIFY(helio_tp1P)
     NULLIFY(swifter_tp1P)
     IF (ntp > 0) THEN
          helio_tp1P => helio_tpA(1)
          swifter_tp1P => helio_tpA(1)%swifter
          NULLIFY(helio_tp1P%prevP)
          NULLIFY(swifter_tp1P%prevP)
          IF (ntp == 1) THEN
               NULLIFY(helio_tp1P%nextP)
               NULLIFY(swifter_tp1P%nextP)
          ELSE
               helio_tp1P%nextP => helio_tpA(2)
               swifter_tp1P%nextP => helio_tpA(2)%swifter
               DO i = 2, ntp - 1
                    helio_tpA(i)%prevP => helio_tpA(i-1)
                    helio_tpA(i)%nextP => helio_tpA(i+1)
                    swifter_tpP => helio_tpA(i)%swifter
                    swifter_tpP%prevP => helio_tpA(i-1)%swifter
                    swifter_tpP%nextP => helio_tpA(i+1)%swifter
               END DO
               helio_tpA(ntp)%prevP => helio_tpA(ntp-1)
               helio_tpP => helio_tpA(ntp)
               NULLIFY(helio_tpP%nextP)
               swifter_tpP => helio_tpA(ntp)%swifter
               swifter_tpP%prevP => helio_tpA(ntp-1)%swifter
               NULLIFY(swifter_tpP%nextP)
          END IF
     END IF

     RETURN

END SUBROUTINE helio_setup
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
