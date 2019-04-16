!**********************************************************************************************************************************
!
!  Unit Name   : tu4_setup
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : tu4
!  Language    : Fortran 90/95
!
!  Description : Set up pointers within TU4 and SWIFTER planet and test particle structure linked-lists
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                tu4_plA      : TU4 planet structure array
!                tu4_tpA      : TU4 test particle structure array
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : tu4_plA      : TU4 planet structure array
!                tu4_tpA      : TU4 test particle structure array
!                tu4_pl1P     : pointer to head of TU4 planet structure linked-list
!                tu4_tp1P     : pointer to head of active TU4 test particle structure linked-list
!                swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!                swifter_tp1P : pointer to head of active SWIFTER test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL tu4_setup(npl, ntp, tu4_plA, tu4_tpA, tu4_pl1P, tu4_tp1P, swifter_pl1P, swifter_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE tu4_setup(npl, ntp, tu4_plA, tu4_tpA, tu4_pl1P, tu4_tp1P, swifter_pl1P, swifter_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_tu4
     USE module_interfaces, EXCEPT_THIS_ONE => tu4_setup
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                          :: npl, ntp
     TYPE(swifter_pl), POINTER                         :: swifter_pl1P
     TYPE(swifter_tp), POINTER                         :: swifter_tp1P
     TYPE(tu4_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: tu4_plA
     TYPE(tu4_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: tu4_tpA
     TYPE(tu4_pl), POINTER                             :: tu4_pl1P
     TYPE(tu4_tp), POINTER                             :: tu4_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(tu4_pl), POINTER     :: tu4_plP
     TYPE(tu4_tp), POINTER     :: tu4_tpP

! Executable code
     tu4_pl1P => tu4_plA(1)
     swifter_pl1P => tu4_plA(1)%swifter
     NULLIFY(tu4_pl1P%prevP)
     NULLIFY(swifter_pl1P%prevP)
     IF (npl == 1) THEN
          NULLIFY(tu4_pl1P%nextP)
          NULLIFY(swifter_pl1P%nextP)
     ELSE
          tu4_pl1P%nextP => tu4_plA(2)
          swifter_pl1P%nextP => tu4_plA(2)%swifter
          DO i = 2, npl - 1
               tu4_plA(i)%prevP => tu4_plA(i-1)
               tu4_plA(i)%nextP => tu4_plA(i+1)
               swifter_plP => tu4_plA(i)%swifter
               swifter_plP%prevP => tu4_plA(i-1)%swifter
               swifter_plP%nextP => tu4_plA(i+1)%swifter
          END DO
          tu4_plA(npl)%prevP => tu4_plA(npl-1)
          tu4_plP => tu4_plA(npl)
          NULLIFY(tu4_plP%nextP)
          swifter_plP => tu4_plA(npl)%swifter
          swifter_plP%prevP => tu4_plA(npl-1)%swifter
          NULLIFY(swifter_plP%nextP)
     END IF
     NULLIFY(tu4_tp1P)
     NULLIFY(swifter_tp1P)
     IF (ntp > 0) THEN
          tu4_tp1P => tu4_tpA(1)
          swifter_tp1P => tu4_tpA(1)%swifter
          NULLIFY(tu4_tp1P%prevP)
          NULLIFY(swifter_tp1P%prevP)
          IF (ntp == 1) THEN
               NULLIFY(tu4_tp1P%nextP)
               NULLIFY(swifter_tp1P%nextP)
          ELSE
               tu4_tp1P%nextP => tu4_tpA(2)
               swifter_tp1P%nextP => tu4_tpA(2)%swifter
               DO i = 2, ntp - 1
                    tu4_tpA(i)%prevP => tu4_tpA(i-1)
                    tu4_tpA(i)%nextP => tu4_tpA(i+1)
                    swifter_tpP => tu4_tpA(i)%swifter
                    swifter_tpP%prevP => tu4_tpA(i-1)%swifter
                    swifter_tpP%nextP => tu4_tpA(i+1)%swifter
               END DO
               tu4_tpA(ntp)%prevP => tu4_tpA(ntp-1)
               tu4_tpP => tu4_tpA(ntp)
               NULLIFY(tu4_tpP%nextP)
               swifter_tpP => tu4_tpA(ntp)%swifter
               swifter_tpP%prevP => tu4_tpA(ntp-1)%swifter
               NULLIFY(swifter_tpP%nextP)
          END IF
     END IF

     RETURN

END SUBROUTINE tu4_setup
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
