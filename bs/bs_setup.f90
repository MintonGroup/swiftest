!**********************************************************************************************************************************
!
!  Unit Name   : bs_setup
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Set up pointers within BS and Swifter planet and test particle structure linked-lists
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                bs_plA       : BS planet structure array
!                bs_tpA       : BS test particle structure array
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : bs_plA       : BS planet structure array
!                bs_tpA       : BS test particle structure array
!                bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL bs_setup(npl, ntp, bs_plA, bs_tpA, bs_pl1P, bs_tp1P, swifter_pl1P, swifter_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE bs_setup(npl, ntp, bs_plA, bs_tpA, bs_pl1P, bs_tp1P, swifter_pl1P, swifter_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_setup
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: npl, ntp
     TYPE(swifter_pl), POINTER                        :: swifter_pl1P
     TYPE(swifter_tp), POINTER                        :: swifter_tp1P
     TYPE(bs_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: bs_plA
     TYPE(bs_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: bs_tpA
     TYPE(bs_pl), POINTER                             :: bs_pl1P
     TYPE(bs_tp), POINTER                             :: bs_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(bs_pl), POINTER      :: bs_plP
     TYPE(bs_tp), POINTER      :: bs_tpP

! Executable code
     bs_pl1P => bs_plA(1)
     swifter_pl1P => bs_plA(1)%swifter
     NULLIFY(bs_pl1P%prevP)
     NULLIFY(swifter_pl1P%prevP)
     IF (npl == 1) THEN
          NULLIFY(bs_pl1P%nextP)
          NULLIFY(swifter_pl1P%nextP)
     ELSE
          bs_pl1P%nextP => bs_plA(2)
          swifter_pl1P%nextP => bs_plA(2)%swifter
          DO i = 2, npl - 1
               bs_plA(i)%prevP => bs_plA(i-1)
               bs_plA(i)%nextP => bs_plA(i+1)
               swifter_plP => bs_plA(i)%swifter
               swifter_plP%prevP => bs_plA(i-1)%swifter
               swifter_plP%nextP => bs_plA(i+1)%swifter
          END DO
          bs_plA(npl)%prevP => bs_plA(npl-1)
          bs_plP => bs_plA(npl)
          NULLIFY(bs_plP%nextP)
          swifter_plP => bs_plA(npl)%swifter
          swifter_plP%prevP => bs_plA(npl-1)%swifter
          NULLIFY(swifter_plP%nextP)
     END IF
     NULLIFY(bs_tp1P)
     NULLIFY(swifter_tp1P)
     IF (ntp > 0) THEN
          bs_tp1P => bs_tpA(1)
          swifter_tp1P => bs_tpA(1)%swifter
          NULLIFY(bs_tp1P%prevP)
          NULLIFY(swifter_tp1P%prevP)
          IF (ntp == 1) THEN
               NULLIFY(bs_tp1P%nextP)
               NULLIFY(swifter_tp1P%nextP)
          ELSE
               bs_tp1P%nextP => bs_tpA(2)
               swifter_tp1P%nextP => bs_tpA(2)%swifter
               DO i = 2, ntp - 1
                    bs_tpA(i)%prevP => bs_tpA(i-1)
                    bs_tpA(i)%nextP => bs_tpA(i+1)
                    swifter_tpP => bs_tpA(i)%swifter
                    swifter_tpP%prevP => bs_tpA(i-1)%swifter
                    swifter_tpP%nextP => bs_tpA(i+1)%swifter
               END DO
               bs_tpA(ntp)%prevP => bs_tpA(ntp-1)
               bs_tpP => bs_tpA(ntp)
               NULLIFY(bs_tpP%nextP)
               swifter_tpP => bs_tpA(ntp)%swifter
               swifter_tpP%prevP => bs_tpA(ntp-1)%swifter
               NULLIFY(swifter_tpP%nextP)
          END IF
     END IF

     RETURN

END SUBROUTINE bs_setup
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
