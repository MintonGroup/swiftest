!**********************************************************************************************************************************
!
!  Unit Name   : whm_setup
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Set up pointers within WHM and SWIFTER planet and test particle structure linked-lists
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                whm_plA      : WHM planet structure array
!                whm_tpA      : WHM test particle structure array
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_plA      : WHM planet structure array
!                whm_tpA      : WHM test particle structure array
!                whm_pl1P     : pointer to head of WHM planet structure linked-list
!                whm_tp1P     : pointer to head of active WHM test particle structure linked-list
!                swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!                swifter_tp1P : pointer to head of active SWIFTER test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_setup(npl, ntp, whm_plA, whm_tpA, whm_pl1P, whm_tp1P, swifter_pl1P, swifter_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE whm_setup(npl, ntp, whm_plA, whm_tpA, whm_pl1P, whm_tp1P, swifter_pl1P, swifter_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_setup
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                          :: npl, ntp
     TYPE(swifter_pl), POINTER                         :: swifter_pl1P
     TYPE(swifter_tp), POINTER                         :: swifter_tp1P
     TYPE(whm_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: whm_plA
     TYPE(whm_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: whm_tpA
     TYPE(whm_pl), POINTER                             :: whm_pl1P
     TYPE(whm_tp), POINTER                             :: whm_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(whm_pl), POINTER     :: whm_plP
     TYPE(whm_tp), POINTER     :: whm_tpP
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     whm_pl1P => whm_plA(1)
     swifter_pl1P => whm_plA(1)%swifter
     NULLIFY(whm_pl1P%prevP)
     NULLIFY(swifter_pl1P%prevP)
     IF (npl == 1) THEN
          NULLIFY(whm_pl1P%nextP)
          NULLIFY(swifter_pl1P%nextP)
     ELSE
          whm_pl1P%nextP => whm_plA(2)
          swifter_pl1P%nextP => whm_plA(2)%swifter
          DO i = 2, npl - 1
               whm_plA(i)%prevP => whm_plA(i-1)
               whm_plA(i)%nextP => whm_plA(i+1)
               swifter_plP => whm_plA(i)%swifter
               swifter_plP%prevP => whm_plA(i-1)%swifter
               swifter_plP%nextP => whm_plA(i+1)%swifter
          END DO
          whm_plA(npl)%prevP => whm_plA(npl-1)
          whm_plP => whm_plA(npl)
          NULLIFY(whm_plP%nextP)
          swifter_plP => whm_plA(npl)%swifter
          swifter_plP%prevP => whm_plA(npl-1)%swifter
          NULLIFY(swifter_plP%nextP)
     END IF
     NULLIFY(whm_tp1P)
     NULLIFY(swifter_tp1P)
     IF (ntp > 0) THEN
          whm_tp1P => whm_tpA(1)
          swifter_tp1P => whm_tpA(1)%swifter
          NULLIFY(whm_tp1P%prevP)
          NULLIFY(swifter_tp1P%prevP)
          IF (ntp == 1) THEN
               NULLIFY(whm_tp1P%nextP)
               NULLIFY(swifter_tp1P%nextP)
          ELSE
               whm_tp1P%nextP => whm_tpA(2)
               swifter_tp1P%nextP => whm_tpA(2)%swifter
               DO i = 2, ntp - 1
                    whm_tpA(i)%prevP => whm_tpA(i-1)
                    whm_tpA(i)%nextP => whm_tpA(i+1)
                    swifter_tpP => whm_tpA(i)%swifter
                    swifter_tpP%prevP => whm_tpA(i-1)%swifter
                    swifter_tpP%nextP => whm_tpA(i+1)%swifter
               END DO
               whm_tpA(ntp)%prevP => whm_tpA(ntp-1)
               whm_tpP => whm_tpA(ntp)
               NULLIFY(whm_tpP%nextP)
               swifter_tpP => whm_tpA(ntp)%swifter
               swifter_tpP%prevP => whm_tpA(ntp-1)%swifter
               NULLIFY(swifter_tpP%nextP)
          END IF
     END IF

     RETURN

END SUBROUTINE whm_setup
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
