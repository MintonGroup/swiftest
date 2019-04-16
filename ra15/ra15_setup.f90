!**********************************************************************************************************************************
!
!  Unit Name   : ra15_setup
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : ra15
!  Language    : Fortran 90/95
!
!  Description : Set up pointers within RA15 and Swifter planet and test particle structure linked-lists
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                ra15_plA     : RA15 planet structure array
!                ra15_tpA     : RA15 test particle structure array
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ra15_plA     : RA15 planet structure array
!                ra15_tpA     : RA15 test particle structure array
!                ra15_pl1P    : pointer to head of RA15 planet structure linked-list
!                ra15_tp1P    : pointer to head of active RA15 test particle structure linked-list
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL ra15_setup(npl, ntp, ra15_plA, ra15_tpA, ra15_pl1P, ra15_tp1P, swifter_pl1P, swifter_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE ra15_setup(npl, ntp, ra15_plA, ra15_tpA, ra15_pl1P, ra15_tp1P, swifter_pl1P, swifter_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_ra15
     USE module_interfaces, EXCEPT_THIS_ONE => ra15_setup
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                           :: npl, ntp
     TYPE(swifter_pl), POINTER                          :: swifter_pl1P
     TYPE(swifter_tp), POINTER                          :: swifter_tp1P
     TYPE(ra15_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: ra15_plA
     TYPE(ra15_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: ra15_tpA
     TYPE(ra15_pl), POINTER                             :: ra15_pl1P
     TYPE(ra15_tp), POINTER                             :: ra15_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(ra15_pl), POINTER    :: ra15_plP
     TYPE(ra15_tp), POINTER    :: ra15_tpP

! Executable code
     ra15_pl1P => ra15_plA(1)
     swifter_pl1P => ra15_plA(1)%swifter
     NULLIFY(ra15_pl1P%prevP)
     NULLIFY(swifter_pl1P%prevP)
     IF (npl == 1) THEN
          NULLIFY(ra15_pl1P%nextP)
          NULLIFY(swifter_pl1P%nextP)
     ELSE
          ra15_pl1P%nextP => ra15_plA(2)
          swifter_pl1P%nextP => ra15_plA(2)%swifter
          DO i = 2, npl - 1
               ra15_plA(i)%prevP => ra15_plA(i-1)
               ra15_plA(i)%nextP => ra15_plA(i+1)
               swifter_plP => ra15_plA(i)%swifter
               swifter_plP%prevP => ra15_plA(i-1)%swifter
               swifter_plP%nextP => ra15_plA(i+1)%swifter
          END DO
          ra15_plA(npl)%prevP => ra15_plA(npl-1)
          ra15_plP => ra15_plA(npl)
          NULLIFY(ra15_plP%nextP)
          swifter_plP => ra15_plA(npl)%swifter
          swifter_plP%prevP => ra15_plA(npl-1)%swifter
          NULLIFY(swifter_plP%nextP)
     END IF
     NULLIFY(ra15_tp1P)
     NULLIFY(swifter_tp1P)
     IF (ntp > 0) THEN
          ra15_tp1P => ra15_tpA(1)
          swifter_tp1P => ra15_tpA(1)%swifter
          NULLIFY(ra15_tp1P%prevP)
          NULLIFY(swifter_tp1P%prevP)
          IF (ntp == 1) THEN
               NULLIFY(ra15_tp1P%nextP)
               NULLIFY(swifter_tp1P%nextP)
          ELSE
               ra15_tp1P%nextP => ra15_tpA(2)
               swifter_tp1P%nextP => ra15_tpA(2)%swifter
               DO i = 2, ntp - 1
                    ra15_tpA(i)%prevP => ra15_tpA(i-1)
                    ra15_tpA(i)%nextP => ra15_tpA(i+1)
                    swifter_tpP => ra15_tpA(i)%swifter
                    swifter_tpP%prevP => ra15_tpA(i-1)%swifter
                    swifter_tpP%nextP => ra15_tpA(i+1)%swifter
               END DO
               ra15_tpA(ntp)%prevP => ra15_tpA(ntp-1)
               ra15_tpP => ra15_tpA(ntp)
               NULLIFY(ra15_tpP%nextP)
               swifter_tpP => ra15_tpA(ntp)%swifter
               swifter_tpP%prevP => ra15_tpA(ntp-1)%swifter
               NULLIFY(swifter_tpP%nextP)
          END IF
     END IF

     RETURN

END SUBROUTINE ra15_setup
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
