!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_setup
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Set up pointers within RMVS, WHM and Swifter planet and test particle structure linked-lists
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                rmvs_plA     : RMVS planet structure array
!                rmvs_tpA     : RMVS test particle structure array
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_plA     : RMVS planet structure array
!                rmvs_tpA     : RMVS test particle structure array
!                rmvs_pl1P    : pointer to head of RMVS planet structure linked-list
!                rmvs_tp1P    : pointer to head of active RMVS test particle structure linked-list
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_setup(npl, ntp, rmvs_plA, rmvs_tpA, rmvs_pl1P, rmvs_tp1P, swifter_pl1P, swifter_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_setup(npl, ntp, rmvs_plA, rmvs_tpA, rmvs_pl1P, rmvs_tp1P, swifter_pl1P, swifter_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_setup
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                           :: npl, ntp
     TYPE(swifter_pl), POINTER                          :: swifter_pl1P
     TYPE(swifter_tp), POINTER                          :: swifter_tp1P
     TYPE(rmvs_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: rmvs_plA
     TYPE(rmvs_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: rmvs_tpA
     TYPE(rmvs_pl), POINTER                             :: rmvs_pl1P
     TYPE(rmvs_tp), POINTER                             :: rmvs_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(whm_pl), POINTER     :: whm_pl1P, whm_plP
     TYPE(whm_tp), POINTER     :: whm_tp1P, whm_tpP
     TYPE(rmvs_pl), POINTER    :: rmvs_plP
     TYPE(rmvs_tp), POINTER    :: rmvs_tpP

! Executable code
     rmvs_pl1P => rmvs_plA(1)
     whm_pl1P => rmvs_plA(1)%whm
     swifter_pl1P => rmvs_plA(1)%whm%swifter
     NULLIFY(rmvs_pl1P%prevP)
     NULLIFY(whm_pl1P%prevP)
     NULLIFY(swifter_pl1P%prevP)
     IF (npl == 1) THEN
          NULLIFY(rmvs_pl1P%nextP)
          NULLIFY(whm_pl1P%nextP)
          NULLIFY(swifter_pl1P%nextP)
     ELSE
          rmvs_pl1P%nextP => rmvs_plA(2)
          whm_pl1P%nextP => rmvs_plA(2)%whm
          swifter_pl1P%nextP => rmvs_plA(2)%whm%swifter
          DO i = 2, npl - 1
               rmvs_plA(i)%prevP => rmvs_plA(i-1)
               rmvs_plA(i)%nextP => rmvs_plA(i+1)
               whm_plP => rmvs_plA(i)%whm
               whm_plP%prevP => rmvs_plA(i-1)%whm
               whm_plP%nextP => rmvs_plA(i+1)%whm
               swifter_plP => rmvs_plA(i)%whm%swifter
               swifter_plP%prevP => rmvs_plA(i-1)%whm%swifter
               swifter_plP%nextP => rmvs_plA(i+1)%whm%swifter
          END DO
          rmvs_plA(npl)%prevP => rmvs_plA(npl-1)
          rmvs_plP => rmvs_plA(npl)
          NULLIFY(rmvs_plP%nextP)
          whm_plP => rmvs_plA(npl)%whm
          whm_plP%prevP => rmvs_plA(npl-1)%whm
          NULLIFY(whm_plP%nextP)
          swifter_plP => rmvs_plA(npl)%whm%swifter
          swifter_plP%prevP => rmvs_plA(npl-1)%whm%swifter
          NULLIFY(swifter_plP%nextP)
     END IF
     NULLIFY(rmvs_tp1P)
     NULLIFY(whm_tp1P)
     NULLIFY(swifter_tp1P)
     IF (ntp > 0) THEN
          rmvs_tp1P => rmvs_tpA(1)
          whm_tp1P => rmvs_tpA(1)%whm
          swifter_tp1P => rmvs_tpA(1)%whm%swifter
          NULLIFY(rmvs_tp1P%prevP)
          NULLIFY(whm_tp1P%prevP)
          NULLIFY(swifter_tp1P%prevP)
          IF (ntp == 1) THEN
               NULLIFY(rmvs_tp1P%nextP)
               NULLIFY(whm_tp1P%nextP)
               NULLIFY(swifter_tp1P%nextP)
          ELSE
               rmvs_tp1P%nextP => rmvs_tpA(2)
               whm_tp1P%nextP => rmvs_tpA(2)%whm
               swifter_tp1P%nextP => rmvs_tpA(2)%whm%swifter
               DO i = 2, ntp - 1
                    rmvs_tpA(i)%prevP => rmvs_tpA(i-1)
                    rmvs_tpA(i)%nextP => rmvs_tpA(i+1)
                    whm_tpP => rmvs_tpA(i)%whm
                    whm_tpP%prevP => rmvs_tpA(i-1)%whm
                    whm_tpP%nextP => rmvs_tpA(i+1)%whm
                    swifter_tpP => rmvs_tpA(i)%whm%swifter
                    swifter_tpP%prevP => rmvs_tpA(i-1)%whm%swifter
                    swifter_tpP%nextP => rmvs_tpA(i+1)%whm%swifter
               END DO
               rmvs_tpA(ntp)%prevP => rmvs_tpA(ntp-1)
               rmvs_tpP => rmvs_tpA(ntp)
               NULLIFY(rmvs_tpP%nextP)
               whm_tpP => rmvs_tpA(ntp)%whm
               whm_tpP%prevP => rmvs_tpA(ntp-1)%whm
               NULLIFY(whm_tpP%nextP)
               swifter_tpP => rmvs_tpA(ntp)%whm%swifter
               swifter_tpP%prevP => rmvs_tpA(ntp-1)%whm%swifter
               NULLIFY(swifter_tpP%nextP)
          END IF
     END IF

     RETURN

END SUBROUTINE rmvs_setup
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
