!**********************************************************************************************************************************
!
!  Unit Name   : util_valid
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Validate planet and test particle ids
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!                swifter_tp1P : pointer to head of active SWIFTER test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL util_valid(npl, ntp, swifter_pl1P, swifter_tp1P)
!
!  Notes       : Subroutine causes program to exit with error if any ids are not unique
!
!**********************************************************************************************************************************
SUBROUTINE util_valid(npl, ntp, swifter_pl1P, swifter_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => util_valid
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl, ntp
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tp1P

! Internals
     INTEGER(I4B)                            :: i
     INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: idarr
     TYPE(swifter_pl), POINTER               :: swifter_plP
     TYPE(swifter_tp), POINTER               :: swifter_tpP

! Executable code
     ALLOCATE(idarr(npl+ntp))
     swifter_plP => swifter_pl1P
     DO i = 1, npl
          idarr(i) = swifter_plP%id
          swifter_plP => swifter_plP%nextP
     END DO
     swifter_tpP => swifter_tp1P
     DO i = npl + 1, npl + ntp
          idarr(i) = swifter_tpP%id
          swifter_tpP => swifter_tpP%nextP
     END DO
     CALL util_sort(idarr)
     DO i = 1, npl + ntp - 1
          IF (idarr(i) == idarr(i+1)) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   More than one body/particle has id = ", idarr(i)
               CALL util_exit(FAILURE)
          END IF
     END DO
     DEALLOCATE(idarr)

     RETURN

END SUBROUTINE util_valid
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
