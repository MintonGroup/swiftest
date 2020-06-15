!**********************************************************************************************************************************
!
!  Unit Name   : util_valid
!  Unit Type   : subroutine
!  Project     : Swiftest
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
SUBROUTINE util_valid(npl, ntp, swiftest_plA, swiftest_tpA)

! Modules
     use swiftest, EXCEPT_THIS_ONE => util_valid
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)        :: npl, ntp
     TYPE(swiftest_pl), INTENT(IN) :: swiftest_plA
     TYPE(swiftest_tp), INTENT(IN) :: swiftest_tpA

! Internals
     INTEGER(I4B)                            :: i
     INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: idarr

! Executable code
     ALLOCATE(idarr(npl+ntp))
     DO i = 1, npl
          idarr(i) = swiftest_plA%name(i)
     END DO
     DO i = 1, ntp
          idarr(npl+i) = swiftest_tpA%name(i)
     END DO
     CALL util_sort(idarr)
     DO i = 1, npl + ntp - 1
          IF (idarr(i) == idarr(i+1)) THEN
               WRITE(*, *) "SWIFTEST Error:"
               WRITE(*, *) "   More than one body/particle has id = ", idarr(i)
               CALL util_exit(FAILURE)
          END IF
     END DO
     DEALLOCATE(idarr)

     RETURN

END SUBROUTINE util_valid
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann (Checked by Jennifer Pouplin & Carlisle Wishard)
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
