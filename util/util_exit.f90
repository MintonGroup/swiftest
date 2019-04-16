!**********************************************************************************************************************************
!
!  Unit Name   : util_exit
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Print termination message and exit program
!
!  Input
!    Arguments : code : exit code
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : termination message
!    File      : none
!
!  Invocation  : CALL util_exit(code)
!
!  Notes       : Adapted from Hal Levison's Swift routine util_exit.f
!
!**********************************************************************************************************************************
SUBROUTINE util_exit(code)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => util_exit
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: code

! Executable code
     IF (code == SUCCESS) THEN
          WRITE(*, 100) VERSION_NUMBER
     ELSE
          WRITE(*, 200) VERSION_NUMBER
     END IF
 100 FORMAT(/, "Normal termination of SWIFTER (Version ", F3.1, ")")
 200 FORMAT(/, "Terminating SWIFTER (Version ", F3.1, ") due to ERROR!!")
     WRITE(*, 300) "------------------------------------------------"
 300 FORMAT(A)

     STOP

END SUBROUTINE util_exit
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
