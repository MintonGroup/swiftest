!**********************************************************************************************************************************
!
!  Unit Name   : util_toupper
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Convert string to uppercase
!
!  Input
!    Arguments : string : string to convert
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : string : converted string
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_toupper(string)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE util_toupper(string)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => util_toupper
     IMPLICIT NONE

! Arguments
     CHARACTER(*), INTENT(INOUT) :: string

! Internals
     INTEGER(I4B) :: i, length, index

! Executable code
     length = LEN(string)
     DO i = 1, length
          index = IACHAR(string(i:i))
          IF ((index >= LOWERCASE_BEGIN) .AND. (index <= LOWERCASE_END)) THEN
               index = index + UPPERCASE_OFFSET
               string(i:i) = ACHAR(index)
          END IF
     END DO

     RETURN

END SUBROUTINE util_toupper
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
