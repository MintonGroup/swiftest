!**********************************************************************************************************************************
!
!  Unit Name   : util_toupper
!  Unit Type   : subroutine
!  Project     : Swiftest
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
     INTEGER(I4B) :: i, length, idx

! Executable code
     length = LEN(string)
     DO i = 1, length
          idx = IACHAR(string(i:i))
          IF ((idx >= LOWERCASE_BEGIN) .AND. (idx <= LOWERCASE_END)) THEN
               idx = idx + UPPERCASE_OFFSET
               string(i:i) = ACHAR(idx)
          END IF
     END DO

     RETURN

END SUBROUTINE util_toupper
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
