!**********************************************************************************************************************************
!
!  Unit Name   : io_get_token
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Return the next token from a character buffer
!
!  Input
!    Arguments : buffer  : input character buffer
!                ilength : buffer length
!                ifirst  : index in character buffer to begin token search
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ifirst  : index in character buffer to begin token search (returned as index of first character of found token)
!                ilast   : index in character buffer of last character of token
!                ierr    : error status (0 = OK, -1 = NO TOKEN FOUND)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL io_get_token(buffer, ilength, ifirst, ilast, ierr)
!
!  Notes       : Here a token is defined as any set of contiguous non-blank characters not beginning with or containing "!"
!                If "!" is present, any remaining part of the buffer including the "!" is ignored
!
!**********************************************************************************************************************************
SUBROUTINE io_get_token(buffer, ilength, ifirst, ilast, ierr)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => io_get_token
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)    :: ilength
     INTEGER(I4B), INTENT(INOUT) :: ifirst
     INTEGER(I4B), INTENT(OUT)   :: ilast, ierr
     CHARACTER(*), INTENT(IN)    :: buffer

! Internals
     INTEGER(I4B) :: i

! Executable code
     IF (ifirst > ilength) THEN
          ilast = ifirst
          ierr = -1
          RETURN
     END IF
     DO i = ifirst, ilength
          IF (buffer(i:i) /= " ") EXIT
     END DO
     IF ((i > ilength) .OR. (buffer(i:i) == "!")) THEN
          ifirst = i
          ilast = i
          ierr = -1
          RETURN
     END IF
     ifirst = i
     DO i = ifirst, ilength
          IF ((buffer(i:i) == " ") .OR. (buffer(i:i) == "!")) EXIT
     END DO
     ilast = i - 1
     ierr = 0

     RETURN

END SUBROUTINE io_get_token
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
