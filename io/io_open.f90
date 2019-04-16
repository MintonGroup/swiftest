!**********************************************************************************************************************************
!
!  Unit Name   : io_open
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Open a file
!
!  Input
!    Arguments : iu        : unit number with which to open file
!                fname     : name of file to open
!                fopenstat : status with which to open file
!                fmt       : format string
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ierr      : open error status (0 = OK, nonzero = ERROR)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL io_open(iu, fname, fopenstat, fmt, ierr)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_open.F
!
!**********************************************************************************************************************************
SUBROUTINE io_open(iu, fname, fopenstat, fmt, ierr)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => io_open
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: iu
     INTEGER(I4B), INTENT(OUT) :: ierr
     CHARACTER(*), INTENT(IN)  :: fname, fopenstat, fmt

! Internals
     CHARACTER(STRMAX) :: fostat

! Executable code
     fostat = TRIM(ADJUSTL(fopenstat))
     CALL util_toupper(fostat)
     IF (fostat == "APPEND") THEN
          OPEN(UNIT = iu, FILE = fname, STATUS = "OLD", POSITION = fostat, FORM = fmt, IOSTAT = ierr)
     ELSE
          OPEN(UNIT = iu, FILE = fname, STATUS = fostat, FORM = fmt, IOSTAT = ierr)
     END IF

     RETURN

END SUBROUTINE io_open
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
