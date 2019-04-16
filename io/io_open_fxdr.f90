!**********************************************************************************************************************************
!
!  Unit Name   : io_open_fxdr
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Open an XDR formatted file using FXDR routines
!
!  Input
!    Arguments : fname     : name of file to open
!                fopenstat : status with which to open file
!                lflag     : logical flag indicating behavior to exhibit on error (.TRUE. = return, .FALSE. = stop)
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : iu        : unit number of open file
!                ierr      : open error status (0 = OK, nonzero = ERROR)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL io_open_fxdr(fname, fopenstat, lflag, iu, ierr)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_open_fxdr.F
!
!**********************************************************************************************************************************
SUBROUTINE io_open_fxdr(fname, fopenstat, lflag, iu, ierr)

! Modules
     USE module_parameters
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_open_fxdr
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)  :: lflag
     INTEGER(I4B), INTENT(OUT) :: iu, ierr
     CHARACTER(*), INTENT(IN)  :: fname
     CHARACTER(1), INTENT(IN)  :: fopenstat

! Executable code
     iu = initxdr(fname, fopenstat, lflag)
     IF (iu >= 0) THEN
          ierr = 0
     ELSE
          ierr = iu
     END IF

     RETURN

END SUBROUTINE io_open_fxdr
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
