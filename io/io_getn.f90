!**********************************************************************************************************************************
!
!  Unit Name   : io_getn
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read the number of planets and test particles from respective input files
!
!  Input
!    Arguments : inplfile : name of input file for planets
!                intpfile : name of input file for test particles
!                in_type  : format of input data files
!                nplmax   : maximum allowed number of planets
!                ntpmax   : maximum allowed number of test particles
!    Terminal  : none
!    File      : npl      : number of planets
!                ntp      : number of active test particles
!
!  Output
!    Arguments : npl      : number of planets
!                ntp      : number of active test particles
!                nplmax   : maximum allowed number of planets
!                ntpmax   : maximum allowed number of test particles
!    Terminal  : error, warning messages
!    File      : none
!
!  Invocation  : CALL io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)
!
!  Notes       : nplmax (ntpmax) is reset to npl (ntp) if the latter exceeds the former
!
!**********************************************************************************************************************************
SUBROUTINE io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)

! Modules
     USE module_parameters
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_getn
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT) :: nplmax, ntpmax
     INTEGER(I4B), INTENT(OUT)   :: npl, ntp
     CHARACTER(*), INTENT(IN)    :: inplfile, intpfile, in_type

! Internals
     INTEGER(I4B), PARAMETER :: LUN = 7
     INTEGER(I4B)            :: iu, ierr

! Executable code
     npl = 0
     IF (in_type == "ASCII") THEN
          CALL io_open(LUN, inplfile, "OLD", "FORMATTED", ierr)
          READ(LUN, *) npl
          CLOSE(UNIT = LUN)
     ELSE
          CALL io_open_fxdr(inplfile, "R", .TRUE., iu, ierr)
          ierr = ixdrint(iu, npl)
          ierr = ixdrclose(iu)
     END IF
     IF (npl < 1) THEN
          WRITE(*, *) "Error: the number of planets, ", npl, ","
          WRITE(*, *) "       must be at least 1"
          CALL util_exit(FAILURE)
     ELSE IF ((npl > nplmax) .AND. .NOT. (nplmax < 0)) THEN
          WRITE(*, *) "Warning: the number of planets, ", npl, ","
          WRITE(*, *) "         exceeds the specified maximum number of planets, ", nplmax
          WRITE(*, *) "         ...resetting nplmax to ", npl
          nplmax = npl
     ELSE IF (nplmax < 0) THEN
          nplmax = npl
     END IF
     ntp = 0
     IF (intpfile /= "") THEN
          IF (in_type == "ASCII") THEN
               CALL io_open(LUN, intpfile, "OLD", "FORMATTED", ierr)
               READ(LUN, *) ntp
               CLOSE(UNIT = LUN)
          ELSE
               CALL io_open_fxdr(intpfile, "R", .TRUE., iu, ierr)
               ierr = ixdrint(iu, ntp)
               ierr = ixdrclose(iu)
          END IF
          IF (ntp < 0) THEN
               WRITE(*, *) "Error: the number of test particles, ", ntp, ","
               WRITE(*, *) "       must be at least 0"
               CALL util_exit(FAILURE)
          ELSE IF ((ntp > ntpmax) .AND. .NOT. (ntpmax < 0)) THEN
               WRITE(*, *) "Warning: the number of test particles, ", ntp, ","
               WRITE(*, *) "         exceeds the specified maximum number of test particles, ", ntpmax
               WRITE(*, *) "         ...resetting ntpmax to ", ntp
               ntpmax = ntp
          ELSE IF (ntpmax < 0) THEN
               ntpmax = ntp
          END IF
     END IF

     RETURN

END SUBROUTINE io_getn
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
