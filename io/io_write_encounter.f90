!**********************************************************************************************************************************
!
!  Unit Name   : io_write_encounter
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write close encounter data to output binary file
!
!  Input
!    Arguments : t              : time
!                id1            : identifier of first planet / particle
!                id2            : identifier of second planet / particle
!                mass1          : mass of first planet / particle
!                mass2          : mass of second planet / particle
!                radius1        : radius of first planet / particle
!                radius2        : radius of second planet / particle
!                xh1            : heliocentric position of first planet / particle
!                xh2            : heliocentric position of second planet / particle
!                vh1            : heliocentric velocity of first planet / particle
!                vh2            : heliocentric velocity of second planet / particle
!                encounter_file : name of output binary file for encounters
!                out_type       : format of output binary file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL io_write_encounter(t, id1, id2, mass1, mass2, radius1, radius2, xh1, xh2, vh1, vh2, encounter_file, out_type)
!
!  Notes       : There is no direct file output from this subroutine
!
!**********************************************************************************************************************************
SUBROUTINE io_write_encounter(t, id1, id2, mass1, mass2, radius1, radius2, xh1, xh2, vh1, vh2, encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_write_encounter
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)              :: id1, id2
     REAL(DP), INTENT(IN)                  :: t, mass1, mass2, radius1, radius2
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: xh1, xh2, vh1, vh2
     CHARACTER(*), INTENT(IN)              :: encounter_file, out_type

! Internals
     LOGICAL(LGT)            :: lxdr
     LOGICAL(LGT), SAVE      :: lfirst = .TRUE.
     INTEGER(I4B), PARAMETER :: LUN = 30
     INTEGER(I4B)            :: ierr
     INTEGER(I4B), SAVE      :: iu = LUN

! Executable code
     lxdr = ((out_type == XDR4_TYPE) .OR. (out_type == XDR8_TYPE))
     IF (lxdr) THEN
          CALL io_open_fxdr(encounter_file, "A", .TRUE., iu, ierr)
          IF ((ierr /= 0) .AND. lfirst) THEN
               CALL io_open_fxdr(encounter_file, "W", .TRUE., iu, ierr)
               lfirst = .FALSE.
          END IF
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Unable to open binary encounter file"
               CALL util_exit(FAILURE)
          END IF
          ierr = ixdrdouble(iu, t)
          IF (ierr < 0) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Unable to write binary file record"
               CALL util_exit(FAILURE)
          END IF
          CALL io_write_line(iu, id1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), XDR8_TYPE, MASS = mass1, RADIUS = radius1)
          CALL io_write_line(iu, id2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), vh2(3), XDR8_TYPE, MASS = mass2, RADIUS = radius2)
          ierr = ixdrclose(iu)
     ELSE
          CALL io_open(iu, encounter_file, "APPEND", "UNFORMATTED", ierr)
          IF ((ierr /= 0) .AND. lfirst) THEN
               CALL io_open(iu, encounter_file, "NEW", "UNFORMATTED", ierr)
               lfirst = .FALSE.
          END IF
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Unable to open binary encounter file"
               CALL util_exit(FAILURE)
          END IF
          WRITE(iu, IOSTAT = ierr) t
          IF (ierr < 0) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Unable to write binary file record"
               CALL util_exit(FAILURE)
          END IF
          CALL io_write_line(iu, id1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), REAL8_TYPE, MASS = mass1, RADIUS = radius1)
          CALL io_write_line(iu, id2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), vh2(3), REAL8_TYPE, MASS = mass2, RADIUS = radius2)
          CLOSE(UNIT = iu, IOSTAT = ierr)
     END IF
     IF (ierr /= 0) THEN
          WRITE(*, *) "SWIFTER Error:"
          WRITE(*, *) "   Unable to close binary encounter file"
          CALL util_exit(FAILURE)
     END IF

     RETURN

END SUBROUTINE io_write_encounter
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
