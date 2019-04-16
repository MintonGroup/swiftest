!**********************************************************************************************************************************
!
!  Unit Name   : io_read_encounter
!  Unit Type   : function
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read close encounter data from input binary file
!
!  Input
!    Arguments : encounter_file : name of input binary file for encounters
!                out_type       : format of input binary file
!    Terminal  : none
!    File      : t              : time
!
!  Output
!    Arguments : t              : time
!                id1            : identifier of first planet / particle
!                id2            : identifier of second planet / particle
!                mass1          : mass of first planet / particle
!                mass2          : mass of second planet / particle
!                xh1            : heliocentric position of first planet / particle
!                xh2            : heliocentric position of second planet / particle
!                vh1            : heliocentric velocity of first planet / particle
!                vh2            : heliocentric velocity of second planet / particle
!    Terminal  : error message
!    File      : none
!
!  Invocation  : istat = io_read_encounter(t, id1, id2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file, out_type)
!
!  Notes       : Other than time t, there is no direct file input from this function
!
!                Function returns read error status (0 = OK, nonzero = ERROR)
!
!**********************************************************************************************************************************
FUNCTION io_read_encounter(t, id1, id2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_read_encounter
     IMPLICIT NONE

! Arguments
     INTEGER(I4B)                           :: io_read_encounter
     INTEGER(I4B), INTENT(OUT)              :: id1, id2
     REAL(DP), INTENT(OUT)                  :: t, mass1, mass2
     REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: xh1, xh2, vh1, vh2
     CHARACTER(*), INTENT(IN)               :: encounter_file, out_type

! Internals
     LOGICAL(LGT)            :: lxdr
     LOGICAL(LGT), SAVE      :: lfirst = .TRUE.
     INTEGER(I4B), PARAMETER :: LUN = 30
     INTEGER(I4B)            :: ierr
     INTEGER(I4B), SAVE      :: iu = LUN

! Executable code
     lxdr = ((out_type == XDR4_TYPE) .OR. (out_type == XDR8_TYPE))
     IF (lfirst) THEN
          IF (lxdr) THEN
               CALL io_open_fxdr(encounter_file, "R", .TRUE., iu, ierr)
          ELSE
               CALL io_open(iu, encounter_file, "OLD", "UNFORMATTED", ierr)
          END IF
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Unable to open binary encounter file"
               CALL util_exit(FAILURE)
          END IF
          lfirst = .FALSE.
     END IF
     IF (lxdr) THEN
          ierr = ixdrdouble(iu, t)
          io_read_encounter = ierr
          IF (ierr /= 0) THEN
               ierr = ixdrclose(iu)
               RETURN
          END IF
          ierr = io_read_line(iu, id1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), XDR8_TYPE, MASS = mass1)
          io_read_encounter = ierr
          IF (ierr /= 0) THEN
               ierr = ixdrclose(iu)
               RETURN
          END IF
          ierr = io_read_line(iu, id2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), vh2(3), XDR8_TYPE, MASS = mass2)
          io_read_encounter = ierr
          IF (ierr /= 0) THEN
               ierr = ixdrclose(iu)
               RETURN
          END IF
     ELSE
          READ(iu, IOSTAT = ierr) t
          io_read_encounter = ierr
          IF (ierr /= 0) THEN
               CLOSE(UNIT = iu, IOSTAT = ierr)
               RETURN
          END IF
          ierr = io_read_line(iu, id1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), REAL8_TYPE, MASS = mass1)
          io_read_encounter = ierr
          IF (ierr /= 0) THEN
               CLOSE(UNIT = iu, IOSTAT = ierr)
               RETURN
          END IF
          ierr = io_read_line(iu, id2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), vh2(3), REAL8_TYPE, MASS = mass2)
          io_read_encounter = ierr
          IF (ierr /= 0) THEN
               CLOSE(UNIT = iu, IOSTAT = ierr)
               RETURN
          END IF
     END IF

     RETURN

END FUNCTION io_read_encounter
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
