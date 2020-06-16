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
!                id1            : identifier of first massive body / particle
!                id2            : identifier of second massive body / particle
!                mass1          : mass of first massive body / particle
!                mass2          : mass of second massive body / particle
!                xh1            : heliocentric position of first massive body / particle
!                xh2            : heliocentric position of second massive body / particle
!                vh1            : heliocentric velocity of first massive body / particle
!                vh2            : heliocentric velocity of second massive body / particle
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
FUNCTION io_read_encounter(t, name1, name2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file, out_type)

! Modules
     USE swiftest, EXCEPT_THIS_ONE => io_read_encounter
     IMPLICIT NONE

! Arguments
     INTEGER(I4B)                           :: io_read_encounter
     INTEGER(I4B), INTENT(OUT)              :: name1, name2
     REAL(DP), INTENT(OUT)                  :: t, mass1, mass2
     REAL(DP), DIMENSION(:), INTENT(OUT)    :: xh1, xh2, vh1, vh2
     CHARACTER(*), INTENT(IN)               :: encounter_file, out_type

! Internals
     LOGICAL(LGT)            :: lxdr
     LOGICAL(LGT), SAVE      :: lfirst = .TRUE.
     INTEGER(I4B), PARAMETER :: LUN = 30
     INTEGER(I4B)            :: ierr
     INTEGER(I4B), SAVE      :: iu = LUN

! Executable code
     IF (lfirst) THEN
         CALL io_open(iu, encounter_file, "OLD", "UNFORMATTED", ierr)
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Unable to open binary encounter file"
               CALL util_exit(FAILURE)
          END IF
          lfirst = .FALSE.
     END IF
    READ(iu, IOSTAT = ierr) t
    io_read_encounter = ierr
    IF (ierr /= 0) THEN
         CLOSE(UNIT = iu, IOSTAT = ierr)
         RETURN
    END IF
    ierr = io_read_line(iu, name1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), REAL8_TYPE, MASS = mass1)
    io_read_encounter = ierr
    IF (ierr /= 0) THEN
         CLOSE(UNIT = iu, IOSTAT = ierr)
         RETURN
    END IF
    ierr = io_read_line(iu, name2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), vh2(3), REAL8_TYPE, MASS = mass2)
    io_read_encounter = ierr
    IF (ierr /= 0) THEN
         CLOSE(UNIT = iu, IOSTAT = ierr)
         RETURN
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
