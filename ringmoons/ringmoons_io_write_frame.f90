!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_io_write_frame
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write a frame (header plus records for each planet and active test particle) to output binary file
!
!  Input
!    Arguments : t            : time
!                out_stat     : open status for output binary file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL ringmoons_io_write_frame(t, ring)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE ringmoons_io_write_frame(t, ring, ring_outfile, out_stat)

! Modules
     USE module_parameters
     USE module_ringmoons
     USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_io_write_frame
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(IN)            :: t
     type(ringmoons_ring),intent(in) :: ring
     CHARACTER(*), INTENT(IN)  :: ring_outfile, out_stat

! Internals
     LOGICAL(LGT), SAVE        :: lfirst = .TRUE.
     INTEGER(I4B), PARAMETER   :: LUN = 88
     INTEGER(I4B)              :: ierr 
     INTEGER(I4B), SAVE        :: iu = LUN

! Executable code
     IF (lfirst) THEN
          IF (out_stat == "APPEND") THEN
              CALL io_open(iu, ring_outfile, out_stat, "UNFORMATTED", ierr)
          ELSE IF (out_stat == "NEW") THEN
              CALL io_open(iu, ring_outfile, out_stat, "UNFORMATTED", ierr)
              IF (ierr /= 0) THEN
                 WRITE(*, *) "RINGMOONS Error:"
                 WRITE(*, *) "   Binary output file already exists"
                 CALL util_exit(FAILURE)
              END IF
          ELSE
              CALL io_open(iu, ring_outfile, "REPLACE", "UNFORMATTED", ierr)
          END IF
          IF (ierr /= 0) THEN
              WRITE(*, *) "RINGMOONS Error:"
              WRITE(*, *) "   Unable to open binary output file"
              CALL util_exit(FAILURE)
          END IF
          WRITE(iu, IOSTAT = ierr) ring%N
          WRITE(iu, IOSTAT = ierr) ring%r_I
          WRITE(iu, IOSTAT = ierr) ring%r_F
          WRITE(iu, IOSTAT = ierr) ring%r_pdisk
          WRITE(iu, IOSTAT = ierr) ring%Gm_pdisk
          WRITE(iu, IOSTAT = ierr) ring%r
          lfirst = .FALSE.
     ELSE
          CALL io_open(iu, ring_outfile, "APPEND", "UNFORMATTED", ierr)
          IF (ierr /= 0) THEN
               WRITE(*, *) "RINGMOONS Error:"
               WRITE(*, *) "   Unable to open binary output file for append"
               CALL util_exit(FAILURE)
          END IF
     END IF
     WRITE(iu, IOSTAT = ierr) t
     WRITE(iu, IOSTAT = ierr) ring%Gsigma
     CLOSE(UNIT = iu, IOSTAT = ierr)
     IF (ierr /= 0) THEN
          WRITE(*, *) "RINGMOONS Error:"
          WRITE(*, *) "   Unable to close binary output file"
          CALL util_exit(FAILURE)
     END IF

     RETURN

END SUBROUTINE ringmoons_io_write_frame
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
