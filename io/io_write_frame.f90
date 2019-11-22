!**********************************************************************************************************************************
!
!  Unit Name   : io_write_frame
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write a frame (header plus records for each planet and active test particle) to output binary file
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                ntp          : number of active test particles
!                swiftest_pl1P : pointer to head of swiftest planet structure linked-list
!                swiftest_tp1P : pointer to head of active swiftest test particle structure linked-list
!                outfile      : name of output binary file
!                out_type     : binary format of output file
!                out_form     : data to write to output file (elements / heliocentric coordinates / filtered elements)
!                out_stat     : open status for output binary file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL io_write_frame(t, npl, ntp, swiftest_pl1P, swiftest_tp1P, outfile, out_type, out_form, out_stat)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_write_frame.F
!
!                There is no direct file output from this subroutine
!
!**********************************************************************************************************************************
SUBROUTINE io_write_frame(t, npl, ntp, swiftest_plA, swiftest_tpA, outfile, out_type, out_form, out_stat)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_write_frame
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)         :: npl, ntp
     REAL(DP), INTENT(IN)             :: t
     CHARACTER(*), INTENT(IN)         :: outfile, out_type, out_form, out_stat
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
     TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA

! Internals
     LOGICAL(LGT)              :: lxdr
     LOGICAL(LGT), SAVE        :: lfirst = .TRUE.
     INTEGER(I4B), PARAMETER   :: LUN = 20
     INTEGER(I4B)              :: i, j, ierr
     INTEGER(I4B), SAVE        :: iu = LUN, iout_form = XV
     REAL(DP)                  :: a, e, inc, capom, omega, capm, mu
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp

! Executable code
     lxdr = ((out_type == XDR4_TYPE) .OR. (out_type == XDR8_TYPE))
     IF (lfirst) THEN
          IF (out_stat == "APPEND") THEN
               IF (lxdr) THEN
                    CALL io_open_fxdr(outfile, "A", .TRUE., iu, ierr)
               ELSE
                    CALL io_open(iu, outfile, out_stat, "UNFORMATTED", ierr)
               END IF
          ELSE IF (out_stat == "NEW") THEN
               IF (lxdr) THEN
                    CALL io_open_fxdr(outfile, "R", .TRUE., iu, ierr)
                    IF (ierr == 0) THEN
                         WRITE(*, *) "SWIFTEST Error:"
                         WRITE(*, *) "   Binary output file already exists"
                         CALL util_exit(FAILURE)
                    END IF
                    CALL io_open_fxdr(outfile, "W", .TRUE., iu, ierr)
               ELSE
                    CALL io_open(iu, outfile, out_stat, "UNFORMATTED", ierr)
                    IF (ierr /= 0) THEN
                         WRITE(*, *) "SWIFTEST Error:"
                         WRITE(*, *) "   Binary output file already exists"
                         CALL util_exit(FAILURE)
                    END IF
               END IF
          ELSE
               IF (lxdr) THEN
                    CALL io_open_fxdr(outfile, "W", .TRUE., iu, ierr)
               ELSE
                    CALL io_open(iu, outfile, "REPLACE", "UNFORMATTED", ierr)
               END IF
          END IF
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTEST Error:"
               WRITE(*, *) "   Unable to open binary output file"
               CALL util_exit(FAILURE)
          END IF
          SELECT CASE (out_form)
               CASE ("EL")
                    iout_form = EL
               CASE ("XV")
                    iout_form = XV
               CASE ("FILT")
                    iout_form = FILT
          END SELECT
          lfirst = .FALSE.
     ELSE
          IF (lxdr) THEN
               CALL io_open_fxdr(outfile, "A", .TRUE., iu, ierr)
          ELSE
               CALL io_open(iu, outfile, "APPEND", "UNFORMATTED", ierr)
          END IF
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTEST Error:"
               WRITE(*, *) "   Unable to open binary output file for append"
               CALL util_exit(FAILURE)
          END IF
     END IF
     CALL io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
     SELECT CASE (iout_form)
          CASE (EL)
               DO i = 2, npl
                    mu = swiftest_plA%mass(1) + swiftest_plA%mass(i)
                    j = swiftest_plA%name(i)
                    CALL orbel_xv2el(swiftest_plA%xh(:,i), swiftest_plA%vh(:,i), mu, a, e, inc, capom, omega, capm)
                    CALL io_write_line(iu, j, a, e, inc, capom, omega, capm, out_type, &
                     MASS = swiftest_plA%mass(i),RADIUS = swiftest_plA%radius(i))
               END DO
               mu = swiftest_plA%mass(1)
               DO i = 1, ntp
                    j = swiftest_tpA%name(i)
                    CALL orbel_xv2el(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i), mu, a, e, inc, capom, omega, capm)
                    CALL io_write_line(iu, j, a, e, inc, capom, omega, capm, out_type)
               END DO
          CASE (XV)
               DO i = 2, npl
                    xtmp(:) = swiftest_plA%xh(:,i)
                    vtmp(:) = swiftest_plA%vh(:,i)
                    j = swiftest_plA%name(i)
                    CALL io_write_line(iu, j, xtmp(1), xtmp(2), xtmp(3), vtmp(1), vtmp(2), vtmp(3), out_type,                     &
                         MASS = swiftest_plA%mass(i), RADIUS = swiftest_plA%radius(i))
               END DO
               DO i = 1, ntp
                    xtmp(:) = swiftest_tpA%xh(:,i)
                    vtmp(:) = swiftest_tpA%vh(:,i)
                    j = swiftest_tpA%name(i)
                    CALL io_write_line(iu, j, xtmp(1), xtmp(2), xtmp(3), vtmp(1), vtmp(2), vtmp(3), out_type)
               END DO
          CASE (FILT)
! DEK - add code here to handle the case for an OUT_FORM = FILT
     END SELECT
     IF (lxdr) THEN
          ierr = ixdrclose(iu)
     ELSE
          CLOSE(UNIT = iu, IOSTAT = ierr)
     END IF
     IF (ierr /= 0) THEN
          WRITE(*, *) "SWIFTEST Error:"
          WRITE(*, *) "   Unable to close binary output file"
          CALL util_exit(FAILURE)
     END IF

     RETURN

END SUBROUTINE io_write_frame
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
