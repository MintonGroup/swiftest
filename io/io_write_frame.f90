!**********************************************************************************************************************************
!
!  Unit Name   : io_write_frame
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write a frame (header plus records for each planet and active test particle) to output binary file
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                ntp          : number of active test particles
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
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
!  Invocation  : CALL io_write_frame(t, npl, ntp, swifter_pl1P, swifter_tp1P, outfile, out_type, out_form, out_stat)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_write_frame.F
!
!                There is no direct file output from this subroutine
!
!**********************************************************************************************************************************
SUBROUTINE io_write_frame(t, npl, ntp, swifter_pl1P, swifter_tp1P, outfile, out_type, out_form, out_stat)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_write_frame
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl, ntp
     REAL(DP), INTENT(IN)      :: t
     CHARACTER(*), INTENT(IN)  :: outfile, out_type, out_form, out_stat
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tp1P

! Internals
     LOGICAL(LGT)              :: lxdr
     LOGICAL(LGT), SAVE        :: lfirst = .TRUE.
     INTEGER(I4B), PARAMETER   :: LUN = 20
     INTEGER(I4B)              :: i, j, ierr
     INTEGER(I4B), SAVE        :: iu = LUN, iout_form = XV
     REAL(DP)                  :: a, e, inc, capom, omega, capm, mu
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP

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
                         WRITE(*, *) "SWIFTER Error:"
                         WRITE(*, *) "   Binary output file already exists"
                         CALL util_exit(FAILURE)
                    END IF
                    CALL io_open_fxdr(outfile, "W", .TRUE., iu, ierr)
               ELSE
                    CALL io_open(iu, outfile, out_stat, "UNFORMATTED", ierr)
                    IF (ierr /= 0) THEN
                         WRITE(*, *) "SWIFTER Error:"
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
               WRITE(*, *) "SWIFTER Error:"
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
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Unable to open binary output file for append"
               CALL util_exit(FAILURE)
          END IF
     END IF
     CALL io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
     SELECT CASE (iout_form)
          CASE (EL)
               swifter_plP => swifter_pl1P
               DO i = 2, npl
                    swifter_plP => swifter_plP%nextP
                    mu = swifter_pl1P%mass + swifter_plP%mass
                    j = swifter_plP%id
                    CALL orbel_xv2el(swifter_plP%xh(:), swifter_plP%vh(:), mu, a, e, inc, capom, omega, capm)
                    CALL io_write_line(iu, j, a, e, inc, capom, omega, capm, out_type, MASS = swifter_plP%mass,                   &
                         RADIUS = swifter_plP%radius)
               END DO
               mu = swifter_pl1P%mass
               swifter_tpP => swifter_tp1P
               DO i = 1, ntp
                    j = swifter_tpP%id
                    CALL orbel_xv2el(swifter_tpP%xh(:), swifter_tpP%vh(:), mu, a, e, inc, capom, omega, capm)
                    CALL io_write_line(iu, j, a, e, inc, capom, omega, capm, out_type)
                    swifter_tpP => swifter_tpP%nextP
               END DO
          CASE (XV)
               swifter_plP => swifter_pl1P
               DO i = 2, npl
                    swifter_plP => swifter_plP%nextP
                    xtmp(:) = swifter_plP%xh(:)
                    vtmp(:) = swifter_plP%vh(:)
                    j = swifter_plP%id
                    CALL io_write_line(iu, j, xtmp(1), xtmp(2), xtmp(3), vtmp(1), vtmp(2), vtmp(3), out_type,                     &
                         MASS = swifter_plP%mass, RADIUS = swifter_plP%radius)
               END DO
               swifter_tpP => swifter_tp1P
               DO i = 1, ntp
                    xtmp(:) = swifter_tpP%xh(:)
                    vtmp(:) = swifter_tpP%vh(:)
                    j = swifter_tpP%id
                    CALL io_write_line(iu, j, xtmp(1), xtmp(2), xtmp(3), vtmp(1), vtmp(2), vtmp(3), out_type)
                    swifter_tpP => swifter_tpP%nextP
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
          WRITE(*, *) "SWIFTER Error:"
          WRITE(*, *) "   Unable to close binary output file"
          CALL util_exit(FAILURE)
     END IF

     RETURN

END SUBROUTINE io_write_frame
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
