!**********************************************************************************************************************************
!
!  Unit Name   : io_init_param
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read in parameters for the integration
!
!  Input
!    Arguments : inparfile
!    Terminal  : none
!    File      : nplmax         : maximum allowed number of planets
!                ntpmax         : maximum allowed number of test particles
!                t0             : integration start time
!                tstop          : integration stop time
!                dt             : time step
!                inplfile       : name of input file for planets
!                intpfile       : name of input file for test particles
!                in_type        : format of input data files
!                istep_out      : number of time steps between binary outputs
!                outfile        : name of output binary file
!                out_type       : binary format of output file
!                out_form       : data to write to output file
!                out_stat       : open status for output binary file
!                j2rp2          : J2 * R**2 for the Sun
!                j4rp4          : J4 * R**4 for the Sun
!                lclose         : logical flag indicating whether to check for planet-test particle encounters
!                rmin           : minimum heliocentric radius for test particle
!                rmax           : maximum heliocentric radius for test particle
!                rmaxu          : maximum unbound heliocentric radius for test particle
!                qmin           : minimum pericenter distance for test particle
!                qmin_coord     : coordinate frame to use for qmin
!                qmin_alo       : minimum semimajor axis for qmin
!                qmin_ahi       : maximum semimajor axis for qmin
!                encounter_file : name of output file for encounters
!                lextra_force   : logical flag indicating whether to use user-supplied accelerations
!                lbig_discard   : logical flag indicating whether to dump planet data with discards
!                lrhill_present : logical flag indicating whether Hill's sphere radii are present in planet data
!                mtiny          : smallest self-gravitating mass (only used for SyMBA)
!  Output
!    Arguments : (same quantities listed above as input from file are passed back to the calling routine as output arguments)
!    Terminal  : status, error messages
!    File      : none
!
!  Invocation  : CALL io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile,
!                                   out_type, out_form, out_stat, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin,
!                                   qmin_coord, qmin_alo, qmin_ahi, encounter_file, lextra_force, lbig_discard, lrhill_present, mtiny)
!
!  Notes       : Adapted from Martin Duncan's Swift routine io_init_param.f
!
!**********************************************************************************************************************************
SUBROUTINE io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile, out_type,     &
     out_form, out_stat, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,               &
     encounter_file, lextra_force, lbig_discard, lrhill_present, mtiny)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => io_init_param
! OPENMP code added by D. Minton
     !$ use omp_lib
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(OUT) :: lclose, lextra_force, lbig_discard, lrhill_present
     INTEGER(I4B), INTENT(OUT) :: nplmax, ntpmax, istep_out
     REAL(DP), INTENT(OUT)     :: t0, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
     REAL(DP), INTENT(OUT),OPTIONAL :: mtiny
     CHARACTER(*), INTENT(IN)  :: inparfile
     CHARACTER(*), INTENT(OUT) :: qmin_coord, encounter_file, inplfile, intpfile, in_type, outfile, out_type, out_form, out_stat

! Internals
     LOGICAL(LGT)            :: t0_set, tstop_set, dt_set
     INTEGER(I4B), PARAMETER :: LUN = 7
     INTEGER(I4B)            :: ierr, ilength, ifirst, ilast
     CHARACTER(STRMAX)       :: line, token

! Executable code
     nplmax = -1
     ntpmax = -1
     t0_set = .FALSE.
     tstop_set = .FALSE.
     dt_set = .FALSE.
     t0 = 0.0_DP
     tstop = 0.0_DP
     dt = 0.0_DP
     inplfile = ""
     intpfile = ""
     in_type = "ASCII"
     istep_out = -1
     outfile = ""
     out_type = REAL8_TYPE
     out_form = "XV"
     out_stat = "NEW"
     j2rp2 = 0.0_DP
     j4rp4 = 0.0_DP
     lclose = .FALSE.
     rmin = -1.0_DP
     rmax = -1.0_DP
     rmaxu = -1.0_DP
     qmin = -1.0_DP
     qmin_coord = "HELIO"
     qmin_alo = -1.0_DP
     qmin_ahi = -1.0_DP
     encounter_file = ""
     lextra_force = .FALSE.
     lbig_discard = .FALSE.
     lrhill_present = .FALSE.
     mtiny = -1.0_DP
     WRITE(*, 100, ADVANCE = "NO") "Parameter data file is "
     WRITE(*, 100) inparfile
     WRITE(*, *) " "
 100 FORMAT(A)
     CALL io_open(LUN, inparfile, "OLD", "FORMATTED", ierr)
     IF (ierr /= 0) THEN
          WRITE(*, 100, ADVANCE = "NO") "Unable to open file "
          WRITE(*, 100) inparfile
          CALL util_exit(FAILURE)
     END IF
     DO
          READ(LUN, 100, IOSTAT = ierr, END = 1) line
          line = ADJUSTL(line)
          ilength = LEN_TRIM(line)
          IF ((ilength /= 0) .AND. (line(1:1) /= "!")) THEN
               ifirst = 1
               CALL io_get_token(line, ilength, ifirst, ilast, ierr)
               token = line(ifirst:ilast)
               CALL util_toupper(token)
               SELECT CASE (token)
                    CASE ("NPLMAX")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) nplmax
                    CASE ("NTPMAX")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) ntpmax
                    CASE ("T0")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) t0
                         t0_set = .TRUE.
                    CASE ("TSTOP")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) tstop
                         tstop_set = .TRUE.
                    CASE ("DT")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) dt
                         dt_set = .TRUE.
                    CASE ("PL_IN")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         inplfile = token
                    CASE ("TP_IN")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         intpfile = token
                    CASE ("IN_TYPE")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         in_type = token
                    CASE ("ISTEP_OUT")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) istep_out
                    CASE ("BIN_OUT")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         outfile = token
                    CASE ("OUT_TYPE")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         out_type = token
                    CASE ("OUT_FORM")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         out_form = token
                    CASE ("OUT_STAT")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         out_stat = token
                    CASE ("J2")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) j2rp2
                    CASE ("J4")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) j4rp4
                    CASE ("CHK_CLOSE")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         IF (token == "YES") lclose = .TRUE.
                    CASE ("CHK_RMIN")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) rmin
                    CASE ("CHK_RMAX")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) rmax
                    CASE ("CHK_EJECT")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) rmaxu
                    CASE ("CHK_QMIN")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) qmin
                    CASE ("CHK_QMIN_COORD")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         qmin_coord = token
                    CASE ("CHK_QMIN_RANGE")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) qmin_alo
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) qmin_ahi
                    CASE ("ENC_OUT")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         encounter_file = token
                    CASE ("EXTRA_FORCE")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         IF (token == "YES") lextra_force = .TRUE.
                    CASE ("BIG_DISCARD")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         IF (token == "YES") lbig_discard = .TRUE.
                    CASE ("RHILL_PRESENT")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         IF (token == "YES") lrhill_present = .TRUE.

                    ! Added by D. Minton
                    CASE ("FRAGMENTATION")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         CALL util_toupper(token)
                         IF (token == "YES") lfragmentation = .TRUE.
                    
                    CASE ("MU2GM")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) MU2GM
                    CASE ("TU2S")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) TU2S
                    CASE ("DU2CM")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         READ(token, *) DU2CM
                    !^^^^^^^^^^^^^^^^^^^^^^
                    CASE ("MTINY")
                         ifirst = ilast + 1
                         CALL io_get_token(line, ilength, ifirst, ilast, ierr)
                         token = line(ifirst:ilast)
                         IF (PRESENT(mtiny)) READ(token, *) mtiny
                    CASE DEFAULT
                         WRITE(*, 100, ADVANCE = "NO") "Unknown parameter -> "
                         WRITE(*, *) token
               END SELECT
          END IF
     END DO
  1  CLOSE(UNIT = LUN)
     WRITE(*, 100, ADVANCE = "NO") "NPLMAX         = "
     WRITE(*, *) nplmax
     WRITE(*, 100, ADVANCE = "NO") "NTPMAX         = "
     WRITE(*, *) ntpmax
     WRITE(*, 100, ADVANCE = "NO") "T0             = "
     WRITE(*, *) t0
     WRITE(*, 100, ADVANCE = "NO") "TSTOP          = "
     WRITE(*, *) tstop
     WRITE(*, 100, ADVANCE = "NO") "DT             = "
     WRITE(*, *) dt
     WRITE(*, 100, ADVANCE = "NO") "PL_IN          = "
     ilength = LEN_TRIM(inplfile)
     WRITE(*, *) inplfile(1:ilength)
     WRITE(*, 100, ADVANCE = "NO") "TP_IN          = "
     ilength = LEN_TRIM(intpfile)
     WRITE(*, *) intpfile(1:ilength)
     WRITE(*, 100, ADVANCE = "NO") "IN_TYPE        = "
     ilength = LEN_TRIM(in_type)
     WRITE(*, *) in_type(1:ilength)
     WRITE(*, 100, ADVANCE = "NO") "ISTEP_OUT      = "
     WRITE(*, *) istep_out
     WRITE(*, 100, ADVANCE = "NO") "BIN_OUT        = "
     ilength = LEN_TRIM(outfile)
     WRITE(*, *) outfile(1:ilength)
     WRITE(*, 100, ADVANCE = "NO") "OUT_TYPE       = "
     ilength = LEN_TRIM(out_type)
     WRITE(*, *) out_type(1:ilength)
     WRITE(*, 100, ADVANCE = "NO") "OUT_FORM       = "
     ilength = LEN_TRIM(out_form)
     WRITE(*, *) out_form(1:ilength)
     WRITE(*, 100, ADVANCE = "NO") "OUT_STAT       = "
     ilength = LEN_TRIM(out_stat)
     WRITE(*, *) out_stat(1:ilength)
     WRITE(*, 100, ADVANCE = "NO") "J2             = "
     WRITE(*, *) j2rp2
     WRITE(*, 100, ADVANCE = "NO") "J4             = "
     WRITE(*, *) j4rp4
     WRITE(*, 100, ADVANCE = "NO") "CHK_CLOSE      = "
     WRITE(*, *) lclose
     WRITE(*, 100, ADVANCE = "NO") "CHK_RMIN       = "
     WRITE(*, *) rmin
     WRITE(*, 100, ADVANCE = "NO") "CHK_RMAX       = "
     WRITE(*, *) rmax
     WRITE(*, 100, ADVANCE = "NO") "CHK_EJECT      = "
     WRITE(*, *) rmaxu
     WRITE(*, 100, ADVANCE = "NO") "CHK_QMIN       = "
     WRITE(*, *) qmin
     WRITE(*, 100, ADVANCE = "NO") "CHK_QMIN_COORD = "
     ilength = LEN_TRIM(qmin_coord)
     WRITE(*, *) qmin_coord(1:ilength)
     WRITE(*, 100, ADVANCE = "NO") "CHK_QMIN_RANGE = "
     WRITE(*, *) qmin_alo, qmin_ahi
     WRITE(*, 100, ADVANCE = "NO") "ENC_OUT        = "
     ilength = LEN_TRIM(encounter_file)
     WRITE(*, *) encounter_file(1:ilength)
     WRITE(*, 100, ADVANCE = "NO") "EXTRA_FORCE    = "
     WRITE(*, *) lextra_force
     WRITE(*, 100, ADVANCE = "NO") "BIG_DISCARD    = "
     WRITE(*, *) lbig_discard
     WRITE(*, 100, ADVANCE = "NO") "RHILL_PRESENT  = "
     WRITE(*, *) lrhill_present
     WRITE(*, *) " "
     ierr = 0
     IF ((.NOT. t0_set) .OR. (.NOT. tstop_set) .OR. (.NOT. dt_set)) ierr = -1
     IF (dt <= 0.0_DP) ierr = -1
     IF (inplfile == "") ierr = -1
     IF ((in_type /= REAL8_TYPE) .AND. (in_type /= "ASCII")) ierr = -1
     IF ((istep_out > 0) .AND. (outfile == "")) ierr = -1
     IF (outfile /= "") THEN
          IF ((out_type /= REAL4_TYPE) .AND. (out_type /= REAL8_TYPE))  ierr = -1
          IF ((out_form /= "EL") .AND. (out_form /= "XV") .AND. (out_form /= "FILT")) ierr = -1
          IF ((out_stat /= "NEW") .AND. (out_stat /= "UNKNOWN") .AND. (out_stat /= "APPEND")) ierr = -1
     END IF
     IF ((j2rp2 == 0.0_DP) .AND. (j4rp4 /= 0.0_DP)) ierr = -1
     IF (qmin > 0.0_DP) THEN
          IF ((qmin_coord /= "HELIO") .AND. (qmin_coord /= "BARY")) ierr = -1
          IF ((qmin_alo <= 0.0_DP) .OR. (qmin_ahi <= 0.0_DP)) ierr = -1
     END IF

     ! Added by D. Minton
     ! The fragmentation model requires the user to set the unit system explicitly.
     WRITE(*, 100, ADVANCE = "NO") "FRAGMENTATION  = "
     write(*,*) lfragmentation
     WRITE(*, 100, ADVANCE = "NO") "MU2GM          = "
     WRITE(*, *) MU2GM
     WRITE(*, 100, ADVANCE = "NO") "TU2S           = "
     WRITE(*, *) TU2S 
     WRITE(*, 100, ADVANCE = "NO") "DU2CM          = "
     WRITE(*, *) DU2CM
     IF ((MU2GM < 0.0_DP) .OR. (TU2S < 0.0_DP) .OR. (DU2CM < 0.0_DP)) ierr = -1
     GU       = GC / (DU2CM**3 / (MU2GM * TU2S**2))

     !Added mtiny to the argument list rather than from the terminal
     IF (PRESENT(mtiny)) THEN
         IF (mtiny < 0.0_DP) THEN
             WRITE(*,*) "MTINY not set or invalid value"
             ierr = -1
         ELSE
            WRITE(*, 100, ADVANCE = "NO") "MTINY          = "
            WRITE(*, *) mtiny    
         END IF
     END IF
            

     IF (ierr < 0) THEN
          WRITE(*, 100) "Input parameter(s) failed check"
          CALL util_exit(FAILURE)
     END IF

     ! OPENMP code added by D. Minton
     ! Define the maximum number of threads
     nthreads = 1                        ! In the *serial* case
     !$ nthreads = omp_get_max_threads() ! In the *parallel* case
     !$ write(*,'(a)')      ' OpenMP parameters:'
     !$ write(*,'(a)')      ' ------------------'
     !$ write(*,'(a,i3,/)') ' Number of threads  = ', nthreads     

     RETURN

END SUBROUTINE io_init_param
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
