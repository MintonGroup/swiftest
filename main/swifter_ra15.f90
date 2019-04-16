!**********************************************************************************************************************************
!
!  Unit Name   : swifter_ra15
!  Unit Type   : program
!  Project     : Swifter
!  Package     : main
!  Language    : Fortran 90/95
!
!  Description : Driver program for the 15th-order RADAU integrator
!
!  Input
!    Arguments : none
!    Terminal  : parameter file name
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : status messages
!    File      : none
!
!  Invocation  : % swifter_ra15
!
!  Notes       : Reference: Everhart, E. 1985. Dynamics of Comets: Their Origin and Evolution, A. Carusi & B. Valsecchi (eds.),
!                           (D. Reidel Publishing Company), 185 - 202.
!
!**********************************************************************************************************************************
PROGRAM swifter_ra15

! Modules
     USE module_parameters
     USE module_swifter
     USE module_ra15
     USE module_random_access
     USE module_interfaces
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT)      :: lclose         ! Check for planet-test particle encounters
     LOGICAL(LGT)      :: lextra_force   ! Use user-supplied force routines
     LOGICAL(LGT)      :: lbig_discard   ! Dump planet data with discards
     LOGICAL(LGT)      :: lrhill_present ! Hill's sphere radius present
     INTEGER(I4B)      :: nplmax         ! Maximum number of planets
     INTEGER(I4B)      :: ntpmax         ! Maximum number of test particles
     INTEGER(I4B)      :: istep_out      ! Time steps between binary outputs
     INTEGER(I4B)      :: istep_dump     ! Time steps between dumps
     REAL(DP)          :: t0             ! Integration start time
     REAL(DP)          :: tstop          ! Integration stop time
     REAL(DP)          :: dt             ! Time step
     REAL(DP)          :: j2rp2          ! J2*R^2 term for central body
     REAL(DP)          :: j4rp4          ! J4*R^4 term for central body
     REAL(DP)          :: rmin           ! Minimum heliocentric radius for test particle
     REAL(DP)          :: rmax           ! Maximum heliocentric radius for test particle
     REAL(DP)          :: rmaxu          ! Maximum unbound heliocentric radius for test particle
     REAL(DP)          :: qmin           ! Minimum pericenter distance for test particle
     REAL(DP)          :: qmin_alo       ! Minimum semimajor axis for qmin
     REAL(DP)          :: qmin_ahi       ! Maximum semimajor axis for qmin
     CHARACTER(STRMAX) :: qmin_coord     ! Coordinate frame to use for qmin
     CHARACTER(STRMAX) :: encounter_file ! Name of output file for encounters
     CHARACTER(STRMAX) :: inplfile       ! Name of input file for planets
     CHARACTER(STRMAX) :: intpfile       ! Name of input file for test particles
     CHARACTER(STRMAX) :: in_type        ! Format of input data files
     CHARACTER(STRMAX) :: outfile        ! Name of output binary file
     CHARACTER(STRMAX) :: out_type       ! Binary format of output file
     CHARACTER(STRMAX) :: out_form       ! Data to write to output file
     CHARACTER(STRMAX) :: out_stat       ! Open status for output binary file

! Internals
     LOGICAL(LGT)                                     :: lfirst
     INTEGER(I4B)                                     :: npl, ntp, ntp0, nsp, iout, idump, iloop
     REAL(DP)                                         :: t, tfrac, tbase, eoffset
     CHARACTER(STRMAX)                                :: inparfile
     TYPE(swifter_pl), POINTER                        :: swifter_pl1P
     TYPE(swifter_tp), POINTER                        :: swifter_tp1P, swifter_tpd1P
     TYPE(ra15_pl), DIMENSION(:), ALLOCATABLE, TARGET :: ra15_plA
     TYPE(ra15_tp), DIMENSION(:), ALLOCATABLE, TARGET :: ra15_tpA
     TYPE(ra15_pl), POINTER                           :: ra15_pl1P
     TYPE(ra15_tp), POINTER                           :: ra15_tp1P, ra15_tpd1P

! Executable code
     CALL util_version
     WRITE(*, 100, ADVANCE = "NO") "Enter name of parameter data file: "
     READ(*, 100) inparfile
 100 FORMAT(A)
     inparfile = TRIM(ADJUSTL(inparfile))
     CALL io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile, out_type,      &
          out_form, out_stat, istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,          &
          encounter_file, lextra_force, lbig_discard, lrhill_present)
     CALL io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)
     ALLOCATE(ra15_plA(nplmax))
     CALL set_point(ra15_plA)
     IF (ntp > 0) THEN
          ALLOCATE(ra15_tpA(ntpmax))
          CALL set_point(ra15_tpA)
     END IF
     CALL ra15_setup(npl, ntp, ra15_plA, ra15_tpA, ra15_pl1P, ra15_tp1P, swifter_pl1P, swifter_tp1P)
     CALL io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, swifter_pl1P)
     CALL io_init_tp(intpfile, in_type, ntp, swifter_tp1P)
     CALL util_valid(npl, ntp, swifter_pl1P, swifter_tp1P)
     lfirst = .TRUE.
     ntp0 = ntp
     t = t0
     tbase = t0
     iloop = 0
     iout = istep_out
     idump = istep_dump
     nsp = 0
     eoffset = 0.0_DP
     NULLIFY(ra15_tpd1P)
     IF (istep_out > 0) CALL io_write_frame(t, npl, ntp, swifter_pl1P, swifter_tp1P, outfile, out_type, out_form, out_stat)
     WRITE(*, *) " *************** MAIN LOOP *************** "
     DO WHILE ((t < tstop) .AND. ((ntp0 == 0) .OR. (ntp > 0)))
          CALL ra15_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, ra15_pl1P, ra15_tp1P, j2rp2, j4rp4, dt)
          iloop = iloop + 1
          IF (iloop == LOOPMAX) THEN
               tbase = tbase + iloop*dt
               iloop = 0
          END IF
          t = tbase + iloop*dt
          CALL ra15_discard(t, npl, ntp, nsp, ra15_pl1P, ra15_tp1P, ra15_tpd1P, dt, rmin, rmax, rmaxu, qmin, qmin_coord,          &
               qmin_alo, qmin_ahi, lclose, lrhill_present)
          IF (nsp > 0) THEN
               swifter_tp1P => ra15_tp1P%swifter
               swifter_tpd1P => ra15_tpd1P%swifter
               CALL io_discard_write(t, npl, nsp, swifter_pl1P, swifter_tpd1P, DISCARD_FILE, lbig_discard)
               nsp = 0
               NULLIFY(ra15_tpd1P)
          END IF
          IF (istep_out > 0) THEN
               iout = iout - 1
               IF (iout == 0) THEN
                    CALL io_write_frame(t, npl, ntp, swifter_pl1P, swifter_tp1P, outfile, out_type, out_form, out_stat)
                    iout = istep_out
               END IF
          END IF
          IF (istep_dump > 0) THEN
               idump = idump - 1
               IF (idump == 0) THEN
                    tfrac = (t - t0)/(tstop - t0)
                    WRITE(*, 200) t, tfrac, npl, ntp
 200                FORMAT(" Time = ", ES12.5, "; fraction done = ", F5.3, "; Number of active pl, tp = ",I5, ", ", I5)
                    CALL io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form,        &
                         istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,               &
                         encounter_file, lextra_force, lbig_discard, lrhill_present)
                    CALL io_dump_pl(npl, swifter_pl1P, lclose, lrhill_present)
                    IF (ntp > 0) CALL io_dump_tp(ntp, swifter_tp1P)
                    idump = istep_dump
               END IF
          END IF
     END DO
     CALL io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form, istep_dump, j2rp2,    &
          j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, encounter_file, lextra_force, lbig_discard,     &
          lrhill_present)
     CALL io_dump_pl(npl, swifter_pl1P, lclose, lrhill_present)
     IF (ntp > 0) CALL io_dump_tp(ntp, swifter_tp1P)
     IF (ALLOCATED(ra15_plA)) DEALLOCATE(ra15_plA)
     IF (ALLOCATED(ra15_tpA)) DEALLOCATE(ra15_tpA)
     CALL util_exit(SUCCESS)

     STOP

END PROGRAM swifter_ra15
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
