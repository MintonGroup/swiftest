!**********************************************************************************************************************************
!
!  Unit Name   : Unit_test
!  Unit Type   : program
!  Project     : Swifter
!  Package     : main
!  Language    : Fortran 90/95
!
!  Description : Unit test for symba_add function
!
!  Input
!    Arguments : none
!    Terminal  : parameter file name
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!                status messages
!    File      : none
!
!  Invocation  : % tool_timetest
!
!  Notes       : 
!
!**********************************************************************************************************************************
PROGRAM unit_test

! Modules
     USE module_parameters
     USE module_swifter
     USE module_symba
     USE module_random_access
     USE module_interfaces
     USE module_timetest
     !Added by D. Minton
     !$ USE omp_lib
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
     LOGICAL(LGT)                                      :: lfirst
     INTEGER(I4B)                                      :: npl, ntp, ntp0, nsppl, nsptp, iout, idump, iloop
     INTEGER(I4B)                                      :: nplplenc, npltpenc, nmergeadd, nmergesub,i,nn
     REAL(DP)                                          :: t, tfrac, tbase, mtiny, ke, pe, te, eoffset
     REAL(DP), DIMENSION(NDIM)                         :: htot
     CHARACTER(STRMAX)                                 :: inparfile
     TYPE(swifter_pl), POINTER                         :: swifter_pl1P
     TYPE(swifter_tp), POINTER                         :: swifter_tp1P
     TYPE(symba_pl), DIMENSION(:), ALLOCATABLE, TARGET :: symba_plA
     TYPE(symba_tp), DIMENSION(:), ALLOCATABLE, TARGET :: symba_tpA
     TYPE(symba_pl), POINTER                           :: symba_pl1P, symba_pld1P
     TYPE(symba_tp), POINTER                           :: symba_tp1P, symba_tpd1P
     TYPE(symba_pl), POINTER                           :: symba_pliP
     TYPE(symba_plplenc), DIMENSION(NENMAX)            :: plplenc_list
     TYPE(symba_pltpenc), DIMENSION(NENMAX)            :: pltpenc_list
     TYPE(symba_merger), DIMENSION(:), ALLOCATABLE     :: mergeadd_list, mergesub_list
     REAL(DP)                                          :: tstart, tend
     REAL(DP),DIMENSION(:),ALLOCATABLE                 :: toriginal,tstatic_arrays,ttype_arrays

     INTEGER(I4B),PARAMETER :: NPLTOT = 12105
     REAL(DP),DIMENSION(NDIM,NPLTOT) :: XX,AA
     REAL(DP),DIMENSION(NPLTOT) :: MM

     CHARACTER(*), PARAMETER :: TEST_FILE = "test_file.out"
     INTEGER(I4B), INTENT(OUT) :: ierr
! Executable code
     CALL util_version
     ! OpenMP code added by D. Minton
     ! Define the maximum number of threads
     nthreads = 1                        ! In the *serial* case
     !$ write(*,*) 'Dynamic thread allocation: ',OMP_get_dynamic()
     !$ nthreads = OMP_get_max_threads() ! In the *parallel* case
     !$ write(*,'(a)')      ' OpenMP parameters:'
     !$ write(*,'(a)')      ' ------------------'
     !$ write(*,'(a,i3,/)') ' Number of threads  = ', nthreads 
     WRITE(*, 100, ADVANCE = "NO") "Enter name of parameter data file: "
     READ(*, 100) inparfile
 100 FORMAT(A)
     inparfile = TRIM(ADJUSTL(inparfile))
     ! Read in the param.in file and get simulation parameters
     CALL io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile, out_type,      &
          out_form, out_stat, istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,          &
          encounter_file, lextra_force, lbig_discard, lrhill_present)
     IF (.NOT. lrhill_present) THEN
          WRITE(*, *) "SWIFTER Error:"
          WRITE(*, *) "   Integrator SyMBA requires planet Hill sphere radii on input"
          CALL util_exit(FAILURE)
     END IF
     ! Read in the total number of bodies from the input files
     CALL io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)

     ! Create arrays of data structures big enough to store the number of bodies we are adding
     ALLOCATE(symba_plA(nplmax), mergeadd_list(nplmax), mergesub_list(nplmax))

     ! Linked-list jiggery pokery that is probably very important
     CALL set_point(symba_plA)

     IF (ntp > 0) THEN
          ALLOCATE(symba_tpA(ntpmax))
          CALL set_point(symba_tpA)
     END IF
     
     ! Converts the array that we allocated into a linked list structure
     CALL symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1P, symba_tp1P, swifter_pl1P, swifter_tp1P)

     ! Reads in initial conditions of all massive bodies from input file and fills the linked list
     CALL io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, swifter_pl1P)
     WRITE(*, 100, ADVANCE = "NO") "Enter the smallest mass to self-gravitate: "
     READ(*, *) mtiny

     ! Reorder linked list by mass 
     CALL symba_reorder_pl(npl, symba_pl1P)
     CALL io_init_tp(intpfile, in_type, ntp, swifter_tp1P)
     CALL util_valid(npl, ntp, swifter_pl1P, swifter_tp1P)
     lfirst = .TRUE.
     ntp0 = ntp
     t = t0
     tbase = t0
     iloop = 0
     iout = istep_out
     idump = istep_dump
     nmergeadd = 0
     nmergesub = 0
     nsppl = 0
     nsptp = 0
     eoffset = 0.0_DP
     NULLIFY(symba_pld1P, symba_tpd1P)
     IF (istep_out > 0) CALL io_write_frame(t, npl, ntp, swifter_pl1P, swifter_tp1P, outfile, out_type, out_form, out_stat)
     WRITE(*, *) " *************** MAIN LOOP *************** "
     DO WHILE ((t < tstop) .AND. ((ntp0 == 0) .OR. (ntp > 0)))
          CALL symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_pl1P, symba_tp1P, j2rp2, j4rp4, dt,    &
               nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset,       &
               mtiny, encounter_file, out_type)
          iloop = iloop + 1
          IF (iloop == LOOPMAX) THEN
               tbase = tbase + iloop*dt
               iloop = 0
          END IF
          t = tbase + iloop*dt
          ! Take the merger info and create fragments
          ! CALL some subroutine that returns the number of fragments and an array of new bodies (swifter_pl type)
          IF (lfragmentation) THEN
               ! CALL symba_fragmentation(t, npl, nplmax, ntp, ntpmax, symba_pl1P, nplplenc, plplenc_list)
               ! update nplmax to add in the new number of bodies
               ! add new bodies into the current body linked list as in CALL symba_setup
               ! reorder bodies (if that is not already going to happen..check the discard subroutines
               CALL symba_add(npl, mergeadd_list, nmergeadd, symba_pl1P, swifter_pl1P, mtiny)
                symba_pliP => symba_pl1P
                DO i = 2, npl
                    symba_pliP => symba_pliP%nextP
                    swifter_pliP => symba_pliP%helio%swifter
                    id1 = swifter_pliP%id
                    m1 = swifter_pliP%mass
                    rad1 = swifter_pliP%radius
                    x1(:) = swifter_pliP%xh(:)
                    v1(:) = swifter_pliP%vh(:) 
                    CALL io_open(40,TEST_FILE, "APPEND", "FORMATTED", ierr)
                    CALL io_write_line(iu, id1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), REAL8_TYPE, MASS = mass, RADIUS = radius)
          END IF

          !
          CALL symba_discard_merge_pl(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
          CALL symba_discard_pl(t, npl, nplmax, nsppl, symba_pl1P, symba_pld1P, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,    &
               qmin_ahi, j2rp2, j4rp4, eoffset)
          CALL symba_discard_tp(t, npl, ntp, nsptp, symba_pl1P, symba_tp1P, symba_tpd1P, dt, rmin, rmax, rmaxu, qmin, qmin_coord, &
               qmin_alo, qmin_ahi, lclose, lrhill_present)
          IF ((nsppl > 0) .OR. (nsptp > 0)) THEN
               swifter_tp1P => symba_tp1P%helio%swifter
               CALL io_discard_write_symba(t, mtiny, npl, nsppl, nsptp, nmergeadd, nmergesub, symba_pl1P, symba_pld1P,            &
                    symba_tpd1P, mergeadd_list, mergesub_list, DISCARD_FILE, lbig_discard)
               nmergeadd = 0
               nmergesub = 0
               nsppl = 0
               nsptp = 0
               NULLIFY(symba_pld1P, symba_tpd1P)
          END IF
          IF (istep_out > 0) THEN
               iout = iout - 1
               IF (iout == 0) THEN
                    CALL io_write_frame(t, npl, ntp, swifter_pl1P, swifter_tp1P, outfile, out_type, out_form, out_stat)
                    iout = istep_out
               END IF
          END IF
     END DO
     IF (ALLOCATED(symba_plA)) DEALLOCATE(symba_plA)
     IF (ALLOCATED(mergeadd_list)) DEALLOCATE(mergeadd_list)
     IF (ALLOCATED(mergesub_list)) DEALLOCATE(mergesub_list)
     IF (ALLOCATED(symba_tpA)) DEALLOCATE(symba_tpA)
     CALL util_exit(SUCCESS)

     STOP

END PROGRAM unit_test
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
