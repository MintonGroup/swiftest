!**********************************************************************************************************************************
!
!  Unit Name   : swiftest_symba_omp
!  Unit Type   : program
!  Project     : Swiftest
!  Package     : main
!  Language    : Fortran 90/95
!
!  Description : Driver program for the Symplectic Massive Body Algorithm with
!  OpenMP parallelization
!
!  Input
!    Arguments : none
!    Terminal  : parameter file name
!                mtiny               : smallest self-gravitating mass
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!                status messages
!    File      : none
!
!  Invocation  : % swiftest_symba
!
!  Notes       : Reference: Duncan, M. J., Levison, H. F. & Lee, M. H. 1998. Astron. J., 116, 2067.
!
!                Adapted from Hal Levison and Martin Duncan's Swift program swift_symba5.f
!                OpenMP parallelization by David Minton
!
!**********************************************************************************************************************************
PROGRAM swiftest_symba_omp

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_symba
     USE module_interfaces
     USE module_swiftestalloc
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
     LOGICAL(LGT)                                               :: lfirst
     INTEGER(I4B)                                               :: npl, ntp, ntp0, nsppl, nsptp, iout, idump, iloop
     INTEGER(I4B)                                               :: nplplenc, npltpenc, nmergeadd, nmergesub
     REAL(DP)                                                   :: t, tfrac, tbase, mtiny, ke, pe, te, eoffset
     REAL(DP), DIMENSION(NDIM)                                  :: htot
     CHARACTER(STRMAX)                                          :: inparfile
     TYPE(symba_pl)                                                   :: symba_plA
     TYPE(symba_tp)                                                   :: symba_tpA
     TYPE(swiftest_pl)                                                :: swiftest_plA
     TYPE(swiftest_tp)                                                :: swiftest_tpA
     TYPE(helio_pl)                                                   :: helio_plA
     TYPE(helio_tp)                                                   :: helio_tpA
     TYPE(symba_plplenc)                                              :: plplenc_list
     TYPE(symba_pltpenc)                                              :: pltpenc_list
     TYPE(symba_merger)                                               :: mergeadd_list, mergesub_list
     REAL(DP), DIMENSION(:,:), allocatable                            :: discard_plA
     REAL(DP), DIMENSION(:,:), allocatable                       :: discard_tpA
     INTEGER(I4B), DIMENSION(:,:), allocatable                   :: discard_plA_id_status
     INTEGER(I4B), DIMENSION(:,:), allocatable                   :: discard_tpA_id_status
     INTEGER(I4B), PARAMETER                                     :: egyiu = 72

! Executable code
     CALL util_version
     nthreads = 1                        ! In the *serial* case
     WRITE(*, 100, ADVANCE = "NO") "Enter name of parameter data file: "
     READ(*, 100) inparfile
 100 FORMAT(A)
     inparfile = TRIM(ADJUSTL(inparfile))
     ! Read in the param.in file and get simulation parameters
     CALL io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile, out_type,      &
          out_form, out_stat, istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,          &
          encounter_file, lextra_force, lbig_discard, lrhill_present)
     IF (.NOT. lrhill_present) THEN
          WRITE(*, *) "SWIFTEST Error:"
          WRITE(*, *) "   Integrator SyMBA requires planet Hill sphere radii on input"
          CALL util_exit(FAILURE)
     END IF
     ! Read in the total number of bodies from the input files
     CALL io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)

     ! Create arrays of data structures big enough to store the number of bodies we are adding
     CALL symba_pl_allocate(symba_plA,nplmax)
     CALL symba_merger_allocate(mergeadd_list,nplmax)
     CALL symba_merger_allocate(mergesub_list,nplmax)
     CALL symba_plplenc_allocate(plplenc_list, nplmax)
     CALL symba_pltpenc_allocate(pltpenc_list, ntpmax)
     ALLOCATE(discard_plA(11,npl))
     ALLOCATE(discard_tpA(11,ntp))
     ALLOCATE(discard_plA_id_status(2,npl))
     ALLOCATE(discard_tpA_id_status(2,ntp))


     IF (ntp > 0) THEN
          CALL symba_tp_allocate(symba_tpA, ntpmax)
     END IF

     ! Reads in initial conditions of all massive bodies from input file and fills the linked list
     CALL io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, symba_plA)
     WRITE(*, 100, ADVANCE = "NO") "Enter the smallest mass to self-gravitate: "
     READ(*, *) mtiny

     ! Reorder linked list by mass 
     CALL symba_reorder_pl(npl, symba_plA)
     CALL io_init_tp(intpfile, in_type, ntp, symba_tpA)
     CALL util_valid(npl, ntp, symba_plA%helio%swiftest, symba_tpA%helio%swiftest)
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
     IF (istep_out > 0) THEN
          CALL io_write_frame(t, npl, ntp, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, outfile, &
          out_type, out_form, out_stat)
     END IF
     IF (out_stat == "OLD") then
        OPEN(UNIT = egyiu, FILE = ENERGY_FILE, FORM = "FORMATTED", STATUS = "OLD", ACTION = "WRITE", POSITION = "APPEND")
     ELSE 
        OPEN(UNIT = egyiu, FILE = ENERGY_FILE, FORM = "FORMATTED", STATUS = "REPLACE", ACTION = "WRITE")
     END IF
 300 FORMAT(7(1X, E23.16))
 310 FORMAT(7(1X, A23))
     WRITE(egyiu,310) "#t","ke","pe","te","htotx","htoty","htotz"
     WRITE(*, *) " *************** MAIN LOOP *************** "
     CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
     WRITE(egyiu,300) t, ke, pe, te, htot
     DO WHILE ((t < tstop) .AND. ((ntp0 == 0) .OR. (ntp > 0)))


          CALL symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, &
           j4rp4, dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
           eoffset, mtiny, encounter_file, out_type)
          iloop = iloop + 1
          IF (iloop == LOOPMAX) THEN
               tbase = tbase + iloop*dt
               iloop = 0
          END IF
          t = tbase + iloop*dt
          ! Take the merger info and create fragments
          ! CALL some subroutine that returns the number of fragments and an array of new bodies (swifter_pl type)                     
          !IF (lfragmentation) THEN
               !CALL symba_fragmentation(t, npl, nplmax, ntp, ntpmax, symba_pl1P, nplplenc, plplenc_list)                               ! CHECK THIS 
               ! update nplmax to add in the new number of bodies
               ! add new bodies into the current body linked list as in CALL symba_setup
               ! reorder bodies (if that is not already going to happen..check the discard subroutines
               
               !CALL symba_add(npl, mergeadd_list, nmergeadd, symba_pl1P, swifter_pl1P, mtiny)                                          ! CHECK THIS 
          !END IF
          ldiscard = .FALSE. 
          ldiscard_tp = .FALSE.
          CALL symba_discard_merge_pl(t, npl, symba_plA, nplplenc, plplenc_list)                                  ! CHECK THIS 
          CALL symba_discard_pl(t, npl, nplmax, nsppl, symba_plA, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,    &    ! CHECK THIS 
               qmin_ahi, j2rp2, j4rp4, eoffset)
          CALL symba_discard_tp(t, npl, ntp, nsptp, symba_plA, symba_tpA, dt, rmin, rmax, rmaxu, qmin, qmin_coord, &    ! CHECK THIS 
               qmin_alo, qmin_ahi, lclose, lrhill_present)
          IF ((ldiscard .eqv. .TRUE.) .or. (ldiscard_tp .eqv. .TRUE.)) THEN
               CALL symba_rearray(t, npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
                    discard_tpA, discard_plA_id_status,discard_tpA_id_status, NPLMAX, j2rp2, j4rp4)
               CALL io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, nmergesub, symba_plA,discard_plA, &
               	discard_tpA, mergeadd_list, mergesub_list, DISCARD_FILE, lbig_discard, discard_plA_id_status, &
               	discard_tpA_id_status) 
               nmergeadd = 0
               nmergesub = 0
               nsppl = 0
               nsptp = 0
               CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
               WRITE(egyiu,300) t, ke, pe, te, htot
          END IF
          IF (istep_out > 0) THEN
               iout = iout - 1
               IF (iout == 0) THEN
                    CALL io_write_frame(t, npl, ntp, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, outfile, out_type, &
                     out_form, out_stat)
                    iout = istep_out
                  CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
                  WRITE(egyiu,300) t, ke, pe, te, htot
               END IF
          END IF
          IF (istep_dump > 0) THEN
               idump = idump - 1
               IF (idump == 0) THEN
                    tfrac = (t - t0)/(tstop - t0)
                    WRITE(*, 200) t, tfrac, npl, ntp
 200                FORMAT(" Time = ", ES12.5, "; fraction done = ", F5.3, "; Number of active pl, tp = ", I5, ", ", I5)
                    CALL io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form,        &
                         istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,               &
                         encounter_file, lextra_force, lbig_discard, lrhill_present)
                    CALL io_dump_pl(npl, symba_plA%helio%swiftest, lclose, lrhill_present)
                    IF (ntp > 0) CALL io_dump_tp(ntp, symba_tpA%helio%swiftest)
                    idump = istep_dump
               END IF
          END IF
          plplenc_list%lvdotr(:) = .FALSE.
          plplenc_list%status(:) = 0
          plplenc_list%level(:) = 0
          plplenc_list%index1(:) = 0
          plplenc_list%index2(:) = 0
          plplenc_list%enc_child(:) = 0 
          plplenc_list%enc_parent(:) = 0

          pltpenc_list%lvdotr(:) = .FALSE.
          pltpenc_list%status(:) = 0
          pltpenc_list%level(:) = 0
          pltpenc_list%indexpl(:) = 0
          pltpenc_list%indextp(:) = 0

          mergeadd_list%name(:) = 0
          mergeadd_list%index_ps(:) = 0
          mergeadd_list%status(:) = 0
          mergeadd_list%ncomp(:) = 0
          mergeadd_list%xh(:,:) = 0
          mergeadd_list%vh(:,:) = 0
          mergeadd_list%mass(:) = 0
          mergeadd_list%radius(:) = 0

          mergesub_list%name(:) = 0
          mergesub_list%index_ps(:) = 0
          mergesub_list%status(:) = 0
          mergesub_list%ncomp(:) = 0
          mergesub_list%xh(:,:) = 0
          mergesub_list%vh(:,:) = 0
          mergesub_list%mass(:) = 0
          mergesub_list%radius(:) = 0

          discard_plA(:,:) = 0
          discard_tpA(:,:) = 0
          discard_plA_id_status(:,:) = 0
          discard_tpA_id_status(:,:) = 0
     END DO
     CALL io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form, istep_dump, j2rp2,    &
          j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, encounter_file, lextra_force, lbig_discard,     &
          lrhill_present)
     CALL io_dump_pl(npl, symba_plA%helio%swiftest, lclose, lrhill_present)
     IF (ntp > 0) CALL io_dump_tp(ntp, symba_tpA%helio%swiftest)

     CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
     WRITE(egyiu,300) t, ke, pe, te, htot
     CLOSE(egyiu)

     CALL symba_pl_deallocate(symba_plA)
     CALL symba_merger_deallocate(mergeadd_list)
     CALL symba_merger_deallocate(mergesub_list)
     CALL symba_plplenc_deallocate(plplenc_list)
     CALL symba_pltpenc_deallocate(pltpenc_list)
     DEALLOCATE(discard_plA)
     DEALLOCATE(discard_plA_id_status)
     IF (ntp > 0) THEN
          CALL symba_tp_deallocate(symba_tpA)
          DEALLOCATE(discard_tpA)
          DEALLOCATE(discard_tpA_id_status)
     END IF

     CALL util_exit(SUCCESS)

     STOP

END PROGRAM swiftest_symba_omp
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
