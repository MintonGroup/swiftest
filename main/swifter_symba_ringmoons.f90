!**********************************************************************************************************************************
!
!  Unit Name   : swifter_symba_ringmoons
!  Unit Type   : program
!  Project     : Swifter
!  Package     : main
!  Language    : Fortran 90/95
!
!  Description : Driver program for the Symplectic Massive Body Algorithm with
!  OpenMP parallelization and RING-MOONS
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
!  Invocation  : % swifter_symba
!
!  Notes       : Reference: Duncan, M. J., Levison, H. F. & Lee, M. H. 1998. Astron. J., 116, 2067.
!
!                Adapted from Hal Levison and Martin Duncan's Swift program swift_symba5.f
!                OpenMP parallelization by David Minton
!
!**********************************************************************************************************************************
PROGRAM swifter_symba_ringmoons

! Modules
     USE module_parameters
     USE module_swifter
     USE module_symba
     USE module_random_access
     USE module_interfaces
     USE module_ringmoons
     USE module_ringmoons_interfaces
     !Added by D. Minton
     !$ USE omp_lib
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT)      :: lclose         ! Check for planet-test particle encounters
     LOGICAL(LGT)      :: lextra_force   ! Use user-supplied force routines
     LOGICAL(LGT)      :: lbig_discard   ! Dump planet data with discards
     LOGICAL(LGT)      :: lrhill_present ! Hill's sphere radius present
     LOGICAL(LGT)      :: lrotation      ! Rotation parameters present
     INTEGER(I4B)      :: nplmax         ! Maximum number of planets
     INTEGER(I4B)      :: ntpmax         ! Maximum number of test particles
     INTEGER(I4B)      :: istep_out      ! Time steps between binary outputs
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
     CHARACTER(STRMAX) :: ring_outfile   ! Name of ringmoons output binary file

! Internals
     LOGICAL(LGT)                                      :: lfirst,lrmfirst
     INTEGER(I4B)                                      :: npl, ntp, ntp0, nsppl, nsptp, iout, iloop, i
     INTEGER(I4B)                                      :: nplplenc, npltpenc, nmergeadd, nmergesub
     REAL(DP)                                          :: t, tfrac, tbase, mtiny, ke, pe, te, eoffset
     REAL(DP), DIMENSION(NDIM)                         :: htot
     CHARACTER(STRMAX)                                 :: inparfile
     TYPE(swifter_pl), POINTER                         :: swifter_pl1P
     TYPE(swifter_tp), POINTER                         :: swifter_tp1P
     TYPE(symba_pl), DIMENSION(:), ALLOCATABLE, TARGET :: symba_plA
     TYPE(symba_tp), DIMENSION(:), ALLOCATABLE, TARGET :: symba_tpA
     TYPE(symba_pl), POINTER                           :: symba_pl1P, symba_pld1P
     TYPE(symba_tp), POINTER                           :: symba_tp1P, symba_tpd1P
     TYPE(symba_plplenc), DIMENSION(NENMAX)            :: plplenc_list
     TYPE(symba_pltpenc), DIMENSION(NENMAX)            :: pltpenc_list
     TYPE(symba_merger), DIMENSION(:), ALLOCATABLE     :: mergeadd_list, mergesub_list
      real(DP)                                         :: Merror,Lerror

     TYPE(ringmoons_ring) :: ring
     TYPE(ringmoons_seeds) :: seeds

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
     CALL io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile, out_type,      &
          out_form, out_stat, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,          &
          encounter_file, lextra_force, lbig_discard, lrhill_present, lrotation, mtiny, ring_outfile)
     IF (.NOT. lrhill_present) THEN
          WRITE(*, *) "SWIFTER Error:"
          WRITE(*, *) "   Integrator SyMBA requires planet Hill sphere radii on input"
          CALL util_exit(FAILURE)
     END IF
     CALL io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)
     ALLOCATE(symba_plA(nplmax), mergeadd_list(nplmax), mergesub_list(nplmax))
     CALL set_point(symba_plA)
     IF (ntp > 0) THEN
          ALLOCATE(symba_tpA(ntpmax))
          CALL set_point(symba_tpA)
     END IF
     CALL symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1P, symba_tp1P, swifter_pl1P, swifter_tp1P)
     CALL io_init_pl(inplfile, in_type, lclose, lrhill_present, lrotation, npl, swifter_pl1P)
     CALL symba_reorder_pl(npl, symba_pl1P)
     CALL io_init_tp(intpfile, in_type, ntp, swifter_tp1P)
     CALL util_valid(npl, ntp, swifter_pl1P, swifter_tp1P)

     !Read in RING-MOONS parameters and data
     call ringmoons_io_init_ring(swifter_pl1P,ring,seeds)
     call ringmoons_ring_construct(swifter_pl1P,ring)
     call ringmoons_viscosity(ring)
     call ringmoons_seed_construct(swifter_pl1P,ring,seeds) 
     DESTRUCTION_EVENT = .false.
     DESTRUCTION_COUNTER = 0

     !Set up integration
     lfirst = .TRUE.
     lrmfirst = .true.
     ntp0 = ntp
     t = t0
     tbase = t0
     iloop = 0
     iout = istep_out
     nmergeadd = 0
     nmergesub = 0
     nsppl = 0
     nsptp = 0
     eoffset = 0.0_DP
     Merror = 0.0_DP
     Lerror = 0.0_DP
     tfrac = 0.0_DP
     NULLIFY(symba_pld1P, symba_tpd1P)
      if (istep_out > 0) then
         CALL io_write_frame(t, npl, ntp, swifter_pl1P, swifter_tp1P, outfile, out_type, out_form, out_stat)
         call ringmoons_io_write_frame(t, ring, seeds, ring_outfile, out_stat = "NEW")
         WRITE(*, 200) t, tfrac, npl, ntp
         WRITE(*,300) Merror,Lerror,count(seeds%active)
      end if
     WRITE(*, *) " *************** MAIN LOOP *************** "
     DO WHILE (t < tstop) 
          CALL ringmoons_step(t,swifter_pl1P,ring,seeds,dt,lrmfirst,Merror,Lerror)
          CALL symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_pl1P, symba_tp1P, j2rp2, j4rp4, dt,    &
               nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset,       &
               mtiny, encounter_file, out_type)
          iloop = iloop + 1
          IF (iloop == LOOPMAX) THEN
               tbase = tbase + iloop*dt
               iloop = 0
          END IF
          t = tbase + iloop*dt
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
                    tfrac = (t - t0) / (tstop - t0)
                    WRITE(*, 200) t, tfrac, npl, ntp
                    WRITE(*,300) Merror,Lerror,count(seeds%active)
 200                FORMAT(" Time = ", ES12.5, "; fraction done = ", F5.3, "; Number of active pl, tp = ", I5, ", ", I5)
 300                FORMAT("dM/M0 = ", ES12.5, ";  dL/L0 = ",ES12.5,       "; Number or active seeds  = ",I5)
                    CALL io_write_frame(t, npl, ntp, swifter_pl1P, swifter_tp1P, outfile, out_type, out_form, out_stat)
                    call ringmoons_io_write_frame(t, ring, seeds, ring_outfile, out_stat = "APPEND")
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

END PROGRAM swifter_symba_ringmoons
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
