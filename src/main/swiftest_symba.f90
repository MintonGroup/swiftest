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
     use omp_lib
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT)      :: lclose         ! Check for planet-test particle encounters
     LOGICAL(LGT)      :: lextra_force   ! Use user-supplied force routines
     LOGICAL(LGT)      :: lbig_discard   ! Dump planet data with discards
     LOGICAL(LGT)      :: lrhill_present ! Hill's sphere radius present
     LOGICAL(LGT)      :: lpython        ! Python flag to output pl_out.dat and tp_out.dat 
     LOGICAL(LGT)      :: lenergy        ! Python flag to output energy.out 
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
     LOGICAL(LGT)                                               :: lfirst, lfrag_add
     INTEGER(I4B)                                               :: npl, ntp, ntp0, nsppl, nsptp, iout, idump, iloop
     INTEGER(I4B)                                               :: nplplenc, npltpenc, nmergeadd, nmergesub, fragmax
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
     INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: k_plpl, k_pltp
     INTEGER(I4B) :: num_plpl_comparisons, num_pltp_comparisons
     REAL(DP) :: start, finish

     INTEGER(I4B), PARAMETER                                    :: egyiu = 72
! Executable code
     CALL CPU_TIME(START)
     CALL util_version
     nthreads = 1                        ! In the *serial* case
     WRITE(*, 100, ADVANCE = "NO") "Enter name of parameter data file: "
     READ(*, 100) inparfile
 100 FORMAT(A)
     inparfile = TRIM(ADJUSTL(inparfile))
     ! Read in the param.in file and get simulation parameters
     CALL io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile, out_type,      &
          out_form, out_stat, istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,          &
          encounter_file, lextra_force, lbig_discard, lrhill_present, mtiny, lpython, lenergy)
     IF (.NOT. lrhill_present) THEN
          WRITE(*, *) "SWIFTEST Error:"
          WRITE(*, *) "   Integrator SyMBA requires planet Hill sphere radii on input"
          CALL util_exit(FAILURE)
     END IF
     ! Read in the total number of bodies from the input files
     CALL io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)

     ! Create arrays of data structures big enough to store the number of bodies we are adding
     CALL symba_pl_allocate(symba_plA,npl)
     CALL symba_merger_allocate(mergeadd_list,npl)
     CALL symba_merger_allocate(mergesub_list,npl)
     CALL symba_plplenc_allocate(plplenc_list, npl)
     CALL symba_pltpenc_allocate(pltpenc_list, ntp)


     IF (ntp > 0) THEN
          CALL symba_tp_allocate(symba_tpA, ntpmax)
     END IF

     ! Reads in initial conditions of all massive bodies from input file
     CALL io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, symba_plA)

     ! Reorder by mass 
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
     fragmax = 0 
     IF (istep_out > 0) THEN
          CALL io_write_frame(t, npl, ntp, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, outfile, &
          out_type, out_form, out_stat)
          IF (lpython) THEN
               call python_io_write_frame_pl(t, symba_plA, npl, out_stat)
               IF (ntp>0) call python_io_write_frame_tp(t, symba_tpA, ntp, out_stat)
          END IF
     END IF
     IF (out_stat == "OLD") then
        OPEN(UNIT = egyiu, FILE = ENERGY_FILE, FORM = "FORMATTED", STATUS = "OLD", ACTION = "WRITE", POSITION = "APPEND")
     END IF
     ELSE 
        OPEN(UNIT = egyiu, FILE = ENERGY_FILE, FORM = "FORMATTED", STATUS = "REPLACE", ACTION = "WRITE")
 300 FORMAT(7(1X, E23.16))
 310 FORMAT(7(1X, A23))
     start = omp_get_wtime()
     ! call cpu_time(start)     
     nplm = count(symba_plA%helio%swiftest%mass>mtiny)
     CALL util_dist_index_plpl(npl, nplm, num_plpl_comparisons, k_plpl)

     CALL util_dist_index_pltp(nplm, ntp, num_pltp_comparisons, k_pltp)
     WRITE(*, *) " *************** MAIN LOOP *************** "
     IF (lenergy) THEN 
          CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
          WRITE(egyiu,300) t, ke, pe, te, htot
     END IF
     CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
     DO WHILE ((t < tstop) .AND. ((ntp0 == 0) .OR. (ntp > 0)))
          if(num_plpl_comparisons > 100000)then
            CALL symba_step_eucl(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, &
                j4rp4, dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
                eoffset, mtiny, encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)
          else
            CALL symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, &
                j4rp4, dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
                eoffset, mtiny, encounter_file, out_type)
          endif
          
          iloop = iloop + 1
          IF (iloop == LOOPMAX) THEN
               tbase = tbase + iloop*dt
               iloop = 0
          END IF
          t = tbase + iloop*dt
          ldiscard = .FALSE. 
          ldiscard_tp = .FALSE.
          CALL symba_discard_merge_pl(t, npl, symba_plA, nplplenc, plplenc_list)                                  ! CHECK THIS 
          CALL symba_discard_pl(t, npl, nplmax, nsppl, symba_plA, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,    &    ! CHECK THIS 
               qmin_ahi, j2rp2, j4rp4, eoffset)
          CALL symba_discard_tp(t, npl, ntp, nsptp, symba_plA, symba_tpA, dt, rmin, rmax, rmaxu, qmin, qmin_coord, &    ! CHECK THIS 
               qmin_alo, qmin_ahi, lclose, lrhill_present)
           ENDIF

            ENDIF
                  pltpenc_list%indextp(1:1+npltpenc) = 0
                  pltpenc_list%indexpl(1:1+npltpenc) = 0
                  pltpenc_list%level(1:1+npltpenc) = 0
                  pltpenc_list%lvdotr(1:1+npltpenc) = .FALSE.
                  pltpenc_list%status(1:1+npltpenc) = 0
            IF(npltpenc > 0)THEN

               plplenc_list%enc_parent(1:1+nplplenc) = 0
               plplenc_list%index2(1:1+nplplenc) = 0
               plplenc_list%enc_child(1:1+nplplenc) = 0 
               plplenc_list%index1(1:1+nplplenc) = 0
               plplenc_list%level(1:1+nplplenc) = 0
               plplenc_list%status(1:1+nplplenc) = 0
               plplenc_list%lvdotr(1:1+nplplenc) = .FALSE.
               
          IF(nplplenc > 0)THEN
          IF ((ldiscard .eqv. .TRUE.) .or. (ldiscard_tp .eqv. .TRUE.) .or. (lfrag_add .eqv. .TRUE.)) THEN
               CALL symba_rearray(t, npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
                    discard_tpA, NPLMAX, j2rp2, j4rp4)
               IF ((ldiscard .eqv. .TRUE.) .or. (ldiscard_tp .eqv. .TRUE.)) THEN
               		CALL io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, nmergesub, symba_plA, &
               			discard_plA, discard_tpA, mergeadd_list, mergesub_list, DISCARD_FILE, lbig_discard) 
               DEALLOCATE(k_plpl)
               nplm = count(symba_plA%helio%swiftest%mass(1:npl)>mtiny)
               CALL util_dist_index_plpl(npl, nplm, num_plpl_comparisons, k_plpl)
               if(ntp>0)then
                    DEALLOCATE(k_pltp)
                    CALL util_dist_index_pltp(nplm, ntp, num_pltp_comparisons, k_pltp)
               endif
               
               mergeadd_list%name(1:1+nmergeadd) = 0
               mergeadd_list%index_ps(1:1+nmergeadd) = 0
               mergeadd_list%ncomp(1:1+nmergeadd) = 0
               mergeadd_list%status(1:1+nmergeadd) = 0
               mergeadd_list%xh(:,1:1+nmergeadd) = 0
               mergeadd_list%vh(:,1:1+nmergeadd) = 0
               mergeadd_list%mass(1:1+nmergeadd) = 0
               mergeadd_list%radius(1:1+nmergeadd) = 0

               mergesub_list%name(1:1+nmergesub) = 0
               mergesub_list%index_ps(1:1+nmergesub) = 0
               mergesub_list%status(1:1+nmergesub) = 0
               mergesub_list%ncomp(1:1+nmergesub) = 0
               mergesub_list%vh(:,1:1+nmergesub) = 0
               mergesub_list%xh(:,1:1+nmergesub) = 0
               mergesub_list%mass(1:1+nmergesub) = 0
               mergesub_list%radius(1:1+nmergesub) = 0

               discard_plA_id_status(:,:) = 0
               discard_plA(:,:) = 0

               if(ntp>0)then

                  discard_tpA_id_status(:,:) = 0
                  discard_tpA(:,:) = 0
               endif


               nmergeadd = 0
               nmergesub = 0
               nsppl = 0
               nsptp = 0
               END IF 
               IF (lenergy) THEN 
                    CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
                    WRITE(egyiu,300) t, ke, pe, te, htot
               END IF
          END IF
          IF (istep_out > 0) THEN
               iout = iout - 1
               IF (iout == 0) THEN
                    CALL io_write_frame(t, npl, ntp, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, outfile, out_type, &
                     out_form, out_stat)
                    iout = istep_out
                    IF (lpython) THEN 
                         call python_io_write_frame_pl(t, symba_plA, npl, out_stat= "APPEND")
                         IF (ntp>0) call python_io_write_frame_tp(t, symba_tpA, ntp, out_stat= "APPEND")
                    END IF 
                    IF (lenergy) THEN 
                         CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
                         WRITE(egyiu,300) t, ke, pe, te, htot
                    END IF
                  CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
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
                         encounter_file, lextra_force, lbig_discard, lrhill_present, mtiny, lpython)
                    CALL io_dump_pl(npl, symba_plA%helio%swiftest, lclose, lrhill_present)
                    IF (ntp > 0) CALL io_dump_tp(ntp, symba_tpA%helio%swiftest)
                    idump = istep_dump
               END IF
          END IF

     END DO
     CALL io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form, istep_dump, j2rp2,    &
          j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, encounter_file, lextra_force, lbig_discard,     &
          lrhill_present, mtiny, lpython)
     CALL io_dump_pl(npl, symba_plA%helio%swiftest, lclose, lrhill_present)
     IF (ntp > 0) CALL io_dump_tp(ntp, symba_tpA%helio%swiftest)
     IF (lenergy) THEN 
          CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
          WRITE(egyiu,300) t, ke, pe, te, htot
          close(egyiu)
     END IF
     
     CALL symba_pl_deallocate(symba_plA)
     CALL symba_merger_deallocate(mergeadd_list)
     CALL symba_merger_deallocate(mergesub_list)
     CALL symba_plplenc_deallocate(plplenc_list)
     CALL symba_pltpenc_deallocate(pltpenc_list)
     IF (ntp > 0) THEN
          CALL symba_tp_deallocate(symba_tpA)
     END IF
     finish = omp_get_wtime()
     ! call cpu_time(finish)     
     print *,'Time: ',finish-start
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
