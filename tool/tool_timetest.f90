!**********************************************************************************************************************************
!
!  Unit Name   : tool_timetest
!  Unit Type   : program
!  Project     : Swifter
!  Package     : main
!  Language    : Fortran 90/95
!
!  Description : Tests of execution time of different methods for looping over SyMBA body lists
!  OpenMP parallelization
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
MODULE module_timetest
   !$ USE omp_lib
   IMPLICIT NONE

   CONTAINS

      SUBROUTINE shadow_array(npl, nplmax, symba_pl1P)
         USE module_parameters
         USE module_symba
         IMPLICIT NONE

         INTEGER(I4B),INTENT(IN)      :: npl,nplmax
         TYPE(symba_pl), POINTER      :: symba_pl1P
         TYPE(swifter_pl), POINTER :: swifter_pliP,swifter_pljP
         TYPE(swifter_pl), POINTER :: swifter_pl1P
         TYPE(helio_pl), POINTER   :: helio_pl1P,helio_pliP,helio_pljP
         TYPE(symba_pl), POINTER   :: symba_pliP, symba_pljP
         INTEGER(I4B)              :: i,j,k

         symba_pliP => symba_pl1P
         helio_pl1P => symba_pl1P%helio
         !Added by D. Minton
         swifter_pl1P => helio_pl1P%swifter
         IF (ALLOCATED(symba_pl1P%symba_plPA)) DEALLOCATE(symba_pl1P%symba_plPA)
         ALLOCATE(symba_pl1P%symba_plPA(npl))
         IF (ALLOCATED(swifter_pl1P%swifter_plPA)) DEALLOCATE(swifter_pl1P%swifter_plPA)
         ALLOCATE(swifter_pl1P%swifter_plPA(npl))
         IF (ALLOCATED(helio_pl1P%helio_plPA)) DEALLOCATE(helio_pl1P%helio_plPA)
         ALLOCATE(helio_pl1P%helio_plPA(npl))
         !^^^^^^^^^^^^^^^^^^
         DO i = 1, npl
             symba_pl1P%symba_plPA(i)%thisP => symba_pliP
             helio_pl1p%helio_plPA(i)%thisP => symba_pliP%helio
             swifter_pl1P%swifter_plPA(i)%thisP => symba_pliP%helio%swifter
             symba_pliP%helio%swifter%vb(:) = 0.0_DP
             symba_pliP => symba_pliP%nextP
             
         END DO
         

        symba_pliP => symba_pl1P
        !$OMP PARALLEL DO DEFAULT(PRIVATE) &
        !$OMP SHARED(npl,symba_pl1P) 
        DO i = 2, npl
             symba_pliP => symba_pl1P%symba_plPA(i)%thisP
             helio_pliP => symba_pliP%helio
             swifter_pliP => helio_pliP%swifter
             helio_pliP%ah(:) = (/9.0_DP, 8.5_DP, 7.0_DP/)
             DO j = i + 1, npl
                  symba_pljP=>symba_pl1P%symba_plPA(j)%thisP
                  helio_pljP => symba_pljP%helio
                  swifter_pljP => helio_pljP%swifter
                  do k = 1, 1
                     swifter_pljP%vb(:) = swifter_pliP%vb(:) + helio_pljP%ah(:)*3.7_DP !Arbitrary calculation
                  end do
             END DO
        END DO
       !$OMP END PARALLEL DO
   
      END SUBROUTINE shadow_array


      SUBROUTINE get_point_method(npl, nplmax, symba_pl1P)
         USE module_parameters
         USE module_symba
         USE module_random_access
         IMPLICIT NONE

         INTEGER(I4B),INTENT(IN)      :: npl,nplmax
         TYPE(symba_pl), POINTER      :: symba_pl1P
         TYPE(swifter_pl), POINTER :: swifter_pliP,swifter_pljP
         TYPE(swifter_pl), POINTER :: swifter_pl1P
         TYPE(helio_pl), POINTER   :: helio_pl1P,helio_pliP,helio_pljP
         TYPE(symba_pl), POINTER   :: symba_pliP, symba_pljP
         INTEGER(I4B)              :: i,j,k

         symba_pliP => symba_pl1P
         helio_pl1P => symba_pl1P%helio
         DO i = 1, npl
             symba_pliP%helio%swifter%vb(:) = 0.0_DP
             symba_pliP => symba_pliP%nextP
         END DO
         

        !$OMP PARALLEL DO DEFAULT(PRIVATE) &
        !$OMP SHARED(npl,symba_pl1P) 
        DO i = 2, npl
             symba_pliP => symba_pl1P
             call get_point(i,symba_pliP)
             helio_pliP => symba_pliP%helio
             swifter_pliP => helio_pliP%swifter
             helio_pliP%ah(:) = (/9.0_DP, 8.5_DP, 7.0_DP/)
             DO j = i + 1, npl
                  symba_pljP => symba_pl1P
                  helio_pljP => symba_pljP%helio
                  call get_point(j,symba_pljP)
                  swifter_pljP => helio_pljP%swifter
                  do k = 1, 1
                     swifter_pliP%vb(:) = swifter_pljP%vb(:) + helio_pljP%ah(:)*3.7_DP !Arbitrary calculation
                  end do
             END DO
        END DO
       !$OMP END PARALLEL DO
   
      END SUBROUTINE get_point_method


      SUBROUTINE forall_method(npl, nplmax, symba_pl1P)
         USE module_parameters
         USE module_symba
         USE module_random_access
         IMPLICIT NONE

         INTEGER(I4B),INTENT(IN)      :: npl,nplmax
         TYPE(symba_pl), POINTER      :: symba_pl1P
         TYPE(swifter_pl), POINTER :: swifter_pliP, swifter_pljP
         TYPE(swifter_pl), POINTER :: swifter_pl1P
         TYPE(helio_pl), POINTER   :: helio_pl1P,helio_pliP,helio_pljP
         TYPE(symba_pl), POINTER   :: symba_pliP, symba_pljP
         INTEGER(I4B)              :: i,j,k

         symba_pliP => symba_pl1P
         helio_pl1P => symba_pl1P%helio
         DO i = 1, npl
             symba_pliP%helio%swifter%vb(:) = 0.0_DP
             symba_pliP => symba_pliP%nextP
         END DO
         

        symba_pliP => symba_pl1P
        symba_pljP => symba_pl1P
        !$OMP WORKSHARE
        FORALL(i=2:npl,j=2:npl,j>i)
            lsymba_pl1P(i)%helio%ah(:) = (/9.0_DP, 8.5_DP, 7.0_DP/)
            lsymba_pl1P(i)%helio%swifter%vb(:) = lsymba_pl1P(j)%helio%swifter%vb(:) + lsymba_pl1P(j)%helio%ah(:)*3.7_DP !Arbitrary calculation
        END FORALL
        !$OMP END WORKSHARE
   
      END SUBROUTINE forall_method


END MODULE module_timetest
       


PROGRAM tool_timetest

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
     INTEGER(I4B)                                      :: nplplenc, npltpenc, nmergeadd, nmergesub,i
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
     REAL(DP)                                          :: tstart, tend
     REAL(DP),DIMENSION(:),ALLOCATABLE                 :: tshadow,tgetpoint,tforall
     

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
     ALLOCATE(tshadow(nthreads))
     ALLOCATE(tgetpoint(nthreads))
     ALLOCATE(tforall(nthreads))
     WRITE(*, 100, ADVANCE = "NO") "Enter name of parameter data file: "
     READ(*, 100) inparfile
 100 FORMAT(A)
     inparfile = TRIM(ADJUSTL(inparfile))
     CALL io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile, out_type,      &
          out_form, out_stat, istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,          &
          encounter_file, lextra_force, lbig_discard, lrhill_present)
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
     CALL io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, swifter_pl1P)
     mtiny = 1e-16_DP
     CALL symba_reorder_pl(npl, symba_pl1P)
     CALL io_init_tp(intpfile, in_type, ntp, swifter_tp1P)
     CALL util_valid(npl, ntp, swifter_pl1P, swifter_tp1P)
     !$ write(*,*) 'nthreads tshadow tgetpoint'
     !$ DO i = 1,nthreads
     !$ CALL OMP_SET_NUM_THREADS(i)
     lfirst = .TRUE.
     ntp0 = ntp
     t = t0
     tbase = t0
     iloop = 0
     NULLIFY(symba_pld1P, symba_tpd1P)

     !$ tstart = omp_get_wtime()
     DO WHILE (t < tstop)
          iloop = iloop + 1
          IF (iloop == LOOPMAX) THEN
               tbase = tbase + iloop*dt
               iloop = 0
          END IF
          t = tbase + iloop*dt
          CALL shadow_array(npl, nplmax, symba_pl1P)
     END DO
     !$ tend = omp_get_wtime()
     !$ tshadow(i) = tend - tstart

     t = t0
     tbase = t0
     iloop = 0

     !$ tstart = omp_get_wtime()
     DO WHILE ((t < tstop) .AND. ((ntp0 == 0) .OR. (ntp > 0)))
          iloop = iloop + 1
          IF (iloop == LOOPMAX) THEN
               tbase = tbase + iloop*dt
               iloop = 0
          END IF
          t = tbase + iloop*dt
          CALL get_point_method(npl, nplmax, symba_pl1P)
     END DO
     !$ tend = omp_get_wtime()
     !$ tgetpoint(i) = tend - tstart


     t = t0
     tbase = t0
     iloop = 0

     !$ tstart = omp_get_wtime()
     DO WHILE ((t < tstop) .AND. ((ntp0 == 0) .OR. (ntp > 0)))
          iloop = iloop + 1
          IF (iloop == LOOPMAX) THEN
               tbase = tbase + iloop*dt
               iloop = 0
          END IF
          t = tbase + iloop*dt
          CALL forall_method(npl, nplmax, symba_pl1P)
     END DO
     !$ tend = omp_get_wtime()
     !$ tforall(i) = tend - tstart




     !$ write(*,*) i,tshadow(i),tgetpoint(i),tforall(i)
     
     !$ end do
     IF (ALLOCATED(symba_plA)) DEALLOCATE(symba_plA)
     IF (ALLOCATED(mergeadd_list)) DEALLOCATE(mergeadd_list)
     IF (ALLOCATED(mergesub_list)) DEALLOCATE(mergesub_list)
     IF (ALLOCATED(symba_tpA)) DEALLOCATE(symba_tpA)
     DEALLOCATE(tshadow,tgetpoint,tforall)     

     STOP

END PROGRAM tool_timetest
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
