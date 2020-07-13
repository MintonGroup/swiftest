program swiftest_driver
   !! author: David A. Minton
   !!
   !! Driver program for the Swiftest integrators. Unlike the earlier Swift and Swifter drivers, in Swiftest all integrators 
   !!    are run from this single program. 
   !!
   !! Adapted from Swifter by David E. Kaufmann's Swifter driver programs swifter_[bs,helio,ra15,rmvs,symba,tu4,whm].f90
   !! Adapted from Hal Levison and Martin Duncan's Swift driver programs
   use swiftest
   implicit none

   class(swiftest_nbody_system), allocatable  :: nbody_system     !! Polymorphic object containing the nbody system to be integrated
   class(swiftest_configuration), allocatable :: config
   integer(I4B)                               :: integrator       !! Integrator type code (see swiftest_globals for symbolic names)
   character(len=:),allocatable               :: config_file_name !! Name of the file containing user-defined configuration parameters
   integer(I4B)                               :: ierr             !! I/O error code 
   logical                                    :: lfirst           !! Flag indicating that this is the first time through the main loop
   integer(I8B)                               :: iloop            !! Loop counter
   integer(I8B)                               :: idump            !! Dump cadence counter
   integer(I8B)                               :: iout             !! Output cadence counter
   integer(I8B), parameter                    :: LOOPMAX = huge(iloop) !! Maximum loop value before resetting 
   real(DP)                                   :: start_wall_time  !! Wall clock time at start of execution
   real(DP)                                   :: finish_wall_time !! Wall clock time when execution has finished
   integer(I4B)                               :: iu               !! Unit number of binary file
   integer(I4B)                               :: ntp, npl
   character(*),parameter :: statusfmt  = '("Time = ", ES12.5, "; fraction done = ", F6.3, "; ' // &
                                             'Number of active pl, tp = ", I5, ", ", I5)'

   !> Define the maximum number of threads
   nthreads = 1            ! In the *serial* case
   !$ nthreads = omp_get_max_threads() ! In the *parallel* case
   !$ write(*,'(a)')   ' OpenMP parameters:'
   !$ write(*,'(a)')   ' ------------------'
   !$ write(*,'(a,i3,/)') ' Number of threads  = ', nthreads  

   ierr = io_get_args(integrator, config_file_name)
   if (ierr /= 0) then
      write(*,*) 'Error reading in arguments from the command line'
      call util_exit(FAILURE)
   end if
   !$ start_wall_time = omp_get_wtime()
   !> Read in the user-defined parameter file and the initial conditions of the system
   call setup_construct_system(nbody_system, config, integrator)
   call config%read_from_file(config_file_name, integrator)
   associate(t          => config%t, &
             t0         => config%t0, &
             dt         => config%dt, &
             tstop      => config%tstop, &
             istep_out  => config%istep_out, &
             istep_dump => config%istep_dump)  
      call nbody_system%initialize(config)
      lfirst = .true.
      t = t0
      iloop = 0
      iout = istep_out
      idump = istep_dump
      if (istep_out > 0) call nbody_system%write_frame(iu, config, t, dt)
      write(*, *) " *************** Main Loop *************** "
      do iloop = 1, LOOPMAX 
         ntp = nbody_system%tp%nbody
         npl = nbody_system%pl%nbody
         t = t0 + iloop * dt
         if (t > tstop) exit 
         !> Step the system forward in time
         call nbody_system%step(config)

         !> Evaluate any discards or mergers
         call nbody_system%discard(config)

         !> If the loop counter is at the output cadence value, append the data file with a single frame
         if (istep_out > 0) then
            iout = iout - 1
            if (iout == 0) then
               call nbody_system%write_frame(iu, config, t, dt)
               iout = istep_out
            end if
         end if

         !> If the loop counter is at the dump cadence value, dump the state of the system to a file in case it needs to be restarted
         if (istep_dump > 0) then
            idump = idump - 1
            if (idump == 0) then
               call nbody_system%dump(config, t, dt, statusfmt)
               idump = istep_dump
            end if
         end if
      end do

      !> Dump the final state of the system to file
      call nbody_system%dump(config, t, dt, statusfmt)
      !$ finish_wall_time = omp_get_wtime()
      !$ write(*,*) 'Time: ', finish_wall_time - start_wall_time
   end associate
   call util_exit(SUCCESS)

   stop

end program swiftest_driver
