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

   type(swiftest_configuration)              :: config           !! Object containing user-defined configuration parameters
   class(swiftest_nbody_system), allocatable :: nbody_system     !! Polymorphic object containing the nbody system to be integrated
   integer(I4B)                              :: integrator       !! Integrator type code (see swiftest_globals for symbolic names)
   character(len=:),allocatable              :: config_file_name !! Name of the file containing user-defined configuration parameters
   integer(I4B)                              :: ierr             !! I/O error code 
   logical                                   :: lfirst           !! Flag indicating that this is the first time through the main loop
   integer(I8B)                              :: iout             !! I/O output cadence counter 
   integer(I8B)                              :: idump            !! Dump cadence counter
   integer(I8B)                              :: iloop            !! Main loop counter
   integer(I8B), parameter                   :: LOOPMAX = huge(iloop) !! Maximum loop value before resetting 
   real(DP)                                  :: t                !! Simulation time
   real(DP)                                  :: dt               !! Simulation step size
   real(DP)                                  :: tfrac            !! Fraction of time remaining in the integration
   real(DP)                                  :: tbase            !! Resets the internal start time if the loop counter reaches its maximum
   real(DP)                                  :: start_wall_time  !! Wall clock time at start of execution
   real(DP)                                  :: finish_wall_time !! Wall clock time when execution has finished
   !character(len=*), parameter               :: fmt_dump = '(" Time = ", es12.5, "; Fraction done = ", f5.3, "; Number of active pl, tp = ", i5, ", ", i5)'


   ierr = io_get_command_line_arguments(integrator, config_file_name)
   if (ierr == 0) then
      !$ start_wall_time = omp_get_wtime()
      !> Read in the user-defined parameter file and the initial conditions of the system
      call config%read_from_file(config_file_name, integrator)
      call nbody_system%construct(config, integrator)
      call nbody_system%initialize(config)

      lfirst = .true.
      t = config%t0
      tbase = config%t0
      iloop = 0
      iout = config%istep_out
      idump = config%istep_dump
      if (istep_out > 0) call nbody_system%write_frame(config, t, dt)
      write(*, *) " *************** Main Loop *************** "
      do while ((t < config%tstop) .and. nbody_system%lkeep_going) 

         !> Step the system forward in time
         call nbody_system%step(config, t, dt)

         !> Advance the loop counter and time value
         iloop = iloop + 1
         if (iloop == LOOPMAX) then  !! Reset loop counter if necessary (this should be rare, as iloop is an 8-byte integer)
            tbase = tbase + iloop * dt
            iloop = 0
         end if
         t = tbase + iloop * dt

         !> Evaluate any discards or mergers
         call nbody_system%discard(config, t, dt)
         if (nbody_system%ldiscard) call nbody_system%write_discard(config, t, ct)

         !> If the loop counter is at the output cadence value, append the data file with a single frame
         if (istep_out > 0) then
            iout = iout - 1
            if (iout == 0) then
               call nbody_system%write_frame(config, t, dt)
               iout = istep_out
            end if
         end if

         !> If the loop counter is at the dump cadence value, dump the state of the system to a file in case it needs to be restarted
         if (istep_dump > 0) then
            idump = idump - 1
            if (idump == 0) then
               tfrac = (t - t0) / (tstop - t0)
               call nbody_system%dump(config, t, dt, tfrac)
               call config%dump(t)
               idump = istep_dump
            end if
         end if
      end do

      !> Dump the final state of the system to file
      call nbody_system%dump(config, t, dt, tfrac)
      call config%dump(t)
      !$ finsih_wall_time = omp_get_wtime()
      !$ write(*,*) 'Time: ', finish_wall_time - start_wall_time
   end if
   call util_exit(SUCCESS)

   stop

end program swiftest_driver
