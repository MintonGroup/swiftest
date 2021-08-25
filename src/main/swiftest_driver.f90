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
   class(swiftest_parameters),   allocatable  :: param            !! Run configuration parameters
   integer(I4B)                               :: integrator       !! Integrator type code (see swiftest_globals for symbolic names)
   character(len=:),allocatable               :: param_file_name  !! Name of the file containing user-defined parameters
   integer(I4B)                               :: ierr             !! I/O error code 
   integer(I8B)                               :: iloop            !! Loop counter
   integer(I8B)                               :: idump            !! Dump cadence counter
   integer(I8B)                               :: iout             !! Output cadence counter
   integer(I8B)                               :: nloops           !! Number of steps to take in the simulation
   real(DP)                                   :: old_t_final = 0.0_DP !! Output time at which writing should start, in order to prevent duplicate lines being written for restarts

   ierr = io_get_args(integrator, param_file_name)
   if (ierr /= 0) then
      write(*,*) 'Error reading in arguments from the command line'
      call util_exit(FAILURE)
   end if
   !> Read in the user-defined parameters file and the initial conditions of the system
   select case(integrator)
   case(symba)
      allocate(symba_parameters :: param)
   case default
      allocate(swiftest_parameters :: param)
   end select
   param%integrator = integrator

   call setup_construct_system(nbody_system, param)
   call param%read_in(param_file_name)

   associate(t          => param%t, &
             t0         => param%t0, &
             dt         => param%dt, &
             tstop      => param%tstop, &
             istep_out  => param%istep_out, &
             istep_dump => param%istep_dump)  

      call nbody_system%initialize(param)
      t = t0
      iloop = 0
      iout = istep_out
      idump = istep_dump
      nloops = ceiling(tstop / dt, kind=I8B)
      ! Prevent duplicate frames from being written if this is a restarted run
      if (param%lrestart) then
         old_t_final = nbody_system%get_old_t_final(param)
      else
         old_t_final = t0
         if (istep_out > 0) call nbody_system%write_frame(param)
         call nbody_system%dump(param)
      end if


      !> Define the maximum number of threads
      nthreads = 1            ! In the *serial* case
      !$ nthreads = omp_get_max_threads() ! In the *parallel* case
      !$ write(*,'(a)')   ' OpenMP parameters:'
      !$ write(*,'(a)')   ' ------------------'
      !$ write(*,'(a,i3,/)') ' Number of threads  = ', nthreads 
      write(*, *) " *************** Main Loop *************** "
      do iloop = 1, nloops
         !> Step the system forward in time
         call nbody_system%step(param, t, dt)

         t = t0 + iloop * dt

         !> Evaluate any discards or collisional outcomes
         call nbody_system%discard(param)

         !> If the loop counter is at the output cadence value, append the data file with a single frame
         if (istep_out > 0) then
            iout = iout - 1
            if (iout == 0) then
               if (t > old_t_final) call nbody_system%write_frame(param)
               iout = istep_out
            end if
         end if

         !> If the loop counter is at the dump cadence value, dump the state of the system to a file in case it needs to be restarted
         if (istep_dump > 0) then
            idump = idump - 1
            if (idump == 0) then
               call nbody_system%dump(param)
               idump = istep_dump
            end if
         end if
      end do
   end associate

   call util_exit(SUCCESS)

   stop
end program swiftest_driver
