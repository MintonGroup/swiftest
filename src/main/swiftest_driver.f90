!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

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
   character(len=:), allocatable              :: display_style    !! Style of the output display {"STANDARD", "COMPACT", "PROGRESS"}). Default is "STANDARD"
   integer(I4B)                               :: ierr             !! I/O error code 
   integer(I8B)                               :: iloop            !! Loop counter
   integer(I8B)                               :: idump            !! Dump cadence counter
   integer(I8B)                               :: iout             !! Output cadence counter
   integer(I8B)                               :: ioutput_t0       !! The output frame counter at time 0
   integer(I8B)                               :: nloops           !! Number of steps to take in the simulation
   real(DP)                                   :: old_t_final = 0.0_DP !! Output time at which writing should start, in order to prevent duplicate lines being written for restarts
   type(walltimer)                            :: integration_timer !! Object used for computing elapsed wall time
   real(DP)                                   :: tfrac
   type(progress_bar)                         :: pbar              !! Object used to print out a progress bar
   character(*), parameter                    :: statusfmt   = '("Time = ", ES12.5, "; fraction done = ", F6.3, ' // & 
                                                                '"; Number of active pl, tp = ", I5, ", ", I5)'
   character(*), parameter                    :: symbastatfmt   = '("Time = ", ES12.5, "; fraction done = ", F6.3, ' // &
                                                                   '"; Number of active plm, pl, tp = ", I5, ", ", I5, ", ", I5)'
   character(*), parameter                    :: pbarfmt = '("Time = ", ES12.5," of ",ES12.5)'
   character(len=64)                          :: pbarmessage


   call io_get_args(integrator, param_file_name, display_style)

   !> Read in the user-defined parameters file and the initial conditions of the system
   select case(integrator)
   case(symba)
      allocate(symba_parameters :: param)
   case default
      allocate(swiftest_parameters :: param)
   end select
   param%integrator = integrator
   call param%set_display(display_style)

   !> Define the maximum number of threads
   nthreads = 1            ! In the *serial* case
   !$ nthreads = omp_get_max_threads() ! In the *parallel* case
   !$ write(param%display_unit,'(a)')   ' OpenMP parameters:'
   !$ write(param%display_unit,'(a)')   ' ------------------'
   !$ write(param%display_unit,'(a,i3,/)') ' Number of threads = ', nthreads 
   !$ if (param%log_output) write(*,'(a,i3)') ' OpenMP: Number of threads = ',nthreads

   call setup_construct_system(nbody_system, param)
   call param%read_in(param_file_name)

   associate(t               => param%t, &
             t0              => param%t0, &
             dt              => param%dt, &
             tstop           => param%tstop, &
             istep_out       => param%istep_out, &
             istep_dump      => param%istep_dump, &
             ioutput         => param%ioutput, &
             display_style   => param%display_style, &
             display_unit    => param%display_unit)

      call nbody_system%initialize(param)
      t = t0
      iloop = 0
      iout = istep_out
      idump = istep_dump
      nloops = ceiling((tstop - t0) / dt, kind=I8B)
      ioutput_t0 = int(t0 / dt / istep_out, kind=I8B)
      ioutput = ioutput_t0
      ! Prevent duplicate frames from being written if this is a restarted run
      if ((param%lrestart) .and. ((param%out_type == REAL8_TYPE) .or. param%out_type == REAL4_TYPE)) then
         old_t_final = nbody_system%get_old_t_final_bin(param)
      else if ((param%lrestart) .and. ((param%out_type == NETCDF_DOUBLE_TYPE) .or. param%out_type == NETCDF_FLOAT_TYPE)) then
         old_t_final = nbody_system%get_old_t_final_netcdf(param)
      else
         old_t_final = t0
         if (param%lenergy) call nbody_system%conservation_report(param, lterminal=.false.) ! This will save the initial values of energy and momentum
         if (istep_out > 0) call nbody_system%write_frame(param)
      end if


      write(display_unit, *) " *************** Main Loop *************** "
      if (param%lrestart .and. param%lenergy) call nbody_system%conservation_report(param, lterminal=.true.)
      if (display_style == "PROGRESS") then
         call pbar%reset(nloops)
         write(pbarmessage,fmt=pbarfmt) t0, tstop
         call pbar%update(1,message=pbarmessage)
      end if
      do iloop = 1, nloops
         !> Step the system forward in time
         call integration_timer%start()
         call nbody_system%step(param, t, dt)
         call integration_timer%stop()

         t = t0 + iloop * dt

         !> Evaluate any discards or collisional outcomes
         call nbody_system%discard(param)
         if (display_style == "PROGRESS") call pbar%update(iloop)

         !> If the loop counter is at the output cadence value, append the data file with a single frame
         if (istep_out > 0) then
            iout = iout - 1
            if (iout == 0) then
               ioutput = ioutput_t0 + iloop / istep_out
               if (t > old_t_final) call nbody_system%write_frame(param)

               tfrac = (param%t - param%t0) / (param%tstop - param%t0)

               select type(pl => nbody_system%pl)
               class is (symba_pl)
                  write(display_unit, symbastatfmt) param%t, tfrac, pl%nplm, pl%nbody, nbody_system%tp%nbody
               class default
                  write(display_unit, statusfmt) param%t, tfrac, pl%nbody, nbody_system%tp%nbody
               end select
               if (param%lenergy) call nbody_system%conservation_report(param, lterminal=.true.)
               call integration_timer%report(message="Integration steps:", unit=display_unit, nsubsteps=istep_out)
               call integration_timer%reset()

               iout = istep_out
               if (display_style == "PROGRESS") then
                  write(pbarmessage,fmt=pbarfmt) t, tstop
                  call pbar%update(1,message=pbarmessage)
               end if
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

   call nbody_system%dealloc()

   call util_exit(SUCCESS)

   stop
end program swiftest_driver
