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
   use symba
   implicit none

   class(swiftest_nbody_system), allocatable      :: system           !! Polymorphic object containing the nbody system to be integrated
   class(swiftest_parameters),   allocatable      :: param             !! Run configuration parameters
   character(len=:), allocatable                  :: integrator        !! Integrator type code (see globals for symbolic names)
   character(len=:), allocatable                  :: param_file_name   !! Name of the file containing user-defined parameters
   character(len=:), allocatable                  :: display_style     !! Style of the output display {"STANDARD", "COMPACT", "PROGRESS"}). Default is "STANDARD"
   integer(I8B)                                   :: istart            !! Starting index for loop counter
   integer(I8B)                                   :: nloops            !! Number of steps to take in the simulation
   integer(I4B)                                   :: iout              !! Output cadence counter
   integer(I4B)                                   :: idump             !! Dump cadence counter
   type(walltimer)                                :: integration_timer !! Object used for computing elapsed wall time
   real(DP)                                       :: tfrac             !! Fraction of total simulation time completed
   type(progress_bar)                             :: pbar              !! Object used to print out a progress bar
   character(*), parameter                        :: statusfmt = '("Time = ", ES12.5, "; fraction done = ", F6.3, ' // & 
                                                                 '"; Number of active pl, tp = ", I6, ", ", I6)'
   character(*), parameter                        :: symbastatfmt = '("Time = ", ES12.5, "; fraction done = ", F6.3, ' // &
                                                                    '"; Number of active plm, pl, tp = ", I6, ", ", I6, ", ", I6)'
   character(*), parameter                        :: pbarfmt = '("Time = ", ES12.5," of ",ES12.5)'
   character(len=64)                              :: pbarmessage

   character(*), parameter                        :: symbacompactfmt = '(";NPLM",ES22.15,$)'
   !type(base_storage(nframes=:)), allocatable :: system_history

   call swiftest_io_get_args(integrator, param_file_name, display_style)

   !> Read in the user-defined parameters file and the initial conditions of the system
   allocate(swiftest_parameters :: param)
   param%integrator = trim(adjustl(integrator))
   call param%set_display(display_style)
   call param%read_in(param_file_name)


   associate(t0       => param%t0, &
      tstart          => param%tstart, &
      dt              => param%dt, &
      tstop           => param%tstop, &
      iloop           => param%iloop, &
      istep_out       => param%istep_out, &
      dump_cadence    => param%dump_cadence, &
      ioutput         => param%ioutput, &
      display_style   => param%display_style, &
      display_unit    => param%display_unit)


      ! Set up loop and output cadence variables
      nloops = ceiling((tstop - t0) / dt, kind=I8B)
      istart =  ceiling((tstart - t0) / dt + 1.0_DP, kind=I8B)
      ioutput = max(int(istart / istep_out, kind=I4B),1)

      ! Set up system storage for intermittent file dumps
      if (dump_cadence == 0) dump_cadence = ceiling(nloops / (1.0_DP * istep_out), kind=I8B)

      ! Construct the main n-body system using the user-input integrator to choose the type of system
      call swiftest_setup_construct_system(system, param)

      !> Define the maximum number of threads
      nthreads = 1            ! In the *serial* case
      !$ nthreads = omp_get_max_threads() ! In the *parallel* case
      !$ write(param%display_unit,'(a)')   ' OpenMP parameters:'
      !$ write(param%display_unit,'(a)')   ' ------------------'
      !$ write(param%display_unit,'(a,i3,/)') ' Number of threads = ', nthreads 
      !$ if (param%log_output) write(*,'(a,i3)') ' OpenMP: Number of threads = ',nthreads

      call system%initialize(param)

      associate (system_history => param%system_history)
         ! If this is a new run, compute energy initial conditions (if energy tracking is turned on) and write the initial conditions to file.
         if (param%lenergy) then
            if (param%lrestart) then
               call system%conservation_report(param, lterminal=.true.)
            else
               call system%conservation_report(param, lterminal=.false.) ! This will save the initial values of energy and momentum
            end if
         end if
         call system_history%take_snapshot(param,system)
         call system%dump(param)

         write(display_unit, *) " *************** Main Loop *************** "

         if (display_style == "PROGRESS") then
            call pbar%reset(nloops)
            write(pbarmessage,fmt=pbarfmt) t0, tstop
            call pbar%update(1,message=pbarmessage)
         else if (display_style == "COMPACT") then
            write(*,*) "SWIFTEST START " // param%integrator
            call system%compact_output(param,integration_timer)
         end if

         iout = 0
         idump = 0
         system%t = tstart
         do iloop = istart, nloops
            !> Step the system forward in time
            call integration_timer%start()
            call system%step(param, system%t, dt)
            call integration_timer%stop()

            system%t = t0 + iloop * dt

            !> Evaluate any discards or collisional outcomes
            call system%discard(param)
            if (display_style == "PROGRESS") call pbar%update(iloop)

            !> If the loop counter is at the output cadence value, append the data file with a single frame
            if (istep_out > 0) then
               iout = iout + 1
               if (iout == istep_out) then
                  iout = 0
                  idump = idump + 1
                  call system_history%take_snapshot(param,system)

                  if (idump == dump_cadence) then
                     idump = 0
                     call system%dump(param)

                  end if

                  tfrac = (system%t - t0) / (tstop - t0)

                  select type(pl => system%pl)
                  class is (symba_pl)
                     write(display_unit, symbastatfmt) system%t, tfrac, pl%nplm, pl%nbody, system%tp%nbody
                  class default
                     write(display_unit, statusfmt) system%t, tfrac, pl%nbody, system%tp%nbody
                  end select
                  if (param%lenergy) call system%conservation_report(param, lterminal=.true.)
                  call integration_timer%report(message="Integration steps:", unit=display_unit, nsubsteps=istep_out)

                  if (display_style == "PROGRESS") then
                     write(pbarmessage,fmt=pbarfmt) system%t, tstop
                     call pbar%update(1,message=pbarmessage)
                  else if (display_style == "COMPACT") then
                     call system%compact_output(param,integration_timer)
                  end if

                  call integration_timer%reset()

               end if
            end if

         end do
         ! Dump any remaining history if it exists
         call system%dump(param)
         call system_history%dump(param)
         if (display_style == "COMPACT") write(*,*) "SWIFTEST STOP" // param%integrator
      end associate
   end associate

   call util_exit(SUCCESS)
end program swiftest_driver
