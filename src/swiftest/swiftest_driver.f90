!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 
submodule(swiftest) s_swiftest_driver
contains

   module subroutine swiftest_driver(integrator, param_file_name, display_style)

      !! author: David A. Minton
      !!
      !! Driver program for the Swiftest integrators. Unlike the earlier Swift and Swifter drivers, in Swiftest all integrators 
      !!    are run from this single program. 
      !!
      !! Adapted from Swifter by David E. Kaufmann's Swifter driver programs swifter_[bs,helio,ra15,rmvs,symba,tu4,whm].f90
      !! Adapted from Hal Levison and Martin Duncan's Swift driver programs
      implicit none

      ! Arguments
      character(len=:), intent(in), allocatable :: integrator      !! Symbolic code of the requested integrator  
      character(len=:), intent(in), allocatable :: param_file_name !! Name of the input parameters file
      character(len=:), intent(in), allocatable :: display_style   !! Style of the output display {"STANDARD", "COMPACT", "PROGRESS"}). Default is "STANDARD"   

      ! Internals
      class(swiftest_nbody_system), allocatable :: nbody_system      !! Polymorphic object containing the nbody system to be integrated
      type(swiftest_parameters)                 :: param             !! Run configuration parameters
      class(swiftest_storage),      allocatable :: system_history    !! Stores the system history between output dumps
      type(walltimer)                           :: integration_timer !! Object used for computing elapsed wall time

      !> Read in the user-defined parameters file and the initial conditions of the nbody_system
      param%integrator = trim(adjustl(integrator))
      param%display_style = trim(adjustl(display_style))
      call param%read_in(param_file_name)
#ifdef COARRAY
      if (.not.param%lcoarray .and. (this_image() /= 1)) stop ! Single image mode
#endif

      associate(t0       => param%t0, &
         tstart          => param%tstart, &
         dt              => param%dt, &
         tstop           => param%tstop, &
         iloop           => param%iloop, &
         istart          => param%istart, &
         iout            => param%iout, &
         idump           => param%idump, &
         nout            => param%nout, &
         istep           => param%istep, &
         nloops          => param%nloops, &
         istep_out       => param%istep_out, &
         fstep_out       => param%fstep_out, &
         ltstretch       => param%ltstretch, &
         dump_cadence    => param%dump_cadence, &
         display_unit    => param%display_unit)

         ! Set up loop and output cadence variables
         nloops = ceiling((tstop - t0) / dt, kind=I8B)
         istart =  ceiling((tstart - t0) / dt + 1.0_DP, kind=I8B)
         iloop = istart - 1
         iout = 0
         idump = 0
         if (ltstretch) then
            nout = floor(log(1.0_DP + iloop * (fstep_out - 1.0_DP) / istep_out) / log(fstep_out)) 
            istep = floor(istep_out * fstep_out**nout, kind=I4B)
         else
            nout = 1
            istep = istep_out
         end if

         ! Set up nbody_system storage for intermittent file dumps
         if (dump_cadence == 0) dump_cadence = int(ceiling(nloops / (1.0_DP * istep_out), kind=I8B), kind=I4B)

         ! Construct the main n-body nbody_system using the user-input integrator to choose the type of nbody_system
         call swiftest_util_setup_construct_system(nbody_system, param)

         !> Define the maximum number of threads
         nthreads = 1            ! In the *serial* case
         !$ nthreads = omp_get_max_threads() ! In the *parallel* case
#ifdef COARRAY
         if (this_image() == 1 .or. param%log_output) then
#endif 
            !$ write(param%display_unit,'(a)')   ' OpenMP parameters:'
            !$ write(param%display_unit,'(a)')   ' ------------------'
            !$ write(param%display_unit,'(a,i3,/)') ' Number of threads = ', nthreads 
            !$ if (param%log_output) write(*,'(a,i3)') ' OpenMP: Number of threads = ',nthreads
#ifdef COARRAY
            if (param%lcoarray) then
               write(param%display_unit,*)   ' Coarray parameters:'
               write(param%display_unit,*)   ' -------------------'
               write(param%display_unit,*) ' Number of images = ', num_images()
               if (param%log_output .and. this_image() == 1) write(*,'(a,i3)') ' Coarray: Number of images = ',num_images()
            else
               write(param%display_unit,*)   ' Coarrays disabled.'
               if (param%log_output) write(*,*)   ' Coarrays disabled.'
            end if
         end if
#endif 
         if (param%log_output) flush(param%display_unit)

#ifdef COARRAY  
         ! The following line lets us read in the input files one image at a time. Letting each image read the input in is faster than broadcasting all of the data
         if (param%lcoarray .and. (this_image() /= 1)) sync images(this_image() - 1)
#endif 
         call nbody_system%initialize(system_history, param)
#ifdef COARRAY  
         if (param%lcoarray .and. (this_image() < num_images())) sync images(this_image() + 1)

         ! Distribute test particles to the various images
         if (param%lcoarray) call nbody_system%coarray_distribute(param)
#endif

         ! If this is a new run, compute energy initial conditions (if energy tracking is turned on) and write the initial conditions to file.
         call nbody_system%display_run_information(param, integration_timer, phase="first")

         if (param%lenergy) then
            if (param%lrestart) then
               call nbody_system%get_t0_values(system_history%nc, param)
            else
               call nbody_system%conservation_report(param, lterminal=.false.) ! This will save the initial values of energy and momentum
            end if
            call nbody_system%conservation_report(param, lterminal=.true.)
         end if

         call system_history%take_snapshot(param,nbody_system)
         call nbody_system%dump(param, system_history)

         do iloop = istart, nloops
            !> Step the nbody_system forward in time
            call integration_timer%start()
            call nbody_system%step(param, nbody_system%t, dt)
            call integration_timer%stop()

            nbody_system%t = t0 + iloop * dt

            !> Evaluate any discards or collisional outcomes
            call nbody_system%discard(param)

            !> If the loop counter is at the output cadence value, append the data file with a single frame
            if (istep_out > 0) then
               iout = iout + 1
               if ((iout == istep) .or. (iloop == nloops)) then
                  iout = 0
                  idump = idump + 1
                  if (ltstretch) then 
                     nout = nout + 1
                     istep = floor(istep_out * fstep_out**nout, kind=I4B)
                  end if

                  call system_history%take_snapshot(param,nbody_system)

                  if (idump == dump_cadence) then
                     idump = 0
                     call nbody_system%dump(param, system_history)
#ifdef COARRAY
                     if (param%lcoarray) call nbody_system%coarray_balance(param)
#endif
                  end if
#ifdef COARRAY
                  if (this_image() == 1 .or. param%log_output) then
#endif
                     call integration_timer%report(message="Integration steps:", unit=display_unit)
#ifdef COARRAY  
                  end if !(this_image() == 1)
#endif
                  call nbody_system%display_run_information(param, integration_timer)
                  call integration_timer%reset()
#ifdef COARRAY
                  if (this_image() == 1 .or. param%log_output) then
#endif
                     if (param%lenergy) call nbody_system%conservation_report(param, lterminal=.true.)
#ifdef COARRAY
                  end if ! (this_image() == 1)
#endif
               end if
            end if

         end do
         ! Dump any remaining history if it exists
         call nbody_system%dump(param, system_history)
         call nbody_system%display_run_information(param, integration_timer, phase="last")
      end associate
   return
   end subroutine swiftest_driver

end submodule s_swiftest_driver
