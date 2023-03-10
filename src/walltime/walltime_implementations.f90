!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(walltime) s_walltime
   use swiftest
contains

   module subroutine walltime_stop(self)
      !! author: David A. Minton
      !!
      !! Pauses the step timer (but not the main timer). 
      implicit none
      ! Arguments
      class(walltimer),           intent(inout) :: self  !! Walltimer object
      ! Internals
      integer(I8B) :: count_delta

      if (self%is_paused) then
         write(*,*) "Wall timer error: Timer is already paused!"
         return
      end if

      call system_clock(self%count_pause)
      self%is_paused = .true.

      self%count_stop_step = self%count_pause

      count_delta = self%count_stop_step - self%count_start_step
      self%wall_step = count_delta / (self%count_rate * 1.0_DP)

      return 
   end subroutine walltime_stop


   module subroutine walltime_report(self, message, unit)
      !! author: David A. Minton
      !!
      !! Prints the elapsed time information to the terminal
      implicit none
      ! Arguments
      class(walltimer),           intent(inout) :: self      !! Walltimer object
      character(len=*),           intent(in)    :: message   !! Message to prepend to the wall time terminal output
      integer(I4B),               intent(in)    :: unit      !! Output file unit for report text to be directed
      ! Internals
      character(len=*), parameter     :: substepfmt   = '" Total wall time: ", es12.5, "; Interval wall time: ", es12.5, ";' //&
                                                        ' Interval wall time/step:  ", es12.5'
      character(len=STRMAX)           :: fmt
      integer(I8B)                    :: count_delta_step, count_delta_main, count_now

      if (.not.self%main_is_started) then
         write(*,*) "Wall timer error: The step finish time cannot be calculated because the timer is not started!"
         return
      end if

      call system_clock(count_now)
      count_delta_main = count_now - self%count_start_main
      self%wall_main = count_delta_main / (self%count_rate * 1.0_DP)

      count_delta_step = self%count_stop_step - self%count_start_step
      self%wall_step = count_delta_step / (self%count_rate * 1.0_DP)
      self%wall_per_substep = self%wall_step / self%nsubsteps

      fmt = '("' //  adjustl(message) // '",' // substepfmt // ')'
      write(unit,trim(adjustl(fmt))) self%wall_main, self%wall_step, self%wall_per_substep

      return
   end subroutine walltime_report


   module subroutine walltime_reset(self)
      !! author: David A. Minton
      !!
      !! Resets the step timer
      implicit none
      ! Arguments
      class(walltimer),           intent(inout) :: self  !! Walltimer object
      ! Internals

      self%is_paused = .false.
      self%wall_step = 0.0_DP
      self%wall_per_substep = 0.0_DP
      self%nsubsteps = 0

      return 
   end subroutine walltime_reset


   module subroutine walltime_start_main(self)
      !! author: David A. Minton
      !!
      !! Resets the clock ticker, settting main_start to the current ticker value
      implicit none
      ! Arguments
      class(walltimer),           intent(inout) :: self  !! Walltimer object
      
      call system_clock(self%count_start_main, self%count_rate, self%count_max)
      self%main_is_started = .true.
      self%wall_main = 0.0_DP

      return
   end subroutine walltime_start_main


   module subroutine walltime_start(self)
      !! author: David A. Minton
      !!
      !! Starts or resumes the step timer
      !!
      implicit none
      ! Arguments
      class(walltimer),           intent(inout) :: self  !! Walltimer object
      ! Internals
      integer(I8B) :: count_resume, count_delta

      if (.not.self%main_is_started) then
         call self%reset()
         call self%start_main()
      end if

      if (self%is_paused) then ! Resume a paused step timer
         call system_clock(count_resume)
         count_delta = count_resume - self%count_pause 
         self%count_pause = 0_I8B 
         self%count_start_step = self%count_start_step + count_delta
         self%is_paused = .false.
      else ! Start a new step timer
         call system_clock(self%count_start_step)
      end if
      self%nsubsteps = self%nsubsteps + 1

      return 
   end subroutine walltime_start

end submodule s_walltime