! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

module walltime
   !! author: David A. Minton
   !!
   !! Classes and methods used to compute elasped wall time
   use globals
   use base
   implicit none
   public

   integer(I4B) :: INTERACTION_TIMER_CADENCE = 1000 !! Minimum number of steps to wait before timing an interaction loop in ADAPTIVE mode
   character(len=*), parameter :: INTERACTION_TIMER_LOG_OUT  = "interaction_timer.log" !! Name of log file for recording results of interaction loop timing
   character(len=*), parameter :: ENCOUNTER_PLPL_TIMER_LOG_OUT  = "encounter_check_plpl_timer.log" !! Name of log file for recording results of encounter check method timing
   character(len=*), parameter :: ENCOUNTER_PLTP_TIMER_LOG_OUT  = "encounter_check_pltp_timer.log" !! Name of log file for recording results of encounter check method timing

   type :: walltimer
      integer(I8B) :: count_rate                 !! Rate at wich the clock ticks
      integer(I8B) :: count_max                  !! Maximum value of the clock ticker
      integer(I8B) :: count_start_main           !! Value of the clock ticker at when the timer is first called
      integer(I8B) :: count_start_step           !! Value of the clock ticker at the start of a timed step
      integer(I8B) :: count_stop_step            !! Value of the clock ticker at the end of a timed step
      integer(I8B) :: count_pause                !! Value of the clock ticker at the end of a timed step
      integer(I4B) :: nsubsteps                  !! Number of substeps in an interval (number of times the timer is turned off and back on again)
      real(DP)     :: wall_step                  !! Value of the step elapsed time
      real(DP)     :: wall_main                  !! Value of the main clock elapsed time
      real(DP)     :: wall_per_substep           !! Value of time per substep 
      logical      :: main_is_started = .false. !! Logical flag indicating whether or not the main timer has been reset or not
      logical      :: is_paused = .false. !! Logical flag indicating whether or not the timer is paused

   contains
      procedure :: reset       => walltime_reset      !! Resets the clock ticker, settting main_start to the current ticker value
      procedure :: start       => walltime_start      !! Starts or resumes the step timer
      procedure :: start_main  => walltime_start_main !! Starts the main timer
      procedure :: stop        => walltime_stop       !! Pauses the step timer
      procedure :: report      => walltime_report     !! Prints the elapsed time information to the terminal
   end type walltimer


   interface
      module subroutine walltime_report(self, message, unit)
         implicit none
         class(walltimer),           intent(inout) :: self      !! Walltimer object
         character(len=*),           intent(in)    :: message   !! Message to prepend to the wall time terminal output
         integer(I4B),               intent(in)    :: unit      !! Output file unit for report text to be directed
      end subroutine walltime_report

      module subroutine walltime_reset(self)
         implicit none
         class(walltimer),           intent(inout) :: self  !! Walltimer object
      end subroutine walltime_reset 

      module subroutine walltime_start(self)
         implicit none
         class(walltimer),           intent(inout) :: self  !! Walltimer object
      end subroutine walltime_start

      module subroutine walltime_start_main(self)
         implicit none
         class(walltimer),           intent(inout) :: self  !! Walltimer object
      end subroutine walltime_start_main

      module subroutine walltime_stop(self)
         implicit none
         class(walltimer),           intent(inout) :: self  !! Walltimer object
      end subroutine walltime_stop
   end interface

end module walltime