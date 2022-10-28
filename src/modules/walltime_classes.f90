!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module walltime_classes
   !! author: David A. Minton
   !!
   !! Classes and methods used to compute elasped wall time
   use swiftest_globals
   use swiftest_classes, only : swiftest_parameters, swiftest_pl
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
      real(DP)     :: wall_step                  !! Value of the step elapsed time
      logical      :: main_is_started = .false. !! Logical flag indicating whether or not the main timer has been reset or not
      logical      :: is_paused = .false. !! Logical flag indicating whether or not the timer is paused

   contains
      procedure :: reset       => walltime_reset      !! Resets the clock ticker, settting main_start to the current ticker value
      procedure :: start       => walltime_start      !! Starts or resumes the step timer
      procedure :: start_main  => walltime_start_main !! Starts the main timer
      procedure :: stop        => walltime_stop       !! Pauses the step timer
      procedure :: report      => walltime_report     !! Prints the elapsed time information to the terminal
   end type walltimer

   type, extends(walltimer) :: interaction_timer
      character(len=STRMAX) :: loopname              !! Stores the name of the loop being timed for logging purposes
      character(len=NAMELEN) :: looptype             !! Stores the type of loop (e.g. INTERACTION or ENCOUNTER)
      integer(I8B) :: max_interactions = huge(1_I8B) !! Stores the number of pl-pl interactions that failed when attempting to flatten (e.g. out of memory). Adapting won't occur if ninteractions > max_interactions
      integer(I8B) :: last_interactions = 0          !! Number of interactions that were computed last time. The timer is only run if there has been a change to the number of interactions
      integer(I4B) :: step_counter = 0               !! Number of steps that have elapsed since the last timed loop
      logical      :: is_on = .false.                !! The loop timer is currently active
      integer(I4B) :: stage = 1                      !! The stage of the loop timing (1 or 2)
      logical      :: stage1_is_advanced            !! Logical flag indicating whether stage1 was done with a flat loop (.true.) or triangular loop (.false.)
      integer(I8B) :: stage1_ninteractions           !! Number of interactions computed during stage 1
      real(DP)     :: stage1_metric                  !! Metric used to judge the performance of a timed loop (e.g. (count_finish_step - count_start_step) / ninteractions)
      real(DP)     :: stage2_metric                  !! Metric used to judge the performance of a timed loop (e.g. (count_finish_step - count_start_step) / ninteractions)
   contains
      procedure :: adapt => walltime_interaction_adapt !! Runs the interaction loop adaptation algorithm on an interaction loop
      procedure :: check => walltime_interaction_check !! Checks whether or not the loop should be timed and starts the timer if the conditions for starting are met
      procedure :: flip  => walltime_interaction_flip_loop_style  !! Flips the interaction loop style from FLAT to TRIANGULAR or vice vers
      procedure :: time_this_loop => walltime_interaction_time_this_loop !! Starts the interaction loop timer
   end type interaction_timer

   interface
      module subroutine walltime_report(self, message, nsubsteps)
         implicit none
         class(walltimer),           intent(inout) :: self      !! Walltimer object
         character(len=*),           intent(in)    :: message   !! Message to prepend to the wall time terminal output
         integer(I4B), optional,     intent(in)    :: nsubsteps !! Number of substeps used to compute the time per step 
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

   interface
      module subroutine walltime_interaction_adapt(self, param, ninteractions, pl)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(interaction_timer),   intent(inout) :: self          !! Interaction loop timer object
         class(swiftest_parameters), intent(inout) :: param         !! Current run configuration parameters
         integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop and to determine if number of interactions has changed since the last timing
         class(swiftest_pl),         intent(inout), optional :: pl            !! Swiftest massive body object
      end subroutine walltime_interaction_adapt

      module function walltime_interaction_check(self, param, ninteractions) result(ltimeit)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(interaction_timer),   intent(inout) :: self    !! Interaction loop timer object
         class(swiftest_parameters), intent(inout) :: param   !! Current run configuration parameters
         integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop and to determine if number of interactions has changed since the last timing
         logical                                   :: ltimeit !! Logical flag indicating whether this loop should be timed or not
      end function walltime_interaction_check

      module subroutine walltime_interaction_flip_loop_style(self, param, pl)
         use swiftest_classes, only : swiftest_parameters, swiftest_pl
         implicit none
         class(interaction_timer),   intent(inout) :: self  !! Interaction loop timer object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
         class(swiftest_pl),         intent(inout), optional :: pl    !! Swiftest massive body object
      end subroutine walltime_interaction_flip_loop_style

      module subroutine walltime_interaction_time_this_loop(self, param, ninteractions, pl)
         use swiftest_classes, only : swiftest_parameters, swiftest_pl
         implicit none
         class(interaction_timer),   intent(inout) :: self          !! Interaction loop timer object
         class(swiftest_parameters), intent(inout) :: param         !! Current run configuration parameters
         integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop)
         class(swiftest_pl),         intent(inout), optional :: pl            !! Swiftest massive body object
      end subroutine walltime_interaction_time_this_loop

   end interface



end module walltime_classes