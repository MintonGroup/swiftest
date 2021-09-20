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

   type :: walltimer
      integer(I8B) :: count_rate                 !! Rate at wich the clock ticks
      integer(I8B) :: count_max                  !! Maximum value of the clock ticker
      integer(I8B) :: count_start_main           !! Value of the clock ticker at when the timer is first called
      integer(I8B) :: count_start_step           !! Value of the clock ticker at the start of a timed step
      integer(I8B) :: count_finish_step          !! Value of the clock ticker at the end of a timed step
      logical      :: lmain_is_started = .false. !! Logical flag indicating whether or not the main timer has been reset or not
   contains
      procedure :: reset  => walltime_reset  !! Resets the clock ticker, settting main_start to the current ticker value
      procedure :: start  => walltime_start  !! Starts the timer, setting step_start to the current ticker value
      procedure :: finish => walltime_finish !! Ends the timer, setting step_finish to the current ticker value and printing the elapsed time information to the terminal
   end type walltimer

   type, extends(walltimer) :: interaction_timer
      integer(I8B) :: max_interactions = huge(1_I8B) !! Stores the number of pl-pl interactions that failed when attempting to flatten (e.g. out of memory). Adapting won't occur if ninteractions > max_interactions
      integer(I8B) :: last_interactions = 0          !! Number of interactions that were computed last time. The timer is only run if there has been a change to the number of interactions
      integer(I4B) :: step_counter = 0               !! Number of steps that have elapsed since the last timed loop
      logical      :: is_on = .false.                !! The loop timer is currently active
      integer(I4B) :: stage = 1                      !! The stage of the loop timing (1 or 2)
      logical      :: stage1_is_flattened            !! Logical flag indicating whether stage1 was done with a flat loop (.true.) or triangular loop (.false.)
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
      module subroutine walltime_finish(self, nsubsteps, message, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(walltimer),           intent(inout) :: self      !! Walltimer object
         integer(I4B),               intent(in)    :: nsubsteps !! Number of substeps used to compute the time per step 
         character(len=*),           intent(in)    :: message   !! Message to prepend to the wall time terminal output
         class(swiftest_parameters), intent(inout) :: param     !! Current run configuration parameters
      end subroutine walltime_finish

      module subroutine walltime_reset(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(walltimer),           intent(inout) :: self  !! Walltimer object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine walltime_reset 

      module subroutine walltime_start(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(walltimer),           intent(inout) :: self  !! Walltimer object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine walltime_start
   end interface

   interface
      module subroutine walltime_interaction_adapt(self, param, pl, ninteractions)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(interaction_timer),   intent(inout) :: self          !! Interaction loop timer object
         class(swiftest_parameters), intent(inout) :: param         !! Current run configuration parameters
         class(swiftest_pl),         intent(inout) :: pl            !! Swiftest massive body object
         integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop and to determine if number of interactions has changed since the last timing
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
         class(swiftest_pl),         intent(inout) :: pl    !! Swiftest massive body object
      end subroutine walltime_interaction_flip_loop_style

      module subroutine walltime_interaction_time_this_loop(self, param, pl, ninteractions)
         use swiftest_classes, only : swiftest_parameters, swiftest_pl
         implicit none
         class(interaction_timer),   intent(inout) :: self          !! Interaction loop timer object
         class(swiftest_parameters), intent(inout) :: param         !! Current run configuration parameters
         class(swiftest_pl),         intent(inout) :: pl            !! Swiftest massive body object
         integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop)
      end subroutine walltime_interaction_time_this_loop

   end interface

   contains

      module subroutine walltime_finish(self, nsubsteps, message, param)
         !! author: David A. Minton
         !!
         !! Ends the timer, setting step_finish to the current ticker value and printing the elapsed time information to the terminal
         implicit none
         ! Arguments
         class(walltimer),           intent(inout) :: self      !! Walltimer object
         integer(I4B),               intent(in)    :: nsubsteps !! Number of substeps used to compute the time per step 
         character(len=*),           intent(in)    :: message   !! Message to prepend to the wall time terminal output
         class(swiftest_parameters), intent(inout) :: param     !! Current run configuration parameters
         ! Internals
         character(len=*), parameter     :: walltimefmt = '" Wall time (s): ", es12.5, "; Wall time/step in this interval (s):  ", es12.5'
         character(len=STRMAX)           :: fmt
         integer(I8B)                    :: count_delta_step, count_delta_main
         real(DP)                        :: wall_main         !! Value of total elapsed time at the end of a timed step
         real(DP)                        :: wall_step         !! Value of elapsed time since the start of a timed step
         real(DP)                        :: wall_per_substep  !! Value of time per substep 

         if (.not.self%lmain_is_started) then
            write(*,*) "Wall timer error: The step finish time cannot be calculated because the timer is not started!"
            return
         end if

         call system_clock(self%count_finish_step)

         count_delta_step = self%count_finish_step - self%count_start_step
         count_delta_main = self%count_finish_step - self%count_start_main
         wall_step = count_delta_step / (self%count_rate * 1.0_DP)
         wall_main = count_delta_main / (self%count_rate * 1.0_DP)
         wall_per_substep = wall_step / nsubsteps

         fmt = '("' //  adjustl(message) // '",' // walltimefmt // ')'
         write(*,trim(adjustl(fmt))) wall_main, wall_per_substep

         call self%start(param)

         return
      end subroutine walltime_finish


      module subroutine walltime_reset(self, param)
         !! author: David A. Minton
         !!
         !! Resets the clock ticker, settting main_start to the current ticker value
         implicit none
         ! Arguments
         class(walltimer),           intent(inout) :: self  !! Walltimer object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
         
         call system_clock(self%count_start_main, self%count_rate, self%count_max)
         self%lmain_is_started = .true.
         call self%start(param)

         return
      end subroutine walltime_reset 


      module subroutine walltime_start(self, param)
         !! author: David A. Minton
         !!
         !! Starts the timer, setting step_start to the current ticker value
         implicit none
         ! Arguments
         class(walltimer),           intent(inout) :: self  !! Walltimer object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters

         if (.not.self%lmain_is_started) then
            write(*,*) "Wall timer error: Cannot start the step time until reset is called at least once!"
            return
         end if

         call system_clock(self%count_start_step)

         return 
      end subroutine walltime_start


      module subroutine walltime_interaction_adapt(self, param, pl, ninteractions)
         !! author: David A. Minton
         !!
         !! Determines which of the two loop styles is fastest and keeps that one
         implicit none
         class(interaction_timer),   intent(inout) :: self          !! Walltimer object
         class(swiftest_parameters), intent(inout) :: param         !! Current run configuration parameters
         class(swiftest_pl),         intent(inout) :: pl            !! Swiftest massive body object
         integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop and to determine if number of interactions has changed since the last timing

         ! Record the elapsed time 
         call system_clock(self%count_finish_step)

         select case(self%stage)
         case(1)
            self%stage1_metric = (self%count_finish_step - self%count_start_step) / real(ninteractions, kind=DP)
         case(2)
            self%stage2_metric = (self%count_finish_step - self%count_start_step) / real(ninteractions, kind=DP)
            self%is_on = .false.
            self%step_counter = 0
            if (self%stage1_metric < self%stage2_metric) call self%flip(param, pl)  ! Go back to the original style, otherwise keep the stage2 style
         end select

         return
      end subroutine walltime_interaction_adapt


      module function walltime_interaction_check(self, param, ninteractions) result(ltimeit)
         !! author: David A. Minton
         !!
         !! Checks whether or not the loop should be timed and starts the timer if the conditions for starting are met
         implicit none
         ! Arguments
         class(interaction_timer),   intent(inout) :: self          !! Walltimer object
         class(swiftest_parameters), intent(inout) :: param         !! Current run configuration parameters
         integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop and to determine if number of interactions has changed since the last timing
         logical                                   :: ltimeit !! Logical flag indicating whether this loop should be timed or not
         ! Internals
         character(len=STRMAX) :: tstring

         if (self%is_on) then ! Entering the second stage of the loop timing. Therefore we will swap the interaction style and time this loop
            self%stage = self%stage + 1
            ltimeit = (self%stage == 2)
         else
            self%step_counter = min(self%step_counter + 1, INTERACTION_TIMER_CADENCE)
            ltimeit = .false.
            if (self%step_counter == INTERACTION_TIMER_CADENCE) then
               ltimeit = (ninteractions /= self%last_interactions)
               if (ltimeit) self%stage = 1
            end if
         end if
         self%is_on = ltimeit

         return
      end function walltime_interaction_check


      module subroutine walltime_interaction_flip_loop_style(self, param, pl)
         !! author: David A. Minton
         !!
         !! Flips the interaction loop style from FLAT to TRIANGULAR or vice versa
         implicit none
         ! Arguments
         class(interaction_timer),   intent(inout) :: self  !! Interaction loop timer object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
         class(swiftest_pl),         intent(inout) :: pl    !! Swiftest massive body object 

         param%lflatten_interactions = .not. param%lflatten_interactions
         if (param%lflatten_interactions) then
            call pl%flatten(param)
         else
            if (allocated(pl%k_plpl)) deallocate(pl%k_plpl)
         end if

         return
      end subroutine walltime_interaction_flip_loop_style


      module subroutine walltime_interaction_time_this_loop(self, param, pl, ninteractions)
         !! author: David A. Minton
         !!
         !! Resets the interaction loop timer, and saves the current value of the array flatten parameter
         implicit none
         ! Arguments
         class(interaction_timer),   intent(inout) :: self          !! Interaction loop timer object
         class(swiftest_parameters), intent(inout) :: param         !! Current run configuration parameters
         class(swiftest_pl),         intent(inout) :: pl            !! Swiftest massive body object
         integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop)

         self%is_on = .true.
         select case(self%stage)
         case(1)
            self%stage1_ninteractions = ninteractions      
            self%stage1_is_flattened = param%lflatten_interactions
         case(2)
            param%lflatten_interactions = self%stage1_is_flattened
            call self%flip(param, pl) 
         case default
            self%stage = 1
         end select
         call self%reset(param)

         return
      end subroutine walltime_interaction_time_this_loop


end module walltime_classes