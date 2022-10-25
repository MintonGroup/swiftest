submodule(walltime_classes) s_walltime
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


   module subroutine walltime_report(self, message, nsubsteps)
      !! author: David A. Minton
      !!
      !! Prints the elapsed time information to the terminal
      implicit none
      ! Arguments
      class(walltimer),           intent(inout) :: self      !! Walltimer object
      character(len=*),           intent(in)    :: message   !! Message to prepend to the wall time terminal output
      integer(I4B), optional,     intent(in)    :: nsubsteps !! Number of substeps used to compute the time per step 
      ! Internals
      character(len=*), parameter     :: nosubstepfmt = '" Total wall time: ", es12.5, "; Interval wall time: ", es12.5 ' 
      character(len=*), parameter     :: substepfmt   = '" Total wall time: ", es12.5, "; Interval wall time: ", es12.5, ";' //&
                                                        'Interval wall time/step:  ", es12.5'
      character(len=STRMAX)           :: fmt
      integer(I8B)                    :: count_delta_step, count_delta_main, count_now
      real(DP)                        :: wall_main         !! Value of total elapsed time at the end of a timed step
      real(DP)                        :: wall_step         !! Value of elapsed time since the start of a timed step
      real(DP)                        :: wall_per_substep  !! Value of time per substep 

      if (.not.self%main_is_started) then
         write(*,*) "Wall timer error: The step finish time cannot be calculated because the timer is not started!"
         return
      end if

      call system_clock(count_now)
      count_delta_main = count_now - self%count_start_main
      count_delta_step = count_now - self%count_start_step
      wall_main = count_delta_main / (self%count_rate * 1.0_DP)
      wall_step = count_delta_step / (self%count_rate * 1.0_DP)
      if (present(nsubsteps)) then
         wall_per_substep = wall_step / nsubsteps
         fmt = '("' //  adjustl(message) // '",' // substepfmt // ')'
         write(*,trim(adjustl(fmt))) wall_main, self%wall_step, wall_per_substep
      else
         fmt = '("' //  adjustl(message) // '",' // nosubstepfmt // ')'
         write(*,trim(adjustl(fmt))) wall_main, self%wall_step
      end if


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

      return 
   end subroutine walltime_start


   module subroutine walltime_interaction_adapt(self, param, ninteractions, pl)
      !! author: David A. Minton
      !!
      !! Determines which of the two loop styles is fastest and keeps that one
      implicit none
      ! Arguments
      class(interaction_timer),   intent(inout) :: self          !! Walltimer object
      class(swiftest_parameters), intent(inout) :: param         !! Current run configuration parameters
      integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop and to determine if number of interactions has changed since the last timing
      class(swiftest_pl),         intent(inout), optional :: pl            !! Swiftest massive body object
      ! Internals
      character(len=STRMAX) :: nstr, cstr, mstr
      character(len=11) :: lstyle, advancedstyle, standardstyle
      character(len=1) :: schar
      logical :: ladvanced_final
      character(len=NAMELEN)  :: logfile

      ! Record the elapsed time 
      call self%stop()

      select case(trim(adjustl(self%looptype)))
      case("INTERACTION")
         write(advancedstyle, *) "FLAT      "
         write(standardstyle, *) "TRIANGULAR"
         write(logfile,*) INTERACTION_TIMER_LOG_OUT
      case("ENCOUNTER_PLPL")
         write(advancedstyle, *) "SORTSWEEP "
         write(standardstyle, *) "TRIANGULAR"
         write(logfile,*) ENCOUNTER_PLPL_TIMER_LOG_OUT
      case("ENCOUNTER_PLTP")
         write(advancedstyle, *) "SORTSWEEP "
         write(standardstyle, *) "TRIANGULAR"
         write(logfile,*) ENCOUNTER_PLTP_TIMER_LOG_OUT
      case default
         write(logfile,*) "unknown_looptimer.log"
      end select

      write(schar,'(I1)') self%stage
      write(nstr,*) ninteractions

      select case(self%stage)
      case(1)
         if (self%stage1_is_advanced) then
            lstyle = advancedstyle
         else
            lstyle = standardstyle
         end if 
         self%stage1_metric = (self%count_stop_step - self%count_start_step) / real(ninteractions, kind=DP)
         write(mstr,*) self%stage1_metric
      case(2)
         if (.not.self%stage1_is_advanced) then
            lstyle = advancedstyle
         else
            lstyle = standardstyle
         end if 

         self%stage2_metric = (self%count_stop_step - self%count_start_step) / real(ninteractions, kind=DP)
         self%is_on = .false.
         self%step_counter = 0
         if (self%stage1_metric < self%stage2_metric) then
            ladvanced_final = self%stage1_is_advanced
            call self%flip(param, pl)  ! Go back to the original style, otherwise keep the stage2 style
         else
            ladvanced_final = .not.self%stage1_is_advanced
         end if
         write(mstr,*) self%stage2_metric
      end select

      write(cstr,*) self%count_stop_step - self%count_start_step

      call io_log_one_message(logfile, adjustl(lstyle) // " " // trim(adjustl(cstr)) // " " // &
                                       trim(adjustl(nstr)) // " " // trim(adjustl(mstr)))

      if (self%stage == 2) then
         if (ladvanced_final) then
            lstyle = advancedstyle
         else
            lstyle = standardstyle
         end if 
         call io_log_one_message(logfile, trim(adjustl(self%loopname)) // &
                                 ": the fastest loop method tested is " // trim(adjustl(lstyle)))
      end if

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
      class(interaction_timer),   intent(inout)           :: self  !! Interaction loop timer object
      class(swiftest_parameters), intent(inout)           :: param !! Current run configuration parameters
      class(swiftest_pl),         intent(inout), optional :: pl    !! Swiftest massive body object 

      select case(trim(adjustl(self%looptype)))
      case("INTERACTION")
         param%lflatten_interactions = .not. param%lflatten_interactions
      case("ENCOUNTER_PLPL")
         param%lencounter_sas_plpl= .not. param%lencounter_sas_plpl
      case("ENCOUNTER_PLTP")
         param%lencounter_sas_pltp= .not. param%lencounter_sas_pltp
      end select

      if (present(pl)) then
         if (param%lflatten_interactions) then
            call pl%flatten(param)
         else
            if (allocated(pl%k_plpl)) deallocate(pl%k_plpl)
         end if
      end if

      return
   end subroutine walltime_interaction_flip_loop_style


   module subroutine walltime_interaction_time_this_loop(self, param, ninteractions, pl)
      !! author: David A. Minton
      !!
      !! Resets the interaction loop timer, and saves the current value of the array flatten parameter
      implicit none
      ! Arguments
      class(interaction_timer),   intent(inout) :: self          !! Interaction loop timer object
      class(swiftest_parameters), intent(inout) :: param         !! Current run configuration parameters
      integer(I8B),               intent(in)    :: ninteractions !! Current number of interactions (used to normalize the timed loop)
      class(swiftest_pl),         intent(inout), optional :: pl            !! Swiftest massive body object
      ! Internals
      character(len=STRMAX) :: tstr
      character(len=1) :: schar
      character(len=NAMELEN) :: logfile

      select case(trim(adjustl(self%looptype)))
      case("INTERACTION")
         write(logfile,*) INTERACTION_TIMER_LOG_OUT
      case("ENCOUNTER_PLPL")
         write(logfile,*) ENCOUNTER_PLPL_TIMER_LOG_OUT
      case("ENCOUNTER_PLTP")
         write(logfile,*) ENCOUNTER_PLTP_TIMER_LOG_OUT
      case default
         write(logfile,*) "unknown_looptimer.log"
      end select

      self%is_on = .true.
      write(tstr,*) param%t
      select case(self%stage)
      case(1)
         self%stage1_ninteractions = ninteractions 
         select case(trim(adjustl(self%looptype)))
         case("INTERACTION")
            self%stage1_is_advanced = param%lflatten_interactions
         case("ENCOUNTER_PLPL")
            self%stage1_is_advanced = param%lencounter_sas_plpl
         case("ENCOUNTER_PLTP")
            self%stage1_is_advanced = param%lencounter_sas_pltp
         end select
         call io_log_one_message(logfile, trim(adjustl(self%loopname)) // ": loop timer turned on at t = " // trim(adjustl(tstr)))
      case(2)
         select case(trim(adjustl(self%looptype)))
         case("INTERACTION")
            param%lflatten_interactions = self%stage1_is_advanced 
         case("ENCOUNTER_PLPL")
            param%lencounter_sas_plpl= self%stage1_is_advanced 
         case("ENCOUNTER_PLTP")
            param%lencounter_sas_pltp= self%stage1_is_advanced 
         end select
         call self%flip(param, pl) 
      case default
         self%stage = 1
      end select

      write(schar,'(I1)') self%stage
      call io_log_one_message(logfile, trim(adjustl(self%loopname)) // ": stage " // schar )

      call self%reset()
      call self%start()

      return
   end subroutine walltime_interaction_time_this_loop

end submodule s_walltime