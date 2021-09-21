module walltime_classes
   !! author: David A. Minton
   !!
   !! Classes and methods used to compute elasped wall time
   use swiftest_globals
   implicit none
   public

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
   end type

   interface
      module subroutine walltime_finish(self, nsubsteps, message)
         implicit none
         class(walltimer), intent(inout) :: self      !! Walltimer object
         integer(I4B),     intent(in)    :: nsubsteps !! Number of substeps used to compute the time per step 
         character(len=*), intent(in)    :: message   !! Message to prepend to the wall time terminal output
      end subroutine walltime_finish

      module subroutine walltime_reset(self)
         implicit none
         class(walltimer), intent(inout) :: self  !! Walltimer object
      end subroutine walltime_reset 

      module subroutine walltime_start(self)
         implicit none
         class(walltimer), intent(inout) :: self  !! Walltimer object
      end subroutine walltime_start
   end interface

   contains

      module subroutine walltime_finish(self, nsubsteps, message)
         !! author: David A. Minton
         !!
         !! Ends the timer, setting step_finish to the current ticker value and printing the elapsed time information to the terminal
         use swiftest_globals
         implicit none
         ! Arguments
         class(walltimer), intent(inout) :: self      !! Walltimer object
         integer(I4B),     intent(in)    :: nsubsteps !! Number of substeps used to compute the time per step 
         character(len=*), intent(in)    :: message   !! Message to prepend to the wall time terminal output
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

         call self%start()

         return
      end subroutine walltime_finish

      module subroutine walltime_reset(self)
         !! author: David A. Minton
         !!
         !! Resets the clock ticker, settting main_start to the current ticker value
         use swiftest_globals
         implicit none
         ! Arguments
         class(walltimer), intent(inout) :: self
         
         call system_clock(self%count_start_main, self%count_rate, self%count_max)
         self%lmain_is_started = .true.
         call self%start()

         return
      end subroutine walltime_reset 


      module subroutine walltime_start(self)
         !! author: David A. Minton
         !!
         !! Starts the timer, setting step_start to the current ticker value
         use swiftest_globals
         implicit none
         ! Arguments
         class(walltimer), intent(inout) :: self

         if (.not.self%lmain_is_started) then
            write(*,*) "Wall timer error: Cannot start the step time until reset is called at least once!"
            return
         end if

         call system_clock(self%count_start_step)

         return 
      end subroutine walltime_start


end module walltime_classes