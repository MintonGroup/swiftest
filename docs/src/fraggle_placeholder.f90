submodule(fraggle_classes) s_fraggle_placeholder
   use swiftest

contains

   !> The following interfaces are placeholders intended to satisfy the required abstract methods given by the parent class
   module subroutine fraggle_placeholder_accel(self, system, param, t, lbeg)
      implicit none
      class(fraggle_fragments),     intent(inout) :: self   !! Fraggle fragment system object 
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      write(*,*) "The type-bound procedure 'accel' is not defined for type fraggle_fragments"
      return
   end subroutine fraggle_placeholder_accel

   module subroutine fraggle_placeholder_kick(self, system, param, t, dt, lbeg)
      implicit none
      class(fraggle_fragments),     intent(inout) :: self   !! Fraggle fragment system object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system objec
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 

      write(*,*) "The type-bound procedure 'kick' is not defined for type fraggle_fragments"
      return
   end subroutine fraggle_placeholder_kick

   module subroutine fraggle_placeholder_step(self, system, param, t, dt)
      implicit none
      class(fraggle_fragments),     intent(inout) :: self   !! Swiftest body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Simulation time
      real(DP),                     intent(in)    :: dt     !! Current stepsize

      write(*,*) "The type-bound procedure 'step' is not defined for type fraggle_fragments"
      return
   end subroutine fraggle_placeholder_step


end submodule s_fraggle_placeholder