submodule(swiftest_classes) s_tides_step_spin
   use swiftest
contains
   module subroutine tides_step_spin_system(self, param, t, dt)
      !! author: Jennifer L.L. Pouplin and David A. Minton
      !!
      !! Integrates the spin equations for central and massive bodies of the system subjected to tides.
      implicit none
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters  
      real(DP),                     intent(in)    :: t     !! Simulation time
      real(DP),                     intent(in)    :: dt    !! Current stepsize
      return
   end subroutine tides_step_spin_system
end submodule s_tides_step_spin