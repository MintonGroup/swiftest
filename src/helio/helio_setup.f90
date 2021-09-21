submodule(helio_classes) s_helio_setup
   use swiftest
contains

   module subroutine helio_setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a Helio nbody system from files, converting all heliocentric quantities to barycentric.
      !!
      implicit none
      ! Arguments
      class(helio_nbody_system),  intent(inout) :: self   !! Helio nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 

      call whm_setup_initialize_system(self, param)
      call self%pl%h2b(self%cb)
      call self%tp%h2b(self%cb)
      call self%pl%sort("mass", ascending=.false.)
      call self%pl%flatten(param)

      return
   end subroutine helio_setup_initialize_system

end submodule s_helio_setup