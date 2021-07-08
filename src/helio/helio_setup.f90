submodule(helio_classes) s_helio_setup
   use swiftest
contains
   module subroutine helio_setup_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a Helio nbody system from files 
      implicit none
      ! Arguments
      class(helio_nbody_system),  intent(inout) :: self   !! Helio system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 

      call io_read_initialize_system(self, param)
      ! Make sure that the discard list gets allocated initially
      call self%tp_discards%setup(self%tp%nbody)
   end subroutine helio_setup_system

   module procedure helio_setup_pl
      !! author: David A. Minton & Carlisle A. Wishard
      !!
      !! Allocate Helio planet structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine helio_setup.f90
      implicit none

      !> Call allocation method for great-grandparent class (we don't need Jacobi variables from WHM/RMVS)
      call setup_pl(self, n) 
      if (n <= 0) return

      allocate(self%ah(NDIM, n))
      self%ah(:,:) = 0.0_DP
      return
   end procedure helio_setup_pl 

   module procedure helio_setup_tp
      !! author: David A. Minton & Carlisle A. Wishard
      !!
      !! Allocate Helio test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine helio_setup.f90
      implicit none

      !> Call allocation method for great-grandparent class 
      call setup_tp(self, n) 
      if (n <= 0) return

      allocate(self%ah(NDIM, n))
      self%ah(:,:) = 0.0_DP

      return
   end procedure helio_setup_tp

end submodule s_helio_setup