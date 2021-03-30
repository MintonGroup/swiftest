submodule(helio_classes) s_helio_setup
contains
   module procedure helio_setup_pl
      !! author: David A. Minton & Carlisle A. Wishard
      !!
      !! Allocate Helio planet structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine helio_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for great-grandparent class (we don't need Jacobi variables from WHM/RMVS)
      call setup_pl(self, n) 
      if (n <= 0) return

      allocate(self%ahi(NDIM, n))
      self%ahi(:,:) = 0.0_DP
      return
   end procedure helio_setup_pl 

   module procedure helio_setup_tp
      !! author: David A. Minton & Carlisle A. Wishard
      !!
      !! Allocate Helio test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine helio_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for great-grandparent class 
      call setup_tp(self, n) 
      if (n <= 0) return

      allocate(self%ahi(NDIM, n))
      self%ahi(:,:) = 0.0_DP

      return
   end procedure helio_setup_tp

end submodule s_helio_setup