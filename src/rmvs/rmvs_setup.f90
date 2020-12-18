submodule(rmvs_classes) s_rmvs_setup
contains
   module procedure rmvs_setup_pl
      !! author: David A. Minton
      !!
      !! Allocate RMVS test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine rmvs_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for parent class
      call whm_setup_pl(self, n) 
      if (n <= 0) return

      allocate(self%nenc(n))
      allocate(self%xpc(NDIM, n))
      allocate(self%tpenc1P(n))

      self%nenc       = 0
      self%xpc(:,:)   = 0.0_DP
      self%tpenc1P(:) = 0

      return
   end procedure rmvs_setup_pl 

   module procedure rmvs_setup_tp
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for parent class
      call whm_setup_tp(self, n) 
      if (n <= 0) return

      allocate(self%lperi(n))
      allocate(self%xpc(NDIM, n))
      allocate(self%vpc(NDIM, n))
      allocate(self%apc(NDIM, n))
      allocate(self%plperP(n))
      allocate(self%plencP(n))
      allocate(self%tpencP(n))

      self%lperi(:)  = .false.
      self%xpc(:,:)  = 0.0_DP
      self%vpc(:,:)  = 0.0_DP
      self%apc(:,:)  = 0.0_DP
      self%plperP(:) = 0
      self%plencP(:) = 0
      self%tpencP(:) = 0

      return
   end procedure rmvs_setup_tp

   module procedure rmvs_setup_system
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      implicit none

      ! Call parent method
      call whm_setup_system(self, config)

   end procedure rmvs_setup_system


end submodule s_rmvs_setup
