submodule(rmvs_classes) s_rmvs_setup
contains
   module subroutine rmvs_setup_pl(self,n)
      !! author: David A. Minton
      !!
      !! Allocate RMVS test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine rmvs_setup.f90
      use swiftest
      implicit none
      ! Arguments
      class(rmvs_pl),                intent(inout) :: self !! RMVS test particle object
      integer,                       intent(in)    :: n    !! Number of test particles to allocate

      !> Call allocation method for parent class
      call whm_setup_pl(self, n) 
      if (n <= 0) return

      allocate(self%nenc(n))
      allocate(self%tpenc1P(n))
      allocate(self%xout(NDIM, n, 0:NTENC))
      allocate(self%vout(NDIM, n, 0:NTENC))
      allocate(self%xin(NDIM, n, 0:NTPHENC))
      allocate(self%vin(NDIM, n, 0:NTPHENC))
      allocate(self%aoblin(NDIM, n, 0:NTPHENC))

      self%nenc          = 0
      self%tpenc1P(:)    = 0
      self%xout(:,:,:)   = 0.0_DP
      self%vout(:,:,:)   = 0.0_DP
      self%xin(:,:,:)    = 0.0_DP
      self%vin(:,:,:)    = 0.0_DP
      self%aoblin(:,:,:) = 0.0_DP

      return
   end subroutine rmvs_setup_pl 

   module subroutine rmvs_setup_tp(self,n)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      use swiftest
      implicit none
      ! Arguments
      class(rmvs_tp),              intent(inout)   :: self !! RMVS test particle object
      integer,                     intent(in)      :: n    !! Number of test particles to allocate

      !> Call allocation method for parent class
      call whm_setup_tp(self, n) 
      if (n <= 0) return

      allocate(self%lperi(n))
      allocate(self%plperP(n))
      allocate(self%plencP(n))
      allocate(self%tpencP(n))

      self%lperi(:)  = .false.
      self%plperP(:) = 0
      self%plencP(:) = 0
      self%tpencP(:) = 0

      return
   end subroutine rmvs_setup_tp

   module subroutine rmvs_setup_system(self, config)
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      use swiftest
      implicit none
      ! Arguments
      class(rmvs_nbody_system),      intent(inout) :: self    !! RMVS system object
      class(swiftest_configuration), intent(inout) :: config  !! Input collection of  configuration parameters 

      ! Call parent method
      call whm_setup_system(self, config)

   end subroutine rmvs_setup_system
   
   module subroutine rmvs_setup_set_beg_end(self, xbeg, xend, vbeg)
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of xbeg, xend, and vbeg
      use swiftest
      implicit none
      ! Arguments
      class(rmvs_tp),                intent(inout) :: self !! RMVS test particle object
      real(DP), dimension(:,:),           optional :: xbeg, xend, vbeg

      if (present(xbeg)) then
         if (allocated(self%xbeg)) deallocate(self%xbeg)
         allocate(self%xbeg, source=xbeg)
      end if
      if (present(xend)) then
         if (allocated(self%xend)) deallocate(self%xend)
         allocate(self%xend, source=xend)
      end if
      if (present(vbeg)) then
         if (allocated(self%vbeg)) deallocate(self%vbeg)
         allocate(self%vbeg, source=vbeg)
      end if

      return

   end subroutine rmvs_setup_set_beg_end

end submodule s_rmvs_setup