submodule (encounter_classes) s_encounter_setup
   use swiftest
contains

   module subroutine encounter_setup_list(self, n)
      !! author: David A. Minton
      !!
      !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      !!
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self !! Swiftest encounter structure
      integer(I4B),              intent(in)    :: n    !! Number of encounters to allocate space for

      if (n < 0) return

      if (allocated(self%lvdotr)) deallocate(self%lvdotr)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%index1)) deallocate(self%index1)
      if (allocated(self%index2)) deallocate(self%index2)
      if (allocated(self%id1)) deallocate(self%id1)
      if (allocated(self%id2)) deallocate(self%id2)
      if (allocated(self%x1)) deallocate(self%x1)
      if (allocated(self%x2)) deallocate(self%x2)
      if (allocated(self%v1)) deallocate(self%v1)
      if (allocated(self%v2)) deallocate(self%v2)
      if (allocated(self%t)) deallocate(self%t)

      self%nenc = n
      if (n == 0) return

      allocate(self%lvdotr(n))
      allocate(self%status(n))
      allocate(self%index1(n))
      allocate(self%index2(n))
      allocate(self%id1(n))
      allocate(self%id2(n))
      allocate(self%x1(NDIM,n))
      allocate(self%x2(NDIM,n))
      allocate(self%v1(NDIM,n))
      allocate(self%v2(NDIM,n))
      allocate(self%t(n))

      self%lvdotr(:) = .false.
      self%status(:) = INACTIVE
      self%index1(:) = 0
      self%index2(:) = 0
      self%id1(:) = 0
      self%id2(:) = 0
      self%x1(:,:) = 0.0_DP
      self%x2(:,:) = 0.0_DP
      self%v1(:,:) = 0.0_DP
      self%v2(:,:) = 0.0_DP
      self%t(:) = 0.0_DP

      return
   end subroutine encounter_setup_list

end submodule s_encounter_setup

 