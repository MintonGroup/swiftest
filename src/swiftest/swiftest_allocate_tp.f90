submodule (swiftest_data_structures) s_swiftest_allocate_tp
contains
   module procedure swiftest_allocate_tp
   !! author: David A. Minton
   !!
   !! Constructor for base Swiftest test particle particle class. Allocates space for 
   !! all particles and initializes all components with a value. 
   use swiftest
   implicit none

   if (self%is_allocated) then
      !write(*,*) 'Swiftest test particle structure already alllocated'
      return
   end if
   !write(*,*) 'Allocating the Swiftest test particle'

   !> Call allocation method for parent class
   call self%swiftest_body%alloc(n)
   if (n <= 0) return

   allocate(self%peri(n))
   allocate(self%atp(n))
   allocate(self%isperi(n))
   allocate(self%xh(NDIM,n))
   allocate(self%vh(NDIM,n))
   allocate(self%xb(NDIM,n))
   allocate(self%vb(NDIM,n))

   self%peri(:) = 0.0_DP
   self%atp(:) = 0.0_DP
   self%isperi(:) = 0.0_DP
   self%xh(:,:) = 0.0_DP
   self%vh(:,:) = 0.0_DP
   self%xb(:,:) = 0.0_DP
   self%vb(:,:) = 0.0_DP

   return
   end procedure swiftest_allocate_tp
end submodule s_swiftest_allocate_tp

