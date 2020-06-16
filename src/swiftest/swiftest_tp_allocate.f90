submodule (swiftest_data_structures) s_swiftest_tp_allocate
contains
   module procedure swiftest_tp_allocate
   !! author: David A. Minton
   !!
   !! Constructor for base Swiftest test particle particle class. Allocates space for 
   !! all particles and initializes all components with a value. 
   implicit none

   class(swiftest_tp), intent(inout)    :: self !! Swiftest test particle object
   integer, intent(in)                  :: n    !! Number of test particles to allocate
  
   call self%swiftest_particle%alloc(n)
   if (n <= 0) return

   if (self%is_allocated) then
      !write(*,*) 'Swiftest test particle structure already alllocated'
      return
   end if
   write(*,*) 'Allocating the Swiftest test particle'

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
   end procedure swiftest_tp_allocate
end submodule s_swiftest_tp_allocate

