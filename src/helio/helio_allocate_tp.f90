submodule (helio) s_helio_allocate_tp
contains
   module procedure helio_allocate_tp
   !! author: David A. Minton
   !!
   !! Constructor for Helio test particle particle class. Allocates space for 
   !! all particles and initializes all components with a value.
   use swiftest
   implicit none

   if (self%is_allocated) then
      !write(*,*) 'Helio test particle structure already allocated'
      return
   end if
   !write(*,*) 'Allocating the Helio test particle structure'
  
   !> Call allocation method for parent class
   call self%swiftest_pl%alloc(n)
   if (n <= 0) return
   allocate(self%ah(NDIM,n))
   allocate(self%ahi(NDIM,n))

   self%ah(:,:) = 0.0_DP
   self%ahi(:,:) = 0.0_DP
   
   return
   end procedure helio_allocate_tp
end submodule s_helio_allocate_tp

