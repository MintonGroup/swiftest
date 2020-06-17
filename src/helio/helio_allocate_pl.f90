submodule (helio_data_structures) s_helio_allocate_pl
contains
   module procedure helio_allocate_pl
   !! author: David A. Minton
   !!
   !! Constructor for Helio massive body class. Allocates space for all particles and
   !! initializes all components with a value. 
   use swiftest
   implicit none

   class(helio_pl), intent(inout) :: self !! Swiftest massive body object
   integer, intent(in)            :: n    !! Number of massive bodies to allocate

   if (self%is_allocated) then
      write(*,*) 'Helio massive body particle structure already allocated'
      return
   end if
   write(*,*) 'Allocating the Helio massive body particle structure'
   
   call self%swiftest_tp%alloc(n)
   if (n <= 0) return
   allocate(self%ah(NDIM,n))
   allocate(self%ahi(NDIM,n))

   self%ah(:,:) = 0.0_DP
   self%ahi(:,:) = 0.0_DP
   return
   end procedure helio_allocate_pl
end submodule s_helio_allocate_pl
