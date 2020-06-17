submodule (swiftest_data_structures) s_swiftest_allocate_pl
contains
   module procedure swiftest_allocate_pl
   !! author: David A. Minton
   !!
   !! Constructor for base Swiftest massive body class. Allocates space for all particles and
   !! initializes all components with a value. 
   implicit none

   class(swiftest_pl), intent(inout)    :: self !! Swiftest massive body object
   integer, intent(in)                  :: n    !! Number of massive bodies to allocate

   call self%swiftest_body%alloc(n)
   if (n <= 0) return

   if (self%is_allocated) then
      !write(*,*) 'Swiftest massive body structure already alllocated'
      return
   end if

   call self%swiftest_tp%alloc(n)

   allocate(self%mass(n))
   allocate(self%radius(n))
   allocate(self%rhill(n))

   self%mass(:) = 0.0_DP
   self%radius(:) = 0.0_DP
   self%rhill(:) = 0.0_DP
   return
   end procedure swiftest_allocate_pl
end submodule s_swiftest_allocate_pl
