submodule (helio) s_helio_allocate_pl
contains
   module procedure helio_allocate_pl
   !! author: David A. Minton
   !!
   !! Constructor for Helio massive body class. Allocates space for all particles and
   !! initializes all components with a value. 
   use swiftest
   implicit none

   if (self%is_allocated) then
      !write(*,*) 'Helio massive body particle structure already allocated'
      return
   end if
   !write(*,*) 'Allocating the Helio massive body particle structure'
  
   !> Call allocation method for parent clas
   call self%helio_tp%alloc(n)
   if (n <= 0) return
   ! No new structures to allocate and initialize

   return
   end procedure helio_allocate_pl
end submodule s_helio_allocate_pl
