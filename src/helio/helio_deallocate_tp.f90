submodule (helio) s_helio_deallocate_tp
contains
   module procedure helio_deallocate_tp
   !! author: David A. Minton
   !!
   !! Finalizer for base helio particle class.
   !! Basic Helio test particle destructor/finalizer
   use swiftest
   implicit none

   if (self%is_allocated) then
      deallocate(self%ah)
      deallocate(self%ahi)
   end if

   return
   end procedure helio_deallocate_tp
end submodule s_helio_deallocate_tp
