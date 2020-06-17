submodule (helio_data_structures) s_helio_tp_deallocate
contains
   module procedure helio_tp_deallocate
   !! author: David A. Minton
   !!
   !! Finalizer for base Helio particle class.
   !! Basic Helio test particle destructor/finalizer
   implicit none

   if (self%is_allocated) then
      deallocate(self%ah)
      deallocate(self%ahi)
   end if
   return

   end procedure helio_tp_deallocate
end submodule s_helio_tp_deallocate
