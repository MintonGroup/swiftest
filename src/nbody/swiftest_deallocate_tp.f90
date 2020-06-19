submodule (nbody_data_structures) s_nbody_deallocate_tp
contains
   module procedure nbody_deallocate_tp
   !! author: David A. Minton
   !!
   !! Finalizer for base Swiftest particle class.
   !! Basic Swiftest test particle destructor/finalizer
   use swiftest
   implicit none

   if (self%is_allocated) then
      deallocate(self%isperi)
      deallocate(self%peri)
      deallocate(self%atp)
      deallocate(self%xh)
      deallocate(self%vh)
      deallocate(self%xb)
      deallocate(self%vb)
   end if
   return
   end procedure nbody_deallocate_tp
end submodule s_nbody_deallocate_tp
