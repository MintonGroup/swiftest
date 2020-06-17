submodule (swiftest_data_structures) s_swiftest_deallocate_tp
contains
   module procedure swiftest_deallocate_tp
   !! author: David A. Minton
   !!
   !! Finalizer for base Swiftest particle class.
   !! Basic Swiftest test particle destructor/finalizer
   implicit none

   type(swiftest_tp), intent(inout)    :: self

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
   end procedure swiftest_deallocate_tp
end submodule s_swiftest_deallocate_tp
