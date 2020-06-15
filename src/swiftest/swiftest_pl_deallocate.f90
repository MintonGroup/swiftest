submodule (swiftest_data_structures) s_swiftest_pl_deallocate
contains
   module procedure swiftest_pl_deallocate
   !! author: David A. Minton
   !!
   !! Finalizer for base Swiftest massive body class. Deallocates all components and sets 
   !! is_allocated flag to false. Mostly this is redundant, so this serves as a placeholder
   !! in case future updates include pointers as part of the class.
   if (self%is_allocated) then
   deallocate(self%mass)
      deallocate(self%radius)
      deallocate(self%rhill)
   end if
   return
   end procedure swiftest_pl_deallocat
end submodule s_swiftest_pl_deallocate
