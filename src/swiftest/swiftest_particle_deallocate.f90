submodule (swiftest_data_structures) s_swiftest_particle_deallocate
contains
   module procedure swiftest_particle_deallocate
   !! author: David A. Minton
   !!
   !! Finalizer for base Swiftest particle class. Deallocates all components and sets 
   !! is_allocated flag to false. Mostly this is redundant, so this serves as a placeholder
   !! in case future updates include pointers as part of the class.
   if (self%is_allocated) then
      deallocate(self%name)
      deallocate(self%status)
      self%is_allocated = .false.
   end if
   return
   end procedure swiftest_particle_deallocate
end submodule s_swiftest_particle_deallocate
