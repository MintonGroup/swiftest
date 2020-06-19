submodule (nbody_data_structures) s_nbody_deallocate_body
contains
   module procedure nbody_deallocate_body
   !! author: David A. Minton
   !!
   !! Finalizer for base Swiftest particle class. Deallocates all components and sets 
   !! is_allocated flag to false. Mostly this is redundant, so this serves as a placeholder
   !! in case future updates include pointers as part of the class.
   use swiftest
   implicit none
   
   if (self%is_allocated) then
      deallocate(self%name)
      deallocate(self%status)
      if (allocated(self%lspill_list)) deallocate(self%lspill_list)
      self%is_allocated = .false.
   end if
   return
   end procedure nbody_deallocate_body
end submodule s_nbody_deallocate_body
