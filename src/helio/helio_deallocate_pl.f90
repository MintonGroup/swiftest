submodule (helio_data_structures) s_helio_deallocate_pl
contains
   module procedure helio_deallocate_pl
   !! author: David A. Minton
   !!
   !! Finalizer for base Helio massive body class. Deallocates all components and sets 
   !! is_allocated flag to false. Mostly this is redundant, so this serves as a placeholder
   !! in case future updates include pointers as part of the class.
   if (self%is_allocated) then
   deallocate(self%mass)
      deallocate(self%radius)
      deallocate(self%rhill)
   end if
   return
   end procedure helio_deallocate_pl
end submodule s_helio_deallocate_pl
