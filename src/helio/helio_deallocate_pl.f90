submodule (helio) s_helio_deallocate_pl
contains
   module procedure helio_deallocate_pl
   !! author: David A. Minton
   !!
   !! Finalizer for base helio massive body class. Deallocates all components and sets 
   !! is_allocated flag to false. Mostly this is redundant, so this serves as a placeholder
   !! in case future updates include pointers as part of the class.
   use swiftest
   implicit none

   return
   end procedure helio_deallocate_pl
end submodule s_helio_deallocate_pl
