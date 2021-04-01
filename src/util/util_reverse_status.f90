submodule (swiftest_classes) s_util_reverse_status
   use swiftest
contains
   module procedure util_reverse_status
      !! author: David A. Minton
      !!
      !! Reverses the active/inactive status of all particles in a structure
      implicit none

      where (self%status(:) == ACTIVE)
         self%status(:) = INACTIVE
      elsewhere (self%status(:) == INACTIVE)
         self%status(:) = ACTIVE
      end where
   end procedure util_reverse_status
end submodule s_util_reverse_status