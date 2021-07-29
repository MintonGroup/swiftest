submodule (swiftest_classes) s_util_reverse_status
   use swiftest
contains

   module subroutine util_reverse_status(self)
      !! author: David A. Minton
      !!
      !! Reverses the active/inactive status of all particles in a structure
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self !! Swiftest body object

      where (self%status(:) == ACTIVE)
         self%status(:) = INACTIVE
      elsewhere (self%status(:) == INACTIVE)
         self%status(:) = ACTIVE
      end where
      self%lmask(:) = self%status(:) == ACTIVE
      
      return
   end subroutine util_reverse_status

end submodule s_util_reverse_status