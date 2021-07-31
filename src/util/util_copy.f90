submodule (swiftest_classes) s_util_copy
   use swiftest
contains

   module subroutine util_copy_into_body(self, source, param, lmask)
      !! author: David A. Minton
      !!
      !! Copies elements from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),            intent(inout) :: self   !! Swiftest body object
      class(swiftest_body),            intent(in)    :: source !! Source object to append
      class(swiftest_parameters),      intent(in)    :: param  !! Current run configuration parameters
      logical, dimension(:), optional, intent(in)    :: lmask  !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B)  :: nnew

      if (present(lmask)) then
         nnew = count(lmask)
      else
         nnew = size(source%status)
      end if
      associate(nold => self%nbody)
         if (nnew > size(self%status)) call self%resize(nnew, param)

      end associate
      return
   end subroutine util_copy_into_body

end submodule s_util_copy