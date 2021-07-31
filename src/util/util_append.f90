submodule (swiftest_classes) s_util_append
   use swiftest
contains

   module subroutine util_append_body(self, source, param, lmask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),            intent(inout) :: self   !! Swiftest body object
      class(swiftest_body),            intent(in)    :: source !! Source object to append
      class(swiftest_parameters),      intent(in)    :: param  !! Current run configuration parameters
      logical, dimension(:), optional, intent(in)    :: lmask  !! Logical mask indicating which elements to append to

      associate(nold => self%nbody, nnew => source%nbody)
         if (nnew > size(self%status)) call self%resize(nnew, param)

      end associate
      return
   end subroutine util_append_body

end submodule s_util_append