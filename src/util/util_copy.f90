submodule (swiftest_classes) s_util_copy
   use swiftest
contains

   module subroutine util_copy_body(self, source, param)
      !! author: David A. Minton
      !!
      !! Non-destructively copy components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self   !! Swiftest body object
      class(swiftest_body),       intent(in)    :: source !! Source object to copy
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters

      return
   end subroutine util_copy_body

end submodule s_util_copy