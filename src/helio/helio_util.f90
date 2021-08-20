submodule(helio_classes) s_helio_eucl
   use swiftest
contains

   module subroutine helio_util_index_eucl_plpl(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper for the indexing method for WHM massive bodies. Sorts the massive bodies by heliocentric distance and then flattens the pl-pl upper triangular matrix
      implicit none
      ! Arguments
      class(helio_pl),            intent(inout) :: self  !! Helio massive body object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters

      call self%sort("mass", ascending=.false.)
      call util_index_eucl_plpl(self, param)

      return
   end subroutine helio_util_index_eucl_plpl

end submodule s_helio_eucl
   
