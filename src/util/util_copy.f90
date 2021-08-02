submodule (swiftest_classes) s_util_copy
   use swiftest
contains

   module subroutine util_copy_into_body(self, source, param, lsource_mask)
      !! author: David A. Minton
      !!
      !! Copies elements from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),            intent(inout) :: self   !! Swiftest body object
      class(swiftest_body),            intent(in)    :: source !! Source object to append
      class(swiftest_parameters),      intent(in)    :: param  !! Current run configuration parameters
      logical, dimension(:), optional, intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B)  :: i,nnew
      logical, dimension(:), allocatable :: lfill_list

      if (present(lsource_mask)) then
         nnew = count(lsource_mask)
      else
         nnew = size(source%status)
      end if
      allocate(lfill_list(size(self%status)))
      lfill_list = .false.
      lfill_list(1:nnew) = .true.
      associate(nold => self%nbody)
         if (nnew > size(self%status)) call self%resize(nnew, param)
         call self%fill(source, lfill_list)
      end associate
      return
   end subroutine util_copy_into_body

end submodule s_util_copy