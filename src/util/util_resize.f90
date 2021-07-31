submodule (swiftest_classes) s_util_resize
   use swiftest
contains

   module subroutine util_resize_body(self, nrequested, param)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self       !! Swiftest body object
      integer(I4B),               intent(in)    :: nrequested !! New size neded
      class(swiftest_parameters), intent(in)    :: param      !! Current run configuration parameters
      ! Internals
      class(swiftest_body), allocatable   :: temp
      integer(I4B)                        :: nold
      logical                             :: lmalloc

      lmalloc = allocated(self%status)
      if (lmalloc) then
         nold = size(self%status)
      else
         nold = 0
      end if
      if (nrequested > nold) then
         if (lmalloc) allocate(temp, source=self)
         call self%setup(nrequested, param)
         if (lmalloc) then
            call self%copy_into(temp, param)
            deallocate(temp)
         end if
      else
         self%status(nrequested+1:nold) = INACTIVE
      end if
      self%nbody = nrequested

      return
   end subroutine util_resize_body

end submodule s_util_resize