submodule (swiftest_classes) s_util_valid
   use swiftest
contains

   module subroutine util_valid_id_system(self, param)
      !! author: David A. Minton
      !!
      !! Validate massive body and test particle ids
      !! Subroutine causes program to exit with error if any ids are not unique
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_valid.f90
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                  :: i
      integer(I4B), dimension(:), allocatable :: idarr

      associate(cb => self%cb, pl => self%pl, npl => self%pl%nbody, tp => self%tp, ntp => self%tp%nbody)
         allocate(idarr(1+npl+ntp))
         idarr(1) = cb%id
         do i = 1, npl
            idarr(1+i) = pl%id(i)
         end do
         do i = 1, ntp
            idarr(1+npl+i) = tp%id(i)
         end do
         call util_sort(idarr)
         do i = 1, npl + ntp 
            if (idarr(i) == idarr(i+1)) then
               write(*, *) "Swiftest error:"
               write(*, *) "   more than one body/particle has id = ", idarr(i)
               call util_exit(FAILURE)
            end if
         end do
         param%maxid = max(param%maxid, maxval(idarr))
      end associate

      return
   end subroutine util_valid_id_system

end submodule s_util_valid
