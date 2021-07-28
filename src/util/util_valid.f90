submodule (swiftest_classes) s_util_valid
   use swiftest
contains
   module subroutine util_valid(pl, tp)
      !! author: David A. Minton
      !!
      !! Validate massive body and test particle ids
      !! Subroutine causes program to exit with error if any ids are not unique
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_valid.f90
      implicit none
      ! Arguments
      class(swiftest_pl), intent(in) :: pl
      class(swiftest_tp), intent(in) :: tp
      ! Internals
      integer(I4B)                  :: i
      integer(I4B), dimension(:), allocatable :: idarr

      associate(npl => pl%nbody, ntp => tp%nbody)
         allocate(idarr(npl+ntp))
         do i = 1, npl
            idarr(i) = pl%id(i)
         end do
         do i = 1, ntp
            idarr(npl+i) = tp%id(i)
         end do
         call util_sort(idarr)
         do i = 1, npl + ntp - 1
            if (idarr(i) == idarr(i+1)) then
               write(*, *) "Swiftest error:"
               write(*, *) "   more than one body/particle has id = ", idarr(i)
               call util_exit(FAILURE)
            end if
         end do
      end associate

      return

   end subroutine util_valid
end submodule s_util_valid
