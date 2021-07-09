submodule (util) s_util_valid
   use swiftest
contains
   module procedure util_valid
   !! author: David A. Minton
   !!
   !! Validate massive body and test particle ids
   !! Subroutine causes program to exit with error if any ids are not unique
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: util_valid.f90
   integer(I4B)                  :: i
   integer(I4B), dimension(:), allocatable :: idarr

! executable code
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
      deallocate(idarr)
   end associate

   return

   end procedure util_valid
end submodule s_util_valid
