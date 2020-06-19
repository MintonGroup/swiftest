submodule (util) s_util_valid
contains
   module procedure util_valid
   !! author: David A. Minton
   !!
   !! Validate massive body and test particle ids
   !! Subroutine causes program to exit with error if any ids are not unique
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: util_valid.f90
   use swiftest
   integer(I4B)                  :: i
   integer(I4B), dimension(:), allocatable :: idarr

! executable code
   allocate(idarr(npl+ntp))
   do i = 1, npl
      idarr(i) = swiftest_plA%name(i)
   end do
   do i = 1, ntp
      idarr(npl+i) = swiftest_tpA%name(i)
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

   return

   end procedure util_valid
end submodule s_util_valid
