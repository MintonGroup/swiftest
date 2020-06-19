submodule (nbody_data_structures) s_io_dump_pl
contains
   module procedure io_dump_pl
   !! author: David A. Minton
   !!
   !! Dump planet data to files
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: io_dump_pl.f90
   !! Adapted from Hal Levison's Swift routine io_dump_pl.f
   use swiftest
   implicit none
   integer(I4B)             :: i, iu, ierr
   integer(I4B), save         :: idx = 1
   integer(I4B),parameter         :: LUN = 7

   open(unit = LUN, file = DUNMP_PL_FILE(idx), form = "unformatted", status = 'replace', iostat = ierr)
   if (ierr /= 0) then
      write(*, *) "Swiftest error:"
      write(*, *) "   Unable to open binary dump file ", trim(DUNMP_PL_FILE(idx))
      call util_exit(FAILURE)
   end if
   write(LUN) npl
   write(LUN) swiftest_plA%name(:)
   write(LUN) swiftest_plA%mass(:)
   if (lrhill_present) write(LUN) swiftest_plA%rhill(:) 
   if (lclose) write(LUN) swiftest_plA%radius(:) 
   write(LUN) swiftest_plA%xh(:,:)
   write(LUN) swiftest_plA%vh(:,:)
   close(LUN)
   idx = idx + 1
   if (idx > 2) idx = 1

   return

 end procedure io_dump_pl
end submodule s_io_dump_pl
