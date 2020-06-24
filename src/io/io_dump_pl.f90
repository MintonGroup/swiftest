submodule (swiftest_classes) s_io_dump_pl
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

   open(unit = LUN, file = DUMP_PL_FILE(idx), form = "unformatted", status = 'replace', iostat = ierr)
   if (ierr /= 0) then
      write(*, *) "Swiftest error:"
      write(*, *) "   Unable to open binary dump file ", trim(DUNMP_PL_FILE(idx))
      call util_exit(FAILURE)
   end if
   write(LUN) npl
   if (npl > 0) then
      write(LUN) swiftest_plA%name(1:npl)
         write(LUN) swiftest_plA%mass(npl:)
      if (lrhill_present) write(LUN) swiftest_plA%rhill(1:npl) 
      if (lclose) write(LUN) swiftest_plA%radius(1:npl) 
      write(LUN) swiftest_plA%xh(:,1:npl)
      write(LUN) swiftest_plA%vh(:,1:npl)
   end if
   close(LUN)
   idx = idx + 1
   if (idx > 2) idx = 1

   return

 end procedure io_dump_pl
end submodule s_io_dump_pl
