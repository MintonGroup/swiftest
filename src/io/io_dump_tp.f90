submodule (swiftest_data_structures) s_io_dump_tp
contains
   module procedure io_dump_tp
   !! author: David A. Minton
   !!
   !! Dump test particle data to files
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: io_dump_tp.f90
   !! Adapted from Hal Levison's Swift routine io_dump_tp.f
   use swiftest
   implicit none
   integer(I4B)                :: i, iu, ierr
   integer(I4B), save            :: idx = 1
   integer(I4B), parameter         :: LUN = 7

   open(unit = LUN, file = dump_tp_file(idx), form = "unformatted", status = 'replace', iostat = ierr)

   if (ierr /= 0) then
      write(*, *) "Swiftest error:"
      write(*, *) "   Unable to open binary dump file ", trim(dump_tp_file(idx))
      call util_exit(failure)
   end if
   write(LUN) ntp
   write(LUN) swiftest_tpa%name(:)
   write(LUN) swiftest_tpa%xh(:,:)
   write(LUN) swiftest_tpa%vh(:,:)
   close(LUN)
   idx = idx + 1
   if (idx > 2) idx = 1

   return

 end procedure io_dump_tp
end submodule s_io_dump_tp
