submodule (symba) s_io_discard_write_symba
contains
   module procedure io_discard_write_symba
   !! author: David A. Minton
   !!
   !! Write out information about discarded and merged planets and test particles in SyMBA
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: io_discard_write_symba.f90
   !! Adapted from Hal Levison's Swift routine io_discard_mass.f and io_discard_merge.f
use swiftest
implicit none
   integer(I4B), parameter   :: lun = 40
   integer(I4B)          :: i, index, ncomp, ierr, nplm

! executable code
   call io_open(lun, fname, "append", "formatted", ierr)
   if (ierr /= 0) then
      call io_open(lun, fname, "new", "formatted", ierr)
      if (ierr /= 0) then
         write(*, *) "swiftest error:"
         write(*, *) "   unable to open discard output file, ", fname
         call util_exit(FAILURE)
      end if
   end if
   write(lun, 100) t, nsppl + nsptp + 2*nmergeadd, lbig_discard 
 100 format(e23.16, 1x, i8, 1x, l1)
   index = 0
   do i = 1, nmergeadd
      write(lun, 200) add, mergeadd_list%name(i), mergeadd_list%status(i)
 200    format(a, 2(1x, i8))
      write(lun, 300) mergeadd_list%xh(:,i)
 300    format(3(e23.16, 1x))
      write(lun, 300) mergeadd_list%vh(:,i)
      ncomp = mergeadd_list%ncomp(i)
      do index = 1, ncomp 
         write(lun, 200) sub, mergesub_list%name(index), mergesub_list%status(index)
         write(lun, 300) mergesub_list%xh(:,index)
         write(lun, 300) mergesub_list%vh(:,index)
         write(lun, 500) mergesub_list%name(index), mergesub_list%mass(index), mergesub_list%radius(index)
      end do
   end do
   do i = 1, nsppl
      if (discard_plA%status(i) /= MERGED) then
         write(lun, 200) sub, discard_plA%name(i), discard_plA%status(i)
         write(lun, 300) discard_plA%xh(1,i),discard_plA%xh(2,i), discard_plA%xh(3,i)
         write(lun, 300) discard_plA%vh(1,i),discard_plA%vh(2,i), discard_plA%vh(3,i)
         write(lun, 500) discard_plA%name(i),discard_plA%mass(i), discard_plA%radius(i)
      end if
   end do
   do i = 1, nsptp
      write(lun, 200) sub, discard_tpA%name(i), discard_tpA%status(i)
      write(lun, 300) discard_tpA%xh(1,i),discard_tpA%xh(2,i), discard_tpA%xh(3,i)
      write(lun, 300) discard_tpA%vh(1,i),discard_tpA%vh(2,i), discard_tpA%vh(3,i)
   end do
   if (lbig_discard) then
      nplm = 0
      do i = 1, npl
         if (symba_plA%mass(i) < mtiny) exit
         nplm = nplm + 1
      end do
      if (nplm > 1) then
         write(lun, 400) nplm
 400       format(i8)
         do i = 2, nplm
            write(lun, 500) symba_plA%name(i), symba_plA%mass(i),& 
             symba_plA%radius(i)
 500          format(i8, 2(1x, e23.16))
            write(lun, 300) symba_plA%xh(:,i)
            write(lun, 300) symba_plA%vh(:,i)
         end do
      end if
   end if
   close(lun)

   return

   end procedure io_discard_write_symba
end submodule s_io_discard_write_symba
