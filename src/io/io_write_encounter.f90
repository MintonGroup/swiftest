submodule (nbody_data_structures) s_io_write_encounter
contains
   module procedure io_write_encounter
   !! author: David A. Minton
   !!
   !! Write close encounter data to output binary files
   !!  There is no direct file output from this subroutine
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: io_write_encounter.f90
   !! Adapted from Hal Levison's Swift routine io_write_encounter.f
use swiftest
implicit none
   logical(lgt)        :: lxdr
   logical(lgt), save    :: lfirst = .true.
   integer(I4B), parameter :: lun = 30
   integer(I4B)        :: ierr
   integer(I4B), save    :: iu = lun

! executable code

   !lxdr = ((out_type == swifter_real4_type) .or. (out_type == swifter_real8_type))
   !if (lxdr) then
   !   call io_open_fxdr(encounter_file, "a", .true., iu, ierr)
   !   if ((ierr /= 0) .and. lfirst) then
   !      call io_open_fxdr(encounter_file, "w", .true., iu, ierr)
   !      lfirst = .false.
   !   end if
   !   if (ierr /= 0) then
   !      write(*, *) "Swiftest Error:"
   !      write(*, *) "   unable to open binary encounter file"
   !      call util_exit(failure)
   !   end if
   !   ierr = ixdrdouble(iu, t)
   !   if (ierr < 0) then
   !      write(*, *) "Swiftest Error:"
   !      write(*, *) "   unable to write binary file record"
   !      call util_exit(failure)
   !   end if
   !   call io_write_line(iu, name1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), swifter_real8_type, mass = mass1, radius = radius1)
   !   call io_write_line(iu, name2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), vh2(3), swifter_real8_type, mass = mass2, radius = radius2)
   !   ierr = ixdrclose(iu)
   !else
      call io_open(iu, encounter_file, "append", "unformatted", ierr)
      if ((ierr /= 0) .and. lfirst) then
         call io_open(iu, encounter_file, "new", "unformatted", ierr)
         lfirst = .false.
      end if
      if (ierr /= 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   unable to open binary encounter file"
         call util_exit(failure)
      end if
      write(iu, iostat = ierr) t
      if (ierr < 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   unable to write binary file record"
         call util_exit(failure)
      end if
      call io_write_line(iu, name1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), real8_type, mass = mass1, radius = radius1)
      call io_write_line(iu, name2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), vh2(3), real8_type, mass = mass2, radius = radius2)
      close(unit = iu, iostat = ierr)
   !end if
   if (ierr /= 0) then
      write(*, *) "Swiftest Error:"
      write(*, *) "   unable to close binary encounter file"
      call util_exit(failure)
   end if

   return

   end procedure io_write_encounter
end submodule s_io_write_encounter
