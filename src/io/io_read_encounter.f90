submodule (swiftest_data_structures) s_io_read_encounter
contains
   module procedure io_read_encounter
   !! author: David A. Minton
   !!
   !! Read close encounter data from input binary files
   !!     Other than time t, there is no direct file input from this function
   !!     Function returns read error status (0 = OK, nonzero = ERROR)
   !! Adapted from David E. Kaufmann's Swifter modules: io_read_encounter.f90
use swiftest
implicit none
   logical(lgt)        :: lxdr
   logical(lgt), save    :: lfirst = .true.
   integer(I4B), parameter :: lun = 30
   integer(I4B)        :: ierr
   integer(I4B), save    :: iu = lun

! executable code
   if (lfirst) then
       call io_open(iu, encounter_file, "old", "unformatted", ierr)
      if (ierr /= 0) then
         write(*, *) "swifter error:"
         write(*, *) "   unable to open binary encounter file"
         call util_exit(failure)
      end if
      lfirst = .false.
   end if
    read(iu, iostat = ierr) t
    io_read_encounter = ierr
    if (ierr /= 0) then
       close(unit = iu, iostat = ierr)
       return
    end if
    ierr = io_read_line(iu, name1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), real8_type, mass = mass1)
    io_read_encounter = ierr
    if (ierr /= 0) then
       close(unit = iu, iostat = ierr)
       return
    end if
    ierr = io_read_line(iu, name2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), vh2(3), real8_type, mass = mass2)
    io_read_encounter = ierr
    if (ierr /= 0) then
       close(unit = iu, iostat = ierr)
       return
    end if

   return

   end procedure io_read_encounter
end submodule s_io_read_encounter
