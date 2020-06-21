submodule (swiftest_classes) s_io_write_hdr
contains   
   module procedure io_write_hdr
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Write frame header to output binary file
   !!
   !! Adapted from David Adapted from David E. Kaufmann's Swifter routine io_write_hdr.f90
   !! Adapted from Hal Levison's Swift routine io_write_hdr.F
   implicit none

   integer(I4B)               :: ierr !! Error code

   select case (out_type)
      case (REAL4_TYPE,SWIFTER_REAL4_TYPE)
         write(iu, iostat = ierr) real(t, kind=SP), npl, ntp, iout_form
         if (ierr < 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Unable to write binary file header"
            call util_exit(FAILURE)
         end if
      case (REAL8_TYPE,SWIFTER_REAL8_TYPE)
         write(iu, iostat = ierr) t, npl, ntp, iout_form
         if (ierr < 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Unable to write binary file header"
            call util_exit(FAILURE)
         end if
   end select

   return

   end procedure io_write_hdr
end submodule s_io_write_hdr
