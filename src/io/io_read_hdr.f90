submodule (swiftest_classes) s_io_read_hdr
contains
   module procedure io_read_hdr
   !! author: David A. Minton
   !!
   !! Read frame header from input binary files
   !!     Function returns read error status (0 = OK, nonzero = ERROR)
   !! Adapted from David E. Kaufmann's Swifter routine: io_read_hdr.f90
   !! Adapted from Hal Levison's Swift routine io_read_hdr.f
   use swiftest
   implicit none
   integer(I4B)         :: ierr
   integer(I4B), dimension(3) :: nn
   real(SP)             :: ttmp

! executable code
   select case (out_type)
      case (real4_type, swifter_real4_type)
         read(iu, iostat = ierr) ttmp, npl, ntp, iout_form
         io_read_hdr = ierr
         if (ierr /= 0) return
         t = ttmp
      case (real8_type, swifter_real8_type)
         read(iu, iostat = ierr) t, npl, ntp, iout_form
         io_read_hdr = ierr
   end select

   return

   end procedure io_read_hdr
end submodule s_io_read_hdr
