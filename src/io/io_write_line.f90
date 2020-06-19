submodule (nbody_data_structures) s_io_write_line
contains
   module procedure io_write_line
   !! author: David A. Minton
   !!
   !! Write a line (record) to output binary files
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: io_write_line.f90
   !! Adapted from Hal Levison's Swift routine io_write_line.f
use swiftest
implicit none
   integer( I4B)       :: ierr
   real(DP), dimension(6) :: dvec
   real(SP), dimension(6) :: svec
   real(SP)         :: smass, sradius
   logical        :: lmass, lradius

! executable code
   dvec(1) = d1; dvec(2) = d2; dvec(3) = d3; dvec(4) = d4; dvec(5) = d5; dvec(6) = d6
   svec(1) = d1; svec(2) = d2; svec(3) = d3; svec(4) = d4; svec(5) = d5; svec(6) = d6
   lmass = present(mass)
   if (lmass) then
      lradius = present(radius)
      if (.not. lradius) then
         write(*, *) "swiftest error:"
         write(*, *) "   io_write_line called with optional mass but without optional radius"
         call util_exit(failure)
      end if
      smass = mass
      sradius = radius
   end if
   select case (out_type)
      case (real4_type, swifter_real4_type)
         if (lmass) then
            write(iu, iostat = ierr) name, smass, sradius, svec
         else
            write(iu, iostat = ierr) name, svec
         end if
         if (ierr < 0) then
            write(*, *) "swiftest error:"
            write(*, *) "   unable to write binary file record"
            call util_exit(failure)
         end if
      case (real8_type, swifter_real8_type)
         if (lmass) then
            write(iu, iostat = ierr) name, mass, radius, dvec
         else
            write(iu, iostat = ierr) name, dvec
         end if
         if (ierr < 0) then
            write(*, *) "swiftest error:"
            write(*, *) "   unable to write binary file record"
            call util_exit(failure)
         end if
   end select

   return

   end procedure io_write_line
end submodule s_io_write_line
