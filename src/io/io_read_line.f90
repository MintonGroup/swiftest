submodule (nbody_data_structures) s_io_read_line
contains
   module procedure io_read_line
   !! author: David A. Minton
   !!
   !! Read a line (record) from input binary files
   !!     Function returns read error status (0 = OK, nonzero = ERROR)
   !! Adapted from David E. Kaufmann's Swifter modules: io_read_line.f90
   !! Adapted from Hal Levison's Swift routine io_read_line.f
use swiftest
implicit none
   logical        :: lmass, lradius
   integer( I4B)       :: ierr
   real(SP)         :: smass, sradius
   real(SP), dimension(6) :: svec
   real(DP), dimension(6) :: dvec

! executable code
   lmass = present(mass)
   if (lmass) then
      lradius = present(radius)
      if (.not. lradius) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   function io_read_line called with optional mass but without optional radius"
         call util_exit(failure)
      end if
   end if
   select case (out_type)
      case (real4_type, swifter_real4_type)
         if (lmass) then
            read(iu, iostat = ierr) name, smass, sradius, svec
         else
            read(iu, iostat = ierr) name, svec
         end if
         io_read_line = ierr
         if (ierr /= 0) return
         if (lmass) mass = smass
         d1 = svec(1); d2 = svec(2); d3 = svec(3); d4 = svec(4); d5 = svec(5); d6 = svec(6)
      case (real8_type, swifter_real8_type)
         if (lmass) then
            read(iu, iostat = ierr) name, mass, radius, dvec
         else
            read(iu, iostat = ierr) name, dvec
         end if
         io_read_line = ierr
         if (ierr /= 0) return
         d1 = dvec(1); d2 = dvec(2); d3 = dvec(3); d4 = dvec(4); d5 = dvec(5); d6 = dvec(6)
   end select

   return

   end procedure io_read_line
end submodule s_io_read_line
