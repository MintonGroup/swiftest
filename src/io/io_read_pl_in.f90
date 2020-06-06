submodule (io) s_io_read_pl_in
contains
   module procedure io_read_pl_in
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Read in massive body data 
   !!
   !! Adapted from David E. Kaufmann's Swifter routine io_init_pl.f90
   !! Adapted from Martin Duncan's Swift routine io_init_pl.f
   use swiftest
   use module_symba
   !use module_helio
   !use module_swiftest
   !use module_fxdr
   !use module_interfaces
   implicit none

   integer(I4B), parameter          :: LUN = 7              !! Unit number of input file
   integer(I4B)                     :: i, iu, inpl

   ierr = 0
   open(unit = LUN, file = param%inplfile, status = 'old', iostat = ierr)
   if (ierr /=  0) then
      write(*,*) 'Error opening massive body initial conditions file ',trim(adjustl(param%inplfile))
      return
   end if

   read(LUN, *, iostat = ierr) inpl
   read(LUN, *, iostat = ierr) swiftest_plA%name(1), swiftest_plA%mass(1)
   swiftest_plA%rhill(1) = 0.0_dp
   swiftest_plA%radius(1) = 0.0_dp
   read(LUN, *, iostat = ierr) swiftest_plA%xh(:,1)
   read(LUN, *, iostat = ierr) swiftest_plA%vh(:,1)
   if (ierr /= 0) then
      write(*,*) 'Error reading central body values in ',trim(adjustl(param%inplfile))
      return
   end if
   do i = 1, NDIM
      if ((swiftest_plA%xh(i,1) /= 0.0_dp) .or. (swiftest_plA%vh(i,1) /= 0.0_dp)) then
         write(*, *) "swiftest error:"
         write(*, *) " input must be in heliocentric coordinates."
         write(*, *) " position/velocity components of body 1 are"
         write(*, *) swiftest_plA%xh(:,1)
         write(*, *) swiftest_plA%vh(:,1)
      end if
   end do
   swiftest_plA%status(1) = active
   do i = 2, npl
      if (param%feature%lrhill_present) then
         read(LUN, *, iostat = ierr) swiftest_plA%name(i), swiftest_plA%mass(i), swiftest_plA%rhill(i)
      else
         read(LUN, *, iostat = ierr) swiftest_plA%name(i), swiftest_plA%mass(i)
         swiftest_plA%rhill(i) = 0.0_dp
      end if
      if (ierr /= 0 ) exit
      if (param%feature%lclose) then
         read(LUN, *, iostat = ierr) swiftest_plA%radius(i)
         if (ierr /= 0 ) exit
      else
          swiftest_plA%radius(i) = 0.0_dp
      end if
      read(LUN, *, iostat = ierr) swiftest_plA%xh(:,i)
      read(LUN, *, iostat = ierr) swiftest_plA%vh(:,i)
      if (ierr /= 0 ) exit
      swiftest_plA%status(i) = active
   end do
   close(unit = LUN)
   if (ierr /= 0 ) then
      write(*,*) 'Error reading in massive body initial conditions from ',trim(adjustl(param%inplfile))
   end if
   return
   end procedure io_read_pl_in

end submodule s_io_read_pl_in

