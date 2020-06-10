submodule (swiftest_data_structures) s_swiftest_read_tp_in
contains
   module procedure swiftest_read_tp_in
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Read in test particle data
   !!
   !! Adapted from David E. Kaufmann's Swifter routine io_init_pl.f90
   !! Adapted from Martin Duncan's Swift routine io_init_tp.f
   implicit none

   integer(I4B), parameter          :: LUN = 7              !! Unit number of input file
   integer(I4B)                     :: i, iu, ierr, ntp

! Executable code
   ierr = 0
   open(unit = LUN, file = param%intpfile, status = 'old', iostat = ierr)
   if (ierr /=  0) then
      write(*,*) 'Error opening test particle initial conditions file ',trim(adjustl(param%intpfile))
      return
   end if
   read(lun, *) ntp
   if (ntp <= 0) return

   call self%alloc(ntp)

   do i = 1, self%nbody
      read(LUN, *) self%name(i)
      read(LUN, *) self%xh(:,i)
      read(LUN, *) self%vh(:,i)
      self%status(i) = ACTIVE
   end do
   close(LUN)

   return
   end procedure swiftest_read_tp_in

end submodule s_swiftest_read_tp_in

