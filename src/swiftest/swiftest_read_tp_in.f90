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

   integer(I4B), parameter  :: LUN = 7              !! Unit number of input file
   integer(I4B)             :: i, iu, ierr, ntp
   logical                  :: is_ascii


! Executable code
   ierr = 0
   is_ascii = (param%in_type == 'ASCII')  
   if (is_ascii) then
      open(unit = LUN, file = param%intpfile, status = 'old', form = 'formatted', iostat = ierr)
   else
      open(unit = LUN, file = param%intpfile, status = 'old', form = 'unformatted', iostat = ierr)
   end if
   if (ierr /=  0) then
      write(*,*) 'Error opening test particle initial conditions file ',trim(adjustl(param%intpfile))
      return
   end if
   if (is_ascii) then
      read(lun, *) ntp
   else
      read(lun) ntp
   end if
   if (ntp <= 0) return

   call self%alloc(ntp)

   if (is_ascii) then
      do i = 1, self%nbody
         read(LUN, *) self%name(i)
         read(LUN, *) self%xh(:,i)
         read(LUN, *) self%vh(:,i)
         self%status(i) = ACTIVE
      end do
   else
      read(LUN) self%name(:)
      read(LUN) self%xh(:,:)
      read(LUN) self%vh(:,:)
      self%status(:) = ACTIVE
   end if
   close(LUN)

   return
   end procedure swiftest_read_tp_in

end submodule s_swiftest_read_tp_in

