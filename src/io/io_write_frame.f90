submodule (io) s_io_write_frame
contains
   module procedure io_write_frame
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Write a frame (header plus records for each plAnet and active test particle) to output binary file
   !! There is no direct file output from this subroutine
   !!
   !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
   !! Adapted from Hal Levison's Swift routine io_write_frame.F
   use module_interfaces
   implicit none
   

   logical                   :: lxdr
   logical, save             :: lfirst = .true.
   integer(I4B), parameter   :: lun = 20
   integer(I4B)              :: i, j, ierr
   integer(I4B), save        :: iu = lun, iout_form = xv
   real(DP)                  :: a, e, inc, capom, omega, capm, mu
   real(DP), dimension(NDIM) :: xtmp, vtmp

   if (lfirst) then
      select case(out_stat)
      case('APPEND')
         open(unit = iu, file = outfile, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr)
      case('NEW')
         open(unit = iu, file = outfile, status = 'NEW', form = 'UNFORMATTED', iostat = ierr)
      case ('REPLACE')
         open(unit = iu, file = outfile, status = 'REPLACE', form = 'UNFORMATTED', iostat = ierr)
      case default
         write(*,*) 'Invalid status code',trim(adjustl(out_stat))
         call util_exit(FAILURE)
      end select
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   binary output file already exists or cannot be accessed"
         call util_exit(FAILURE)
      end if

      select case (out_form)
      case ("EL")
         iout_form = EL
      case ("XV")
         iout_form = XV
      end select
      lfirst = .false.
   else
      open(unit = iu, file = outfile, status = 'OLD', position =  'APPEND', form = 'UNFORMATTED', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   unable to open binary output file for APPEND"
         call util_exit(FAILURE)
      end if
   end if

   call io_write_hdr(iu, t, swiftest_plA%npl, swiftest_tpA%ntp, iout_form, out_type)
   select case (iout_form)
   case (EL)
      do i = 2, swiftest_plA%npl
         mu = swiftest_plA%mass(1) + swiftest_plA%mass(i)
         j = swiftest_plA%name(i)
         call orbel_xv2el(swiftest_plA%xh(:,i), swiftest_plA%vh(:,i), mu, a, e, inc, capom, omega, capm)
         call io_write_line(iu, j, a, e, inc, capom, omega, capm, out_type, &
         mass = swiftest_plA%mass(i),radius = swiftest_plA%radius(i))
      end do
      mu = swiftest_plA%mass(1)
      do i = 1, swiftest_tpA%ntp
         j = swiftest_tpA%name(i)
         call orbel_xv2el(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i), mu, a, e, inc, capom, omega, capm)
         call io_write_line(iu, j, a, e, inc, capom, omega, capm, out_type)
      end do
   case (XV)
      do i = 2, swiftest_plA%npl
         xtmp(:) = swiftest_plA%xh(:,i)
         vtmp(:) = swiftest_plA%vh(:,i)
         j = swiftest_plA%name(i)
         call io_write_line(iu, j, xtmp(1), xtmp(2), xtmp(3), vtmp(1), vtmp(2), vtmp(3), out_type,                     &
         mass = swiftest_plA%mass(i), radius = swiftest_plA%radius(i))
      end do
      do i = 1, swiftest_tpA%ntp
         xtmp(:) = swiftest_tpA%xh(:,i)
         vtmp(:) = swiftest_tpA%vh(:,i)
         j = swiftest_tpA%name(i)
         call io_write_line(iu, j, xtmp(1), xtmp(2), xtmp(3), vtmp(1), vtmp(2), vtmp(3), out_type)
      end do
   end select

   close(unit = iu, iostat = ierr)
   if (ierr /= 0) then
      write(*, *) "Swiftest error:"
      write(*, *) "   unable to close binary output file"
      call util_exit(FAILURE)
   end if

   return

   end procedure io_write_frame
end submodule s_io_write_frame
