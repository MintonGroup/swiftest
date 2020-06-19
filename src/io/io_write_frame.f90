submodule (nbody_data_structures) s_io_write_frame
contains
   module procedure io_write_frame
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Write a frame (header plus records for each massive body and active test particle) to output binary file
   !! There is no direct file output from this subroutine
   !!
   !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
   !! Adapted from Hal Levison's Swift routine io_write_frame.F
   use swiftest
   implicit none
   

   logical, save             :: lfirst = .true.
   integer(I4B), parameter   :: lun = 20
   integer(I4B)              :: i, j, ierr
   integer(I4B), save        :: iu = lun, iout_form = XV
   real(DP),dimension(:),allocatable :: a, e, inc, capom, omega, capm
   real(DP), dimension(NDIM) :: xtmp, vtmp

   if (lfirst) then
      select case(config%out_stat)
      case('APPEND')
         open(unit = iu, file = config%outfile, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr)
      case('NEW')
         open(unit = iu, file = config%outfile, status = 'NEW', form = 'UNFORMATTED', iostat = ierr)
      case ('REPLACE')
         open(unit = iu, file = config%outfile, status = 'REPLACE', form = 'UNFORMATTED', iostat = ierr)
      case default
         write(*,*) 'Invalid status code',trim(adjustl(config%out_stat))
         call util_exit(FAILURE)
      end select
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   binary output file already exists or cannot be accessed"
         call util_exit(FAILURE)
      end if

      select case (config%out_form)
      case ("EL")
         iout_form = EL
      case ("XV")
         iout_form = XV
      end select
      lfirst = .false.
   else
      open(unit = iu, file = config%outfile, status = 'OLD', position =  'APPEND', form = 'UNFORMATTED', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   unable to open binary output file for APPEND"
         call util_exit(FAILURE)
      end if
   end if

   call io_write_hdr(iu, t, swiftest_plA%nbody, swiftest_tpA%nbody, iout_form, config%out_type)
   select case (iout_form)
   case (EL)
      associate (npl => swiftest_plA%npl, mu => swiftest_plA%mu_vec)
         mu(2:npl) = swiftest_plA%mass(1) + swiftest_plA%mass(2:npl)
         call orbel_xv2el(mu(2:npl), swiftest_plA%xh(1,2:npl), &
                                     swiftest_plA%xh(2,2:npl), &
                                     swiftest_plA%xh(3,2:npl), &
                                     swiftest_plA%vh(1,2:npl), &
                                     swiftest_plA%vh(2,2:npl), &
                                     swiftest_plA%vh(3,2:npl), &
                                     a, e, inc, capom, omega, capm)
      end associate
         !call io_write_line(iu, j, a, e, inc, capom, omega, capm, config%out_type, &
         mass = swiftest_plA%mass(i),radius = swiftest_plA%radius(i)
      end do
      mu = swiftest_plA%mass(1)
      do i = 1, swiftest_tpA%nbody
         j = swiftest_tpA%name(i)
         call orbel_xv2el(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i), mu, a, e, inc, capom, omega, capm)
         !call io_write_line(iu, j, a, e, inc, capom, omega, capm, config%out_type)
      end do
   case (XV)
      do i = 2, swiftest_plA%nbody
         xtmp(:) = swiftest_plA%xh(:,i)
         vtmp(:) = swiftest_plA%vh(:,i)
         j = swiftest_plA%name(i)
         call io_write_line(iu, j, xtmp(1), xtmp(2), xtmp(3), vtmp(1), vtmp(2), vtmp(3), config%out_type,                     &
         mass = swiftest_plA%mass(i), radius = swiftest_plA%radius(i))
      end do
      do i = 1, swiftest_tpA%nbody
         xtmp(:) = swiftest_tpA%xh(:,i)
         vtmp(:) = swiftest_tpA%vh(:,i)
         j = swiftest_tpA%name(i)
         call io_write_line(iu, j, xtmp(1), xtmp(2), xtmp(3), vtmp(1), vtmp(2), vtmp(3), config%out_type)
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
