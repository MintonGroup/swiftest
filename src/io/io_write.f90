submodule (swiftest_classes) io_write
   !! author: David A. Minton
   !! 
   !! This submodule contains implementations of the following procedures:
   !!    io_write_frame
   !!    io_write_config_in
   !!    io_config_reader
   !!    io_write_cb_in
   !!    io_write_pl_in
   !!    io_write_tp_in
   !!    io_write_line_swifter
contains
   module procedure io_write_frame_system
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
      integer(I4B)              :: i, j, ierr
      real(DP),dimension(:),allocatable :: a, e, inc, capom, omega, capm
      real(DP), dimension(NDIM) :: xtmp, vtmp

      iu = BINUNIT
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
            write(*, *) "   Binary output file already exists or cannot be accessed"
            call util_exit(FAILURE)
         end if
      else
         open(unit = iu, file = config%outfile, status = 'OLD', position =  'APPEND', form = 'UNFORMATTED', iostat = ierr)
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Unable to open binary output file for APPEND"
            call util_exit(FAILURE)
         end if
      end if
      call io_write_hdr(iu, t, self%pl%nbody, self%tp%nbody, config%out_form, config%out_type)
      if (config%out_form == EL) then ! Do an orbital element conversion prior to writing out the frame, as we have access to the central body here
         call self%pl%xv2el(self%cb)
         call self%tp%xv2el(self%cb)
      end if
      call self%cb%write_frame(iu, config, t, dt)
      call self%pl%write_frame(iu, config, t, dt)
      call self%tp%write_frame(iu, config, t, dt)

      return
   end procedure io_write_frame_system

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
            write(iu, iostat = ierr) real(t, kind=SP), npl, ntp, out_form
            if (ierr /= 0) then
               write(*, *) "Swiftest error:"
               write(*, *) "   Unable to write binary file header"
               call util_exit(FAILURE)
            end if
         case (REAL8_TYPE,SWIFTER_REAL8_TYPE)
            write(iu, iostat = ierr) t, npl, ntp, out_form
            if (ierr /= 0) then
               write(*, *) "Swiftest error:"
               write(*, *) "   Unable to write binary file header"
               call util_exit(FAILURE)
            end if
      end select
   
      return
   
   end procedure io_write_hdr

   module procedure io_write_frame_cb
      !! author: David A. Minton
      !!
      !! Write a frame of output of central body data to the binary output file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      use swiftest
      implicit none

      write(iu) self%mass
      write(iu) self%radius
      write(iu) self%j2rp2 
      write(iu) self%j4rp4 
      if (config%lrotation) then
         write(iu) self%Ip(:)
         write(iu) self%rot(:)
      end if
      if (config%ltides) then
         write(iu) self%k2
         write(iu) self%Q
      end if

      return
   end procedure io_write_frame_cb

   module procedure io_write_frame_body
      !! author: David A. Minton
      !!
      !! Write a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      use swiftest
      implicit none

      associate(n => self%nbody)
         select case (config%out_form)
         case (EL) 
            write(iu) self%a(1:n)
            write(iu) self%e(1:n)
            write(iu) self%inc(1:n)
            write(iu) self%capom(:)
            write(iu) self%omega(:)
            write(iu) self%capm(:)
         case (XV)
            write(iu) self%xh(1,1:n)
            write(iu) self%xh(2,1:n)
            write(iu) self%xh(3,1:n)
            write(iu) self%vh(1,1:n)
            write(iu) self%vh(2,1:n)
            write(iu) self%vh(3,1:n)
         end select
         select type(self)  
         class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
            write(iu) self%mass(1:n)
            write(iu) self%radius(1:n)
            if (config%lrotation) then
               write(iu) self%Ip(1,1:n)
               write(iu) self%Ip(2,1:n)
               write(iu) self%Ip(3,1:n)
               write(iu) self%rot(1,1:n)
               write(iu) self%rot(2,1:n)
               write(iu) self%rot(3,1:n)
            end if
            if (config%ltides) then
               write(iu) self%k2(1:n)
               write(iu) self%Q(1:n)
            end if
         end select
      end associate

      return
   end procedure io_write_frame_body

   module procedure io_write_encounter
      !! author: David A. Minton
      !!
      !! Write close encounter data to output binary files
      !!  There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_write_encounter.f90
      !! Adapted from Hal Levison's Swift routine io_write_encounter.f
      use swiftest
      implicit none
      logical         :: lxdr
      logical , save    :: lfirst = .true.
      integer(I4B), parameter :: lun = 30
      integer(I4B)        :: ierr
      integer(I4B), save    :: iu = lun

      open(unit = iu, file = encounter_file, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED')
      if ((ierr /= 0) .and. lfirst) then
         open(unit = iu, file = encounter_file, status = 'NEW', form = 'UNFORMATTED')
      end if
      if (ierr /= 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to open binary encounter file"
         call util_exit(FAILURE)
      end if
      lfirst = .false.
      write(iu, iostat = ierr) t
      if (ierr < 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to write binary file record"
         call util_exit(FAILURE)
      end if
      write(iu) name1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), mass1, radius1
      write(iu) name2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), mass2, radius2
      close(unit = iu, iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to close binary encounter file"
         call util_exit(FAILURE)
      end if

      return

   end procedure io_write_encounter

end submodule io_write
