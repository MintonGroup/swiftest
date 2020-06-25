submodule (swiftest_classes) s_io_read_config_in
contains
   module procedure io_read_config_in
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in parameters for the integration
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_config.f90
      !! Adapted from Martin Duncan's Swift routine io_init_config.f
      !$ use omp_lib
      implicit none

      integer(I4B), parameter :: LUN = 7                 !! Unit number of input file
      integer(I4B)            :: ierr = 0                !! Input error code
      character(STRMAX)       :: error_message           !! Error message in UDIO procedure

      ! Read in name of parameter file
      write(*, *) 'Configuration data file is ', trim(adjustl(config_file_name))
      write(*, *) ' '
      100 format(A)
      open(unit = LUN, file = config_file_name, status = 'old', iostat = ierr)
      if (ierr /= 0) then
         write(*,*) 'Swiftest error: ', ierr
         write(*,*) '   Unable to open file ',trim(adjustl(config_file_name))
         call util_exit(FAILURE)
      end if

      config%integrator = integrator

      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    as the newline characters are ignored in the input file when compiled in ifort.

      !read(LUN,'(DT)', iostat= ierr, iomsg = error_message) config
      call config_reader(LUN,iotype="none",v_list=(/0/),iostat=ierr,iomsg=error_message)
      if (ierr /= 0) then
         write(*,*) 'Swiftest error reading ', trim(adjustl(config_file_name))
         write(*,*) ierr,trim(adjustl(error_message))
         call util_exit(FAILURE)
      end if

      close(LUN)

      ierr = 0
      select case(integrator)
      case(SYMBA)
         config%lmtiny = .true. ! SyMBA requires the mtiny value
         if (.not. config%lrhill_present) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Integrator SyMBA requires massive body Hill sphere radii on input"
            ierr = -1 
         end if
      case default
         ierr = 0
      end select

      if (ierr /= 0) then
         write(*, 100) "Input parameter(s) failed check"
         call util_exit(FAILURE)
      end if

      !> Define the maximum number of threads
      nthreads = 1            ! In the *serial* case
      !$ nthreads = omp_get_max_threads() ! In the *parallel* case
      !$ write(*,'(a)')   ' OpenMP parameters:'
      !$ write(*,'(a)')   ' ------------------'
      !$ write(*,'(a,i3,/)') ' Number of threads  = ', nthreads   

      return 
   end procedure io_read_config_in

   module procedure io_read_cb_in
      !! author: David A. Minton
      !!
      !! Read in central body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f
      implicit none

      integer(I4B), parameter :: LUN = 7              !! Unit number of input file
      integer(I4B)            :: i, iu, ierr, npl
      logical                 :: is_ascii 

      ierr = 0
      is_ascii = (config%in_type == 'ASCII') 
      if (is_ascii) then
         open(unit = LUN, file = config%incbfile, status = 'old', form = 'formatted', iostat = ierr)
      else
         open(unit = LUN, file = config%incbfile, status = 'old', form = 'unformatted', iostat = ierr)
      end if
      if (ierr /=  0) then
         write(*,*) 'Error opening massive body initial conditions file ',trim(adjustl(config%inplfile))
         call util_exit(FAILURE)
      end if

      return
   end procedure io_read_cb_in

   module procedure io_read_pl_in
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in massive body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f
      implicit none

      integer(I4B), parameter :: LUN = 7              !! Unit number of input file
      integer(I4B)            :: i, iu, ierr, npl
      logical                 :: is_ascii 

      ierr = 0
      is_ascii = (config%in_type == 'ASCII') 
      if (is_ascii) then
         open(unit = LUN, file = config%inplfile, status = 'old', form = 'formatted', iostat = ierr)
      else
         open(unit = LUN, file = config%inplfile, status = 'old', form = 'unformatted', iostat = ierr)
      end if
      if (ierr /=  0) then
         write(*,*) 'Error opening massive body initial conditions file ',trim(adjustl(config%inplfile))
         return
      end if

      if (is_ascii) then
         read(LUN, *, iostat = ierr) npl
      else
         read(LUN, iostat = ierr) npl
      end if
      if (npl <= 0) return
      call self%alloc(npl)

      if (is_ascii) then
         read(LUN, *, iostat = ierr) self%name(1), self%mass(1)
         self%rhill(1) = 0.0_dp
         self%radius(1) = 0.0_dp
         read(LUN, *, iostat = ierr) self%xh(:,1)
         read(LUN, *, iostat = ierr) self%vh(:,1)
         if (ierr /= 0) then
            write(*,*) 'Error reading central body values in ',trim(adjustl(config%inplfile))
            return
         end if
         do i = 1, NDIM
            if ((self%xh(i,1) /= 0.0_dp) .or. (self%vh(i,1) /= 0.0_dp)) then
               write(*, *) "Swiftest error:"
               write(*, *) " Input must be in heliocentric coordinates."
               write(*, *) " position/velocity components of body 1 are"
               write(*, *) self%xh(:,1)
               write(*, *) self%vh(:,1)
            end if
         end do
         self%status(1) = ACTIVE
         do i = 2, self%nbody
            if (config%lrhill_present) then
               read(LUN, *, iostat = ierr) self%name(i), self%mass(i), self%rhill(i)
            else
               read(LUN, *, iostat = ierr) self%name(i), self%mass(i)
               self%rhill(i) = 0.0_dp
            end if
            if (ierr /= 0 ) exit
            if (config%lclose) then
               read(LUN, *, iostat = ierr) self%radius(i)
               if (ierr /= 0 ) exit
            else
               self%radius(i) = 0.0_dp
            end if
            read(LUN, *, iostat = ierr) self%xh(:,i)
            read(LUN, *, iostat = ierr) self%vh(:,i)
            if (ierr /= 0 ) exit
            self%status(i) = ACTIVE
         end do
      else
         read(LUN, iostat = ierr) self%name(:)
         read(LUN, iostat = ierr) self%mass(:)
         if (config%lrhill_present) then
            read(LUN, iostat = ierr) self%rhill(:)
         else
            self%rhill(:) = 0.0_dp
         end if
         self%status(:) = ACTIVE
         if (config%lclose) then
            read(LUN, iostat = ierr) self%radius(:)
         else
            self%radius(:) = 0.0_dp
         end if
         read(LUN, iostat = ierr) self%xh(:,:)
         read(LUN, iostat = ierr) self%vh(:,:)
         self%status(:) = ACTIVE
      end if
      close(unit = LUN)
      if (ierr /= 0 ) then
         write(*,*) 'Error reading in massive body initial conditions from ',trim(adjustl(config%inplfile))
         call util_exit(FAILURE)
      end if

      return
   end procedure io_read_pl_in

   module procedure io_read_tp_in
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in test particle data
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine io_init_tp.f
      implicit none
   
      integer(I4B), parameter  :: LUN = 7              !! Unit number of input file
      integer(I4B)             :: i, iu, ierr, ntp
      logical                  :: is_ascii
   
      ierr = 0
      is_ascii = (config%in_type == 'ASCII')  
      if (is_ascii) then
         open(unit = LUN, file = config%intpfile, status = 'old', form = 'formatted', iostat = ierr)
      else
         open(unit = LUN, file = config%intpfile, status = 'old', form = 'unformatted', iostat = ierr)
      end if
      if (ierr /=  0) then
         write(*,*) 'Error opening test particle initial conditions file ',trim(adjustl(config%intpfile))
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
   end procedure io_read_tp_in
end submodule s_io_read_config_in
