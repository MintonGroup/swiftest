submodule (symba_classes) s_symba_io
   use swiftest
contains
   module subroutine symba_io_dump_particle_info(self, param, msg) 
      !! author: David A. Minton
      !!
      !! Dumps the particle information data to a file
      implicit none
      class(symba_particle_info), intent(inout) :: self  !! Swiftest base object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      character(*), optional,     intent(in)    :: msg   !! Message to display with dump operation
   end subroutine symba_io_dump_particle_info

   module subroutine symba_io_initialize_particle_info(self, param) 
      !! author: David A. Minton
      !!
      !! Initializes a particle info data structure, either starting a new one or reading one in 
      !! from a file if it is a restarted run
      implicit none
      class(symba_particle_info), intent(inout) :: self  !! SyMBA particle info object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
   end subroutine symba_io_initialize_particle_info

   module subroutine symba_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in parameters specific to the SyMBA integrator, then calls the base io_param_reader.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_param.f90
      !! Adapted from Martin Duncan's Swift routine io_init_param.f
      implicit none
      ! Arguments
      class(symba_parameters), intent(inout) :: self       !! Collection of parameters
      integer,                 intent(in)    :: unit       !! File unit number
      character(len=*),        intent(in)    :: iotype     !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                           !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
      integer,                 intent(in)    :: v_list(:)  !! The first element passes the integrator code to the reader
      integer,                 intent(out)   :: iostat     !! IO status code
      character(len=*),        intent(inout) :: iomsg      !! Message to pass if iostat /= 0
      ! internals
      integer(I4B)                   :: ilength, ifirst, ilast  !! Variables used to parse input file
      character(STRMAX)              :: line                    !! Line of the input file
      character (len=:), allocatable :: line_trim,param_name, param_value !! Strings used to parse the param file
      integer(I4B) :: nseeds, nseeds_from_file, i
      logical                 :: seed_set = .false.      !! Is the random seed set in the input file?
      character(len=*),parameter    :: linefmt = '(A)'

      associate(param => self)
         call io_param_reader(param, unit, iotype, v_list, iostat, iomsg) 

         call random_seed(size = nseeds)
         if (allocated(param%seed)) deallocate(param%seed)
         allocate(param%seed(nseeds))
         do
            read(unit = unit, fmt = linefmt, iostat = iostat, end = 1) line
            line_trim = trim(adjustl(line))
            ilength = len(line_trim)
            if ((ilength /= 0)) then 
               ifirst = 1
               ! Read the pair of tokens. The first one is the parameter name, the second is the value.
               param_name = io_get_token(line_trim, ifirst, ilast, iostat)
               if (param_name == '') cycle ! No parameter name (usually because this line is commented out)
               call util_toupper(param_name)
               ifirst = ilast + 1
               param_value = io_get_token(line_trim, ifirst, ilast, iostat)
               select case (param_name)
               case ("FRAGMENTATION")
                  call util_toupper(param_value)
                  if (param_value == "YES" .or. param_value == "T") self%lfragmentation = .true.
               case ("ROTATION")
                  call util_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') self%lrotation = .true. 
               case ("TIDES")
                  call util_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') self%ltides = .true. 
               case ("MTINY")
                  read(param_value, *) param%mtiny
               case("SEED")
                  read(param_value, *) nseeds_from_file
                  ! Because the number of seeds can vary between compilers/systems, we need to make sure we can handle cases in which the input file has a different
                  ! number of seeds than the current system. If the number of seeds in the file is smaller than required, we will use them as a source to fill in the missing elements.
                  ! If the number of seeds in the file is larger than required, we will truncate the seed array.
                  if (nseeds_from_file > nseeds) then
                     nseeds = nseeds_from_file
                     deallocate(param%seed)
                     allocate(param%seed(nseeds))
                     do i = 1, nseeds
                        ifirst = ilast + 1
                        param_value = io_get_token(line, ifirst, ilast, iostat) 
                        read(param_value, *) param%seed(i)
                     end do
                  else ! Seed array in file is too small
                     do i = 1, nseeds_from_file
                        ifirst = ilast + 1
                        param_value = io_get_token(line, ifirst, ilast, iostat) 
                        read(param_value, *) param%seed(i)
                     end do
                     param%seed(nseeds_from_file+1:nseeds) = [(param%seed(1) - param%seed(nseeds_from_file) + i, i=nseeds_from_file+1, nseeds)]
                  end if
                  seed_set = .true.
               end select
            end if
         end do
         1 continue

         write(*,*) "FRAGMENTATION = ", param%lfragmentation
         if (param%lfragmentation) then
            if (seed_set) then
               call random_seed(put = param%seed)
            else
               call random_seed(get = param%seed)
            end if
            write(*,*) "SEED: N,VAL    = ",size(param%seed), param%seed(:)
         end if
         write(*,*) "ROTATION      = ", param%lrotation
         write(*,*) "TIDES         = ", param%ltides

         if (self%mtiny < 0.0_DP) then
            write(iomsg,*) "MTINY invalid or not set: ", self%mtiny
            iostat = -1
            return
         else
            write(*,*) "MTINY          = ", self%mtiny   
         end if

         if (.not.self%lclose) then
            write(iomsg,*) 'This integrator requires CHK_CLOSE to be enabled.'
            iostat = -1
            return
         end if
      end associate
      return

   end subroutine symba_io_param_reader

   module subroutine symba_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
      !! author: David A. Minton
      !!
      !! Dump integration parameters specific to SyMBA to file and then call the base io_param_writer method.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_param.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_param.f
      implicit none
      ! Arguments
      class(symba_parameters),intent(in)    :: self      !! Collection of SyMBA parameters
      integer,                intent(in)    :: unit      !! File unit number
      character(len=*),       intent(in)    :: iotype    !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                         !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
      integer,                intent(in)    :: v_list(:) !! Not used in this procedure
      integer,                intent(out)   :: iostat    !! IO status code
      character(len=*),       intent(inout) :: iomsg     !! Message to pass if iostat /= 0
      ! Internals
      character(*),parameter :: Ifmt  = '(I0)'         !! Format label for integer values
      character(*),parameter :: Rfmt  = '(ES25.17)'    !! Format label for real values 
      character(*),parameter :: Rarrfmt  = '(3(ES25.17,1X))'    !! Format label for real values 
      character(*),parameter :: Lfmt  = '(L1)'         !! Format label for logical values 
      character(len=*), parameter :: Afmt = '(A25,1X,64(:,A25,1X))'
      character(256)          :: param_name, param_value
      type character_array
         character(25) :: value
      end type character_array
      type(character_array), dimension(:), allocatable :: param_array
      integer(I4B) :: i

      associate(param => self)
         call io_param_writer(param, unit, iotype, v_list, iostat, iomsg) 

         ! Special handling is required for writing the random number seed array as its size is not known until runtime
         ! For the "SEED" parameter line, the first value will be the size of the seed array and the rest will be the seed array elements
         write(param_name, Afmt) "PARTICLE_FILE"; write(param_value, Afmt) trim(adjustl(param%particle_file)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "MTINY"; write(param_value, Rfmt) param%mtiny; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "ROTATION"; write(param_value, Lfmt)  param%lrotation; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "TIDES"; write(param_value, Lfmt)  param%ltides; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "FRAGMENTATION"; write(param_value, Lfmt)  param%lfragmentation; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         if (param%lfragmentation) then
            write(param_name, Afmt) "SEED"
            if (allocated(param_array)) deallocate(param_array)
            allocate(param_array(0:size(param%seed)))
            write(param_array(0)%value, Ifmt) size(param%seed)
            do i = 1, size(param%seed)
               write(param_array(i)%value, Ifmt) param%seed(i)
            end do
            write(unit, Afmt, advance='no') adjustl(param_name), adjustl(param_array(0)%value)
            do i = 1, size(param%seed)
               if (i < size(param%seed)) then
                  write(unit, Afmt, advance='no') adjustl(param_array(i)%value)
               else
                  write(unit, Afmt) adjustl(param_array(i)%value)
               end if
            end do
         end if

         iostat = 0
      end associate

      return

   end subroutine symba_io_param_writer

   module subroutine symba_io_read_frame_cb(self, iu, param, form, ierr)
      !! author: David A. Minton
      !!
      !! Reads a frame of output of central body data to the binary output file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_read_frame.f90
      !! Adapted from Hal Levison's Swift routine io_read_frame.F
      implicit none
      ! Arguments
      class(symba_cb),            intent(inout) :: self     !! Swiftest central body object
      integer(I4B),               intent(inout) :: iu       !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(inout) :: param   !! Current run configuration parameters 
      character(*),               intent(in)    :: form     !! Input format code ("XV" or "EL")
      integer(I4B),               intent(out)   :: ierr     !! Error code

      call io_read_frame_cb(self, iu, param, form, ierr)
      select type(param)
      class is (symba_parameters)
         if (param%lrotation) then
            read(iu, iostat = ierr) self%Ip(:)
            read(iu, iostat = ierr) self%rot(:)
         end if
         if (param%ltides) then
            read(iu, iostat = ierr) self%k2
            read(iu, iostat = ierr) self%Q
         end if
      end select
      if (ierr /=0) then
         write(*,*) 'Error reading SyMBA central body data'
         call util_exit(FAILURE)
      end if
      return
   end subroutine symba_io_read_frame_cb

   module subroutine symba_io_read_frame_pl(self, iu, param, form, ierr)
      !! author: David A. Minton
      !!
      !! Reads a frame of output of a SyMBA massive body object 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_read_frame.f90
      !! Adapted from Hal Levison's Swift routine io_read_frame.F
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout) :: self    !! Swiftest particle object
      integer(I4B),               intent(inout) :: iu      !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      character(*),               intent(in)    :: form    !! Input format code ("XV" or "EL")
      integer(I4B),               intent(out)   :: ierr    !! Error code

      call io_read_frame_body(self, iu, param, form, ierr)
      select type(param)
      class is (symba_parameters)
         associate(pl => self, npl => self%nbody)
            if (param%lrotation) then
               read(iu, iostat = ierr) pl%rot(1, 1:npl)
               read(iu, iostat = ierr) pl%rot(2, 1:npl)
               read(iu, iostat = ierr) pl%rot(3, 1:npl)
               read(iu, iostat = ierr) pl%Ip(1, 1:npl)
               read(iu, iostat = ierr) pl%Ip(2, 1:npl)
               read(iu, iostat = ierr) pl%Ip(3, 1:npl)
            end if
            if (param%ltides) then
               read(iu, iostat = ierr) pl%k2(1:npl)
               read(iu, iostat = ierr) pl%Q(1:npl)
            end if
         end associate
      end select

      if (ierr /=0) then
         write(*,*) 'Error reading SyMBA massive body body data'
         call util_exit(FAILURE)
      end if
      return
   end subroutine symba_io_read_frame_pl

   module subroutine symba_io_read_frame_info(self, iu, param, form, ierr)
      !! author: David A. Minton
      !!
      !! Reads a single frame of a particle info data from a file.
      implicit none
      class(symba_particle_info), intent(inout) :: self  !! SyMBA particle info object
      integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      character(*),               intent(in)    :: form  !! Input format code ("XV" or "EL")
      integer(I4B),               intent(out)   :: ierr  !! Error code

      ierr = 0
   end subroutine symba_io_read_frame_info

   module subroutine symba_io_read_cb_in(self, param) 
      !! author: David A. Minton
      !!
      !! Reads in central body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f
      implicit none
      ! Arguments
      class(symba_cb),            intent(inout) :: self
      class(swiftest_parameters), intent(inout) :: param
      ! Internals
      integer(I4B), parameter :: LUN = 7              !! Unit number of input file
      integer(I4B)            :: iu = LUN
      integer(I4B)            :: ierr
      logical                 :: is_ascii 
      real(DP)                :: t
      real(QP)                :: val

      select type(param)
      class is (symba_parameters)
         ierr = 0
         is_ascii = (param%in_type == 'ASCII') 
         if (is_ascii) then
            open(unit = iu, file = param%incbfile, status = 'old', form = 'FORMATTED', iostat = ierr)
            read(iu, *, iostat = ierr) val 
            self%Gmass = real(val, kind=DP)
            self%mass = real(val / param%GU, kind=DP)
            read(iu, *, iostat = ierr) self%radius
            read(iu, *, iostat = ierr) self%j2rp2
            read(iu, *, iostat = ierr) self%j4rp4
            if (param%lrotation) then
               read(iu, *, iostat = ierr) self%Ip(:)
               read(iu, *, iostat = ierr) self%rot(:)
            end if
            if (param%ltides) then
               read(iu, *, iostat = ierr) self%k2
               read(iu, *, iostat = ierr) self%Q
            end if
         else
            open(unit = iu, file = param%incbfile, status = 'old', form = 'UNFORMATTED', iostat = ierr)
            call self%read_frame(iu, param, XV, ierr)
         end if
         close(iu)
         if (ierr /=  0) then
            write(*,*) 'Error opening massive body initial conditions file ',trim(adjustl(param%incbfile))
            call util_exit(FAILURE)
         end if
         if (self%j2rp2 /= 0.0_DP) param%loblatecb = .true.
      end select

      return
   end subroutine symba_io_read_cb_in

   module subroutine symba_io_read_pl_in(self, param) 
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in either test particle or massive body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90 and swiftest_init_tp.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f and swiftest_init_tp.f
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout) :: self   !! Swiftest particle object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B), parameter       :: LUN = 7              !! Unit number of input file
      integer(I4B)                  :: iu = LUN
      integer(I4B)                  :: i, ierr, nbody
      logical                       :: is_ascii
      character(len=:), allocatable :: infile
      real(DP)                      :: t
      real(QP)                      :: val

      select type(param)
      class is (symba_parameters)
         ierr = 0
         is_ascii = (param%in_type == 'ASCII') 
         select case(param%in_type)
         case(ASCII_TYPE)
            open(unit = iu, file = infile, status = 'old', form = 'FORMATTED', iostat = ierr)
            read(iu, *, iostat = ierr) nbody
            call self%setup(nbody)
            if (nbody > 0) then
               do i = 1, nbody
                  read(iu, *, iostat = ierr) self%id(i), val 
                  self%mass(i) = real(val / param%GU, kind=DP)
                  self%Gmass(i) = real(val, kind=DP)
                  read(iu, *, iostat = ierr) self%radius(i)
                  if (param%lrotation) then
                     read(iu, iostat = ierr) self%Ip(:, i)
                     read(iu, iostat = ierr) self%rot(:, i)
                  end if
                  if (param%ltides) then
                     read(iu, iostat = ierr) self%k2(i)
                     read(iu, iostat = ierr) self%Q(i)
                  end if
                  if (ierr /= 0 ) exit
                  read(iu, *, iostat = ierr) self%xh(1, i), self%xh(2, i), self%xh(3, i)
                  read(iu, *, iostat = ierr) self%vh(1, i), self%vh(2, i), self%vh(3, i)
                  if (ierr /= 0 ) exit
                  self%status(i) = ACTIVE
               end do
            end if
         case (REAL4_TYPE, REAL8_TYPE)  
            open(unit = iu, file = infile, status = 'old', form = 'UNFORMATTED', iostat = ierr)
            read(iu, iostat = ierr) nbody
            call self%setup(nbody)
            if (nbody > 0) then
               call self%read_frame(iu, param, XV, ierr)
               self%status(:) = ACTIVE
            end if
         case default
            write(*,*) trim(adjustl(param%in_type)) // ' is an unrecognized file type'
            ierr = -1
         end select
         close(iu)
         if (ierr /= 0 ) then
            write(*,*) 'Error reading in initial conditions from ',trim(adjustl(infile))
            call util_exit(FAILURE)
         end if
      end select

      return
   end subroutine symba_io_read_pl_in

   module subroutine symba_io_write_frame_cb(self, iu, param)
      !! author: David A. Minton
      !!
      !! Writes a single frame of a SyMBA pl file 
      implicit none
      class(symba_cb),            intent(in)    :: self  !! SyMBA massive body object
      integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 

      call io_write_frame_cb(self, iu, param)
      select type(param)
      class is (symba_parameters)
         if (param%lrotation) then
            write(iu) self%rot(1)
            write(iu) self%rot(2)
            write(iu) self%rot(3)
            write(iu) self%Ip(1)
            write(iu) self%Ip(2)
            write(iu) self%Ip(3)
         end if
         if (param%ltides) then
            write(iu) self%k2
            write(iu) self%Q
         end if
      end select
   end subroutine symba_io_write_frame_cb

   module subroutine symba_io_write_frame_pl(self, iu, param)
      !! author: David A. Minton
      !!
      !! Writes a single frame of a SyMBA pl file 
      implicit none
      class(symba_pl),            intent(in)    :: self  !! SyMBA massive body object
      integer(I4B),               intent(inout) :: iu    !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 

      call io_write_frame_body(self, iu, param)
      select type(param)
      class is (symba_parameters)
         associate(pl => self, npl => self%nbody)
            if (param%lrotation) then
               write(iu) pl%rot(1, 1:npl)
               write(iu) pl%rot(2, 1:npl)
               write(iu) pl%rot(3, 1:npl)
               write(iu) pl%Ip(1, 1:npl)
               write(iu) pl%Ip(2, 1:npl)
               write(iu) pl%Ip(3, 1:npl)
            end if
            if (param%ltides) then
               write(iu) pl%k2(1:npl)
               write(iu) pl%Q(1:npl)
            end if
         end associate
      end select
   end subroutine symba_io_write_frame_pl

end submodule s_symba_io

