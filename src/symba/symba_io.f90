submodule (symba_classes) s_symba_io
   use swiftest
contains

   module subroutine symba_io_dump_particle_info(self, param, msg) 
      !! author: David A. Minton
      !!
      !! Dumps the particle information data to a file
      implicit none
      class(symba_particle_info), intent(inout) :: self  !! Swiftest base object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
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
         rewind(unit)
         do
            read(unit = unit, fmt = linefmt, iostat = iostat, end = 1) line
            line_trim = trim(adjustl(line))
            ilength = len(line_trim)
            if ((ilength /= 0)) then 
               ifirst = 1
               ! Read the pair of tokens. The first one is the parameter name, the second is the value.
               param_name = io_get_token(line_trim, ifirst, ilast, iostat)
               if (param_name == '') cycle ! No parameter name (usually because this line is commented out)
               call io_toupper(param_name)
               ifirst = ilast + 1
               param_value = io_get_token(line_trim, ifirst, ilast, iostat)
               select case (param_name)
               case ("FRAGMENTATION")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == "T") self%lfragmentation = .true.
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

      iostat = 0

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


   module subroutine symba_io_write_discard(self, param)
      implicit none
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      ! Internals
      integer(I4B), parameter   :: LUN = 40
      integer(I4B)          :: iadd, isub, j, ierr, nsub, nadd
      logical, save :: lfirst = .true. 
      real(DP), dimension(:,:), allocatable :: vh
      character(*), parameter :: HDRFMT    = '(E23.16, 1X, I8, 1X, L1)'
      character(*), parameter :: NAMEFMT   = '(A, 2(1X, I8))'
      character(*), parameter :: VECFMT    = '(3(E23.16, 1X))'
      character(*), parameter :: NPLFMT    = '(I8)'
      character(*), parameter :: PLNAMEFMT = '(I8, 2(1X, E23.16))'
      class(swiftest_body), allocatable :: pltemp

      associate(pl => self%pl, npl => self%pl%nbody, mergesub_list => self%mergesub_list, mergeadd_list => self%mergeadd_list)
         if (self%tp_discards%nbody > 0) call io_write_discard(self, param)

         if (mergesub_list%nbody == 0) return
         select case(param%out_stat)
         case('APPEND')
            open(unit = LUN, file = param%discard_out, status = 'OLD', position = 'APPEND', form = 'FORMATTED', iostat = ierr)
         case('NEW', 'REPLACE', 'UNKNOWN')
            open(unit = LUN, file = param%discard_out, status = param%out_stat, form = 'FORMATTED', iostat = ierr)
         case default
            write(*,*) 'Invalid status code for OUT_STAT: ',trim(adjustl(param%out_stat))
            call util_exit(FAILURE)
         end select
         lfirst = .false.
         if (param%lgr) then
            call mergesub_list%pv2v(param) 
            call mergeadd_list%pv2v(param) 
         end if

         write(LUN, HDRFMT) param%t, mergesub_list%nbody, param%lbig_discard
         iadd = 1
         isub = 1
         do while (iadd <= mergeadd_list%nbody)
            nadd = mergeadd_list%ncomp(iadd)
            nsub = mergesub_list%ncomp(isub)
            do j = 1, nadd
               if (iadd <= mergeadd_list%nbody) then
                  write(LUN, NAMEFMT) ADD, mergesub_list%id(iadd), mergesub_list%status(iadd)
                  write(LUN, VECFMT) mergeadd_list%xh(1, iadd), mergeadd_list%xh(2, iadd), mergeadd_list%xh(3, iadd)
                  write(LUN, VECFMT) mergeadd_list%vh(1, iadd), mergeadd_list%vh(2, iadd), mergeadd_list%vh(3, iadd)
               else 
                  exit
               end if
               iadd = iadd + 1
            end do
            do j = 1, nsub
               if (isub <= mergesub_list%nbody) then
                  write(LUN, NAMEFMT) SUB, mergesub_list%id(isub), mergesub_list%status(isub)
                  write(LUN, VECFMT) mergesub_list%xh(1, isub), mergesub_list%xh(2, isub), mergesub_list%xh(3, isub)
                  write(LUN, VECFMT) mergesub_list%vh(1, isub), mergesub_list%vh(2, isub), mergesub_list%vh(3, isub)
               else
                  exit
               end if
               isub = isub + 1
            end do
         end do

         close(LUN)
      end associate

      return
   end subroutine symba_io_write_discard


   module subroutine symba_io_write_frame_info(self, iu, param)
      implicit none
      class(symba_particle_info), intent(in)    :: self  !! SyMBA particle info object
      integer(I4B),               intent(inout) :: iu      !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(in)    :: param   !! Current run configuration parameters 
   end subroutine symba_io_write_frame_info 

end submodule s_symba_io

