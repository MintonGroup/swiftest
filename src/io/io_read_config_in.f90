submodule (swiftest_data_structures) s_io_read_config_in
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
   write(*, *) 'Configuration data file is ', trim(adjustl(inparfile))
   write(*, *) ' '
   100 format(A)
   open(unit = LUN, file = inparfile, status = 'old', iostat = ierr)
   if (ierr /= 0) then
      write(*,*) 'Swiftest error: ', ierr
      write(*,*) '   Unable to open file ',trim(adjustl(inparfile))
      call util_exit(FAILURE)
   end if

   !! todo: Currently this procedure does not work in user-defined derived-type input mode 
   !!    as the newline characters are ignored in the input file when compiled in ifort.

   !read(LUN,'(DT)', iostat= ierr, iomsg = error_message) config
   call config%udio_reader(LUN,iotype="none",v_list=(/0/),iostat=ierr,iomsg=error_message)
   if (ierr /= 0) then
      write(*,*) 'Swiftest error reading ', trim(adjustl(inparfile))
      write(*,*) ierr,trim(adjustl(error_message))
      call util_exit(FAILURE)
   end if

   close(LUN)

   ierr = 0

   if (ierr < 0) then
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

end submodule s_io_read_config_in
