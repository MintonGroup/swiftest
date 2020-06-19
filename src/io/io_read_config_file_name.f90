submodule (nbody_data_structures) s_io_read_config_file_name

contains

   module procedure io_read_config_file_name
      !! author: David A. Minton
      !!
      !! Reads in the name of the configuration file. It first tries the command line arguments, and if
      !!    none are present, it will prompt the user to input the file name manually
      use swiftest
      implicit none

      character(len=STRMAX)) :: arg
      integer :: ierr,i,narg
      character(len=*),parameter    :: linefmt = '(A)'

      narg = command_argument_count() !
      if (narg > 0) then
         call get_command_argument(1, arg, status=ierr)
         if (ierr == 0) then
            config_file_name = trim(adjustl(arg))
         else
            write(*,*) 'Invalid command line argument passed to Swiftest'
            call util_exit(FAILURE)
         end if
      else
         write(*, linefmt, advance = "no") "Enter name of configuration data file: "
         read(*, linefmt) arg
         config_file_name = trim(adjustl(arg))
      end if
   end procedure io_read_config_file_name

end submodule s_io_read_config_file_name
