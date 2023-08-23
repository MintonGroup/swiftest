!! Copyright 2023 - David Minton
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

program main
   !! author: David A. Minton
   !!
   !! This is used as a wrapper for the swiftest_driver subroutine so that the driver can be used as either
   !! a function call or a standalone executable program.
   use swiftest
   implicit none

   character(len=:), allocatable             :: integrator        !! Integrator type code (see globals for symbolic names)
   character(len=:), allocatable             :: param_file_name   !! Name of the file containing user-defined parameters
   character(len=:), allocatable             :: display_style     !! Style of the output display {"STANDARD", "COMPACT", "PROGRESS"}). Default is "STANDARD"

   ! Read in command line arguments
   call swiftest_io_get_args(integrator, param_file_name, display_style, from_cli=.true.)

   ! Execute the driver
   call swiftest_driver(integrator, param_file_name, display_style)


#ifdef COARRAY
   if (this_image() == 1) then
#endif
      call base_util_exit(SUCCESS)
#ifdef COARRAY
   end if ! (this_image() == 1) 
#endif

end program main