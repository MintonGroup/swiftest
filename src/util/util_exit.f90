!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_exit
   use swiftest
contains

   module subroutine util_exit(code)
      !! author: David A. Minton
      !!
      !! Print termination message and exit program
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_exit.f90
      !! Adapted from Hal Levison's Swift routine util_exit.f
      implicit none
      ! Arguments
      integer(I4B), intent(in) :: code
      ! Internals
      character(*), parameter :: BAR = '("------------------------------------------------")'

      select case(code)
      case(SUCCESS)
         write(*, SUCCESS_MSG) VERSION_NUMBER
         write(*, BAR)
      case(USAGE) 
         write(*, USAGE_MSG)
      case(HELP)
         write(*, HELP_MSG)
      case default
         write(*, FAIL_MSG) VERSION_NUMBER
         write(*, BAR)
         error stop
      end select

      stop

   end subroutine util_exit
   
end submodule s_util_exit
