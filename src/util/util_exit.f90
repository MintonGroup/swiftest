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
      end select

      stop

   end subroutine util_exit
   
end submodule s_util_exit
