submodule (util) s_util_exit
contains
   module procedure util_exit
   !! author: David A. Minton
   !!
   !! Print termination message and exit program
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: util_exit.f90
   !! Adapted from Hal Levison's Swift routine util_exit.f
   use swiftest
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

   end procedure util_exit
end submodule s_util_exit
