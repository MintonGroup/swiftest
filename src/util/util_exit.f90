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

   if (code == SUCCESS) then
      write(*, 100) VERSION_NUMBER
   else
      write(*, 200) VERSION_NUMBER
   end if
 100 format(/, "Normal termination of swiftest (version ", f3.1, ")")
 200 format(/, "Terminating Swiftest (version ", f3.1, ") due to error!!")
   write(*, 300) "------------------------------------------------"
 300 format(a)

   stop

   end procedure util_exit
end submodule s_util_exit
