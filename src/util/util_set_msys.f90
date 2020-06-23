submodule (swiftest_classes) s_util_set_msys
contains
   module procedure util_set_msys
      !! author: David A. Minton
      !!
      !! Calculates the total system mass for a generic Swiftest system
      use swiftest
      implicit none

      self%msys = sum(swiftest_plA%mass(1:swiftest_plA%nbody))

      return
   end procedure util_set_msys

end submodule s_util_set_msys

