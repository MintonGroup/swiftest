submodule (swiftest_classes) s_setup_set_msys
contains
   module procedure setup_set_msys
      !! author: David A. Minton
      !!
      !! Calculates the total system mass for a generic Swiftest system
      use swiftest
      implicit none
   
      self%msys = self%cb%mass + sum(self%pl%mass(1:self%pl%nbody))
   
      return
      end procedure setup_set_msys


end submodule s_setup_set_msys
