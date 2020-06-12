submodule (swiftest_data_structures) s_swiftest_set_mu
contains
   module procedure swiftest_set_mu
   !! author: David A. Minton
   !!
   !! Saves the central body mass term mu in vector form to be used in elemental functions
   implicit none

   self%mu(:) = mu

   return
   end procedure swiftest_set_mu

end submodule s_swiftest_set_mu

