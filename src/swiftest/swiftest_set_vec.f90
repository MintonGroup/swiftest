submodule (swiftest_data_structures) s_swiftest_set_vec
contains
   module procedure swiftest_set_vec
   !! author: David A. Minton
   !!
   !! Converts certain scalar values to arrays so that they can be used in elemental functions
   implicit none

   self%mu_vec(:) = mu
   self%dt_vec(:) = dt

   return
   end procedure swiftest_set_vec

end submodule s_swiftest_set_vec

