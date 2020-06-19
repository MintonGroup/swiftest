submodule (nbody_data_structures) s_nbody_set_vec
contains
   module procedure nbody_set_vec
   !! author: David A. Minton
   !!
   !! Converts certain scalar values to arrays so that they can be used in elemental functions
   use swiftest
   implicit none

   self%mu_vec(:) = mu
   self%dt_vec(:) = dt

   return
   end procedure nbody_set_vec

end submodule s_nbody_set_vec

