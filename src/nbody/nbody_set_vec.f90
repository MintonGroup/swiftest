submodule (nbody_data_structures) s_nbody_set_vec
contains

   module procedure nbody_set_vec_pl
      !! author: David A. Minton
      !!
      !! Converts certain scalar values to arrays so that they can be used in elemental functions
      use swiftest
      implicit none

      self%mu_vec(:) = self%mu + self%mass(:)
      self%dt_vec(:) = dt

      return
   end procedure nbody_set_vec_pl

   module procedure nbody_set_vec_tp
      !! author: David A. Minton
      !!
      !! Converts certain scalar values to arrays so that they can be used in elemental functions
      use swiftest
      implicit none

      self%mu_vec(:) = self%mu
      self%dt_vec(:) = dt

      return
   end procedure nbody_set_vec_tp


end submodule s_nbody_set_vec

