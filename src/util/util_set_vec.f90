submodule (swiftest_classes) s_util_set_vec
contains

   module procedure util_set_vec_dt
      !! author: David A. Minton
      !!
      !! Converts scalar dt value into vector for use in elemental functions
      use swiftest
      implicit none

      self%dt_vec(:) = dt

      return
   end procedure util_set_vec_dt

   module procedure util_set_vec_mu_pl
      !! author: David A. Minton
      !!
      !! Computes G * (M + m) for each massive body
      use swiftest
      implicit none

      self%mu_vec(:) = cb%Gmass + self%Gmass(:)

      return
   end procedure util_set_vec_mu_pl


   module procedure util_set_vec_mu_tp
      !! author: David A. Minton
      !!
      !! Converts certain scalar values to arrays so that they can be used in elemental functions
      use swiftest
      implicit none

      self%mu_vec(:) = self%mu
      self%dt_vec(:) = dt

      return
   end procedure util_set_vec_mu_tp

end submodule s_util_set_vec

