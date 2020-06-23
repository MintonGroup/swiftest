submodule(swiftest_class) s_nbody_cb

contains
   module procedure nbody_initialize_cb
      use swiftest
      implicit none
      return
   end procedure nbody_initialize_cb

   module procedure nbody_step_cb
      use swiftest
      implicit none
      return
   end procedure nbody_step_cb
end submodule s_nbody_cb
