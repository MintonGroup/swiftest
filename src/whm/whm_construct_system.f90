submodule (whm_classes) s_whm_construct_system

contains

   module procedure whm_construct_system
      !! author: David A. Minton
      !!
      !! Constructor for a WHM nbody system. Creates the nbody system object based on the user-input integrator
      !! 
      use swiftest
      implicit none

      ! self, config, integrator
      allocate(whm_pl :: self%pl)
      allocate(whm_tp :: self%tp)
      allocate(whm_tp :: self%tp_discards)
  

      end procedure whm_construct_system

   end submodule s_whm_construct_system
