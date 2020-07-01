submodule (whm_classes) s_whm_construct_system

contains

   module procedure whm_construct_system
      !! author: David A. Minton
      !!
      !! Constructor for a WHM nbody system. Creates the nbody system object with all of the WHM typed component object.
      !! 
      use swiftest
      implicit none
      
      allocate(whm_central_body :: self%cb)
      allocate(whm_pl :: self%pl)
      allocate(whm_tp :: self%tp)
      allocate(whm_tp :: self%tp_discards)
  
      return 
   end procedure whm_construct_system

end submodule s_whm_construct_system
