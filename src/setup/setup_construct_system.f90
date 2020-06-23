submodule (swiftest_classes) s_setup_construct_system

contains

   module procedure setup_construct_system
      !! author: David A. Minton
      !!
      !! Constructor for a Swiftest nbody system. Creates the nbody system object based on the user-input integrator
      !! 
      use swiftest
      implicit none

      ! self, config, integrator

      select case(integrator)
      case (WHM)
         allocate(whm_system :: self)
      case default
         write(*,*) 'Integrator not yet enabled'
         call util_exit(FAILURE)
      end select
      call self%construct(config, integrator)

      end procedure setup_construct_system

   end submodule s_setup_construct_system
