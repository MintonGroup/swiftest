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
      case (BS)
         write(*,*) 'Bulirsch-Stoer integrator not yet enabled'
      case (HELIO)
         write(*,*) 'Democratic Heliocentric integrator not yet enabled'
      case (RA15)
         write(*,*) 'Radau integrator not yet enabled'
      case (TU4)
         write(*,*) 'TU4 integrator not yet enabled'
      case (WHM)
         allocate(whm_system :: self)
      case (RMVS)
         write(*,*) 'RMVS integrator not yet enabled'
      case (SYMBA)
         write(*,*) 'SyMBA integrator not yet enabled'
      case (RINGMOONS)
         write(*,*) 'RINGMOONS-SyMBA integrator not yet enabled'
      case default
         write(*,*) 'Unkown integrator',integrator
         call util_exit(FAILURE)
      end select
      call self%construct(config, integrator)

      end procedure setup_construct_system

   end submodule s_setup_construct_system

   integer(I4B), parameter :: BS                 = 2
   integer(I4B), parameter :: HELIO              = 3
   integer(I4B), parameter :: RA15               = 4
   integer(I4B), parameter :: TU4                = 5
   integer(I4B), parameter :: WHM                = 6
   integer(I4B), parameter :: RMVS               = 7
   integer(I4B), parameter :: SYMBA              = 8
   integer(I4B), parameter :: RINGMOONS          = 9
