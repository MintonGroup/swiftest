submodule (swiftest_classes) s_io_initialize_system
contains
   module procedure io_initialize_system
   !! author: David A. Minton
   !!
   !! Wrapper method to initialize a basic Swiftest nbody system from files
   !!
   implicit none

   call self%cb%initialize(config)
   call self%pl%initialize(config)
   call self%tp%initialize(config)

end procedure io_initialize_system

end submodule s_io_initialize_system

