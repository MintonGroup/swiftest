submodule (swiftest_data_structures) s_swiftest_allocate_body
contains
   module procedure swiftest_allocate_body
   !! author: David A. Minton
   !!
   !! Constructor for base Swiftest particle class. Allocates space for all particles and
   !! initializes all components with a value. Also sets the is_allocated flag to true.
   use swiftest
   implicit none

   self%nbody = n
   if (n <= 0) return

   if (self%is_allocated) then
      !write(*,*) 'Swiftest particle structure already alllocated'
      return
   end if
   !write(*,*) 'Allocating the basic Swiftest particle'
   allocate(self%name(n))
   allocate(self%status(n))
   allocate(self%mu_vec(n))
   allocate(self%dt_vec(n))
   allocate(self%lspill_list(n))

   self%name(:) = 0
   self%status(:) = 0
   self%is_allocated = .true.
   self%lspill = .false.
   self%ldiscard = .false.


   return
   end procedure swiftest_allocate_body
end submodule s_swiftest_allocate_body
