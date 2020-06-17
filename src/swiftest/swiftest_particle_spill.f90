submodule (swiftest_data_structures) s_swiftest_particle_spill
contains
   module procedure swiftest_particle_spill
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest particle structure from active list to discard list
   use swiftest    
   integer(I4B) :: nspill, np

   np = self%nbody
   if (.not. self%lspill) then  ! Only calculate the number of spilled particles if this method is called directly. It won't recompute if this method
                                ! is called from higher up in the class heirarchy 
      self%lspill_list(1:np) = (self%status(1:np) /= ACTIVE) ! Use automatic allocation to allocate this logical flag
      nspill = count(self%lspill_list(1:np))
      call discard%alloc(nspill) ! Create the discard object
   end if

   ! Pack the discarded bodies into the discard object
   discard%name(:)   = pack(self%name(1:np),   self%lspill_list(1:np))
   discard%status(:) = pack(self%status(1:np), self%lspill_list(1:np))
   discard%mu_vec(:) = pack(self%mu_vec(1:np), self%lspill_list(1:np))
   discard%dt_vec(:) = pack(self%dt_vec(1:np), self%lspill_list(1:np))

   ! Pack the kept bodies back into the original object
   self%name(:)   = pack(self%name(1:np),   .not. self%lspill_list(1:np))
   self%status(:) = pack(self%status(1:np), .not. self%lspill_list(1:np))
   self%mu_vec(:) = pack(self%mu_vec(1:np), .not. self%lspill_list(1:np))
   self%dt_vec(:) = pack(self%dt_vec(1:np), .not. self%lspill_list(1:np))

   ! This is the base class, so will be the last to be called. Therefore we must reset spill flag for the next discard operation.
   self%lspill = .false.
   self%nbody = self%nbody - discard%nbody
   
   end procedure swiftest_particle_spill
end submodule s_swiftest_particle_spill



    


