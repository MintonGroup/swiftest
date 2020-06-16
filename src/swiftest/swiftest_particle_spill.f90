submodule (swiftest_data_structures) s_swiftest_particle_spill
contains
   module procedure swiftest_particle_spill
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest particle structure from active list to discard list
   use swiftest    
   integer(I4B) :: nspill

   if (.not. self%lspill) then  ! Only calculate the number of spilled particles if this method is called directly. It won't recompute if this method
                                ! is called from higher up in the class heirarchy 
      self%lspill_list = (self%status(:) /= ACTIVE) ! Use automatic allocation to allocate this logical flag
      nspill = count(self%lspill_list)
      self%nbody = self%nbody - nspill
      call discard%alloc(nspill) ! Create the discard object
   end if

   ! Pack the discarded bodies into the discard object
   discard%name(:)   = pack(self%name(:),   self%lspill_list)
   discard%status(:) = pack(self%status(:), self%lspill_list)
   discard%mu_vec(:) = pack(self%mu_vec(:), self%lspill_list)
   discard%dt_vec(:) = pack(self%dt_vec(:), self%lspill_list)

   ! Pack the kept bodies back into the original object
   self%name(:)   = pack(self%name(:),   .not. self%lspill_list)
   self%status(:) = pack(self%status(:), .not. self%lspill_list)
   self%mu_vec(:) = pack(self%mu_vec(:), .not. self%lspill_list)
   self%dt_vec(:) = pack(self%dt_vec(:), .not. self%lspill_list)

   ! This is the base class, so will be the last to be called. Therefore we must reset spill flag for the next discard operation.
   self%lspill = .false.
   
   end procedure swiftest_particle_spill
end submodule s_swiftest_particle_spill



    


