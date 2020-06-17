submodule (swiftest_data_structures) s_swiftest_spill_body
contains
   module procedure swiftest_spill_body
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest particle structure from active list to discard list
   use swiftest    
   implicit none

   associate(np => self%nbody)
      self%lspill_list(1:np) = (self%status(1:np) /= ACTIVE) ! Use automatic allocation to allocate this logical flag
      self%nspill = count(self%lspill_list(1:np))

      select type(self)
      class is (swiftest_tp)
         call swiftest_spill_tp(self, discard)
      class is (swiftest_pl)
         call swiftest_spill_pl(self, discard)
      class is (helio_tp)
         call helio_spill_tp(self, discard)
      class is (helio_pl)
         call helio_spill_pl(self, discard)
      class is (symba_tp)
         call symba_spill_tp(self, discard)
      class is (symba_pl)
         call symba_spill_pl(self, discard)
      class default
         write(*,*) 'Unrecognized type of body'
      end select

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
      self%nbody = np - self%nspill
   end associate
   
   end procedure swiftest_spill_body
end submodule s_swiftest_spill_body



    


