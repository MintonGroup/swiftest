submodule (swiftest_data_structures) s_swiftest_tp_spill
contains
   module procedure swiftest_tp_spill
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
   use swiftest
   integer(I4B) :: nspill, ntp

   ntp = self%nbody
   if (.not. self%lspill) then  ! Only calculate the number of spilled particles if this method is called directly. It won't recompute if this method
                                ! is called from higher up in the class heirarchy 
      self%lspill_list(1:ntp) = (self%status(:) /= ACTIVE) ! Use automatic allocation to allocate this logical flag
      nspill = count(self%lspill_list)
      call discard%alloc(nspill) ! Create the discard object
      self%lspill = .true.
    end if


   ! Pack the discarded bodies into the discard object
   discard%peri(:)   = pack(self%peri(1:ntp),   self%lspill_list(1:ntp))
   discard%atp(:)    = pack(self%atp(1:ntp),    self%lspill_list(1:ntp))
   discard%isperi(:) = pack(self%isperi(1:ntp), self%lspill_list(1:ntp))
   discard%xh(:)     = pack(self%xh(1:ntp),     self%lspill_list(1:ntp))
   discard%vh(:)     = pack(self%vh(1:ntp),     self%lspill_list(1:ntp))
   discard%xb(:)     = pack(self%xb(1:ntp),     self%lspill_list(1:ntp))
   discard%vb(:)     = pack(self%vb(1:ntp),     self%lspill_list(1:ntp))

   ! Pack the kept bodies back into the original object
   self%peri(:)   = pack(self%peri(1:ntp),   .not. self%lspill_list(1:ntp))
   self%atp(:)    = pack(self%atp(1:ntp),    .not. self%lspill_list(1:ntp))
   self%isperi(:) = pack(self%isperi(1:ntp), .not. self%lspill_list(1:ntp))
   self%xh(:)     = pack(self%xh(1:ntp),     .not. self%lspill_list(1:ntp))
   self%vh(:)     = pack(self%vh(1:ntp),     .not. self%lspill_list(1:ntp))
   self%xb(:)     = pack(self%xb(1:ntp),     .not. self%lspill_list(1:ntp))
   self%vb(:)     = pack(self%vb(1:ntp),     .not. self%lspill_list(1:ntp)))

   ! Call the spill method for the parent class (the base class in this case)
   call self%swiftest_particle%spill(discard)

   end procedure swiftest_tp_spill
end submodule s_swiftest_tp_spill



    


