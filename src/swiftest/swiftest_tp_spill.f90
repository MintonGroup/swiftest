submodule (swiftest_data_structures) s_swiftest_tp_spill
contains
   module procedure swiftest_tp_spill
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
   use swiftest
   integer(I4B) :: nspill

   if (.not. self%lspill) then  ! Only calculate the number of spilled particles if this method is called directly. It won't recompute if this method
                                ! is called from higher up in the class heirarchy 
      self%lspill_list = (self%status(:) /= ACTIVE) ! Use automatic allocation to allocate this logical flag
      nspill = count(self%lspill_list)
      self%nbody = self%nbody - nspill
      call discard%alloc(nspill) ! Create the discard object
      self%lspill = .true.
    end if

   ! Call the spill method for the parent class (the base class in this case)
   call self%swiftest_particle%spill(discard)

   ! Pack the discarded bodies into the discard object
   discard%peri(:)   = pack(self%peri(:),   self%lspill_list)
   discard%atp(:)    = pack(self%atp(:),    self%lspill_list)
   discard%isperi(:) = pack(self%isperi(:), self%lspill_list)
   discard%xh(:)     = pack(self%xh(:),     self%lspill_list)
   discard%vh(:)     = pack(self%vh(:),     self%lspill_list)
   discard%xb(:)     = pack(self%xb(:),     self%lspill_list)
   discard%vb(:)     = pack(self%vb(:),     self%lspill_list)

   ! Pack the kept bodies back into the original object
   self%peri(:)   = pack(self%peri(:),   .not. self%lspill_list)
   self%atp(:)    = pack(self%atp(:),    .not. self%lspill_list)
   self%isperi(:) = pack(self%isperi(:), .not. self%lspill_list)
   self%xh(:)     = pack(self%xh(:),     .not. self%lspill_list)
   self%vh(:)     = pack(self%vh(:),     .not. self%lspill_list)
   self%xb(:)     = pack(self%xb(:),     .not. self%lspill_list)
   self%vb(:)     = pack(self%vb(:),     .not. self%lspill_list))

   end procedure swiftest_tp_spill
end submodule s_swiftest_tp_spill



    


