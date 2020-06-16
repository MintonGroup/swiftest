submodule (swiftest_data_structures) s_swiftest_tp_spill
contains
   module procedure swiftest_tp_spill
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest massive body particle structure from active list to discard list
   use swiftest
   integer(I4B) :: nspill

   if (.not. self%lspill) then  ! Only calculate the number of spilled particles if this method is called directly. It won't recompute if this method
                                ! is called from higher up in the class heirarchy 
      self%ldiscard = (self%status(:) /= ACTIVE) ! Use automatic allocation to allocate this logical flag
      nspill = count(self%ldiscard)
      self%nbody = self%nbody - nspill
      call discard%alloc(nspill) ! Create the discard object
      self%lspill = .true.
    end if

   ! Call the spill method for the parent class 
   call self%swiftest_tp%spill(discard)

   ! Pack the discarded bodies into the discard object
   discard%mass(:)   = pack(self%mass(:),   self%ldiscard)
   discard%radius(:) = pack(self%radius(:), self%ldiscard)
   discard%rhill(:)  = pack(self%rhill(:),  self%ldiscard)

   ! Pack the kept bodies back into the original object
   self%mass(:)   = pack(self%mass(:),   .not. self%ldiscard)
   self%radius(:) = pack(self%radius(:), .not. self%ldiscard)
   self%rhill(:)  = pack(self%rhill(:),  .not. self%ldiscard)

   return
   end procedure swiftest_tp_spill
end submodule s_swiftest_tp_spill



    


