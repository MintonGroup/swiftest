submodule (swiftest_data_structures) s_swiftest_spill_tp
contains
   module procedure swiftest_spill_pl
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest massive body particle structure from active list to discard list
   use swiftest
   integer(I4B) :: nspill

   npl = self%nbody
   if (.not. self%lspill) then  ! Only calculate the number of spilled particles if this method is called directly. It won't recompute if this method
                                ! is called from higher up in the class heirarchy 
      self%lspill_list(1:npl) = (self%status(1:npl) /= ACTIVE) ! Use automatic allocation to allocate this logical flag
      nspill = count(self%lspill_list(1:npl))
      call discard%alloc(nspill) ! Create the discard object
      self%lspill = .true.
    end if


   ! Pack the discarded bodies into the discard object
   discard%mass(:)   = pack(self%mass(1:npl),   self%lspill_list(1:npl))
   discard%radius(:) = pack(self%radius(1:npl), self%lspill_list(1:npl))
   discard%rhill(:)  = pack(self%rhill(1:npl),  self%lspill_list(1:npl))

   ! Pack the kept bodies back into the original object
   self%mass(:)   = pack(self%mass(1:npl),   .not. self%lspill_list(1:npl))
   self%radius(:) = pack(self%radius(1:npl), .not. self%lspill_list(1:npl))
   self%rhill(:)  = pack(self%rhill(1:npl),  .not. self%lspill_list(1:npl))

   ! Call the spill method for the parent class 
   call self%swiftest_tp%spill(discard)

   return
   end procedure swiftest_spill_pl
end submodule s_swiftest_spill_pl



    


