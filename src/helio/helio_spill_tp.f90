submodule (helio) s_helio_spill_tp
contains
   module procedure helio_spill_tp
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
   call self%swiftest_pl%spill(discard)

   ! Pack the discarded bodies into the discard object
   discard%ah(:)  = pack(self%ah(:),  self%lspill_list)
   discard%ahi(:) = pack(self%ahi(:), self%lspill_list)

   ! Pack the kept bodies back into the original object
   self%ah(:) = pack(self%ah(:),  .not. self%lspill_list)
   self%ahi(:)= pack(self%ahi(:), .not. self%lspill_list)

   end procedure helio_spill_tp
end submodule s_helio_spill_tp



    


