submodule (swiftest_data_structures) s_helio_spill_tp
contains
   module procedure helio_spill_tp
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
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

   ! Call the spill method for the parent class (the base class in this case)
   call self%swiftest_pl%spill(discard)

   real(DP), dimension(:,:), allocatable :: ah  !! Total heliocentric acceleration
   real(DP), dimension(:,:), allocatable :: ahi !! Heliocentric acceleration due to interactions

   ! Pack the discarded bodies into the discard object
   discard%ah(:)  = pack(self%ah(:),  self%ldiscard)
   discard%ahi(:) = pack(self%ahi(:), self%ldiscard)

   ! Pack the kept bodies back into the original object
   self%ah(:) = pack(self%ah(:),  .not. self%ldiscard)
   self%ahi(:)= pack(self%ahi(:), .not. self%ldiscard)

   end procedure helio_spill_tp
end submodule s_helio_spill_tp



    


