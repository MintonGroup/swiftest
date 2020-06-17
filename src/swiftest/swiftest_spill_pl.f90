submodule (swiftest_data_structures) s_swiftest_spill_pl
contains
   module procedure swiftest_spill_pl
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest massive body particle structure from active list to discard list
   use swiftest
   implicit none

   associate(npl => self%nbody, nspill => self%nspill)
      if (.not. self%lspill) then
         call discard%alloc(nspill) ! Create the discard object for this type
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
      call swiftest_spill_tp(self,discard)
   end associate

   return
   end procedure swiftest_spill_pl
end submodule s_swiftest_spill_pl



    


