submodule (symba) s_symba_spill_tp
contains
   module procedure symba_spill_tp
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) symba test particle structure from active list to discard list
   use swiftest
   integer(I4B) :: nspill, ntp

   ntp = self%nbody
   nspill = self%nspill
   if (.not. self%lspill) then
      call discard%alloc(nspill) ! Create the discard object for this type
      self%lspill = .true.
   end if

   ! Pack the discarded bodies into the discard object
   discard%nplenc(:)  = pack(self%nplenc(1:ntp),  self%lspill_list(1:ntp))
   discard%levelg(:) = pack(self%levelg(1:ntp), self%lspill_list(1:ntp))
   discard%levelm(:) = pack(self%levelm(1:ntp), self%lspill_list(1:ntp))

   ! Pack the kept bodies back into the original object
   self%nplenc(:) = pack(self%nplenc(1:ntp),  .not. self%lspill_list(1:ntp))
   self%levelg(:)= pack(self%levelg(1:ntp), .not. self%lspill_list(1:ntp))
   self%levelm(:)= pack(self%levelm(1:ntp), .not. self%lspill_list(1:ntp))

   ! Call the spill method for the parent class 
   call helio_spill_pl(self,discard)

   return

   end procedure symba_spill_tp
end submodule s_symba_spill_tp