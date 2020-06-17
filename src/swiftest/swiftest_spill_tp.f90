submodule (swiftest_data_structures) s_swiftest_spill_tp
contains
   module procedure swiftest_spill_tp
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
   use swiftest
   integer(I4B) :: nspill, ntp

   ntp = self%nbody
   nspill = self%nspill
   if (.not. self%lspill) then
      call discard%alloc(nspill) ! Create the discard object for this type
      self%lspill = .true.
   end if

   ! Pack the discarded bodies into the discard object
   discard%peri(:)   = pack(self%peri(1:ntp),   self%lspill_list(1:ntp))
   discard%atp(:)    = pack(self%atp(1:ntp),    self%lspill_list(1:ntp))
   discard%isperi(:) = pack(self%isperi(1:ntp), self%lspill_list(1:ntp))

   ! Pack the kept bodies back into the original object
   self%peri(:)   = pack(self%peri(1:ntp),   .not. self%lspill_list(1:ntp))
   self%atp(:)    = pack(self%atp(1:ntp),    .not. self%lspill_list(1:ntp))
   self%isperi(:) = pack(self%isperi(1:ntp), .not. self%lspill_list(1:ntp))

   do concurrent (i = 1:NDIM)
      discard%xh(i,:) = pack(self%xh(i,1:ntp), self%lspill_list(1:ntp))
      discard%vh(i,:) = pack(self%vh(i,1:ntp), self%lspill_list(1:ntp))
      discard%xb(i,:) = pack(self%xb(i,1:ntp), self%lspill_list(1:ntp))
      discard%vb(i,:) = pack(self%vb(i,1:ntp), self%lspill_list(1:ntp))
      self%xh(i,:)    = pack(self%xh(i,1:ntp), .not. self%lspill_list(1:ntp))
      self%vh(i,:)    = pack(self%vh(i,1:ntp), .not. self%lspill_list(1:ntp))
      self%xb(i,:)    = pack(self%xb(i,1:ntp), .not. self%lspill_list(1:ntp))
      self%vb(i,:)    = pack(self%vb(i,1:ntp), .not. self%lspill_list(1:ntp)))
   end do

   end procedure swiftest_spill_tp
end submodule s_swiftest_spill_tp



    


