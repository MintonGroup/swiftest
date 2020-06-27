submodule (symba) s_symba_spill_pl
contains
   module procedure symba_spill_pl
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) symba massive body structure from active list to discard list
   use swiftest
   integer(I4B) :: i,nspill, npl

   npl = self%nbody
   nspill = self%nspill
   if (.not. self%lspill) then
      call discard%alloc(nspill) ! Create the discard object for this type
      self%lspill = .true.
   end if

   ! Pack the discarded bodies into the discard object
   discard%lmerged(:) = pack(self%lmerged(1:npl), self%lspill_list(1:npl))
   discard%nplenc(:)  = pack(self%nplenc(1:npl),  self%lspill_list(1:npl))
   discard%nchild(:) = pack(self%nchild(1:npl), self%lspill_list(1:npl))
   discard%index_parent(:) = pack(self%index_parent(1:npl), self%lspill_list(1:npl))

   ! Pack the kept bodies back into the original object
   self%lmerged(:)= pack(self%lmerged(1:npl), .not. self%lspill_list(1:npl))
   self%nplenc(:) = pack(self%nplenc(1:npl),  .not. self%lspill_list(1:npl))
   self%nchild(:)= pack(self%nchild(1:npl), .not. self%lspill_list(1:npl))
   self%index_parent(:)= pack(self%index_parent(1:npl), .not. self%lspill_list(1:npl))

   do concurrent (i = 1:NDIM)
      discard%index_child(:, i) = pack(self%index_child(1:npl, i), self%lspill_list(1:npl))
      self%index_child(:, i)= pack(self%index_child(1:npl, i), .not. self%lspill_list(1:npl))
   end do

   ! Call the spill method for the parent class 
   call symba_spill_tp(self,discard)

   return

   end procedure symba_spill_pl
end submodule s_symba_spill_pl



    


