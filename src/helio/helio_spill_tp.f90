submodule (helio) s_helio_spill_tp
contains
   module procedure helio_spill_tp
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
   use swiftest
   implicit none

   integer(I4B) :: nspill, ntp, i

   ntp = self%nbody
   nspill = self%nspill
   if (.not. self%lspill) then
      call discard%alloc(nspill) ! Create the discard object for this type
      self%lspill = .true.
   end if

   do concurrent (i = 1:NDIM)
      ! Pack the discarded bodies into the discard object
      discard%ah(i,:)  = pack(self%ah(i,1:ntp),  self%lspill_list(1:ntp))
      discard%ahi(i,:) = pack(self%ahi(i,1:ntp), self%lspill_list(1:ntp))

      ! Pack the kept bodies back into the original object
      self%ah(i,:)  = pack(self%ah(i,1:ntp),  .not. self%lspill_list(1:ntp))
      self%ahi(i,:) = pack(self%ahi(i,1:ntp), .not. self%lspill_list(1:ntp))
   end do

   ! Call the spill method for the parent class (the base class in this case)
   call nbody_spill_pl(self,discard)

   return

   end procedure helio_spill_tp
end submodule s_helio_spill_tp



    


