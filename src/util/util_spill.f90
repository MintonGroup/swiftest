submodule (swiftest_classes) s_util_spill
   use swiftest
contains

   module subroutine util_spill_arr_char_string(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of type character strings
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,               dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                                          intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (.not.allocated(discards)) then
         allocate(discards(nspill))
      else if (size(discards) /= nspill) then
         deallocate(discards)
         allocate(discards(nspill))
      end if

      discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
      if (ldestructive) then
         if (nkeep > 0) then
            keeps(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine util_spill_arr_char_string

   module subroutine util_spill_arr_DP(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of type DP
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      real(DP), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,  dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
      logical,                             intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (.not.allocated(discards)) then
         allocate(discards(nspill))
      else if (size(discards) /= nspill) then
         deallocate(discards)
         allocate(discards(nspill))
      end if

      discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
      if (ldestructive) then
         if (nkeep > 0) then
            keeps(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine util_spill_arr_DP

   module subroutine util_spill_arr_DPvec(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of DP vectors with shape (NDIM, n)
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      real(DP), dimension(:,:), allocatable, intent(inout) :: discards     !! Array discards
      logical,  dimension(:),                intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: i, nspill, nkeep, nlist

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (.not.allocated(discards)) then
         allocate(discards(NDIM, nspill))
      else if (size(discards, dim=2) /= nspill) then
         deallocate(discards)
         allocate(discards(NDIM, nspill))
      end if

      do i = 1, NDIM
         discards(i,:) = pack(keeps(i,1:nlist), lspill_list(1:nlist))
      end do
      if (ldestructive) then
         if (nkeep > 0) then
            do i = 1, NDIM
               keeps(i,:) = pack(keeps(i,1:nlist), .not. lspill_list(1:nlist))
            end do
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine util_spill_arr_DPvec

   module subroutine util_spill_arr_I4B(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of type I4B
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      integer(I4B), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,      dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                                 intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (.not.allocated(discards)) then
         allocate(discards(nspill))
      else if (size(discards) /= nspill) then
         deallocate(discards)
         allocate(discards(nspill))
      end if

      discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
      if (ldestructive) then
         if (nkeep > 0) then
            keeps(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine util_spill_arr_I4B

   module subroutine util_spill_arr_logical(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of logicals
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      logical, dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical, dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                            intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or no
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (.not.allocated(discards)) then
         allocate(discards(nspill))
      else if (size(discards) /= nspill) then
         deallocate(discards)
         allocate(discards(nspill))
      end if

      discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
      if (ldestructive) then
         if (nkeep > 0) then
            keeps(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine util_spill_arr_logical


   module subroutine util_spill_body(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest generic particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_body),  intent(inout) :: self         !! Swiftest generic body object
      class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      ! Internals

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         call util_spill(keeps%id, discards%id, lspill_list, ldestructive)
         call util_spill(keeps%name, discards%name, lspill_list, ldestructive)
         call util_spill(keeps%status, discards%status, lspill_list, ldestructive)
         call util_spill(keeps%lmask, discards%lmask, lspill_list, ldestructive)
         call util_spill(keeps%ldiscard, discards%ldiscard, lspill_list, ldestructive)
         call util_spill(keeps%mu, discards%mu, lspill_list, ldestructive)
         call util_spill(keeps%xh, discards%xh, lspill_list, ldestructive)
         call util_spill(keeps%vh, discards%vh, lspill_list, ldestructive)
         call util_spill(keeps%xb, discards%xb, lspill_list, ldestructive)
         call util_spill(keeps%vb, discards%vb, lspill_list, ldestructive)
         call util_spill(keeps%ah, discards%ah, lspill_list, ldestructive)
         call util_spill(keeps%aobl, discards%aobl, lspill_list, ldestructive)
         call util_spill(keeps%agr, discards%agr, lspill_list, ldestructive)
         call util_spill(keeps%atide, discards%atide, lspill_list, ldestructive)
         call util_spill(keeps%a, discards%a, lspill_list, ldestructive)
         call util_spill(keeps%e, discards%e, lspill_list, ldestructive)
         call util_spill(keeps%inc, discards%inc, lspill_list, ldestructive)
         call util_spill(keeps%capom, discards%capom, lspill_list, ldestructive)
         call util_spill(keeps%omega, discards%omega, lspill_list, ldestructive)
         call util_spill(keeps%capm, discards%capm, lspill_list, ldestructive)

         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nbody values for both the keeps and discareds
         discards%nbody = count(lspill_list(:))
         keeps%nbody = keeps%nbody - discards%nbody
         if (keeps%nbody > size(keeps%status)) keeps%status(keeps%nbody+1:size(keeps%status)) = INACTIVE
      end associate
     
      return
   end subroutine util_spill_body

   module subroutine util_spill_encounter(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest encounter structure from active list to discard list
      implicit none
      ! Arguments
      class(swiftest_encounter), intent(inout) :: self         !! Swiftest encounter list 
      class(swiftest_encounter), intent(inout) :: discards     !! Discarded object 
      logical, dimension(:),     intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                   intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
  
      associate(keeps => self)
         call util_spill(keeps%lvdotr, discards%lvdotr, lspill_list, ldestructive)
         call util_spill(keeps%status, discards%status, lspill_list, ldestructive)
         call util_spill(keeps%index1, discards%index1, lspill_list, ldestructive)
         call util_spill(keeps%index2, discards%index2, lspill_list, ldestructive)
         call util_spill(keeps%x1, discards%x1, lspill_list, ldestructive)
         call util_spill(keeps%x2, discards%x2, lspill_list, ldestructive)
         call util_spill(keeps%v1, discards%v1, lspill_list, ldestructive)
         call util_spill(keeps%v2, discards%v2, lspill_list, ldestructive)
         call util_spill(keeps%t, discards%t, lspill_list, ldestructive)

         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nenc values for both the keeps and discareds
         discards%nenc = count(lspill_list(:))
         keeps%nenc = count(.not.lspill_list(:)) 
         if (keeps%nenc > size(keeps%status)) keeps%status(keeps%nenc+1:size(keeps%status)) = INACTIVE
      end associate
   
      return
   end subroutine util_spill_encounter


   module subroutine util_spill_pl(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest massive body structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_pl),    intent(inout) :: self        !! Swiftest massive body object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list

      associate(keeps => self)
         select type (discards) ! The standard requires us to select the type of both arguments in order to access all the components
         class is (swiftest_pl)
            !> Spill components specific to the massive body class
            call util_spill(keeps%mass, discards%mass, lspill_list, ldestructive)
            call util_spill(keeps%Gmass, discards%Gmass, lspill_list, ldestructive)
            call util_spill(keeps%rhill, discards%rhill, lspill_list, ldestructive)
            call util_spill(keeps%radius, discards%radius, lspill_list, ldestructive)
            call util_spill(keeps%density, discards%density, lspill_list, ldestructive)
            call util_spill(keeps%k2, discards%k2, lspill_list, ldestructive)
            call util_spill(keeps%Q, discards%Q, lspill_list, ldestructive)
            call util_spill(keeps%tlag, discards%tlag, lspill_list, ldestructive)
            call util_spill(keeps%xbeg, discards%xbeg, lspill_list, ldestructive)
            call util_spill(keeps%vbeg, discards%vbeg, lspill_list, ldestructive)
            call util_spill(keeps%Ip, discards%Ip, lspill_list, ldestructive)
            call util_spill(keeps%rot, discards%rot, lspill_list, ldestructive)

            call util_spill_body(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on swiftest_pl'
         end select
      end associate

      return
   end subroutine util_spill_pl   


   module subroutine util_spill_tp(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_tp),    intent(inout) :: self        !! Swiftest test particle object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discardse
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list

      associate(keeps => self, ntp => self%nbody)
         select type(discards)
         class is (swiftest_tp)
            !> Spill components specific to the test particle class
            call util_spill(keeps%isperi, discards%isperi, lspill_list, ldestructive)
            call util_spill(keeps%peri, discards%peri, lspill_list, ldestructive)
            call util_spill(keeps%atp, discards%atp, lspill_list, ldestructive)

            call util_spill_body(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on swiftest_tp'
         end select
      end associate

      return
   end subroutine util_spill_tp

end submodule s_util_spill