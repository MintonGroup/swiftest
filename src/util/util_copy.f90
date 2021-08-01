submodule (swiftest_classes) s_util_copy
   use swiftest
contains

   module subroutine util_copy_into_body(self, source, param, lsource_mask)
      !! author: David A. Minton
      !!
      !! Copies elements from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),            intent(inout) :: self   !! Swiftest body object
      class(swiftest_body),            intent(in)    :: source !! Source object to append
      class(swiftest_parameters),      intent(in)    :: param  !! Current run configuration parameters
      logical, dimension(:), optional, intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B)  :: i,nnew
      logical, dimension(:), allocatable :: lfill_list

      if (present(lsource_mask)) then
         nnew = count(lsource_mask)
      else
         nnew = size(source%status)
      end if
      allocate(lfill_list(size(self%status)))
      lfill_list = .false.
      lfill_list(1:nnew) = .true.
      associate(nold => self%nbody)
         if (nnew > size(self%status)) call self%resize(nnew, param)
         call self%fill(source, lfill_list)
      end associate
      return
   end subroutine util_copy_into_body


   module subroutine util_copy_spill_body(self, discards, lspill_list)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest generic particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self       !! Swiftest generic body object
      class(swiftest_body), intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)   :: lspill_list !! Logical array of bodies to spill into the discards
      ! Internals
      integer(I4B) :: i

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         discards%id(:)       = pack(keeps%id(:),     lspill_list(:))
         discards%name(:)     = pack(keeps%name(:),   lspill_list(:))
         discards%status(:)   = pack(keeps%status(:), lspill_list(:))
         discards%mu(:)       = pack(keeps%mu(:),     lspill_list(:))
         discards%lmask(:)    = pack(keeps%lmask(:),  lspill_list(:))
         do i = 1, NDIM
            discards%xh(i, :) = pack(keeps%xh(i, :),     lspill_list(:))
            discards%vh(i, :) = pack(keeps%vh(i, :),     lspill_list(:))
            discards%xb(i, :) = pack(keeps%xb(i, :),     lspill_list(:))
            discards%vb(i, :) = pack(keeps%vb(i, :),     lspill_list(:))
            discards%ah(i, :) = pack(keeps%ah(i, :),     lspill_list(:))
         end do

         if (allocated(keeps%a))     discards%a(:)        = pack(keeps%a(:),      lspill_list(:))
         if (allocated(keeps%e))     discards%e(:)        = pack(keeps%e(:),      lspill_list(:))
         if (allocated(keeps%capom)) discards%capom(:)    = pack(keeps%capom(:),  lspill_list(:))
         if (allocated(keeps%omega)) discards%omega(:)    = pack(keeps%omega(:),  lspill_list(:))
         if (allocated(keeps%capm))  discards%capm(:)     = pack(keeps%capm(:),   lspill_list(:))

         if (allocated(keeps%aobl)) then
            do i = 1, NDIM
               discards%aobl(i, :) = pack(keeps%aobl(i, :), lspill_list(:))
            end do
         end if
         if (allocated(keeps%agr)) then
            do i = 1, NDIM
               discards%agr(i, :)   = pack(keeps%agr(i, :), lspill_list(:))
            end do
         end if
         if (allocated(keeps%atide)) then
            do i = 1, NDIM
               discards%atide(i, :) = pack(keeps%atide(i, :), lspill_list(:))
            end do
         end if

         if (count(.not.lspill_list(:)) > 0) then 
            keeps%id(:)         = pack(keeps%id(:),     .not. lspill_list(:))
            keeps%name(:)       = pack(keeps%name(:),   .not. lspill_list(:))
            keeps%status(:)     = pack(keeps%status(:), .not. lspill_list(:))
            keeps%mu(:)         = pack(keeps%mu(:),     .not. lspill_list(:))
            keeps%lmask(:)      = pack(keeps%lmask(:),  .not. lspill_list(:))

            do i = 1, NDIM
               keeps%xh(i, :)    = pack(keeps%xh(i, :),   .not. lspill_list(:))
               keeps%vh(i, :)    = pack(keeps%vh(i, :),   .not. lspill_list(:))
               keeps%xb(i, :)    = pack(keeps%xb(i, :),   .not. lspill_list(:))
               keeps%vb(i, :)    = pack(keeps%vb(i, :),   .not. lspill_list(:))
               keeps%ah(i, :)    = pack(keeps%ah(i, :),   .not. lspill_list(:))
            end do

            if (allocated(keeps%a))     keeps%a(:)          = pack(keeps%a(:),      .not. lspill_list(:))
            if (allocated(keeps%e))     keeps%e(:)          = pack(keeps%e(:),      .not. lspill_list(:))
            if (allocated(keeps%inc))   keeps%inc(:)        = pack(keeps%inc(:),    .not. lspill_list(:))
            if (allocated(keeps%capom)) keeps%capom(:)      = pack(keeps%capom(:),  .not. lspill_list(:))
            if (allocated(keeps%omega)) keeps%omega(:)      = pack(keeps%omega(:),  .not. lspill_list(:))
            if (allocated(keeps%capm))  keeps%capm(:)       = pack(keeps%capm(:),   .not. lspill_list(:))

            if (allocated(keeps%aobl)) then
               do i = 1, NDIM
                  keeps%aobl(i, :)  = pack(keeps%aobl(i, :), .not. lspill_list(:))
               end do
            end if

            if (allocated(keeps%agr)) then
               do i = 1, NDIM
                  keeps%agr(i, :)   = pack(keeps%agr(i, :),  .not. lspill_list(:))
               end do
            end if

            if (allocated(keeps%atide)) then
               do i = 1, NDIM
                  keeps%atide(i, :)  = pack(keeps%atide(i, :), .not. lspill_list(:))
               end do
            end if

         end if
         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nbody values for both the keeps and discareds
         discards%nbody = count(lspill_list(:))
         keeps%nbody = count(.not.lspill_list(:)) 
         if (allocated(keeps%ldiscard)) deallocate(keeps%ldiscard)
         if (allocated(discards%ldiscard)) deallocate(discards%ldiscard)
         allocate(keeps%ldiscard(keeps%nbody))
         allocate(discards%ldiscard(discards%nbody))
         keeps%ldiscard = .false.
         discards%ldiscard = .true.
      end associate
     
      return
   end subroutine util_copy_spill_body


   module subroutine util_copy_spill_pl(self, discards, lspill_list)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest massive body structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_pl),    intent(inout) :: self        !! Swiftest massive body object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      ! Internals
      integer(I4B) :: i

      associate(keeps => self)

         select type (discards) ! The standard requires us to select the type of both arguments in order to access all the components
         class is (swiftest_pl)
         !> Spill components specific to the massive body class
            discards%mass(:)     = pack(keeps%mass(:),    lspill_list(:))
            discards%Gmass(:)    = pack(keeps%Gmass(:),   lspill_list(:))
            discards%rhill(:)    = pack(keeps%rhill(:),   lspill_list(:))

            if (allocated(keeps%radius))  discards%radius(:)  = pack(keeps%radius(:),  lspill_list(:))
            if (allocated(keeps%density)) discards%density(:) = pack(keeps%density(:), lspill_list(:))
            if (allocated(keeps%k2))      discards%k2(:)      = pack(keeps%k2(:),      lspill_list(:))
            if (allocated(keeps%Q))       discards%Q(:)       = pack(keeps%Q(:),       lspill_list(:))
            if (allocated(keeps%tlag))    discards%tlag(:)    = pack(keeps%tlag(:),    lspill_list(:))

            if (allocated(keeps%xbeg)) then
               do i = 1, NDIM
                  discards%xbeg(i, :)  = pack(keeps%xbeg(i, :),     lspill_list(:))
               end do
            end if

            if (allocated(keeps%xend)) then
               do i = 1, NDIM
                  discards%xend(i, :)  = pack(keeps%xend(i, :),     lspill_list(:))
               end do
            end if

            if (allocated(keeps%vbeg)) then
               do i = 1, NDIM
                  discards%vbeg(i, :)  = pack(keeps%vbeg(i, :),     lspill_list(:))
               end do
            end if

            if (allocated(keeps%Ip)) then
               do i = 1, NDIM
                  discards%Ip(i, :)  = pack(keeps%Ip(i, :),     lspill_list(:))
               end do
            end if

            if (allocated(keeps%rot)) then
               do i = 1, NDIM
                  discards%rot(i, :)  = pack(keeps%rot(i, :),     lspill_list(:))
               end do
            end if

            if (count(.not.lspill_list(:))  > 0) then 
               keeps%mass(:)        = pack(keeps%mass(:),    .not. lspill_list(:))
               keeps%Gmass(:)       = pack(keeps%Gmass(:),   .not. lspill_list(:))
               keeps%rhill(:)       = pack(keeps%rhill(:),   .not. lspill_list(:))
               if (allocated(keeps%radius))  keeps%radius(:)      = pack(keeps%radius(:),  .not. lspill_list(:))
               if (allocated(keeps%density)) keeps%density(:)     = pack(keeps%density(:), .not. lspill_list(:))
               if (allocated(keeps%k2))      keeps%k2(:)          = pack(keeps%k2(:),      .not. lspill_list(:))
               if (allocated(keeps%Q))       keeps%Q(:)           = pack(keeps%Q(:),       .not. lspill_list(:))
               if (allocated(keeps%tlag))    keeps%tlag(:)        = pack(keeps%tlag(:),    .not. lspill_list(:))

               if (allocated(keeps%xbeg)) then
                  do i = 1, NDIM
                     keeps%xbeg(i,:)     = pack(keeps%xbeg(i,:),   .not. lspill_list(:))
                  end do
               end if

               if (allocated(keeps%xend)) then
                  do i = 1, NDIM
                     keeps%xend(i,:)     = pack(keeps%xend(i,:),   .not. lspill_list(:))
                  end do
               end if

               if (allocated(keeps%vbeg)) then
                  do i = 1, NDIM
                     keeps%vbeg(i,:)     = pack(keeps%vbeg(i,:),   .not. lspill_list(:))
                  end do
               end if

               if (allocated(keeps%Ip)) then
                  do i = 1, NDIM
                     keeps%Ip(i,:)       = pack(keeps%Ip(i,:),   .not. lspill_list(:))
                  end do
               end if

               if (allocated(keeps%rot)) then
                  do i = 1, NDIM
                     keeps%rot(i,:)       = pack(keeps%rot(i,:),   .not. lspill_list(:))
                  end do
               end if

            end if

            call util_copy_spill_body(keeps, discards, lspill_list)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on swiftest_pl'
         end select
      end associate

      return
   end subroutine util_copy_spill_pl   


   module subroutine util_copy_spill_tp(self, discards, lspill_list)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_tp),    intent(inout) :: self        !! Swiftest test particle object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discardse

      associate(keeps => self, ntp => self%nbody)
         select type(discards)
         class is (swiftest_tp)
            !> Spill components specific to the test particle class
            discards%isperi(:) = pack(keeps%isperi(:),       lspill_list(:))
            discards%peri(:)   = pack(keeps%peri(:),         lspill_list(:))
            discards%atp(:)    = pack(keeps%atp(:),          lspill_list(:))
            if (count(.not.lspill_list(:))  > 0) then 
               keeps%atp(:)       = pack(keeps%atp(:),    .not. lspill_list(:))
               keeps%peri(:)      = pack(keeps%peri(:),   .not. lspill_list(:))
               keeps%isperi(:)    = pack(keeps%isperi(:), .not. lspill_list(:))
            end if
            call util_copy_spill_body(keeps, discards, lspill_list)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on swiftest_tp'
         end select
      end associate

      return
   end subroutine util_copy_spill_tp

end submodule s_util_copy