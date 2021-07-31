submodule (swiftest_classes) s_util_copy
   use swiftest
contains

   module subroutine util_copy_fill_body(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest generic particle structure into an old one. 
      !! This is the inverse of a fill operation.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self       !! Swiftest generic body object
      class(swiftest_body), intent(in)    :: inserts     !! Inserted object 
      logical, dimension(:), intent(in)   :: lfill_list  !! Logical array of bodies to merge into the keeps
      ! internals
      integer(I4B) :: i

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         keeps%id(:)     = unpack(keeps%id(:), .not.lfill_list(:), keeps%id(:))
         keeps%id(:)     = unpack(inserts%id(:), lfill_list(:), keeps%id(:))

         keeps%name(:)     = unpack(keeps%name(:), .not.lfill_list(:), keeps%name(:))
         keeps%name(:)     = unpack(inserts%name(:), lfill_list(:), keeps%name(:))

         keeps%status(:)   = unpack(keeps%status(:), .not.lfill_list(:), keeps%status(:))
         keeps%status(:)   = unpack(inserts%status(:), lfill_list(:), keeps%status(:))

         keeps%ldiscard(:) = unpack(keeps%ldiscard(:), .not.lfill_list(:), keeps%ldiscard(:))
         keeps%ldiscard(:) = unpack(inserts%ldiscard(:), lfill_list(:), keeps%ldiscard(:))

         keeps%mu(:)    = unpack(keeps%mu(:),   .not.lfill_list(:), keeps%mu(:))
         keeps%mu(:)    = unpack(inserts%mu(:),   lfill_list(:), keeps%mu(:))

         keeps%lmask(:) = unpack(keeps%lmask(:), .not.lfill_list(:), keeps%ldiscard(:))
         keeps%lmask(:) = unpack(inserts%lmask(:), lfill_list(:), keeps%ldiscard(:))

         do i = 1, NDIM
            keeps%xh(i, :)    = unpack(keeps%xh(i, :), .not.lfill_list(:), keeps%xh(i, :))
            keeps%xh(i, :)    = unpack(inserts%xh(i, :), lfill_list(:), keeps%xh(i, :))

            keeps%vh(i, :)    = unpack(keeps%vh(i, :), .not.lfill_list(:), keeps%vh(i, :))
            keeps%vh(i, :)    = unpack(inserts%vh(i, :), lfill_list(:), keeps%vh(i, :))

            keeps%xb(i, :)    = unpack(keeps%xb(i, :), .not.lfill_list(:), keeps%xb(i, :))
            keeps%xb(i, :)    = unpack(inserts%xb(i, :), lfill_list(:), keeps%xb(i, :))

            keeps%vb(i, :)    = unpack(keeps%vb(i, :), .not.lfill_list(:), keeps%vb(i, :))
            keeps%vb(i, :)    = unpack(inserts%vb(i, :), lfill_list(:), keeps%vb(i, :))
            
            keeps%ah(i, :)    = unpack(keeps%ah(i, :), .not.lfill_list(:), keeps%ah(i, :))
            keeps%ah(i, :)    = unpack(inserts%ah(i, :), lfill_list(:), keeps%ah(i, :))
         end do

         if (allocated(keeps%aobl)) then
            do i = 1, NDIM
               keeps%aobl(i, :)  = unpack(keeps%aobl(i, :), .not.lfill_list(:), keeps%aobl(i, :))
               keeps%aobl(i, :)  = unpack(inserts%aobl(i, :), lfill_list(:), keeps%aobl(i, :))
            end do
         end if

         if (allocated(keeps%agr)) then
            do i = 1, NDIM
               keeps%agr(i, :)  = unpack(keeps%agr(i, :), .not.lfill_list(:), keeps%agr(i, :))
               keeps%agr(i, :)  = unpack(inserts%agr(i, :), lfill_list(:), keeps%agr(i, :))
            end do
         end if

         if (allocated(keeps%atide)) then
            do i = 1, NDIM
               keeps%atide(i, :)  = unpack(keeps%atide(i, :), .not.lfill_list(:), keeps%atide(i, :))
               keeps%atide(i, :)  = unpack(inserts%atide(i, :), lfill_list(:), keeps%atide(i, :))
            end do
         end if
        
         if (allocated(keeps%a)) then
            keeps%a(:)     = unpack(keeps%a(:),    .not.lfill_list(:), keeps%a(:))
            keeps%a(:)     = unpack(inserts%a(:),    lfill_list(:), keeps%a(:))
         end if
           
         if (allocated(keeps%e)) then
            keeps%e(:)     = unpack(keeps%e(:),    .not.lfill_list(:), keeps%e(:))
            keeps%e(:)     = unpack(inserts%e(:),    lfill_list(:), keeps%e(:))
         end if
           
         if (allocated(keeps%inc)) then
            keeps%inc(:)   = unpack(keeps%inc(:),  .not.lfill_list(:), keeps%inc(:))
            keeps%inc(:)   = unpack(inserts%inc(:),  lfill_list(:), keeps%inc(:))
         end if
           
         if (allocated(keeps%capom)) then
            keeps%capom(:) = unpack(keeps%capom(:),.not.lfill_list(:), keeps%capom(:))
            keeps%capom(:) = unpack(inserts%capom(:),lfill_list(:), keeps%capom(:))
         end if
           
         if (allocated(keeps%omega)) then
            keeps%omega(:) = unpack(keeps%omega(:),.not.lfill_list(:), keeps%omega(:))
            keeps%omega(:) = unpack(inserts%omega(:),lfill_list(:), keeps%omega(:))
         end if
           
         if (allocated(keeps%capm)) then
            keeps%capm(:)  = unpack(keeps%capm(:), .not.lfill_list(:), keeps%capm(:))
            keeps%capm(:)  = unpack(inserts%capm(:), lfill_list(:), keeps%capm(:))
         end if
            
         ! This is the base class, so will be the last to be called in the cascade. 
         keeps%nbody = size(keeps%id(:))
      end associate
     
      return
   end subroutine util_copy_fill_body


   module subroutine util_copy_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest massive body structure into an old one. 
      !! This is the inverse of a fill operation.
      implicit none
      ! Arguments
      class(swiftest_pl),    intent(inout) :: self       !! Swiftest massive body object
      class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! Internals
      integer(I4B) :: i

      associate(keeps => self)

      select type (inserts) ! The standard requires us to select the type of both arguments in order to access all the components
         class is (swiftest_pl)
         !> Spill components specific to the massive body class
            keeps%mass(:)     = unpack(keeps%mass(:),.not.lfill_list(:), keeps%mass(:))
            keeps%mass(:)     = unpack(inserts%mass(:),lfill_list(:), keeps%mass(:))
         
            keeps%Gmass(:)    = unpack(keeps%Gmass(:),.not.lfill_list(:), keeps%Gmass(:))
            keeps%Gmass(:)    = unpack(inserts%Gmass(:),lfill_list(:), keeps%Gmass(:))
         
            keeps%rhill(:)    = unpack(keeps%rhill(:),.not.lfill_list(:), keeps%rhill(:))
            keeps%rhill(:)    = unpack(inserts%rhill(:),lfill_list(:), keeps%rhill(:))
        
            if (allocated(keeps%radius) .and. allocated(inserts%radius)) then
               keeps%radius(:)   = unpack(keeps%radius(:),.not.lfill_list(:), keeps%radius(:))
               keeps%radius(:)   = unpack(inserts%radius(:),lfill_list(:), keeps%radius(:))
            end if
        
            if (allocated(keeps%density) .and. allocated(inserts%density)) then
               keeps%density(:)  = unpack(keeps%density(:),.not.lfill_list(:), keeps%density(:))
               keeps%density(:)  = unpack(inserts%density(:),lfill_list(:), keeps%density(:))
            end if

            if (allocated(keeps%k2) .and. allocated(inserts%k2)) then
               keeps%k2(:)  = unpack(keeps%k2(:),.not.lfill_list(:), keeps%k2(:))
               keeps%k2(:)  = unpack(inserts%k2(:),lfill_list(:), keeps%k2(:))
            end if

            if (allocated(keeps%Q) .and. allocated(inserts%Q)) then
               keeps%Q(:)  = unpack(keeps%Q(:),.not.lfill_list(:), keeps%Q(:))
               keeps%Q(:)  = unpack(inserts%Q(:),lfill_list(:), keeps%Q(:))
            end if

            if (allocated(keeps%tlag) .and. allocated(inserts%tlag)) then
               keeps%tlag(:)  = unpack(keeps%tlag(:),.not.lfill_list(:), keeps%tlag(:))
               keeps%tlag(:)  = unpack(inserts%tlag(:),lfill_list(:), keeps%tlag(:))
            end if

            if (allocated(keeps%xbeg) .and. allocated(inserts%xbeg)) then
               do i = 1, NDIM
                  keeps%xbeg(i, :)    = unpack(keeps%xbeg(i, :), .not.lfill_list(:), keeps%xbeg(i, :))
                  keeps%xbeg(i, :)    = unpack(inserts%xbeg(i, :), lfill_list(:), keeps%xbeg(i, :))
               end do
            end if

            if (allocated(keeps%xend) .and. allocated(inserts%xend)) then
               do i = 1, NDIM
                  keeps%xend(i, :)    = unpack(keeps%xend(i, :), .not.lfill_list(:), keeps%xend(i, :))
                  keeps%xend(i, :)    = unpack(inserts%xend(i, :), lfill_list(:), keeps%xend(i, :))
               end do
            end if

            if (allocated(keeps%vbeg) .and. allocated(inserts%vbeg)) then
               do i = 1, NDIM
                  keeps%vbeg(i, :)    = unpack(keeps%vbeg(i, :), .not.lfill_list(:), keeps%vbeg(i, :))
                  keeps%vbeg(i, :)    = unpack(inserts%vbeg(i, :), lfill_list(:), keeps%vbeg(i, :))
               end do
            end if

            if (allocated(keeps%Ip) .and. allocated(inserts%Ip)) then
               do i = 1, NDIM
                  keeps%Ip(i, :)    = unpack(keeps%Ip(i, :), .not.lfill_list(:), keeps%Ip(i, :))
                  keeps%Ip(i, :)    = unpack(inserts%Ip(i, :), lfill_list(:), keeps%Ip(i, :))
               end do
            end if

            if (allocated(keeps%rot) .and. allocated(inserts%rot)) then
               do i = 1, NDIM
                  keeps%rot(i, :)    = unpack(keeps%rot(i, :), .not.lfill_list(:), keeps%rot(i, :))
                  keeps%rot(i, :)    = unpack(inserts%rot(i, :), lfill_list(:), keeps%rot(i, :))
               end do
            end if

            keeps%ldiscard(:) = unpack(inserts%ldiscard(:), lfill_list(:), keeps%ldiscard(:))
         
            call util_copy_fill_body(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on swiftest_pl'
         end select
      end associate

      return
   end subroutine util_copy_fill_pl


   module subroutine util_copy_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      implicit none
      ! Arguments
      class(swiftest_tp),    intent(inout) :: self       !! Swiftest test particle object
      class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (swiftest_tp)
         !> Spill components specific to the test particle class
            keeps%isperi(:) = unpack(keeps%isperi(:), .not.lfill_list(:), keeps%isperi(:))
            keeps%isperi(:) = unpack(inserts%isperi(:), lfill_list(:), keeps%isperi(:))
         
            keeps%peri(:)   = unpack(keeps%peri(:),   .not.lfill_list(:), keeps%peri(:))
            keeps%peri(:)   = unpack(inserts%peri(:),   lfill_list(:), keeps%peri(:))
         
            keeps%atp(:)    = unpack(keeps%atp(:),    .not.lfill_list(:), keeps%atp(:))
            keeps%atp(:)    = unpack(inserts%atp(:),    lfill_list(:), keeps%atp(:))
         
            call util_copy_fill_body(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on swiftest_tp'
         end select
      end associate

      return
   end subroutine util_copy_fill_tp


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