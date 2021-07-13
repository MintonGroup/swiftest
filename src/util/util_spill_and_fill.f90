submodule (swiftest_classes) s_util_spill_and_fill
   use swiftest
contains
   module subroutine util_spill_body(self, discards, lspill_list)
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
         discards%a(:)        = pack(keeps%a(:),      lspill_list(:))
         discards%e(:)        = pack(keeps%e(:),      lspill_list(:))
         discards%capom(:)    = pack(keeps%capom(:),  lspill_list(:))
         discards%omega(:)    = pack(keeps%omega(:),  lspill_list(:))
         discards%capm(:)     = pack(keeps%capm(:),   lspill_list(:))
         discards%mu(:)       = pack(keeps%mu(:),     lspill_list(:))
         do i = 1, NDIM
            discards%xh(i, :) = pack(keeps%xh(i, :),     lspill_list(:))
            discards%vh(i, :) = pack(keeps%vh(i, :),     lspill_list(:))
            discards%xb(i, :) = pack(keeps%xb(i, :),     lspill_list(:))
            discards%vb(i, :) = pack(keeps%vb(i, :),     lspill_list(:))
            discards%ah(i, :) = pack(keeps%ah(i, :),     lspill_list(:))
            discards%aobl(i, :) = pack(keeps%aobl(i, :), lspill_list(:))
            discards%agr(i, :)   = pack(keeps%agr(i, :), lspill_list(:))
         end do
         if (count(.not.lspill_list(:))  > 0) then 
            keeps%id(:)         = pack(keeps%id(:),     .not. lspill_list(:))
            keeps%name(:)       = pack(keeps%name(:),   .not. lspill_list(:))
            keeps%status(:)     = pack(keeps%status(:), .not. lspill_list(:))
            keeps%a(:)          = pack(keeps%a(:),      .not. lspill_list(:))
            keeps%e(:)          = pack(keeps%e(:),      .not. lspill_list(:))
            keeps%inc(:)        = pack(keeps%inc(:),    .not. lspill_list(:))
            keeps%capom(:)      = pack(keeps%capom(:),  .not. lspill_list(:))
            keeps%omega(:)      = pack(keeps%omega(:),  .not. lspill_list(:))
            keeps%capm(:)       = pack(keeps%capm(:),   .not. lspill_list(:))
            keeps%mu(:)         = pack(keeps%mu(:),     .not. lspill_list(:))
            do i = 1, NDIM
               keeps%xh(i, :)    = pack(keeps%xh(i, :),   .not. lspill_list(:))
               keeps%vh(i, :)    = pack(keeps%vh(i, :),   .not. lspill_list(:))
               keeps%xb(i, :)    = pack(keeps%xb(i, :),   .not. lspill_list(:))
               keeps%vb(i, :)    = pack(keeps%vb(i, :),   .not. lspill_list(:))
               keeps%ah(i, :)    = pack(keeps%ah(i, :),   .not. lspill_list(:))
               keeps%aobl(i, :)  = pack(keeps%aobl(i, :), .not. lspill_list(:))
               keeps%agr(i, :)   = pack(keeps%agr(i, :),  .not. lspill_list(:))
            end do
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
      
   end subroutine util_spill_body

   module subroutine util_fill_body(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest generic particle structure into an old one. 
      !! This is the inverse of a fill operation.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self       !! Swiftest generic body object
      class(swiftest_body), intent(inout) :: inserts     !! Insertted object 
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
            
            keeps%aobl(i, :)  = unpack(keeps%aobl(i, :), .not.lfill_list(:), keeps%aobl(i, :))
            keeps%aobl(i, :)  = unpack(inserts%aobl(i, :), lfill_list(:), keeps%aobl(i, :))

            keeps%agr(i, :)  = unpack(keeps%agr(i, :), .not.lfill_list(:), keeps%agr(i, :))
            keeps%agr(i, :)  = unpack(inserts%agr(i, :), lfill_list(:), keeps%agr(i, :))
         end do
         
         keeps%a(:)     = unpack(keeps%a(:),    .not.lfill_list(:), keeps%a(:))
         keeps%a(:)     = unpack(inserts%a(:),    lfill_list(:), keeps%a(:))
            
         keeps%e(:)     = unpack(keeps%e(:),    .not.lfill_list(:), keeps%e(:))
         keeps%e(:)     = unpack(inserts%e(:),    lfill_list(:), keeps%e(:))
            
         keeps%inc(:)   = unpack(keeps%inc(:),  .not.lfill_list(:), keeps%inc(:))
         keeps%inc(:)   = unpack(inserts%inc(:),  lfill_list(:), keeps%inc(:))
            
         keeps%capom(:) = unpack(keeps%capom(:),.not.lfill_list(:), keeps%capom(:))
         keeps%capom(:) = unpack(inserts%capom(:),lfill_list(:), keeps%capom(:))
            
         keeps%omega(:) = unpack(keeps%omega(:),.not.lfill_list(:), keeps%omega(:))
         keeps%omega(:) = unpack(inserts%omega(:),lfill_list(:), keeps%omega(:))
            
         keeps%capm(:)  = unpack(keeps%capm(:), .not.lfill_list(:), keeps%capm(:))
         keeps%capm(:)  = unpack(inserts%capm(:), lfill_list(:), keeps%capm(:))
            
         keeps%mu(:)    = unpack(keeps%mu(:),   .not.lfill_list(:), keeps%mu(:))
         keeps%mu(:)    = unpack(inserts%mu(:),   lfill_list(:), keeps%mu(:))
            

         ! This is the base class, so will be the last to be called in the cascade. 
         keeps%nbody = size(keeps%id(:))
      end associate
      
   end subroutine util_fill_body

   module subroutine util_spill_pl(self, discards, lspill_list)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest massive body structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_pl),    intent(inout) :: self        !! Swiftest massive body object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discardse
      ! Internals
      integer(I4B) :: i

      associate(keeps => self)

         select type (discards) ! The standard requires us to select the type of both arguments in order to access all the components
         class is (swiftest_pl)
         !> Spill components specific to the massive body class
            discards%mass(:)     = pack(keeps%mass(:),    lspill_list(:))
            discards%Gmass(:)    = pack(keeps%Gmass(:),   lspill_list(:))
            discards%rhill(:)    = pack(keeps%rhill(:),   lspill_list(:))
            discards%radius(:)   = pack(keeps%radius(:),  lspill_list(:))
            discards%density(:)  = pack(keeps%density(:), lspill_list(:))
            discards%k2(:)       = pack(keeps%k2(:),       lspill_list(:))
            discards%Q(:)        = pack(keeps%Q(:),        lspill_list(:))
            discards%tlag(:)      = pack(keeps%tlag(:),      lspill_list(:))
            do i = 1, NDIM
               discards%Ip(i, :) = pack(keeps%Ip(i, :),     lspill_list(:))
               discards%Ip(i, :)  = pack(keeps%Ip(i, :),     lspill_list(:))
            end do
            if (count(.not.lspill_list(:))  > 0) then 
               keeps%mass(:)        = pack(keeps%mass(:),    .not. lspill_list(:))
               keeps%Gmass(:)       = pack(keeps%Gmass(:),   .not. lspill_list(:))
               keeps%rhill(:)       = pack(keeps%rhill(:),   .not. lspill_list(:))
               keeps%radius(:)      = pack(keeps%radius(:),  .not. lspill_list(:))
               keeps%density(:)     = pack(keeps%density(:), .not. lspill_list(:))
            end if

            call util_spill_body(keeps, discards, lspill_list)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on swiftest_pl'
         end select
      end associate
      return
   end subroutine util_spill_pl

   module subroutine util_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest massive body structure into an old one. 
      !! This is the inverse of a fill operation.
      implicit none
      ! Arguments
      class(swiftest_pl),    intent(inout) :: self       !! Swiftest massive body object
      class(swiftest_body),  intent(inout) :: inserts    !! Swiftest body object to be inserted
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
         
            keeps%radius(:)   = unpack(keeps%radius(:),.not.lfill_list(:), keeps%radius(:))
            keeps%radius(:)   = unpack(inserts%radius(:),lfill_list(:), keeps%radius(:))
         
            keeps%density(:)  = unpack(keeps%density(:),.not.lfill_list(:), keeps%density(:))
            keeps%density(:)  = unpack(inserts%density(:),lfill_list(:), keeps%density(:))

            keeps%k2(:)  = unpack(keeps%k2(:),.not.lfill_list(:), keeps%k2(:))
            keeps%k2(:)  = unpack(inserts%k2(:),lfill_list(:), keeps%k2(:))

            keeps%Q(:)  = unpack(keeps%Q(:),.not.lfill_list(:), keeps%Q(:))
            keeps%Q(:)  = unpack(inserts%Q(:),lfill_list(:), keeps%Q(:))

            keeps%tlag(:)  = unpack(keeps%tlag(:),.not.lfill_list(:), keeps%tlag(:))
            keeps%tlag(:)  = unpack(inserts%tlag(:),lfill_list(:), keeps%tlag(:))

            do i = 1, NDIM
               keeps%Ip(i, :)    = unpack(keeps%Ip(i, :), .not.lfill_list(:), keeps%Ip(i, :))
               keeps%Ip(i, :)    = unpack(inserts%Ip(i, :), lfill_list(:), keeps%Ip(i, :))

               keeps%Ip(i, :)    = unpack(keeps%Ip(i, :), .not.lfill_list(:), keeps%Ip(i, :))
               keeps%rot(i, :)    = unpack(inserts%rot(i, :), lfill_list(:), keeps%rot(i, :))
            end do
            keeps%ldiscard(:) = unpack(inserts%ldiscard(:), lfill_list(:), keeps%ldiscard(:))
         
            call util_fill_body(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on swiftest_pl'
         end select
      end associate
      return
   end subroutine util_fill_pl

   module subroutine util_spill_tp(self, discards, lspill_list)
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
         call util_spill_body(keeps, discards, lspill_list)
      class default
         write(*,*) 'Error! spill method called for incompatible return type on swiftest_tp'
      end select
      end associate
      return
   end subroutine util_spill_tp

   module subroutine util_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      implicit none
      ! Arguments
      class(swiftest_tp),    intent(inout) :: self       !! Swiftest test particle object
      class(swiftest_body),  intent(inout) :: inserts    !! Swiftest body object to be inserted
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
      
         call util_fill_body(keeps, inserts, lfill_list)
      class default
         write(*,*) 'Error! fill method called for incompatible return type on swiftest_tp'
      end select
      end associate
      return
   end subroutine util_fill_tp

end submodule s_util_spill_and_fill



      


