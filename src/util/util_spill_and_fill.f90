submodule (swiftest_classes) s_util_spill_and_fill
!! This submodule contains the methods spill, spill_body, spill_pl, and spill_tp for basic swiftest particles
contains
   module procedure util_spill_body
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest generic particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      use swiftest    
      implicit none
      integer(I4B) :: i

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self, n => self%nbody)
         discards%name(:)     = pack(keeps%name(1:n),          lspill_list(1:n))
         discards%status(:)   = pack(keeps%status(1:n),        lspill_list(1:n))
         discards%ldiscard(:) = pack(keeps%ldiscard(1:n),      lspill_list(1:n))
         discards%a(:)        = pack(keeps%a(1:n),             lspill_list(1:n))
         discards%e(:)        = pack(keeps%e(1:n),             lspill_list(1:n))
         discards%capom(:)    = pack(keeps%capom(1:n),         lspill_list(1:n))
         discards%omega(:)    = pack(keeps%omega(1:n),         lspill_list(1:n))
         discards%capm(:)     = pack(keeps%capm(1:n),          lspill_list(1:n))
         discards%mu(:)       = pack(keeps%mu(1:n),            lspill_list(1:n))
         do i = 1, NDIM
            discards%xh(i, :) = pack(keeps%xh(i, 1:n),          lspill_list(1:n))
            discards%vh(i, :) = pack(keeps%vh(i, 1:n),          lspill_list(1:n))
            discards%xb(i, :) = pack(keeps%xb(i, 1:n),          lspill_list(1:n))
            discards%vb(i, :) = pack(keeps%vb(i, 1:n),          lspill_list(1:n))
            discards%ah(i, :) = pack(keeps%ah(i, 1:n),          lspill_list(1:n))
            discards%aobl(i, :) = pack(keeps%aobl(i, 1:n),          lspill_list(1:n))
         end do
         if (count(.not.lspill_list(1:n))  > 0) then 
            keeps%name(:)       = pack(keeps%name(1:n),     .not. lspill_list(1:n))
            keeps%status(:)     = pack(keeps%status(1:n),   .not. lspill_list(1:n))
            keeps%ldiscard(:)   = pack(keeps%ldiscard(1:n), .not. lspill_list(1:n))
            keeps%a(:)          = pack(keeps%a(1:n),        .not. lspill_list(1:n))
            keeps%e(:)          = pack(keeps%e(1:n),        .not. lspill_list(1:n))
            keeps%inc(:)        = pack(keeps%inc(1:n),      .not. lspill_list(1:n))
            keeps%capom(:)      = pack(keeps%capom(1:n),    .not. lspill_list(1:n))
            keeps%omega(:)      = pack(keeps%omega(1:n),    .not. lspill_list(1:n))
            keeps%capm(:)       = pack(keeps%capm(1:n),     .not. lspill_list(1:n))
            keeps%mu(:)         = pack(keeps%mu(1:n),       .not. lspill_list(1:n))
            do i = 1, NDIM
               keeps%xh(i, :)    = pack(keeps%xh(i, 1:n),    .not. lspill_list(1:n))
               keeps%vh(i, :)    = pack(keeps%vh(i, 1:n),    .not. lspill_list(1:n))
               keeps%xb(i, :)    = pack(keeps%xb(i, 1:n),    .not. lspill_list(1:n))
               keeps%vb(i, :)    = pack(keeps%vb(i, 1:n),    .not. lspill_list(1:n))
               keeps%ah(i, :)    = pack(keeps%ah(i, 1:n),    .not. lspill_list(1:n))
               keeps%aobl(i, :)    = pack(keeps%aobl(i, 1:n),    .not. lspill_list(1:n))
            end do
         end if
         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nbody values for both the keeps and discareds
         discards%nbody = count(lspill_list(1:n))
         keeps%nbody = count(.not.lspill_list(1:n)) 

      end associate
      
   end procedure util_spill_body

   module procedure util_fill_body
      !! author: David A. Minton
      !!
      !! Insert new Swiftest generic particle structure into an old one. 
      !! This is the inverse of a fill operation.
      use swiftest    
      implicit none
      integer(I4B) :: i

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         keeps%name(:)     = unpack(inserts%name(:), lfill_list(:), keeps%name(:))
         keeps%status(:)   = unpack(inserts%status(:), lfill_list(:), keeps%status(:))
         keeps%ldiscard(:) = unpack(inserts%ldiscard(:), lfill_list(:), keeps%ldiscard(:))

         do i = 1, NDIM
            keeps%xh(i, :)    = unpack(inserts%xh(i, :), lfill_list(:), keeps%xh(i, :))
            keeps%vh(i, :)    = unpack(inserts%vh(i, :), lfill_list(:), keeps%vh(i, :))
            keeps%xb(i, :)    = unpack(inserts%xb(i, :), lfill_list(:), keeps%xb(i, :))
            keeps%vb(i, :)    = unpack(inserts%vb(i, :), lfill_list(:), keeps%vb(i, :))
            keeps%ah(i, :)    = unpack(inserts%ah(i, :), lfill_list(:), keeps%ah(i, :))
            keeps%aobl(i, :)  = unpack(inserts%aobl(i, :), lfill_list(:), keeps%aobl(i, :))
         end do
         
         keeps%a(:)     = unpack(inserts%a(:),    lfill_list(:), keeps%a(:))
         keeps%e(:)     = unpack(inserts%e(:),    lfill_list(:), keeps%e(:))
         keeps%inc(:)   = unpack(inserts%inc(:),  lfill_list(:), keeps%inc(:))
         keeps%capom(:) = unpack(inserts%capom(:),lfill_list(:), keeps%capom(:))
         keeps%omega(:) = unpack(inserts%omega(:),lfill_list(:), keeps%omega(:))
         keeps%capm(:)  = unpack(inserts%capm(:), lfill_list(:), keeps%capm(:))
         keeps%mu(:)    = unpack(inserts%mu(:),   lfill_list(:), keeps%mu(:))

         ! This is the base class, so will be the last to be called in the cascade. 
         keeps%nbody = count(keeps%status(:) == ACTIVE) 
      end associate
      
   end procedure util_fill_body

   module procedure util_spill_pl
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest massive body structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      use swiftest
      implicit none

      integer(I4B) :: i

      associate(keeps => self, npl => self%nbody)

         select type (discards) ! The standard requires us to select the type of both arguments in order to access all the components
            class is (swiftest_pl)
            !> Spill components specific to the massive body class
               discards%mass(:)     = pack(keeps%mass(1:npl),          lspill_list(1:npl))
               discards%Gmass(:)    = pack(keeps%Gmass(1:npl),         lspill_list(1:npl))
               discards%rhill(:)    = pack(keeps%rhill(1:npl),         lspill_list(1:npl))
               discards%radius(:)   = pack(keeps%radius(1:npl),        lspill_list(1:npl))
               discards%density(:)  = pack(keeps%density(1:npl),       lspill_list(1:npl))
               discards%k2(:)       = pack(keeps%k2(1:npl),            lspill_list(1:npl))
               discards%Q(:)        = pack(keeps%Q(1:npl),              lspill_list(1:npl))
               do i = 1, NDIM
                  discards%Ip(i, :)  = pack(keeps%Ip(i, 1:npl),          lspill_list(1:npl))
                  discards%rot(i, :) = pack(keeps%rot(i, 1:npl),         lspill_list(1:npl))
               end do
               if (count(.not.lspill_list(1:npl))  > 0) then 
                  keeps%mass(:)        = pack(keeps%mass(1:npl),    .not. lspill_list(1:npl))
                  keeps%Gmass(:)       = pack(keeps%Gmass(1:npl),   .not. lspill_list(1:npl))
                  keeps%rhill(:)       = pack(keeps%rhill(1:npl),   .not. lspill_list(1:npl))
                  keeps%radius(:)      = pack(keeps%radius(1:npl),  .not. lspill_list(1:npl))
                  keeps%density(:)     = pack(keeps%density(1:npl), .not. lspill_list(1:npl))
                  keeps%k2(:)          = pack(keeps%k2(1:npl),      .not. lspill_list(1:npl))
                  keeps%Q(:)           = pack(keeps%Q(1:npl),        .not. lspill_list(1:npl))
                  do i = 1, NDIM
                     keeps%Ip(i, :)     = pack(keeps%Ip(i, 1:npl),    .not. lspill_list(1:npl))
                     keeps%rot(i, :)    = pack(keeps%rot(i, 1:npl),   .not. lspill_list(1:npl))
                  end do
               end if

               call util_spill_body(keeps, discards, lspill_list)
            class default
               write(*,*) 'Error! spill method called for incompatible return type on swiftest_pl'
            end select
         end associate
         return
   
   end procedure util_spill_pl

   module procedure util_fill_pl
      !! author: David A. Minton
      !!
      !! Insert new Swiftest massive body structure into an old one. 
      !! This is the inverse of a fill operation.
      use swiftest
      implicit none

      integer(I4B) :: i

      associate(keeps => self)

         select type (inserts) ! The standard requires us to select the type of both arguments in order to access all the components
            class is (swiftest_pl)
            !> Spill components specific to the massive body class
               keeps%mass(:)     = unpack(inserts%mass(:),lfill_list(:), keeps%mass(:))
               keeps%Gmass(:)    = unpack(inserts%Gmass(:),lfill_list(:), keeps%Gmass(:))
               keeps%rhill(:)    = unpack(inserts%rhill(:),lfill_list(:), keeps%rhill(:))
               keeps%radius(:)   = unpack(inserts%radius(:),lfill_list(:), keeps%radius(:))
               keeps%density(:)  = unpack(inserts%density(:),lfill_list(:), keeps%density(:))
               do i = 1, NDIM
                  keeps%Ip(i, :)   = unpack(inserts%Ip(i, :), lfill_list(:), keeps%Ip(i, :))
                  keeps%rot(i, :)  = unpack(inserts%rot(i, :), lfill_list(:), keeps%rot(i, :))
               end do
               keeps%k2(:)       = unpack(inserts%k2(:), lfill_list(:), keeps%k2(:))
               keeps%Q(:)        = unpack(inserts%Q(:), lfill_list(:), keeps%Q(:))
               call util_fill_body(keeps, inserts, lfill_list)
            class default
               write(*,*) 'Error! fill method called for incompatible return type on swiftest_pl'
            end select
         end associate
         return
   
   end procedure util_fill_pl

   module procedure util_spill_tp
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      use swiftest
      implicit none

      associate(keeps => self, ntp => self%nbody)
      select type(discards)
         class is (swiftest_tp)
         !> Spill components specific to the test particle class
            discards%isperi(:) = pack(keeps%isperi(1:ntp),       lspill_list(1:ntp))
            discards%peri(:)   = pack(keeps%peri(1:ntp),         lspill_list(1:ntp))
            discards%atp(:)    = pack(keeps%atp(1:ntp),          lspill_list(1:ntp))
            if (count(.not.lspill_list(1:ntp))  > 0) then 
               keeps%atp(:)       = pack(keeps%atp(1:ntp),    .not. lspill_list(1:ntp))
               keeps%peri(:)      = pack(keeps%peri(1:ntp),   .not. lspill_list(1:ntp))
               keeps%isperi(:)    = pack(keeps%isperi(1:ntp), .not. lspill_list(1:ntp))
            end if
            call util_spill_body(keeps, discards, lspill_list)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on swiftest_tp'
         end select
      end associate
      return
      end procedure util_spill_tp

      module procedure util_fill_tp
         !! author: David A. Minton
         !!
         !! Insert new Swiftest test particle structure into an old one. 
         !! This is the inverse of a fill operation.
         use swiftest
         implicit none
   
         associate(keeps => self)
         select type(inserts)
            class is (swiftest_tp)
            !> Spill components specific to the test particle class
               keeps%isperi(:) = unpack(inserts%isperi(:), lfill_list(:), keeps%isperi(:))
               keeps%peri(:)   = unpack(inserts%peri(:),   lfill_list(:), keeps%peri(:))
               keeps%atp(:)    = unpack(inserts%atp(:),    lfill_list(:), keeps%atp(:))
               call util_fill_body(keeps, inserts, lfill_list)
            class default
               write(*,*) 'Error! fill method called for incompatible return type on swiftest_tp'
            end select
         end associate
         return
         end procedure util_fill_tp

end submodule s_util_spill_and_fill



      


