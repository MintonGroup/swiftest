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
         discards%name(:)    = pack(keeps%name(1:n),          lspill_list(1:n))
         keeps%name(:)       = pack(keeps%name(1:n),    .not. lspill_list(1:n))

         discards%status(:)  = pack(keeps%status(1:n),        lspill_list(1:n))
         keeps%status(:)     = pack(keeps%status(1:n),  .not. lspill_list(1:n))

         discards%ldiscard(:)  = pack(keeps%ldiscard(1:n),        lspill_list(1:n))
         keeps%ldiscard(:)     = pack(keeps%ldiscard(1:n),  .not. lspill_list(1:n))

         !do concurrent (i = 1:NDIM)
         do i = 1, NDIM
            discards%xh(i, :) = pack(keeps%xh(i, 1:n),          lspill_list(1:n))
            keeps%xh(i, :)    = pack(keeps%xh(i, 1:n),    .not. lspill_list(1:n))

            discards%vh(i, :) = pack(keeps%vh(i, 1:n),          lspill_list(1:n))
            keeps%vh(i, :)    = pack(keeps%vh(i, 1:n),    .not. lspill_list(1:n))

            discards%xb(i, :) = pack(keeps%xb(i, 1:n),          lspill_list(1:n))
            keeps%xb(i, :)    = pack(keeps%xb(i, 1:n),    .not. lspill_list(1:n))

            discards%vb(i, :) = pack(keeps%vb(i, 1:n),          lspill_list(1:n))
            keeps%vb(i, :)    = pack(keeps%vb(i, 1:n),    .not. lspill_list(1:n))

            discards%ah(i, :) = pack(keeps%ah(i, 1:n),          lspill_list(1:n))
            keeps%ah(i, :)    = pack(keeps%ah(i, 1:n),    .not. lspill_list(1:n))

            discards%aobl(i, :) = pack(keeps%aobl(i, 1:n),          lspill_list(1:n))
            keeps%aobl(i, :)    = pack(keeps%aobl(i, 1:n),    .not. lspill_list(1:n))
         end do
         
         discards%a(:)       = pack(keeps%a(1:n),             lspill_list(1:n))
         keeps%a(:)          = pack(keeps%a(1:n),       .not. lspill_list(1:n))
         
         discards%e(:)       = pack(keeps%e(1:n),             lspill_list(1:n))
         keeps%e(:)          = pack(keeps%e(1:n),       .not. lspill_list(1:n))
         
         discards%inc(:)     = pack(keeps%inc(1:n),           lspill_list(1:n))
         keeps%inc(:)        = pack(keeps%inc(1:n),     .not. lspill_list(1:n))
         
         discards%capom(:)   = pack(keeps%capom(1:n),         lspill_list(1:n))
         keeps%capom(:)      = pack(keeps%capom(1:n),   .not. lspill_list(1:n))
         
         discards%omega(:)   = pack(keeps%omega(1:n),         lspill_list(1:n))
         keeps%omega(:)      = pack(keeps%omega(1:n),   .not. lspill_list(1:n))
         
         discards%capm(:)    = pack(keeps%capm(1:n),          lspill_list(1:n))
         keeps%capm(:)       = pack(keeps%capm(1:n),    .not. lspill_list(1:n))

         discards%mu(:) = pack(keeps%mu(1:n),         lspill_list(1:n))
         keeps%mu(:)    = pack(keeps%mu(1:n),   .not. lspill_list(1:n))

         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nbody values for both the keeps and discareds
         discards%nbody = count(lspill_list(1:n))
         keeps%nbody = n - discards%nbody
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
         keeps%name(:)     = merge(inserts%name(:), keeps%name(:), lfill_list(:))
         keeps%status(:)   = merge(inserts%status(:), keeps%status(:),   lfill_list(:))
         keeps%ldiscard(:) = merge(inserts%ldiscard(:), keeps%ldiscard(:),   lfill_list(:))

         do i = 1, NDIM
            keeps%xh(i, :)    = merge(inserts%xh(i, :), keeps%xh(i, :),     lfill_list(:))
            keeps%vh(i, :)    = merge(inserts%vh(i, :), keeps%vh(i, :),     lfill_list(:))
            keeps%xb(i, :)    = merge(inserts%xb(i, :), keeps%xb(i, :),     lfill_list(:))
            keeps%vb(i, :)    = merge(inserts%vb(i, :), keeps%vb(i, :),     lfill_list(:))
            keeps%ah(i, :)    = merge(inserts%ah(i, :), keeps%ah(i, :),     lfill_list(:))
            keeps%aobl(i, :)  = merge(inserts%aobl(i, :), keeps%aobl(i, :),     lfill_list(:))
         end do
         
         keeps%a(:)     = merge(inserts%a(:),     keeps%a(:),     lfill_list(:))
         keeps%e(:)     = merge(inserts%e(:),     keeps%e(:),     lfill_list(:))
         keeps%inc(:)   = merge(inserts%inc(:),   keeps%inc(:),   lfill_list(:))
         keeps%capom(:) = merge(inserts%capom(:), keeps%capom(:), lfill_list(:))
         keeps%omega(:) = merge(inserts%omega(:), keeps%omega(:), lfill_list(:))
         keeps%capm(:)  = merge(inserts%capm(:),  keeps%capm(:),  lfill_list(:))
         keeps%mu(:)    = merge(inserts%mu(:),    keeps%mu(:),    lfill_list(:))

         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nbody values for both the keeps and discareds
         keeps%nbody = count(lfill_list(:))
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
               keeps%mass(:)        = pack(keeps%mass(1:npl),    .not. lspill_list(1:npl))
   
               discards%Gmass(:)    = pack(keeps%Gmass(1:npl),         lspill_list(1:npl))
               keeps%Gmass(:)       = pack(keeps%Gmass(1:npl),   .not. lspill_list(1:npl))
   
               discards%rhill(:)    = pack(keeps%rhill(1:npl),         lspill_list(1:npl))
               keeps%rhill(:)       = pack(keeps%rhill(1:npl),   .not. lspill_list(1:npl))
   
               discards%radius(:)   = pack(keeps%radius(1:npl),        lspill_list(1:npl))
               keeps%radius(:)      = pack(keeps%radius(1:npl),  .not. lspill_list(1:npl))
   
               discards%density(:)  = pack(keeps%density(1:npl),       lspill_list(1:npl))
               keeps%density(:)     = pack(keeps%density(1:npl), .not. lspill_list(1:npl))
   
               !do concurrent (i = 1:NDIM)
               do i = 1, NDIM
                  discards%Ip(i, :)  = pack(keeps%Ip(i, 1:npl),          lspill_list(1:npl))
                  keeps%Ip(i, :)     = pack(keeps%Ip(i, 1:npl),    .not. lspill_list(1:npl))
   
                  discards%rot(i, :) = pack(keeps%rot(i, 1:npl),         lspill_list(1:npl))
                  keeps%rot(i, :)    = pack(keeps%rot(i, 1:npl),   .not. lspill_list(1:npl))
               end do
               
               discards%k2(:)       = pack(keeps%k2(1:npl),            lspill_list(1:npl))
               keeps%k2(:)          = pack(keeps%k2(1:npl),      .not. lspill_list(1:npl))
               
               discards%Q(:)        = pack(keeps%Q(1:npl),              lspill_list(1:npl))
               keeps%Q(:)           = pack(keeps%Q(1:npl),        .not. lspill_list(1:npl))

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
               keeps%mass(:)     = merge(inserts%mass(:), keeps%mass(:),     lfill_list(:))
               keeps%Gmass(:)    = merge(inserts%Gmass(:), keeps%Gmass(:),    lfill_list(:))
               keeps%rhill(:)    = merge(inserts%rhill(:), keeps%rhill(:),    lfill_list(:))
               keeps%radius(:)   = merge(inserts%radius(:), keeps%radius(:),   lfill_list(:))
               keeps%density(:)  = merge(inserts%density(:), keeps%density(:),  lfill_list(:))
               do i = 1, NDIM
                  keeps%Ip(i, :)   = merge(inserts%Ip(i, :), keeps%Ip(i, :),     lfill_list(:))
                  keeps%rot(i, :)  = merge(inserts%rot(i, :), keeps%rot(i, :),    lfill_list(:))
               end do
               keeps%k2(:)       = merge(inserts%k2(:), keeps%k2(:),       lfill_list(:))
               keeps%Q(:)        = merge(inserts%Q(:), keeps%Q(:),         lfill_list(:))
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
            keeps%isperi(:)    = pack(keeps%isperi(1:ntp), .not. lspill_list(1:ntp))

            discards%peri(:)   = pack(keeps%peri(1:ntp),         lspill_list(1:ntp))
            keeps%peri(:)      = pack(keeps%peri(1:ntp),   .not. lspill_list(1:ntp))

            discards%atp(:)    = pack(keeps%atp(1:ntp),          lspill_list(1:ntp))
            keeps%atp(:)       = pack(keeps%atp(1:ntp),    .not. lspill_list(1:ntp))
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
               keeps%isperi(:) = merge(inserts%isperi(:), keeps%isperi(:), lfill_list(:))
               keeps%peri(:)   = merge(inserts%peri(:),   keeps%peri(:),   lfill_list(:))
               keeps%atp(:)    = merge(inserts%atp(:),    keeps%atp(:),    lfill_list(:))
               call util_fill_body(keeps, inserts, lfill_list)
            class default
               write(*,*) 'Error! fill method called for incompatible return type on swiftest_tp'
            end select
         end associate
         return
         end procedure util_fill_tp

end submodule s_util_spill_and_fill



      


