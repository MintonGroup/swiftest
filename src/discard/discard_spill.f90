submodule (swiftest_classes) s_discard_spill
!! This submodule contains the methods spill, spill_body, spill_pl, and spill_tp for basic swiftest particles
contains
   module procedure discard_spill_body
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest particle structure from active list to discard list
      use swiftest    
      implicit none
      integer(I4B) :: i

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      associate(n => keeps%nbody)
         discards%name(:)    = pack(keeps%name(1:n),          lspill_list(1:n))
         keeps%name(:)       = pack(keeps%name(1:n),    .not. lspill_list(1:n))

         discards%status(:)  = pack(keeps%status(1:n),        lspill_list(1:n))
         keeps%status(:)     = pack(keeps%status(1:n),  .not. lspill_list(1:n))

         do concurrent (i = 1:NDIM)
            discards%xh(:, i) = pack(keeps%xh(1:n, i),          lspill_list(1:n))
            keeps%xh(:, i)    = pack(keeps%xh(1:n, i),    .not. lspill_list(1:n))

            discards%vh(:, i) = pack(keeps%vh(1:n, i),          lspill_list(1:n))
            keeps%vh(:, i)    = pack(keeps%vh(1:n, i),    .not. lspill_list(1:n))

            discards%xb(:, i) = pack(keeps%xb(1:n, i),          lspill_list(1:n))
            keeps%xb(:, i)    = pack(keeps%xb(1:n, i),    .not. lspill_list(1:n))

            discards%vb(:, i) = pack(keeps%vb(1:n, i),          lspill_list(1:n))
            keeps%vb(:, i)    = pack(keeps%vb(1:n, i),    .not. lspill_list(1:n))
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

         discards%mu_vec(:) = pack(keeps%mu_vec(1:n),         lspill_list(1:n))
         keeps%mu_vec(:)    = pack(keeps%mu_vec(1:n),   .not. lspill_list(1:n))

         discards%dt_vec(:) = pack(keeps%dt_vec(1:n),         lspill_list(1:n))
         keeps%dt_vec(:)    = pack(keeps%dt_vec(1:n),   .not. lspill_list(1:n))

         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nbody values for both the keeps and discareds
         discards%nbody = count(lspill_list(1:n))
         keeps%nbody = n - discards%nbody
      end associate
      
   end procedure discard_spill_body

   module procedure discard_spill_pl
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest massive body particle structure from active list to discard list
      use swiftest
      implicit none
      integer(I4B) :: i

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      associate(npl => keeps%nbody)
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

         do i = 1, NDIM
            discards%Ip(:, i)  = pack(keeps%Ip(1:npl, i),          lspill_list(1:npl))
            keeps%Ip(:, i)     = pack(keeps%Ip(1:npl, i),    .not. lspill_list(1:npl))

            discards%rot(:, i) = pack(keeps%rot(1:npl, i),         lspill_list(1:npl))
            keeps%rot(:, i)    = pack(keeps%rot(1:npl, i),   .not. lspill_list(1:npl))
         end do
         
         discards%k2(:)       = pack(keeps%k2(1:npl),            lspill_list(1:npl))
         keeps%k2(:)          = pack(keeps%k2(1:npl),      .not. lspill_list(1:npl))
         
         discards%Q(:)        = pack(keeps%Q(1:npl),              lspill_list(1:npl))
         keeps%Q(:)           = pack(keeps%Q(1:npl),        .not. lspill_list(1:npl))
         
      end associate

      return
   end procedure discard_spill_pl

   module procedure discard_spill_tp
      !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
      use swiftest
      implicit none

      integer(I4B) :: i

      associate(ntp => keeps%nbody)

         discards%isperi(:) = pack(keeps%isperi(1:ntp),       lspill_list(1:ntp))
         keeps%isperi(:)    = pack(keeps%isperi(1:ntp), .not. lspill_list(1:ntp))

         discards%peri(:)   = pack(keeps%peri(1:ntp),         lspill_list(1:ntp))
         keeps%peri(:)      = pack(keeps%peri(1:ntp),   .not. lspill_list(1:ntp))

         discards%atp(:)    = pack(keeps%atp(1:ntp),          lspill_list(1:ntp))
         keeps%atp(:)       = pack(keeps%atp(1:ntp),    .not. lspill_list(1:ntp))

      end associate
      return
   end procedure discard_spill_tp
end submodule s_discard_spill



      


