submodule(whm_classes) s_whm_getacch
contains
   module procedure whm_getacch_pl  ! (self, cb, config, t)
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of planets
   !!
   !! Adapted from Hal Levison's Swift routine getacch.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_getacch.f90
   use swiftest
   implicit none
   integer(I4B)                     :: i
   logical, save :: lmalloc = .false.  
   real(DP), dimension(:), allocatable, save :: r2, fac
   real(DP), dimension(:), allocatable, save    :: irh, irj, ir3h, ir3j
   real(DP), dimension(NDIM) :: ah0

   associate(pl => self, npl => self%nbody, Gmpl => self%Gmass, xj => self%xj, xh => self%xh, &
             ah => self%ah, ah1 => self%ah1, ah2 => self%ah2, ah3 => self%ah3, &
             j2rp2 => cb%j2rp2, j4rp4 => cb%j4rp4, aobl => self%aobl, aobl0 => cb%aobl)

      if (.not. allocated(r2)) allocate(r2(npl))
      if (.not. allocated(irj)) allocate(irj(npl))
      if (.not. allocated(ir3j)) allocate(ir3j(npl))
      if (.not. allocated(irh)) allocate(irh(npl))
      if (.not. allocated(ir3h)) allocate(ir3h(npl))
      if (.not. allocated(fac)) allocate(fac(npl))
      r2(:) = xj(1:npl, :) .dot. xj(1:npl, :)
      irj(:)= 1.0_DP / sqrt(r2(1:npl))
      ir3j(:) = irj(1:npl) / r2(1:npl)

      r2(:) = xh(1:npl, :) .dot. xh(1:npl, :)
      irh(:)= 1.0_DP / sqrt(r2(1:npl))
      ir3h(:) = irh(1:npl) / r2(1:npl)

      ah0(:) = 0.0_DP

      fac(2:npl) = Gmpl(2:npl) * ir3h(2:npl)
      do concurrent (i = 1:NDIM)
         ah0(i) = - sum(fac(2:npl) * xh(2:npl, i))
      end do
      call whm_getacch_ah1(cb, pl, ir3h, ir3j)
      call whm_getacch_ah2(cb, pl, ir3j)
      call whm_getacch_ah3(pl)
      
      do concurrent (i = 1:NDIM)
         ah(1:npl, i) = ah0(i) + ah1(1:npl, i) + ah2(1:npl, i) + ah3(1:npl, i)
      end do

      if (j2rp2 /= 0.0_DP) then
         call  self%obl_acc(cb, irh, xh)

         do concurrent (i = 1:NDIM)
            ah(1:npl, i) = ah(1:npl, i) + aobl(1:npl, i) - aobl0(i)
         end do
      end if
      if (config%lextra_force) call pl%user_getacch(cb, config, t)
      if (config%lgr) call pl%gr_getacch(cb, config) 

   end associate
   return

   contains
      pure subroutine whm_getacch_ah1(cb, pl, ir3h, ir3j)
         !! author: David A. Minton
         !!
         !! Compute first term heliocentric accelerations of planets
         !!
         !! Adapted from Hal Levison's Swift routine getacch_ah1.f
         !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah1.f90
         use swiftest
         implicit none
         class(whm_central_body), intent(in) :: cb
         class(whm_pl), intent(inout)        :: pl
         real(DP), dimension(:), intent(in) :: ir3h, ir3j

         integer(I4B)              :: i
         real(DP), dimension(NDIM) :: ah1h, ah1j

         associate(npl => pl%nbody, msun => cb%Gmass, xh => pl%xh, xj => pl%xj, ah1 => pl%ah1)
            ah1(1:npl, :) = 0.0_DP
            do concurrent (i = 2:npl)
               ah1j(:) = xj(i, :) * ir3j(i)
               ah1h(:) = xh(i, :) * ir3h(i)
               ah1(i,:) = msun * (ah1j(:) - ah1h(:))
            end do
         end associate
      
         return
      
      end subroutine whm_getacch_ah1

      pure subroutine whm_getacch_ah2(cb, pl, ir3j) 
         !! author: David A. Minton
         !!
         !! Compute second term heliocentric accelerations of planets
         !!
         !! Adapted from Hal Levison's Swift routine getacch_ah2.f
         !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah2.f90
         use swiftest
         implicit none

         class(whm_central_body), intent(in)    :: cb
         class(whm_pl),           intent(inout) :: pl
         real(DP), dimension(:),  intent(in)    :: ir3j
         integer(I4B)                           :: i
         real(DP)                               :: etaj, fac
     
         associate(npl => pl%nbody, Gmsun => cb%Gmass, xh => pl%xh, xj => pl%xj, ah2 => pl%ah2, Gmpl => pl%Gmass)
            ah2(1:npl, :) = 0.0_DP
            etaj = Gmsun
            do i = 2, npl
               etaj = etaj + Gmpl(i - 1)
               fac = Gmpl(i) * Gmsun * ir3j(i) / etaj
               ah2(i, :) = ah2(i - 1, :) + fac * xj(i, :)
            end do
         end associate
      
         return
      end subroutine whm_getacch_ah2

      pure subroutine whm_getacch_ah3(pl)
         !! author: David A. Minton
         !!
         !! Compute direct cross (third) term heliocentric accelerations of planets
         !!
         !! Adapted from Hal Levison's Swift routine getacch_ah3.f
         !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
         use swiftest
         implicit none

         class(whm_pl),           intent(inout) :: pl
         integer(I4B)                           :: i, j
         real(DP)                               :: rji2, irij3, faci, facj
         real(DP), dimension(NDIM)              :: dx
     
         associate(npl => pl%nbody, xh => pl%xh, ah3 => pl%ah3, Gmpl => pl%Gmass) 
            ah3(1:npl,:) = 0.0_DP

            do i = 1, npl - 1
               do j = i + 1, npl
                  dx(:) = xh(j, :) - xh(i, :)
                  rji2  = dx(:) .dot. dx(:) 
                  irij3 = 1.0_DP / (rji2 * sqrt(rji2))
                  faci = Gmpl(i) * irij3
                  facj = Gmpl(j) * irij3
                  ah3(i, :) = ah3(i, :) + facj * dx(:)
                  ah3(j, :) = ah3(j, :) - faci * dx(:)
               end do
            end do
         end associate
      
         return
      end subroutine whm_getacch_ah3

   end procedure whm_getacch_pl

   module procedure whm_getacch_tp !((self, cb, pl, config, t, xh)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from Hal Levison's Swift routine getacch_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_tp.f90
      use swiftest
      implicit none
      integer(I4B)                                 :: i
      real(DP), dimension(:), allocatable, save    :: r2, r2t, fac, mu
      real(DP), dimension(:), allocatable, save    :: irh, ir3h
      real(DP), dimension(:), allocatable, save    :: irht
      real(DP), dimension(:, :), allocatable, save :: aobl
      real(DP), dimension(:, :), allocatable, save :: xht, aoblt
      real(DP), dimension(NDIM)                    :: ah0
   
      associate(tp => self, ntp => self%nbody, npl => pl%nbody, aht => self%ah, &
                status => self%status, xht => self%xh, Gmpl => pl%Gmass, &
                j2rp2 => cb%j2rp2, j4rp4 => cb%j4rp4, aobl => self%aobl, aobl0 => cb%aobl)
         if(.not.allocated(r2)) allocate(r2(npl))
         if(.not.allocated(irh)) allocate(irh(npl))
         if (.not.allocated(ir3h)) allocate(ir3h(npl))
         if (.not. allocated(irht)) allocate(irht(ntp))
         if (.not. allocated(r2t)) allocate(r2t(ntp))
         r2(:) = xh(1:npl, :) .dot. xh(1:npl, :)
         irh(:)= 1.0_DP / sqrt(r2(1:npl))
         ir3h(:) = irh(1:npl) / r2(1:npl)

         do concurrent (i = 1:ntp, status(i) == ACTIVE)
            r2t(i) = xht(i, :) .dot. xht(i, :) 
            irht(i) = 1.0_DP / sqrt(r2t(i))
         end do

         ah0(:) = 0.0_DP
         if(.not.allocated(fac)) allocate(fac(npl))
         fac(:) = Gmpl(1:npl) * ir3h(1:npl)
         do concurrent (i = 1:NDIM)
            ah0(i) = - sum(fac(1:npl) * xh(1:npl, i))
         end do
         call whm_getacch_ah3_tp(cb, pl, tp, xh)
         do concurrent (i = 1:ntp, status(i) == ACTIVE)
            aht(i, :) = aht(i, :) + ah0(:)
         end do
         if (j2rp2 /= 0.0_DP) then
            call pl%obl_acc(cb, irh, xh) 
            call tp%obl_acc(cb, irht, xht)
            do concurrent (i = 1:ntp, status(i) == ACTIVE)
               aht(i, :) = aht(i, :) + tp%aobl(i, :) - aobl0(:)
            end do
         end if
         if (config%lextra_force) call tp%user_getacch(cb, config, t)
         if (config%lgr) call tp%gr_getacch(cb, config) 
      end associate
      return
   
      contains

         pure subroutine whm_getacch_ah3_tp(cb, pl, tp, xh) 
            !! author: David A. Minton
            !!
            !! Compute direct cross (third) term heliocentric accelerations of test particles
            !!
            !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
            !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
            use swiftest
            implicit none
            class(whm_central_body), intent(in) :: cb 
            class(whm_pl), intent(in) :: pl 
            class(whm_tp), intent(inout) :: tp
            real(DP), dimension(:,:), intent(in) :: xh
            integer(I4B)          :: i, j
            real(DP)            :: rji2, irij3, fac
            real(DP), dimension(NDIM) :: dx, acc

            associate(ntp => tp%nbody, npl => pl%nbody, msun => cb%Gmass,  Gmpl => pl%Gmass, &
                      xht => tp%xh, aht => tp%ah)
        
               aht(:,:) = 0.0_DP
               do i = 1, ntp
                  do j = 1, npl
                     dx(:) = xht(i, :) - xh(j, :)
                     rji2 = dx(:) .dot. dx(:) 
                     irij3 = 1.0_DP / (rji2 * sqrt(rji2))
                     fac = Gmpl(j) * irij3
                     aht(i, :) = aht(i, :) - fac * dx(:)
                  end do
               end do
            end associate
            return
         end subroutine whm_getacch_ah3_tp
   end procedure whm_getacch_tp
end submodule s_whm_getacch
