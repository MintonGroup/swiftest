!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_swiftest_obl
contains
   module subroutine swiftest_obl_acc(n, GMcb, j2rp2, j4rp4, rh, lmask, aobl, GMpl, aoblcb)
      !! author: David A. Minton
      !!
      !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      !!      Returned values do not include monopole term or terms higher than J4
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_acc.f90 and obl_acc_tp.f90
      !! Adapted from Hal Levison's Swift routine obl_acc.f and obl_acc_tp.f
      implicit none
      ! Arguments
      integer(I4B),             intent(in)            :: n      !! Number of bodies
      real(DP),                 intent(in)            :: GMcb   !! Central body G*Mass
      real(DP),                 intent(in)            :: j2rp2  !! J2 * R**2 for the central body
      real(DP),                 intent(in)            :: j4rp4  !! J4 * R**4 for the central body
      real(DP), dimension(:,:), intent(in)            :: rh     !! Heliocentric positions of bodies
      logical,  dimension(:),   intent(in)            :: lmask  !! Logical mask of bodies to compute aobl
      real(DP), dimension(:,:), intent(out)           :: aobl   !! Barycentric acceleration of bodies due to central body oblateness
      real(DP), dimension(:),   intent(in),  optional :: GMpl   !! Masses of input bodies if they are not test particles
      real(DP), dimension(:),   intent(out), optional :: aoblcb !! Barycentric acceleration of central body (only needed if input bodies are massive)
      ! Internals
      integer(I4B) :: i
      real(DP)     :: r2, irh, rinv2, t0, t1, t2, t3, fac1, fac2

      if (n == 0) return

      aobl(:,:) = 0.0_DP
#ifdef DOCONLOC
      do concurrent(i = 1:n, lmask(i)) shared(lmask,rh,aobl) local(r2,irh,rinv2,t0,t1,t2,t3,fac1,fac2)
#else
      do concurrent(i = 1:n, lmask(i))
#endif
         r2 = dot_product(rh(:, i), rh(:, i))
         irh = 1.0_DP / sqrt(r2)
         rinv2 = irh**2
         t0 = -GMcb * rinv2 * rinv2 * irh
         t1 = 1.5_DP * j2rp2
         t2 = rh(3, i) * rh(3, i) * rinv2
         t3 = 1.875_DP * j4rp4 * rinv2
         fac1 = t0 * (t1 - t3 - (5 * t1 - (14.0_DP - 21.0_DP * t2) * t3) * t2)
         fac2 = 2 * t0 * (t1 - (2.0_DP - (14.0_DP * t2 / 3.0_DP)) * t3)
         aobl(:, i) = fac1 * rh(:, i)
         aobl(3, i) = fac2 * rh(3, i) + aobl(3, i)
      end do

      if (present(GMpl) .and. present(aoblcb)) then
         aoblcb(:) = 0.0_DP
         do i = n, 1, -1
            if (lmask(i)) aoblcb(:) = aoblcb(:) - GMpl(i) * aobl(:, i) / GMcb
         end do
      end if

      return

   end subroutine swiftest_obl_acc


   module subroutine swiftest_obl_acc_pl(self, nbody_system)
      !! author: David A. Minton
      !!
      !! Compute the barycentric accelerations of massive bodies due to the oblateness of the central body
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_acc.f90 and obl_acc_tp.f90
      !! Adapted from Hal Levison's Swift routine obl_acc.f and obl_acc_tp.f
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      ! Internals
      integer(I4B) :: i, npl

      if (self%nbody == 0) return

      associate(pl => self, cb => nbody_system%cb)
         npl = self%nbody
         call swiftest_obl_acc(npl, cb%Gmass, cb%j2rp2, cb%j4rp4, pl%rh, pl%lmask, pl%aobl, pl%Gmass, cb%aobl)

#ifdef DOCONLOC
         do concurrent(i = 1:npl, pl%lmask(i)) shared(cb,pl)
#else
         do concurrent(i = 1:npl, pl%lmask(i))
#endif
            pl%ah(:, i) = pl%ah(:, i) + pl%aobl(:, i) - cb%aobl(:)
         end do
      end associate

      return

   end subroutine swiftest_obl_acc_pl


   module subroutine swiftest_obl_acc_tp(self, nbody_system)
      !! author: David A. Minton
      !!
      !! Compute the barycentric accelerations of massive bodies due to the oblateness of the central body
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_acc.f90 and obl_acc_tp.f90
      !! Adapted from Hal Levison's Swift routine obl_acc.f and obl_acc_tp.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      ! Internals
      real(DP), dimension(NDIM)                   :: aoblcb
      integer(I4B) :: i, ntp

      if (self%nbody == 0) return

      associate(tp => self, cb => nbody_system%cb)
         ntp = self%nbody
         call swiftest_obl_acc(ntp, cb%Gmass, cb%j2rp2, cb%j4rp4, tp%rh, tp%lmask, tp%aobl)
         if (nbody_system%lbeg) then
            aoblcb = cb%aoblbeg
         else
            aoblcb = cb%aoblend
         end if

#ifdef DOCONLOC
         do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp,aoblcb)
#else
         do concurrent(i = 1:ntp, tp%lmask(i))
#endif
            tp%ah(:, i) = tp%ah(:, i) + tp%aobl(:, i) - aoblcb(:)
         end do

      end associate
      return

   end subroutine swiftest_obl_acc_tp


   module subroutine swiftest_obl_pot_system(self) 
      !! author: David A. Minton
      !!
      !! Compute the contribution to the total gravitational potential due solely to the oblateness of the central body
      !!    Returned value does not include monopole term or terms higher than J4
      !!
      !!    Reference: MacMillan, W. D. 1958. The Theory of the Potential, (Dover Publications), 363.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_pot.f90 
      !! Adapted from Hal Levison's Swift routine obl_pot.f 
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout)  :: self   !! Swiftest nbody system object
      ! Internals
      integer(I4B) :: i, npl
      real(DP), dimension(self%pl%nbody) :: oblpot_arr

      associate(nbody_system => self, pl => self%pl, cb => self%cb)
         npl = self%pl%nbody
         if (.not. any(pl%lmask(1:npl))) return
#ifdef DOCONLOC
         do concurrent (i = 1:npl, pl%lmask(i)) shared(cb,pl,oblpot_arr)
#else
         do concurrent (i = 1:npl, pl%lmask(i))
#endif
            oblpot_arr(i) = swiftest_obl_pot_one(cb%Gmass, pl%Gmass(i), cb%j2rp2, cb%j4rp4, pl%rh(3,i), 1.0_DP / norm2(pl%rh(:,i)))
         end do
         nbody_system%oblpot = sum(oblpot_arr, pl%lmask(1:npl))
      end associate
         
      return
   end subroutine swiftest_obl_pot_system


   elemental function swiftest_obl_pot_one(GMcb, GMpl, j2rp2, j4rp4, zh, irh) result(oblpot)
      !! author: David A. Minton
      !!
      !! Compute the contribution to the total gravitational potential due solely to the oblateness of the central body from a single massive body
      !!    Returned value does not include monopole term or terms higher than J4
      !!
      !!    Reference: MacMillan, W. D. 1958. The Theory of the Potential, (Dover Publications), 363.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_pot.f90 
      !! Adapted from Hal Levison's Swift routine obl_pot.f 
      implicit none
      ! Arguments
      real(DP),     intent(in)  :: GMcb   !! G*mass of the central body
      real(DP),     intent(in)  :: GMpl   !! G*mass of the massive body
      real(DP),     intent(in)  :: j2rp2  !! J_2 / R**2 of the central body
      real(DP),     intent(in)  :: j4rp4  !! J_2 / R**4 of the central body
      real(DP),     intent(in)  :: zh     !! z-component of the heliocentric distance vector of the massive body
      real(DP),     intent(in)  :: irh    !! Inverse of the heliocentric distance magnitude of the massive body
      ! Result
      real(DP)                  :: oblpot !! Gravitational potential
         
      ! Internals
      real(DP)                  :: rinv2, t0, t1, t2, t3, p2, p4
         
      rinv2 = irh**2
      t0 = GMcb * GMpl * rinv2 * irh
      t1 = j2rp2
      t2 = zh**2 * rinv2
      t3 = j4rp4 * rinv2
      p2 = 0.5_DP * (3 * t2 - 1.0_DP)
      p4 = 0.125_DP * ((35 * t2 - 30.0_DP) * t2 + 3.0_DP)
      oblpot = t0 * (t1 * p2 + t3 * p4)
         
      return
   end function swiftest_obl_pot_one

end submodule s_swiftest_obl
