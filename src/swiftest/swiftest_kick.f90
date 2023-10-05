!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(swiftest) s_swiftest_kick
contains
   module subroutine swiftest_kick_getacch_int_pl(self, param)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of massive bodies
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f90
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
      class(swiftest_parameters), intent(inout) :: param !! Current swiftest run configuration parameters
      ! Internals
#ifdef PROFILE
      type(walltimer), save :: timer 
#endif

      if (param%lflatten_interactions) then
         if (param%lclose) then
            call swiftest_kick_getacch_int_all(self%nbody, self%nplpl, self%k_plpl, self%rh, self%Gmass, self%radius, self%ah)
         else
            call swiftest_kick_getacch_int_all(self%nbody, self%nplpl, self%k_plpl, self%rh, self%Gmass, self%ah)
         end if
      else
         if (param%lclose) then
            call swiftest_kick_getacch_int_all(self%nbody, self%nbody, self%rh, self%Gmass, self%radius, self%ah)
         else
            call swiftest_kick_getacch_int_all(self%nbody, self%nbody, self%rh, self%Gmass, self%ah)
         end if
      end if

      return
   end subroutine swiftest_kick_getacch_int_pl


   module subroutine swiftest_kick_getacch_int_tp(self, param, GMpl, rhp, npl)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of test particles by massive bodies
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int_tp.f90
      implicit none
      ! Arguments
      class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
      class(swiftest_parameters), intent(inout) :: param !! Current swiftest run configuration parameters
      real(DP), dimension(:),     intent(in)    :: GMpl  !! Massive body masses
      real(DP), dimension(:,:),   intent(in)    :: rhp   !! Massive body position vectors
      integer(I4B),               intent(in)    :: npl   !! Number of active massive bodies

      if ((self%nbody == 0) .or. (npl == 0)) return

      call swiftest_kick_getacch_int_all_tp(self%nbody, npl, self%rh, rhp, GMpl, self%lmask, self%ah)
      
      return
   end subroutine swiftest_kick_getacch_int_tp


   module subroutine swiftest_kick_getacch_int_all_flat_rad_pl(npl, nplpl, k_plpl, r, Gmass, radius, acc)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations for massive bodies, with parallelization.
      !! This is the flattened (single loop) version that uses the k_plpl interaction pair index array
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f9
      implicit none
      integer(I4B),                 intent(in)             :: npl    !! Number of massive bodies
      integer(I8B),                 intent(in)             :: nplpl  !! Number of massive body interactions to compute
      integer(I4B), dimension(:,:), intent(in)             :: k_plpl !! Array of interaction pair indices (flattened upper triangular matrix)
      real(DP),     dimension(:,:), intent(in)             :: r      !! Position vector array
      real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
      real(DP),     dimension(:),   intent(in)             :: radius !! Array of massive body radii
      real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      ! Internals
      integer(I8B)                      :: k
      real(DP), dimension(NDIM,npl) :: ahi, ahj
      integer(I4B) :: i, j
      real(DP)     :: rji2, rlim2
      real(DP)     :: rx, ry, rz

      ahi(:,:) = 0.0_DP
      ahj(:,:) = 0.0_DP

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(nplpl, k_plpl, r, Gmass, radius) &
      !$omp lastprivate(i, j, rji2, rlim2, rx, ry, rz) &
      !$omp reduction(+:ahi,ahj) 
      do k = 1_I8B, nplpl
         i = k_plpl(1, k)
         j = k_plpl(2, k)
         rx = r(1, j) - r(1, i) 
         ry = r(2, j) - r(2, i) 
         rz = r(3, j) - r(3, i) 
         rji2 = rx**2 + ry**2 + rz**2
         rlim2 = (radius(i) + radius(j))**2
         if (rji2 > rlim2) call swiftest_kick_getacch_int_one_pl(rji2, rx, ry, rz, Gmass(i), Gmass(j), &
                                 ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
      end do
      !$omp end parallel do 

      acc(:,:) = acc(:,:) + ahi(:,:) + ahj(:,:)

      return
   end subroutine swiftest_kick_getacch_int_all_flat_rad_pl


   module subroutine swiftest_kick_getacch_int_all_flat_norad_pl(npl, nplpl, k_plpl, r, Gmass, acc)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations for massive bodies, with parallelization.
      !! This is the flattened (single loop) version that uses the k_plpl interaction pair index array
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f9
      implicit none
      integer(I4B),                 intent(in)             :: npl    !! Number of massive bodies
      integer(I8B),                 intent(in)             :: nplpl  !! Number of massive body interactions to compute
      integer(I4B), dimension(:,:), intent(in)             :: k_plpl !! Array of interaction pair indices (flattened upper triangular matrix)
      real(DP),     dimension(:,:), intent(in)             :: r      !! Position vector array
      real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
      real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      ! Internals
      integer(I8B)                      :: k
      real(DP), dimension(NDIM,npl) :: ahi, ahj
      integer(I4B) :: i, j
      real(DP)     :: rji2
      real(DP)     :: rx, ry, rz

      ahi(:,:) = 0.0_DP
      ahj(:,:) = 0.0_DP

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(nplpl, k_plpl, r, Gmass) &
      !$omp lastprivate(i, j, rji2, rx, ry, rz) &
      !$omp reduction(+:ahi,ahj) 
      do k = 1_I8B, nplpl
         i = k_plpl(1, k)
         j = k_plpl(2, k)
         rx = r(1, j) - r(1, i) 
         ry = r(2, j) - r(2, i) 
         rz = r(3, j) - r(3, i) 
         rji2 = rx**2 + ry**2 + rz**2
         call swiftest_kick_getacch_int_one_pl(rji2, rx, ry, rz, Gmass(i), Gmass(j), &
                                       ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
      end do
      !$omp end parallel do
     
      acc(:,:) = acc(:,:) + ahi(:,:) + ahj(:,:)

      return
   end subroutine swiftest_kick_getacch_int_all_flat_norad_pl


   module subroutine swiftest_kick_getacch_int_all_tri_rad_pl(npl, nplm, r, Gmass, radius, acc)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations for massive bodies, with parallelization.
      !! This is the upper triangular matrix (double loop) version.
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f9
      implicit none
      integer(I4B),                 intent(in)             :: npl    !! Total number of massive bodies
      integer(I4B),                 intent(in)             :: nplm   !! Number of fully interacting massive bodies
      real(DP),     dimension(:,:), intent(in)             :: r      !! Position vector array
      real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
      real(DP),     dimension(:),   intent(in)             :: radius !! Array of massive body radii
      real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      ! Internals
      integer(I4B) :: i, j, nplt
      real(DP)     :: rji2, rlim2, fac, rx, ry, rz
      real(DP), dimension(NDIM,npl) :: ahi, ahj
      logical :: lmtiny

      nplt = npl - nplm
      lmtiny = (nplt > nplm)

      if (lmtiny) then
         ahi(:,:) = 0.0_DP
         ahj(:,:) = 0.0_DP
         !$omp parallel do default(private) schedule(static)&
         !$omp shared(npl, nplm, r, Gmass, radius) &
         !$omp reduction(+:ahi,ahj)
         do i = 1, nplm
#ifdef DOCONLOC
            do concurrent(j = i+1:npl) shared(i,r,radius,ahi,ahj,Gmass) local(rx,ry,rz,rji2,rlim2)
#else
            do concurrent(j = i+1:npl)
#endif
               rx = r(1, j) - r(1, i) 
               ry = r(2, j) - r(2, i) 
               rz = r(3, j) - r(3, i) 
               rji2 = rx**2 + ry**2 + rz**2
               rlim2 = (radius(i) + radius(j))**2
               if (rji2 > rlim2) call swiftest_kick_getacch_int_one_pl(rji2, rx, ry, rz, Gmass(i), Gmass(j), &
                                          ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
            end do
         end do
         !$omp end parallel do
#ifdef DOCONLOC
         do concurrent(i = 1:npl) shared(acc,ahi,ahj)
#else
         do concurrent(i = 1:npl)
#endif
            acc(:,i) = acc(:,i) + ahi(:,i) + ahj(:,i)
         end do
      else 
         !$omp parallel do default(private) schedule(static)&
         !$omp shared(npl, nplm, r, Gmass, radius, acc)
         do i = 1, nplm
#ifdef DOCONLOC
            do concurrent(j = 1:npl, i/=j) shared(i,r,radius,Gmass,acc) local(rx,ry,rz,rji2,rlim2,fac)
#else
            do concurrent(j = 1:npl, i/=j)
#endif
               rx = r(1,j) - r(1,i)
               ry = r(2,j) - r(2,i)
               rz = r(3,j) - r(3,i)
               rji2 = rx**2 + ry**2 + rz**2
               rlim2 = (radius(i) + radius(j))**2
               if (rji2 > rlim2)  then
                  fac = Gmass(j) / (rji2 * sqrt(rji2))
                  acc(1,i) = acc(1,i) + fac * rx
                  acc(2,i) = acc(2,i) + fac * ry
                  acc(3,i) = acc(3,i) + fac * rz
               end if
            end do
         end do
         !$omp end parallel do

         if (nplt > 0) then
            !$omp parallel do default(private) schedule(static)&
            !$omp shared(npl, nplm, r, Gmass, radius, acc)
            do i = nplm+1,npl
#ifdef DOCONLOC
               do concurrent(j = 1:nplm) shared(i,r,radius,Gmass,acc) local(rx,ry,rz,rji2,rlim2,fac)
#else
               do concurrent(j = 1:nplm)
#endif
                  rx = r(1,j) - r(1,i)
                  ry = r(2,j) - r(2,i)
                  rz = r(3,j) - r(3,i)
                  rji2 = rx**2 + ry**2 + rz**2
                  rlim2 = (radius(i) + radius(j))**2
                  if (rji2 > rlim2)  then
                     fac = Gmass(j) / (rji2 * sqrt(rji2))
                     acc(1,i) = acc(1,i) + fac * rx
                     acc(2,i) = acc(2,i) + fac * ry
                     acc(3,i) = acc(3,i) + fac * rz
                  end if
               end do
            end do
            !$omp end parallel do
         end if

      end if


      return
   end subroutine swiftest_kick_getacch_int_all_tri_rad_pl


   module subroutine swiftest_kick_getacch_int_all_tri_norad_pl(npl, nplm, r, Gmass, acc)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations for massive bodies, with parallelization.
      !! This is the upper triangular matrix (double loop) version.
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f9
      implicit none
      integer(I4B),                 intent(in)             :: npl    !! Total number of massive bodies
      integer(I4B),                 intent(in)             :: nplm   !! Number of fully interacting massive bodies
      real(DP),     dimension(:,:), intent(in)             :: r      !! Position vector array
      real(DP),     dimension(:),   intent(in)             :: Gmass  !! Array of massive body G*mass
      real(DP),     dimension(:,:), intent(inout)          :: acc    !! Acceleration vector array 
      ! Internals
      integer(I4B) :: i, j, nplt
      real(DP)     :: rji2, fac, rx, ry, rz
      real(DP), dimension(NDIM,npl) :: ahi, ahj
      logical :: lmtiny

      nplt = npl - nplm
      lmtiny = (nplt > nplm)

      if (lmtiny) then
         ahi(:,:) = 0.0_DP
         ahj(:,:) = 0.0_DP
         !$omp parallel do default(private) schedule(static)&
         !$omp shared(npl, nplm, r, Gmass) &
         !$omp reduction(+:ahi,ahj)
         do i = 1, nplm
#ifdef DOCONLOC
            do concurrent(j = i+1:npl) shared(i,r,Gmass,ahi,ahj) local(rx,ry,rz,rji2)
#else
            do concurrent(j = i+1:npl)
#endif
               rx = r(1, j) - r(1, i) 
               ry = r(2, j) - r(2, i) 
               rz = r(3, j) - r(3, i) 
               rji2 = rx**2 + ry**2 + rz**2
               call swiftest_kick_getacch_int_one_pl(rji2, rx, ry, rz, Gmass(i), Gmass(j), &
                                          ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
            end do
         end do
         !$omp end parallel do
#ifdef DOCONLOC
         do concurrent(i = 1:npl) shared(acc,ahi,ahj)
#else
         do concurrent(i = 1:npl)
#endif
            acc(:,i) = acc(:,i) + ahi(:,i) + ahj(:,i)
         end do
      else 
         !$omp parallel do default(private) schedule(static)&
         !$omp shared(npl, nplm, r, Gmass, acc)
         do i = 1, nplm
#ifdef DOCONLOC
            do concurrent(j = 1:npl, j/=i) shared(i,r,Gmass, acc) local(rx,ry,rz,rji2,fac)
#else
            do concurrent(j = 1:npl, j/=i)
#endif
               rx = r(1,j) - r(1,i)
               ry = r(2,j) - r(2,i)
               rz = r(3,j) - r(3,i)
               rji2 = rx**2 + ry**2 + rz**2
               fac = Gmass(j) / (rji2 * sqrt(rji2))
               acc(1,i) = acc(1,i) + fac * rx
               acc(2,i) = acc(2,i) + fac * ry
               acc(3,i) = acc(3,i) + fac * rz
            end do
         end do
         !$omp end parallel do

         if (nplt > 0) then
            !$omp parallel do default(private) schedule(static)&
            !$omp shared(npl, nplm, r, Gmass, acc)
            do i = nplm+1,npl
#ifdef DOCONLOC
               do concurrent(j = 1:nplm) shared(i,r,Gmass,acc) local(rx,ry,rz,rji2,fac)
#else
               do concurrent(j = 1:nplm)
#endif
                  rx = r(1,j) - r(1,i)
                  ry = r(2,j) - r(2,i)
                  rz = r(3,j) - r(3,i)
                  rji2 = rx**2 + ry**2 + rz**2
                  fac = Gmass(j) / (rji2 * sqrt(rji2))
                  acc(1,i) = acc(1,i) + fac * rx
                  acc(2,i) = acc(2,i) + fac * ry
                  acc(3,i) = acc(3,i) + fac * rz
               end do
            end do
            !$omp end parallel do
         end if

      end if

      return
   end subroutine swiftest_kick_getacch_int_all_tri_norad_pl


   module subroutine swiftest_kick_getacch_int_all_tp(ntp, npl, rtp, rpl, GMpl, lmask, acc)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of test particles by massive bodies with parallelisim
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int_tp.f99
      implicit none
      integer(I4B),                 intent(in)    :: ntp    !! Number of test particles
      integer(I4B),                 intent(in)    :: npl    !! Number of massive bodies
      real(DP),     dimension(:,:), intent(in)    :: rtp    !! Test particle position vector array
      real(DP),     dimension(:,:), intent(in)    :: rpl    !! Massive body particle position vector array
      real(DP),     dimension(:),   intent(in)    :: GMpl   !! Array of massive body G*mass
      logical,      dimension(:),   intent(in)    :: lmask  !! Logical mask indicating which test particles should be computed
      real(DP),     dimension(:,:), intent(inout) :: acc    !! Acceleration vector array 
      ! Internals
      real(DP)     :: rji2
      real(DP)     :: rx, ry, rz
      integer(I4B) :: i, j

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(npl, ntp, lmask, rtp, rpl, GMpl) &
      !$omp reduction(+:acc)
      do i = 1, ntp
         if (lmask(i)) then
#ifdef DOCONLOC
            do concurrent (j = 1:npl) shared(rtp,rpl,GMpl,acc) local(rx,ry,rz,rji2)
#else
            do concurrent (j = 1:npl)
#endif
               rx = rtp(1, i) - rpl(1, j)
               ry = rtp(2, i) - rpl(2, j)
               rz = rtp(3, i) - rpl(3, j)
               rji2 = rx**2 + ry**2 + rz**2
               call swiftest_kick_getacch_int_one_tp(rji2, rx, ry, rz, GMpl(j), acc(1,i), acc(2,i), acc(3,i))
            end do
         end if
      end do
      !$omp end parallel do
      
      return
   end subroutine swiftest_kick_getacch_int_all_tp


   pure module subroutine swiftest_kick_getacch_int_one_pl(rji2, xr, yr, zr, Gmi, Gmj, axi, ayi, azi, axj, ayj, azj)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations for a single pair of massive bodies
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int.f9
      implicit none
      real(DP), intent(in)  :: rji2            !! Square of distance between the two bodies
      real(DP), intent(in)  :: xr, yr, zr      !! Distances between the two bodies in x, y, and z directions
      real(DP), intent(in)  :: Gmi             !! G*mass of body i
      real(DP), intent(in)  :: Gmj             !! G*mass of body j
      real(DP), intent(inout) :: axi, ayi, azi !! Acceleration vector components of body i
      real(DP), intent(inout) :: axj, ayj, azj !! Acceleration vector components of body j
      ! Internals
      real(DP) :: faci, facj, irij3

      irij3 = 1.0_DP / (rji2 * sqrt(rji2))
      faci = Gmi * irij3
      facj = Gmj * irij3
      axi = axi + facj * xr
      ayi = ayi + facj * yr
      azi = azi + facj * zr
      axj = axj - faci * xr
      ayj = ayj - faci * yr
      azj = azj - faci * zr

      return
   end subroutine swiftest_kick_getacch_int_one_pl


   pure module subroutine swiftest_kick_getacch_int_one_tp(rji2, xr, yr, zr, GMpl, ax, ay, az)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of a single test particle massive body pair.
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah3.f90 and helio_kick_getacch_int_tp.f90
      implicit none
      real(DP), intent(in)  :: rji2         !! Square of distance between the test particle and massive body
      real(DP), intent(in)  :: xr, yr, zr   !! Distances between the two bodies in x, y, and z directions
      real(DP), intent(in)  :: Gmpl         !! G*mass of massive body
      real(DP), intent(inout) :: ax, ay, az !! Acceleration vector components of test particle
      ! Internals
      real(DP) :: fac

      fac = GMpl / (rji2*sqrt(rji2))
      ax = ax - fac * xr
      ay = ay - fac * yr
      az = az - fac * zr

      return
   end subroutine swiftest_kick_getacch_int_one_tp

end submodule s_swiftest_kick
