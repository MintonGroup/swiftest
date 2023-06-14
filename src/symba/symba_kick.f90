!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(symba) s_symba_kick
   use swiftest
contains

   module subroutine symba_kick_getacch_int_pl(self, param)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of massive bodies, with no mutual interactions between bodies below GMTINY
      !!
      !! Adapted from Hal Levison's Swift routine symba5_helio_getacch.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_kick_getacch_int.f90
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
      class(swiftest_parameters), intent(inout) :: param !! Current Swiftest run configuration parameter

      if (param%lflatten_interactions) then
         call swiftest_kick_getacch_int_all(self%nbody, self%nplplm, self%k_plpl, self%rh, self%Gmass, self%radius, self%ah)
      else
         call swiftest_kick_getacch_int_all(self%nbody, self%nplm, self%rh, self%Gmass, self%radius, self%ah)
      end if

      return
   end subroutine symba_kick_getacch_int_pl


   module subroutine symba_kick_getacch_pl(self, nbody_system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine symba_kick_getacch.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick_getacch.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I8B)              :: nplplenc
      real(DP), dimension(NDIM,self%nbody) :: ah_enc
      integer(I4B), dimension(:,:), allocatable :: k_plpl_enc

      if (self%nbody == 0) return
      select type(nbody_system)
      class is (symba_nbody_system)
         associate(pl => self, npl => self%nbody, nplm => self%nplm, plpl_encounter => nbody_system%plpl_encounter, radius => self%radius)
            ! Apply kicks to all bodies (including those in the encounter list)
            call helio_kick_getacch_pl(pl, nbody_system, param, t, lbeg)
            if (plpl_encounter%nenc > 0) then 
               ! Remove kicks from bodies involved currently in the encounter list, as these are dealt with separately.
               ah_enc(:,:) = 0.0_DP
               nplplenc = int(plpl_encounter%nenc, kind=I8B)
               allocate(k_plpl_enc(2,nplplenc))
               k_plpl_enc(1,1:nplplenc) = plpl_encounter%index1(1:nplplenc)
               k_plpl_enc(2,1:nplplenc) = plpl_encounter%index2(1:nplplenc)
               call swiftest_kick_getacch_int_all(npl, nplplenc, k_plpl_enc, pl%rh, pl%Gmass, pl%radius, ah_enc)
               pl%ah(:,1:npl) = pl%ah(:,1:npl) - ah_enc(:,1:npl)
            end if

         end associate
      end select

      return
   end subroutine symba_kick_getacch_pl


   module subroutine symba_kick_getacch_tp(self, nbody_system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from David E. Kaufmann's Swifter routine symba_kick_getacch_tp.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick_getacch.f
      implicit none
      ! Arguments
      class(symba_tp),              intent(inout) :: self   !! SyMBA test particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)              :: i, j
      integer(I8B)              :: k
      real(DP)                  :: rjj, fac
      real(DP), dimension(NDIM) :: dx

      if (self%nbody == 0) return
      select type(nbody_system)
      class is (symba_nbody_system)
         associate(tp => self, cb => nbody_system%cb, pl => nbody_system%pl, &
            pltp_encounter => nbody_system%pltp_encounter, npltpenc => nbody_system%pltp_encounter%nenc)
            call helio_kick_getacch_tp(tp, nbody_system, param, t, lbeg)
            ! Remove accelerations from encountering pairs
            do k = 1, npltpenc
               i = pltp_encounter%index1(k)
               j = pltp_encounter%index2(k)
               if (tp%lmask(j)) then
                  if (lbeg) then
                     dx(:) = tp%rh(:,j) - pl%rbeg(:,i)
                  else
                     dx(:) = tp%rh(:,j) - pl%rend(:,i)
                  end if
                  rjj = dot_product(dx(:), dx(:))
                  fac = pl%Gmass(i) / (rjj * sqrt(rjj))
                  tp%ah(:,j) = tp%ah(:,j) + fac * dx(:)
               end if
            end do
         end associate
      end select
      return
   end subroutine symba_kick_getacch_tp


   module subroutine symba_kick_list_plpl(self, nbody_system, dt, irec, sgn)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of massive bodies within SyMBA recursion.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_kick.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick.f
      implicit none
      ! Arguments
      class(symba_list_plpl),    intent(in)    :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system), intent(inout) :: nbody_system !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt     !! step size
      integer(I4B),              intent(in)    :: irec   !! Current recursion level
      integer(I4B),              intent(in)    :: sgn    !! sign to be applied to acceleration
      ! Internals
      integer(I4B)              :: i, j, irm1, irecl, ngood
      integer(I8B)              :: k
      real(DP)                  :: r, rr, ri, ris, rim1, r2, ir3, fac, faci, facj
      real(DP), dimension(NDIM) :: dx
      logical, dimension(:), allocatable :: lgoodlevel
      integer(I8B), dimension(:), allocatable :: good_idx

      if (self%nenc == 0) return

      select type(pl => nbody_system%pl)
      class is (symba_pl)
         associate(npl => pl%nbody, nenc => self%nenc)
            if (npl == 0)  return
            pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE
            allocate(lgoodlevel(nenc))

            irm1 = irec - 1

            if (sgn < 0) then
               irecl = irec - 1
            else
               irecl = irec
            end if

            do k = 1, nenc
               i = self%index1(k)
               j = self%index2(k)
               lgoodlevel(k) = (pl%levelg(i) >= irm1) .and. (pl%levelg(j) >= irm1)
               lgoodlevel(k) = (self%status(k) == ACTIVE) .and. lgoodlevel(k)
            end do
            ngood = count(lgoodlevel(:))
            if (ngood > 0_I8B) then
               allocate(good_idx(ngood))
               good_idx(:) = pack([(k, k = 1_I8B, nenc)], lgoodlevel(:))

#ifdef DOCONLOC
               do concurrent (k = 1:ngood) shared(self,pl,good_idx) local(i,j)
#else
               do concurrent (k = 1:ngood)
#endif
                  i = self%index1(good_idx(k))
                  j = self%index2(good_idx(k))
                  pl%ah(:,i) = 0.0_DP
                  pl%ah(:,j) = 0.0_DP
               end do

               do k = 1, ngood
                  i = self%index1(good_idx(k))
                  j = self%index2(good_idx(k))
                  ri = ((pl%rhill(i)  + pl%rhill(j))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                  rim1 = ri * (RSHELL**2)
                  dx(:) = pl%rh(:,j) - pl%rh(:,i)

                  r2 = dot_product(dx(:), dx(:))
                  if (r2 < rim1) then
                     fac = 0.0_DP
                     lgoodlevel(good_idx(k)) = .false.
                     cycle
                  end if
                  if (r2 < ri) then
                     ris = sqrt(ri)
                     r = sqrt(r2)
                     rr = (ris - r) / (ris * (1.0_DP - RSHELL))
                     fac = (r2**(-1.5_DP)) * (1.0_DP - 3 * (rr**2) + 2 * (rr**3))
                  else
                     ir3 = 1.0_DP / (r2 * sqrt(r2))
                     fac = ir3
                  end if
                  faci = fac * pl%Gmass(i)
                  facj = fac * pl%Gmass(j)
                  pl%ah(:, i) = pl%ah(:, i) + facj * dx(:)
                  pl%ah(:, j) = pl%ah(:, j) - faci * dx(:)
               end do
               ngood = count(lgoodlevel(:))
               if (ngood == 0_I8B) return
               good_idx(1:ngood) = pack([(k, k = 1_I8B, nenc)], lgoodlevel(:))

               do k = 1, ngood
                  i = self%index1(good_idx(k))
                  j = self%index2(good_idx(k))
                  pl%vb(:,i) = pl%vb(:,i) + sgn * dt * pl%ah(:,i)
                  pl%vb(:,j) = pl%vb(:,j) + sgn * dt * pl%ah(:,j)
                  pl%ah(:,i) = 0.0_DP
                  pl%ah(:,j) = 0.0_DP
               end do

            end if
         end associate
      end select
      
      return
   end subroutine symba_kick_list_plpl


   module subroutine symba_kick_list_pltp(self, nbody_system, dt, irec, sgn)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of ACTIVE test particles within SyMBA recursion.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_kick.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick.f
      implicit none
      ! Arguments
      class(symba_list_pltp),    intent(in)    :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system), intent(inout) :: nbody_system !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt     !! step size
      integer(I4B),              intent(in)    :: irec   !! Current recursion level
      integer(I4B),              intent(in)    :: sgn    !! sign to be applied to acceleration
      ! Internals
      integer(I4B)              :: i, j, irm1, irecl, ngood
      integer(I8B)              :: k
      real(DP)                  :: r, rr, ri, ris, rim1, r2, ir3, fac, faci
      real(DP), dimension(NDIM) :: dx
      logical, dimension(:), allocatable :: lgoodlevel
      integer(I8B), dimension(:), allocatable :: good_idx

      if (self%nenc == 0) return

      select type(pl => nbody_system%pl)
      class is (symba_pl)
      select type(tp => nbody_system%tp)
      class is (symba_tp)
         associate(npl => pl%nbody, ntp => tp%nbody, nenc => self%nenc)
            if ((npl == 0) .or. (ntp == 0)) return
            pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE
            tp%lmask(1:ntp) = tp%status(1:ntp) /= INACTIVE
            allocate(lgoodlevel(nenc))

            irm1 = irec - 1

            if (sgn < 0) then
               irecl = irec - 1
            else
               irecl = irec
            end if

            do k = 1, nenc
               i = self%index1(k)
               j = self%index2(k)
               lgoodlevel(k) = (pl%levelg(i) >= irm1) .and. (tp%levelg(j) >= irm1)
               lgoodlevel(k) = (self%status(k) == ACTIVE) .and. lgoodlevel(k)
            end do

            ngood = count(lgoodlevel(:))

            if (ngood > 0_I8B) then
               allocate(good_idx(ngood))
               good_idx(:) = pack([(k, k = 1_I8B, nenc)], lgoodlevel(:))

#ifdef DOCONLOC
               do concurrent (k = 1_I8B:ngood) shared(self,tp,good_idx) local(j)
#else
               do concurrent (k = 1_I8B:ngood)
#endif
                  j = self%index2(good_idx(k))
                  tp%ah(:,j) = 0.0_DP
               end do

               do k = 1, ngood
                  i = self%index1(good_idx(k))
                  j = self%index2(good_idx(k))

                  ri = ((pl%rhill(i))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                  rim1 = ri * (RSHELL**2)
                  dx(:) = tp%rh(:,j) - pl%rh(:,i)
                  r2 = dot_product(dx(:), dx(:))
                  if (r2 < rim1) then
                     fac = 0.0_DP
                     lgoodlevel(good_idx(k)) = .false.
                     cycle
                  end if
                  if (r2 < ri) then
                     ris = sqrt(ri)
                     r = sqrt(r2)
                     rr = (ris - r) / (ris * (1.0_DP - RSHELL))
                     fac = (r2**(-1.5_DP)) * (1.0_DP - 3 * (rr**2) + 2 * (rr**3))
                  else
                     ir3 = 1.0_DP / (r2 * sqrt(r2))
                     fac = ir3
                  end if
                  faci = fac * pl%Gmass(i)

                  tp%ah(:, j) = tp%ah(:, j) - faci * dx(:)
               end do
               ngood = count(lgoodlevel(:))
               if (ngood == 0_I8B) return
               good_idx(1:ngood) = pack([(k, k = 1_I8B, nenc)], lgoodlevel(:))

               do k = 1, ngood
                  j = self%index2(good_idx(k))
                  tp%vb(:,j) = tp%vb(:,j) + sgn * dt * tp%ah(:,j)
                  tp%ah(:,j) = 0.0_DP
               end do
            end if
         end associate
      end select
      end select
      
      return
   end subroutine symba_kick_list_pltp


end submodule s_symba_kick