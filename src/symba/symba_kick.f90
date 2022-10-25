submodule(symba_classes) s_symba_kick
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
      class(swiftest_parameters), intent(inout) :: param !! Current swiftest run configuration parameter
      ! Internals
      type(interaction_timer), save :: itimer
      logical, save :: lfirst = .true.

      if (param%ladaptive_interactions) then
         if (self%nplplm > 0) then
            if (lfirst) then
               write(itimer%loopname, *)  "symba_kick_getacch_int_pl"
               write(itimer%looptype, *)  "INTERACTION"
               call itimer%time_this_loop(param, self%nplplm, self)
               lfirst = .false.
            else
               if (itimer%check(param, self%nplplm)) call itimer%time_this_loop(param, self%nplplm, self)
            end if
         else
            param%lflatten_interactions = .false.
         end if
      end if

      if (param%lflatten_interactions) then
         call kick_getacch_int_all_flat_pl(self%nbody, self%nplplm, self%k_plpl, self%xh, self%Gmass, self%radius, self%ah)
      else
         call kick_getacch_int_all_triangular_pl(self%nbody, self%nplm, self%xh, self%Gmass, self%radius, self%ah)
      end if

      if (param%ladaptive_interactions .and. self%nplplm > 0) then 
         if (itimer%is_on) call itimer%adapt(param, self%nplplm, self)
      end if

      return
   end subroutine symba_kick_getacch_int_pl


   module subroutine symba_kick_getacch_pl(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine symba_kick_getacch.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick_getacch.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I8B)              :: nplplenc
      real(DP), dimension(NDIM,self%nbody) :: ah_enc
      integer(I4B), dimension(:,:), allocatable :: k_plpl_enc

      if (self%nbody == 0) return
      select type(system)
      class is (symba_nbody_system)
         associate(pl => self, npl => self%nbody, nplm => self%nplm, plplenc_list => system%plplenc_list, radius => self%radius)
            ! Apply kicks to all bodies (including those in the encounter list)
            call helio_kick_getacch_pl(pl, system, param, t, lbeg)
            if (plplenc_list%nenc > 0) then 
               ! Remove kicks from bodies involved currently in the encounter list, as these are dealt with separately.
               ah_enc(:,:) = 0.0_DP
               nplplenc = int(plplenc_list%nenc, kind=I8B)
               allocate(k_plpl_enc(2,nplplenc))
               k_plpl_enc(1,1:nplplenc) = plplenc_list%index1(1:nplplenc)
               k_plpl_enc(2,1:nplplenc) = plplenc_list%index2(1:nplplenc)
               call kick_getacch_int_all_flat_pl(npl, nplplenc, k_plpl_enc, pl%xh, pl%Gmass, pl%radius, ah_enc)
               pl%ah(:,1:npl) = pl%ah(:,1:npl) - ah_enc(:,1:npl)
            end if

         end associate
      end select

      return
   end subroutine symba_kick_getacch_pl


   module subroutine symba_kick_getacch_tp(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from David E. Kaufmann's Swifter routine symba_kick_getacch_tp.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick_getacch.f
      implicit none
      ! Arguments
      class(symba_tp),              intent(inout) :: self   !! SyMBA test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)              :: i, j, k
      real(DP)                  :: rjj, fac
      real(DP), dimension(NDIM) :: dx

      if (self%nbody == 0) return
      select type(system)
      class is (symba_nbody_system)
         associate(tp => self, cb => system%cb, pl => system%pl, &
            pltpenc_list => system%pltpenc_list, npltpenc => system%pltpenc_list%nenc)
            call helio_kick_getacch_tp(tp, system, param, t, lbeg)
            ! Remove accelerations from encountering pairs
            do k = 1, npltpenc
               i = pltpenc_list%index1(k)
               j = pltpenc_list%index2(k)
               if (tp%lmask(j)) then
                  if (lbeg) then
                     dx(:) = tp%xh(:,j) - pl%xbeg(:,i)
                  else
                     dx(:) = tp%xh(:,j) - pl%xend(:,i)
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


   module subroutine symba_kick_encounter(self, system, dt, irec, sgn)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of massive bodies and ACTIVE test particles within SyMBA recursion.
      !! Note: This method works for the polymorphic symba_pltpenc and symba_plplenc types
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_kick.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick.f
      implicit none
      ! Arguments
      class(symba_encounter),    intent(in)    :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt     !! step size
      integer(I4B),              intent(in)    :: irec   !! Current recursion level
      integer(I4B),              intent(in)    :: sgn    !! sign to be applied to acceleration
      ! Internals
      integer(I4B)              :: i, j, irm1, irecl, ngood
      integer(I8B)              :: k
      real(DP)                  :: r, rr, ri, ris, rim1, r2, ir3, fac, faci, facj
      real(DP), dimension(NDIM) :: dx
      logical                   :: isplpl
      logical, dimension(:), allocatable :: lgoodlevel
      integer(I4B), dimension(:), allocatable :: good_idx

      if (self%nenc == 0) return

      select type(self)
      class is (symba_plplenc)
         isplpl = .true.
      class is (symba_pltpenc)
         isplpl = .false.
      end select
      select type(pl => system%pl)
      class is (symba_pl)
         select type(tp => system%tp)
         class is (symba_tp)
            associate(npl => pl%nbody, ntp => tp%nbody, nenc => self%nenc)
               if (npl == 0)  return
               pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE
               if (.not. isplpl) then
                  if (ntp == 0) return
                  tp%lmask(1:ntp) = tp%status(1:ntp) /= INACTIVE
               end if
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
                  if (isplpl) then
                     lgoodlevel(k) = (pl%levelg(i) >= irm1) .and. (pl%levelg(j) >= irm1)
                  else
                     lgoodlevel(k) = (pl%levelg(i) >= irm1) .and. (tp%levelg(j) >= irm1)
                  end if
                  lgoodlevel(k) = (self%status(k) == ACTIVE) .and. lgoodlevel(k)
               end do
               ngood = count(lgoodlevel(:))
               if (ngood > 0_I8B) then
                  allocate(good_idx(ngood))
                  good_idx(:) = pack([(i, i = 1, nenc)], lgoodlevel(:))

                  if (isplpl) then
                     do concurrent (k = 1:ngood)
                        i = self%index1(good_idx(k))
                        j = self%index2(good_idx(k))
                        pl%ah(:,i) = 0.0_DP
                        pl%ah(:,j) = 0.0_DP
                     end do
                  else
                     do concurrent (k = 1_I8B:ngood)
                        j = self%index2(good_idx(k))
                        tp%ah(:,j) = 0.0_DP
                     end do
                  end if

                  do k = 1, ngood
                     i = self%index1(good_idx(k))
                     j = self%index2(good_idx(k))
                     if (isplpl) then
                        ri = ((pl%rhill(i)  + pl%rhill(j))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                        rim1 = ri * (RSHELL**2)
                        dx(:) = pl%xh(:,j) - pl%xh(:,i)
                     else
                        ri = ((pl%rhill(i))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                        rim1 = ri * (RSHELL**2)
                        dx(:) = tp%xh(:,j) - pl%xh(:,i)
                     end if
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
                     if (isplpl) then
                        facj = fac * pl%Gmass(j)
                        pl%ah(:, i) = pl%ah(:, i) + facj * dx(:)
                        pl%ah(:, j) = pl%ah(:, j) - faci * dx(:)
                     else
                        tp%ah(:, j) = tp%ah(:, j) - faci * dx(:)
                     end if
                  end do
                  ngood = count(lgoodlevel(:))
                  if (ngood == 0_I8B) return
                  good_idx(1:ngood) = pack([(i, i = 1, nenc)], lgoodlevel(:))

                  if (isplpl) then
                     do k = 1, ngood
                        i = self%index1(good_idx(k))
                        j = self%index2(good_idx(k))
                        pl%vb(:,i) = pl%vb(:,i) + sgn * dt * pl%ah(:,i)
                        pl%vb(:,j) = pl%vb(:,j) + sgn * dt * pl%ah(:,j)
                        pl%ah(:,i) = 0.0_DP
                        pl%ah(:,j) = 0.0_DP
                     end do
                  else
                     do k = 1, ngood
                        j = self%index2(good_idx(k))
                        tp%vb(:,j) = tp%vb(:,j) + sgn * dt * tp%ah(:,j)
                        tp%ah(:,j) = 0.0_DP
                     end do
                  end if
               end if
            end associate
         end select
      end select
      
      return
   end subroutine symba_kick_encounter

end submodule s_symba_kick