submodule(symba_classes) s_symba_kick
   use swiftest
contains

   module subroutine symba_kick_getacch_int_pl(self)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of massive bodies, with no mutual interactions between bodies below GMTINY
      !!
      !! Adapted from Hal Levison's Swift routine symba5_helio_getacch.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_kick_getacch_int.f90
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self
      ! Internals
      integer(I8B)                      :: k, nplplm
      real(DP)                          :: rji2, rlim2
      real(DP)                          :: dx, dy, dz
      integer(I4B) :: i, j
      real(DP), dimension(:,:), pointer :: ah, xh
      real(DP), dimension(NDIM,self%nbody) :: ahi, ahj
      integer(I4B), dimension(:,:), pointer :: k_plpl
      logical, dimension(:), pointer :: lmask
      real(DP), dimension(:), pointer :: Gmass

      associate(ah => self%ah, xh => self%xh, k_plpl => self%k_plpl, lmask => self%lmask, Gmass => self%Gmass)
         nplplm = self%nplplm
         ahi(:,:) = 0.0_DP
         ahj(:,:) = 0.0_DP
         !$omp parallel do default(shared)&
         !$omp private(k, i, j, dx, dy, dz, rji2)  &
         !$omp reduction(+:ahi) &
         !$omp reduction(-:ahj) 
         do k = 1_I8B, nplplm
            i = k_plpl(1,k)
            j = k_plpl(2,k)
            dx = xh(1, j) - xh(1, i)
            dy = xh(2, j) - xh(2, i)
            dz = xh(3, j) - xh(3, i)
            rji2 = dx**2 + dy**2 + dz**2
            if (lmask(i) .and. lmask(j)) call kick_getacch_int_one_pl(rji2, dx, dy, dz, Gmass(i), Gmass(j), ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
         end do
         !$omp end parallel do
         ah(:,:) = ah(:,:) + ahi(:,:) + ahj(:,:)
      end associate

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
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)              :: i, j
      integer(I8B)              :: k, nplplenc
      real(DP)                  :: rji2, dx, dy, dz
      real(DP), dimension(NDIM,self%nbody) :: ahi, ahj
      class(symba_plplenc), pointer :: plplenc_list

      if (self%nbody == 0) return
      select type(system)
      class is (symba_nbody_system)
         associate(pl => self, xh => self%xh, ah => self%ah, Gmass => self%Gmass, plplenc_list => system%plplenc_list)
            call helio_kick_getacch_pl(pl, system, param, t, lbeg)
            ! Remove accelerations from encountering pairs
            nplplenc = int(plplenc_list%nenc, kind=I8B)
            ahi(:,:) = 0.0_DP
            ahj(:,:) = 0.0_DP
            !$omp parallel do default(shared)&
            !$omp private(k, i, j, dx, dy, dz, rji2)  &
            !$omp reduction(+:ahi) &
            !$omp reduction(-:ahj) 
            do k = 1_I8B, nplplenc
               i = plplenc_list%index1(k)
               j = plplenc_list%index2(k)
               dx = xh(1, j) - xh(1, i)
               dy = xh(2, j) - xh(2, i)
               dz = xh(3, j) - xh(3, i)
               rji2 = dx**2 + dy**2 + dz**2
               call kick_getacch_int_one_pl(rji2, dx, dy, dz, Gmass(i), Gmass(j), ahi(1,i), ahi(2,i), ahi(3,i), ahj(1,j), ahj(2,j), ahj(3,j))
            end do
            !$omp end parallel do
            ah(:,:) = ah(:,:) - ahi(:,:) - ahj(:,:)
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
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)              :: k
      real(DP)                  :: rji2, fac, rlim2
      real(DP), dimension(NDIM) :: dx

      if (self%nbody == 0) return
      select type(system)
      class is (symba_nbody_system)
         associate(tp => self, cb => system%cb, pl => system%pl, pltpenc_list => system%pltpenc_list, npltpenc => system%pltpenc_list%nenc)
            call helio_kick_getacch_tp(tp, system, param, t, lbeg)
            ! Remove accelerations from encountering pairs
            do k = 1, npltpenc
               associate(i => pltpenc_list%index1(k), j => pltpenc_list%index2(k))
                  if (tp%lmask(j)) THEN
                     if (lbeg) then
                        dx(:) = tp%xh(:,j) - pl%xbeg(:,i)
                     else
                        dx(:) = tp%xh(:,j) - pl%xend(:,i)
                     end if
                     rji2 = dot_product(dx(:), dx(:))
                     fac = pl%Gmass(i) / (rji2 * sqrt(rji2))
                     tp%ah(:,j) = tp%ah(:,j) + fac * dx(:)
                  end IF
               end associate
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
      integer(I4B)              :: i, irm1, irecl, ngood
      real(DP)                  :: r, rr, ri, ris, rim1, r2, ir3, fac, faci, facj
      real(DP), dimension(NDIM) :: dx
      logical                   :: isplpl
      logical, dimension(self%nenc) :: lgoodlevel
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
            associate(ind1 => self%index1, ind2 => self%index2, npl => pl%nbody, ntp => tp%nbody, nenc => self%nenc)
               if (npl > 0) pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE
               if (ntp > 0) tp%lmask(1:ntp) = tp%status(1:ntp) /= INACTIVE

               irm1 = irec - 1

               if (sgn < 0) then
                  irecl = irec - 1
               else
                  irecl = irec
               end if

               do i = 1, nenc
                  if (isplpl) then
                     lgoodlevel(i) = (pl%levelg(ind1(i)) >= irm1) .and. (pl%levelg(ind2(i)) >= irm1)
                  else
                     lgoodlevel(i) = (pl%levelg(ind1(i)) >= irm1) .and. (tp%levelg(ind2(i)) >= irm1)
                  end if
                  lgoodlevel(i) = (self%status(i) == ACTIVE) .and. lgoodlevel(i)
               end do
               ngood = count(lgoodlevel(:))
               if (ngood > 0) then
                  allocate(good_idx(ngood))
                  good_idx(:) = pack([(i, i = 1, nenc)], lgoodlevel(:))

                  if (isplpl) then
                     do i = 1, ngood
                        associate(i1 => ind1(good_idx(i)), i2 => ind2(good_idx(i)))
                           pl%ah(:,i1) = 0.0_DP
                           pl%ah(:,i2) = 0.0_DP
                        end associate
                     end do
                  else
                     do i = 1, ngood
                        associate(i2 => ind2(good_idx(i)))
                           tp%ah(:,i2) = 0.0_DP
                        end associate
                     end do
                  end if

                  do i = 1, ngood
                     associate(k => good_idx(i), i1 => ind1(good_idx(i)), i2 => ind2(good_idx(i)))
                        if (isplpl) then
                           ri = ((pl%rhill(i1)  + pl%rhill(i2))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                           rim1 = ri * (RSHELL**2)
                           dx(:) = pl%xh(:,i2) - pl%xh(:,i1)
                        else
                           ri = ((pl%rhill(i1))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                           rim1 = ri * (RSHELL**2)
                           dx(:) = tp%xh(:,i2) - pl%xh(:,i1)
                        end if
                        r2 = dot_product(dx(:), dx(:))
                        if (r2 < rim1) then
                           fac = 0.0_DP
                           lgoodlevel(k) = .false.
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
                        faci = fac * pl%Gmass(i1)
                        if (isplpl) then
                           facj = fac * pl%Gmass(i2)
                           pl%ah(:, i1) = pl%ah(:, i1) + facj * dx(:)
                           pl%ah(:, i2) = pl%ah(:, i2) - faci * dx(:)
                        else
                           tp%ah(:, i2) = tp%ah(:, i2) - faci * dx(:)
                        end if
                     end associate
                  end do
                  ngood = count(lgoodlevel(:))
                  if (ngood == 0) return
                  good_idx(1:ngood) = pack([(i, i = 1, nenc)], lgoodlevel(:))

                  if (isplpl) then
                     do i = 1, ngood
                        associate(i1 => ind1(good_idx(i)), i2 => ind2(good_idx(i)))
                           pl%vb(:,i1) = pl%vb(:,i1) + sgn * dt * pl%ah(:,i1)
                           pl%vb(:,i2) = pl%vb(:,i2) + sgn * dt * pl%ah(:,i2)
                           pl%ah(:,i1) = 0.0_DP
                           pl%ah(:,i2) = 0.0_DP
                        end associate
                     end do
                  else
                     do i = 1, ngood
                        associate(i2 => ind2(good_idx(i)))
                           tp%vb(:,i2) = tp%vb(:,i2) + sgn * dt * tp%ah(:,i2)
                           tp%ah(:,i2) = 0.0_DP
                        end associate
                     end do
                  end if
               end if
            end associate
         end select
      end select
      
      return
   end subroutine symba_kick_encounter

end submodule s_symba_kick