submodule(symba_classes) s_symba_kick
   use swiftest
contains

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
      integer(I4B)              :: k
      real(DP)                  :: irij3, rji2, rlim2, faci, facj
      real(DP), dimension(NDIM) :: dx

      if (self%nbody == 0) return
      select type(system)
      class is (symba_nbody_system)
         associate(pl => self, cb => system%cb, plplenc_list => system%plplenc_list, nplplenc => system%plplenc_list%nenc)
            call helio_kick_getacch_pl(pl, system, param, t, lbeg)
            ! Remove accelerations from encountering pairs
            do k = 1, nplplenc
               associate(i => plplenc_list%index1(k), j => plplenc_list%index2(k))
                  dx(:) = pl%xh(:, j) - pl%xh(:, i)
                  rji2 = dot_product(dx(:), dx(:))
                  irij3 = 1.0_DP / (rji2 * sqrt(rji2))
                  faci = pl%Gmass(i) * irij3
                  facj = pl%Gmass(j) * irij3
                  pl%ah(:, i) = pl%ah(:, i) - facj * dx(:)
                  pl%ah(:, j) = pl%ah(:, j) + faci * dx(:)
               end associate
            end do
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


   module subroutine symba_kick_pltpenc(self, system, dt, irec, sgn)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of massive bodies and ACTIVE test particles within SyMBA recursion.
      !! Note: This method works for the polymorphic symba_pltpenc and symba_plplenc types
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_kick.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick.f
      implicit none
      ! Arguments
      class(symba_pltpenc),      intent(in)   :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
      real(DP),                  intent(in)   :: dt    !! step size
      integer(I4B),              intent(in)   :: irec   !! Current recursion level
      integer(I4B),              intent(in)   :: sgn   !! sign to be applied to acceleration
      ! Internals
      integer(I4B)              :: k, irm1, irecl
      real(DP)                  :: r, rr, ri, ris, rim1, r2, ir3, fac, faci, facj
      real(DP), dimension(NDIM) :: dx
      logical                   :: isplpl, lgoodlevel

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
            associate(ind1 => self%index1, ind2 => self%index2)
               if (pl%nbody > 0) pl%lmask(:) = pl%status(:) /= INACTIVE
               if (tp%nbody > 0) tp%lmask(:) = tp%status(:) /= INACTIVE

               irm1 = irec - 1
               if (sgn < 0) then
                  irecl = irec - 1
               else
                  irecl = irec
               end if
               if (isplpl) then
                  pl%ah(:,ind1(1:self%nenc)) = 0.0_DP
                  pl%ah(:,ind2(1:self%nenc)) = 0.0_DP
               else
                  tp%ah(:,ind2(1:self%nenc)) = 0.0_DP
               end if
               do k = 1, self%nenc
                  if (isplpl) then
                     lgoodlevel = (pl%levelg(ind1(k)) >= irm1) .and. (pl%levelg(ind2(k)) >= irm1)
                  else
                     lgoodlevel = (pl%levelg(ind1(k)) >= irm1) .and. (tp%levelg(ind2(k)) >= irm1)
                  end if
                  if ((self%status(k) /= INACTIVE) .and. lgoodlevel) then
                     if (isplpl) then
                        ri = ((pl%rhill(ind1(k))  + pl%rhill(ind2(k)))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                        rim1 = ri * (RSHELL**2)
                        dx(:) = pl%xh(:,ind2(k)) - pl%xh(:,ind1(k))
                     else
                        ri = ((pl%rhill(ind1(k)))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                        rim1 = ri * (RSHELL**2)
                        dx(:) = tp%xh(:,ind2(k)) - pl%xh(:,ind1(k))
                     end if
                     r2 = dot_product(dx(:), dx(:))
                     if (r2 < rim1) then
                        fac = 0.0_DP
                     else if (r2 < ri) then
                        ris = sqrt(ri)
                        r = sqrt(r2)
                        rr = (ris - r) / (ris * (1.0_DP - RSHELL))
                        fac = (r2**(-1.5_DP)) * (1.0_DP - 3 * (rr**2) + 2 * (rr**3))
                     else
                        ir3 = 1.0_DP / (r2 * sqrt(r2))
                        fac = ir3
                     end if
                     faci = fac * pl%Gmass(ind1(k))
                     if (isplpl) then
                        facj = fac * pl%Gmass(ind2(k))
                        pl%ah(:,ind1(k)) = pl%ah(:,ind1(k)) + facj * dx(:)
                        pl%ah(:,ind2(k)) = pl%ah(:,ind2(k)) - faci * dx(:)
                     else
                        tp%ah(:,ind2(k)) = tp%ah(:,ind2(k)) - faci * dx(:)
                     end if
                  end if
               end do
               if (isplpl) then
                  do k = 1, self%nenc
                     pl%vb(:,ind1(k)) = pl%vb(:,ind1(k)) + sgn * dt * pl%ah(:,ind1(k))
                     pl%vb(:,ind2(k)) = pl%vb(:,ind2(k)) + sgn * dt * pl%ah(:,ind2(k))
                     pl%ah(:,ind1(k)) = 0.0_DP
                     pl%ah(:,ind1(k)) = 0.0_DP
                  end do
               else
                  do k = 1, self%nenc
                     tp%vb(:,ind2(k)) = tp%vb(:,ind2(k)) + sgn * dt * tp%ah(:,ind2(k))
                     tp%ah(:,ind2(k)) = 0.0_DP
                  end do
               end if
            end associate
         end select
      end select
      
      return
   end subroutine symba_kick_pltpenc

end submodule s_symba_kick