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
   logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
   ! Internals
   integer(I4B)              :: k
   real(DP)                  :: irij3, rji2, rlim2, faci, facj
   real(DP), dimension(NDIM) :: dx

   select type(system)
   class is (symba_nbody_system)
      associate(pl => self, cb => system%cb, plplenc_list => system%plplenc_list, nplplenc => system%plplenc_list%nenc)
         ! Remove accelerations from encountering pairs
         do k = 1, nplplenc
            associate(i => plplenc_list%index1(k), j => plplenc_list%index2(k))
               dx(:) = pl%xh(:, j) - pl%xh(:, i)
               rji2 = dot_product(dx(:), dx(:))
               rlim2 = (pl%radius(i) + pl%radius(j))**2
               if (rji2 > rlim2) then
                  irij3 = 1.0_DP / (rji2 * sqrt(rji2))
                  faci = pl%Gmass(i) * irij3
                  facj = pl%Gmass(j) * irij3
                  pl%ah(:, i) = pl%ah(:, i) - facj * dx(:)
                  pl%ah(:, j) = pl%ah(:, j) + faci * dx(:)
               end if
            end associate
         end do
         call helio_kick_getacch_pl(pl, system, param, t, lbeg)
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
      logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)              :: k
      real(DP)                  :: rji2, fac, rlim2
      real(DP), dimension(NDIM) :: dx
  
      select type(system)
      class is (symba_nbody_system)
         associate(tp => self, cb => system%cb, pl => system%pl, pltpenc_list => system%pltpenc_list, npltpenc => system%pltpenc_list%nenc)
            ! Remove accelerations from encountering pairs
            do k = 1, npltpenc
               associate(i => pltpenc_list%index1(k), j => pltpenc_list%index2(k))
                  if (tp%status(j) == ACTIVE) THEN
                     dx(:) = tp%xh(:,j) - pl%xh(:,i)
                     rji2 = dot_product(dx(:), dx(:))
                     rlim2 = (pl%radius(i))**2
                     if (rji2 > rlim2) then
                        fac = pl%Gmass(i) / (rji2 * sqrt(rji2))
                        tp%ah(:,j) = tp%ah(:,j) + fac * dx(:)
                     end if
                  end IF
               end associate
            end do
            call helio_kick_getacch_tp(tp, system, param, t, lbeg)
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
            irm1 = irec - 1
            if (sgn < 0) then
               irecl = irec - 1
            else
               irecl = irec
            end if
            do k = 1, self%nenc
               associate(i => self%index1(k), j => self%index2(k))
                  if (isplpl) then
                     pl%ah(:,i) = 0.0_DP
                     pl%ah(:,j) = 0.0_DP
                  else
                     tp%ah(:,j) = 0.0_DP
                  end if
                  if (isplpl) then
                     lgoodlevel = (pl%levelg(i) >= irm1) .and. (pl%levelg(j) >= irm1)
                  else
                     lgoodlevel = (pl%levelg(i) >= irm1) .and. (tp%levelg(j) >= irm1)
                  end if
                  if ((self%status(i) == ACTIVE) .and. lgoodlevel) then
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
                     else if (r2 < ri) then
                        ris = sqrt(ri)
                        r = sqrt(r2)
                        rr = (ris - r) / (ris * (1.0_DP - RSHELL))
                        fac = (r2**(-1.5_DP)) * (1.0_DP - 3 * (rr**2) + 2 * (rr**3))
                     else
                        ir3 = 1.0_DP / (r2 * sqrt(r2))
                        fac = ir3
                     end if
                     faci = fac * pl%mass(i)
                     if (isplpl) then
                        facj = fac * pl%mass(j)
                        pl%ah(:,i) = pl%ah(:,i) + facj*dx(:)
                        pl%ah(:,j) = pl%ah(:,j) - faci*dx(:)
                     else
                        tp%ah(:,j) = tp%ah(:,j) - faci*dx(:)
                     end if
                  end if
               end associate
            end do
            if (isplpl) then
               do k = 1, self%nenc
                  associate(i => self%index1(k), j => self%index2(k))
                     pl%vb(:,i) = pl%vb(:,i) + sgn * dt * pl%ah(:,i)
                     pl%vb(:,j) = pl%vb(:,j) + sgn * dt * pl%ah(:,j)
                     pl%ah(:,i) = 0.0_DP
                     pl%ah(:,j) = 0.0_DP
                  end associate
               end do
            else
               where(tp%status(self%index2(1:self%nenc)) == ACTIVE)
                  tp%vb(1,self%index2(:)) = tp%vb(1,self%index2(:)) + sgn * dt * tp%ah(1,self%index2(:))
                  tp%vb(2,self%index2(:)) = tp%vb(2,self%index2(:)) + sgn * dt * tp%ah(2,self%index2(:))
                  tp%vb(3,self%index2(:)) = tp%vb(3,self%index2(:)) + sgn * dt * tp%ah(3,self%index2(:))
               end where
               tp%ah(:,self%index2(1:self%nenc)) = 0.0_DP
            end if
         end select
      end select
      return
   end subroutine symba_kick_pltpenc

end submodule s_symba_kick