submodule(symba_classes) s_symba_kick
   use swiftest
contains

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
      class(symba_pltpenc),     intent(in)   :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
      real(DP),              intent(in)   :: dt    !! step size
      integer(I4B),           intent(in)   :: irec   !! Current recursion level
      integer(I4B),           intent(in)   :: sgn   !! sign to be applied to acceleration
      ! Internals
      integer(I4B)           :: i, irm1, irecl
      real(DP)              :: r, rr, ri, ris, rim1, r2, ir3, fac, faci, facj
      real(DP), dimension(NDIM) :: dx
      logical               :: isplpl, lgoodlevel

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
            do i = 1, self%nenc
               associate(index_i => self%index1(i), index_j => self%index2(i))
                  if (isplpl) then
                     pl%ah(:,index_i) = 0.0_DP
                     pl%ah(:,index_j) = 0.0_DP
                  else
                     tp%ah(:,index_j) = 0.0_DP
                  end if
                  if (isplpl) then
                     lgoodlevel = (pl%levelg(index_i) >= irm1) .and. (pl%levelg(index_j) >= irm1)
                  else
                     lgoodlevel = (pl%levelg(index_i) >= irm1) .and. (tp%levelg(index_j) >= irm1)
                  end if
                  if ((self%status(i) == ACTIVE) .and. lgoodlevel) then
                     if (isplpl) then
                        ri = ((pl%rhill(index_i)  + pl%rhill(index_j))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                        rim1 = ri * (RSHELL**2)
                        dx(:) = pl%xh(:,index_j) - pl%xh(:,index_i)
                     else
                        ri = ((pl%rhill(index_i))**2) * (RHSCALE**2) * (RSHELL**(2*irecl))
                        rim1 = ri * (RSHELL**2)
                        dx(:) = tp%xh(:,index_j) - pl%xh(:,index_i)
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
                     faci = fac * pl%mass(index_i)
                     if (isplpl) then
                        facj = fac * pl%mass(index_j)
                        pl%ah(:,index_i) = pl%ah(:,index_i) + facj*dx(:)
                        pl%ah(:,index_j) = pl%ah(:,index_j) - faci*dx(:)
                     else
                        tp%ah(:,index_j) = tp%ah(:,index_j) - faci*dx(:)
                     end if
                  end if
               end associate
            end do
            if (isplpl) then
               do i = 1, self%nenc
                  associate(index_i => self%index1(i), index_j => self%index2(i))
                     pl%vb(:,index_i) = pl%vb(:,index_i) + sgn * dt * pl%ah(:,index_i)
                     pl%vb(:,index_j) = pl%vb(:,index_j) + sgn * dt * pl%ah(:,index_j)
                     pl%ah(:,index_i) = 0.0_DP
                     pl%ah(:,index_j) = 0.0_DP
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