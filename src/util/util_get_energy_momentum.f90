!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_get_energy_momentum
   use swiftest
contains   
   module subroutine util_get_energy_momentum_system(self, param)
      !! author: David A. Minton
      !!
      !! Compute total system angular momentum vector and kinetic, potential and total system energy
      !!  
      !! Adapted from David E. Kaufmann Swifter routine symba_energy_eucl.f90
      !!  
      !! Adapted from Martin Duncan's Swift routine anal_energy.f
      implicit none
      class(swiftest_nbody_system), intent(inout) :: self     !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param    !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i
      real(DP) :: kecb, kespincb
      real(DP), dimension(self%pl%nbody) :: kepl, kespinpl
      real(DP), dimension(self%pl%nbody) :: Lplorbitx, Lplorbity, Lplorbitz
      real(DP), dimension(self%pl%nbody) :: Lplspinx, Lplspiny, Lplspinz
      real(DP), dimension(NDIM) :: Lcborbit, Lcbspin
      real(DP) :: hx, hy, hz

      associate(system => self, pl => self%pl, npl => self%pl%nbody, cb => self%cb)
         system%Lorbit(:) = 0.0_DP
         system%Lspin(:) = 0.0_DP
         system%Ltot(:) = 0.0_DP
         system%ke_orbit = 0.0_DP
         system%ke_spin = 0.0_DP

         kepl(:) = 0.0_DP
         Lplorbitx(:) = 0.0_DP
         Lplorbity(:) = 0.0_DP
         Lplorbitz(:) = 0.0_DP
         Lplspinx(:) = 0.0_DP
         Lplspiny(:) = 0.0_DP
         Lplspinz(:) = 0.0_DP

         pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE

         system%GMtot = cb%Gmass + sum(pl%Gmass(1:npl), pl%lmask(1:npl)) 
         kecb = cb%mass * dot_product(cb%vb(:), cb%vb(:))
         Lcborbit(:) = cb%mass * (cb%rb(:) .cross. cb%vb(:))

         do concurrent (i = 1:npl, pl%lmask(i))
            hx = pl%rb(2,i) * pl%vb(3,i) - pl%rb(3,i) * pl%vb(2,i)
            hy = pl%rb(3,i) * pl%vb(1,i) - pl%rb(1,i) * pl%vb(3,i)
            hz = pl%rb(1,i) * pl%vb(2,i) - pl%rb(2,i) * pl%vb(1,i)

            ! Angular momentum from orbit 
            Lplorbitx(i) = pl%mass(i) * hx
            Lplorbity(i) = pl%mass(i) * hy
            Lplorbitz(i) = pl%mass(i) * hz

            ! Kinetic energy from orbit
            kepl(i) = pl%mass(i) * dot_product(pl%vb(:,i), pl%vb(:,i)) 
         end do

         if (param%lrotation) then
            kespincb = cb%mass * cb%Ip(3) * cb%radius**2 * dot_product(cb%rot(:), cb%rot(:))

            ! For simplicity, we always assume that the rotation pole is the 3rd principal axis
            Lcbspin(:) = cb%Ip(3) * cb%mass * cb%radius**2 * cb%rot(:)

            do concurrent (i = 1:npl, pl%lmask(i))
               ! Currently we assume that the rotation pole is the 3rd principal axis
               ! Angular momentum from spin
               Lplspinx(i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(1,i)
               Lplspiny(i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(2,i)  
               Lplspinz(i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(3,i)  

               ! Kinetic energy from spin
               kespinpl(i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * dot_product(pl%rot(:,i), pl%rot(:,i))
            end do
         else
            kespincb = 0.0_DP
            kespinpl(:) = 0.0_DP
         end if
  
         if (param%lflatten_interactions) then
            call util_get_energy_potential_flat(npl, pl%nplpl, pl%k_plpl, pl%lmask, cb%Gmass, pl%Gmass, pl%mass, pl%rb, system%pe)
         else
            call util_get_energy_potential_triangular(npl, pl%lmask, cb%Gmass, pl%Gmass, pl%mass, pl%rb, system%pe)
         end if

         ! Potential energy from the oblateness term
         if (param%loblatecb) then
            call system%obl_pot()
            system%pe = system%pe + system%oblpot
         end if

         system%ke_orbit = 0.5_DP * (kecb + sum(kepl(1:npl), pl%lmask(1:npl)))
         if (param%lrotation) system%ke_spin = 0.5_DP * (kespincb + sum(kespinpl(1:npl), pl%lmask(1:npl)))
   
         system%Lorbit(1) = Lcborbit(1) + sum(Lplorbitx(1:npl), pl%lmask(1:npl)) 
         system%Lorbit(2) = Lcborbit(2) + sum(Lplorbity(1:npl), pl%lmask(1:npl)) 
         system%Lorbit(3) = Lcborbit(3) + sum(Lplorbitz(1:npl), pl%lmask(1:npl)) 
  
         if (param%lrotation) then
            system%Lspin(1) = Lcbspin(1) + sum(Lplspinx(1:npl), pl%lmask(1:npl)) 
            system%Lspin(2) = Lcbspin(2) + sum(Lplspiny(1:npl), pl%lmask(1:npl)) 
            system%Lspin(3) = Lcbspin(3) + sum(Lplspinz(1:npl), pl%lmask(1:npl)) 
         end if

         system%te = system%ke_orbit + system%ke_spin + system%pe
         system%Ltot(:) = system%Lorbit(:) + system%Lspin(:)
      end associate

      return
   end subroutine util_get_energy_momentum_system


   subroutine util_get_energy_potential_flat(npl, nplpl, k_plpl, lmask, GMcb, Gmass, mass, rb, pe)
      !! author: David A. Minton
      !!
      !! Compute total system potential energy
      implicit none
      ! Arguments
      integer(I4B),                 intent(in)  :: npl
      integer(I8B),                 intent(in)  :: nplpl
      integer(I4B), dimension(:,:), intent(in)  :: k_plpl
      logical,      dimension(:),   intent(in)  :: lmask
      real(DP),                     intent(in)  :: GMcb
      real(DP),     dimension(:),   intent(in)  :: Gmass
      real(DP),     dimension(:),   intent(in)  :: mass
      real(DP),     dimension(:,:), intent(in)  :: rb
      real(DP),                     intent(out) :: pe
      ! Internals
      integer(I4B) :: i, j
      integer(I8B) :: k
      real(DP), dimension(npl) :: pecb
      real(DP), dimension(nplpl) :: pepl 
      logical, dimension(nplpl) :: lstatpl

      ! Do the central body potential energy component first
      where(.not. lmask(1:npl))
         pecb(1:npl) = 0.0_DP
      end where

      do concurrent(i = 1:npl, lmask(i))
         pecb(i) = -GMcb * mass(i) / norm2(rb(:,i)) 
      end do

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(k_plpl, rb, mass, Gmass, pepl, lstatpl, lmask) &
      !$omp firstprivate(nplpl)
      do k = 1, nplpl
         i = k_plpl(1,k)
         j = k_plpl(2,k)
         lstatpl(k) = (lmask(i) .and. lmask(j))
         if (lstatpl(k)) then
            pepl(k) = -(Gmass(i) * mass(j)) / norm2(rb(:, i) - rb(:, j))
         else
            pepl(k) = 0.0_DP
         end if
      end do
      !$omp end parallel do 

      pe = sum(pepl(:), lstatpl(:)) + sum(pecb(1:npl), lmask(1:npl))

      return
   end subroutine util_get_energy_potential_flat


   subroutine util_get_energy_potential_triangular(npl, lmask, GMcb, Gmass, mass, rb, pe)
      !! author: David A. Minton
      !!
      !! Compute total system potential energy
      implicit none
      ! Arguments
      integer(I4B),                 intent(in)  :: npl
      logical,      dimension(:),   intent(in)  :: lmask
      real(DP),                     intent(in)  :: GMcb
      real(DP),     dimension(:),   intent(in)  :: Gmass
      real(DP),     dimension(:),   intent(in)  :: mass
      real(DP),     dimension(:,:), intent(in)  :: rb
      real(DP),                     intent(out) :: pe
      ! Internals
      integer(I4B) :: i, j
      real(DP), dimension(npl) :: pecb, pepl

      ! Do the central body potential energy component first
      where(.not. lmask(1:npl))
         pecb(1:npl) = 0.0_DP
      end where

      do concurrent(i = 1:npl, lmask(i))
         pecb(i) = -GMcb * mass(i) / norm2(rb(:,i)) 
      end do

      pe = 0.0_DP
      !$omp parallel do default(private) schedule(static)&
      !$omp shared(lmask, Gmass, mass, rb) &
      !$omp firstprivate(npl) &
      !$omp reduction(+:pe) 
      do i = 1, npl
         if (lmask(i)) then
            do concurrent(j = i+1:npl, lmask(i) .and. lmask(j))
               pepl(j) = - (Gmass(i) * mass(j)) / norm2(rb(:, i) - rb(:, j))
            end do
            pe = pe + sum(pepl(i+1:npl), lmask(i+1:npl))
         end if
      end do
      !$omp end parallel do
      pe = pe + sum(pecb(1:npl), lmask(1:npl))

      return
   end subroutine util_get_energy_potential_triangular

end submodule s_util_get_energy_momentum
