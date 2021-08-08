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
      integer(I4B) :: i, j
      integer(I8B) :: k
      real(DP) :: rmag, v2, rot2, oblpot, hx, hy, hz, hsx, hsy, hsz
      real(DP), dimension(self%pl%nbody) :: irh, kepl, kespinpl, pecb
      real(DP), dimension(self%pl%nbody) :: Lplorbitx, Lplorbity, Lplorbitz
      real(DP), dimension(self%pl%nbody) :: Lplspinx, Lplspiny, Lplspinz
      real(DP), dimension(self%pl%nplpl) :: pepl 
      logical, dimension(self%pl%nplpl) :: lstatpl
      logical, dimension(self%pl%nbody) :: lstatus

      associate(system => self, pl => self%pl, npl => self%pl%nbody, cb => self%cb)
         system%Lorbit(:) = 0.0_DP
         system%Lspin(:) = 0.0_DP
         system%ke_orbit = 0.0_DP
         system%ke_spin = 0.0_DP

         kepl(:) = 0.0_DP
         Lplorbitx(:) = 0.0_DP
         Lplorbity(:) = 0.0_DP
         Lplorbitz(:) = 0.0_DP
         Lplspinx(:) = 0.0_DP
         Lplspiny(:) = 0.0_DP
         Lplspinz(:) = 0.0_DP
         lstatus(1:npl) = pl%status(1:npl) /= INACTIVE
         call pl%h2b(cb)
         !!$omp simd private(v2, rot2, hx, hy, hz)
         do i = 1, npl
            v2 = dot_product(pl%vb(:,i), pl%vb(:,i))
            hx = pl%xb(2,i) * pl%vb(3,i) - pl%xb(3,i) * pl%vb(2,i)
            hy = pl%xb(3,i) * pl%vb(1,i) - pl%xb(1,i) * pl%vb(3,i)
            hz = pl%xb(1,i) * pl%vb(2,i) - pl%xb(2,i) * pl%vb(1,i)

            ! Angular momentum from orbit 
            Lplorbitx(i) = pl%mass(i) * hx
            Lplorbity(i) = pl%mass(i) * hy
            Lplorbitz(i) = pl%mass(i) * hz
   
            ! Kinetic energy from orbit and spin
            kepl(i) = pl%mass(i) * v2
         end do

         if (param%lrotation) then
            do i = 1, npl
               rot2 = dot_product(pl%rot(:,i), pl%rot(:,i))
               ! For simplicity, we always assume that the rotation pole is the 3rd principal axis
               hsx = pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(1,i) 
               hsy = pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(2,i) 
               hsz = pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(3,i) 

               ! Angular momentum from spin
               Lplspinx(i) = pl%mass(i) * hsx
               Lplspiny(i) = pl%mass(i) * hsy
               Lplspinz(i) = pl%mass(i) * hsz
               kespinpl(i) = pl%mass(i) * pl%Ip(3, i) * pl%radius(i)**2 * rot2
            end do
         else
            kespinpl(:) = 0.0_DP
         end if
   
         ! Do the central body potential energy component first
         !$omp simd 
         do i = 1, npl
            associate(px => pl%xb(1,i), py => pl%xb(2,i), pz => pl%xb(3,i))
               pecb(i) = -cb%Gmass * pl%mass(i) / sqrt(px**2 + py**2 + pz**2)
            end associate
         end do
   
         ! Do the potential energy between pairs of massive bodies
         do k = 1, pl%nplpl
            associate(ik => pl%k_plpl(1, k), jk => pl%k_plpl(2, k))
               pepl(k) = -pl%Gmass(ik) * pl%mass(jk) / norm2(pl%xb(:, jk) - pl%xb(:, ik)) 
               lstatpl(k) = (lstatus(ik) .and. lstatus(jk))
            end associate
         end do
   
         system%ke_orbit = 0.5_DP * sum(kepl(1:npl), lstatus(:))
         if (param%lrotation) system%ke_spin = 0.5_DP * sum(kespinpl(1:npl), lstatus(:))
   
         system%pe = sum(pepl(:), lstatpl(:)) + sum(pecb(1:npl), lstatus(1:npl))
   
         ! Potential energy from the oblateness term
         if (param%loblatecb) then
            !$omp simd 
            do i = 1, npl
               irh(i) = 1.0_DP / norm2(pl%xh(:,i))
            end do
            call obl_pot(npl, cb%Gmass, pl%mass, cb%j2rp2, cb%j4rp4, pl%xh, irh, oblpot)
            system%pe = system%pe + oblpot
         end if
   
         system%Lorbit(1) = sum(Lplorbitx(1:npl), lstatus(1:npl)) 
         system%Lorbit(2) = sum(Lplorbity(1:npl), lstatus(1:npl)) 
         system%Lorbit(3) = sum(Lplorbitz(1:npl), lstatus(1:npl)) 
   
         system%Lspin(1) = sum(Lplspinx(1:npl), lstatus(1:npl)) 
         system%Lspin(2) = sum(Lplspiny(1:npl), lstatus(1:npl)) 
         system%Lspin(3) = sum(Lplspinz(1:npl), lstatus(1:npl)) 
      end associate

      return
   end subroutine util_get_energy_momentum_system

end submodule s_util_get_energy_momentum
