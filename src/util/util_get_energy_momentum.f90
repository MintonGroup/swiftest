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
      real(DP) :: kecb, kespincb, Lcborbitx, Lcborbity, Lcborbitz, Lcbspinx, Lcbspiny, Lcbspinz
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

         kecb = cb%mass * dot_product(cb%vb(:), cb%vb(:))
         hx = cb%xb(2) * cb%vb(3) - cb%xb(3) * cb%vb(2)
         hy = cb%xb(3) * cb%vb(1) - cb%xb(1) * cb%vb(3)
         hz = cb%xb(1) * cb%vb(2) - cb%xb(2) * cb%vb(1) 
         Lcborbitx = cb%mass * hx
         Lcborbity = cb%mass * hy
         Lcborbitz = cb%mass * hz
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
            kespincb = cb%mass * cb%Ip(3) * cb%radius**2 * dot_product(cb%rot(:), cb%rot(:))
            ! For simplicity, we always assume that the rotation pole is the 3rd principal axis
            hsx = cb%Ip(3) * cb%radius**2 * cb%rot(1) 
            hsy = cb%Ip(3) * cb%radius**2 * cb%rot(2) 
            hsz = cb%Ip(3) * cb%radius**2 * cb%rot(3) 

            ! Angular momentum from spin
            Lcbspinx = cb%mass * hsx
            Lcbspiny = cb%mass * hsy
            Lcbspinz = cb%mass * hsz

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
            kespincb = 0.0_DP
            kespinpl(:) = 0.0_DP
         end if
   
         ! Do the central body potential energy component first
         !$omp simd 
         associate(px => pl%xb(1,:), py => pl%xb(2,:), pz => pl%xb(3,:))
            do concurrent(i = 1:npl, lstatus(i))
               pecb(i) = -cb%Gmass * pl%mass(i) / sqrt(px(i)**2 + py(i)**2 + pz(i)**2)
            end do
         end associate
   
         ! Do the potential energy between pairs of massive bodies
         associate(indi => pl%k_plpl(1, :), indj => pl%k_plpl(2, :))
            do concurrent (k = 1:pl%nplpl)
               lstatpl(k) = (lstatus(indi(k)) .and. lstatus(indj(k)))
            end do

            do concurrent (k = 1:pl%nplpl, lstatpl(k))
               pepl(k) = -pl%Gmass(indi(k)) * pl%mass(indj(k)) / norm2(pl%xb(:, indi(k)) - pl%xb(:, indj(k))) 
            end do
         end associate
   
         system%pe = sum(pepl(:), lstatpl(:)) + sum(pecb(1:npl), lstatus(1:npl))

         system%ke_orbit = 0.5_DP * (kecb + sum(kepl(1:npl), lstatus(:)))
         if (param%lrotation) system%ke_spin = 0.5_DP * (kespincb + sum(kespinpl(1:npl), lstatus(:)))
   
         ! Potential energy from the oblateness term
         if (param%loblatecb) then
            !$omp simd 
            do concurrent(i = 1:npl, lstatus(i))
               irh(i) = 1.0_DP / norm2(pl%xh(:,i))
            end do
            call obl_pot(npl, cb%Gmass, pl%mass, cb%j2rp2, cb%j4rp4, pl%xh, irh, oblpot)
            system%pe = system%pe + oblpot
         end if
   
         system%Lorbit(1) = Lcborbitx + sum(Lplorbitx(1:npl), lstatus(1:npl)) 
         system%Lorbit(2) = Lcborbity + sum(Lplorbity(1:npl), lstatus(1:npl)) 
         system%Lorbit(3) = Lcborbitz + sum(Lplorbitz(1:npl), lstatus(1:npl)) 
  
         if (param%lrotation) then
            system%Lspin(1) = Lcbspinx + sum(Lplspinx(1:npl), lstatus(1:npl)) 
            system%Lspin(2) = Lcbspiny + sum(Lplspiny(1:npl), lstatus(1:npl)) 
            system%Lspin(3) = Lcbspinz + sum(Lplspinz(1:npl), lstatus(1:npl)) 
         end if
      end associate

      return
   end subroutine util_get_energy_momentum_system

end submodule s_util_get_energy_momentum
