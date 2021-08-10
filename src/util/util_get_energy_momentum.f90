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
      real(DP) :: oblpot, kecb, kespincb
      real(DP), dimension(self%pl%nbody) :: irh, kepl, kespinpl, pecb
      real(DP), dimension(self%pl%nbody) :: Lplorbitx, Lplorbity, Lplorbitz
      real(DP), dimension(self%pl%nbody) :: Lplspinx, Lplspiny, Lplspinz
      real(DP), dimension(self%pl%nplpl) :: pepl 
      real(DP), dimension(NDIM) :: Lcborbit, Lcbspin
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
         Lcborbit(:) = cb%mass * cb%xb(:) .cross. cb%vb(:)

         do concurrent (i = 1:npl, lstatus(i))
            block ! We use a block construct to prevent generating temporary arrays for local variables
               real(DP) :: v2, hx, hy, hz
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
            end block
         end do

         if (param%lrotation) then
            kespincb = cb%mass * cb%Ip(3) * cb%radius**2 * dot_product(cb%rot(:), cb%rot(:))

            ! For simplicity, we always assume that the rotation pole is the 3rd principal axis
            Lcbspin(:) = cb%Ip(3) * cb%mass * cb%radius**2 * cb%rot(:)

            do concurrent (i = 1:npl, lstatus(i))
               block 
                  real(DP) :: rot2, hsx, hsy, hsz

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
               end block
            end do
         else
            kespincb = 0.0_DP
            kespinpl(:) = 0.0_DP
         end if
   
         ! Do the central body potential energy component first
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
            do concurrent(i = 1:npl, lstatus(i))
               irh(i) = 1.0_DP / norm2(pl%xh(:,i))
            end do
            call obl_pot(npl, cb%Gmass, pl%mass, cb%j2rp2, cb%j4rp4, pl%xh, irh, oblpot)
            system%pe = system%pe + oblpot
         end if
   
         system%Lorbit(1) = Lcborbit(1) + sum(Lplorbitx(1:npl), lstatus(1:npl)) 
         system%Lorbit(2) = Lcborbit(2) + sum(Lplorbity(1:npl), lstatus(1:npl)) 
         system%Lorbit(3) = Lcborbit(3) + sum(Lplorbitz(1:npl), lstatus(1:npl)) 
  
         if (param%lrotation) then
            system%Lspin(1) = Lcbspin(1) + sum(Lplspinx(1:npl), lstatus(1:npl)) 
            system%Lspin(2) = Lcbspin(2) + sum(Lplspiny(1:npl), lstatus(1:npl)) 
            system%Lspin(3) = Lcbspin(3) + sum(Lplspinz(1:npl), lstatus(1:npl)) 
         end if
      end associate

      return
   end subroutine util_get_energy_momentum_system

end submodule s_util_get_energy_momentum
