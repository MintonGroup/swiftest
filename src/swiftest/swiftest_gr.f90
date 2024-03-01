! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(swiftest) s_swiftest_gr
contains

   pure module subroutine swiftest_gr_kick_getaccb_ns_body(self, nbody_system, param) 
      !! author: David A. Minton
      !!
      !! Add relativistic correction acceleration for non-symplectic integrators.
      !! Based on Quinn et al. (1991) eq. 5
      !!
      !! Reference:  
      !!    Quinn, T.R., Tremaine, S., Duncan, M., 1991. A three million year integration of the earth’s orbit. 
      !!       AJ 101, 2287–2305. https://doi.org/10.1086/115850
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_kick_getaccb_ns.f90
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   
         !! Swiftest generic body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system 
         !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  
         !! Current run configuration parameters 
      ! Internals
      real(DP)                  :: rmag, rdotv, vmag2
      integer(I4B)              :: i

      associate(n => self%nbody, cb => nbody_system%cb, inv_c2 => param%inv_c2)
         if (n == 0) return
         do i = 1, n
            rmag = norm2(self%rh(:,i))
            vmag2 = dot_product(self%vh(:,i), self%vh(:,i))
            rdotv = dot_product(self%rh(:,i), self%vh(:,i))
            self%agr(:, i) = self%mu * inv_c2 / rmag**3 * ((4 * self%mu(i) / rmag - vmag2) &
                           * self%rh(:,i) + 4 * rdotv * self%vh(:,i))
         end do

         select type(self)
         class is (swiftest_pl)
            cb%agr(:) = 0.0_DP
            do i = n, 1, -1
               cb%agr(:) = cb%agr(:) - self%Gmass(i) * self%agr(:, i) / cb%Gmass
            end do
         end select
      end associate 

      return
   end subroutine swiftest_gr_kick_getaccb_ns_body


   pure module subroutine swiftest_gr_kick_getacch(mu, x, lmask, n, inv_c2, agr) 
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of massive bodies
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_kick_getacch.f90
      implicit none
      ! Arguments
      real(DP), dimension(:),     intent(in)    :: mu     
         !! Gravitational constant
      real(DP), dimension(:,:),   intent(in)    :: x      
         !! Position vectors
      logical,  dimension(:),     intent(in)    :: lmask  
         !! Logical mask indicating which bodies to compute
      integer(I4B),               intent(in)    :: n      
         !! Total number of bodies
      real(DP),                   intent(in)    :: inv_c2 
         !! Inverse speed of light squared: 1 / c**2
      real(DP), dimension(:,:),   intent(out)   :: agr    
         !! Accelerations
      ! Internals
      integer(I4B)                              :: i
      real(DP)                                  :: beta, rjmag4
     
      agr(:,:) = 0.0_DP
#ifdef DOCONLOC
      do concurrent (i = 1:n, lmask(i)) shared(lmask,x,mu,agr,inv_c2) local(rjmag4,beta)
#else
      do concurrent (i = 1:n, lmask(i))
#endif
         rjmag4 = (dot_product(x(:, i), x(:, i)))**2
         beta = -mu(i)**2 * inv_c2 
         agr(:, i) = 2 * beta * x(:, i) / rjmag4
      end do

      return
   end subroutine swiftest_gr_kick_getacch


   pure elemental module subroutine swiftest_gr_p4_pos_kick(inv_c2, rx, ry, rz, vx, vy, vz, dt)
      !! author: David A. Minton
      !!
      !! Position kick due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Reference:
      !!    Saha, P., Tremaine, S., 1994. Long-term planetary integration with individual time steps. 
      !!       AJ 108, 1962–1969. https://doi.org/10.1086/117210
      !!
      !! Adapted from David A. Minton's Swifter routine gr_whm_p4.f90  
      implicit none
      ! Arguments
      real(DP), intent(in)    :: inv_c2     
         !! One over speed of light squared (1/c**2)
      real(DP), intent(inout) :: rx, ry, rz  
         !! Position vector
      real(DP), intent(in)    :: vx, vy, vz 
         !! Velocity vector
      real(DP), intent(in)    :: dt         
         !! Step size
      ! Internals
      real(DP) :: drx, dry, drz
      real(DP) :: vmag2

      vmag2 = vx*vx + vy*vy + vz*vz
      drx = -2 * inv_c2 * vmag2 * vx
      dry = -2 * inv_c2 * vmag2 * vy
      drz = -2 * inv_c2 * vmag2 * vz
      rx = rx + drx * dt
      ry = ry + dry * dt
      rz = rz + drz * dt

      return
   end subroutine swiftest_gr_p4_pos_kick


   pure module subroutine swiftest_gr_pseudovel2vel(param, mu, rh, pv, vh) 
      !! author: David A. Minton
      !!
      !! Converts the relativistic pseudovelocity back into a veliocentric velocity
      !!    Based on Saha & Tremaine (1994) Eq. 32
      !! 
      !! Reference:
      !!    Saha, P., Tremaine, S., 1994. Long-term planetary integration with individual time steps. 
      !!       AJ 108, 1962–1969. https://doi.org/10.1086/117210
      !!
      !! Adapted from David A. Minton's Swifter routine gr_pseudovel2vel.f90 
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(in)  :: param 
         !! Current run configuration parameters 
      real(DP),                   intent(in)  :: mu    
         !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
      real(DP), dimension(:),     intent(in)  :: rh    
         !! Heliocentric position vector 
      real(DP), dimension(:),     intent(in)  :: pv    
         !! Pseudovelocity velocity vector - see Saha & Tremain (1994), eq. (32)
      real(DP), dimension(:),     intent(out) :: vh    
         !! Heliocentric velocity vector 
      ! Internals
      real(DP) :: vmag2, rmag, grterm
   
      associate(inv_c2 => param%inv_c2)
         vmag2 = dot_product(pv(:), pv(:))
         rmag  = norm2(rh(:))
         grterm = 1.0_DP - inv_c2 * (0.5_DP * vmag2 + 3 * mu / rmag)
         vh(:) = pv(:) * grterm
      end associate

      return
   end subroutine swiftest_gr_pseudovel2vel


   pure module subroutine swiftest_gr_pv2vh_body(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from pseudovelocity to heliocentric velocity for swiftest bodies
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self  !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      ! Internals
      integer(I4B)                              :: i
      real(DP), dimension(:,:), allocatable     :: vh    !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(n => self%nbody)
         if (n == 0) return
         allocate(vh, mold = self%vh)
         do i = 1, n
            call swiftest_gr_pseudovel2vel(param, self%mu(i), self%rh(:, i), self%vh(:, i), vh(:, i))
         end do
         call move_alloc(vh, self%vh)
      end associate

      return
   end subroutine swiftest_gr_pv2vh_body


   pure module subroutine swiftest_gr_vel2pseudovel(param, mu, rh, vh, pv)
      !! author: David A. Minton
      !!
      !! Converts the heliocentric velocity into a pseudovelocity with relativistic corrections. 
      !! Uses Newton-Raphson method with direct inversion of the Jacobian (yeah, it's slow, but 
      !! this is only done once per run).
      !!
      !! Reference:
      !!    Saha, P., Tremaine, S., 1994. Long-term planetary integration with individual time steps. 
      !!       AJ 108, 1962–1969. https://doi.org/10.1086/117210
      !!
      !! Adapted from David A. Minton's Swifter routine gr_vel2pseudovel.f90
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(in)  :: param 
         !! Current run configuration parameters 
      real(DP),                   intent(in)  :: mu    
         !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
      real(DP), dimension(:),     intent(in)  :: rh    
         !! Heliocentric position vector 
      real(DP), dimension(:),     intent(in)  :: vh    
         !! Heliocentric velocity vector 
      real(DP), dimension(:),     intent(out) :: pv    
         !! Pseudovelocity vector - see Saha & Tremain (1994), eq. (32)
      ! Internals
      real(DP)                       :: v2, G, pv2, rterm, det
      real(DP), dimension(NDIM,NDIM) :: J,Jinv
      real(DP), dimension(NDIM)      :: F
      integer(I4B)                   :: n,i,k
      integer(I4B), parameter        :: MAXITER = 50
      real(DP),parameter             :: TOL = 1.0e-12_DP

      associate(inv_c2 => param%inv_c2)
         pv(1:NDIM) = vh(1:NDIM) ! Initial guess
         rterm = 3 * mu / norm2(rh(:))
         v2 = dot_product(vh(:), vh(:))
         if (v2 < TINY(1.0_DP)) then
            pv(:) = 0.0_DP
            return
         end if
         do n = 1, MAXITER
            pv2 = dot_product(pv(:), pv(:))
            G = 1.0_DP - inv_c2 * (0.5_DP * pv2 + rterm)
            F(:) = pv(:) * G - vh(:)
            if (abs(sum(F) / v2 ) < TOL) exit ! Root found

            ! Calculate the Jacobian
            do k = 1, NDIM
               do i = 1, NDIM
                  if (i == k) then
                     J(i,k) = G - inv_c2 * pv(k)
                  else
                     J(i,k) = -inv_c2 * pv(k)
                  end if
               end do
            end do

            ! Inverse of the Jacobian
            det = J(1,1) * (J(3,3) * J(2,2) - J(3,2) * J(2,3))
            det = det - J(2,1) * (J(3,3) * J(1,2)-J(3,2) * J(1,3))
            det = det + J(3,1) * (J(2,3) * J(1,2)-J(2,2) * J(1,3))

            Jinv(1,1) =   J(3,3) * J(2,2) - J(3,2) * J(2,3)
            Jinv(1,2) = -(J(3,3) * J(1,2) - J(3,2) * J(1,3))
            Jinv(1,3) =   J(2,3) * J(1,2) - J(2,2) * J(1,3)

            Jinv(2,1) = -(J(3,3) * J(2,1) - J(3,1) * J(2,3))
            Jinv(2,2) =   J(3,3) * J(1,1) - J(3,1) * J(1,3)
            Jinv(2,3) = -(J(2,3) * J(1,1) - J(2,1) * J(1,3))

            Jinv(3,1) =   J(3,2) * J(2,1) - J(3,1) * J(2,2)
            Jinv(3,2) = -(J(3,2) * J(1,1) - J(3,1) * J(1,2))
            Jinv(3,3) =   J(2,2) * J(1,1) - J(2,1) * J(1,2)

            Jinv = Jinv * det

            pv(1) = pv(1) - dot_product(Jinv(1,:), F(:))
            pv(2) = pv(2) - dot_product(Jinv(2,:), F(:))
            pv(3) = pv(3) - dot_product(Jinv(3,:), F(:))
         end do 
      end associate
   
      return
   end subroutine swiftest_gr_vel2pseudovel


   pure module subroutine swiftest_gr_vh2pv_body(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from heliocentric velocity to pseudovelocity for Swiftest bodies
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self  
         !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param 
         !! Current run configuration parameters 
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(:,:), allocatable        :: pv    
         !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(n => self%nbody)
         if (n == 0) return
         allocate(pv, mold = self%vh)
         do i = 1, n
            call swiftest_gr_vel2pseudovel(param, self%mu(i), self%rh(:, i), self%vh(:, i), pv(:, i))
         end do
         call move_alloc(pv, self%vh)
      end associate

      return
   end subroutine swiftest_gr_vh2pv_body

end submodule s_swiftest_gr