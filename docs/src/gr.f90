submodule(swiftest_classes) s_gr
   use swiftest
contains

   module pure subroutine gr_kick_getaccb_ns_body(self, system, param) 
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
      class(swiftest_body),         intent(inout) :: self   !! Swiftest generic body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      ! Internals
      real(DP)                  :: rmag, rdotv, vmag2
      integer(I4B)              :: i

      associate(n => self%nbody, cb => system%cb, inv_c2 => param%inv_c2)
         if (n == 0) return
         do i = 1, n
            rmag = norm2(self%xh(:,i))
            vmag2 = dot_product(self%vh(:,i), self%vh(:,i))
            rdotv = dot_product(self%xh(:,i), self%vh(:,i))
            self%agr(:, i) = self%mu * inv_c2 / rmag**3 * ((4 * self%mu(i) / rmag - vmag2) * self%xh(:,i) + 4 * rdotv * self%vh(:,i))
         end do

         select type(self)
         class is (swiftest_pl)
            do i = 1, NDIM
               cb%agr(i) = -sum(self%Gmass(1:n) * self%agr(1:n, i) / cb%Gmass)
            end do
         end select
      end associate 

      return
   end subroutine gr_kick_getaccb_ns_body


   module subroutine gr_kick_getacch(mu, x, lmask, n, inv_c2, agr) 
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of massive bodies
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_kick_getacch.f90
      implicit none
      ! Arguments
      real(DP), dimension(:),     intent(in)    :: mu     !! Gravitational constant
      real(DP), dimension(:,:),   intent(in)    :: x      !! Position vectors
      logical,  dimension(:),     intent(in)    :: lmask  !! Logical mask indicating which bodies to compute
      integer(I4B),               intent(in)    :: n      !! Total number of bodies
      real(DP),                   intent(in)    :: inv_c2 !! Inverse speed of light squared: 1 / c**2
      real(DP), dimension(:,:),   intent(out)   :: agr    !! Accelerations
      ! Internals
      integer(I4B)                              :: i
      real(DP)                                  :: beta, rjmag4
     
      agr(:,:) = 0.0_DP
      do concurrent (i = 1:n, lmask(i))
         rjmag4 = (dot_product(x(:, i), x(:, i)))**2
         beta = -mu(i)**2 * inv_c2 
         agr(:, i) = 2 * beta * x(:, i) / rjmag4
      end do

      return
   end subroutine gr_kick_getacch


   module pure subroutine gr_p4_pos_kick(param, x, v, dt)
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
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      real(DP), dimension(:),     intent(inout) :: x     !! Position vector
      real(DP), dimension(:),     intent(in)    :: v     !! Velocity vector
      real(DP),                   intent(in)    :: dt    !! Step size
      ! Internals
      real(DP), dimension(NDIM)             :: dr
      real(DP)                              :: vmag2

      vmag2 = dot_product(v(:), v(:)) 
      dr(:) = -2 * param%inv_c2 * vmag2 * v(:)
      x(:) = x(:) + dr(:) * dt

      return
   end subroutine gr_p4_pos_kick


   module pure subroutine gr_pseudovel2vel(param, mu, xh, pv, vh) 
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
      class(swiftest_parameters), intent(in)  :: param !! Current run configuration parameters 
      real(DP),                   intent(in)  :: mu    !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
      real(DP), dimension(:),     intent(in)  :: xh    !! Heliocentric position vector 
      real(DP), dimension(:),     intent(in)  :: pv    !! Pseudovelocity velocity vector - see Saha & Tremain (1994), eq. (32)
      real(DP), dimension(:),     intent(out) :: vh    !! Heliocentric velocity vector 
      ! Internals
      real(DP) :: vmag2, rmag, grterm
   
      associate(inv_c2 => param%inv_c2)
         vmag2 = dot_product(pv(:), pv(:))
         rmag  = norm2(xh(:))
         grterm = 1.0_DP - inv_c2 * (0.5_DP * vmag2 + 3 * mu / rmag)
         vh(:) = pv(:) * grterm
      end associate

      return
   end subroutine gr_pseudovel2vel


   module pure subroutine gr_pv2vh_body(self, param)
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
            call gr_pseudovel2vel(param, self%mu(i), self%xh(:, i), self%vh(:, i), vh(:, i))
         end do
         call move_alloc(vh, self%vh)
      end associate

      return
   end subroutine gr_pv2vh_body


   module pure subroutine gr_vel2pseudovel(param, mu, xh, vh, pv)
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
      class(swiftest_parameters), intent(in)  :: param !! Current run configuration parameters 
      real(DP),                   intent(in)  :: mu    !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
      real(DP), dimension(:),     intent(in)  :: xh    !! Heliocentric position vector 
      real(DP), dimension(:),     intent(in)  :: vh    !! Heliocentric velocity vector 
      real(DP), dimension(:),     intent(out) :: pv    !! Pseudovelocity vector - see Saha & Tremain (1994), eq. (32)
      ! Internals
      real(DP)                       :: v2, G, pv2, rterm, det
      real(DP), dimension(NDIM,NDIM) :: J,Jinv
      real(DP), dimension(NDIM)      :: F
      integer(I4B)                   :: n,i,k
      integer(I4B), parameter        :: MAXITER = 50
      real(DP),parameter             :: TOL = 1.0e-12_DP

      associate(inv_c2 => param%inv_c2)
         pv(1:NDIM) = vh(1:NDIM) ! Initial guess
         rterm = 3 * mu / norm2(xh(:))
         v2 = dot_product(vh(:), vh(:))
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

            do i = 1, NDIM
               pv(i) = pv(i) - dot_product(Jinv(i,:), F(:))
            end do
         end do 
      end associate
   
      return
   end subroutine gr_vel2pseudovel


   module pure subroutine gr_vh2pv_body(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from heliocentric velocity to pseudovelocity for Swiftest bodies
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self  !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(:,:), allocatable        :: pv !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(n => self%nbody)
         if (n == 0) return
         allocate(pv, mold = self%vh)
         do i = 1, n
            call gr_vel2pseudovel(param, self%mu(i), self%xh(:, i), self%vh(:, i), pv(:, i))
         end do
         call move_alloc(pv, self%vh)
      end associate

      return
   end subroutine gr_vh2pv_body

end submodule s_gr