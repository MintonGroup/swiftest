submodule(swiftest_classes) gr_implementations
contains
   module procedure gr_getacch_body ! (self, cb, config
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of planet
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_getacch.f90
      use swiftest
      implicit none

      integer(I4B)                   :: i
      real(DP), dimension(NDIM)      :: suma
      real(DP), dimension(:, :), allocatable :: aj
      real(DP), dimension(:), allocatable :: rjmag4, beta
      

      associate(n => self%nbody, msun => cb%Gmass, mu => self%mu, c2 => config%inv_c2, ah => self%ah)
         if (n == 0) return
         allocate(rjmag4(n))
         allocate(beta(n))
         select type(self)
         class is (whm_pl)
            allocate(aj, mold = ah)
            rjmag4(:) = (self%xj(1:n, :) .dot. self%xj(1:n, :))**2
            beta(:)   = - mu(1:n)**2 * c2 
            do concurrent (i = 1:NDIM)
               aj(:, i) = beta(1:n) / rjmag4(1:n) * self%xj(1:n, i)
            end do
            suma(:) = 0.0_DP
            ah(1, :) = ah(1, :) + aj(2, :)
            do i = 2, n
               suma(:) = suma(:) + self%Gmass(i) * aj(i, :)
               ah(i, :) = ah(i, :) + aj(i, :) + suma(:)
            end do
         class is (swiftest_tp)
            rjmag4(:) = (self%xh(1:n, :) .dot. self%xh(1:n, :))**2
            beta(:) = - mu(1:n)**2 * c2 
            ah(1, :) = ah(1, :) + aj(2, :)
            do i = 2, n
               ah(i, :) = ah(i, :) + aj(i, :) 
            end do
         end select
      end associate
      return
   end procedure gr_getacch_body

   module procedure gr_pv2vh_body 
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from pseudovelocity to heliocentric velocity for all bodies
      !! in a Swiftest body object
      use swiftest
      implicit none

      integer(I4B) :: i
      
      associate(n => self%nbody, xh => self%xh, pv => self%vh, mu => self%mu, status => self%status)
         if (n == 0) return
         do concurrent(i = 1:n, status(i) == ACTIVE)
            call gr_pseudovel2vel(config, mu(i), xh(i, :), pv(i, :), vh(i, :))
         end do
      end associate
      return
   end procedure gr_pv2vh_body

   module procedure gr_vh2pv_body 
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from heliocentric velocity to pseudovelocity for all bodies
      !! in a Swiftest body object
      use swiftest
      implicit none
      integer(I4B) :: i
      
      associate(n => self%nbody, xh => self%xh, vh => self%vh, mu => self%mu, status => self%status)
         if (n == 0) return
         do concurrent(i = 1:n, status(i) == ACTIVE)
            call gr_vel2pseudovel(config, mu(i), xh(i, :), vh(i, :), pv(i, :))
         end do
      end associate

   end procedure gr_vh2pv_body

   module procedure gr_p4_body
      !! author: David A. Minton
      !!
      !! Position kick due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      use swiftest
      implicit none

      associate(n => self%nbody)
         if (n == 0) return
         select type(self)
         class is (whm_pl)
            associate(xj => self%xj, vj => self%vj)
               call p4_func(xj(1:n, :), vj(1:n, :), dt, config%inv_c2)
            end associate
         class is (swiftest_tp)
            associate(xh => self%xh, vh => self%vh)
               call p4_func(xh(1:n, :), vh(1:n, :), dt, config%inv_c2)
            end associate
         end select
      end associate
 
     return
   end procedure gr_p4_body

   pure subroutine p4_func(x, v, dt, c2)
      use swiftest
      implicit none
      real(DP), dimension(:,:), intent(inout) :: x
      real(DP), dimension(:,:), intent(in)    :: v
      real(DP),                 intent(in)    :: dt, c2 
      real(DP), dimension(:,:), allocatable :: dr
      real(DP), dimension(:),   allocatable :: vmag2
      integer(I4B) :: i

      allocate(vmag2, mold = v(:,1))
      allocate(dr, mold = x)
      vmag2(:) = v(:, :) .dot. v(:, :)
      do concurrent (i = 1:NDIM)
         dr(:, i) = - c2 * vmag2(:) * v(:,i)
         x(:, i)  =   x(:, i) + dr(:, i) * dt
      end do

      return
   end subroutine p4_func  

   pure subroutine gr_vel2pseudovel(config, mu, xh, vh, pv)
      !! author: David A. Minton
      !!
      !! Converts the heliocentric velocity into a pseudovelocity with relativistic corrections. 
      !! Uses Newton-Raphson method with direct inversion of the Jacobian (yeah, it's slow, but 
      !! this is only done once per run).
      !!
      !! Adapted from David A. Minton's Swifter routine gr_vel2pseudovel.f90
      use swiftest
      implicit none
      class(swiftest_configuration), intent(in)  :: config !! Input collection of user-defined configuration parameters 
      real(DP),                      intent(in)  :: mu     !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
      real(DP), dimension(:),        intent(in)  :: xh     !! Heliocentric position vector 
      real(DP), dimension(:),        intent(in)  :: vh     !! Heliocentric velocity vector 
      real(DP), dimension(:),        intent(out) :: pv     !! Pseudovelocity vector - see Saha & Tremain (1994), eq. (32)

      real(DP)                       :: v2, G, pv2, rterm, det
      real(DP), dimension(NDIM,NDIM) :: J,Jinv
      real(DP), dimension(NDIM)      :: F
      integer(I4B)                   :: n,i,k
      integer(I4B), parameter        :: MAXITER = 50
      real(DP),parameter             :: TOL = 1.0e-12_DP

      associate (c2 => config%inv_c2)
         pv(1:NDIM) = vh(1:NDIM) ! Initial guess
         rterm = 3 * mu / (.mag. xh(:)) 
         v2 = vh(:) .dot. vh(:)
         do n = 1, MAXITER
            pv2 = pv(:) .dot. pv(:) 
            G = 1.0_DP - c2 * (0.5_DP * pv2 + rterm)
            F(:) = pv(:) * G - vh(:)
            if (abs(sum(F) / v2 ) < TOL) exit ! Root found
   
            ! Calculate the Jacobian
            do concurrent (k = 1:NDIM)
               do concurrent (i = 1:NDIM)
                  if (i == k) then
                     J(i,k) = G - c2 * pv(k)
                  else
                     J(i,k) = - c2 * pv(k)
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
               pv(i) = pv(i) - dot_product(Jinv(i,:) ,F(:))
            end DO
         end do 
   
      end associate
      return
   end subroutine gr_vel2pseudovel

   pure subroutine gr_pseudovel2vel(config, mu, xh, pv, vh) 
      !! author: David A. Minton
      !!
      !! Converts the relativistic pseudovelocity back into a veliocentric velocity
      !!    Based on Saha & Tremaine (1994) Eq. 32
      !!
      !! Adapted from David A. Minton's Swifter routine gr_pseudovel2vel.f90 
      use swiftest
      implicit none
      class(swiftest_configuration), intent(in)  :: config !! Input collection of user-defined configuration parameters 
      real(DP),                      intent(in)  :: mu     !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body
      real(DP), dimension(:),        intent(in)  :: xh     !! Heliocentric position vector 
      real(DP), dimension(:),        intent(in)  :: pv     !! Pseudovelocity velocity vector - see Saha & Tremain (1994), eq. (32)
      real(DP), dimension(:),        intent(out) :: vh     !! Heliocentric velocity vector 

      real(DP) :: vmag2, rmag, grterm
   
      associate(c2 => config%inv_c2)
         vmag2 = pv(:) .dot. pv(:) 
         rmag  = .mag. xh(:) 
         grterm = 1.0_DP - c2 * (0.5_DP * vmag2 + 3 * mu / rmag)
         vh(:) = pv(:) * grterm
      end associate
      return
   end subroutine gr_pseudovel2vel

   module procedure gr_getaccb_ns_body
      !! author: David A. Minton
      !!
      !! Add relativistic correction acceleration for non-symplectic integrators
      !!    Based on Quinn, Tremaine, & Duncan 1998
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_getaccb_ns.f90
      use swiftest
      implicit none

      real(DP), dimension(NDIM) :: xh, vh
      real(DP)                  :: rmag, rdotv, vmag2
      integer(I4B)              :: i

      associate(n => self%nbody, msun => cb%Gmass, vbsun => cb%vb, xbsun => cb%xb, mu => self%mu, c2 => config%inv_c2, &
                xb => self%xb, vb => self%vb)
         if (n == 0) return
         do i = 1, n
            xh(:) = xb(i, :) - xbsun(:)
            vh(:) = vb(i, :) - vbsun(:)
            rmag = .mag. xh 
            vmag2 = vh .dot. vh
            rdotv = xh .dot. vh
            agr(i, :) =  mu * c2 / rmag**3 * ((4 * mu(i) / rmag - vmag2) * xh(:) + 4 * rdotv * vh(:))
         end do

         agr0 =  0.0_DP
         select type(self)
         class is (swiftest_pl)
            do concurrent(i = 1:NDIM)
               agr0(i) = -sum(self%Gmass(1:n) * agr(1:n, i) / msun)
            end do
         end select

      end associate 

      return

   end procedure gr_getaccb_ns_body

end submodule gr_implementations