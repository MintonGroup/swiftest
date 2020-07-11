submodule(whm_classes) whm_gr_implementations
contains
   module procedure whm_gr_getacch_pl 
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of massive bodies
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_getacch.f90
      use swiftest
      implicit none

      integer(I4B)                   :: i
      real(DP), dimension(NDIM)      :: suma
      real(DP), dimension(:, :), allocatable :: aj
      real(DP), dimension(:), allocatable :: rjmag4, beta
      
      associate(n => self%nbody, msun => cb%Gmass, mu => self%muj, c2 => config%inv_c2, &
         ah => self%ah, xj => self%xj, GMpl => self%Gmass, eta => self%eta)
         if (n == 0) return
         allocate(rjmag4(n))
         allocate(beta(n))
         allocate(aj, mold = ah)
         rjmag4(:) = (xj(:, 1:n) .dot. xj(:, 1:n))**2
         beta(:)   = - mu(1:n)**2 * c2 
         do concurrent (i = 1:n)
            aj(:, i) = beta(i) * xj(:, i) / rjmag4(i)
         end do
         suma(:) = 0.0_DP
         ah(:, 1) = ah(:, 1) + aj(:, 1)
         do i = 2, n
            suma(:) = suma(:) + GMpl(i) * aj(:, i) / eta(i)
            ah(:, i) = ah(:, i) + aj(:, i) + suma(:)
         end do
      end associate
      return
   end procedure whm_gr_getacch_pl

   module procedure whm_gr_getacch_tp
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of test particles
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_getacch.f90
      use swiftest
      implicit none

      integer(I4B)                   :: i
      real(DP) :: rjmag4, beta
      
      associate(n => self%nbody, msun => cb%Gmass, mu => self%mu,& 
         c2 => config%inv_c2, ah => self%ah, xh => self%xh, status => self%status)
         if (n == 0) return
         do concurrent (i = 1:n, status(i) == active)
            rjmag4 = (xh(i, 1:n) .dot. xh(i, 1:n))**2
            beta = - mu(i)**2 * c2 
            ah(:, i) = ah(:, i) + beta * xh(:, i) / rjmag4
         end do
      end associate
      return
   end procedure whm_gr_getacch_tp

   module procedure whm_gr_p4_pl
      !! author: David A. Minton
      !!
      !! Position kick to massive bodies due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      use swiftest
      implicit none
      integer(I4B) :: i

      associate(n => self%nbody, xj => self%xj, vj => self%vj, status => self%status, c2 => config%inv_c2)
         if (n == 0) return
         do concurrent (i = 1:n, status(i) == ACTIVE)
            call p4_func(xj(:, i), vj(:, i), dt, c2)
         end do
      end associate
 
     return
   end procedure whm_gr_p4_pl

   module procedure whm_gr_p4_tp
      !! author: David A. Minton
      !!
      !! Position kick to test particles due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      use swiftest
      implicit none
      integer(I4B) :: i

      associate(n => self%nbody, xh => self%xh, vh => self%vh, status => self%status, c2 => config%inv_c2)
         if (n == 0) return
         do concurrent (i = 1:n, status(i) == ACTIVE)
            call p4_func(xh(:, i), vh(:, i), dt, c2)
         end do
      end associate
 
     return
   end procedure whm_gr_p4_tp

   module procedure whm_gr_pv2vh_pl
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from pseudovelocity to heliocentric velocity for massive bodies
      !! in a WHM object
      use swiftest
      implicit none

      integer(I4B) :: i
      real(DP), dimension(:), allocatable :: mu
      real(DP), dimension(:,:), allocatable :: vh !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(n => self%nbody, xh => self%xh, pv => self%vh, status => self%status, mu => self%muj)
         if (n == 0) return
         allocate(vh, mold = pv)
         do concurrent(i = 1:n, status(i) == ACTIVE)
            call gr_pseudovel2vel(config, mu(i), xh(:, i), pv(:, i), vh(:, i))
            pv(:, i) = vh(:, i)
         end do
         deallocate(vh)
      end associate
      return
   end procedure whm_gr_pv2vh_pl

   module procedure whm_gr_pv2vh_tp
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from pseudovelocity to heliocentric velocity for test particles bodies
      !! in a WHM object
      use swiftest
      implicit none

      integer(I4B) :: i
      real(DP), dimension(:), allocatable :: mu
      real(DP), dimension(:,:), allocatable :: vh !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(n => self%nbody, xh => self%xh, pv => self%vh, status => self%status, mu => self%mu)
         if (n == 0) return
         allocate(vh, mold = pv)
         do concurrent(i = 1:n, status(i) == ACTIVE)
            call gr_pseudovel2vel(config, mu(i), xh(:, i), pv(:, i), vh(:, i))
            pv(:, i) = vh(:, i)
         end do
         deallocate(vh)
      end associate
      return
   end procedure whm_gr_pv2vh_tp

   module procedure whm_gr_vh2pv_pl
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from heliocentric velocity to pseudovelocity for massive bodies
      !! in a WHM object
      use swiftest
      implicit none
      integer(I4B) :: i
      real(DP), dimension(:), allocatable :: mu
      real(DP), dimension(:,:), allocatable :: pv !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(n => self%nbody, xh => self%xh, vh => self%vh, status => self%status, mu => self%muj)
         if (n == 0) return
         allocate(pv, mold = vh)
         do concurrent(i = 1:n, status(i) == ACTIVE)
            call gr_vel2pseudovel(config, mu(i), xh(:, i), vh(:, i), pv(:, i))
            vh(:, i) = pv(:, i)
         end do
         deallocate(pv)
      end associate
      return
   end procedure whm_gr_vh2pv_pl

   module procedure whm_gr_vh2pv_tp
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from heliocentric velocity to pseudovelocity for teset particles
      !! in a WHM object
      use swiftest
      implicit none
      integer(I4B) :: i
      real(DP), dimension(:), allocatable :: mu
      real(DP), dimension(:,:), allocatable :: pv !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(n => self%nbody, xh => self%xh, vh => self%vh, status => self%status, mu => self%mu)
         if (n == 0) return
         allocate(pv, mold = vh)
         do concurrent(i = 1:n, status(i) == ACTIVE)
            call gr_vel2pseudovel(config, mu(i), xh(:, i), vh(:, i), pv(:, i))
            vh(:, i) = pv(:, i)
         end do
         deallocate(pv)
      end associate
      return
   end procedure whm_gr_vh2pv_tp

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
            !do concurrent (k = 1:NDIM)
            !   do concurrent (i = 1:NDIM)
            do k = 1, NDIM
               do i = 1, NDIM
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
            end do
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

   pure subroutine p4_func(x, v, dt, c2)
      use swiftest
      implicit none
      real(DP), dimension(:), intent(inout) :: x
      real(DP), dimension(:), intent(in)    :: v
      real(DP),               intent(in)    :: dt, c2 
      real(DP), dimension(NDIM)             :: dr
      real(DP)                              :: vmag2
      integer(I4B) :: i

      vmag2 = v .dot. v 
      dr(:) = -c2 * vmag2 * v(:)
      x(:) = x(:) + dr(:) * dt

      return
   end subroutine p4_func  

end submodule whm_gr_implementations