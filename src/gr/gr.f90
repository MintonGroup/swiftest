submodule(swiftest_classes) gr_implementations
contains

   module procedure gr_pv2vh_body 
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from pseudovelocity to heliocentric velocity for all bodies
      !! in a Swiftest body object
      use swiftest
      implicit none

      integer(I4B) :: i
      real(DP), dimension(:), allocatable :: mu
      
      associate(n => self%nbody, xh => self%xh, pv => self%vh, status => self%status)
         if (n == 0) return
         select type(self)
         class is (whm_pl)
            allocate(mu, source = self%muj)
         class is (whm_tp)
            allocate(mu, source = self%mu)
         end select 
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
      real(DP), dimension(:), allocatable :: mu
      
      associate(n => self%nbody, xh => self%xh, vh => self%vh, status => self%status)
         if (n == 0) return
         select type(self)
         class is (whm_pl)
            allocate(mu, source = self%muj)
         class is (whm_tp)
            allocate(mu, source = self%mu)
         end select 
         do concurrent(i = 1:n, status(i) == ACTIVE)
            call gr_vel2pseudovel(config, mu(i), xh(i, :), vh(i, :), pv(i, :))
         end do
      end associate

   end procedure gr_vh2pv_body


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