submodule(swiftest_classes) gr_implementations
contains
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
         select type(self)
         class is (whm_pl)
            associate(xj => self%xj, vj => self%vj)
               call p4_func(xj(1:n, :), vj(1:n, :))
            end associate
         class is (whm_tp)
            associate(xh => self%xh, vh => self%vh)
               call p4_func(xh(1:n, :), vh(1:n, :))
            end associate
         end select
      end associate
 
     return
     contains
         pure subroutine p4_func(x,v)
            implicit none
            real(DP), dimension(:,:), intent(inout) :: x
            real(DP), dimension(:,:), intent(in)    :: v
            real(DP), dimension(:,:), allocatable :: dr, vjmag2
            real(DP), dimension(:),   allocatable :: vmag2
            integer(I4B) :: i

            vmag2(:) = v(:, :) .dot. v(:, :)
            do concurrent (i = 1:NDIM)
               dr(:, i) = - config%inv_c2 * vmag2(:) * v(:,i)
               x(:, i)  =   x(:, i) + dr(:, i) * dt
            end do

            return
         end subroutine p4_func  

   end procedure gr_p4_body

   module procedure gr_vh2pv_body !(self, config, dt)
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from heliocentric velocity to pseudovelocity for all bodies
      !! in a Swiftest body object
      implicit none
      integer(I4B) :: i
      
      associate(n => self%nbody, xh => self%xh, vh => self%vh, mu => self%mu_vec, status => self%status)
         do concurrent(i = 1:n, status(i) == ACTIVE)
            call gr_vel2pseudovel(mu(i), xh(i, :), vh(i, :))
         end do
      end associate

      contains
         pure subroutine gr_vel2pseudovel(mu, xh, vh) !(xh, pvh, mu, c2, vh)
            !! author: David A. Minton
            !!
            !! Converts the heliocentric velocity into a pseudovelocity with relativistic corrections. 
            !! Uses Newton-Raphson method with direct inversion of the Jacobian (yeah, it's slow, but 
            !! this is only done once per run).
            !!
            !! Adapted from David A. Minton's Swifter routine gr_vel2pseudovel.f90 
            use swiftest
            implicit none
            real(DP),               intent(in)    :: mu
            real(DP), dimension(:), intent(in)    :: xh
            real(DP), dimension(:), intent(inout) :: vh

            real(DP) :: vmag2, rmag, grterm
         
            real(DP)                       :: v2, G, pv2, rterm, det
            real(DP), dimension(NDIM,NDIM) :: J,Jinv
            real(DP), dimension(NDIM)      :: F
            integer(I4B)                   :: n,i,k
            integer(I4B), parameter        :: MAXITER = 50
            real(DP),parameter             :: TOL = 1.0e-12_DP
            real(DP), dimension(NDIM)      :: pvh
      
            associate (c2 => config%inv_c2)
               pvh(1:NDIM) = vh(1:NDIM) ! Initial guess
               rterm = 3 * mu / (.mag. xh(:)) 
               v2 = vh(:) .dot. vh(:)
               do n = 1, MAXITER
                  pv2 = pvh(:) .dot. pvh(:) 
                  G = 1.0_DP - c2 * (0.5_DP * pv2 + rterm)
                  F(:) = pvh(:) * G - vh(:)
                  if (abs(sum(F) / v2 ) < TOL) exit ! Root found
         
                  ! Calculate the Jacobian
                  do concurrent (k = 1:NDIM)
                     do concurrent (i = 1:NDIM)
                        if (i == k) then
                           J(i,k) = G - c2 * pvh(k)
                        else
                           J(i,k) = - c2 * pvh(k)
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
         
                  do concurrent (i = 1:NDIM)
                     pvh(i) = pvh(i) - Jinv(:, i) .dot. F(:) 
                  end do
               end do 
         
               vh(:) = pvh(:)
            end associate
            return
         end subroutine gr_vel2pseudovel
   end procedure gr_vh2pv_body

   module procedure gr_pv2vh_body !(self, config, dt)
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from pseudovelocity to heliocentric velocity for all bodies
      !! in a Swiftest body object
      use swiftest
      implicit none

      integer(I4B) :: i
      
      associate(n => self%nbody, xh => self%xh, vh => self%vh, mu => self%mu_vec, status => self%status)
         do concurrent(i = 1:n, status(i) == ACTIVE)
            call gr_pseudovel2vel(mu(i), xh(i, :), vh(i, :))
         end do
      end associate

      contains
         pure subroutine gr_pseudovel2vel(mu, xh, vh) !(xh, pvh, mu, c2, vh)
            !! author: David A. Minton
            !!
            !! Converts the relativistic pseudovelocity back into a veliocentric velocity
            !!    Based on Saha & Tremaine (1994) Eq. 32
            !!
            !! Adapted from David A. Minton's Swifter routine gr_pseudovel2vel.f90 
            use swiftest
            implicit none

            real(DP),               intent(in) :: mu
            real(DP), dimension(:), intent(in) :: xh 
            real(DP), dimension(:), intent(inout) :: vh
            real(DP) :: vmag2, rmag, grterm
        
            associate(c2 => config%inv_c2)
               vmag2 = vh(:) .dot. vh(:) 
               rmag  = .mag. xh(:) 
               grterm = 1.0_DP - c2 * (0.5_DP * vmag2 + 3 * mu / rmag)
               vh(:) = vh(:) * grterm
            end associate
            return
         end subroutine gr_pseudovel2vel

   end procedure gr_pv2vh_body

end submodule gr_implementations