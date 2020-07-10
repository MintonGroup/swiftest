submodule(swiftest_classes) gr_implementations
contains

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
            xh(:) = xb(:, i) - xbsun(:)
            vh(:) = vb(:, i) - vbsun(:)
            rmag = .mag. xh 
            vmag2 = vh .dot. vh
            rdotv = xh .dot. vh
            agr(:, i) =  mu * c2 / rmag**3 * ((4 * mu(i) / rmag - vmag2) * xh(:) + 4 * rdotv * vh(:))
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