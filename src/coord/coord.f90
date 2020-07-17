submodule (swiftest_classes) coord
contains
   module procedure coord_h2b_body
      !! author: David A. Minton
      !!
      !! Convert from heliocentric to barycentric coordinates (position and velocity)
      !! Polymorphic method that accepts either test particles or massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b.f90 and coord_h2b_tp.f90
      !! Adapted from Hal Levison's Swift routine coord_h2b.f and coord_h2b_tp.f
      use swiftest
      integer(I4B)  :: i
      real(DP)      :: msys

      associate(n => self%nbody, xbcb => cb%xb, vbcb => cb%vb, status => self%status, Mcb => cb%Gmass, &
         xb => self%xb, xh => self%xh, vb => self%vb, vh => self%vh)

         select type(self)
         class is (swiftest_pl)
            associate(Mpl => self%Gmass)
               msys = Mcb + sum(Mpl(1:n), status(1:n) == ACTIVE)
               do i = 1, NDIM
                  xbcb(i) = -sum(Mpl(1:n) * xh(i, 1:n), status(1:n) == ACTIVE) / msys
                  vbcb(i) = -sum(Mpl(1:n) * vh(i, 1:n), status(1:n) == ACTIVE) / msys
               end do
            end associate
         end select

         do concurrent(i = 1:n, status(i) == ACTIVE) !shared(n, status, xb, xh, xbcb, vb, vh, xbcb, vbcb)
            xb(:, i) = xh(:, i) + xbcb(:)
            vb(:, i) = vh(:, i) + vbcb(:)
         end do
      end associate

      return
   end procedure coord_h2b_body

   module procedure coord_b2h_body
      !! author: David A. Minton
      !!
      !! Convert from barycentric to heliocentric coordinates (position and velocity)
      !! Polymorphic method that accepts either test particles or massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b.f90 and coord_h2b_tp.f90
      !! Adapted from Hal Levison's Swift routine coord_h2b.f and coord_h2b_tp.f
      use swiftest
      integer(I4B)          :: i

      associate(n => self%nbody, xbcb => cb%xb, vbcb => cb%vb, status => self%status, & 
         xb => self%xb, xh => self%xh, vb => self%vb, vh => self%vh)
         do concurrent(i = 1:n, status(i) == ACTIVE) !shared(n, status, xb, xbh, vb, vh, xbcb, vbcb)
            xh(:, i) = xb(:, i) - xbcb(:)
            vh(:, i) = vb(:, i) - vbcb(:)
         end do
      end associate

      return
   end procedure coord_b2h_body

   module procedure coord_vb2vh_body
      !! author: David A. Minton
      !!
      !! Convert from barycentric to heliocentric coordinates (velocity only)
      !! Polymorphic method that accepts either test particles or massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vb2vh.f90 and coord_vb2vh_tp.f90
      !! Adapted from Hal Levison's Swift routine coord_vb2vh.f and coord_vb2vh_tp.f
      use swiftest
      implicit none
      integer(I4B)              :: i
      real(DP), dimension(NDIM) :: vtmp
   
      vtmp(:) = 0.0_DP
      associate(n => self%nbody, vbcb => cb%vb, status => self%status, & 
         xh => self%xh, vb => self%vb, vh => self%vh, Mcb => cb%Gmass)
         select type(self)
         class is (swiftest_pl)
            associate(Mpl => self%Gmass)
               do i = 1, NDIM
                  vbcb(i) = -sum(Mpl(1:n) * vb(i, 1:n), status(1:n) == ACTIVE) / Mcb
               end do
            end associate
         end select
         do concurrent(i = 1:n, status(i) == ACTIVE) !shared(n, status, vh, vb, vbcb)
            vh(:, i) = vb(:, i) - vbcb(:)
         end do
      end associate
   
      return
   end procedure coord_vb2vh_body 

   module procedure coord_vh2vb_body
      !! author: David A. Minton
      !!
      !! Convert from heliocentric to barycentric coordinates (velocity only)
      !! Polymorphic method that accepts either test particles or massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b.f90 and coord_h2b_tp.f90
      !! Adapted from Hal Levison's Swift routine coord_h2b.f and coord_h2b_tp.f
      use swiftest
      integer(I4B)  :: i
      real(DP)      :: msys

      associate(n => self%nbody, vbcb => cb%vb, status => self%status, & 
         vb => self%vb, vh => self%vh, Mcb => cb%Gmass)
         select type(self)
         class is (swiftest_pl)
            associate(Mpl => self%Gmass)
               msys = Mcb + sum(Mpl(1:n))
               do i = 1, NDIM
                  vbcb(i) = -sum(Mpl(1:n) * vh(i, 1:n), status(1:n) == ACTIVE) / msys
               end do
            end associate
         end select
         do concurrent(i = 1:n, status(i) == ACTIVE) !shared(n, status, vb, vh, vbcb)
            vb(:, i) = vh(:, i) + vbcb(:)
         end do
      end associate

      return
   end procedure coord_vh2vb_body
end submodule coord
