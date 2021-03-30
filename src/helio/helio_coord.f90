submodule (helio_classes) s_helio_coord
contains

   module subroutine helio_coord_h2b_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from heliocentric to barycentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b.f90 
      !! Adapted from Hal Levison's Swift routine coord_h2b.f 
      use swiftest
      implicit none
      ! Arguments
      class(helio_pl),    intent(inout) :: self !! Helio massive body object
      class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)  :: i
      real(DP)      :: msys

      associate(npl => self%nbody, xbcb => cb%xb, vbcb => cb%vb, Mcb => cb%Gmass, &
         xb => self%xb, xh => self%xh, vb => self%vb, vh => self%vh, Mpl => self%Gmass)

         msys = Mcb + sum(Mpl(1:npl))
         do i = 1, NDIM
            xbcb(i) = -sum(Mpl(1:npl) * xh(i, 1:npl)) / Mcb
            vbcb(i) = -sum(Mpl(1:npl) * vh(i, 1:npl)) / Mcb
            xb(i, 1:npl) = xh(i, 1:npl) + xbcb(i)
            vb(i, 1:npl) = vh(i, 1:npl) + vbcb(i)
         end do
      end associate

      return
   end subroutine helio_coord_h2b_pl

   module subroutine helio_coord_h2b_tp(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from heliocentric to barycentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b_tp.f90 
      !! Adapted from Hal Levison's Swift routine coord_h2b_tp.f 
      use swiftest
      implicit none
      ! Arguments
      class(helio_tp),    intent(inout) :: self !! Helio test particle object
      class(swiftest_cb), intent(in) :: cb   !! Swiftest central body object

      associate(ntp => self%nbody, xbcb => cb%xb, vbcb => cb%vb, status => self%status, &
               xb => self%xb, xh => self%xh, vb => self%vb, vh => self%vh)

         where(status(1:ntp) == ACTIVE)
            xb(1, 1:ntp) = xh(1, 1:ntp) + xbcb(1)
            xb(2, 1:ntp) = xh(2, 1:ntp) + xbcb(2)
            xb(3, 1:ntp) = xh(3, 1:ntp) + xbcb(3)

            vb(1, 1:ntp) = vh(1, 1:ntp) + vbcb(1)
            vb(2, 1:ntp) = vh(2, 1:ntp) + vbcb(2)
            vb(3, 1:ntp) = vh(3, 1:ntp) + vbcb(3)
         end where
      end associate

      return
   end subroutine helio_coord_h2b_tp

   module subroutine helio_coord_b2h_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from barycentric to heliocentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_b2h.f90 
      !! Adapted from Hal Levison's Swift routine coord_b2h.f 
      use swiftest
      implicit none
      ! Arguments
      class(helio_pl),     intent(inout) :: self !! Helio massive body object
      class(swiftest_cb),  intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)          :: i

      associate(npl => self%nbody, xbcb => cb%xb, vbcb => cb%vb, xb => self%xb, xh => self%xh, &
                vb => self%vb, vh => self%vh)
         do i = 1, NDIM
            xh(i, 1:npl) = xb(i, 1:npl) - xbcb(i)
            vh(i, 1:npl) = vb(i, 1:npl) - vbcb(i)
         end do
      end associate

      return
   end subroutine helio_coord_b2h_pl

   module subroutine helio_coord_b2h_tp(self, cb)
      !! author: David A. Minton
      !!
      !! Convert test particles from barycentric to heliocentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_b2h_tp.f90 
      !! Adapted from Hal Levison's Swift routine coord_b2h_tp.f 
      use swiftest
      implicit none
      ! Arguments
      class(helio_tp),     intent(inout) :: self !! Helio massive body object
      class(swiftest_cb),  intent(in)    :: cb   !! Swiftest central body object

      associate(ntp => self%nbody, xbcb => cb%xb, vbcb => cb%vb, xb => self%xb, xh => self%xh, &
                vb => self%vb, vh => self%vh, status => self%status)
         where(status(1:ntp) == ACTIVE)
            xh(1, 1:ntp) = xb(1, 1:ntp) - xbcb(1)
            xh(2, 1:ntp) = xb(2, 1:ntp) - xbcb(2)
            xh(3, 1:ntp) = xb(3, 1:ntp) - xbcb(3)

            vh(1, 1:ntp) = vb(1, 1:ntp) - vbcb(1)
            vh(2, 1:ntp) = vb(2, 1:ntp) - vbcb(2)
            vh(3, 1:ntp) = vb(3, 1:ntp) - vbcb(3)
         end where
      end associate

      return
   end subroutine helio_coord_b2h_tp

   module subroutine helio_coord_vb2vh_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from barycentric to heliocentric coordinates (velocity only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vb2vh.f90 
      !! Adapted from Hal Levison's Swift routine coord_vb2vh.f 
      use swiftest
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self !! Helio massive body object
      class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)              :: i

      associate(npl => self%nbody, vbcb => cb%vb, xh => self%xh, vb => self%vb, &
                vh => self%vh, Mcb => cb%Gmass, Mpl => self%Gmass)
         do i = 1, NDIM
            vbcb(i) = -sum(Mpl(1:npl) * vb(i, 1:npl)) / Mcb
            vh(i, 1:npl) = vb(i, 1:npl) - vbcb(i)
         end do
      end associate

      return
   end subroutine helio_coord_vb2vh_pl

   module subroutine helio_coord_vb2vh_tp(self, vbcb)
      !! author: David A. Minton
      !!
      !! Convert test particles from barycentric to heliocentric coordinates (velocity only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vb2vh_tp.f90 
      !! Adapted from Hal Levison's Swift routine coord_vb2h_tp.f 
      use swiftest
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self !! Helio massive body object
      real(DP), dimension(:),       intent(in)    :: vbcb  !! Barycentric velocity of the central body

      associate(ntp => self%nbody, vb => self%vb, vh => self%vh, status => self%status)
         where (status(1:ntp) == ACTIVE)
            vh(1, 1:ntp) = vb(1, 1:ntp) - vbcb(1)
            vh(2, 1:ntp) = vb(2, 1:ntp) - vbcb(2)
            vh(3, 1:ntp) = vb(3, 1:ntp) - vbcb(3)
         end where
      end associate

      return
   end subroutine helio_coord_vb2vh_tp

   module subroutine helio_coord_vh2vb_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from heliocentric to barycentric coordinates (velocity only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vh2vb.f90 
      !! Adapted from Hal Levison's Swift routine coord_vh2b.f 
      use swiftest
      implicit none
      ! Arguments
      class(helio_pl),        intent(inout) :: self !! Helio massive body object
      class(swiftest_cb),     intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)  :: i
      real(DP)      :: msys

      associate(npl => self%nbody, vbcb => cb%vb, vb => self%vb, vh => self%vh, &
                Mcb => cb%Gmass, Mpl => self%Gmass)
         msys = Mcb + sum(Mpl(1:npl))
         do i = 1, NDIM
            vbcb(i) = -sum(Mpl(1:npl) * vh(i, 1:npl)) / msys
            vb(i, 1:npl) = vh(i, 1:npl) + vbcb(i)
         end do
      end associate

      return
   end subroutine helio_coord_vh2vb_pl

   module subroutine helio_coord_vh2vb_tp(self, vbcb)
      !! author: David A. Minton
      !!
      !! Convert test particles from heliocentric to barycentric coordinates (velocity only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vh2vb_tp.f90
      !! Adapted from Hal Levison's Swift routine coord_vh2b_tp.f 
      use swiftest
      implicit none
      ! Arguments
      class(helio_tp),        intent(inout) :: self !! Helio massive body object
      real(DP), dimension(:), intent(in)    :: vbcb !! Barycentric velocity of the central body

      associate(ntp => self%nbody, vb => self%vb, vh => self%vh, status => self%status)
         where (status(1:ntp) == ACTIVE)
            vb(1, 1:ntp) = vh(1, 1:ntp) + vbcb(1)
            vb(2, 1:ntp) = vh(2, 1:ntp) + vbcb(2)
            vb(3, 1:ntp) = vh(3, 1:ntp) + vbcb(3)
         end where
      end associate

      return
   end subroutine helio_coord_vh2vb_tp
end submodule s_helio_coord

