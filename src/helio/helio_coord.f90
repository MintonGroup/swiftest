submodule (helio_classes) s_helio_coord
   use swiftest
contains
   module subroutine helio_coord_vb2vh_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from barycentric to heliocentric coordinates (velocity only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vb2vh.f90 
      !! Adapted from Hal Levison's Swift routine coord_vb2vh.f 
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

