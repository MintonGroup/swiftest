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

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         do i = 1, NDIM
            cb%vb(i) = -sum(pl%Gmass(1:npl) * pl%vb(i, 1:npl)) / cb%Gmass
            pl%vh(i, 1:npl) = pl%vb(i, 1:npl) - cb%vb(i)
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

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         where (tp%lmask(1:ntp))
            tp%vh(1, 1:ntp) = tp%vb(1, 1:ntp) - vbcb(1)
            tp%vh(2, 1:ntp) = tp%vb(2, 1:ntp) - vbcb(2)
            tp%vh(3, 1:ntp) = tp%vb(3, 1:ntp) - vbcb(3)
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

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         msys = cb%Gmass + sum(pl%Gmass(1:npl))
         do i = 1, NDIM
            cb%vb(i) = -sum(pl%Gmass(1:npl) * pl%vh(i, 1:npl)) / msys
            pl%vb(i, 1:npl) = pl%vh(i, 1:npl) + cb%vb(i)
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

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         where (tp%lmask(1:ntp))
            tp%vb(1, 1:ntp) = tp%vh(1, 1:ntp) + vbcb(1)
            tp%vb(2, 1:ntp) = tp%vh(2, 1:ntp) + vbcb(2)
            tp%vb(3, 1:ntp) = tp%vh(3, 1:ntp) + vbcb(3)
         end where
      end associate

      return
   end subroutine helio_coord_vh2vb_tp
   
end submodule s_helio_coord

