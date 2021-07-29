submodule(swiftest_classes) s_util_coord
   use swiftest
contains

   module subroutine util_coord_h2b_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from heliocentric to barycentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b.f90 
      !! Adapted from Hal Levison's Swift routine coord_h2b.f 
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)  :: i
      real(DP)      :: msys
      real(DP), dimension(NDIM) :: xtmp, vtmp

      associate(pl => self, npl => self%nbody)
         msys = cb%Gmass
         xtmp(:) = 0.0_DP
         vtmp(:) = 0.0_DP
         do i = 1, npl
            msys = msys + pl%Gmass(i)
            xtmp(:) = xtmp(:) + pl%Gmass(i) * pl%xh(:,i)
            vtmp(:) = vtmp(:) + pl%Gmass(i) * pl%vh(:,i)
         end do
         cb%xb(:) = -xtmp(:) / msys
         cb%vb(:) = -vtmp(:) / msys
         do i = 1, npl
            pl%xb(:,i) = pl%xh(:,i) + cb%xb(:)
            pl%vb(:,i) = pl%vh(:,i) + cb%vb(:)
         end do
      end associate

      return
   end subroutine util_coord_h2b_pl


   module subroutine util_coord_h2b_tp(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from heliocentric to barycentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b_tp.f90 
      !! Adapted from Hal Levison's Swift routine coord_h2b_tp.f 
      implicit none
      ! Arguments
      class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
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
   end subroutine util_coord_h2b_tp


   module subroutine util_coord_b2h_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from barycentric to heliocentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_b2h.f90 
      !! Adapted from Hal Levison's Swift routine coord_b2h.f 
      implicit none
      ! Arguments
      class(swiftest_pl),     intent(inout) :: self !! Swiftest massive body object
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
   end subroutine util_coord_b2h_pl


   module subroutine util_coord_b2h_tp(self, cb)
      !! author: David A. Minton
      !!
      !! Convert test particles from barycentric to heliocentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_b2h_tp.f90 
      !! Adapted from Hal Levison's Swift routine coord_b2h_tp.f 
      implicit none
      ! Arguments
      class(swiftest_tp),     intent(inout) :: self !! Swiftest massive body object
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
   end subroutine util_coord_b2h_tp
   
end submodule s_util_coord