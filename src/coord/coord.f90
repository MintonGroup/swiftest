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

      associate(n => self%nbody)
         select type(self)
         class is (swiftest_pl)
            msys = cb%mass + sum(self%mass(1:n), self%status(i) == ACTIVE)
            do i = 1, NDIM
               cb%xb(i) = -sum(self%mass(1:n) * self%xh(i, 1:n), self%status(i) == ACTIVE) / msys
               cb%vb(i) = -sum(self%mass(1:n) * self%vh(i, 1:n), self%status(i) == ACTIVE) / msys
            end do
         end select
         do concurrent(i = 1:n, self%status(i) == ACTIVE)
            self%xb(:, i) = self%xh(:, i) + cb%xb(:)
            self%vb(:, i) = self%vh(:, i) + cb%vb(:)
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

      associate(n => self%nbody)
         do concurrent(i = 1:n, self%status(i) == ACTIVE)
            self%xh(:, i) = self%xb(:, i) - cb%xb(:)
            self%vh(:, i) = self%vb(:, i) - cb%vb(:)
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
      associate(n => self%nbody)
         select type(self)
         class is (swiftest_pl)
            do i = 1, NDIM
               cb%vb(i) = -sum(self%mass(1:n) * self%vb(i, 1:n), self%status(1:n) == ACTIVE) / cb%mass
            end do
         end select
         do concurrent(i = 1:n, self%status(i) == ACTIVE)
            self%vh(:, i) = self%vb(:, i) - cb%vb(:)
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

      associate(n => self%nbody)
         select type(self)
         class is (swiftest_pl)
            msys = cb%mass + sum(self%mass(1:n))
            do i = 1, NDIM
               cb%vb(i) = -sum(self%mass(1:n) * self%vh(i, 1:n), self%status(1:n) == ACTIVE) / msys
            end do
         end select
         do concurrent(i = 1:n, self%status(i) == ACTIVE)
            self%vb(:, i) = self%vh(:, i) + cb%vb(:)
         end do
      end associate

      return
   end procedure coord_vh2vb_body
end submodule coord
