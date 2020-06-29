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
            msys = cb%mass + sum(self%mass(1:n))
            do i = 1, NDIM
               cb%xb(i) = -sum(self%mass(1:n) * self%xh(1:n, i)) / msys
               cb%vb(i) = -sum(self%mass(1:n) * self%vh(1:n, i)) / msys
            end do
         end select
         where (self%status(1:n) == ACTIVE)
            self%xb(1:n, 1) = self%xh(1:n, 1) + cb%xb(1)
            self%xb(1:n, 2) = self%xh(1:n, 2) + cb%xb(2)
            self%xb(1:n, 3) = self%xh(1:n, 3) + cb%xb(3)
            self%vb(1:n, 1) = self%vh(1:n, 1) + cb%vb(1)
            self%vb(1:n, 2) = self%vh(1:n, 2) + cb%vb(2)
            self%vb(1:n, 3) = self%vh(1:n, 3) + cb%vb(3)
         end where
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
         where (self%status(1:n) == ACTIVE)
            self%xh(1:n, 1) = self%xb(1:n, 1) - cb%xb(1)
            self%xh(1:n, 2) = self%xb(1:n, 2) - cb%xb(2)
            self%xh(1:n, 3) = self%xb(1:n, 3) - cb%xb(3)
            self%vh(1:n, 1) = self%vb(1:n, 1) - cb%vb(1)
            self%vh(1:n, 2) = self%vb(1:n, 2) - cb%vb(2)
            self%vh(1:n, 3) = self%vb(1:n, 3) - cb%vb(3)
         end where
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
               cb%vb(i) = -sum(self%mass(1:n) * self%vb(:, i)) / cb%mass
            end do
         end select
         where (self%status(1:n) == ACTIVE)
            self%vh(1:n, 1) = self%vb(1:n, 1) - cb%vb(1)
            self%vh(1:n, 2) = self%vb(1:n, 2) - cb%vb(2)
            self%vh(1:n, 3) = self%vb(1:n, 3) - cb%vb(3)
         end where
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
               cb%vb(i) = -sum(self%mass(1:n) * self%vh(1:n, i)) / msys
            end do
         end select
         where (self%status(1:n) == ACTIVE)
            self%vb(1:n, 1) = self%vh(1:n, 1) + cb%vb(1)
            self%vb(1:n, 2) = self%vh(1:n, 2) + cb%vb(2)
            self%vb(1:n, 3) = self%vh(1:n, 3) + cb%vb(3)
         end where
      end associate

      return
   end procedure coord_vh2vb_body
end submodule coord
