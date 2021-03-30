submodule (whm_classes) whm_coord_implementations
contains
   module subroutine whm_coord_h2j_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert from heliocentric to Jacobi coordinates, massive bodies only
      !!
      !! Uses pre-computed eta rather than computing it each time
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2j.f90 
      !!
      !! Adapted from Hal Levison's Swift routine coord_h2j.f 
      use swiftest
      implicit none
      ! Arguments
      class(whm_pl),         intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),    intent(inout) :: cb     !! Swiftest central body particle data structuree
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: sumx, sumv, cap, capv

      associate(npl => self%nbody, GMpl => self%Gmass, eta => self%eta, xh => self%xh, vh => self%vh, &
                xj => self%xj, vj => self%vj)
         if (npl == 0) return
         xj(:, 1) = xh(:, 1)
         vj(:, 1) = vh(:, 1)
         sumx(:) = 0.0_DP
         sumv(:) = 0.0_DP
         do i = 2, npl
            sumx(:) = sumx(:) + GMpl(i - 1) * xh(:, i - 1)
            sumv(:) = sumv(:) + GMpl(i - 1) * vh(:, i - 1)
            cap(:) = sumx(:) / eta(i - 1)
            capv(:) = sumv(:) / eta(i - 1)
            xj(:, i) = xh(:, i) - cap(:)
            vj(:, i) = vh(:, i) - capv(:)
         end do
      end associate
    
      return
   end subroutine whm_coord_h2j_pl

   module subroutine whm_coord_j2h_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert from Jacobi to heliocentric coordinates, massive bodies only.
      !! 
      !! Uses pre-computed eta rather than computing it each time
      !!
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_j2h.f90 
      !!
      !! Adapted from Hal Levison's Swift routine coord_j2h.f 
      use swiftest
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),                 intent(inout) :: cb     !! Swiftest central body particle data structuree
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: sumx, sumv

      associate(npl => self%nbody, GMpl => self%Gmass, eta => self%eta, xh => self%xh, vh => self%vh, &
                xj => self%xj, vj => self%vj)
         if (npl == 0) return
         xh(:, 1) = xj(:, 1)
         vh(:, 1) = vj(:, 1)
         sumx(:) = 0.0_DP
         sumv(:) = 0.0_DP
         do i = 2, npl 
            sumx(:) = sumx(:) + GMpl(i - 1) * xj(:, i - 1) / eta(i - 1)
            sumv(:) = sumv(:) + GMpl(i - 1) * vj(:, i - 1) / eta(i - 1)
            xh(:, i) = xj(:, i) + sumx(:)
            vh(:, i) = vj(:, i) + sumv(:)
         end do
      end associate
    
      return
   end subroutine whm_coord_j2h_pl

   module subroutine whm_coord_vh2vj_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert from heliocentric to Jadcobi coordinates, massive body velocities only
      !! 
      !! Uses pre-computed eta rather than computing it each time
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vh2vj.f90 
      !!
      !! Adapted from Hal Levison's Swift routine coord_vh2vj.f 
      use swiftest
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),                 intent(inout) :: cb     !! Swiftest central body particle data structuree
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: sumv, capv

      associate(npl => self%nbody, GMpl => self%Gmass, vh => self%vh, vj => self%vj, eta => self%eta)
         if (npl == 0) return
         vj(:, 1) = vh(:, 1)
         sumv(:) = 0.0_DP
         do i = 2, npl
            sumv(:) = sumv(:) + GMpl(i - 1) * vh(:, i - 1)
            capv(:) = sumv(:) / eta(i - 1)
            vj(:, i) = vh(:, i) - capv(:)
         end do
      end associate
    
      return
   end subroutine whm_coord_vh2vj_pl
end submodule whm_coord_implementations

