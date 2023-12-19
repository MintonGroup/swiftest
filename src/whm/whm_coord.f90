! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (whm) s_whm_coord
   use swiftest
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
      implicit none
      ! Arguments
      class(whm_pl),         intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),    intent(inout) :: cb     !! Swiftest central body particle data structuree
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: sumx, sumv, cap, capv

      if (self%nbody == 0) return

      associate(npl => self%nbody, GMpl => self%Gmass, eta => self%eta, rh => self%rh, vh => self%vh, &
                xj => self%xj, vj => self%vj)
         xj(:, 1) = rh(:, 1)
         vj(:, 1) = vh(:, 1)
         sumx(:) = 0.0_DP
         sumv(:) = 0.0_DP
         do i = 2, npl
            sumx(:) = sumx(:) + GMpl(i - 1) * rh(:, i - 1)
            sumv(:) = sumv(:) + GMpl(i - 1) * vh(:, i - 1)
            cap(:) = sumx(:) / eta(i - 1)
            capv(:) = sumv(:) / eta(i - 1)
            xj(:, i) = rh(:, i) - cap(:)
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
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),                 intent(inout) :: cb     !! Swiftest central body particle data structuree
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: sumx, sumv

      if (self%nbody == 0) return

      associate(npl => self%nbody, GMpl => self%Gmass, eta => self%eta, rh => self%rh, vh => self%vh, &
                xj => self%xj, vj => self%vj)
         rh(:, 1) = xj(:, 1)
         vh(:, 1) = vj(:, 1)
         sumx(:) = 0.0_DP
         sumv(:) = 0.0_DP
         do i = 2, npl 
            sumx(:) = sumx(:) + GMpl(i - 1) * xj(:, i - 1) / eta(i - 1)
            sumv(:) = sumv(:) + GMpl(i - 1) * vj(:, i - 1) / eta(i - 1)
            rh(:, i) = xj(:, i) + sumx(:)
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
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),                 intent(inout) :: cb     !! Swiftest central body particle data structuree
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: sumv, capv

      if (self%nbody == 0) return

      associate(npl => self%nbody, GMpl => self%Gmass, vh => self%vh, vj => self%vj, eta => self%eta)
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

end submodule s_whm_coord

