!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_swiftest_obl
contains

   pure function matinv3(A) result(B)
   !! Performs a direct calculation of the inverse of a 3×3 matrix.
   !!
   !! from https://fortranwiki.org/fortran/show/Matrix+inversion
   !!

      complex(wp), intent(in) :: A(3,3)   !! Matrix
      complex(wp)             :: B(3,3)   !! Inverse matrix
      complex(wp)             :: detinv

      ! Calculate the inverse determinant of the matrix
      detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
               - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
               + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

      ! Calculate the inverse of the matrix
      B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
      B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
      B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
      B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
      B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
      B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
      B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
      B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
      B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
   end function

   module subroutine swiftest_obl_rot_matrix(self, nbody_system, rot_matrix, rot_matrix_inv)
      !! author: Kaustub P. Anand
      !! 
      !! Generate a rotation matrix and its inverse to rotate the coordinate frame to align the rotation axis along the z axis for correct spin calculation
      !! 

      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   !! Swiftest body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      real(DP), dimension(NDIM, NDIM), intent(inout) :: rot_matrix ! rotation matrix and its inverse
      real(DP), dimension(NDIM, NDIM), intent(inout) :: rot_matrix_inv ! inverse of the rotation matrix

      ! Internals
      real(DP)     :: theta ! angle to rotate it through
      real(DP), dimension(NDIM) :: u ! unit vector about which we rotate
      
      rot_matrix(:, :) = 0.0_DP
      rot_matrix_inv(:, :) = 0.0_DP

      if (self%nbody == 0) return

      associate(n => self%nbody, cb => nbody_system%cb)
         if (cb%rot(0) == 0 .and. cb%rot(1) == 0) then
            do (i = 1, NDIM)
                  rot_matrix(i, i) = 1.0
                  rot_matrix_inv(i, i) = 1.0
            end do

            return ! rotation axis is about the z-axis, no need to change
         end if
         
         u(:) = .unit. (cb%rot(:) .cross. (0, 0, 1))
         theta = acos(dot_product((.unit. cb%rot(:)), (0, 0, 1)))
         
         S_matrix = ((0, -u(3), u(2)), (u(3), 0, -u(1)), (-u(2), u(1), 0)) ! skew-symmetric matrix
         ! assuming NDIM = 3
         ! CHECK for a general formula for the skew-symmetric matrix

         do i = 1, NDIM
            do j = 1, NDIM
               if (i == j) then
                  rot_matrix(i, j) = rot_matrix(i, j) + cos(theta) ! identity matrix
                  continue
               end if

               rot_matrix(i, j) = rot_matrix(i, j) + u(i) * u(j) * (1 - cos(theta)) + S_matrix(i, j) * sin(theta) ! Skew-symmetric matrix + Tensor product matrix

            end do
         end do

         rot_matrix_inv = matinv3(rot_matrix)

      end associate

      return
      end subroutine swiftest_obl_rot_matrix
   
   module subroutine swiftest_obl_acc_body(self, nbody_system)
      !! author: David A. Minton
      !!
      !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      !!      Returned values do not include monopole term or terms higher than J4
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_acc.f90 and obl_acc_tp.f90
      !! Adapted from Hal Levison's Swift routine obl_acc.f and obl_acc_tp.f
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   !! Swiftest body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      ! Internals
      integer(I4B) :: i
      real(DP)     :: r2, irh, rinv2, t0, t1, t2, t3, fac1, fac2
      real(DP), dimension(NDIM)       :: rh_transformed ! rotated position vector
      real(DP), dimension(NDIM, NDIM) :: rot_matrix, rot_matrix_inv ! rotation matrix and its inverse

      if (self%nbody == 0) return

      ! generate the rotation matrix
      call swiftest_obl_rot_matrix(self, nbody_system, rot_matrix, rot_matrix_inv)

      associate(n => self%nbody, cb => nbody_system%cb)
         self%aobl(:,:) = 0.0_DP
         do concurrent(i = 1:n, self%lmask(i))
            ! rotate the position vectors
            rh_transformed = MATMUL(self%rh(:, i), self%rot_matrix) ! 1x3 vector * 3x3 matrix
            r2 = dot_product(rh_transformed, rh_transformed)
            irh = 1.0_DP / sqrt(r2)
            rinv2 = irh**2
            t0 = -cb%Gmass * rinv2 * rinv2 * irh
            t1 = 1.5_DP * cb%j2rp2
            t2 = rh_transformed(3) * rh_transformed(3) * rinv2
            t3 = 1.875_DP * cb%j4rp4 * rinv2
            fac1 = t0 * (t1 - t3 - (5 * t1 - (14.0_DP - 21.0_DP * t2) * t3) * t2)
            fac2 = 2 * t0 * (t1 - (2.0_DP - (14.0_DP * t2 / 3.0_DP)) * t3)
            self%aobl(:, i) = fac1 * rh_transformed(:)
            self%aobl(3, i) = fac2 * rh_transformed(3) + self%aobl(3, i)

            ! rotate the acceleration and position vectors back to the original coordinate frame
            self%aobl(:, i) = MATMUL(self%aobl(:, i), self%rot_matrix_inv)

         end do
      end associate
      return

   end subroutine swiftest_obl_acc_body


   module subroutine swiftest_obl_acc_pl(self, nbody_system)
      !! author: David A. Minton
      !!
      !! Compute the barycentric accelerations of massive bodies due to the oblateness of the central body
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_acc.f90 and obl_acc_tp.f90
      !! Adapted from Hal Levison's Swift routine obl_acc.f and obl_acc_tp.f
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      ! Internals
      integer(I4B) :: i

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody, cb => nbody_system%cb)
         call swiftest_obl_acc_body(pl, nbody_system)
         cb%aobl(:) = 0.0_DP
         do i = npl, 1, -1
            if (pl%lmask(i)) cb%aobl(:) = cb%aobl(:) - pl%Gmass(i) * pl%aobl(:, i) / cb%Gmass
         end do

         do concurrent(i = 1:npl, pl%lmask(i))
            pl%ah(:, i) = pl%ah(:, i) + pl%aobl(:, i) - cb%aobl(:)
         end do
      end associate

      return

   end subroutine swiftest_obl_acc_pl


   module subroutine swiftest_obl_acc_tp(self, nbody_system)
      !! author: David A. Minton
      !!
      !! Compute the barycentric accelerations of massive bodies due to the oblateness of the central body
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_acc.f90 and obl_acc_tp.f90
      !! Adapted from Hal Levison's Swift routine obl_acc.f and obl_acc_tp.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      ! Internals
      real(DP), dimension(NDIM)                   :: aoblcb
      integer(I4B) :: i

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, cb => nbody_system%cb)
         call swiftest_obl_acc_body(tp, nbody_system)
         if (nbody_system%lbeg) then
            aoblcb = cb%aoblbeg
         else
            aoblcb = cb%aoblend
         end if

         do concurrent(i = 1:ntp, tp%lmask(i))
            tp%ah(:, i) = tp%ah(:, i) + tp%aobl(:, i) - aoblcb(:)
         end do

      end associate
      return

   end subroutine swiftest_obl_acc_tp


   module subroutine swiftest_obl_pot_system(self) 
      !! author: David A. Minton
      !!
      !! Compute the contribution to the total gravitational potential due solely to the oblateness of the central body
      !!    Returned value does not include monopole term or terms higher than J4
      !!
      !!    Reference: MacMillan, W. D. 1958. The Theory of the Potential, (Dover Publications), 363.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_pot.f90 
      !! Adapted from Hal Levison's Swift routine obl_pot.f 
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout)  :: self   !! Swiftest nbody system object
      ! Internals
      integer(I4B) :: i
      real(DP), dimension(self%pl%nbody) :: oblpot_arr

      associate(nbody_system => self, pl => self%pl, npl => self%pl%nbody, cb => self%cb)
         if (.not. any(pl%lmask(1:npl))) return
         do concurrent (i = 1:npl, pl%lmask(i))
            oblpot_arr(i) = swiftest_obl_pot_one(cb%Gmass, pl%Gmass(i), cb%j2rp2, cb%j4rp4, pl%rh(3,i), 1.0_DP / norm2(pl%rh(:,i)))
         end do
         nbody_system%oblpot = sum(oblpot_arr, pl%lmask(1:npl))
      end associate
         
      return
   end subroutine swiftest_obl_pot_system


   elemental function swiftest_obl_pot_one(GMcb, GMpl, j2rp2, j4rp4, zh, irh) result(oblpot)
      !! author: David A. Minton
      !!
      !! Compute the contribution to the total gravitational potential due solely to the oblateness of the central body from a single massive body
      !!    Returned value does not include monopole term or terms higher than J4
      !!
      !!    Reference: MacMillan, W. D. 1958. The Theory of the Potential, (Dover Publications), 363.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_pot.f90 
      !! Adapted from Hal Levison's Swift routine obl_pot.f 
      implicit none
      ! Arguments
      real(DP),     intent(in)  :: GMcb   !! G*mass of the central body
      real(DP),     intent(in)  :: GMpl   !! G*mass of the massive body
      real(DP),     intent(in)  :: j2rp2  !! J_2 / R**2 of the central body
      real(DP),     intent(in)  :: j4rp4  !! J_2 / R**4 of the central body
      real(DP),     intent(in)  :: zh     !! z-component of the heliocentric distance vector of the massive body
      real(DP),     intent(in)  :: irh    !! Inverse of the heliocentric distance magnitude of the massive body
      ! Result
      real(DP)                  :: oblpot !! Gravitational potential
         
      ! Internals
      real(DP)                  :: rinv2, t0, t1, t2, t3, p2, p4
         
      rinv2 = irh**2
      t0 = GMcb * GMpl * rinv2 * irh
      t1 = j2rp2
      t2 = zh**2 * rinv2
      t3 = j4rp4 * rinv2
      p2 = 0.5_DP * (3 * t2 - 1.0_DP)
      p4 = 0.125_DP * ((35 * t2 - 30.0_DP) * t2 + 3.0_DP)
      oblpot = t0 * (t1 * p2 + t3 * p4)
         
      return
   end function swiftest_obl_pot_one

end submodule s_swiftest_obl
