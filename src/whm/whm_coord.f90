submodule (whm_classes) whm_coord_implementations
contains
   module procedure whm_coord_h2j_pl
      !! author: David A. Minton
      !!
      !! Convert from heliocentric to Jacobi coordinates, massive bodies only
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2j.f90 
      !! Adapted from Hal Levison's Swift routine coord_h2j.f 
      use swiftest
      implicit none

      integer(I4B)              :: i
      real(DP)                  :: eta
      real(DP), dimension(NDIM) :: sumx, sumv, cap, capv

      associate(npl => self%nbody)
         eta = cb%Gmass 
         cb%eta = eta
         do i = 1, npl
            eta = eta + self%Gmass(i)
            self%eta(i) = eta
         end do
         cb%xj(:) = 0.0_DP
         cb%vj(:) = 0.0_DP
         if (npl > 0) then
            self%xj(1, :) = self%xh(1, :)
            self%vj(1, :) = self%vh(1, :)
            sumx(:) = 0.0_DP
            sumv(:) = 0.0_DP
         end if
         do i = 2, npl
            sumx(:) = sumx(:) + self%Gmass(i - 1) * self%xh(i - 1, :)
            sumv(:) = sumv(:) + self%Gmass(i - 1) * self%vh(i - 1, :)
            cap(:) = sumx(:) / self%eta(i - 1)
            capv(:) = sumv(:) / self%eta(i - 1)
            self%xj(i, :) = self%xh(i, :) - cap(:)
            self%vj(i, :) = self%vh(i, :) - capv(:)
         end do
      end associate
    
      return
   end procedure whm_coord_h2j_pl

   module procedure whm_coord_j2h_pl
      !! author: David A. Minton
      !!
      !! Convert from Jacobi to heliocentric coordinates, massive bodies only
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_j2h.f90 
      !! Adapted from Hal Levison's Swift routine coord_j2h.f 
      use swiftest
      implicit none
      integer(I4B)              :: i
      real(DP)                  :: eta
      real(DP), dimension(NDIM) :: sumx, sumv
     
      associate(npl => self%nbody)
         eta = cb%Gmass
         cb%eta = eta
         do i = 1, npl
            eta = eta + self%Gmass(i)
            self%eta(i) = eta
         end do
         if (npl > 0) then
            self%xh(1, :) = self%xj(1, :)
            self%vh(1, :) = self%vj(1,:)
            sumx(:) = 0.0_DP
            sumv(:) = 0.0_DP
         end if
         do i = 2, npl 
            sumx(:) = sumx(:) + self%Gmass(i) * self%xj(i - 1, :) / self%eta(i - 1)
            sumv(:) = sumv(:) + self%Gmass(i) * self%vj(i - 1, :) / self%eta(i -1)
            self%xh(i, :) = self%xj(i, :) + sumx(:)
            self%vh(i, :) = self%vj(i, :) + sumv(:)
         end do
      end associate
    
      return
   end procedure whm_coord_j2h_pl

   module procedure whm_coord_vh2vj_pl
      !! author: David A. Minton
      !!
      !! Convert from heliocentric to Jadcobi coordinates, massive body velocities only
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vh2vj.f90 
      !! Adapted from Hal Levison's Swift routine coord_vh2vj.f 
      use swiftest
      implicit none
      integer(I4B)              :: i
      real(DP)                  :: eta
      real(DP), dimension(NDIM) :: sumv, capv

      associate(npl => self%nbody)
         eta = cb%Gmass 
         cb%eta = eta
         do i = 1, npl
            eta = eta + self%Gmass(i)
            self%eta(i) = eta
         end do
         cb%vj(:) = 0.0_DP
         if (npl > 0) then
            self%vj(1, :) = self%vh(1, :)
            sumv(:) = 0.0_DP
         end if
         do i = 2, npl
            sumv(:) = sumv(:) + self%Gmass(i - 1) * self%vh(i - 1, :)
            capv(:) = sumv(:) / self%eta(i - 1)
            self%vj(i, :) = self%vh(i, 1) - capv(:)
         end do
      end associate
    
      return
   end procedure whm_coord_vh2vj_pl
end submodule whm_coord_implementations

