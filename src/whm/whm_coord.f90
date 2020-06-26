submodule (whm_classes) whm_coord
contains
    module procedure coord_h2j_pl
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
            eta = cb%mass 
            cb%eta = eta
            do i = 1, npl
                eta = eta + self%mass(i)
                self%eta(i) = eta
            end do
            cb%xj(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
            cb%vj(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
            if (npl > 1) then
                self%xj(:,1) = self%xh(:,1)
                self%vj(:,1) = self%vh(:,1)
                sumx(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
                sumv(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
            end if
            do i = 1, npl - 1
                sumx(:) = sumx(:) + self%mass(i) * self%xh(:, i)
                sumv(:) = sumv(:) + self%mass(i) * self%vh(:, i)
                cap(:) = sumx(:) / self%eta(i)
                capv(:) = sumv(:) / self%eta(i)
                self%xj(:, i + 1)    = self%xh(:, i + 1) - cap(:)
                self%vj(:, i + 1) = self%vh(:, i + 1) - capv(:)
            end do
        end associate
    
        return
    end procedure coord_h2j_pl

    module procedure coord_j2h_pl
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
            eta = cb%mass
            cb%eta = eta
            do i = 1, npl
                eta = eta + self%mass(i)
                self%eta(i) = eta
            end DO
            if (npl > 1) then
                self%xh(:,1) = self%xj(:,1)
                self%vh(:,1) = self%vj(:,1)
                sumx(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
                sumv(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
            end if
            do i = 1, npl - 1
                sumx(:) = sumx(:) + self%mass(i) * self%xj(:, i) / self%eta(i)
                sumv(:) = sumv(:) + self%mass(i) * self%vj(:, i) / self%eta(i)
                self%xh(:, i + 1) = self%xj(:, i + 1) + sumx(:)
                self%vh(:, i + 1) = self%vj(:, i + 1) + sumv(:)
            end do
        end associate
    
        return
    end procedure coord_j2h_pl

    module procedure coord_vh2vj_pl
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
            eta = cb%mass 
            cb%eta = eta
            do i = 1, npl
                eta = eta + self%mass(i)
                self%eta(i) = eta
            end do
            cb%vj(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
            if (npl > 1) then
                self%vj(:,1) = self%vh(:,1)
                sumv(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
            end if
            do i = 1, npl - 1
                sumv(:) = sumv(:) + self%mass(i) * self%vh(:, i)
                capv(:) = sumv(:) / self%eta(i)
                self%vj(:, i + 1) = self%vh(:, i + 1) - capv(:)
            end do
        end associate
    
        return
    end procedure coord_vh2vj_pl
end submodule whm_coord

