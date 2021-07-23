submodule (swiftest_classes) s_eucl
   use swiftest
contains
   module subroutine eucl_dist_index_plpl(self)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      class(swiftest_pl),             intent(inout) :: self  !! Swiftest massive body objec
      ! Internals
      integer(I8B) :: i, j, counter, npl

      npl = int(self%nbody, kind=I8B)
      associate(num_comparisons => self%num_comparisons)
         num_comparisons = (npl * (npl - 1) / 2) ! number of entries in a strict lower triangle, nplm x npl, minus first column
         if (allocated(self%k_eucl)) deallocate(self%k_eucl) ! Reset the index array if it's been set previously
         if (allocated(self%irij3)) deallocate(self%irij3)  
         allocate(self%k_eucl(2, num_comparisons))
         allocate(self%irij3(num_comparisons))
         do i = 1, npl
            counter = (i - 1_I8B) * npl - i * (i - 1_I8B) / 2_I8B + 1_I8B
            do j = i + 1_I8B, npl
               self%k_eucl(1, counter) = i
               self%k_eucl(2, counter) = j
               counter = counter + 1_I8B
            end do
         end do
      end associate
      return

   end subroutine eucl_dist_index_plpl

   module subroutine eucl_dist_index_pltp(self, pl)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751      implicit none
      class(swiftest_tp),             intent(inout) :: self  !! Swiftest test particle object
      class(swiftest_pl),             intent(inout) :: pl    !! Swiftest massive body object
   end subroutine eucl_dist_index_pltp

   module subroutine eucl_irij3_plpl(self)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Efficient parallel loop-blocking algrorithm for evaluating the Euclidean distance matrix for planet-planet
      implicit none
      ! Arguments
      class(swiftest_pl),             intent(inout) :: self  !! Swiftest massive body object
      ! Internals
      integer(I4B) :: k, i, j
      real(DP), dimension(NDIM) :: dx
      real(DP) :: rji2

      associate(k_eucl => self%k_eucl, xh => self%xh, irij3 => self%irij3, nk => self%num_comparisons)
         do k = 1, nk
            i = k_eucl(1, k)
            j = k_eucl(2, k)
            dx(:) = xh(:, j) - xh(:, i)
            rji2  = dot_product(dx(:), dx(:))
            irij3(k) = 1.0_DP / (rji2 * sqrt(rji2))
         end do
      end associate

      return
   end subroutine eucl_irij3_plpl

 end submodule s_eucl
