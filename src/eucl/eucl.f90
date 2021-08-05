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
      associate(nplpl => self%nplpl)
         nplpl = (npl * (npl - 1) / 2) ! number of entries in a strict lower triangle, nplm x npl, minus first column
         if (allocated(self%k_plpl)) deallocate(self%k_plpl) ! Reset the index array if it's been set previously
         allocate(self%k_plpl(2, nplpl))
         do i = 1, npl
            counter = (i - 1_I8B) * npl - i * (i - 1_I8B) / 2_I8B + 1_I8B
            do j = i + 1_I8B, npl
               self%k_plpl(1, counter) = i
               self%k_plpl(2, counter) = j
               counter = counter + 1_I8B
            end do
         end do
      end associate

      return
   end subroutine eucl_dist_index_plpl

 end submodule s_eucl
