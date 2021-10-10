submodule (swiftest_classes) s_util_index
   use swiftest
contains

   module pure subroutine util_flatten_eucl_ij_to_k(n, i, j, k)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix for pl-pl interactions.
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      integer(I4B), intent(in)  :: n !! Number of bodies
      integer(I4B), intent(in)  :: i !! Index of the ith body
      integer(I4B), intent(in)  :: j !! Index of the jth body
      integer(I8B), intent(out) :: k !! Index of the flattened matrix
      ! Internals
      integer(I8B) :: i8, j8, n8
     
      i8 = int(i, kind=I8B)
      j8 = int(j, kind=I8B)
      n8 = int(n, kind=I8B)
      k = (i8 - 1_I8B) * n8 - i8 * (i8 - 1_I8B) / 2_I8B + (j8 - i8)

      return
   end subroutine util_flatten_eucl_ij_to_k


   module pure subroutine util_flatten_eucl_k_to_ij(n, k, i, j)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns k index into i,j indices for use in the Euclidean distance matrix for pl-pl interactions.
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      integer(I4B), intent(in)  :: n !! Number of bodies
      integer(I8B), intent(in)  :: k !! Index of the flattened matrix
      integer(I4B), intent(out) :: i !! Index of the ith body
      integer(I4B), intent(out) :: j !! Index of the jth body
      ! Internals
      integer(I8B) :: kp, p, i8, j8, n8

      n8 = int(n, kind=I8B)
    
      kp = n8 * (n8 - 1_I8B) / 2_I8B - k
      p = floor((sqrt(1._DP + 8_I8B * kp) - 1_I8B) / 2_I8B)
      i8 = n8 - 1_I8B - p
      j8 = k - (n8 - 1_I8B) * (n8 - 2_I8B) / 2_I8B + p * (p + 1_I8B) / 2_I8B + 1_I8B

      i = int(i8, kind=I4B)
      j = int(j8, kind=I4B)

      return
   end subroutine util_flatten_eucl_k_to_ij


   module subroutine util_flatten_eucl_plpl(self, param)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix for pl-pl interactions for a Swiftest massive body object
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i, j, npl, err
      integer(I8B) :: k

      npl = int(self%nbody, kind=I8B)
      associate(nplpl => self%nplpl)
         nplpl = (npl * (npl - 1) / 2) ! number of entries in a strict lower triangle, npl x npl
         if (param%lflatten_interactions) then
            if (allocated(self%k_plpl)) deallocate(self%k_plpl) ! Reset the index array if it's been set previously
            allocate(self%k_plpl(2, nplpl), stat=err)
            if (err /=0) then ! An error occurred trying to allocate this big array. This probably means it's too big to fit in memory, and so we will force the run back into triangular mode
               param%lflatten_interactions = .false.
            else
               do concurrent (i=1:npl, j=1:npl, j>i)
                  call util_flatten_eucl_ij_to_k(npl, i, j, k)
                  self%k_plpl(1, k) = i
                  self%k_plpl(2, k) = j
               end do
            end if
         end if
      end associate

      return
   end subroutine util_flatten_eucl_plpl


   module subroutine util_flatten_eucl_pltp(self, pl, param)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix for pl-tp interactions
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
      class(swiftest_pl),         intent(in)    :: pl    !! Swiftest massive body object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I8B) :: i, j, counter, npl, ntp

      ntp = int(self%nbody, kind=I8B)
      npl = int(pl%nbody, kind=I8B)
      associate(npltp => self%npltp)
         npltp = npl * ntp
         if (allocated(self%k_pltp)) deallocate(self%k_pltp) ! Reset the index array if it's been set previously
         allocate(self%k_pltp(2, npltp))
         do i = 1_I8B, npl
            counter = (i - 1_I8B) * npl + 1_I8B
            do j = 1_I8B,  ntp
               self%k_pltp(1, counter) = i
               self%k_pltp(2, counter) = j
               counter = counter + 1_I8B
            end do
         end do
      end associate

      return
   end subroutine util_flatten_eucl_pltp

end submodule s_util_index
