submodule(whm_classes) s_whm_setup
contains
   module procedure whm_setup_pl
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for parent class
      call self%swiftest_pl%alloc(n)
      if (n <= 0) return

      real(DP), dimension(:),   allocatable :: eta     ! Jacobi mass
      real(DP), dimension(:,:), allocatable :: xj      ! Jacobi position
      real(DP), dimension(:,:), allocatable :: vj      ! Jacobi velocity
      real(DP), dimension(:,:), allocatable :: ah1     ! First term of heliocentric acceleration
      real(DP), dimension(:,:), allocatable :: ah2     ! Second term of heliocentric acceleration
      real(DP), dimension(:,:), allocatable :: ah3     ! Third term of heliocentric acceleration
      real(DP), dimension(:,:), allocatable :: ah      ! Total heliocentric acceleration
      
      allocate(self%eta(n))
      allocate(self%xj(n, NDIM))
      allocate(self%vj(n, NDIM))
      allocate(self%ah1(n, NDIM))
      allocate(self%ah2(n, NDIM))
      allocate(self%ah3(n, NDIM))
      allocate(self%ah(n, NDIM))

      self%eta(:)   = 0.0_DP
      self%xj(:,:)  = 0.0_DP
      self%vj(:,:)  = 0.0_DP
      self%ah1(:,:) = 0.0_DP
      self%ah2(:,:) = 0.0_DP
      self%ah3(:,:) = 0.0_DP
      self%ah(:,:)  = 0.0_DP

      return
   end procedure whm_setup_pl 

   module procedure whm_setup_tp
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for parent class
      call self%swiftest_tp%alloc(n)
      if (n <= 0) return
      allocate(self%ah(n, NDIM))

      selt%ah(:,:) = 0.0_DP

      return
   end procedure whm_setup_tp
end submodule s_whm_setup
