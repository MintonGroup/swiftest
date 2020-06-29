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
      call setup_pl(self, n) 
      if (n <= 0) return

      allocate(self%eta(n))
      allocate(self%xj(n, NDIM))
      allocate(self%vj(n, NDIM))
      allocate(self%ah1(n, NDIM))
      allocate(self%ah2(n, NDIM))
      allocate(self%ah3(n, NDIM))

      self%eta(:)   = 0.0_DP
      self%xj(:,:)  = 0.0_DP
      self%vj(:,:)  = 0.0_DP
      self%ah1(:,:) = 0.0_DP
      self%ah2(:,:) = 0.0_DP
      self%ah3(:,:) = 0.0_DP

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
      call setup_tp(self, n) 
      if (n <= 0) return

      return
   end procedure whm_setup_tp

   module procedure whm_setup_set_eta
      !! author: David A. Minton
      !!
      !! Sets the Jacobi mass value eta for all massive bodies
      implicit none
      integer(I4B) :: i

      associate(npl => self%nbody, Gmsun => cb%Gmass, Gmpl => self%Gmass, mu => self%mu_vec, eta => self%eta)
         if (npl > 0) then
            eta(1) = Gmsun + Gmpl(1)
            mu(1) = eta(1) 
         end if
         do i = 2, npl
            eta(i) = eta(i - 1) + Gmpl(i)
            mu(i) = Gmsun * eta(i) / eta(i - 1)
         end do
      end associate

   end procedure whm_setup_set_eta

   module procedure whm_setup_system
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      implicit none

      call io_read_initialize_system(self, config)

   end procedure whm_setup_system


end submodule s_whm_setup
