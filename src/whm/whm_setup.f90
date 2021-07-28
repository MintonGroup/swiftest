submodule(whm_classes) s_whm_setup
   use swiftest
contains
   module subroutine whm_setup_pl(self,n)
      !! author: David A. Minton
      !!
      !! Allocate WHM planet structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self !! Swiftest test particle object
      integer(I4B),                     intent(in) :: n    !! Number of test particles to allocate
      !> Call allocation method for parent class
      call setup_pl(self, n) 
      if (n <= 0) return

      allocate(self%eta(n))
      allocate(self%muj(n))
      allocate(self%xj(NDIM, n))
      allocate(self%vj(NDIM, n))
      allocate(self%ir3j(n))

      self%eta(:)   = 0.0_DP
      self%muj(:)   = 0.0_DP
      self%xj(:,:)  = 0.0_DP
      self%vj(:,:)  = 0.0_DP
      self%ir3j(:) = 0.0_DP

      return
   end subroutine whm_setup_pl 

   module subroutine whm_setup_tp(self,n)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(whm_tp),                 intent(inout) :: self   !! WHM test particle data structure
      integer,                       intent(in)    :: n      !! Number of test particles to allocate
      !> Call allocation method for parent class
      call setup_tp(self, n) 
      if (n <= 0) return

      return
   end subroutine whm_setup_tp

   module subroutine whm_util_set_mu_eta_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Sets the Jacobi mass value eta for all massive bodies
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self    !! Swiftest system object
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structure
      ! Internals
      integer(I4B)                                 :: i

      associate(pl => self, npl => self%nbody)
         if (npl == 0) return
         call util_set_mu_pl(pl, cb)
         pl%eta(1) = cb%Gmass + pl%Gmass(1)
         pl%muj(1) = pl%eta(1)
         do i = 2, npl
            pl%eta(i) = pl%eta(i - 1) + pl%Gmass(i)
            pl%muj(i) = cb%Gmass * pl%eta(i) / pl%eta(i - 1)
         end do
      end associate

   end subroutine whm_util_set_mu_eta_pl

   module subroutine whm_setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a WHM nbody system from files
      !!
      implicit none
      ! Arguments
      class(whm_nbody_system),    intent(inout) :: self    !! Swiftest system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 

      call setup_initialize_system(self, param)
      ! Make sure that the discard list gets allocated initially
      call self%tp_discards%setup(self%tp%nbody)
      call self%pl%set_mu(self%cb)
      call self%tp%set_mu(self%cb)
      if (param%lgr) then
         call self%pl%v2pv(param)
         call self%tp%v2pv(param)
      end if

   end subroutine whm_setup_initialize_system

end submodule s_whm_setup