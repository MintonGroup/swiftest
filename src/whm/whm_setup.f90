submodule(whm_classes) s_whm_setup
   use swiftest
contains

   module subroutine whm_setup_pl(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate WHM planet structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(whm_pl),             intent(inout) :: self  !! Swiftest test particle object
      integer(I4B),              intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class
      call setup_pl(self, n, param) 
      if (n < 0) return

      if (allocated(self%eta)) deallocate(self%eta)
      if (allocated(self%muj)) deallocate(self%muj)
      if (allocated(self%xj)) deallocate(self%xj)
      if (allocated(self%vj)) deallocate(self%vj)
      if (allocated(self%ir3j)) deallocate(self%ir3j)

      if (n == 0) return

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


   module subroutine whm_util_set_mu_eta_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Sets the Jacobi mass value eta for all massive bodies
      implicit none
      ! Arguments
      class(whm_pl),      intent(inout) :: self   !! WHM system object
      class(swiftest_cb), intent(inout) :: cb     !! Swiftest central body object
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

      return
   end subroutine whm_util_set_mu_eta_pl


   module subroutine whm_setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a WHM nbody system from files
      !!
      implicit none
      ! Arguments
      class(whm_nbody_system),    intent(inout) :: self   !! WHM nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 

      call setup_initialize_system(self, param)
      ! First we need to make sure that the massive bodies are sorted by heliocentric distance before computing jacobies
      call util_set_ir3h(self%pl)

      ! Make sure that the discard list gets allocated initially
      call self%tp_discards%setup(0, param)
      call self%pl%set_mu(self%cb)
      call self%tp%set_mu(self%cb)
      if (param%lgr) then
         call self%pl%v2pv(param)
         call self%tp%v2pv(param)
      end if

      return
   end subroutine whm_setup_initialize_system

end submodule s_whm_setup