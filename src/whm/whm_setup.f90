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

   module subroutine whm_setup_set_mu_eta_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Sets the Jacobi mass value eta for all massive bodies
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self    !! Swiftest system object
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structure
      ! Internals
      integer(I4B)                                 :: i

      associate(pl => self, npl => self%nbody,  GMpl => self%Gmass, muj => self%muj, &
                eta => self%eta, GMcb => cb%Gmass)
         if (npl == 0) return
         call setup_set_mu_pl(pl, cb)
         eta(1) = GMcb + GMpl(1)
         muj(1) = eta(1)
         do i = 2, npl
            eta(i) = eta(i - 1) + GMpl(i)
            muj(i) = GMcb * eta(i) / eta(i - 1)
         end do
      end associate

   end subroutine whm_setup_set_mu_eta_pl

   module subroutine whm_setup_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a WHM nbody system from files
      !!
      implicit none
      ! Arguments
      class(whm_nbody_system),       intent(inout) :: self    !! Swiftest system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters of on parameters 
      call io_read_initialize_system(self, param)
      ! Make sure that the discard list gets allocated initially
      call self%tp_discards%setup(self%tp%nbody)

      if (self%pl%nbody > 0) then
         select type(pl => self%pl)
         class is (whm_pl)
            call pl%set_mu(self%cb)
            if (param%lgr) call pl%gr_vh2pv(param)
            !call pl%eucl_index()
         end select
      end if

      if (self%tp%nbody > 0) then
         select type(tp => self%tp)
         class is (whm_tp)
            call tp%set_mu(self%cb)
            if (param%lgr) call tp%gr_vh2pv(param)
         end select
      end if

   end subroutine whm_setup_system

   module subroutine whm_setup_set_ir3j(self)
      !! author: David A. Minton
      !!
      !! Sets the inverse Jacobi and heliocentric radii cubed (1/rj**3 and 1/rh**3)
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self    !! WHM massive body object
      ! Internals
      integer(I4B)                                 :: i
      real(DP)                                     :: r2, ir

      if (self%nbody > 0) then
         do i = 1, self%nbody
            r2 = dot_product(self%xh(:, i), self%xh(:, i))
            ir = 1.0_DP / sqrt(r2)
            self%ir3h(i) = ir / r2
            r2 = dot_product(self%xj(:, i), self%xj(:, i))
            ir = 1.0_DP / sqrt(r2)
            self%ir3j(i) = ir / r2
         end do
      end if
   end subroutine whm_setup_set_ir3j

   module subroutine whm_setup_set_beg_end(self, xbeg, xend, vbeg)
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of xbeg and xend
      implicit none
      ! Arguments
      class(whm_nbody_system),  intent(inout)           :: self !! WHM nbody system object
      real(DP), dimension(:,:), intent(in),    optional :: xbeg, xend
      real(DP), dimension(:,:), intent(in),    optional :: vbeg ! vbeg is an unused variable to keep this method forward compatible with RMVS

      if (present(xbeg)) then
         if (allocated(self%xbeg)) deallocate(self%xbeg)
         allocate(self%xbeg, source=xbeg)
      end if
      if (present(xend)) then
         if (allocated(self%xend)) deallocate(self%xend)
         allocate(self%xend, source=xend)
      end if

      return

   end subroutine whm_setup_set_beg_end

end submodule s_whm_setup