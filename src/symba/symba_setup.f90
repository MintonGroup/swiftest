submodule(symba_classes) s_symba_setup
   use swiftest
contains

   module subroutine symba_setup_initialize_particle_info_system(self, param) 
      !! author: David A. Minton
      !!
      !! Initializes a new particle information data structure with initial conditions recorded
      implicit none
      ! Argumets
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i

      select type(cbinfo => self%cb%info)
      class is (symba_particle_info)
         cbinfo%origin_type = "Initial conditions"
         cbinfo%origin_time = param%t0
         cbinfo%origin_xh(:) = 0.0_DP
         cbinfo%origin_vh(:) = 0.0_DP
      end select

      select type(plinfo => self%pl%info)
      class is (symba_particle_info)
         do i = 1, self%pl%nbody
            plinfo(i)%origin_type = "Initial conditions"
            plinfo(i)%origin_time = param%t0
            plinfo(i)%origin_xh(:) = self%pl%xh(:,i)
            plinfo(i)%origin_vh(:) = self%pl%vh(:,i)
         end do
      end select

      select type(tpinfo => self%tp%info)
      class is (symba_particle_info)
         do i = 1, self%tp%nbody
            tpinfo(i)%origin_type = "Initial conditions"
            tpinfo(i)%origin_time = param%t0
            tpinfo(i)%origin_xh(:) = self%tp%xh(:,i)
            tpinfo(i)%origin_vh(:) = self%tp%vh(:,i)
         end do
      end select
      call setup_initialize_particle_info_system(self, param)

      return
   end subroutine symba_setup_initialize_particle_info_system


   module subroutine symba_setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize an SyMBA nbody system from files and sets up the planetocentric structures.
      !! This subroutine will also sort the massive bodies in descending order by mass
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self    !! SyMBA system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i, j

      ! Call parent method
      associate(system => self)
         call helio_setup_initialize_system(system, param)
         call system%pltpenc_list%setup(0)
         call system%plplenc_list%setup(0)
         call system%plplcollision_list%setup(0)
      end associate

      return
   end subroutine symba_setup_initialize_system


   module subroutine symba_setup_merger(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate SyMBA test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine symba_setup.f90
      implicit none
      ! Arguments
      class(symba_merger),        intent(inout) :: self  !! SyMBA merger list object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter
      ! Internals
      integer(I4B) :: i

      !> Call allocation method for parent class. In this case, helio_pl does not have its own setup method so we use the base method for swiftest_pl
      call symba_setup_pl(self, n, param) 
      if (n <= 0) return

      if (allocated(self%ncomp)) deallocate(self%ncomp)
      allocate(self%ncomp(n))
      self%ncomp(:) = 0

      return
   end subroutine symba_setup_merger


   module subroutine symba_setup_pl(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate SyMBA test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine symba_setup.f90
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter
      ! Internals
      integer(I4B) :: i

      !> Call allocation method for parent class. In this case, helio_pl does not have its own setup method so we use the base method for swiftest_pl
      call setup_pl(self, n, param) 
      if (n <= 0) return


      if (allocated(self%info)) deallocate(self%info)
      if (allocated(self%lcollision)) deallocate(self%lcollision)
      if (allocated(self%lencounter)) deallocate(self%lencounter)
      if (allocated(self%lmtiny)) deallocate(self%lmtiny)
      if (allocated(self%nplenc)) deallocate(self%nplenc)
      if (allocated(self%ntpenc)) deallocate(self%ntpenc)
      if (allocated(self%levelg)) deallocate(self%levelg)
      if (allocated(self%levelm)) deallocate(self%levelm)
      if (allocated(self%isperi)) deallocate(self%isperi)
      if (allocated(self%peri)) deallocate(self%peri)
      if (allocated(self%atp)) deallocate(self%atp)
      if (allocated(self%kin)) deallocate(self%kin)
      if (allocated(self%info)) deallocate(self%info)

      allocate(symba_particle_info :: self%info(n))
      allocate(self%lcollision(n))
      allocate(self%lencounter(n))
      allocate(self%lmtiny(n))
      allocate(self%nplenc(n))
      allocate(self%ntpenc(n))
      allocate(self%levelg(n))
      allocate(self%levelm(n))
      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))
      allocate(self%kin(n))

      self%lcollision(:) = .false.
      self%lencounter(:) = .false.
      self%lmtiny(:) = .false.
      self%nplenc(:) = 0
      self%ntpenc(:) = 0
      self%levelg(:) = -1
      self%levelm(:) = -1
      self%isperi(:) = 0
      self%peri(:) = 0.0_DP
      self%atp(:) = 0.0_DP
      self%kin(:)%nchild = 0
      self%kin(:)%parent = [(i, i=1, n)]
      return
   end subroutine symba_setup_pl


   module subroutine symba_setup_encounter(self, n)
      !! author: David A. Minton
      !!
      !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      !!
      implicit none
      ! Arguments
      class(symba_encounter), intent(inout) :: self !! SyMBA pl-tp encounter structure
      integer(I4B),         intent(in)    :: n    !! Number of encounters to allocate space for

      call setup_encounter(self, n)
      if (n == 0) return

      if (allocated(self%level)) deallocate(self%level)
      allocate(self%level(n))

      self%level(:) = -1

      return
   end subroutine symba_setup_encounter


   module subroutine symba_setup_tp(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(symba_tp),            intent(inout) :: self  !! SyMBA test particle object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class. In this case, helio_tp does not have its own setup method so we use the base method for swiftest_tp
      call setup_tp(self, n, param) 
      if (n <= 0) return

      if (allocated(self%nplenc)) deallocate(self%nplenc)
      if (allocated(self%levelg)) deallocate(self%levelg)
      if (allocated(self%levelm)) deallocate(self%levelm)
      if (allocated(self%info)) deallocate(self%info)

      allocate(self%nplenc(n))
      allocate(self%levelg(n))
      allocate(self%levelm(n))
      allocate(self%info(n))

      self%nplenc(:) = 0
      self%levelg(:) = -1
      self%levelm(:) = -1
      
      return
   end subroutine symba_setup_tp

end submodule s_symba_setup
