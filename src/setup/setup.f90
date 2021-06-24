submodule (swiftest_classes) s_setup
   use swiftest
contains
   module subroutine setup_construct_system(system, param)
      !! author: David A. Minton
      !!
      !! Constructor for a Swiftest nbody system. Creates the nbody system object based on the user-input integrator
      !! 
      implicit none
      ! Arguments
      class(swiftest_nbody_system),  allocatable,  intent(inout) :: system     !! Swiftest system object
      type(swiftest_parameters),                intent(in)    :: param     !! Swiftest parameters parameters

      select case(param%integrator)
      case (BS)
         write(*,*) 'Bulirsch-Stoer integrator not yet enabled'
      case (HELIO)
         write(*,*) 'Democratic Heliocentric integrator not yet enabled'
         allocate(helio_nbody_system :: system)
         select type(system)
         class is (helio_nbody_system)
            allocate(helio_cb :: system%cb)
            allocate(helio_pl :: system%pl)
            allocate(helio_tp :: system%tp)
            allocate(helio_tp :: system%tp_discards)
         end select
      case (RA15)
         write(*,*) 'Radau integrator not yet enabled'
      case (TU4)
         write(*,*) 'TU4 integrator not yet enabled'
      case (WHM)
         allocate(whm_nbody_system :: system)
         select type(system)
         class is (whm_nbody_system)
            allocate(whm_cb :: system%cb)
            allocate(whm_pl :: system%pl)
            allocate(whm_tp :: system%tp)
            allocate(whm_tp :: system%tp_discards)
         end select
      case (RMVS)
         allocate(rmvs_nbody_system :: system)
         select type(system)
         class is (rmvs_nbody_system)
            allocate(rmvs_cb :: system%cb)
            allocate(rmvs_pl :: system%pl)
            allocate(rmvs_tp :: system%tp)
            allocate(rmvs_tp :: system%tp_discards)
         end select
      case (SYMBA)
         write(*,*) 'SyMBA integrator not yet enabled'
      case (RINGMOONS)
         write(*,*) 'RINGMOONS-SyMBA integrator not yet enabled'
      case default
         write(*,*) 'Unkown integrator',param%integrator
         call util_exit(FAILURE)
      end select

      return
   end subroutine setup_construct_system

   module procedure setup_body
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest particle class. Allocates space for all particles and
      !! initializes all components with a value.
      !! Note: Timing tests indicate that (NDIM, n) is more efficient than (NDIM, n) 
      implicit none

      self%nbody = n
      if (n <= 0) return

      !write(*,*) 'Allocating the basic Swiftest particle'
      allocate(self%name(n))
      allocate(self%status(n))
      allocate(self%ldiscard(n))
      allocate(self%xh(NDIM, n))
      allocate(self%vh(NDIM, n))
      allocate(self%xb(NDIM, n))
      allocate(self%vb(NDIM, n))
      allocate(self%ah(NDIM, n))
      allocate(self%aobl(NDIM, n))
      allocate(self%ir3h(n))
      allocate(self%a(n))
      allocate(self%e(n))
      allocate(self%inc(n))
      allocate(self%capom(n))
      allocate(self%omega(n))
      allocate(self%capm(n))
      allocate(self%mu(n))

      self%name(:)   = 0
      self%status(:) = INACTIVE
      self%ldiscard(:) = .false.
      self%xh(:,:)   = 0.0_DP
      self%vh(:,:)   = 0.0_DP
      self%xb(:,:)   = 0.0_DP
      self%vb(:,:)   = 0.0_DP
      self%ah(:,:)   = 0.0_DP
      self%aobl(:,:) = 0.0_DP
      self%ir3h(:)   = 0.0_DP
      self%a(:)      = 0.0_DP
      self%e(:)      = 0.0_DP
      self%inc(:)    = 0.0_DP
      self%capom(:)  = 0.0_DP
      self%omega(:)  = 0.0_DP
      self%capm(:)   = 0.0_DP
      self%a(:)      = 0.0_DP
      self%mu(:)     = 0.0_DP

      return
   end procedure setup_body

   module procedure setup_pl
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest massive body class. Allocates space for all particles and
      !! initializes all components with a value. 
      implicit none

      !> Call allocation method for parent class
      !> The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call setup_body(self, n)
      if (n <= 0) return 

      allocate(self%mass(n))
      allocate(self%Gmass(n))
      allocate(self%rhill(n))
      allocate(self%radius(n))
      allocate(self%density(n))
      allocate(self%Ip(NDIM, n))
      allocate(self%rot(NDIM, n))
      allocate(self%k2(n))
      allocate(self%Q(n))

      self%mass(:) = 0.0_DP
      self%Gmass(:) = 0.0_DP
      self%rhill(:) = 0.0_DP
      self%radius(:) = 0.0_DP
      self%density(:) = 0.0_DP
      self%Ip(:,:) = 0.0_DP
      self%rot(:,:) = 0.0_DP
      self%k2(:) = 0.0_DP
      self%Q(:) = 0.0_DP
      self%num_comparisons = 0   
      return
   end procedure setup_pl

   module procedure setup_tp
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest test particle particle class. Allocates space for 
      !! all particles and initializes all components with a value. 
      implicit none

      !> Call allocation method for parent class
      !> The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call setup_body(self, n)
      if (n <= 0) return

      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))

      self%isperi(:) = 0
      self%peri(:)   = 0.0_DP
      self%atp(:)    = 0.0_DP

      return
   end procedure setup_tp

   module procedure setup_set_msys
      !! author: David A. Minton
      !!
      !! Sets the value of msys and the vector mass quantities based on the total mass of the system
      self%msys = self%cb%mass + sum(self%pl%mass(1:self%pl%nbody))

      return
   end procedure setup_set_msys

   module procedure setup_set_mu_pl
      !! author: David A. Minton
      !!
      !! Computes G * (M + m) for each massive body
      implicit none

      if (self%nbody > 0) self%mu(:) = cb%Gmass + self%Gmass(:)

      return
   end procedure setup_set_mu_pl

   module procedure setup_set_mu_tp
      !! author: David A. Minton
      !!
      !! Converts certain scalar values to arrays so that they can be used in elemental functions
      implicit none

      if (self%nbody > 0) self%mu(:) = cb%Gmass

      return
   end procedure setup_set_mu_tp

   module procedure setup_set_rhill
      !! author: David A. Minton
      !!
      !! Sets the value of the Hill's radius
      implicit none


      if (self%nbody > 0) then
         call self%xv2el(cb) 
         self%rhill(:) = self%a(:) * (self%Gmass(:) / cb%Gmass / 3)**THIRD 
      end if

      return
   end procedure setup_set_rhill

   module procedure setup_set_ir3h  
      !! author: David A. Minton
      !!
      !! Sets the inverse heliocentric radius term (1/rh**3) for all bodies in a structure
      implicit none

      integer(I4B) :: i
      real(DP) :: r2, irh

      if (self%nbody > 0) then

         do i = 1, self%nbody
            r2 = dot_product(self%xh(:, i), self%xh(:, i))
            irh = 1.0_DP / sqrt(r2)
            self%ir3h(i) = irh / r2
         end do
      end if

      return
   end procedure setup_set_ir3h

end submodule s_setup
