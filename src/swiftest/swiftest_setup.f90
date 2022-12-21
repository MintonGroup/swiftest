!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_setup
   use whm
   use rmvs
   use helio
   use symba
contains

   module subroutine swiftest_setup_construct_system(system, param)
      !! author: David A. Minton
      !!
      !! Constructor for a Swiftest nbody system. Creates the nbody system object based on the user-input integrator
      !! 
      implicit none
      ! Arguments
      class(swiftest_nbody_system), allocatable, intent(inout) :: system !! Swiftest system object
      class(swiftest_parameters),                intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      type(encounter_storage)  :: encounter_history
      type(collision_storage)  :: collision_history

      select type(param)
      class is (swiftest_parameters)
         allocate(swiftest_storage(param%dump_cadence) :: param%system_history)
         allocate(swiftest_netcdf_parameters :: param%system_history%nc)
         call param%system_history%reset()

         select case(param%integrator)
         case (INT_BS)
            write(*,*) 'Bulirsch-Stoer integrator not yet enabled'
         case (INT_HELIO)
            allocate(helio_nbody_system :: system)
            select type(system)
            class is (helio_nbody_system)
               allocate(helio_cb :: system%cb)
               allocate(helio_pl :: system%pl)
               allocate(helio_tp :: system%tp)
               allocate(helio_tp :: system%tp_discards)
            end select
         case (INT_RA15)
            write(*,*) 'Radau integrator not yet enabled'
         case (INT_TU4)
            write(*,*) 'INT_TU4 integrator not yet enabled'
         case (INT_WHM)
            allocate(whm_nbody_system :: system)
            select type(system)
            class is (whm_nbody_system)
               allocate(whm_cb :: system%cb)
               allocate(whm_pl :: system%pl)
               allocate(whm_tp :: system%tp)
               allocate(whm_tp :: system%tp_discards)
            end select
         case (INT_RMVS)
            allocate(rmvs_nbody_system :: system)
            select type(system)
            class is (rmvs_nbody_system)
               allocate(rmvs_cb :: system%cb)
               allocate(rmvs_pl :: system%pl)
               allocate(rmvs_tp :: system%tp)
               allocate(rmvs_tp :: system%tp_discards)
            end select
         case (INT_SYMBA)
            allocate(symba_nbody_system :: system)
            select type(system)
            class is (symba_nbody_system)
               allocate(symba_cb :: system%cb)
               allocate(symba_pl :: system%pl)
               allocate(symba_tp :: system%tp)

               allocate(symba_tp :: system%tp_discards)
               allocate(symba_pl :: system%pl_adds)
               allocate(symba_pl :: system%pl_discards)

               allocate(collision_list_pltp :: system%pltp_encounter)
               allocate(collision_list_plpl :: system%plpl_encounter)
               allocate(collision_list_plpl :: system%plpl_collision)

               if (param%lfragmentation) then
                  allocate(fraggle_system :: system%collision_system)
                  call system%collision_system%setup(system)
               end if

               if (param%lenc_save_trajectory .or. param%lenc_save_closest) then
                  allocate(encounter_netcdf_parameters :: encounter_history%nc)
                  call encounter_history%reset()
                  select type(nc => encounter_history%nc)
                  class is (encounter_netcdf_parameters)
                     nc%file_number = param%iloop / param%dump_cadence
                  end select
                  allocate(system%encounter_history, source=encounter_history)
               end if
               
               allocate(collision_netcdf_parameters :: collision_history%nc)
               call collision_history%reset()
               select type(nc => collision_history%nc)
               class is (collision_netcdf_parameters)
                  nc%file_number = param%iloop / param%dump_cadence
               end select
               allocate(system%collision_history, source=collision_history)

            end select
         case (INT_RINGMOONS)
            write(*,*) 'RINGMOONS-SyMBA integrator not yet enabled'
         case default
            write(*,*) 'Unkown integrator',param%integrator
            call util_exit(FAILURE)
         end select

         allocate(swiftest_particle_info :: system%cb%info)
      end select

      return
   end subroutine swiftest_setup_construct_system


   module subroutine swiftest_setup_initialize_particle_info_system(self, param)
      !! author: David A. Minton
      !!
      !! Setup up particle information metadata from initial conditions
      !
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i

      associate(pl => self%pl, npl => self%pl%nbody, tp => self%tp, ntp => self%tp%nbody)

         if (.not. allocated(self%cb%info)) allocate(swiftest_particle_info :: self%cb%info)

         call self%cb%info%set_value(particle_type=CB_TYPE_NAME, status="ACTIVE", origin_type="Initial conditions", &
                                origin_time=param%t0, origin_rh=[0.0_DP, 0.0_DP, 0.0_DP], origin_vh=[0.0_DP, 0.0_DP, 0.0_DP])
         do i = 1, self%pl%nbody
            call pl%info(i)%set_value(particle_type=PL_TYPE_NAME, status="ACTIVE", origin_type="Initial conditions", &
                                       origin_time=param%t0, origin_rh=self%pl%rh(:,i), origin_vh=self%pl%vh(:,i))
         end do
         do i = 1, self%tp%nbody
            call tp%info(i)%set_value(particle_type=TP_TYPE_NAME, status="ACTIVE", origin_type="Initial conditions", &
                                      origin_time=param%t0, origin_rh=self%tp%rh(:,i), origin_vh=self%tp%vh(:,i))
         end do

      end associate

      return
   end subroutine swiftest_setup_initialize_particle_info_system


   module subroutine swiftest_setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters

      associate(system => self, cb => self%cb, pl => self%pl, tp => self%tp)

         call system%read_in(param)
         call system%validate_ids(param)
         call system%set_msys()
         call pl%set_mu(cb) 
         call tp%set_mu(cb) 
         if (param%in_form == "EL") then
            call pl%el2xv(cb)
            call tp%el2xv(cb)
         end if
         call pl%flatten(param)
         if (.not.param%lrhill_present) call pl%set_rhill(cb)
         pl%lfirst = param%lfirstkick
         tp%lfirst = param%lfirstkick

         if (.not.param%lrestart) then
            call system%init_particle_info(param)
         end if
      end associate

      return
   end subroutine swiftest_setup_initialize_system


   module subroutine swiftest_setup_body(self, n, param)
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest particle class. Allocates space for all particles and
      !! initializes all components with a value.
      !! Note: Timing tests indicate that (NDIM, n) is more efficient than (NDIM, n) 
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self  !! Swiftest generic body object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter
      ! Internals
      integer(I4B) :: i

      if (n < 0) return

      self%lfirst = .true.

      call self%dealloc()

      self%nbody = n
      if (n == 0) return

      allocate(swiftest_particle_info :: self%info(n))
      allocate(self%id(n))
      allocate(self%status(n))
      allocate(self%ldiscard(n))
      allocate(self%lmask(n))
      allocate(self%mu(n))
      allocate(self%rh(NDIM, n))
      allocate(self%vh(NDIM, n))
      allocate(self%rb(NDIM, n))
      allocate(self%vb(NDIM, n))
      allocate(self%ah(NDIM, n))
      allocate(self%ir3h(n))
      allocate(self%aobl(NDIM, n))
      if (param%lclose) then
         allocate(self%lcollision(n))
         allocate(self%lencounter(n))
         self%lcollision(:) = .false.
         self%lencounter(:) = .false.
      end if

      self%id(:) = 0
      do i = 1, n
         call self%info(i)%set_value(&
            name = "UNNAMED", &
            particle_type = "UNKNOWN", &
            status = "INACTIVE", & 
            origin_type = "UNKNOWN", &
            collision_id = 0, &
            origin_time = -huge(1.0_DP), & 
            origin_rh = [0.0_DP, 0.0_DP, 0.0_DP], &
            origin_vh = [0.0_DP, 0.0_DP, 0.0_DP], &
            discard_time = huge(1.0_DP), & 
            discard_rh = [0.0_DP, 0.0_DP, 0.0_DP], &
            discard_vh = [0.0_DP, 0.0_DP, 0.0_DP], &
            discard_body_id = -1  &
         )
      end do

      self%status(:) = INACTIVE
      self%ldiscard(:) = .false.
      self%lmask(:)  = .false.
      self%mu(:)     = 0.0_DP
      self%rh(:,:)   = 0.0_DP
      self%vh(:,:)   = 0.0_DP
      self%rb(:,:)   = 0.0_DP
      self%vb(:,:)   = 0.0_DP
      self%ah(:,:)   = 0.0_DP
      self%ir3h(:)   = 0.0_DP
      self%aobl(:,:) = 0.0_DP

      if (param%ltides) then
         allocate(self%atide(NDIM, n))
         self%atide(:,:) = 0.0_DP
      end if
      if (param%lgr) then
         allocate(self%agr(NDIM, n))
         self%agr(:,:) = 0.0_DP
      end if

      return
   end subroutine swiftest_setup_body


   module subroutine swiftest_setup_pl(self, n, param)
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest massive body class. Allocates space for all particles and
      !! initializes all components with a value. 
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class
      !> The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call swiftest_setup_body(self, n, param)
      if (n == 0) return

      allocate(self%mass(n))
      allocate(self%Gmass(n))
      allocate(self%rhill(n))
      allocate(self%renc(n))

      self%mass(:) = 0.0_DP
      self%Gmass(:) = 0.0_DP
      self%rhill(:) = 0.0_DP
      self%renc(:) = 0.0_DP

      self%nplpl = 0   

      if (param%lclose) then
         allocate(self%nplenc(n))
         allocate(self%ntpenc(n))
         allocate(self%radius(n))
         allocate(self%density(n))

         self%nplenc(:) = 0
         self%ntpenc(:) = 0
         self%radius(:) = 0.0_DP
         self%density(:) = 1.0_DP

      end if

      if (param%lmtiny_pl) then
         allocate(self%lmtiny(n))
         self%lmtiny(:) = .false.
      end if

      if (param%lrotation) then
         allocate(self%rot(NDIM, n))
         allocate(self%Ip(NDIM, n))
         self%rot(:,:) = 0.0_DP
         self%Ip(:,:) = 0.0_DP
      end if

      if (param%ltides) then
         allocate(self%k2(n))
         allocate(self%Q(n))
         allocate(self%tlag(n))
         self%k2(:) = 0.0_DP
         self%Q(:) = 0.0_DP
         self%tlag(:) = 0.0_DP
      end if
      
      return
   end subroutine swiftest_setup_pl
   

   module subroutine swiftest_setup_tp(self, n, param)
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest test particle particle class. Allocates space for 
      !! all particles and initializes all components with a value. 
      implicit none
      ! Arguments
      class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class
      !> The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call swiftest_setup_body(self, n, param)
      if (n == 0) return

      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))
      allocate(self%nplenc(n))

      self%isperi(:) = 0
      self%peri(:)   = 0.0_DP
      self%atp(:)    = 0.0_DP
      self%nplenc(:) = 0

      return
   end subroutine swiftest_setup_tp

end submodule s_setup
