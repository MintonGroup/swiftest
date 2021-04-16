submodule(rmvs_classes) s_rmvs_setup
   use swiftest
contains
   module subroutine rmvs_setup_pl(self,n)
      !! author: David A. Minton
      !!
      !! Allocate RMVS test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine rmvs_setup.f90
      implicit none
      ! Arguments
      class(rmvs_base_pl),           intent(inout) :: self !! RMVS test particle object
      integer(I4B),                  intent(in)    :: n    !! Number of test particles to allocate
      ! Internals
      integer(I4B)                                 :: i,j
      !type(swiftest_configuration)                 :: encounter_config

      !> Call allocation method for parent class
      associate(system => self)
         call whm_setup_pl(system, n) 
         if (n <= 0) return

         if (.not.system%lplanetocentric) then
            allocate(self%nenc(n))
            self%nenc(:)         = 0

            ! Set up inner and outer planet interpolation vector storage containers
            allocate(self%outer(0:NTENC))
            do i = 0, NTENC
               allocate(self%outer(i)%x(NDIM, n))
               allocate(self%outer(i)%v(NDIM, n))
            end do
            allocate(self%inner(0:NTPHENC))
            do i = 0, NTPHENC
               allocate(self%inner(i)%x(NDIM, n))
               allocate(self%inner(i)%v(NDIM, n))
               allocate(self%inner(i)%aobl(NDIM, n))
            end do
         end if
      end associate
      return
   end subroutine rmvs_setup_pl 

   module subroutine rmvs_setup_tp(self,n)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(rmvs_tp),              intent(inout)   :: self !! RMVS test particle object
      integer,                     intent(in)      :: n    !! Number of test particles to allocate

      !> Call allocation method for parent class
      call whm_setup_tp(self, n) 
      if (n <= 0) return

      allocate(self%lperi(n))
      allocate(self%plperP(n))
      allocate(self%plencP(n))
      if (self%lplanetocentric) then
         allocate(self%xheliocentric(NDIM, n))
      end if

      self%lperi(:)  = .false.

      return
   end subroutine rmvs_setup_tp

   module subroutine rmvs_setup_system(self, config)
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files.
      !! 
      !! We currently rearrange the pl order to keep it consistent with the way Swifter does it 
      !! In Swifter, the central body occupies the first position in the pl list, and during
      !! encounters, the encountering planet is skipped in loops. In Swiftest, we instantiate an
      !! RMVS nbody system object attached to each pl to store planetocentric versions of the system
      !! to use during close encounters. 
      implicit none
      ! Arguments
      class(rmvs_nbody_system),      intent(inout) :: self    !! RMVS system object
      class(swiftest_configuration), intent(inout) :: config  !! Input collection of  configuration parameters 
      ! Internals
      integer(I4B) :: i, j
      ! Call parent method
      call whm_setup_system(self, config)

      ! Set up the pl-tp planetocentric encounter structures for pl and cb. The planetocentric tp structures are 
      ! generated as necessary during close encounter steps.
      select type(pl => self%pl)
      class is(rmvs_pl)
      select type(cb => self%cb)
      class is (rmvs_cb)
      select type (tp => self%tp)
      class is (rmvs_tp)
         tp%cb_heliocentric = cb
         pl%lplanetocentric = .false.
         tp%lplanetocentric = .false.
         cb%lplanetocentric = .false.
         associate(npl => pl%nbody)
            allocate(pl%planetocentric(npl))
            do i = 1, npl
               allocate(pl%planetocentric(i)%cb, source=cb)
               allocate(rmvs_pl :: pl%planetocentric(i)%pl)
               associate(plenci => pl%planetocentric(i)%pl, cbenci => pl%planetocentric(i)%cb)
                  cbenci%lplanetocentric = .true.
                  plenci%lplanetocentric = .true.
                  call plenci%setup(npl)
                  plenci%status(:) = ACTIVE

                  ! plind stores the heliocentric index value of a planetocentric planet
                  ! e.g. Consider an encounter with planet 3.  
                  ! Then the following will be the values of plind:
                  ! pl%planetocentric(3)%pl%plind(1) = 0 (central body - never used)  
                  ! pl%planetocentric(3)%pl%plind(2) = 1  
                  ! pl%planetocentric(3)%pl%plind(3) = 2
                  ! pl%planetocentric(3)%pl%plind(4) = 4
                  ! pl%planetocentric(3)%pl%plind(5) = 5
                  ! etc.  
                  allocate(plenci%plind(npl))
                  plenci%plind(1:npl) = [(j,j=1,npl)] 
                  plenci%plind(2:npl) = pack(plenci%plind(1:npl), plenci%plind(1:npl) /= i)
                  plenci%plind(1)     = 0
                  plenci%Gmass(1)     = cb%Gmass
                  plenci%Gmass(2:npl) = pl%Gmass(plenci%plind(2:npl))
                  cbenci%Gmass        = pl%Gmass(i)
               end associate
            end do
         end associate
      end select
      end select
      end select
   
   end subroutine rmvs_setup_system
   
   module subroutine rmvs_setup_set_beg_end(self, xbeg, xend, vbeg)
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of xbeg, xend, and vbeg
      implicit none
      ! Arguments
      class(rmvs_tp),           intent(inout)          :: self !! RMVS test particle object
      real(DP), dimension(:,:), intent(in),   optional :: xbeg, xend, vbeg

      if (present(xbeg)) then
         if (allocated(self%xbeg)) deallocate(self%xbeg)
         allocate(self%xbeg, source=xbeg)
      end if
      if (present(xend)) then
         if (allocated(self%xend)) deallocate(self%xend)
         allocate(self%xend, source=xend)
      end if
      if (present(vbeg)) then
         if (allocated(self%vbeg)) deallocate(self%vbeg)
         allocate(self%vbeg, source=vbeg)
      end if

      return

   end subroutine rmvs_setup_set_beg_end

end submodule s_rmvs_setup
