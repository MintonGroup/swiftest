!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(rmvs_classes) s_rmvs_setup
   use swiftest
contains

   module subroutine rmvs_setup_pl(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate RMVS test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine rmvs_setup.f90
      implicit none
      ! Arguments
      class(rmvs_pl),            intent(inout) :: self  !! RMVS test particle object
      integer(I4B),              intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter
      ! Internals
      integer(I4B) :: i

      !> Call allocation method for parent class
      associate(pl => self)
         call whm_setup_pl(pl, n, param) 
         if (n == 0) return

         allocate(pl%outer(0:NTENC))
         allocate(pl%inner(0:NTPHENC))
         if (.not.pl%lplanetocentric) then
            allocate(pl%nenc(n))
            pl%nenc(:) = 0
            ! Set up inner and outer planet interpolation vector storage containers
            do i = 0, NTENC
               allocate(pl%outer(i)%x(NDIM, n))
               allocate(pl%outer(i)%v(NDIM, n))
               pl%outer(i)%x(:,:) = 0.0_DP
               pl%outer(i)%v(:,:) = 0.0_DP
            end do
            do i = 0, NTPHENC
               allocate(pl%inner(i)%x(NDIM, n))
               allocate(pl%inner(i)%v(NDIM, n))
               allocate(pl%inner(i)%aobl(NDIM, n))
               pl%inner(i)%x(:,:) = 0.0_DP
               pl%inner(i)%v(:,:) = 0.0_DP
               pl%inner(i)%aobl(:,:) = 0.0_DP
            end do
            ! if (param%ltides) then
            !    do i = 0, NTPHENC
            !       allocate(pl%inner(i)%atide(NDIM, n))
            !       pl%inner(i)%atide(:,:) = 0.0_DP
            !    end do
            ! end if
         end if
      end associate
      return
   end subroutine rmvs_setup_pl 


   module subroutine rmvs_setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize an RMVS nbody system from files and sets up the planetocentric structures.
      !! 
      !! We currently rearrange the pl order to keep it consistent with the way Swifter does it 
      !! In Swifter, the central body occupies the first position in the pl list, and during
      !! encounters, the encountering planet is skipped in loops. In Swiftest, we instantiate an
      !! RMVS nbody system object attached to each pl to store planetocentric versions of the system
      !! to use during close encounters. 
      implicit none
      ! Arguments
      class(rmvs_nbody_system),   intent(inout) :: self    !! RMVS system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i, j

      ! Call parent method
      call whm_setup_initialize_system(self, param)

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
                  pl%planetocentric(:)%lplanetocentric = .true.
                  do i = 1, npl
                     allocate(pl%planetocentric(i)%cb, source=cb)
                     allocate(rmvs_pl :: pl%planetocentric(i)%pl)
                     select type(cbenci => pl%planetocentric(i)%cb)
                     class is (rmvs_cb)
                        select type(plenci => pl%planetocentric(i)%pl)
                        class is (rmvs_pl)
                           cbenci%lplanetocentric = .true.
                           plenci%lplanetocentric = .true.
                           call plenci%setup(npl, param)
                           plenci%status(:) = ACTIVE
                           plenci%lmask(:) = .true.
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
                        end select
                     end select
                  end do
               end associate
            end select
         end select
      end select
      return
   end subroutine rmvs_setup_initialize_system


   module subroutine rmvs_setup_tp(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(rmvs_tp),             intent(inout) :: self  !! RMVS test particle object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class. In this case, whm does not have its own setup method, so we use the base method for swiftest_tp
      call setup_tp(self, n, param) 
      if (n <= 0) return

      allocate(self%lperi(n))
      allocate(self%plperP(n))
      allocate(self%plencP(n))

      if (self%lplanetocentric) allocate(self%xheliocentric(NDIM, n))

      self%lperi(:)  = .false.

      return
   end subroutine rmvs_setup_tp

end submodule s_rmvs_setup
