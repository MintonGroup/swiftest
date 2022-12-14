!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(swiftest_classes) s_util_set
   !! author: David A. Minton
   !! This submodule contains a collection of setter method implementations
   use swiftest
contains

   module subroutine util_set_beg_end_pl(self, rbeg, xend, vbeg)
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of rbeg, xend, and vbeg
      implicit none
      ! Arguments
      class(swiftest_pl),       intent(inout)          :: self !! Swiftest massive body object
      real(DP), dimension(:,:), intent(in),   optional :: rbeg, xend, vbeg

      if (present(rbeg)) then
         if (allocated(self%rbeg)) deallocate(self%rbeg)
         allocate(self%rbeg, source=rbeg)
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
   end subroutine util_set_beg_end_pl


   module subroutine util_set_ir3h(self)
      !! author: David A. Minton
      !!
      !! Sets the inverse heliocentric radius term (1/rh**3) for all bodies in a structure
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self !! Swiftest generic body object
      ! Internals
      integer(I4B) :: i
      real(DP) :: r2, irh

      if (self%nbody > 0) then

         do i = 1, self%nbody
            r2 = dot_product(self%rh(:, i), self%rh(:, i))
            irh = 1.0_DP / sqrt(r2)
            self%ir3h(i) = irh / r2
         end do
      end if

      return
   end subroutine util_set_ir3h


   module subroutine util_set_msys(self)
      !! author: David A. Minton
      !!
      !! Sets the value of msys and the vector mass quantities based on the total mass of the system
      implicit none
      ! Arguments
      class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest nobdy system object

      self%Gmtot = self%cb%Gmass + sum(self%pl%Gmass(1:self%pl%nbody), self%pl%status(1:self%pl%nbody) /= INACTIVE)

      return
   end subroutine util_set_msys


   module subroutine util_set_mu_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Computes G * (M + m) for each massive body
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object

      if (self%nbody > 0) self%mu(1:self%nbody) = cb%Gmass + self%Gmass(1:self%nbody)

      return
   end subroutine util_set_mu_pl


   module subroutine util_set_mu_tp(self, cb)
      !! author: David A. Minton
      !!
      !! Converts certain scalar values to arrays so that they can be used in elemental functions
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: self !! Swiftest test particle object
      class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object

      if (self%nbody == 0) return
      self%mu(1:self%nbody) = cb%Gmass

      return
   end subroutine util_set_mu_tp

   module subroutine util_set_particle_info(self, name, particle_type, status, origin_type, origin_time, collision_id, origin_rh,&
                                            origin_vh, discard_time, discard_rh, discard_vh, discard_body_id)
      !! author: David A. Minton
      !!
      !! Sets one or more values of the particle information metadata object
      implicit none
      ! Arguments
      class(swiftest_particle_info), intent(inout)           :: self
      character(len=*),              intent(in),    optional :: name            !! Non-unique name
      character(len=*),              intent(in),    optional :: particle_type   !! String containing a description of the particle type (e.g. Central Body, Massive Body, Test Particle)
      character(len=*),              intent(in),    optional :: status          !! Particle status description: ACTIVE, MERGED, FRAGMENTED, etc.
      character(len=*),              intent(in),    optional :: origin_type     !! String containing a description of the origin of the particle (e.g. Initial Conditions, Supercatastrophic, Disruption, etc.)
      real(DP),                      intent(in),    optional :: origin_time     !! The time of the particle's formation
      integer(I4B),                  intent(in),    optional :: collision_id    !! The ID fo the collision that formed the particle
      real(DP), dimension(:),        intent(in),    optional :: origin_rh       !! The heliocentric distance vector at the time of the particle's formation
      real(DP), dimension(:),        intent(in),    optional :: origin_vh       !! The heliocentric velocity vector at the time of the particle's formation
      real(DP),                      intent(in),    optional :: discard_time    !! The time of the particle's discard
      real(DP), dimension(:),        intent(in),    optional :: discard_rh      !! The heliocentric distance vector at the time of the particle's discard
      real(DP), dimension(:),        intent(in),    optional :: discard_vh      !! The heliocentric velocity vector at the time of the particle's discard
      integer(I4B),                  intent(in),    optional :: discard_body_id !! The id of the other body involved in the discard (0 if no other body involved)
      ! Internals
      character(len=NAMELEN) :: lenstr
      character(len=:), allocatable :: fmtlabel

      write(lenstr, *) NAMELEN
      fmtlabel = "(A" // trim(adjustl(lenstr)) // ")"

      if (present(name)) then
         write(self%name, fmtlabel) trim(adjustl(name))
      end if
      if (present(particle_type)) then
         write(self%particle_type, fmtlabel) trim(adjustl(particle_type))
      end if 
      if (present(status)) then
         write(self%status, fmtlabel) trim(adjustl(status))
      end if
      if (present(origin_type)) then
         write(self%origin_type, fmtlabel) trim(adjustl(origin_type))
      end if
      if (present(origin_time)) then
         self%origin_time = origin_time
      end if
      if (present(collision_id)) then
         self%collision_id = collision_id
      end if
      if (present(origin_rh)) then
         self%origin_rh(:) = origin_rh(:)
      end if
      if (present(origin_vh)) then
         self%origin_vh(:) = origin_vh(:)
      end if
      if (present(discard_time)) then
         self%discard_time = discard_time
      end if
      if (present(discard_rh)) then
         self%discard_rh(:) = discard_rh(:)
      end if
      if (present(discard_vh)) then
         self%discard_vh(:) = discard_vh(:)
      end if
      if (present(discard_body_id)) then
         self%discard_body_id = discard_body_id
      end if

      return
   end subroutine util_set_particle_info


   module subroutine util_set_renc_I4B(self, scale)
      !! author: David A. Minton
      !!
      !! Sets the critical radius for encounter given an input scale factor
      !!
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      integer(I4B),       intent(in)    :: scale !! Input scale factor (multiplier of Hill's sphere size)

      associate(pl => self, npl => self%nbody)
         pl%renc(1:npl) = pl%rhill(1:npl) * scale
      end associate

      return
   end subroutine util_set_renc_I4B


   module subroutine util_set_renc_DP(self, scale)
      !! author: David A. Minton
      !!
      !! Sets the critical radius for encounter given an input scale factor
      !!
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      real(DP),           intent(in)    :: scale !! Input scale factor (multiplier of Hill's sphere size)

      associate(pl => self, npl => self%nbody)
         pl%renc(1:npl) = pl%rhill(1:npl) * scale
      end associate

      return
   end subroutine util_set_renc_DP


   module subroutine util_set_rhill(self,cb)
      !! author: David A. Minton
      !!
      !! Sets the value of the Hill's radius
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object

      if (self%nbody == 0) return

      call self%xv2el(cb) 
      self%rhill(1:self%nbody) = self%a(1:self%nbody) * (self%Gmass(1:self%nbody) / cb%Gmass / 3)**THIRD 

      return
   end subroutine util_set_rhill


   module subroutine util_set_rhill_approximate(self,cb)
      !! author: David A. Minton
      !!
      !! Sets the approximate value of the Hill's radius using the heliocentric radius instead of computing the semimajor axis
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      real(DP), dimension(:), allocatable :: rh

      if (self%nbody == 0) return

      rh(1:self%nbody) = .mag. self%rh(:,1:self%nbody)
      self%rhill(1:self%nbody) = rh(1:self%nbody) * (self%Gmass(1:self%nbody) / cb%Gmass / 3)**THIRD 

      return
   end subroutine util_set_rhill_approximate

end submodule s_util_set