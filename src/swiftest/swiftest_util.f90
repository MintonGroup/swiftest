! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_swiftest_util
   use whm
   use rmvs
   use helio
   use symba
   use fraggle
contains


   module subroutine swiftest_util_append_arr_info(arr, source, nold, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of particle information type onto another. If the destination array is not allocated, or is not big 
      !! enough, this will allocate space for it.
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr  
         !! Destination array 
      type(swiftest_particle_info), dimension(:), allocatable, intent(in) :: source 
         !! Array to append 
      integer(I4B), intent(in), optional :: nold 
         !! Extent of original array. If passed the source array will begin at arr(nold+1). Otherwise, the size of arr will be used.
      logical, dimension(:), intent(in), optional :: lsource_mask 
         !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew, nsrc, nend_orig, i
      integer(I4B), dimension(:), allocatable :: idx

      if (.not.allocated(source)) return

      if (present(lsource_mask)) then
         nsrc = count(lsource_mask(:))
      else
         nsrc = size(source)
      end if
      if (nsrc == 0) return

      if (.not.allocated(arr)) then
         nend_orig = 0
         allocate(arr(nsrc))
      else
         if (present(nold)) then
            nend_orig = nold
         else
            nend_orig = size(arr)
         end if
         call util_resize(arr, nend_orig + nsrc)
      end if
      nnew = nend_orig + nsrc

      allocate(idx(nsrc))
      if (present(lsource_mask)) then
         idx = pack([(i, i = 1, size(lsource_mask))], lsource_mask(:))
      else
         idx = [(i, i = 1,nsrc)]
      end  if

      call swiftest_util_copy_particle_info_arr(source(:), arr(nend_orig+1:nnew), idx)

      return
   end subroutine swiftest_util_append_arr_info


   module subroutine swiftest_util_append_arr_kin(arr, source, nold, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of kinship type onto another. If the destination array is not allocated, or is not big enough, this 
      !! will allocate space for it.
      implicit none
      ! Arguments
      type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: arr 
         !! Destination array 
      type(swiftest_kinship), dimension(:), allocatable, intent(in) :: source 
         !! Array to append 
      integer(I4B), intent(in), optional :: nold 
         !! Extent of original array. If passed the source array will begin at arr(nold+1). Otherwise, the size of arr will be used.
      logical, dimension(:), intent(in), optional :: lsource_mask 
         !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew, nsrc, nend_orig

      if (.not.allocated(source)) return

      if (present(lsource_mask)) then
         nsrc = count(lsource_mask(:))
      else
         nsrc = size(source)
      end if
      if (nsrc == 0) return

      if (.not.allocated(arr)) then
         nend_orig = 0
         allocate(arr(nsrc))
      else
         if (present(nold)) then
            nend_orig = nold
         else
            nend_orig = size(arr)
         end if
         call util_resize(arr, nend_orig + nsrc)
      end if
      nnew = nend_orig + nsrc

      if (present(lsource_mask)) then
         arr(nend_orig + 1:nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))
      else
         arr(nend_orig + 1:nnew) = source(1:nsrc)
      end if

      return
   end subroutine swiftest_util_append_arr_kin


   module subroutine swiftest_util_append_body(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),  intent(inout) :: self  
         !! Swiftest body object
      class(swiftest_body),  intent(in) :: source  
         !! Source object to append
      logical, dimension(:), intent(in) :: lsource_mask 
         !! Logical mask indicating which elements to append to

      call util_append(self%id, source%id, lsource_mask=lsource_mask)
      call util_append(self%info, source%info, lsource_mask=lsource_mask)
      call util_append(self%lmask, source%lmask, lsource_mask=lsource_mask)
      call util_append(self%status, source%status, lsource_mask=lsource_mask)
      call util_append(self%ldiscard, source%ldiscard, lsource_mask=lsource_mask)
      call util_append(self%lencounter, source%lencounter, lsource_mask=lsource_mask)
      call util_append(self%lcollision, source%lcollision, lsource_mask=lsource_mask)
      call util_append(self%mu, source%mu, lsource_mask=lsource_mask)
      call util_append(self%rh, source%rh, lsource_mask=lsource_mask)
      call util_append(self%vh, source%vh, lsource_mask=lsource_mask)
      call util_append(self%rb, source%rb, lsource_mask=lsource_mask)
      call util_append(self%vb, source%vb, lsource_mask=lsource_mask)
      call util_append(self%ah, source%ah, lsource_mask=lsource_mask)
      call util_append(self%aobl, source%aobl, lsource_mask=lsource_mask)
      call util_append(self%atide, source%atide, lsource_mask=lsource_mask)
      call util_append(self%agr, source%agr, lsource_mask=lsource_mask)
      call util_append(self%ir3h, source%ir3h, lsource_mask=lsource_mask)
      call util_append(self%isperi, source%isperi, lsource_mask=lsource_mask)
      call util_append(self%peri, source%peri, lsource_mask=lsource_mask)
      call util_append(self%atp, source%atp, lsource_mask=lsource_mask)
      call util_append(self%a, source%a, lsource_mask=lsource_mask)
      call util_append(self%e, source%e, lsource_mask=lsource_mask)
      call util_append(self%inc, source%inc, lsource_mask=lsource_mask)
      call util_append(self%capom, source%capom, lsource_mask=lsource_mask)
      call util_append(self%omega, source%omega, lsource_mask=lsource_mask)
      call util_append(self%capm, source%capm, lsource_mask=lsource_mask)

      self%nbody = self%nbody + count(lsource_mask(:))

      return
   end subroutine swiftest_util_append_body


   module subroutine swiftest_util_append_pl(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self 
         !! Swiftest massive body object
      class(swiftest_body),  intent(in) :: source 
         !! Source object to append
      logical, dimension(:), intent(in) :: lsource_mask 
         !! Logical mask indicating which elements to append to

      select type(source)
      class is (swiftest_pl)
         call util_append(self%mass, source%mass, lsource_mask=lsource_mask)
         call util_append(self%Gmass, source%Gmass, lsource_mask=lsource_mask)
         call util_append(self%rhill, source%rhill, lsource_mask=lsource_mask)
         call util_append(self%renc, source%renc, lsource_mask=lsource_mask)
         call util_append(self%radius, source%radius, lsource_mask=lsource_mask)
         call util_append(self%density, source%density, lsource_mask=lsource_mask)
         call util_append(self%rbeg, source%rbeg, lsource_mask=lsource_mask)
         call util_append(self%rend, source%rend, lsource_mask=lsource_mask)
         call util_append(self%vbeg, source%vbeg, lsource_mask=lsource_mask)
         call util_append(self%Ip, source%Ip, lsource_mask=lsource_mask)
         call util_append(self%rot, source%rot, lsource_mask=lsource_mask)
         call util_append(self%k2, source%k2, lsource_mask=lsource_mask)
         call util_append(self%Q, source%Q, lsource_mask=lsource_mask)
         call util_append(self%tlag, source%tlag, lsource_mask=lsource_mask)
         call util_append(self%kin, source%kin, lsource_mask=lsource_mask)
         call util_append(self%lmtiny, source%lmtiny, lsource_mask=lsource_mask)
         call util_append(self%nplenc, source%nplenc, lsource_mask=lsource_mask)
         call util_append(self%ntpenc, source%ntpenc, lsource_mask=lsource_mask)

         if (allocated(self%k_plpl)) deallocate(self%k_plpl)

         call swiftest_util_append_body(self, source, lsource_mask)
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class swiftest_pl or its descendents"
         call base_util_exit(FAILURE)
      end select

      return
   end subroutine swiftest_util_append_pl


   module subroutine swiftest_util_append_tp(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_tp), intent(inout) :: self 
         !! Swiftest test particle object
      class(swiftest_body), intent(in) :: source  
         !! Source object to append
      logical, dimension(:), intent(in):: lsource_mask 
         !! Logical mask indicating which elements to append to

      select type(source)
      class is (swiftest_tp)
         call util_append(self%nplenc, source%nplenc, lsource_mask=lsource_mask)

         call swiftest_util_append_body(self, source, lsource_mask)
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class swiftest_tp or its descendents"
         call base_util_exit(FAILURE)
      end select

      return
   end subroutine swiftest_util_append_tp


   module subroutine swiftest_util_coord_h2b_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from heliocentric to barycentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b.f90 
      !! Adapted from Hal Levison's Swift routine coord_h2b.f 
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)  :: i
      real(DP)      :: Gmtot
      real(DP), dimension(NDIM) :: xtmp, vtmp

      if (self%nbody == 0) return
      associate(pl => self, npl => self%nbody)
         Gmtot = cb%Gmass
         xtmp(:) = 0.0_DP
         vtmp(:) = 0.0_DP
         do i = 1, npl
            if (pl%status(i) == INACTIVE) cycle
            Gmtot = Gmtot + pl%Gmass(i)
            xtmp(:) = xtmp(:) + pl%Gmass(i) * pl%rh(:,i)
            vtmp(:) = vtmp(:) + pl%Gmass(i) * pl%vh(:,i)
         end do
         cb%rb(:) = -xtmp(:) / Gmtot
         cb%vb(:) = -vtmp(:) / Gmtot
         do i = 1, npl
            if (pl%status(i) == INACTIVE) cycle
            pl%rb(:,i) = pl%rh(:,i) + cb%rb(:)
            pl%vb(:,i) = pl%vh(:,i) + cb%vb(:)
         end do
      end associate

      return
   end subroutine swiftest_util_coord_h2b_pl


   module subroutine swiftest_util_coord_h2b_tp(self, cb)
      !! author: David A. Minton
      !!
      !! Convert test particles from heliocentric to barycentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b_tp.f90 
      !! Adapted from Hal Levison's Swift routine coord_h2b_tp.f 
      implicit none
      ! Arguments
      class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
      class(swiftest_cb), intent(in) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B) :: i, ntp

      if (self%nbody == 0) return
      associate(tp => self)
         ntp = self%nbody
#ifdef DOCONLOC
         do concurrent (i = 1:ntp, tp%status(i) /= INACTIVE) shared(cb,tp)
#else
         do concurrent (i = 1:ntp, tp%status(i) /= INACTIVE)
#endif
            tp%rb(:, i) = tp%rh(:, i) + cb%rb(:)
            tp%vb(:, i) = tp%vh(:, i) + cb%vb(:)
         end do
      end associate

      return
   end subroutine swiftest_util_coord_h2b_tp


   module subroutine swiftest_util_coord_b2h_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from barycentric to heliocentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_b2h.f90 
      !! Adapted from Hal Levison's Swift routine coord_b2h.f 
      implicit none
      ! Arguments
      class(swiftest_pl),     intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb),  intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)          :: i, npl

      if (self%nbody == 0) return

      associate(pl => self)
         npl = self%nbody
#ifdef DOCONLOC
         do concurrent (i = 1:npl, pl%status(i) /= INACTIVE) shared(cb,pl)
#else
         do concurrent (i = 1:npl, pl%status(i) /= INACTIVE)
#endif
            pl%rh(:, i) = pl%rb(:, i) - cb%rb(:)
            pl%vh(:, i) = pl%vb(:, i) - cb%vb(:)
         end do
      end associate

      return
   end subroutine swiftest_util_coord_b2h_pl


   module subroutine swiftest_util_coord_b2h_tp(self, cb)
      !! author: David A. Minton
      !!
      !! Convert test particles from barycentric to heliocentric coordinates (position and velocity)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_b2h_tp.f90 
      !! Adapted from Hal Levison's Swift routine coord_b2h_tp.f 
      implicit none
      ! Arguments
      class(swiftest_tp),     intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb),  intent(in)    :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B) :: i, ntp

      if (self%nbody == 0) return

      associate(tp => self)
         ntp = self%nbody
#ifdef DOCONLOC
         do concurrent(i = 1:ntp, tp%status(i) /= INACTIVE) shared(cb,tp)
#else
         do concurrent(i = 1:ntp, tp%status(i) /= INACTIVE)
#endif
            tp%rh(:, i) = tp%rb(:, i) - cb%rb(:)
            tp%vh(:, i) = tp%vb(:, i) - cb%vb(:)
         end do
      end associate

      return
   end subroutine swiftest_util_coord_b2h_tp


   module subroutine swiftest_util_coord_vb2vh_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from barycentric to heliocentric coordinates (velocity only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vb2vh.f90 
      !! Adapted from Hal Levison's Swift routine coord_vb2vh.f 
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)              :: i, npl

      if (self%nbody == 0) return

      associate(pl => self)
         npl = self%nbody
         cb%vb(:) = 0.0_DP
         do i = npl, 1, -1
            if (pl%status(i) /= INACTIVE) cb%vb(:) = cb%vb(:) - pl%Gmass(i) * pl%vb(:, i) / cb%Gmass
         end do
#ifdef DOCONLOC
         do concurrent(i = 1:npl) shared(cb,pl)
#else
         do concurrent(i = 1:npl)
#endif
            pl%vh(:, i) = pl%vb(:, i) - cb%vb(:)
         end do
      end associate

      return
   end subroutine swiftest_util_coord_vb2vh_pl


   module subroutine swiftest_util_coord_vb2vh_tp(self, vbcb)
      !! author: David A. Minton
      !!
      !! Convert test particles from barycentric to heliocentric coordinates (velocity only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vb2vh_tp.f90 
      !! Adapted from Hal Levison's Swift routine coord_vb2h_tp.f 
      implicit none
      ! Arguments
      class(swiftest_tp),     intent(inout) :: self !! Swiftest test particle object
      real(DP), dimension(:), intent(in)    :: vbcb !! Barycentric velocity of the central body

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         where (tp%lmask(1:ntp))
            tp%vh(1, 1:ntp) = tp%vb(1, 1:ntp) - vbcb(1)
            tp%vh(2, 1:ntp) = tp%vb(2, 1:ntp) - vbcb(2)
            tp%vh(3, 1:ntp) = tp%vb(3, 1:ntp) - vbcb(3)
         end where
      end associate

      return
   end subroutine swiftest_util_coord_vb2vh_tp


   module subroutine swiftest_util_coord_vh2vb_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert massive bodies from heliocentric to barycentric coordinates (velocity only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vh2vb.f90 
      !! Adapted from Hal Levison's Swift routine coord_vh2b.f 
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)  :: i, npl
      real(DP)      :: Gmtot

      if (self%nbody == 0) return

      associate(pl => self)
         npl = self%nbody
         Gmtot = cb%Gmass + sum(pl%Gmass(1:npl))
         cb%vb(:) = 0.0_DP
         do i = 1, npl
            cb%vb(:) = cb%vb(:) - pl%Gmass(i) * pl%vh(:, i) 
         end do
         cb%vb(:) = cb%vb(:) / Gmtot
#ifdef DOCONLOC
         do concurrent(i = 1:npl) shared(cb,pl)
#else
         do concurrent(i = 1:npl)
#endif
            pl%vb(:, i) = pl%vh(:, i) + cb%vb(:)
         end do
      end associate

      return
   end subroutine swiftest_util_coord_vh2vb_pl


   module subroutine swiftest_util_coord_vh2vb_tp(self, vbcb)
      !! author: David A. Minton
      !!
      !! Convert test particles from heliocentric to barycentric coordinates (velocity only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_vh2vb_tp.f90
      !! Adapted from Hal Levison's Swift routine coord_vh2b_tp.f 
      implicit none
      ! Arguments
      class(swiftest_tp),     intent(inout) :: self !! Swiftest test particle object
      real(DP), dimension(:), intent(in)    :: vbcb !! Barycentric velocity of the central body

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         where (tp%lmask(1:ntp))
            tp%vb(1, 1:ntp) = tp%vh(1, 1:ntp) + vbcb(1)
            tp%vb(2, 1:ntp) = tp%vh(2, 1:ntp) + vbcb(2)
            tp%vb(3, 1:ntp) = tp%vh(3, 1:ntp) + vbcb(3)
         end where
      end associate

      return
   end subroutine swiftest_util_coord_vh2vb_tp
   

   module subroutine swiftest_util_coord_rh2rb_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Convert position vectors of massive bodies from heliocentric to barycentric coordinates (position only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b.f90 
      !! Adapted from Hal Levison's Swift routine coord_h2b.f 
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)  :: i
      real(DP)      :: Gmtot
      real(DP), dimension(NDIM) :: xtmp

      if (self%nbody == 0) return
      associate(pl => self, npl => self%nbody)
         Gmtot = cb%Gmass
         xtmp(:) = 0.0_DP
         do i = 1, npl
            if (pl%status(i) == INACTIVE) cycle
            Gmtot = Gmtot + pl%Gmass(i)
            xtmp(:) = xtmp(:) + pl%Gmass(i) * pl%rh(:,i)
         end do
         cb%rb(:) = -xtmp(:) / Gmtot
         do i = 1, npl
            if (pl%status(i) == INACTIVE) cycle
            pl%rb(:,i) = pl%rh(:,i) + cb%rb(:)
         end do
      end associate

      return
   end subroutine swiftest_util_coord_rh2rb_pl


   module subroutine swiftest_util_coord_rh2rb_tp(self, cb)
      !! author: David A. Minton
      !!
      !! Convert test particles from heliocentric to barycentric coordinates (position only)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine coord_h2b_tp.f90 
      !! Adapted from Hal Levison's Swift routine coord_h2b_tp.f 
      implicit none
      ! Arguments
      class(swiftest_tp), intent(inout) :: self !! Swiftest test particle object
      class(swiftest_cb), intent(in) :: cb      !! Swiftest central body object
      ! Internals
      integer(I4B) :: i, ntp

      if (self%nbody == 0) return
      associate(tp => self)
         ntp = self%nbody
#ifdef DOCONLOC
         do concurrent (i = 1:ntp, tp%status(i) /= INACTIVE) shared(cb,tp)
#else
         do concurrent (i = 1:ntp, tp%status(i) /= INACTIVE)
#endif
            tp%rb(:, i) = tp%rh(:, i) + cb%rb(:)
         end do
      end associate

      return
   end subroutine swiftest_util_coord_rh2rb_tp


   module subroutine swiftest_util_copy_particle_info(self, source)
      !! author: David A. Minton
      !!
      !! Copies one set of information object components into another, component-by-component
      implicit none
      class(swiftest_particle_info),  intent(inout) :: self
      class(swiftest_particle_info),  intent(in)    :: source

      call self%set_value(&
         name = source%name, &
         particle_type = source%particle_type, &
         status = source%status, & 
         origin_type = source%origin_type, &
         origin_time = source%origin_time, & 
         collision_id = source%collision_id, &
         origin_rh = source%origin_rh(:), &
         origin_vh = source%origin_vh(:), &
         discard_time = source%discard_time, & 
         discard_rh = source%discard_rh(:), &
         discard_vh = source%discard_vh(:), &
         discard_body_id = source%discard_body_id &
      )

      return
   end subroutine swiftest_util_copy_particle_info


   module subroutine swiftest_util_copy_particle_info_arr(source, dest, idx)
      !! author: David A. Minton
      !!
      !! Copies contents from an array of one particle information objects to another.
      implicit none
      class(swiftest_particle_info), dimension(:), intent(in)             :: source !! Source object to copy into
      class(swiftest_particle_info), dimension(:), intent(inout)          :: dest   !! Swiftest body object with particle metadata 
                                                                                    !! information object
      integer(I4B),                  dimension(:), intent(in),   optional :: idx    !! Optional array of indices to draw the source 
                                                                                    !! object
      ! Internals
      integer(I4B) :: i, j, n, nsource, ndest

      if (size(source) == 0) return

      if (present(idx)) then
         n = size(idx)
      else
         n = size(source)
      end if

      nsource = size(source)
      ndest = size(dest)

      if ((n == 0) .or. (n > ndest) .or. (n > nsource)) then
         write(*,*) 'Particle info copy operation failed. n, nsource, ndest: ',n, nsource, ndest
         return
      end if

      do i = 1, n
         if (present(idx)) then
            j = idx(i)
         else
            j = i
         end if 
         call dest(i)%copy(source(j))
      end do
      
      return
   end subroutine swiftest_util_copy_particle_info_arr


   module subroutine swiftest_util_dealloc_body(self)
      !! author: David A. Minton
      !!
      !! Finalize the Swiftest body object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_body),  intent(inout) :: self

      self%lfirst = .true.
      self%nbody = 0
      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%info)) deallocate(self%info)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%lmask)) deallocate(self%lmask)
      if (allocated(self%ldiscard)) deallocate(self%ldiscard)
      if (allocated(self%lcollision)) deallocate(self%lcollision)
      if (allocated(self%lencounter)) deallocate(self%lencounter)
      if (allocated(self%mu)) deallocate(self%mu)
      if (allocated(self%rh)) deallocate(self%rh)
      if (allocated(self%vh)) deallocate(self%vh)
      if (allocated(self%rb)) deallocate(self%rb)
      if (allocated(self%vb)) deallocate(self%vb)
      if (allocated(self%ah)) deallocate(self%ah)
      if (allocated(self%aobl)) deallocate(self%aobl)
      if (allocated(self%agr)) deallocate(self%agr)
      if (allocated(self%atide)) deallocate(self%atide)
      if (allocated(self%ir3h)) deallocate(self%ir3h)
      if (allocated(self%isperi)) deallocate(self%isperi)
      if (allocated(self%peri)) deallocate(self%peri)
      if (allocated(self%atp)) deallocate(self%atp)
      if (allocated(self%a)) deallocate(self%a)
      if (allocated(self%e)) deallocate(self%e)
      if (allocated(self%inc)) deallocate(self%inc)
      if (allocated(self%capom)) deallocate(self%capom)
      if (allocated(self%omega)) deallocate(self%omega)
      if (allocated(self%capm)) deallocate(self%capm)

      return
   end subroutine swiftest_util_dealloc_body


   module subroutine swiftest_util_dealloc_cb(self)
      !! author: David A. Minton
      !!
      !! Finalize the Swiftest central body object - deallocates all allocatables
      implicit none
      ! Arguments
      class(swiftest_cb), intent(inout) :: self !! Swiftest central body object

      if (allocated(self%info)) deallocate(self%info)

      self%id       = 0      
      self%mass     = 0.0_DP 
      self%Gmass    = 0.0_DP 
      self%radius   = 0.0_DP 
      self%density  = 1.0_DP 
      self%j2rp2    = 0.0_DP 
      self%j4rp4    = 0.0_DP 
      self%aobl     = 0.0_DP 
      self%atide    = 0.0_DP 
      self%aoblbeg  = 0.0_DP 
      self%aoblend  = 0.0_DP 
      self%atidebeg = 0.0_DP 
      self%atideend = 0.0_DP 
      self%rb       = 0.0_DP 
      self%vb       = 0.0_DP 
      self%agr      = 0.0_DP 
      self%Ip       = 0.0_DP 
      self%rot      = 0.0_DP 
      self%k2       = 0.0_DP 
      self%Q        = 0.0_DP 
      self%tlag     = 0.0_DP 
      self%L0       = 0.0_DP 
      self%dL       = 0.0_DP 
      self%GM0      = 0.0_DP 
      self%dGM      = 0.0_DP 
      self%R0       = 0.0_DP 
      self%dR       = 0.0_DP 

      return
   end subroutine swiftest_util_dealloc_cb


   module subroutine swiftest_util_dealloc_kin(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Arguments
      class(swiftest_kinship),  intent(inout) :: self !! Swiftest kinship object

      if (allocated(self%child)) deallocate(self%child)

      return
   end subroutine swiftest_util_dealloc_kin


   module subroutine swiftest_util_dealloc_pl(self)
      !! author: David A. Minton
      !!
      !! Finalize the Swiftest massive body object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_pl),  intent(inout) :: self !! Swiftest massive body object
      ! Internals
      integer(I4B) :: i

      if (allocated(self%mass)) deallocate(self%mass)
      if (allocated(self%Gmass)) deallocate(self%Gmass)
      if (allocated(self%rhill)) deallocate(self%rhill)
      if (allocated(self%renc)) deallocate(self%renc)
      if (allocated(self%radius)) deallocate(self%radius)
      if (allocated(self%density)) deallocate(self%density)
      if (allocated(self%rbeg)) deallocate(self%rbeg)
      if (allocated(self%rend)) deallocate(self%rend)
      if (allocated(self%vbeg)) deallocate(self%vbeg)
      if (allocated(self%Ip)) deallocate(self%Ip)
      if (allocated(self%rot)) deallocate(self%rot)
      if (allocated(self%k2)) deallocate(self%k2)
      if (allocated(self%Q)) deallocate(self%Q)
      if (allocated(self%tlag)) deallocate(self%tlag)
      if (allocated(self%k_plpl)) deallocate(self%k_plpl)
      if (allocated(self%lmtiny)) deallocate(self%lmtiny)
      if (allocated(self%nplenc)) deallocate(self%nplenc)
      if (allocated(self%ntpenc)) deallocate(self%ntpenc)


      if (allocated(self%kin)) then
         do i = 1, self%nbody
            call self%kin(i)%dealloc()
         end do
         deallocate(self%kin)
      end if

      call swiftest_util_dealloc_body(self)

      return
   end subroutine swiftest_util_dealloc_pl


   module subroutine swiftest_util_dealloc_storage(self)
      !! author: David A. Minton
      !!
      !! Resets a storage object by deallocating all items and resetting the frame counter to 0
      use base, only : base_util_dealloc_storage
      implicit none
      ! Arguments
      class(swiftest_storage), intent(inout) :: self !! Swiftest storage object

      if (allocated(self%nc)) deallocate(self%nc)
      call base_util_dealloc_storage(self)

      return
   end subroutine swiftest_util_dealloc_storage


   module subroutine swiftest_util_dealloc_system(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatables and resets all values to defaults. Acts as a base for a finalizer
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self

      if (allocated(self%cb)) deallocate(self%cb)
      if (allocated(self%pl)) deallocate(self%pl)
      if (allocated(self%tp)) deallocate(self%tp)
      if (allocated(self%tp_discards)) deallocate(self%tp_discards)
      if (allocated(self%pl_discards)) deallocate(self%pl_discards)
      if (allocated(self%pl_adds)) deallocate(self%pl_adds)
      if (allocated(self%tp_adds)) deallocate(self%tp_adds)
      if (allocated(self%pltp_encounter)) deallocate(self%pltp_encounter)
      if (allocated(self%plpl_encounter)) deallocate(self%plpl_encounter)
      if (allocated(self%plpl_collision)) deallocate(self%plpl_collision)
      if (allocated(self%pltp_collision)) deallocate(self%pltp_collision)
      if (allocated(self%collider)) deallocate(self%collider)
      if (allocated(self%encounter_history)) deallocate(self%encounter_history)
      if (allocated(self%collision_history)) deallocate(self%collision_history)

      self%t = -1.0_DP            
      self%GMtot = 0.0_DP         
      self%ke_orbit = 0.0_DP      
      self%ke_spin = 0.0_DP       
      self%pe = 0.0_DP            
      self%be = 0.0_DP            
      self%te = 0.0_DP            
      self%oblpot = 0.0_DP        
      self%L_orbit = 0.0_DP        
      self%L_spin = 0.0_DP         
      self%L_total = 0.0_DP          
      self%ke_orbit_orig = 0.0_DP 
      self%ke_spin_orig = 0.0_DP  
      self%pe_orig = 0.0_DP       
      self%be_orig = 0.0_DP       
      self%E_orbit_orig = 0.0_DP   
      self%GMtot_orig = 0.0_DP    
      self%L_total_orig = 0.0_DP     
      self%L_orbit_orig = 0.0_DP   
      self%L_spin_orig = 0.0_DP    
      self%L_escape = 0.0_DP       
      self%GMescape = 0.0_DP      
      self%E_collisions = 0.0_DP   
      self%E_untracked = 0.0_DP    

      self%ke_orbit_error    = 0.0_DP
      self%ke_spin_error     = 0.0_DP
      self%pe_error          = 0.0_DP
      self%be_error          = 0.0_DP
      self%E_orbit_error     = 0.0_DP
      self%Ecoll_error       = 0.0_DP
      self%E_untracked_error = 0.0_DP
      self%te_error          = 0.0_DP
      self%L_orbit_error     = 0.0_DP
      self%L_spin_error      = 0.0_DP
      self%L_escape_error    = 0.0_DP
      self%L_total_error     = 0.0_DP
      self%Mtot_error        = 0.0_DP
      self%Mescape_error     = 0.0_DP

      return
   end subroutine swiftest_util_dealloc_system


   module subroutine swiftest_util_dealloc_tp(self)
      !! author: David A. Minton
      !!
      !! Finalize the Swiftest test particle object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_tp),  intent(inout) :: self !! Swiftest test particle object

      if (allocated(self%k_pltp)) deallocate(self%k_pltp)
      if (allocated(self%nplenc)) deallocate(self%nplenc)

      call swiftest_util_dealloc_body(self)

      return
   end subroutine swiftest_util_dealloc_tp


   module subroutine swiftest_util_fill_arr_info(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of particle origin information types
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      type(swiftest_particle_info), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,                      dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into
                                                                                           !! the keeps
      ! Internals
      integer(I4B), dimension(:), allocatable  :: insert_idx
      integer(I4B) :: i, nkeep, ninsert

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      nkeep = size(keeps)
      ninsert = count(lfill_list)

      allocate(insert_idx(ninsert))

      insert_idx(:) = pack([(i, i = 1, nkeep)], lfill_list)
      call swiftest_util_copy_particle_info_arr(inserts, keeps, insert_idx)

      return
   end subroutine swiftest_util_fill_arr_info


   module subroutine swiftest_util_fill_arr_kin(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of particle kinship types
      !! This is the inverse of a spill operation   
      implicit none
      ! Arguments
      type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      type(swiftest_kinship), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,             dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))
   
      return
   end subroutine swiftest_util_fill_arr_kin


   module subroutine swiftest_util_fill_body(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest generic particle structure into an old one. 
      !! This is the inverse of a spill operation.
      implicit none
      ! Arguments
      class(swiftest_body),  intent(inout) :: self       !! Swiftest generic body object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !! Fill all the common components
      associate(keeps => self)
         call util_fill(keeps%id,         inserts%id,         lfill_list)
         call util_fill(keeps%info,       inserts%info,       lfill_list)
         call util_fill(keeps%lmask,      inserts%lmask,      lfill_list)
         call util_fill(keeps%status,     inserts%status,     lfill_list)
         call util_fill(keeps%ldiscard,   inserts%ldiscard,   lfill_list)
         call util_fill(keeps%lcollision, inserts%lcollision, lfill_list)
         call util_fill(keeps%lencounter, inserts%lencounter, lfill_list)
         call util_fill(keeps%mu,         inserts%mu,         lfill_list)
         call util_fill(keeps%rh,         inserts%rh,         lfill_list)
         call util_fill(keeps%vh,         inserts%vh,         lfill_list)
         call util_fill(keeps%rb,         inserts%rb,         lfill_list)
         call util_fill(keeps%vb,         inserts%vb,         lfill_list)
         call util_fill(keeps%ah,         inserts%ah,         lfill_list)
         call util_fill(keeps%aobl,       inserts%aobl,       lfill_list)
         call util_fill(keeps%agr,        inserts%agr,        lfill_list)
         call util_fill(keeps%atide,      inserts%atide,      lfill_list)
         call util_fill(keeps%ir3h,       inserts%ir3h,       lfill_list)
         call util_fill(keeps%isperi,     inserts%isperi,     lfill_list)
         call util_fill(keeps%peri,       inserts%peri,       lfill_list)
         call util_fill(keeps%atp,        inserts%atp,        lfill_list)
         call util_fill(keeps%a,          inserts%a,          lfill_list)
         call util_fill(keeps%e,          inserts%e,          lfill_list)
         call util_fill(keeps%inc,        inserts%inc,        lfill_list)
         call util_fill(keeps%capom,      inserts%capom,      lfill_list)
         call util_fill(keeps%omega,      inserts%omega,      lfill_list)
         call util_fill(keeps%capm,       inserts%capm,       lfill_list)
           
         ! This is the base class, so will be the last to be called in the cascade. 
         keeps%nbody = size(keeps%id(:))
      end associate
     
      return
   end subroutine swiftest_util_fill_body


   module subroutine swiftest_util_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest massive body structure into an old one. 
      !! This is the inverse of a spill operation.
      implicit none
      ! Arguments
      class(swiftest_pl),    intent(inout) :: self       !! Swiftest massive body object
      class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)

      select type (inserts) ! The standard requires us to select the type of both arguments in order to access all the components
         class is (swiftest_pl)
            !! Fill components specific to the massive body class
            call util_fill(keeps%mass,    inserts%mass,    lfill_list)
            call util_fill(keeps%Gmass,   inserts%Gmass,   lfill_list)
            call util_fill(keeps%rhill,   inserts%rhill,   lfill_list)
            call util_fill(keeps%renc,    inserts%renc,    lfill_list)
            call util_fill(keeps%radius,  inserts%radius,  lfill_list)
            call util_fill(keeps%density, inserts%density, lfill_list)
            call util_fill(keeps%rbeg,    inserts%rbeg,    lfill_list)
            call util_fill(keeps%rend,    inserts%rend,    lfill_list)
            call util_fill(keeps%vbeg,    inserts%vbeg,    lfill_list)
            call util_fill(keeps%Ip,      inserts%Ip,      lfill_list)
            call util_fill(keeps%rot,     inserts%rot,     lfill_list)
            call util_fill(keeps%k2,      inserts%k2,      lfill_list)
            call util_fill(keeps%Q,       inserts%Q,       lfill_list)
            call util_fill(keeps%tlag,    inserts%tlag,    lfill_list)
            call util_fill(keeps%kin,     inserts%kin,     lfill_list)
            call util_fill(keeps%nplenc,  inserts%nplenc,  lfill_list)
            call util_fill(keeps%ntpenc,  inserts%ntpenc,  lfill_list)

            if (allocated(keeps%k_plpl)) deallocate(keeps%k_plpl)
            
            call swiftest_util_fill_body(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on swiftest_pl'
         end select
      end associate

      return
   end subroutine swiftest_util_fill_pl


   module subroutine swiftest_util_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      implicit none
      ! Arguments
      class(swiftest_tp),    intent(inout) :: self       !! Swiftest test particle object
      class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (swiftest_tp)
            !! Spill components specific to the test particle class
            call util_fill(keeps%nplenc,  inserts%nplenc,  lfill_list)

            call swiftest_util_fill_body(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on swiftest_tp'
         end select
      end associate

      return
   end subroutine swiftest_util_fill_tp


   pure module subroutine swiftest_util_flatten_eucl_ij_to_k(n, i, j, k)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix for pl-pl interactions.
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      integer(I4B), intent(in)  :: n !! Number of bodies
      integer(I4B), intent(in)  :: i !! Index of the ith body
      integer(I4B), intent(in)  :: j !! Index of the jth body
      integer(I8B), intent(out) :: k !! Index of the flattened matrix
      ! Internals
      integer(I8B) :: i8, j8, n8
     
      i8 = int(i, kind=I8B)
      j8 = int(j, kind=I8B)
      n8 = int(n, kind=I8B)
      k = (i8 - 1_I8B) * n8 - i8 * (i8 - 1_I8B) / 2_I8B + (j8 - i8)

      return
   end subroutine swiftest_util_flatten_eucl_ij_to_k


   pure module subroutine swiftest_util_flatten_eucl_k_to_ij(n, k, i, j)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns k index into i,j indices for use in the Euclidean distance matrix for pl-pl interactions.
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      integer(I4B), intent(in)  :: n !! Number of bodies
      integer(I8B), intent(in)  :: k !! Index of the flattened matrix
      integer(I4B), intent(out) :: i !! Index of the ith body
      integer(I4B), intent(out) :: j !! Index of the jth body
      ! Internals
      integer(I8B) :: kp, p, i8, j8, n8

      n8 = int(n, kind=I8B)
    
      kp = n8 * (n8 - 1_I8B) / 2_I8B - k
      p = floor((sqrt(1._DP + 8_I8B * kp) - 1_I8B) / 2_I8B)
      i8 = n8 - 1_I8B - p
      j8 = k - (n8 - 1_I8B) * (n8 - 2_I8B) / 2_I8B + p * (p + 1_I8B) / 2_I8B + 1_I8B

      i = int(i8, kind=I4B)
      j = int(j8, kind=I4B)

      return
   end subroutine swiftest_util_flatten_eucl_k_to_ij


   module subroutine swiftest_util_flatten_eucl_plpl(self, param)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix for pl-pl interactions for a Swiftest massive body 
      !! object
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: err, i, j, npl
      integer(I8B) :: k, npl8

      associate(nplpl => self%nplpl)
         npl = self%nbody
         npl8 = int(npl, kind=I8B)
         nplpl = npl8 * (npl8 - 1_I8B) / 2_I8B ! number of entries in a strict lower triangle, npl x npl
         if (param%lflatten_interactions) then
            if (allocated(self%k_plpl)) deallocate(self%k_plpl) ! Reset the index array if it's been set previously
            allocate(self%k_plpl(2, nplpl), stat=err)
            if (err /=0) then ! An error occurred trying to allocate this big array. This probably means it's too big to fit in 
                              ! memory, and so we will force the run back into triangular mode
               param%lflatten_interactions = .false.
            else
#ifdef DOCONLOC
               do concurrent (i=1:npl, j=1:npl, j>i) shared(self) local(k)
#else
               do concurrent (i=1:npl, j=1:npl, j>i)
#endif
                  call swiftest_util_flatten_eucl_ij_to_k(npl, i, j, k)
                  self%k_plpl(1, k) = i
                  self%k_plpl(2, k) = j
               end do
            end if
         end if
      end associate

      return
   end subroutine swiftest_util_flatten_eucl_plpl


   module subroutine swiftest_util_flatten_eucl_pltp(self, pl, param)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix for pl-tp interactions
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
      class(swiftest_pl),         intent(in)    :: pl    !! Swiftest massive body object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i, j
      integer(I8B) :: counter, npl8, ntp8

      associate(ntp => self%nbody, npl => pl%nbody, npltp => self%npltp)
         npl8 = int(npl, kind=I8B)
         ntp8 = int(ntp, kind=I8B)
         npltp = npl8 * ntp8 
         if (allocated(self%k_pltp)) deallocate(self%k_pltp) ! Reset the index array if it's been set previously
         allocate(self%k_pltp(2, npltp))
         counter = 1_I8B
         do i = 1, npl
            do j = 1,  ntp
               self%k_pltp(1, counter) = i
               self%k_pltp(2, counter) = j
               counter = counter + 1_I8B
            end do
         end do
      end associate

      return
   end subroutine swiftest_util_flatten_eucl_pltp


   module subroutine swiftest_util_get_energy_and_momentum_system(self, param)
      !! author: David A. Minton
      !!
      !! Compute total nbody_system angular momentum vector and kinetic, potential and total nbody_system energy
      !!  
      !! Adapted from David E. Kaufmann Swifter routine symba_energy_eucl.f90
      !!  
      !! Adapted from Martin Duncan's Swift routine anal_energy.f
      implicit none
      class(swiftest_nbody_system), intent(inout) :: self     !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param    !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i,j, npl
      real(DP) :: kecb, kespincb
      real(DP), dimension(self%pl%nbody) :: kepl, kespinpl
      real(DP), dimension(NDIM,self%pl%nbody) :: Lplorbit
      real(DP), dimension(NDIM,self%pl%nbody) :: Lplspin
      real(DP), dimension(NDIM) :: Lcborbit, Lcbspin
      real(DP), dimension(NDIM) :: h

      associate(nbody_system => self, pl => self%pl, cb => self%cb)
         npl = self%pl%nbody
         nbody_system%L_orbit(:) = 0.0_DP
         nbody_system%L_spin(:) = 0.0_DP
         nbody_system%L_total(:) = 0.0_DP
         nbody_system%ke_orbit = 0.0_DP
         nbody_system%ke_spin = 0.0_DP

         nbody_system%GMtot = cb%Gmass
         if (npl > 0) then
            kepl(:) = 0.0_DP
            Lplorbit(:,:) = 0.0_DP
            Lplspin(:,:) = 0.0_DP
            pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE
            nbody_system%GMtot = nbody_system%GMtot + sum(pl%Gmass(1:npl), pl%lmask(1:npl)) 
         end if
            
         kecb = cb%mass * dot_product(cb%vb(:), cb%vb(:))
         nbody_system%be_cb = -3*cb%Gmass * cb%mass / (5 * cb%radius) 
         Lcborbit(:) = cb%mass * (cb%rb(:) .cross. cb%vb(:))
         if (npl > 0) then

#ifdef DOCONLOC
            do concurrent (i = 1:npl, pl%lmask(i)) shared(pl,Lplorbit,kepl,npl) local(h) 
#else
            do concurrent (i = 1:npl, pl%lmask(i))
#endif
               h(1) = pl%rb(2,i) * pl%vb(3,i) - pl%rb(3,i) * pl%vb(2,i)
               h(2) = pl%rb(3,i) * pl%vb(1,i) - pl%rb(1,i) * pl%vb(3,i)
               h(3) = pl%rb(1,i) * pl%vb(2,i) - pl%rb(2,i) * pl%vb(1,i)
   
            ! Angular momentum from orbit 
               Lplorbit(:,i) = pl%mass(i) * h(:)

               ! Kinetic energy from orbit
               kepl(i) = pl%mass(i) * dot_product(pl%vb(:,i), pl%vb(:,i)) 
            end do
         end if

         if (param%lrotation) then
            kespincb = cb%mass * cb%Ip(3) * cb%radius**2 * dot_product(cb%rot(:), cb%rot(:))

            ! For simplicity, we always assume that the rotation pole is the 3rd principal axis
            Lcbspin(:) = cb%Ip(3) * cb%mass * cb%radius**2 * cb%rot(:)

            if (npl > 0) then
#ifdef DOCONLOC
               do concurrent (i = 1:npl, pl%lmask(i)) shared(pl,Lplspin,kespinpl)
#else
               do concurrent (i = 1:npl, pl%lmask(i))
#endif
                  ! Currently we assume that the rotation pole is the 3rd principal axis
                  ! Angular momentum from spin
                  Lplspin(:,i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(:,i)

                  ! Kinetic energy from spin
                  kespinpl(i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * dot_product(pl%rot(:,i), pl%rot(:,i))
               end do

               nbody_system%ke_spin = 0.5_DP * (kespincb + sum(kespinpl(1:npl), pl%lmask(1:npl)))
            else
               nbody_system%ke_spin = 0.5_DP * kespincb
            end if

            if (npl > 0) then
#ifdef DOCONLOC
               do concurrent (j = 1:NDIM) shared(nbody_system,pl,Lplspin,Lcbspin)
#else
               do concurrent (j = 1:NDIM)
#endif
                  nbody_system%L_spin(j) = Lcbspin(j) + sum(Lplspin(j,1:npl), pl%lmask(1:npl))
               end do
            else
               nbody_system%L_spin(:) = Lcbspin(:)
            end if
         else
            nbody_system%ke_spin = 0.0_DP
            nbody_system%L_spin(:) = 0.0_DP
         end if
 
         if (npl > 0) then
            if (param%lflatten_interactions) then
               call swiftest_util_get_potential_energy(npl, pl%nplpl, pl%k_plpl, pl%lmask, cb%Gmass, pl%Gmass, pl%mass, pl%rb, &
                                                      nbody_system%pe)
            else
               call swiftest_util_get_potential_energy(npl, pl%lmask, cb%Gmass, pl%Gmass, pl%mass, pl%rb, nbody_system%pe)
            end if
         end if

         ! Potential energy from the oblateness term
         if (param%lnon_spherical_cb) then
            call nbody_system%obl_pot()
            nbody_system%pe = nbody_system%pe + nbody_system%oblpot
         end if

         if (npl > 0) then
            nbody_system%ke_orbit = 0.5_DP * (kecb + sum(kepl(1:npl), pl%lmask(1:npl)))
#ifdef DOCONLOC
            do concurrent (j = 1:NDIM) shared(nbody_system,pl,Lcborbit,Lplorbit,npl)
#else  
            do concurrent (j = 1:NDIM)
#endif
               nbody_system%L_orbit(j) = Lcborbit(j) + sum(Lplorbit(j,1:npl), pl%lmask(1:npl)) 
            end do
         else
            nbody_system%ke_orbit = 0.5_DP * kecb
            nbody_system%L_orbit(:) = Lcborbit(:)
         end if

         if ((param%lclose .and. (npl > 0))) then
            nbody_system%be = sum(-3*pl%Gmass(1:npl)*pl%mass(1:npl)/(5*pl%radius(1:npl)), pl%lmask(1:npl)) 
         else
            nbody_system%be = 0.0_DP
         end if
         nbody_system%te = nbody_system%ke_orbit + nbody_system%ke_spin + nbody_system%pe + nbody_system%be 
         nbody_system%L_total(:) = nbody_system%L_orbit(:) + nbody_system%L_spin(:)
      end associate

      return
   end subroutine swiftest_util_get_energy_and_momentum_system


   module subroutine swiftest_util_get_potential_energy_flat(npl, nplpl, k_plpl, lmask, GMcb, Gmass, mass, rb, pe)
      !! author: David A. Minton
      !!
      !! Compute total nbody_system potential energy
      implicit none
      ! Arguments
      integer(I4B),                 intent(in)  :: npl
      integer(I8B),                 intent(in)  :: nplpl
      integer(I4B), dimension(:,:), intent(in)  :: k_plpl
      logical,      dimension(:),   intent(in)  :: lmask
      real(DP),                     intent(in)  :: GMcb
      real(DP),     dimension(:),   intent(in)  :: Gmass
      real(DP),     dimension(:),   intent(in)  :: mass
      real(DP),     dimension(:,:), intent(in)  :: rb
      real(DP),                     intent(out) :: pe
      ! Internals
      integer(I4B) :: i, j
      integer(I8B) :: k
      real(DP), dimension(npl) :: pecb
      real(DP), dimension(nplpl) :: pepl 
      logical, dimension(nplpl) :: lstatpl

      ! Do the central body potential energy component first
      where(.not. lmask(1:npl))
         pecb(1:npl) = 0.0_DP
      end where

#ifdef DOCONLOC
      do concurrent(i = 1:npl, lmask(i)) shared(lmask,pecb,GMcb,mass,rb)
#else
      do concurrent(i = 1:npl, lmask(i))
#endif
         pecb(i) = -GMcb * mass(i) / norm2(rb(:,i)) 
      end do

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(k_plpl, rb, mass, Gmass, pepl, lstatpl, lmask) &
      !$omp firstprivate(nplpl)
      do k = 1, nplpl
         i = k_plpl(1,k)
         j = k_plpl(2,k)
         lstatpl(k) = (lmask(i) .and. lmask(j))
         if (lstatpl(k)) then
            pepl(k) = -(Gmass(i) * mass(j)) / norm2(rb(:, i) - rb(:, j))
         else
            pepl(k) = 0.0_DP
         end if
      end do
      !$omp end parallel do 

      pe = sum(pepl(:), lstatpl(:)) + sum(pecb(1:npl), lmask(1:npl))

      return
   end subroutine swiftest_util_get_potential_energy_flat


   module subroutine swiftest_util_get_potential_energy_triangular(npl, lmask, GMcb, Gmass, mass, rb, pe)
      !! author: David A. Minton
      !!
      !! Compute total nbody_system potential energy
      implicit none
      ! Arguments
      integer(I4B),                 intent(in)  :: npl
      logical,      dimension(:),   intent(in)  :: lmask
      real(DP),                     intent(in)  :: GMcb
      real(DP),     dimension(:),   intent(in)  :: Gmass
      real(DP),     dimension(:),   intent(in)  :: mass
      real(DP),     dimension(:,:), intent(in)  :: rb
      real(DP),                     intent(out) :: pe
      ! Internals
      integer(I4B) :: i, j
      real(DP), dimension(npl) :: pecb, pepl

      ! Do the central body potential energy component first
      where(.not. lmask(1:npl))
         pecb(1:npl) = 0.0_DP
      end where

#ifdef DOCONLOC
      do concurrent(i = 1:npl, lmask(i)) shared(lmask, pecb, GMcb, mass, rb, lmask)
#else
      do concurrent(i = 1:npl, lmask(i))
#endif
         pecb(i) = -GMcb * mass(i) / norm2(rb(:,i)) 
      end do

      pe = 0.0_DP
      !$omp parallel do default(private) schedule(static)&
      !$omp shared(lmask, Gmass, mass, rb) &
      !$omp firstprivate(npl) &
      !$omp reduction(+:pe) 
      do i = 1, npl
         if (lmask(i)) then
#ifdef DOCONLOC
            do concurrent(j = i+1:npl, lmask(i) .and. lmask(j)) shared(lmask, pepl, rb, mass, Gmass, lmask) 
#else
            do concurrent(j = i+1:npl, lmask(i) .and. lmask(j))
#endif
               pepl(j) = - (Gmass(i) * mass(j)) / norm2(rb(:, i) - rb(:, j))
            end do
            pe = pe + sum(pepl(i+1:npl), lmask(i+1:npl))
         end if
      end do
      !$omp end parallel do
      pe = pe + sum(pecb(1:npl), lmask(1:npl))

      return
   end subroutine swiftest_util_get_potential_energy_triangular


   module subroutine swiftest_util_get_idvalues_system(self, idvals)
      !! author: David A. Minton
      !!
      !! Returns an array of all id values saved in this snapshot
      implicit none
      ! Arguments
      class(swiftest_nbody_system),            intent(in)  :: self   !! Encounter snapshot object
      integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values saved in this snapshot
      ! Internals
      integer(I4B) :: npl, ntp

      if (allocated(self%pl)) then
         npl = self%pl%nbody
      else
         npl = 0
      end if 
      if (allocated(self%tp)) then
         ntp = self%tp%nbody
      else
         ntp = 0
      end if

      allocate(idvals(1 + npl+ntp))

      idvals(1) = self%cb%id
      if (npl > 0) idvals(2:npl+1) = self%pl%id(:)
      if (ntp > 0) idvals(npl+2:npl+ntp+1) = self%tp%id(:)

      return

   end subroutine swiftest_util_get_idvalues_system


   module subroutine swiftest_util_get_vals_storage(self, idvals, tvals)
      !! author: David A. Minton
      !!
      !! Gets the id values in a storage object, regardless of whether it is encounter of collision
      ! Argument
      class(swiftest_storage),              intent(in)  :: self   !! Swiftest storage object
      integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values in all snapshots
      real(DP),     dimension(:), allocatable, intent(out) :: tvals  !! Array of all time values in all snapshots
      ! Internals
      integer(I4B) :: i, n, nlo, nhi, ntotal
      integer(I4B), dimension(:), allocatable :: itmp

      associate(storage => self, nsnaps => self%iframe)

         allocate(tvals(nsnaps))
         tvals(:) = 0.0_DP

         ! First pass to get total number of ids
         ntotal = 0
         do i = 1, nsnaps
            if (allocated(storage%frame(i)%item)) then
               select type(snapshot => storage%frame(i)%item)
               class is (swiftest_nbody_system)
                  tvals(i) = snapshot%t
                  call snapshot%get_idvals(itmp)
                  if (allocated(itmp)) then
                     n = size(itmp)
                     ntotal = ntotal + n
                  end if
               end select
            end if
         end do

         allocate(idvals(ntotal))
         nlo = 1
         ! Second pass to store all ids get all of the ids stored
         do i = 1, nsnaps
            if (allocated(storage%frame(i)%item)) then
               select type(snapshot => storage%frame(i)%item)
               class is (swiftest_nbody_system)
                  tvals(i) = snapshot%t
                  call snapshot%get_idvals(itmp)
                  if (allocated(itmp)) then
                     n = size(itmp)
                     nhi = nlo + n - 1
                     idvals(nlo:nhi) = itmp(1:n)
                     nlo = nhi + 1 
                  end if
               end select
            end if
         end do

      end associate 
      return
   end subroutine swiftest_util_get_vals_storage


   module subroutine swiftest_util_index_array(ind_arr, n)
      !! author: David A. Minton
      !!
      !! Creates or resizes an index array of size n where ind_arr = [1, 2, ... n].
      !! This subroutine swiftest_assumes that if ind_arr is already allocated, it is a pre-existing index array of a different size
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind_arr !! Index array. Input is a pre-existing index array where 
                                                                        !! n /= size(ind_arr). Output is a new index array 
                                                                        !! ind_arr = [1, 2, ... n]
      integer(I4B),                            intent(in)    :: n       !! The new size of the index array
      ! Internals
      integer(I4B) :: nold, i
      integer(I4B), dimension(:), allocatable :: itmp

      if (allocated(ind_arr)) then
         nold = size(ind_arr)
         if (nold == n) return ! Nothing to do, so go home
      else
         nold = 0
      end if

      allocate(itmp(n))
      if (n >= nold) then
         if (nold > 0) itmp(1:nold) = ind_arr(1:nold)
         itmp(nold+1:n) = [(i, i = nold + 1, n)]
         call move_alloc(itmp, ind_arr)
      else
         itmp(1:n) = ind_arr(1:n)
         call move_alloc(itmp, ind_arr)
      end if

      return
   end subroutine swiftest_util_index_array


   module subroutine swiftest_util_index_map_storage(self)
      !! author: David A. Minton
      !!
      !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      implicit none
      ! Arguments
      class(swiftest_storage), intent(inout) :: self  !! Swiftest storage object
      ! Internals
      integer(I4B), dimension(:), allocatable :: idvals
      real(DP), dimension(:), allocatable :: tvals
 
      call swiftest_util_get_vals_storage(self, idvals, tvals)

      call util_unique(idvals,self%idvals,self%idmap)
      self%nid = size(self%idvals)

      call util_unique(tvals,self%tvals,self%tmap)
      self%nt = size(self%tvals)

      return
   end subroutine swiftest_util_index_map_storage


   module subroutine swiftest_util_make_impactors_pl(self, idx)
      !! author: David A. Minton
      !!
      !! This is a simple wrapper function that is used to make a type-bound procedure using a subroutine whose interface is in the
      !! collision module, which must be defined first
      implicit none
      class(swiftest_pl),         intent(inout) :: self  !! Massive body object
      integer(I4B), dimension(:), intent(in)    :: idx !! Array holding the indices of the two bodies involved in the collision)

      call collision_resolve_make_impactors_pl(self, idx)

      return
   end subroutine swiftest_util_make_impactors_pl


   module subroutine swiftest_util_peri(n,m, r, v, atp, q, isperi)
      !! author: David A. Minton
      !!
      !! Helper function that does the pericenter passage computation for any body
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_peri.f90
      !! Adapted from Hal Levison's Swift routine util_mass_peri.f
      implicit none
      ! Arguments
      integer(I4B),                 intent(in)    :: n      !! Number of bodies
      real(DP),     dimension(:),   intent(in)    :: m      !! Mass term (mu for HELIO coordinates, and Gmtot for BARY)
      real(DP),     dimension(:,:), intent(in)    :: r      !! Position vectors (rh for HELIO coordinates, rb for BARY)
      real(DP),     dimension(:,:), intent(in)    :: v      !! Position vectors (vh for HELIO coordinates, rb for BARY)
      real(DP),     dimension(:),   intent(out)   :: atp    !! Semimajor axis 
      real(DP),     dimension(:),   intent(out)   :: q      !! Periapsis
      integer(I4B), dimension(:),   intent(inout) :: isperi !! Periapsis passage flag
      ! Internals
      integer(I4B) :: i
      real(DP), dimension(n) :: e !! Temporary, just to make use of the xv2aeq subroutine
      real(DP) :: vdotr

      do i = 1,n
         vdotr = dot_product(r(:,i),v(:,i))
         if (isperi(i) == -1) then
            if (vdotr >= 0.0) then
               isperi(i) = 0
               call swiftest_orbel_xv2aeq(m(i),r(1,i),r(2,i),r(3,i),v(1,i),v(2,i),v(3,i),atp(i),e(i),q(i))
            end if
         else
            if (vdotr > 0.0) then
               isperi(i) = -1
            else
               isperi(i) = 1
            end if
         end if
      end do

      return
   end subroutine swiftest_util_peri


   module subroutine swiftest_util_peri_body(self, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Determine nbody_system pericenter passages for bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_peri.f90
      !! Adapted from Hal Levison's Swift routine util_mass_peri.f
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i

      if (self%nbody == 0) return

      select type(self)
      class is (swiftest_pl)
         if (self%lfirst) self%isperi(:) = 0
      end select

      if (param%qmin_coord == "HELIO") then
         call swiftest_util_peri(self%nbody, self%mu, self%rh, self%vh, self%atp, self%peri, self%isperi)
      else 
         call swiftest_util_peri(self%nbody, [(nbody_system%Gmtot,i=1,self%nbody)], self%rb, self%vb, self%atp, self%peri, &
                                 self%isperi)
      end if

      return
   end subroutine swiftest_util_peri_body


   module subroutine swiftest_util_rearray_pl(self, nbody_system, param)
      !! Author: the Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Clean up the massive body structures to remove discarded bodies and add new bodies
      use symba
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      class(swiftest_pl), allocatable :: tmp !! The discarded body list.
      integer(I4B) :: i, npl, nadd, idnew1, idnew2, idold1, idold2
      integer(I8B) :: k, nenc_old, nencmin
      logical, dimension(:), allocatable :: lmask
      class(encounter_list), allocatable :: plplenc_old
      logical :: lencounter

      associate(pl => self, tp => nbody_system%tp, cb => nbody_system%cb, pl_adds => nbody_system%pl_adds)

         npl = pl%nbody
         nadd = pl_adds%nbody
         if (npl == 0) return
         ! Deallocate any temporary variables
         if (allocated(pl%rbeg)) deallocate(pl%rbeg)
         if (allocated(pl%rend)) deallocate(pl%rend)

         ! Remove the discards and destroy the list, as the nbody_system already tracks pl_discards elsewhere
         allocate(lmask(npl))
         lmask(1:npl) = pl%ldiscard(1:npl)
         if (count(lmask(:)) > 0) then
            allocate(tmp, mold=self)
            call pl%spill(tmp, lspill_list=lmask, ldestructive=.true.)
            npl = pl%nbody
            call tmp%setup(0,param)
            deallocate(tmp)
            deallocate(lmask)
         end if

         ! Add in any new bodies
         if (nadd > 0) then
            ! Append the adds to the main pl object
            call pl%append(pl_adds, lsource_mask=[(.true., i=1, nadd)])
            npl = pl%nbody
         end if

         if (npl == 0) then
            if (param%lmtiny_pl) pl%nplm = 0
            ! There are no more massive bodies. Reset the encounter lists and move on
            if (allocated(nbody_system%plpl_encounter)) call nbody_system%plpl_encounter%setup(0_I8B)
            if (allocated(nbody_system%pltp_encounter)) call nbody_system%pltp_encounter%setup(0_I8B)
            return
         end if

         ! Reset all of the status flags for the remaining bodies
         pl%status(1:npl) = ACTIVE
         do i = 1, npl
            call pl%info(i)%set_value(status="ACTIVE")
         end do
         pl%ldiscard(1:npl) = .false.
         pl%lcollision(1:npl) = .false.
         pl%lmask(1:npl) = .true.

         if (param%lmtiny_pl) then
            pl%lmtiny(1:npl) = pl%Gmass(1:npl) < param%GMTINY
            where(pl%lmtiny(1:npl))
               pl%info(1:npl)%particle_type = PL_TINY_TYPE_NAME 
            elsewhere
               pl%info(1:npl)%particle_type = PL_TYPE_NAME 
            end where
            pl%nplm = count(.not.pl%lmtiny(1:npl))
         end if

         ! Reindex the new list of bodies 
         select type(pl)
         class is (helio_pl)
            call pl%sort("mass", ascending=.false.)
         class is (whm_pl) 
            call pl%sort("ir3h", ascending=.false.)
         end select
         call pl%flatten(param)

         call pl%set_rhill(cb)

         ! Reset the kinship trackers
         call pl%reset_kinship([(i, i=1, npl)])

         if (allocated(nbody_system%plpl_encounter)) then
            ! Store the original plplenc list so we don't remove any of the original encounters
            nenc_old = nbody_system%plpl_encounter%nenc
            if (nenc_old > 0_I8B) then 
               allocate(plplenc_old, source=nbody_system%plpl_encounter)
               call plplenc_old%copy(nbody_system%plpl_encounter)
            end if

            ! Re-build the encounter list
            ! Be sure to get the level info if this is a SyMBA nbody_system
            select type(nbody_system)
            class is (symba_nbody_system)
            select type(pl)
            class is (symba_pl)
            select type(tp)
            class is (symba_tp)
               lencounter = pl%encounter_check(param, nbody_system, param%dt, nbody_system%irec) 
               if (tp%nbody > 0) then
                  lencounter = tp%encounter_check(param, nbody_system, param%dt, nbody_system%irec)
               end if
            end select
            end select
            end select

            ! Re-index the encounter list as the index values may have changed
            if (nenc_old > 0_I8B) then
               nencmin = min(nbody_system%plpl_encounter%nenc, plplenc_old%nenc) 
               nbody_system%plpl_encounter%nenc = nencmin
               do k = 1_I8B, nencmin
                  idnew1 = nbody_system%plpl_encounter%id1(k)
                  idnew2 = nbody_system%plpl_encounter%id2(k)
                  idold1 = plplenc_old%id1(k)
                  idold2 = plplenc_old%id2(k)
                  if ((idnew1 == idold1) .and. (idnew2 == idold2)) then
                     ! This is an encounter we already know about, so save the old information
                     nbody_system%plpl_encounter%lvdotr(k) = plplenc_old%lvdotr(k) 
                     nbody_system%plpl_encounter%lclosest(k) = plplenc_old%lclosest(k) 
                     nbody_system%plpl_encounter%status(k) = plplenc_old%status(k) 
                     nbody_system%plpl_encounter%r1(:,k) = plplenc_old%r1(:,k)
                     nbody_system%plpl_encounter%r2(:,k) = plplenc_old%r2(:,k)
                     nbody_system%plpl_encounter%v1(:,k) = plplenc_old%v1(:,k)
                     nbody_system%plpl_encounter%v2(:,k) = plplenc_old%v2(:,k)
                     nbody_system%plpl_encounter%tcollision(k) = plplenc_old%tcollision(k)
                     nbody_system%plpl_encounter%level(k) = plplenc_old%level(k)
                  else if (((idnew1 == idold2) .and. (idnew2 == idold1))) then
                     ! This is an encounter we already know about, but with the order reversed, so save the old information
                     nbody_system%plpl_encounter%lvdotr(k) = plplenc_old%lvdotr(k) 
                     nbody_system%plpl_encounter%lclosest(k) = plplenc_old%lclosest(k) 
                     nbody_system%plpl_encounter%status(k) = plplenc_old%status(k) 
                     nbody_system%plpl_encounter%r1(:,k) = plplenc_old%r2(:,k)
                     nbody_system%plpl_encounter%r2(:,k) = plplenc_old%r1(:,k)
                     nbody_system%plpl_encounter%v1(:,k) = plplenc_old%v2(:,k)
                     nbody_system%plpl_encounter%v2(:,k) = plplenc_old%v1(:,k)
                     nbody_system%plpl_encounter%tcollision(k) = plplenc_old%tcollision(k)
                     nbody_system%plpl_encounter%level(k) = plplenc_old%level(k)
                  end if
                  nbody_system%plpl_encounter%index1(k) = findloc(pl%id(1:npl), nbody_system%plpl_encounter%id1(k), dim=1)
                  nbody_system%plpl_encounter%index2(k) = findloc(pl%id(1:npl), nbody_system%plpl_encounter%id2(k), dim=1)
               end do
               if (allocated(lmask)) deallocate(lmask)
               allocate(lmask(nencmin))
               nenc_old = nencmin
               if (any(nbody_system%plpl_encounter%index1(1:nencmin) == 0) .or. &
                  any(nbody_system%plpl_encounter%index2(1:nencmin) == 0)) then
                  lmask(:) = nbody_system%plpl_encounter%index1(1:nencmin) /= 0 .and. &
                             nbody_system%plpl_encounter%index2(1:nencmin) /= 0
               else
                  return
               end if
               nencmin = count(lmask(:))
               nbody_system%plpl_encounter%nenc = nencmin
               if (nencmin > 0_I8B) then
                  nbody_system%plpl_encounter%index1(1:nencmin) = pack(nbody_system%plpl_encounter%index1(1:nenc_old), &
                                                                       lmask(1:nenc_old))
                  nbody_system%plpl_encounter%index2(1:nencmin) = pack(nbody_system%plpl_encounter%index2(1:nenc_old), &
                                                                       lmask(1:nenc_old))
                  nbody_system%plpl_encounter%id1(1:nencmin) = pack(nbody_system%plpl_encounter%id1(1:nenc_old), lmask(1:nenc_old))
                  nbody_system%plpl_encounter%id2(1:nencmin) = pack(nbody_system%plpl_encounter%id2(1:nenc_old), lmask(1:nenc_old))
                  nbody_system%plpl_encounter%lvdotr(1:nencmin) = pack(nbody_system%plpl_encounter%lvdotr(1:nenc_old), &
                                                                       lmask(1:nenc_old))
                  nbody_system%plpl_encounter%lclosest(1:nencmin) = pack(nbody_system%plpl_encounter%lclosest(1:nenc_old), &
                                                                         lmask(1:nenc_old))
                  nbody_system%plpl_encounter%status(1:nencmin) = pack(nbody_system%plpl_encounter%status(1:nenc_old), &
                                                                       lmask(1:nenc_old))
                  nbody_system%plpl_encounter%tcollision(1:nencmin) = pack(nbody_system%plpl_encounter%tcollision(1:nenc_old), &
                                                                           lmask(1:nenc_old))
                  nbody_system%plpl_encounter%level(1:nencmin) = pack(nbody_system%plpl_encounter%level(1:nenc_old), &
                                                                      lmask(1:nenc_old))
                  do i = 1, NDIM
                     nbody_system%plpl_encounter%r1(i, 1:nencmin) = pack(nbody_system%plpl_encounter%r1(i, 1:nenc_old), &
                                                                         lmask(1:nenc_old))
                     nbody_system%plpl_encounter%r2(i, 1:nencmin) = pack(nbody_system%plpl_encounter%r2(i, 1:nenc_old), &
                                                                         lmask(1:nenc_old))
                     nbody_system%plpl_encounter%v1(i, 1:nencmin) = pack(nbody_system%plpl_encounter%v1(i, 1:nenc_old), &
                                                                         lmask(1:nenc_old))
                     nbody_system%plpl_encounter%v2(i, 1:nencmin) = pack(nbody_system%plpl_encounter%v2(i, 1:nenc_old), &
                                                                         lmask(1:nenc_old))
                  end do
               end if
            end if
         end if
      end associate

      return
   end subroutine swiftest_util_rearray_pl


   module subroutine swiftest_util_rearray_tp(self, nbody_system, param)
      !! Author: David A. Minton
      !!
      !! Clean up the test particle structures to remove discarded bodies
      use symba
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      class(swiftest_tp), allocatable :: tmp !! The discarded body list.
      integer(I4B) :: i, ntp, npl
      integer(I8B) :: k, nenc
      logical, dimension(:), allocatable :: lmask
      logical :: lencounter

      associate(tp => self, pl => nbody_system%pl, cb => nbody_system%cb, pl_adds => nbody_system%pl_adds)

         ntp = tp%nbody
         if (ntp == 0) return
         npl = pl%nbody

         ! Remove the discards and destroy the list, as the nbody_system already tracks tp_discards elsewhere
         allocate(lmask(ntp))
         lmask(1:ntp) = tp%ldiscard(1:ntp)
         if (count(lmask(:)) > 0) then
            allocate(tmp, mold=self)
            call tp%spill(tmp, lspill_list=lmask, ldestructive=.true.)
            ntp = tp%nbody
            call tmp%setup(0,param)
            deallocate(tmp)
            deallocate(lmask)
         end if
         ntp = tp%nbody
         if (ntp == 0) then
            ! There are no more test particles. Reset the encounter list and move on
            if (allocated(nbody_system%pltp_encounter)) call nbody_system%pltp_encounter%setup(0_I8B)
            return
         end if

         ! Reset all of the status flags for the remaining bodies
         tp%status(1:ntp) = ACTIVE
         do i = 1, ntp
            call tp%info(i)%set_value(status="ACTIVE")
         end do
         tp%ldiscard(1:ntp) = .false.
         tp%lcollision(1:ntp) = .false.
         tp%lmask(1:ntp) = .true.

         if (allocated(nbody_system%pltp_encounter)) then
            ! Index values may have changed, so re-index the encounter list
            nenc = nbody_system%pltp_encounter%nenc
            do k = 1_I8B, nenc
               nbody_system%pltp_encounter%index1(k) = findloc(pl%id(1:npl), nbody_system%pltp_encounter%id1(k), dim=1)
               nbody_system%pltp_encounter%index2(k) = findloc(tp%id(1:ntp), nbody_system%pltp_encounter%id2(k), dim=1)
            end do

         end if

      end associate

      return
   end subroutine swiftest_util_rearray_tp


   module subroutine swiftest_util_rescale_system(self, param, mscale, dscale, tscale)
      !! author: David A. Minton
      !!
      !! Rescales an nbody system to a new set of units. Inputs are the multipliers on the mass (mscale), distance (dscale), and 
      !! time units (tscale). Rescales all united quantities in the nbody_system, as well as the mass conversion factors, 
      !! gravitational constant, and Einstein's constant in the parameter object.
      implicit none
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters. Returns with new values of the 
                                                            !! scale vactors and GU
      real(DP),                     intent(in)    :: mscale, dscale, tscale !! Scale factors for mass, distance, and time units 
      ! Internals
      real(DP) :: vscale

      param%MU2KG = param%MU2KG * mscale
      param%DU2M = param%DU2M * dscale
      param%TU2S = param%TU2S * tscale

      ! Calculate the G for the nbody_system units
      param%GU = GC / (param%DU2M**3 / (param%MU2KG * param%TU2S**2))

      if (param%lgr) then
         ! Calculate the inverse speed of light in the nbody_system units
         param%inv_c2 = einsteinC * param%TU2S / param%DU2M
         param%inv_c2 = (param%inv_c2)**(-2)
      end if

      vscale = dscale / tscale

      associate(cb => self%cb, pl => self%pl, npl => self%pl%nbody, tp => self%tp, ntp => self%tp%nbody)

         cb%mass = cb%mass / mscale
         cb%Gmass = param%GU * cb%mass 
         cb%radius = cb%radius / dscale
         cb%rb(:) = cb%rb(:) / dscale
         cb%vb(:) = cb%vb(:) / vscale
         cb%rot(:) = cb%rot(:) * tscale
         pl%mass(1:npl) = pl%mass(1:npl) / mscale
         pl%Gmass(1:npl) = param%GU * pl%mass(1:npl) 
         pl%radius(1:npl) = pl%radius(1:npl) / dscale
         pl%rh(:,1:npl) = pl%rh(:,1:npl) / dscale
         pl%vh(:,1:npl) = pl%vh(:,1:npl) / vscale
         pl%rb(:,1:npl) = pl%rb(:,1:npl) / dscale
         pl%vb(:,1:npl) = pl%vb(:,1:npl) / vscale
         pl%rot(:,1:npl) = pl%rot(:,1:npl) * tscale

      end associate


      return
   end subroutine swiftest_util_rescale_system 


   module subroutine swiftest_util_reset_kinship_pl(self, idx)
      !! author: David A. Minton
      !! 
      !! Resets the kinship status of bodies.
      !!
      implicit none
      class(swiftest_pl),         intent(inout) :: self !! SyMBA massive body object
      integer(I4B), dimension(:), intent(in)    :: idx  !! Index array of bodies to reset
      ! Internals
      integer(I4B) :: i, j

      self%kin(idx(:))%parent = idx(:)
      self%kin(idx(:))%nchild = 0
      do j = 1, size(idx(:))
         i = idx(j)
         if (allocated(self%kin(i)%child)) deallocate(self%kin(i)%child)
      end do

      return
   end subroutine swiftest_util_reset_kinship_pl


   module subroutine swiftest_util_resize_arr_info(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. Array will only be resized if has previously been allocated. 
      !! Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                         intent(in)    :: nnew !! New size
      ! Internals
      type(swiftest_particle_info), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already 
                                                                     !! allocated
      integer(I4B) :: nold !! Old size

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr)
      else
         nold = 0
      end if

      if (nnew == nold) return
      
      allocate(tmp(nnew))
      if (nnew > nold) then
         call swiftest_util_copy_particle_info_arr(arr(1:nold), tmp(1:nold))
      else
         call swiftest_util_copy_particle_info_arr(arr(1:nnew), tmp(1:nnew))
      end if

      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_resize_arr_info

  
   module subroutine swiftest_util_resize_arr_kin(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. Array will only be resized if has previously been allocated. 
      !! Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                      intent(in)    :: nnew !! New size
      ! Internals
      type(swiftest_kinship), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already 
                                                               !! allocated
      integer(I4B) :: nold !! Old size

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr)
      else
         nold = 0
      end if

      allocate(tmp(nnew))
      if (nnew > nold) then
         tmp(1:nold) = arr(1:nold)
      else
         tmp(1:nnew) = arr(1:nnew)
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_resize_arr_kin


   module subroutine swiftest_util_resize_body(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self  !! Swiftest body object
      integer(I4B),         intent(in)    :: nnew  !! New size neded

      call util_resize(self%info, nnew)
      call util_resize(self%id, nnew)
      call util_resize(self%status, nnew)
      call util_resize(self%lcollision, nnew)
      call util_resize(self%lencounter, nnew)
      call util_resize(self%ldiscard, nnew)
      call util_resize(self%lmask, nnew)
      call util_resize(self%mu, nnew)
      call util_resize(self%rh, nnew)
      call util_resize(self%vh, nnew)
      call util_resize(self%rb, nnew)
      call util_resize(self%vb, nnew)
      call util_resize(self%ah, nnew)
      call util_resize(self%aobl, nnew)
      call util_resize(self%atide, nnew)
      call util_resize(self%agr, nnew)
      call util_resize(self%ir3h, nnew)
      call util_resize(self%a, nnew)
      call util_resize(self%e, nnew)
      call util_resize(self%inc, nnew)
      call util_resize(self%capom, nnew)
      call util_resize(self%omega, nnew)
      call util_resize(self%capm, nnew)
      self%nbody = count(self%status(1:nnew) /= INACTIVE)

      return
   end subroutine swiftest_util_resize_body


   module subroutine swiftest_util_resize_pl(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest massive body against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self  !! Swiftest massive body object
      integer(I4B),       intent(in)    :: nnew  !! New size neded

      call swiftest_util_resize_body(self, nnew)

      call util_resize(self%mass, nnew)
      call util_resize(self%Gmass, nnew)
      call util_resize(self%rhill, nnew)
      call util_resize(self%renc, nnew)
      call util_resize(self%radius, nnew)
      call util_resize(self%rbeg, nnew)
      call util_resize(self%rend, nnew)
      call util_resize(self%vbeg, nnew)
      call util_resize(self%density, nnew)
      call util_resize(self%Ip, nnew)
      call util_resize(self%rot, nnew)
      call util_resize(self%k2, nnew)
      call util_resize(self%Q, nnew)
      call util_resize(self%tlag, nnew)
      call util_resize(self%kin, nnew)
      call util_resize(self%lmtiny, nnew)
      call util_resize(self%nplenc, nnew)
      call util_resize(self%ntpenc, nnew)



      if (allocated(self%k_plpl)) deallocate(self%k_plpl)

      return
   end subroutine swiftest_util_resize_pl


   module subroutine swiftest_util_resize_tp(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest test particle against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_tp), intent(inout) :: self  !! Swiftest test particle object
      integer(I4B),       intent(in)    :: nnew  !! New size neded

      call swiftest_util_resize_body(self, nnew)

      call util_resize(self%nplenc, nnew)
      call util_resize(self%isperi, nnew)
      call util_resize(self%peri, nnew)
      call util_resize(self%atp, nnew)

      return
   end subroutine swiftest_util_resize_tp


   module subroutine swiftest_util_set_beg_end_pl(self, rbeg, rend, vbeg)
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of rbeg, rend, and vbeg
      implicit none
      ! Arguments
      class(swiftest_pl),       intent(inout)          :: self !! Swiftest massive body object
      real(DP), dimension(:,:), intent(in),   optional :: rbeg, rend, vbeg

      if (present(rbeg)) then
         if (allocated(self%rbeg)) deallocate(self%rbeg)
         allocate(self%rbeg, source=rbeg)
      end if
      if (present(rend)) then
         if (allocated(self%rend)) deallocate(self%rend)
         allocate(self%rend, source=rend)
      end if
      if (present(vbeg)) then
         if (allocated(self%vbeg)) deallocate(self%vbeg)
         allocate(self%vbeg, source=vbeg)
      end if

      return
   end subroutine swiftest_util_set_beg_end_pl


   module subroutine swiftest_util_set_ir3h(self)
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
   end subroutine swiftest_util_set_ir3h


   module subroutine swiftest_util_set_msys(self)
      !! author: David A. Minton
      !!
      !! Sets the value of msys and the vector mass quantities based on the total mass of the nbody_system
      implicit none
      ! Arguments
      class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest nobdy nbody_system object

      self%Gmtot = self%cb%Gmass
      if (self%pl%nbody > 0) self%Gmtot = self%Gmtot + sum(self%pl%Gmass(1:self%pl%nbody), &
                                                           self%pl%status(1:self%pl%nbody) /= INACTIVE)

      return
   end subroutine swiftest_util_set_msys


   module subroutine swiftest_util_set_mu_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Computes G * (M + m) for each massive body
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object

      if (self%nbody > 0) self%mu(1:self%nbody) = cb%Gmass + self%Gmass(1:self%nbody)

      return
   end subroutine swiftest_util_set_mu_pl


   module subroutine swiftest_util_set_mu_tp(self, cb)
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
   end subroutine swiftest_util_set_mu_tp


   module subroutine swiftest_util_set_particle_info(self, name, particle_type, status, origin_type, origin_time, collision_id, &
                                            origin_rh, origin_vh, discard_time, discard_rh, discard_vh, discard_body_id)
      !! author: David A. Minton
      !!
      !! Sets one or more values of the particle information metadata object
      implicit none
      ! Arguments
      class(swiftest_particle_info), intent(inout)           :: self
      character(len=*),              intent(in),    optional :: name            !! Non-unique name
      character(len=*),              intent(in),    optional :: particle_type   !! String containing a description of the particle 
                                                                                !!  type (Central Body, Massive Body, Test Particle)
      character(len=*),              intent(in),    optional :: status          !! Particle status description: ACTIVE, MERGED, 
                                                                                !!  FRAGMENTED, etc.
      character(len=*),              intent(in),    optional :: origin_type     !! String containing a description of the origin of
                                                                                !!  the particle (e.g. Initial Conditions, 
                                                                                !!  Supercatastrophic, Disruption, etc.)
      real(DP),                      intent(in),    optional :: origin_time     !! The time of the particle's formation
      integer(I4B),                  intent(in),    optional :: collision_id    !! The ID fo the collision that formed the particle
      real(DP), dimension(:),        intent(in),    optional :: origin_rh       !! The heliocentric distance vector at the time of 
                                                                                !!  the particle's formation
      real(DP), dimension(:),        intent(in),    optional :: origin_vh       !! The heliocentric velocity vector at the time of 
                                                                                !!  the particle's formation
      real(DP),                      intent(in),    optional :: discard_time    !! The time of the particle's discard
      real(DP), dimension(:),        intent(in),    optional :: discard_rh      !! The heliocentric distance vector at the time of 
                                                                                !!  the particle's discard
      real(DP), dimension(:),        intent(in),    optional :: discard_vh      !! The heliocentric velocity vector at the time of 
                                                                                !! the particle's discard
      integer(I4B),                  intent(in),    optional :: discard_body_id !! The id of the other body involved in the discard
                                                                                !! (0 if no other body involved)
      ! Internals
      character(len=NAMELEN) :: lenstr
      character(len=:), allocatable :: fmtlabel

      write(lenstr, *) NAMELEN
      fmtlabel = "(A" // trim(adjustl(lenstr)) // ")"

      if (present(name)) then
         write(self%name, fmtlabel) name
         call swiftest_io_remove_nul_char(self%name)
      end if
      if (present(particle_type)) then
         write(self%particle_type, fmtlabel) particle_type
         call swiftest_io_remove_nul_char(self%particle_type)
      end if 
      if (present(status)) then
         write(self%status, fmtlabel) status
         call swiftest_io_remove_nul_char(self%status)
      end if
      if (present(origin_type)) then
         write(self%origin_type, fmtlabel) origin_type
         call swiftest_io_remove_nul_char(self%origin_type)
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
   end subroutine swiftest_util_set_particle_info


   module subroutine swiftest_util_set_renc_I4B(self, scale)
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
   end subroutine swiftest_util_set_renc_I4B


   module subroutine swiftest_util_set_renc_DP(self, scale)
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
   end subroutine swiftest_util_set_renc_DP


   module subroutine swiftest_util_set_rhill(self,cb)
      !! author: David A. Minton
      !!
      !! Sets the value of the Hill's radius
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self !! Swiftest massive body object
      class(swiftest_cb), intent(inout) :: cb   !! Swiftest central body object

      if (self%nbody == 0) return
      if (cb%Gmass <= tiny(1.0_DP)) return

      call self%xv2el(cb) 
      where(self%e(1:self%nbody) < 1.0_DP)
         self%rhill(1:self%nbody) = self%a(1:self%nbody) * (self%Gmass(1:self%nbody) / cb%Gmass / 3)**THIRD 
      elsewhere
         self%rhill(1:self%nbody) = (.mag.self%rh(:,1:self%nbody)) * (self%Gmass(1:self%nbody) / cb%Gmass / 3)**THIRD 
      end where

      return
   end subroutine swiftest_util_set_rhill


   module subroutine swiftest_util_set_rhill_approximate(self,cb)
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
   end subroutine swiftest_util_set_rhill_approximate

   module subroutine swiftest_util_setup_construct_system(nbody_system, param)

      !! author: David A. Minton
      !!
      !! Constructor for a Swiftest nbody system. Creates the nbody system object based on the user-input integrator
      !! 
      implicit none
      ! Arguments
      class(swiftest_nbody_system), allocatable, intent(inout) :: nbody_system !! Swiftest nbody_system object
      class(swiftest_parameters),                intent(inout) :: param        !! Current run configuration parameters
      select case(param%integrator)
      case (INT_BS)
         write(*,*) 'Bulirsch-Stoer integrator not yet enabled'
       case (INT_HELIO)
         allocate(helio_nbody_system :: nbody_system)
         select type(nbody_system)
         class is (helio_nbody_system)
            allocate(helio_cb :: nbody_system%cb)
            allocate(helio_pl :: nbody_system%pl)
            allocate(helio_tp :: nbody_system%tp)
            allocate(helio_pl :: nbody_system%pl_discards)
            allocate(helio_tp :: nbody_system%tp_discards)
         end select
         param%collision_model = "MERGE"
      case (INT_RA15)
         write(*,*) 'Radau integrator not yet enabled'
      case (INT_TU4)
         write(*,*) 'INT_TU4 integrator not yet enabled'
      case (INT_WHM)
         allocate(whm_nbody_system :: nbody_system)
         select type(nbody_system)
         class is (whm_nbody_system)
            allocate(whm_cb :: nbody_system%cb)
            allocate(whm_pl :: nbody_system%pl)
            allocate(whm_tp :: nbody_system%tp)
            allocate(whm_pl :: nbody_system%pl_discards)
            allocate(whm_tp :: nbody_system%tp_discards)
         end select
         param%collision_model = "MERGE"
      case (INT_RMVS)
         allocate(rmvs_nbody_system :: nbody_system)
         select type(nbody_system)
         class is (rmvs_nbody_system)
            allocate(rmvs_cb :: nbody_system%cb)
            allocate(rmvs_pl :: nbody_system%pl)
            allocate(rmvs_tp :: nbody_system%tp)
            allocate(rmvs_pl :: nbody_system%pl_discards)
            allocate(rmvs_tp :: nbody_system%tp_discards)
         end select
         param%collision_model = "MERGE"
      case (INT_SYMBA)
         allocate(symba_nbody_system :: nbody_system)
         select type(nbody_system)
         class is (symba_nbody_system)
            allocate(symba_cb :: nbody_system%cb)
            allocate(symba_pl :: nbody_system%pl)
            allocate(symba_tp :: nbody_system%tp)

            allocate(symba_tp :: nbody_system%tp_discards)
            allocate(symba_pl :: nbody_system%pl_adds)
            allocate(symba_pl :: nbody_system%pl_discards)

            allocate(symba_list_pltp     :: nbody_system%pltp_encounter)
            allocate(symba_list_plpl     :: nbody_system%plpl_encounter)
            allocate(collision_list_plpl :: nbody_system%plpl_collision)
            allocate(collision_list_pltp :: nbody_system%pltp_collision)
         end select
      case (INT_RINGMOONS)
         write(*,*) 'RINGMOONS-SyMBA integrator not yet enabled'
      case default
         write(*,*) 'Unkown integrator',param%integrator
         call base_util_exit(FAILURE,param%display_unit)
      end select
      nbody_system%lfirst_io = .true.
      nbody_system%lfirst_peri = .true.

      allocate(swiftest_particle_info :: nbody_system%cb%info)

      nbody_system%t = param%tstart

      return
   end subroutine swiftest_util_setup_construct_system


   module subroutine swiftest_util_setup_initialize_particle_info_system(self, param)
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
   end subroutine swiftest_util_setup_initialize_particle_info_system


   module subroutine swiftest_util_setup_initialize_system(self, system_history, param)
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system),              intent(inout) :: self           !! Swiftest nbody_system object
      class(swiftest_storage),      allocatable, intent(inout) :: system_history !! Stores the system history between output dumps
      class(swiftest_parameters),                intent(inout) :: param          !! Current run configuration parameters
      ! Internals
      type(encounter_storage)  :: encounter_history
      type(collision_storage)  :: collision_history

      call encounter_history%setup(4096)
      call collision_history%setup(4096)

      if (allocated(system_history)) then
         call system_history%dealloc()
         deallocate(system_history)
      end if
      allocate(swiftest_storage :: system_history)
      call system_history%setup(param%dump_cadence)
      allocate(swiftest_netcdf_parameters :: system_history%nc)

      associate(nbody_system => self, cb => self%cb, pl => self%pl, tp => self%tp, nc => system_history%nc)
         call nbody_system%read_in(nc, param)
         call nbody_system%validate_ids(param)
         call nbody_system%set_msys()
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
            call nbody_system%init_particle_info(param)
         end if

         ! Write initial conditions to file
         nc%file_name = param%outfile
         call nbody_system%initialize_output_file(nc, param) 
         call nc%close()

         allocate(collision_basic :: nbody_system%collider)
         call nbody_system%collider%setup(nbody_system, param)

         if (param%lenc_save_trajectory .or. param%lenc_save_closest) then
            allocate(encounter_netcdf_parameters :: encounter_history%nc)
            select type(nc => encounter_history%nc)
            class is (encounter_netcdf_parameters)
               nc%file_name = ENCOUNTER_OUTFILE
               if (.not.param%lrestart) then
                  call nc%initialize(param)
                  call nc%close()
               end if
            end select
            allocate(nbody_system%encounter_history, source=encounter_history)
         end if
        
         allocate(collision_netcdf_parameters :: collision_history%nc)
         select type(nc => collision_history%nc)
         class is (collision_netcdf_parameters)
            nc%file_name = COLLISION_OUTFILE
            if (param%lrestart) then
               call nc%open(param) ! This will find the nc%max_idslot variable
            else
               call nc%initialize(param)
            end if
            call nc%close()
            nbody_system%collider%maxid_collision = nc%max_idslot
         end select

         allocate(nbody_system%collision_history, source=collision_history)        

      end associate

      return
   end subroutine swiftest_util_setup_initialize_system


   module subroutine swiftest_util_setup_body(self, n, param)
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest particle class. Allocates space for all particles and
      !! initializes all components with a value.
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
      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))
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
      self%isperi(:) = 1
      self%peri(:)   = 0.0_DP
      self%atp(:)    = 0.0_DP

      if (param%lnon_spherical_cb) then
         allocate(self%aobl(NDIM, n))
         self%aobl(:,:) = 0.0_DP
      end if
      if (param%ltides) then
         allocate(self%atide(NDIM, n))
         self%atide(:,:) = 0.0_DP
      end if
      if (param%lgr) then
         allocate(self%agr(NDIM, n))
         self%agr(:,:) = 0.0_DP
      end if

      return
   end subroutine swiftest_util_setup_body


   module subroutine swiftest_util_setup_pl(self, n, param)
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest massive body class. Allocates space for all particles and
      !! initializes all components with a value. 
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter
      ! Internals
      integer(I4B) :: i

      !! Call allocation method for parent class
      !! The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call swiftest_util_setup_body(self, n, param)
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
         allocate(self%kin(n))

         self%nplenc(:) = 0
         self%ntpenc(:) = 0
         self%radius(:) = 0.0_DP
         self%density(:) = 1.0_DP
         call self%reset_kinship([(i, i=1, n)])
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
   end subroutine swiftest_util_setup_pl
   

   module subroutine swiftest_util_setup_tp(self, n, param)
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest test particle particle class. Allocates space for 
      !! all particles and initializes all components with a value. 
      implicit none
      ! Arguments
      class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !! Call allocation method for parent class
      !! The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call swiftest_util_setup_body(self, n, param)
      if (n == 0) return

      allocate(self%nplenc(n))

      self%npltp = 0_I8B
      self%nplenc(:) = 0

      return
   end subroutine swiftest_util_setup_tp
 

   module subroutine swiftest_util_snapshot_system(self, param, nbody_system, t, arg)
      !! author: David A. Minton
      !!
      !! Takes a snapshot of the nbody_system for later file storage
      implicit none
      ! Arguments
      class(swiftest_storage),      intent(inout)        :: self         !! Swiftest storage object
      class(swiftest_parameters),   intent(inout)        :: param        !! Current run configuration parameters
      class(swiftest_nbody_system), intent(inout)        :: nbody_system !! Swiftest nbody system object to store
      real(DP),                     intent(in), optional :: t            !! Time of snapshot if different from nbody_system time
      character(*),                 intent(in), optional :: arg          !! Optional argument (needed for extended storage type used
                                                                         !!  in collision snapshots)
      ! Internals
      class(swiftest_nbody_system), allocatable :: snapshot

      ! To allow for runs to be restarted in a bit-identical way, we'll need to run the same coordinate conversion routines we would
      !  run upon restarting
      select type(pl => nbody_system%pl)
      class is (whm_pl)
         call pl%h2j(nbody_system%cb)
      class is (helio_pl)
         call pl%vh2vb(nbody_system%cb)
      end select

      select type(tp => nbody_system%tp)
      class is (helio_tp)
      select type(cb => nbody_system%cb)
      class is (helio_cb)
         call tp%vh2vb(vbcb = -cb%ptbeg)
      end select
      end select

      if (.not.param%lrhill_present) call nbody_system%pl%set_rhill(nbody_system%cb)

      ! Take a minimal snapshot wihout all of the extra storage objects
      allocate(snapshot, mold=nbody_system)
      allocate(snapshot%cb, source=nbody_system%cb )
      allocate(snapshot%pl, source=nbody_system%pl )
      allocate(snapshot%tp, source=nbody_system%tp )

      snapshot%t                 = nbody_system%t
      snapshot%GMtot             = nbody_system%GMtot
      snapshot%ke_orbit          = nbody_system%ke_orbit
      snapshot%ke_spin           = nbody_system%ke_spin
      snapshot%pe                = nbody_system%pe
      snapshot%be                = nbody_system%be
      snapshot%te                = nbody_system%te
      snapshot%oblpot            = nbody_system%oblpot
      snapshot%L_orbit           = nbody_system%L_orbit
      snapshot%L_spin            = nbody_system%L_spin
      snapshot%L_total           = nbody_system%L_total
      snapshot%ke_orbit_orig     = nbody_system%ke_orbit_orig
      snapshot%ke_spin_orig      = nbody_system%ke_spin_orig
      snapshot%pe_orig           = nbody_system%pe_orig
      snapshot%be_orig           = nbody_system%be_orig
      snapshot%E_orbit_orig      = nbody_system%E_orbit_orig
      snapshot%GMtot_orig        = nbody_system%GMtot_orig
      snapshot%L_total_orig      = nbody_system%L_total_orig
      snapshot%L_orbit_orig      = nbody_system%L_orbit_orig
      snapshot%L_spin_orig       = nbody_system%L_spin_orig
      snapshot%L_escape          = nbody_system%L_escape
      snapshot%GMescape          = nbody_system%GMescape
      snapshot%E_collisions      = nbody_system%E_collisions
      snapshot%E_untracked       = nbody_system%E_untracked
      snapshot%ke_orbit_error    = nbody_system%ke_orbit_error   
      snapshot%ke_spin_error     = nbody_system%ke_spin_error    
      snapshot%pe_error          = nbody_system%pe_error         
      snapshot%be_error          = nbody_system%be_error         
      snapshot%E_orbit_error     = nbody_system%E_orbit_error    
      snapshot%Ecoll_error       = nbody_system%Ecoll_error      
      snapshot%E_untracked_error = nbody_system%E_untracked_error
      snapshot%te_error          = nbody_system%te_error         
      snapshot%L_orbit_error     = nbody_system%L_orbit_error    
      snapshot%L_spin_error      = nbody_system%L_spin_error     
      snapshot%L_escape_error    = nbody_system%L_escape_error   
      snapshot%L_total_error     = nbody_system%L_total_error    
      snapshot%Mtot_error        = nbody_system%Mtot_error       
      snapshot%Mescape_error     = nbody_system%Mescape_error    
      snapshot%lbeg              = nbody_system%lbeg


      ! Store a snapshot of the nbody_system for posterity
      call base_util_snapshot_save(self, snapshot)
      self%nt = self%iframe
      self%nid = self%nid + 1 ! Central body
      if (allocated(nbody_system%pl)) self%nid = self%nid + nbody_system%pl%nbody
      if (allocated(nbody_system%tp)) self%nid = self%nid + nbody_system%tp%nbody
       
      return
   end subroutine swiftest_util_snapshot_system


   module subroutine swiftest_util_sort_body(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a Swiftest body structure in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self      !! Swiftest body object
      character(*),         intent(in)    :: sortby    !! Sorting attribute
      logical,              intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending 
                                                       !!  or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(body => self, n => self%nbody)
         select case(sortby)
         case("id")
            call util_sort(direction * body%id(1:n), ind)
         case("status")
            call util_sort(direction * body%status(1:n), ind)
         case("ir3h")
            call util_sort(direction * body%ir3h(1:n), ind)
         case("a")
            call util_sort(direction * body%a(1:n), ind)
         case("e")
            call util_sort(direction * body%e(1:n), ind)
         case("inc")
            call util_sort(direction * body%inc(1:n), ind)
         case("capom")
            call util_sort(direction * body%capom(1:n), ind)
         case("mu")
            call util_sort(direction * body%mu(1:n), ind)
         case("peri")
            call util_sort(direction * body%peri(1:n), ind)
         case("atp")
            call util_sort(direction * body%atp(1:n), ind)
         case("info","lfirst","nbody","ldiscard","lcollision","lencounter","rh","vh","rb","vb","ah","aobl","atide","agr","isperi")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not found!'
            return
         end select

         call body%rearrange(ind)

      end associate

      return
   end subroutine swiftest_util_sort_body


   module subroutine swiftest_util_sort_pl(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a Swiftest massive body object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self      !! Swiftest massive body object
      character(*),       intent(in)    :: sortby    !! Sorting attribute
      logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or
                                                     !!  descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(pl => self, npl => self%nbody)
         select case(sortby)
         case("Gmass","mass")
            call util_sort(direction * pl%Gmass(1:npl), ind)
         case("rhill")
            call util_sort(direction * pl%rhill(1:npl), ind)
         case("renc")
            call util_sort(direction * pl%renc(1:npl), ind)
         case("radius")
            call util_sort(direction * pl%radius(1:npl), ind)
         case("density")
            call util_sort(direction * pl%density(1:npl), ind)
         case("k2")
            call util_sort(direction * pl%k2(1:npl), ind)
         case("Q")
            call util_sort(direction * pl%Q(1:npl), ind)
         case("tlag")
            call util_sort(direction * pl%tlag(1:npl), ind)
         case("nplenc")
            call util_sort(direction * pl%nplenc(1:npl), ind)
         case("ntpenc")
            call util_sort(direction * pl%ntpenc(1:npl), ind)
         case("lmtiny", "nplm", "nplplm", "kin", "rbeg", "rend", "vbeg", "Ip", "rot", "k_plpl", "nplpl")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class
            call swiftest_util_sort_body(pl, sortby, ascending)
            return
         end select

         call pl%rearrange(ind)

      end associate

      return
   end subroutine swiftest_util_sort_pl


   module subroutine swiftest_util_sort_tp(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a Swiftest test particle object  in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(swiftest_tp), intent(inout) :: self      !! Swiftest test particle object
      character(*),       intent(in)    :: sortby    !! Sorting attribute
      logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or
                                                     !!  descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(tp => self, ntp => self%nbody)
         select case(sortby)
         case("nplenc")
            call util_sort(direction * tp%nplenc(1:ntp), ind)
         case default ! Look for components in the parent class
            call swiftest_util_sort_body(tp, sortby, ascending)
            return
         end select

         call tp%rearrange(ind)

      end associate

      return
   end subroutine swiftest_util_sort_tp


   module subroutine swiftest_util_sort_rearrange_body(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange Swiftest body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(swiftest_body),               intent(inout) :: self !! Swiftest body object
      integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body 
                                                                !!  (should contain all 1:n index values in the desired order)

      associate(n => self%nbody)
         call util_sort_rearrange(self%id,       ind, n)
         call util_sort_rearrange(self%lmask,    ind, n)
         call util_sort_rearrange(self%info,     ind, n)
         call util_sort_rearrange(self%status,   ind, n)
         call util_sort_rearrange(self%ldiscard, ind, n)
         call util_sort_rearrange(self%lcollision, ind, n)
         call util_sort_rearrange(self%lencounter, ind, n)
         call util_sort_rearrange(self%rh,       ind, n)
         call util_sort_rearrange(self%vh,       ind, n)
         call util_sort_rearrange(self%rb,       ind, n)
         call util_sort_rearrange(self%vb,       ind, n)
         call util_sort_rearrange(self%ah,       ind, n)
         call util_sort_rearrange(self%aobl,     ind, n)
         call util_sort_rearrange(self%agr,      ind, n)
         call util_sort_rearrange(self%atide,    ind, n)
         call util_sort_rearrange(self%ir3h,     ind, n)
         call util_sort_rearrange(self%isperi,   ind, n)
         call util_sort_rearrange(self%peri,     ind, n)
         call util_sort_rearrange(self%atp,      ind, n)
         call util_sort_rearrange(self%mu,       ind, n)
         call util_sort_rearrange(self%a,        ind, n)
         call util_sort_rearrange(self%e,        ind, n)
         call util_sort_rearrange(self%inc,      ind, n)
         call util_sort_rearrange(self%capom,    ind, n)
         call util_sort_rearrange(self%omega,    ind, n)
         call util_sort_rearrange(self%capm,     ind, n)
      end associate

      return
   end subroutine swiftest_util_sort_rearrange_body


   module subroutine swiftest_util_sort_rearrange_arr_info(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of particle information type in-place from an index list.
      implicit none
      ! Arguments
      type(swiftest_particle_info),  dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B),                  dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                                             intent(in)    :: n   !! Number of elements in arr & ind to rearrange
      ! Internals
      type(swiftest_particle_info),  dimension(:), allocatable :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)

      call swiftest_util_copy_particle_info_arr(arr, tmp, ind)
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_sort_rearrange_arr_info
 

   pure module subroutine swiftest_util_sort_rearrange_arr_kin(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of particle kinship type in-place from an index list.
      implicit none
      ! Arguments
      type(swiftest_kinship),  dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B),            dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                                       intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      type(swiftest_kinship),  dimension(:), allocatable :: tmp !! Temporary copy of array used during rearrange operation
      integer(I4B) :: i,j

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, source=arr)
      tmp(1:n) = arr(ind(1:n))

      do i = 1, n
         do j = 1, tmp(i)%nchild
            tmp(i)%child(j) = ind(tmp(i)%child(j))
         end do
      end do

      call move_alloc(tmp, arr)
      return
   end subroutine swiftest_util_sort_rearrange_arr_kin


   module subroutine swiftest_util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange Swiftest massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      class(swiftest_pl),               intent(inout) :: self !! Swiftest massive body object
      integer(I4B),       dimension(:), intent(in)    :: ind  !! Index array used to restructure the body 
                                                              !!  (should contain all 1:n index values in the desired order)

      associate(pl => self, npl => self%nbody)
         call util_sort_rearrange(pl%mass,    ind, npl)
         call util_sort_rearrange(pl%Gmass,   ind, npl)
         call util_sort_rearrange(pl%rhill,   ind, npl)
         call util_sort_rearrange(pl%renc,    ind, npl)
         call util_sort_rearrange(pl%radius,  ind, npl)
         call util_sort_rearrange(pl%density, ind, npl)
         call util_sort_rearrange(pl%rbeg,    ind, npl)
         call util_sort_rearrange(pl%vbeg,    ind, npl)
         call util_sort_rearrange(pl%Ip,      ind, npl)
         call util_sort_rearrange(pl%rot,     ind, npl)
         call util_sort_rearrange(pl%k2,      ind, npl)
         call util_sort_rearrange(pl%Q,       ind, npl)
         call util_sort_rearrange(pl%tlag,    ind, npl)
         call util_sort_rearrange(pl%kin,        ind, npl)
         call util_sort_rearrange(pl%lmtiny,     ind, npl)
         call util_sort_rearrange(pl%nplenc,     ind, npl)
         call util_sort_rearrange(pl%ntpenc,     ind, npl)

         if (allocated(pl%k_plpl)) deallocate(pl%k_plpl)

         call swiftest_util_sort_rearrange_body(pl, ind)
      end associate

      return
   end subroutine swiftest_util_sort_rearrange_pl


   module subroutine swiftest_util_sort_rearrange_tp(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange Swiftest massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(swiftest_tp),                 intent(inout) :: self !! Swiftest test particle object
      integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body 
                                                                !!  (should contain all 1:n index values in the desired order)

      associate(tp => self, ntp => self%nbody)
         call util_sort_rearrange(tp%nplenc,  ind, ntp)

         if (allocated(tp%k_pltp)) deallocate(tp%k_pltp)

         call swiftest_util_sort_rearrange_body(tp, ind)
      end associate

      return
   end subroutine swiftest_util_sort_rearrange_tp


   module subroutine swiftest_util_spill_arr_info(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of particle origin information types
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,                       dimension(:),              intent(in)   :: lspill_list  !! Logical array of bodies to spill 
                                                                                             !!  into the discardss
      logical,                                                  intent(in)   :: ldestructive !! Logical flag indicating whether or
                                                                                             !!  not this operation should alter the
                                                                                             !!  keeps array or not
      ! Internals
      integer(I4B) :: i, nspill, nkeep, nlist
      integer(I4B), dimension(:), allocatable :: idx
      type(swiftest_particle_info), dimension(:), allocatable :: tmp

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (size(keeps) < nkeep) return
      if (.not.allocated(discards)) then
         allocate(discards(nspill))
      else if (size(discards) /= nspill) then
         deallocate(discards)
         allocate(discards(nspill))
      end if

      allocate(idx(nspill))
      idx(:) = pack([(i, i = 1, nlist)], lspill_list)
      call swiftest_util_copy_particle_info_arr(keeps, discards, idx)
      if (ldestructive) then
         if (nkeep > 0) then
            deallocate(idx)
            allocate(idx(nkeep))
            allocate(tmp(nkeep))
            idx(:) = pack([(i, i = 1, nlist)], .not. lspill_list)
            call swiftest_util_copy_particle_info_arr(keeps, tmp, idx)
            call move_alloc(tmp, keeps)
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine swiftest_util_spill_arr_info


   module subroutine swiftest_util_spill_arr_kin(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of particle kinships
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,                dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the
                                                                                       !!  discardss
      logical,                                           intent(in)    :: ldestructive !! Logical flag indicating whether or not
                                                                                       !!  this operation should alter the keeps
                                                                                       !!  array or not
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist
      type(swiftest_kinship), dimension(:), allocatable :: tmp

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (size(keeps) < nkeep) return
      if (.not.allocated(discards)) then
         allocate(discards(nspill))
      else if (size(discards) /= nspill) then
         deallocate(discards)
         allocate(discards(nspill))
      end if

      discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
      if (ldestructive) then
         if (nkeep > 0) then
            allocate(tmp(nkeep))
            tmp(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
            call move_alloc(tmp, keeps)
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine swiftest_util_spill_arr_kin


   module subroutine swiftest_util_spill_body(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest generic particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_body),  intent(inout) :: self         !! Swiftest generic body object
      class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter
                                                           !!  body by removing the discard list
      ! Internals
      integer(I4B) :: nbody_old

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !! Spill all the common components
      associate(keeps => self)
         
         call util_spill(keeps%id,         discards%id,         lspill_list, ldestructive)
         call util_spill(keeps%info,       discards%info,       lspill_list, ldestructive)
         call util_spill(keeps%lmask,      discards%lmask,      lspill_list, ldestructive)
         call util_spill(keeps%status,     discards%status,     lspill_list, ldestructive)
         call util_spill(keeps%ldiscard,   discards%ldiscard,   lspill_list, ldestructive)
         call util_spill(keeps%lcollision, discards%lcollision, lspill_list, ldestructive)
         call util_spill(keeps%lencounter, discards%lencounter, lspill_list, ldestructive)
         call util_spill(keeps%mu,         discards%mu,         lspill_list, ldestructive)
         call util_spill(keeps%rh,         discards%rh,         lspill_list, ldestructive)
         call util_spill(keeps%vh,         discards%vh,         lspill_list, ldestructive)
         call util_spill(keeps%rb,         discards%rb,         lspill_list, ldestructive)
         call util_spill(keeps%vb,         discards%vb,         lspill_list, ldestructive)
         call util_spill(keeps%ah,         discards%ah,         lspill_list, ldestructive)
         call util_spill(keeps%aobl,       discards%aobl,       lspill_list, ldestructive)
         call util_spill(keeps%agr,        discards%agr,        lspill_list, ldestructive)
         call util_spill(keeps%atide,      discards%atide,      lspill_list, ldestructive)
         call util_spill(keeps%ir3h,       discards%ir3h,       lspill_list, ldestructive)
         call util_spill(keeps%isperi,     discards%isperi,     lspill_list, ldestructive)
         call util_spill(keeps%peri,       discards%peri,       lspill_list, ldestructive)
         call util_spill(keeps%atp,        discards%atp,        lspill_list, ldestructive)
         call util_spill(keeps%a,          discards%a,          lspill_list, ldestructive)
         call util_spill(keeps%e,          discards%e,          lspill_list, ldestructive)
         call util_spill(keeps%inc,        discards%inc,        lspill_list, ldestructive)
         call util_spill(keeps%capom,      discards%capom,      lspill_list, ldestructive)
         call util_spill(keeps%omega,      discards%omega,      lspill_list, ldestructive)
         call util_spill(keeps%capm,       discards%capm,       lspill_list, ldestructive)

         nbody_old = keeps%nbody

         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nbody values for both the keeps and discareds
         discards%nbody = count(lspill_list(1:nbody_old))
         if (ldestructive) keeps%nbody = nbody_old- discards%nbody
      end associate
     
      return
   end subroutine swiftest_util_spill_body


   module subroutine swiftest_util_spill_pl(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest massive body structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_pl),    intent(inout) :: self        !! Swiftest massive body object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter
                                                           !!  body by removing the discard list

      associate(keeps => self)
         select type (discards) !The standard requires us to select the type of both arguments in order to access all the components
         class is (swiftest_pl)
            !! Spill components specific to the massive body class
            call util_spill(keeps%mass,    discards%mass,    lspill_list, ldestructive)
            call util_spill(keeps%Gmass,   discards%Gmass,   lspill_list, ldestructive)
            call util_spill(keeps%rhill,   discards%rhill,   lspill_list, ldestructive)
            call util_spill(keeps%renc,    discards%renc,    lspill_list, ldestructive)
            call util_spill(keeps%radius,  discards%radius,  lspill_list, ldestructive)
            call util_spill(keeps%density, discards%density, lspill_list, ldestructive)
            call util_spill(keeps%rbeg,    discards%rbeg,    lspill_list, ldestructive)
            call util_spill(keeps%rend,    discards%rend,    lspill_list, ldestructive)
            call util_spill(keeps%vbeg,    discards%vbeg,    lspill_list, ldestructive)
            call util_spill(keeps%Ip,      discards%Ip,      lspill_list, ldestructive)
            call util_spill(keeps%rot,     discards%rot,     lspill_list, ldestructive)
            call util_spill(keeps%k2,      discards%k2,      lspill_list, ldestructive)
            call util_spill(keeps%Q,       discards%Q,       lspill_list, ldestructive)
            call util_spill(keeps%tlag,    discards%tlag,    lspill_list, ldestructive)
            call util_spill(keeps%kin,     discards%kin,     lspill_list, ldestructive)
            call util_spill(keeps%lmtiny,  discards%lmtiny,  lspill_list, ldestructive)
            call util_spill(keeps%nplenc,  discards%nplenc,  lspill_list, ldestructive)
            call util_spill(keeps%ntpenc,  discards%ntpenc,  lspill_list, ldestructive)

            if (ldestructive .and. allocated(keeps%k_plpl)) deallocate(keeps%k_plpl)

            call swiftest_util_spill_body(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on swiftest_pl'
         end select
      end associate

      return
   end subroutine swiftest_util_spill_pl   


   module subroutine swiftest_util_spill_tp(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(swiftest_tp),    intent(inout) :: self         !! Swiftest test particle object
      class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardse
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter
                                                           !!  body by removing the discard list

      associate(keeps => self, ntp => self%nbody)
         select type(discards)
         class is (swiftest_tp)
            !! Spill components specific to the test particle class
            call util_spill(keeps%nplenc,  discards%nplenc,  lspill_list, ldestructive)
            call swiftest_util_spill_body(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on swiftest_tp'
         end select
      end associate

      return
   end subroutine swiftest_util_spill_tp


   module subroutine swiftest_util_valid_id_system(self, param)
      !! author: David A. Minton
      !!
      !! Validate massive body and test particle ids
      !! If non-unique values detected, it will replace them
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_valid.f90
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i, nid
      integer(I4B), dimension(:), allocatable :: idarr, unique_idarr, idmap

      associate(cb => self%cb, pl => self%pl, npl => self%pl%nbody, tp => self%tp, ntp => self%tp%nbody, maxid => self%maxid)
         nid = 1 + npl+ ntp
         allocate(idarr(nid))
         ! Central body should always be id=0
         cb%id = 0
         idarr(1) = cb%id
         do i = 1, npl
            idarr(1+i) = pl%id(i)
         end do
         do i = 1, ntp
            idarr(1+npl+i) = tp%id(i)
         end do
         maxid = maxval(idarr)

         ! Check to see if the ids are unique
         call util_unique(idarr, unique_idarr, idmap)
         if (size(unique_idarr) == nid) return ! All id values are unique

         ! Fix any duplicate id values and update the maxid
         call util_sort(idmap)
         do i = 2, size(idmap)
            if (idmap(i) == idmap(i-1)) then
               maxid = maxid + 1
               if (i < 1 + npl) then
                  pl%id(i - 1) = maxid 
               else
                  tp%id(i - 1 - npl) = maxid
               end if
            end if
         end do

      end associate

      return
   end subroutine swiftest_util_valid_id_system


   module subroutine swiftest_util_version()
      !! author: David A. Minton
      !!
      !! Print program version information to terminale
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_version.f90
      implicit none
      write(*, 200) VERSION
      200 format(/, "************* Swiftest: Version ", f3.1, " *************",/,  &
            "Based off of Swifter:",/,                                             &
            "Authors:",/,                                                          &
            "    The Purdue University Swiftest Development team ",/,              &
            "    Lead by David A. Minton ",/,                                      &
            "    Carlisle Wishard, Jennifer Pouplin, Jacob Elliott, Dana Singh.",/,&
            "Please address comments and questions to:",/,                         &
            "    David A. Minton",/,                                               &
            "    Department Earth, Atmospheric, & Planetary Sciences ",/,          &
            "    Purdue University",/,                                             &
            "    550 Stadium Mall Drive",/,                                        &
            "    West Lafayette, Indiana 47907", /,                                &
            "    765-494-3292 ",/,                                                 &
            "    daminton@purdue.edu",/,                                           &
            "Special thanks to Hal Levison, Martin Duncan, and David Kaufmann",/,  &
            "for the original SWIFTER and SWIFT codes that made this possible.",/, &
            "************************************************", /)


      return
   end subroutine swiftest_util_version

end submodule s_swiftest_util