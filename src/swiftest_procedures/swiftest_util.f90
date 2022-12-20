!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_util
contains

   module subroutine swiftest_util_append_arr_char_string(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of character string type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      character(len=STRMAX), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                                     intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,               dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine swiftest_util_append_arr_char_string


   module subroutine swiftest_util_append_arr_DP(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      real(DP), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                        intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,  dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine swiftest_util_append_arr_DP


   module subroutine swiftest_util_append_arr_DPvec(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision vector type of size (NDIM, n) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: arr          !! Destination array 
      real(DP), dimension(:,:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                          intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,  dimension(:),                intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(NDIM,nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(1, nold + 1:nold + nnew) = pack(source(1,1:nsrc), lsource_mask(1:nsrc))
      arr(2, nold + 1:nold + nnew) = pack(source(2,1:nsrc), lsource_mask(1:nsrc))
      arr(3, nold + 1:nold + nnew) = pack(source(3,1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine swiftest_util_append_arr_DPvec


   module subroutine swiftest_util_append_arr_I4B(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of integer(I4B) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      integer(I4B), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                            intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,      dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine swiftest_util_append_arr_I4B


   module subroutine swiftest_util_append_arr_info(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of particle information type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      type(swiftest_particle_info), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                                            intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,                       dimension(:),             intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew, i
      integer(I4B), dimension(:), allocatable :: idx

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      allocate(idx(nnew))

      idx = pack([(i, i = 1, nsrc)], lsource_mask(1:nsrc))

      call util_copy_particle_info_arr(source(1:nsrc), arr(nold+1:nold+nnew), idx)

      return
   end subroutine swiftest_util_append_arr_info


   module subroutine swiftest_util_append_arr_logical(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of logical type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      logical, dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                       intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical, dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine swiftest_util_append_arr_logical


   module subroutine swiftest_util_append_body(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),  intent(inout) :: self         !! Swiftest body object
      class(swiftest_body),  intent(in)    :: source       !! Source object to append
      logical, dimension(:), intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nold, nsrc, nnew

      nold = self%nbody
      nsrc = source%nbody
      nnew = count(lsource_mask(1:nsrc))

      call util_append(self%info, source%info, nold, nsrc, lsource_mask)
      call util_append(self%id, source%id, nold, nsrc, lsource_mask)
      call util_append(self%status, source%status, nold, nsrc, lsource_mask)
      call util_append(self%ldiscard, source%ldiscard, nold, nsrc, lsource_mask)
      call util_append(self%lmask, source%lmask, nold, nsrc, lsource_mask)
      call util_append(self%mu, source%mu, nold, nsrc, lsource_mask)
      call util_append(self%rh, source%rh, nold, nsrc, lsource_mask)
      call util_append(self%vh, source%vh, nold, nsrc, lsource_mask)
      call util_append(self%rb, source%rb, nold, nsrc, lsource_mask)
      call util_append(self%vb, source%vb, nold, nsrc, lsource_mask)
      call util_append(self%ah, source%ah, nold, nsrc, lsource_mask)
      call util_append(self%aobl, source%aobl, nold, nsrc, lsource_mask)
      call util_append(self%atide, source%atide, nold, nsrc, lsource_mask)
      call util_append(self%agr, source%agr, nold, nsrc, lsource_mask)
      call util_append(self%ir3h, source%ir3h, nold, nsrc, lsource_mask)
      call util_append(self%a, source%a, nold, nsrc, lsource_mask)
      call util_append(self%e, source%e, nold, nsrc, lsource_mask)
      call util_append(self%inc, source%inc, nold, nsrc, lsource_mask)
      call util_append(self%capom, source%capom, nold, nsrc, lsource_mask)
      call util_append(self%omega, source%omega, nold, nsrc, lsource_mask)
      call util_append(self%capm, source%capm, nold, nsrc, lsource_mask)

      self%nbody = nold + nnew

      return
   end subroutine swiftest_util_append_body


   module subroutine swiftest_util_append_pl(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_pl),              intent(inout) :: self         !! Swiftest massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (swiftest_pl)
         associate(nold => self%nbody, nsrc => source%nbody)
            call util_append(self%mass, source%mass, nold, nsrc, lsource_mask)
            call util_append(self%Gmass, source%Gmass, nold, nsrc, lsource_mask)
            call util_append(self%rhill, source%rhill, nold, nsrc, lsource_mask)
            call util_append(self%renc, source%renc, nold, nsrc, lsource_mask)
            call util_append(self%radius, source%radius, nold, nsrc, lsource_mask)
            call util_append(self%rbeg, source%rbeg, nold, nsrc, lsource_mask)
            call util_append(self%rend, source%rend, nold, nsrc, lsource_mask)
            call util_append(self%vbeg, source%vbeg, nold, nsrc, lsource_mask)
            call util_append(self%density, source%density, nold, nsrc, lsource_mask)
            call util_append(self%Ip, source%Ip, nold, nsrc, lsource_mask)
            call util_append(self%rot, source%rot, nold, nsrc, lsource_mask)
            call util_append(self%k2, source%k2, nold, nsrc, lsource_mask)
            call util_append(self%Q, source%Q, nold, nsrc, lsource_mask)
            call util_append(self%tlag, source%tlag, nold, nsrc, lsource_mask)

            if (allocated(self%k_plpl)) deallocate(self%k_plpl)

            call util_append_body(self, source, lsource_mask)
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class swiftest_pl or its descendents"
         call util_exit(FAILURE)
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
      class(swiftest_tp),              intent(inout) :: self         !! Swiftest test particle object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (swiftest_tp)
         associate(nold => self%nbody, nsrc => source%nbody)
            call util_append(self%isperi, source%isperi, nold, nsrc, lsource_mask)
            call util_append(self%peri, source%peri, nold, nsrc, lsource_mask)
            call util_append(self%atp, source%atp, nold, nsrc, lsource_mask)

            call util_append_body(self, source, lsource_mask)
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class swiftest_tp or its descendents"
         call util_exit(FAILURE)
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
      integer(I4B) :: i

      if (self%nbody == 0) return
      associate(tp => self, ntp => self%nbody)
         do concurrent (i = 1:ntp, tp%status(i) /= INACTIVE)
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
      integer(I4B)          :: i

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         do concurrent (i = 1:npl, pl%status(i) /= INACTIVE)
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
      integer(I4B) :: i

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         do concurrent(i = 1:ntp, tp%status(i) /= INACTIVE)
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
      integer(I4B)              :: i

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         cb%vb(:) = 0.0_DP
         do i = npl, 1, -1
            cb%vb(:) = cb%vb(:) - pl%Gmass(i) * pl%vb(:, i) / cb%Gmass
         end do
         do concurrent(i = 1:npl)
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
      integer(I4B)  :: i
      real(DP)      :: Gmtot

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         Gmtot = cb%Gmass + sum(pl%Gmass(1:npl))
         cb%vb(:) = 0.0_DP
         do i = 1, npl
            cb%vb(:) = cb%vb(:) - pl%Gmass(i) * pl%vh(:, i) 
         end do
         cb%vb(:) = cb%vb(:) / Gmtot
         do concurrent(i = 1:npl)
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
      integer(I4B) :: i

      if (self%nbody == 0) return
      associate(tp => self, ntp => self%nbody)
         do concurrent (i = 1:ntp, tp%status(i) /= INACTIVE)
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
      class(swiftest_particle_info), dimension(:), intent(inout)          :: dest   !! Swiftest body object with particle metadata information object
      integer(I4B),                  dimension(:), intent(in),   optional :: idx    !! Optional array of indices to draw the source object
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
      !! Finalize the swiftest body object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_body),  intent(inout) :: self

      if (allocated(self%info)) deallocate(self%info)
      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%ldiscard)) deallocate(self%ldiscard)
      if (allocated(self%lcollision)) deallocate(self%lcollision)
      if (allocated(self%lencounter)) deallocate(self%lencounter)
      if (allocated(self%lmask)) deallocate(self%lmask)
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
      if (allocated(self%a)) deallocate(self%a)
      if (allocated(self%e)) deallocate(self%e)
      if (allocated(self%e)) deallocate(self%e)
      if (allocated(self%inc)) deallocate(self%inc)
      if (allocated(self%capom)) deallocate(self%capom)
      if (allocated(self%omega)) deallocate(self%omega)
      if (allocated(self%capm)) deallocate(self%capm)

      return
   end subroutine swiftest_util_dealloc_body


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
      !! Finalize the swiftest massive body object - deallocates all allocatables
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
      if (allocated(self%rot)) deallocate(self%rot)
      if (allocated(self%Ip)) deallocate(self%Ip)
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

      call util_dealloc_body(self)

      return
   end subroutine swiftest_util_dealloc_pl


   module subroutine swiftest_util_dealloc_tp(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest test particle object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_tp),  intent(inout) :: self !! Swiftest test particle object

      if (allocated(self%nplenc)) deallocate(self%nplenc)
      if (allocated(self%isperi)) deallocate(self%isperi)
      if (allocated(self%peri)) deallocate(self%peri)
      if (allocated(self%atp)) deallocate(self%atp)
      if (allocated(self%k_pltp)) deallocate(self%k_pltp)

      call util_dealloc_body(self)

      return
   end subroutine swiftest_util_dealloc_tp


   module subroutine swiftest_util_exit(code)
      !! author: David A. Minton
      !!
      !! Print termination message and exit program
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_exit.f90
      !! Adapted from Hal Levison's Swift routine util_exit.f
      implicit none
      ! Arguments
      integer(I4B), intent(in) :: code
      ! Internals
      character(*), parameter :: BAR = '("------------------------------------------------")'
      character(*), parameter :: SUCCESS_MSG = '(/, "Normal termination of Swiftest (version ", f3.1, ")")'
      character(*), parameter :: FAIL_MSG = '(/, "Terminating Swiftest (version ", f3.1, ") due to error!!")'
      character(*), parameter :: USAGE_MSG = '("Usage: swiftest [bs|helio|ra15|rmvs|symba|tu4|whm] <paramfile> [standard|compact|progress|NONE]")'
      character(*), parameter :: HELP_MSG  = USAGE_MSG

      select case(code)
      case(SUCCESS)
         write(*, SUCCESS_MSG) VERSION_NUMBER
         write(*, BAR)
      case(USAGE) 
         write(*, USAGE_MSG)
      case(HELP)
         write(*, HELP_MSG)
      case default
         write(*, FAIL_MSG) VERSION_NUMBER
         write(*, BAR)
         error stop
      end select

      stop

   end subroutine swiftest_util_exit


   module subroutine swiftest_util_fill_arr_char_string(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of type character strings
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      character(len=STRMAX), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,               dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))

      return
   end subroutine swiftest_util_fill_arr_char_string


   module subroutine swiftest_util_fill_arr_DP(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of type DP
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      real(DP), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,  dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))

      return
   end subroutine swiftest_util_fill_arr_DP

   module subroutine swiftest_util_fill_arr_DPvec(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of DP vectors with shape (NDIM, n)
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      real(DP), dimension(:,:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,  dimension(:),                intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! Internals
      integer(I4B) :: i

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      do i = 1, NDIM
         keeps(i,:) = unpack(keeps(i,:),   .not.lfill_list(:), keeps(i,:))
         keeps(i,:) = unpack(inserts(i,:),      lfill_list(:), keeps(i,:))
      end do

      return
   end subroutine swiftest_util_fill_arr_DPvec

   module subroutine swiftest_util_fill_arr_I4B(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of type I4B
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      integer(I4B), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,      dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))

      return
   end subroutine swiftest_util_fill_arr_I4B


   module subroutine swiftest_util_fill_arr_info(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of particle origin information types
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      type(swiftest_particle_info), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,                      dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! Internals
      integer(I4B), dimension(:), allocatable  :: insert_idx
      integer(I4B) :: i, nkeep, ninsert

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      nkeep = size(keeps)
      ninsert = count(lfill_list)

      allocate(insert_idx(ninsert))

      insert_idx(:) = pack([(i, i = 1, nkeep)], lfill_list)
      call util_copy_particle_info_arr(inserts, keeps, insert_idx)

      return
   end subroutine swiftest_util_fill_arr_info


   module subroutine swiftest_util_fill_arr_logical(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of logicals
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      logical, dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical, dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))

      return
   end subroutine swiftest_util_fill_arr_logical


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
      !> Fill all the common components
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
            !> Fill components specific to the massive body class
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
            
            call util_fill_body(keeps, inserts, lfill_list)
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
            !> Spill components specific to the test particle class
            call util_fill(keeps%nplenc,  inserts%nplenc,  lfill_list)

            call util_fill_body(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on swiftest_tp'
         end select
      end associate

      return
   end subroutine swiftest_util_fill_tp


   module subroutine swiftest_util_final_storage(self)
      !! author: David A. Minton
      !!
      !! Finalizer for the storage data type
      implicit none
      ! Arguments
      type(swiftest_storage(*)) :: self
      ! Internals
      integer(I4B) :: i

      do i = 1, self%nframes
         if (allocated(self%frame(i)%item)) deallocate(self%frame(i)%item)
      end do

      return
   end subroutine swiftest_util_final_storage


   module subroutine swiftest_util_final_system(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest nbody system object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_nbody_system),  intent(inout) :: self !! Swiftest nbody system object

      if (allocated(self%cb)) deallocate(self%cb)
      if (allocated(self%pl)) deallocate(self%pl)
      if (allocated(self%tp)) deallocate(self%tp)
      if (allocated(self%tp_discards)) deallocate(self%tp_discards)
      if (allocated(self%pl_discards)) deallocate(self%pl_discards)

      return
   end subroutine swiftest_util_final_system


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
      !! Turns i,j indices into k index for use in the Euclidean distance matrix for pl-pl interactions for a Swiftest massive body object
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
      integer(I4B) :: i, j, err
      integer(I8B) :: k, npl

      npl = int(self%nbody, kind=I8B)
      associate(nplpl => self%nplpl)
         nplpl = npl * (npl - 1_I8B) / 2_I8B ! number of entries in a strict lower triangle, npl x npl
         if (param%lflatten_interactions) then
            if (allocated(self%k_plpl)) deallocate(self%k_plpl) ! Reset the index array if it's been set previously
            allocate(self%k_plpl(2, nplpl), stat=err)
            if (err /=0) then ! An error occurred trying to allocate this big array. This probably means it's too big to fit in memory, and so we will force the run back into triangular mode
               param%lflatten_interactions = .false.
            else
               do concurrent (i=1:npl, j=1:npl, j>i)
                  call util_flatten_eucl_ij_to_k(self%nbody, i, j, k)
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
      integer(I8B) :: i, j, counter, npl, ntp

      ntp = int(self%nbody, kind=I8B)
      npl = int(pl%nbody, kind=I8B)
      associate(npltp => self%npltp)
         npltp = npl * ntp
         if (allocated(self%k_pltp)) deallocate(self%k_pltp) ! Reset the index array if it's been set previously
         allocate(self%k_pltp(2, npltp))
         do i = 1_I8B, npl
            counter = (i - 1_I8B) * npl + 1_I8B
            do j = 1_I8B,  ntp
               self%k_pltp(1, counter) = i
               self%k_pltp(2, counter) = j
               counter = counter + 1_I8B
            end do
         end do
      end associate

      return
   end subroutine swiftest_util_flatten_eucl_pltp


   module subroutine swiftest_util_get_energy_momentum_system(self, param)
      !! author: David A. Minton
      !!
      !! Compute total system angular momentum vector and kinetic, potential and total system energy
      !!  
      !! Adapted from David E. Kaufmann Swifter routine symba_energy_eucl.f90
      !!  
      !! Adapted from Martin Duncan's Swift routine anal_energy.f
      implicit none
      class(swiftest_nbody_system), intent(inout) :: self     !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param    !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i
      real(DP) :: kecb, kespincb
      real(DP), dimension(self%pl%nbody) :: kepl, kespinpl
      real(DP), dimension(self%pl%nbody) :: Lplorbitx, Lplorbity, Lplorbitz
      real(DP), dimension(self%pl%nbody) :: Lplspinx, Lplspiny, Lplspinz
      real(DP), dimension(NDIM) :: Lcborbit, Lcbspin
      real(DP) :: hx, hy, hz

      associate(system => self, pl => self%pl, npl => self%pl%nbody, cb => self%cb)
         system%Lorbit(:) = 0.0_DP
         system%Lspin(:) = 0.0_DP
         system%Ltot(:) = 0.0_DP
         system%ke_orbit = 0.0_DP
         system%ke_spin = 0.0_DP

         kepl(:) = 0.0_DP
         Lplorbitx(:) = 0.0_DP
         Lplorbity(:) = 0.0_DP
         Lplorbitz(:) = 0.0_DP
         Lplspinx(:) = 0.0_DP
         Lplspiny(:) = 0.0_DP
         Lplspinz(:) = 0.0_DP

         pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE

         system%GMtot = cb%Gmass + sum(pl%Gmass(1:npl), pl%lmask(1:npl)) 
         kecb = cb%mass * dot_product(cb%vb(:), cb%vb(:))
         Lcborbit(:) = cb%mass * (cb%rb(:) .cross. cb%vb(:))

         do concurrent (i = 1:npl, pl%lmask(i))
            hx = pl%rb(2,i) * pl%vb(3,i) - pl%rb(3,i) * pl%vb(2,i)
            hy = pl%rb(3,i) * pl%vb(1,i) - pl%rb(1,i) * pl%vb(3,i)
            hz = pl%rb(1,i) * pl%vb(2,i) - pl%rb(2,i) * pl%vb(1,i)

            ! Angular momentum from orbit 
            Lplorbitx(i) = pl%mass(i) * hx
            Lplorbity(i) = pl%mass(i) * hy
            Lplorbitz(i) = pl%mass(i) * hz

            ! Kinetic energy from orbit
            kepl(i) = pl%mass(i) * dot_product(pl%vb(:,i), pl%vb(:,i)) 
         end do

         if (param%lrotation) then
            kespincb = cb%mass * cb%Ip(3) * cb%radius**2 * dot_product(cb%rot(:), cb%rot(:))

            ! For simplicity, we always assume that the rotation pole is the 3rd principal axis
            Lcbspin(:) = cb%Ip(3) * cb%mass * cb%radius**2 * cb%rot(:)

            do concurrent (i = 1:npl, pl%lmask(i))
               ! Currently we assume that the rotation pole is the 3rd principal axis
               ! Angular momentum from spin
               Lplspinx(i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(1,i)
               Lplspiny(i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(2,i)  
               Lplspinz(i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * pl%rot(3,i)  

               ! Kinetic energy from spin
               kespinpl(i) = pl%mass(i) * pl%Ip(3,i) * pl%radius(i)**2 * dot_product(pl%rot(:,i), pl%rot(:,i))
            end do
         else
            kespincb = 0.0_DP
            kespinpl(:) = 0.0_DP
         end if
  
         if (param%lflatten_interactions) then
            call util_get_energy_potential_flat(npl, pl%nplpl, pl%k_plpl, pl%lmask, cb%Gmass, pl%Gmass, pl%mass, pl%rb, system%pe)
         else
            call util_get_energy_potential_triangular(npl, pl%lmask, cb%Gmass, pl%Gmass, pl%mass, pl%rb, system%pe)
         end if

         ! Potential energy from the oblateness term
         if (param%loblatecb) then
            call system%obl_pot()
            system%pe = system%pe + system%oblpot
         end if

         system%ke_orbit = 0.5_DP * (kecb + sum(kepl(1:npl), pl%lmask(1:npl)))
         if (param%lrotation) system%ke_spin = 0.5_DP * (kespincb + sum(kespinpl(1:npl), pl%lmask(1:npl)))
   
         system%Lorbit(1) = Lcborbit(1) + sum(Lplorbitx(1:npl), pl%lmask(1:npl)) 
         system%Lorbit(2) = Lcborbit(2) + sum(Lplorbity(1:npl), pl%lmask(1:npl)) 
         system%Lorbit(3) = Lcborbit(3) + sum(Lplorbitz(1:npl), pl%lmask(1:npl)) 
  
         if (param%lrotation) then
            system%Lspin(1) = Lcbspin(1) + sum(Lplspinx(1:npl), pl%lmask(1:npl)) 
            system%Lspin(2) = Lcbspin(2) + sum(Lplspiny(1:npl), pl%lmask(1:npl)) 
            system%Lspin(3) = Lcbspin(3) + sum(Lplspinz(1:npl), pl%lmask(1:npl)) 
         end if

         system%te = system%ke_orbit + system%ke_spin + system%pe
         system%Ltot(:) = system%Lorbit(:) + system%Lspin(:)
      end associate

      return
   end subroutine swiftest_util_get_energy_momentum_system


   subroutine swiftest_util_get_energy_potential_flat(npl, nplpl, k_plpl, lmask, GMcb, Gmass, mass, rb, pe)
      !! author: David A. Minton
      !!
      !! Compute total system potential energy
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

      do concurrent(i = 1:npl, lmask(i))
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
   end subroutine swiftest_util_get_energy_potential_flat


   subroutine swiftest_util_get_energy_potential_triangular(npl, lmask, GMcb, Gmass, mass, rb, pe)
      !! author: David A. Minton
      !!
      !! Compute total system potential energy
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

      do concurrent(i = 1:npl, lmask(i))
         pecb(i) = -GMcb * mass(i) / norm2(rb(:,i)) 
      end do

      pe = 0.0_DP
      !$omp parallel do default(private) schedule(static)&
      !$omp shared(lmask, Gmass, mass, rb) &
      !$omp firstprivate(npl) &
      !$omp reduction(+:pe) 
      do i = 1, npl
         if (lmask(i)) then
            do concurrent(j = i+1:npl, lmask(i) .and. lmask(j))
               pepl(j) = - (Gmass(i) * mass(j)) / norm2(rb(:, i) - rb(:, j))
            end do
            pe = pe + sum(pepl(i+1:npl), lmask(i+1:npl))
         end if
      end do
      !$omp end parallel do
      pe = pe + sum(pecb(1:npl), lmask(1:npl))

      return
   end subroutine swiftest_util_get_energy_potential_triangular


   module subroutine swiftest_util_index_array(ind_arr, n)
      !! author: David A. Minton
      !!
      !! Creates or resizes an index array of size n where ind_arr = [1, 2, ... n].
      !! This subroutine swiftest_assumes that if ind_arr is already allocated, it is a pre-existing index array of a different size.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind_arr !! Index array. Input is a pre-existing index array where n /= size(ind_arr). Output is a new index array ind_arr = [1, 2, ... n]
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
      class(swiftest_storage(*)),              intent(in)  :: self   !! Swiftest storage object
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


   module subroutine swiftest_util_index_map_storage(self)
      !! author: David A. Minton
      !!
      !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      implicit none
      ! Arguments
      class(swiftest_storage(*)), intent(inout) :: self  !! Swiftest storage object
      ! Internals
      integer(I4B), dimension(:), allocatable :: idvals
      real(DP), dimension(:), allocatable :: tvals
 
      call util_get_vals_storage(self, idvals, tvals)

      call util_unique(idvals,self%idvals,self%idmap)
      self%nid = size(self%idvals)

      call util_unique(tvals,self%tvals,self%tmap)
      self%nt = size(self%tvals)

      return
   end subroutine swiftest_util_index_map_storage

   module subroutine swiftest_util_minimize_bfgs(f, N, x0, eps, maxloop, lerr, x1)
      !! author: David A. Minton
      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
      !! This function implements the Broyden-Fletcher-Goldfarb-Shanno method to determine the minimum of a function of N variables.  
      !! It recieves as input:
      !!   f%eval(x) : lambda function object containing the objective function as the eval metho
      !!   N       : Number of variables of function f
      !!   x0      : Initial starting value of x
      !!   eps     : Accuracy of 1 - dimensional minimization at each step
      !!   maxloop : Maximum number of loops to attempt to find a solution
      !! The outputs include
      !!   lerr :  Returns .true. if it could not find the minimum
      !! Returns
      !!   x1   :  Final minimum (all 0 if none found)
      !!   0 = No miniumum found
      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      integer(I4B),           intent(in)    :: N
      class(lambda_obj),      intent(inout) :: f
      real(DP), dimension(:), intent(in)    :: x0
      real(DP),               intent(in)    :: eps
      integer(I4B),           intent(in)    :: maxloop
      logical,                intent(out)   :: lerr
      ! Result
      real(DP), dimension(:), intent(out), allocatable :: x1
      ! Internals
      integer(I4B) ::  i, j, k, l, conv
      real(DP), parameter     :: graddelta = 1e-4_DP !! Delta x for gradient calculations
      real(DP), dimension(N) :: S               !! Direction vectors 
      real(DP), dimension(N,N) :: H             !! Approximated inverse Hessian matrix 
      real(DP), dimension(N) :: grad1           !! gradient of f 
      real(DP), dimension(N) :: grad0           !! old value of gradient 
      real(DP) :: astar                         !! 1D minimized value 
      real(DP), dimension(N) :: y, P
      real(DP), dimension(N,N) :: PP, PyH, HyP
      real(DP), save :: yHy, Py
      type(ieee_status_type) :: original_fpe_status
      logical, dimension(:), allocatable :: fpe_flag 

      call ieee_get_status(original_fpe_status) ! Save the original floating point exception status
      call ieee_set_flag(ieee_all, .false.) ! Set all flags to quiet
      allocate(fpe_flag(size(ieee_usual)))

      lerr = .false.
      allocate(x1, source=x0)
      ! Initialize approximate Hessian with the identity matrix (i.e. begin with method of steepest descent) 
      ! Get initial gradient and initialize arrays for updated values of gradient and x
      H(:,:) = reshape([((0._DP, i=1, j-1), 1._DP, (0._DP, i=j+1, N), j=1, N)], [N,N])  
      grad0 = gradf(f, N, x0(:), graddelta, lerr)
      if (lerr) then
         call ieee_set_status(original_fpe_status)
         return
      end if
      grad1(:) = grad0(:)
      do i = 1, maxloop 
         !check for convergence
         conv = count(abs(grad1(:)) > eps)
         if (conv == 0) exit 
         S(:) = -matmul(H(:,:), grad1(:))
         astar = minimize1D(f, x1, S, N, graddelta, lerr)
         if (lerr) exit
         ! Get new x values 
         P(:) = astar * S(:) 
         x1(:) = x1(:) + P(:)
         ! Calculate new gradient
         grad0(:) = grad1(:)
         grad1 = gradf(f, N, x1, graddelta, lerr)
         y(:) = grad1(:) - grad0(:)
         Py = sum(P(:) * y(:))
         ! set up factors for H matrix update 
         yHy = 0._DP
         !$omp do simd schedule(static)&
         !$omp firstprivate(N, y, H) &
         !$omp reduction(+:yHy)
         do k = 1, N 
            do j = 1, N
               yHy = yHy + y(j) * H(j,k) * y(k)
            end do
         end do
         !$omp end do simd
         ! prevent divide by zero (convergence) 
         if (abs(Py) < tiny(Py)) exit
         ! set up update 
         PyH(:,:) = 0._DP
         HyP(:,:) = 0._DP
         !$omp parallel do default(private) schedule(static)&
         !$omp shared(N, PP, P, y, H) &
         !$omp reduction(+:PyH, HyP)
         do k = 1, N 
            do j = 1, N
               PP(j, k) = P(j) * P(k)
               do l = 1, N
                  PyH(j, k) = PyH(j, k) + P(j) * y(l) * H(l,k)
                  HyP(j, k) = HyP(j, k) + P(k) * y(l) * H(j,l)
               end do
            end do
         end do
         !$omp end parallel do 
         ! update H matrix 
         H(:,:) = H(:,:) + ((1._DP - yHy / Py) * PP(:,:) - PyH(:,:) - HyP(:,:)) / Py
         ! Normalize to prevent it from blowing up if it takes many iterations to find a solution
         H(:,:) = H(:,:) / norm2(H(:,:))
         ! Stop everything if there are any exceptions to allow the routine to fail gracefully
         call ieee_get_flag(ieee_usual, fpe_flag)
         if (any(fpe_flag)) exit 
         if (i == maxloop) then
            lerr = .true.
         end if
      end do
      call ieee_get_flag(ieee_usual, fpe_flag)
      lerr = lerr .or. any(fpe_flag)  
      call ieee_set_status(original_fpe_status)

      return 

      contains

         function gradf(f, N, x1, dx, lerr) result(grad)
            !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
            !! Purpose:  Estimates the gradient of a function using a central difference
            !! approximation
            !! Inputs:
            !!   f%eval(x) : lambda function object containing the objective function as the eval metho
            !!   N    :  number of variables N
            !!   x1   :  x value array
            !!   dx   :  step size to use when calculating derivatives
            !! Outputs: 
            !!   lerr : .true. if an error occurred. Otherwise returns .false.
            !! Returns
            !!   grad :  N sized array containing estimated gradient of f at x1
            !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
            implicit none
            ! Arguments
            integer(I4B),           intent(in)    :: N
            class(lambda_obj),      intent(inout) :: f
            real(DP), dimension(:), intent(in)    :: x1
            real(DP),               intent(in)    :: dx
            logical,                intent(out)   :: lerr
            ! Result
            real(DP), dimension(N)                :: grad
            ! Internals
            integer(I4B) :: i, j
            real(DP), dimension(N) :: xp, xm
            real(DP) :: fp, fm
            logical :: lerrp, lerrm

            do i = 1, N
               do j = 1, N
                  if (j == i) then
                     xp(j) = x1(j) + dx
                     xm(j) = x1(j) - dx
                  else
                     xp(j) = x1(j)
                     xm(j) = x1(j)
                  end if
               end do
               select type (f)
               class is (lambda_obj_err)
                  fp = f%eval(xp)
                  lerrp = f%lerr
                  fm = f%eval(xm)
                  lerrm = f%lerr
                  lerr = lerrp .or. lerrm
               class is (lambda_obj)
                  fp = f%eval(xp)
                  fm = f%eval(xm)
                  lerr = .false.
               end select
               grad(i) = (fp - fm) / (2 * dx)
               if (lerr) return
            end do
            return 
         end function gradf


         function minimize1D(f, x0, S, N, eps, lerr) result(astar)
            !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
            !! This program find the minimum of a function of N variables in a single direction
            !! S using in sequence:
            !!    1.  A Bracketing method
            !!    2.  The golden section method
            !!    3.  A quadratic polynomial fit
            !! Inputs
            !!   f%eval(x) : lambda function object containing the objective function as the eval metho
            !!   x0   :  Array of size N of initial x values
            !!   S    :  Array of size N that determines the direction of minimization
            !!   N    :  Number of variables of function f
            !!   eps  :  Accuracy of 1 - dimensional minimization at each step
            !! Output
            !!   lerr : .true. if an error occurred. Otherwise returns .false.
            !! Returns
            !!   astar      :  Final minimum along direction S
            !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
            implicit none
            ! Arguments
            integer(I4B),           intent(in)  :: N
            class(lambda_obj),      intent(inout)  :: f
            real(DP), dimension(:), intent(in)  :: x0, S
            real(DP),               intent(in)  :: eps
            logical,                intent(out) :: lerr
            ! Result
            real(DP)                            :: astar
            ! Internals
            integer(I4B) :: num = 0
            real(DP), parameter :: step = 0.7_DP     !! Bracketing method step size   
            real(DP), parameter :: gam = 1.2_DP      !! Bracketing method expansion parameter   
            real(DP), parameter :: greduce = 0.2_DP  !! Golden section method reduction factor   
            real(DP), parameter :: greduce2 = 0.1_DP ! Secondary golden section method reduction factor   
            real(DP) :: alo, ahi                     !! High and low values for 1 - D minimization routines   
            real(DP), parameter :: a0 = epsilon(1.0_DP)       !! Initial guess of alpha   
         
            alo = a0
            call bracket(f, x0, S, N, gam, step, alo, ahi, lerr)
            if (lerr) then
               !write(*,*) "BFGS bracketing step failed!"
               !write(*,*) "alo: ",alo, "ahi: ", ahi
               return 
            end if
            if (abs(alo - ahi) < eps) then
               astar = alo
               lerr = .false.
               return 
            end if
            call golden(f, x0, S, N, greduce, alo, ahi, lerr)
            if (lerr) then
               !write(*,*) "BFGS golden section step failed!"
               return 
            end if
            if (abs(alo - ahi) < eps) then
               astar = alo
               lerr = .false.
               return 
            end if
            call quadfit(f, x0, S, N, eps, alo, ahi, lerr)
            if (lerr) then
               !write(*,*) "BFGS quadfit failed!"
               return 
            end if
            if (abs(alo - ahi) < eps) then
               astar = alo
               lerr = .false.
               return 
            end if 
            ! Quadratic fit method won't converge, so finish off with another golden section   
            call golden(f, x0, S, N, greduce2, alo, ahi, lerr)
            if (.not. lerr) astar = (alo + ahi) / 2.0_DP
            return 
         end function minimize1D


         function n2one(f, x0, S, N, a, lerr) result(fnew)
            implicit none
            ! Arguments
            integer(I4B),           intent(in) :: N
            class(lambda_obj),      intent(inout) :: f
            real(DP), dimension(:), intent(in) :: x0, S
            real(DP),               intent(in) :: a
            logical,                intent(out) :: lerr

            ! Return
            real(DP) :: fnew
            ! Internals
            real(DP), dimension(N) :: xnew
            integer(I4B) :: i
            
            xnew(:) = x0(:) + a * S(:)
            fnew = f%eval(xnew(:))
            select type(f)
            class is (lambda_obj_err)
               lerr = f%lerr
            class is (lambda_obj)
               lerr = .false.
            end select
            return 
         end function n2one

         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   

         subroutine swiftest_bracket(f, x0, S, N, gam, step, lo, hi, lerr)
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
            !! This subroutine swiftest_brackets the minimum.  It recieves as input:
            !!   f%eval(x) : lambda function object containing the objective function as the eval metho
            !!   x0   :  Array of size N of initial x values
            !!   S    :  Array of size N that determines the direction of minimization
            !!   gam  :  expansion parameter
            !!   step :  step size
            !!   lo   :  initial guess of lo bracket value
            !! The outputs include
            !!   lo   :  lo bracket
            !!   hi   :  hi bracket
            !!   lerr : .true. if an error occurred. Otherwise returns .false.
            !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
            implicit none
            ! Arguments
            integer(I4B),           intent(in)    :: N
            class(lambda_obj),      intent(inout) :: f
            real(DP), dimension(:), intent(in)    :: x0, S
            real(DP),               intent(in)    :: gam, step
            real(DP),               intent(inout) :: lo
            real(DP),               intent(out)   :: hi
            logical,                intent(out)   :: lerr
            ! Internals
            real(DP) :: a0, a1, a2, atmp, da
            real(DP) :: f0, f1, f2
            integer(I4B) :: i, j
            integer(I4B), parameter :: MAXLOOP = 100 ! maximum number of loops before method is determined to have failed   
            real(DP), parameter :: eps = epsilon(lo) ! small number precision to test floating point equality   

            ! set up initial bracket points   
            a0 =  lo
            da = step
            a1 = a0 + da
            a2 = a0 + 2 * da
            f0 = n2one(f, x0, S, N, a0, lerr)
            if (lerr) return
            f1 = n2one(f, x0, S, N, a1, lerr)
            if (lerr) return
            f2 = n2one(f, x0, S, N, a2, lerr)
            if (lerr) return
            ! loop over bracket method until either min is bracketed method fails   
            do i = 1, MAXLOOP 
               if ((f0 > f1) .and. (f1 < f2)) then  ! Minimum was found   
                  lo = a0
                  hi = a2
                  return 
               else if ((f0 >= f1) .and. (f1 > f2)) then ! Function appears to decrease   
                  da = da * gam
                  atmp = a2 + da
                  a0 = a1
                  a1 = a2
                  a2 = atmp
                  f0 = f1
                  f1 = f2
                  f2 = n2one(f, x0, S, N, a2, lerr)
               else if ((f0 < f1) .and. (f1 <= f2)) then ! Function appears to increase   
                  da = da * gam
                  atmp = a0 - da
                  a2 = a1
                  a1 = a0
                  a0 = atmp
                  f2 = f1
                  f0 = n2one(f, x0, S, N, a0, lerr)
               else if ((f0 < f1) .and. (f1 > f2)) then ! We are at a peak. Pick the direction that descends the fastest
                  da = da * gam
                  if (f2 > f0) then ! LHS is lower than RHS
                     atmp = a2 + da
                     a0 = a1
                     a1 = a2
                     a2 = atmp
                     f0 = f1
                     f1 = f2
                     f2 = n2one(f, x0, S, N, a2, lerr)
                  else ! RHS is lower than LHS
                     atmp = a0 - da
                     a2 = a1
                     a1 = a0
                     a0 = atmp
                     f2 = f1
                     f1 = f2
                     f0 = n2one(f, x0, S, N, a0, lerr)
                  end if
               else if ((f0 > f1) .and. (abs(f2 - f1) <= eps)) then ! Decrasging but RHS equal   
                  da = da * gam
                  atmp = a2 + da
                  a2 = atmp
                  f2 = n2one(f, x0, S, N, a2, lerr)
               else if ((abs(f0 - f1) < eps) .and. (f1 < f2)) then ! Increasing but LHS equal   
                  da = da * gam
                  atmp = a0 - da
                  a0 = atmp
                  f0 = n2one(f, x0, S, N, a0, lerr)
               else  ! all values equal. Expand in either direction and try again
                  a0 = a0 - da
                  a2 = a2 + da
                  f0 = n2one(f, x0, S, N, a0, lerr)
                  if (lerr) exit ! An error occurred while evaluating the function
                  f2 = n2one(f, x0, S, N, a2, lerr)
               end if
               if (lerr) exit ! An error occurred while evaluating the function
            end do
            lerr = .true.
            return ! no minimum found   
         end subroutine swiftest_bracket

         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   

         subroutine swiftest_golden(f, x0, S, N, eps, lo, hi, lerr) 
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
            !! This function uses the golden section method to reduce the starting interval lo, hi by some amount sigma.  
            !! It recieves as input:
            !!   f%eval(x) : lambda function object containing the objective function as the eval metho
            !!   x0   :  Array of size N of initial x values
            !!   S    :  Array of size N that determines the direction of minimization
            !!   gam  :  expansion parameter
            !!   eps  :  reduction interval in range (0 < sigma < 1) such that:
            !!             hi(new) - lo(new) = eps * (hi(old) - lo(old))
            !!   lo   :  initial guess of lo bracket value
            !! The outputs include
            !!   lo   :  lo bracket
            !!   hi   :  hi bracket
            !!   lerr : .true. if an error occurred. Otherwise returns .false.
            !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
            implicit none
            ! Arguments
            integer(I4B),           intent(in)    :: N
            class(lambda_obj),      intent(inout) :: f
            real(DP), dimension(:), intent(in)    :: x0, S
            real(DP),               intent(in)    :: eps
            real(DP),               intent(inout) :: lo
            real(DP),               intent(out)   :: hi
            logical,                intent(out)   :: lerr
            ! Internals 
            real(DP), parameter :: tau = 0.5_DP * (sqrt(5.0_DP) - 1.0_DP)  ! Golden section constant   
            integer(I4B), parameter :: MAXLOOP = 40 ! maximum number of loops before method is determined to have failed (unlikely, but could occur if no minimum exists between lo and hi)   
            real(DP) :: i0 ! Initial interval value   
            real(DP) :: a1, a2
            real(DP) :: f1, f2
            integer(I4B) :: i, j

            i0 =  hi - lo
            a1 =  hi - tau * i0
            a2 =  lo + tau * i0
            f1 = n2one(f, x0, S, N, a1, lerr)
            if (lerr) return
            f2 = n2one(f, x0, S, N, a2, lerr)
            if (lerr) return
            do i = 1, MAXLOOP 
               if (abs((hi - lo) / i0) <= eps) return ! interval reduced to input amount   
               if (f2 > f1) then
                  hi = a2
                  a2 = a1
                  f2 = f1
                  a1 = hi - tau * (hi - lo)
                  f1 = n2one(f, x0, S, N, a1, lerr)
               else 
                  lo = a1
                  a1 = a2
                  f2 = f1
                  a2 = hi - (1.0_DP - tau) * (hi - lo)
                  f2 = n2one(f, x0, S, N, a2, lerr)
               end if
               if (lerr) exit
            end do
            lerr = .true.
            return ! search took too many iterations - no minimum found   
         end subroutine swiftest_golden

         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   

         subroutine swiftest_quadfit(f, x0, S, N, eps, lo, hi, lerr) 
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
            !! This function uses a quadratic polynomial fit to locate the minimum of a function
            !! to some accuracy eps.  It recieves as input:
            !!   f%eval(x) : lambda function object containing the objective function as the eval metho
            !!   lo    :  low bracket value
            !!   hi    :  high bracket value
            !!   eps   :  desired accuracy of final minimum location
            !! The outputs include
            !!   lo   :  final minimum location
            !!   hi   :  final minimum location
            !! Notes: Uses the ieee_exceptions intrinsic module to allow for graceful failure due to floating point exceptions, which won't terminate the run.
            !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
            implicit none
            ! Arguments
            integer(I4B),           intent(in)    :: N
            class(lambda_obj),      intent(inout) :: f
            real(DP), dimension(:), intent(in)    :: x0, S
            real(DP),               intent(in)    :: eps
            real(DP),               intent(inout) :: lo
            real(DP),               intent(out)   :: hi
            logical,                intent(out)   :: lerr
            ! Internals 
            integer(I4B), parameter :: MAXLOOP = 20 ! maximum number of loops before method is determined to have failed.   
            real(DP) :: a1, a2, a3, astar   ! three points for the polynomial fit and polynomial minimum   
            real(DP) :: f1, f2, f3, fstar   ! three function values for the polynomial and polynomial minimum   
            real(DP), dimension(3) :: row_1, row_2, row_3, rhs, soln        ! matrix for 3 equation solver (gaussian elimination)   
            real(DP), dimension(3,3) :: lhs
            real(DP) :: d1, d2, d3, aold, denom, errval
            integer(I4B) :: i

            lerr = .false.
            ! Get initial a1, a2, a3 values   
            a1 =  lo
            a2 =  lo + 0.5_DP * (hi - lo)
            a3 =  hi
            aold = a1
            astar = a2
            f1 = n2one(f, x0, S, N, a1, lerr)
            if (lerr) return
            f2 = n2one(f, x0, S, N, a2, lerr)
            if (lerr) return
            f3 = n2one(f, x0, S, N, a3, lerr)
            if (lerr) return
            do i = 1, MAXLOOP 
               ! check to see if convergence is reached and exit   
               errval = abs((astar - aold) / astar)
               call ieee_get_flag(ieee_usual, fpe_flag)
               if (any(fpe_flag)) then
                  !write(*,*) 'quadfit fpe'
                  !write(*,*) 'aold : ',aold
                  !write(*,*) 'astar: ',astar
                  lerr = .true.
                  exit
               end if
               if (errval < eps) then
                  lo = astar
                  hi = astar
                  exit
               end if
               ! Set up system for gaussian elimination equation solver   
               row_1 = [1.0_DP, a1, a1**2]
               row_2 = [1.0_DP, a2, a2**2]
               row_3 = [1.0_DP, a3, a3**2]
               rhs = [f1, f2, f3]
               lhs(1, :) = row_1
               lhs(2, :) = row_2
               lhs(3, :) = row_3
               ! Solve system of equations   
               soln(:) = util_solve_linear_system(lhs, rhs, 3, lerr)
               call ieee_set_flag(ieee_all, .false.) ! Set all flags back to quiet
               call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
               if (lerr) then
                  !write(*,*) 'quadfit fpe:'
                  !write(*,*) 'util_solve_linear_system failed'
                  exit
               end if
               aold = astar
               if (soln(2) == soln(3)) then ! Handles the case where they are both 0. 0/0 is an unhandled exception
                  astar = -0.5_DP
               else
                  astar =  -soln(2) / (2 * soln(3))
               end if
               call ieee_get_flag(ieee_usual, fpe_flag)
               if (any(fpe_flag)) then
                  !write(*,*) 'quadfit fpe'
                  !write(*,*) 'soln(2:3): ',soln(2:3)
                  !write(*,*) 'a1, a2, a3'
                  !write(*,*) a1, a2, a3
                  !write(*,*) 'f1, f2, f3'
                  !write(*,*) f1, f2, f3
                  lerr = .true.
                  exit
               end if
               fstar = n2one(f, x0, S, N, astar, lerr)
               if (lerr) exit
               ! keep the three closest a values to astar and discard the fourth  
               d1 = abs(a1 - astar)
               d2 = abs(a2 - astar)
               d3 = abs(a3 - astar)

               if (d1 > d2) then
                  if (d1 > d3) then
                     f1 = fstar
                     a1 = astar
                  else if (d3 > d2) then
                     f3 = fstar
                     a3 = astar
                  end if
               else 
                  if (d2 > d3) then
                     f2 = fstar
                     a2 = astar
                  else if (d3 > d1) then
                     f3 = fstar
                     a3 = astar
                  end if
               end if
            end do
            if (lerr) return
            lo = a1
            hi = a3
            return 
         end subroutine swiftest_quadfit

   end subroutine swiftest_util_minimize_bfgs

   subroutine swiftest_util_peri(n,m, r, v, atp, q, isperi)
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

      do concurrent(i = 1:n)
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


   module subroutine swiftest_util_peri_body(self, system, param)
      !! author: David A. Minton
      !!
      !! Determine system pericenter passages for bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_peri.f90
      !! Adapted from Hal Levison's Swift routine util_mass_peri.f
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i

      select type(self)
      class is (swiftest_pl)
         if (self%lfirst) self%isperi(:) = 0
      end select

      if (param%qmin_coord == "HELIO") then
         call swiftest_util_peri(self%nbody, self%mu, self%rh, self%vh, self%atp, self%peri, self%isperi)
      else 
         call swiftest_util_peri(self%nbody, [(system%Gmtot,i=1,self%nbody)], self%rb, self%vb, self%atp, self%peri, self%isperi)
      end if

      return
   end subroutine swiftest_util_peri_body


   module subroutine swiftest_util_rearray_pl(self, system, param)
      !! Author: the Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Clean up the massive body structures to remove discarded bodies and add new bodies
      use symba
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      class(swiftest_pl), allocatable :: tmp !! The discarded body list.
      integer(I4B) :: i, k, npl, nadd, nencmin, nenc_old, idnew1, idnew2, idold1, idold2
      logical, dimension(:), allocatable :: lmask, ldump_mask
      class(collision_list_plpl), allocatable :: plplenc_old
      logical :: lencounter
      integer(I4B), dimension(:), allocatable :: levelg_orig_pl, levelm_orig_pl, levelg_orig_tp, levelm_orig_tp
      integer(I4B), dimension(:), allocatable :: nplenc_orig_pl, nplenc_orig_tp, ntpenc_orig_pl

      associate(pl => self, pl_adds => system%pl_adds)

         npl = pl%nbody
         nadd = pl_adds%nbody
         if (npl == 0) return
         ! Deallocate any temporary variables
         if (allocated(pl%rbeg)) deallocate(pl%rbeg)
         if (allocated(pl%rend)) deallocate(pl%rend)

         ! Remove the discards and destroy the list, as the system already tracks pl_discards elsewhere
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

         ! Store the original plplenc list so we don't remove any of the original encounters
         nenc_old = system%plpl_encounter%nenc
         if (nenc_old > 0) then 
            allocate(plplenc_old, source=system%plpl_encounter)
            call plplenc_old%copy(system%plpl_encounter)
         end if

         ! Add in any new bodies
         if (nadd > 0) then
            ! Append the adds to the main pl object
            call pl%append(pl_adds, lsource_mask=[(.true., i=1, nadd)])

            allocate(ldump_mask(npl+nadd))  ! This mask is used only to append the original Fortran binary particle.dat file with new bodies. This is ignored for NetCDF output
            ldump_mask(1:npl) = .false.
            ldump_mask(npl+1:npl+nadd) = pl%status(npl+1:npl+nadd) == NEW_PARTICLE
            npl = pl%nbody
         else
            allocate(ldump_mask(npl))
            ldump_mask(:) = .false.
         end if

         ! Reset all of the status flags for this body
         pl%status(1:npl) = ACTIVE
         do i = 1, npl
            call pl%info(i)%set_value(status="ACTIVE")
         end do
         pl%ldiscard(1:npl) = .false.
         pl%lcollision(1:npl) = .false.
         pl%lmask(1:npl) = .true.

         pl%lmtiny(1:npl) = pl%Gmass(1:npl) < param%GMTINY
         where(pl%lmtiny(1:npl))
            pl%info(1:npl)%particle_type = PL_TINY_TYPE_NAME 
         elsewhere
            pl%info(1:npl)%particle_type = PL_TYPE_NAME 
         end where

         call pl%write_info(param%system_history%nc, param)
         deallocate(ldump_mask)

         ! Reindex the new list of bodies 
         call pl%sort("mass", ascending=.false.)
         call pl%flatten(param)

         ! Reset the kinship trackers
         call pl%reset_kinship([(i, i=1, npl)])

         ! Re-build the zero-level encounter list, being sure to save the original level information for all bodies
         allocate(nplenc_orig_pl, source=pl%nplenc)
         select type(pl)
         class is (symba_pl)
            allocate(levelg_orig_pl, source=pl%levelg)
            allocate(levelm_orig_pl, source=pl%levelm)
            call move_alloc(levelg_orig_pl, pl%levelg)
            call move_alloc(levelm_orig_pl, pl%levelm)
         end select
         lencounter = pl%encounter_check(param, system, param%dt, system%irec) 
         if (system%tp%nbody > 0) then
            allocate(ntpenc_orig_pl, source=pl%ntpenc)
            allocate(nplenc_orig_tp, source=tp%nplenc)
            select type(tp => system%tp)
            class is (symba_tp)
               allocate(levelg_orig_tp, source=tp%levelg)
               allocate(levelm_orig_tp, source=tp%levelm)
               lencounter = tp%encounter_check(param, system, param%dt, system%irec)
               call move_alloc(levelg_orig_tp, tp%levelg)
               call move_alloc(levelm_orig_tp, tp%levelm)
               call move_alloc(nplenc_orig_tp, tp%nplenc)
               call move_alloc(ntpenc_orig_pl, pl%ntpenc)
            end select
         end if

         call move_alloc(nplenc_orig_pl, pl%nplenc)

         ! Re-index the encounter list as the index values may have changed
         if (nenc_old > 0) then
            nencmin = min(system%plpl_encounter%nenc, plplenc_old%nenc) 
            system%plpl_encounter%nenc = nencmin
            do k = 1, nencmin
               idnew1 = system%plpl_encounter%id1(k)
               idnew2 = system%plpl_encounter%id2(k)
               idold1 = plplenc_old%id1(k)
               idold2 = plplenc_old%id2(k)
               if ((idnew1 == idold1) .and. (idnew2 == idold2)) then
                  ! This is an encounter we already know about, so save the old information
                  system%plpl_encounter%lvdotr(k) = plplenc_old%lvdotr(k) 
                  system%plpl_encounter%lclosest(k) = plplenc_old%lclosest(k) 
                  system%plpl_encounter%status(k) = plplenc_old%status(k) 
                  system%plpl_encounter%r1(:,k) = plplenc_old%r1(:,k)
                  system%plpl_encounter%r2(:,k) = plplenc_old%r2(:,k)
                  system%plpl_encounter%v1(:,k) = plplenc_old%v1(:,k)
                  system%plpl_encounter%v2(:,k) = plplenc_old%v2(:,k)
                  system%plpl_encounter%tcollision(k) = plplenc_old%tcollision(k)
                  system%plpl_encounter%level(k) = plplenc_old%level(k)
               else if (((idnew1 == idold2) .and. (idnew2 == idold1))) then
                  ! This is an encounter we already know about, but with the order reversed, so save the old information
                  system%plpl_encounter%lvdotr(k) = plplenc_old%lvdotr(k) 
                  system%plpl_encounter%lclosest(k) = plplenc_old%lclosest(k) 
                  system%plpl_encounter%status(k) = plplenc_old%status(k) 
                  system%plpl_encounter%r1(:,k) = plplenc_old%r2(:,k)
                  system%plpl_encounter%r2(:,k) = plplenc_old%r1(:,k)
                  system%plpl_encounter%v1(:,k) = plplenc_old%v2(:,k)
                  system%plpl_encounter%v2(:,k) = plplenc_old%v1(:,k)
                  system%plpl_encounter%tcollision(k) = plplenc_old%tcollision(k)
                  system%plpl_encounter%level(k) = plplenc_old%level(k)
               end if
               system%plpl_encounter%index1(k) = findloc(pl%id(1:npl), system%plpl_encounter%id1(k), dim=1)
               system%plpl_encounter%index2(k) = findloc(pl%id(1:npl), system%plpl_encounter%id2(k), dim=1)
            end do
            if (allocated(lmask)) deallocate(lmask)
            allocate(lmask(nencmin))
            nenc_old = nencmin
            if (any(system%plpl_encounter%index1(1:nencmin) == 0) .or. any(system%plpl_encounter%index2(1:nencmin) == 0)) then
               lmask(:) = system%plpl_encounter%index1(1:nencmin) /= 0 .and. system%plpl_encounter%index2(1:nencmin) /= 0
            else
               return
            end if
            nencmin = count(lmask(:))
            system%plpl_encounter%nenc = nencmin
            if (nencmin > 0) then
               system%plpl_encounter%index1(1:nencmin) = pack(system%plpl_encounter%index1(1:nenc_old), lmask(1:nenc_old))
               system%plpl_encounter%index2(1:nencmin) = pack(system%plpl_encounter%index2(1:nenc_old), lmask(1:nenc_old))
               system%plpl_encounter%id1(1:nencmin) = pack(system%plpl_encounter%id1(1:nenc_old), lmask(1:nenc_old))
               system%plpl_encounter%id2(1:nencmin) = pack(system%plpl_encounter%id2(1:nenc_old), lmask(1:nenc_old))
               system%plpl_encounter%lvdotr(1:nencmin) = pack(system%plpl_encounter%lvdotr(1:nenc_old), lmask(1:nenc_old))
               system%plpl_encounter%lclosest(1:nencmin) = pack(system%plpl_encounter%lclosest(1:nenc_old), lmask(1:nenc_old))
               system%plpl_encounter%status(1:nencmin) = pack(system%plpl_encounter%status(1:nenc_old), lmask(1:nenc_old))
               system%plpl_encounter%tcollision(1:nencmin) = pack(system%plpl_encounter%tcollision(1:nenc_old), lmask(1:nenc_old))
               system%plpl_encounter%level(1:nencmin) = pack(system%plpl_encounter%level(1:nenc_old), lmask(1:nenc_old))
               do i = 1, NDIM
                  system%plpl_encounter%r1(i, 1:nencmin) = pack(system%plpl_encounter%r1(i, 1:nenc_old), lmask(1:nenc_old))
                  system%plpl_encounter%r2(i, 1:nencmin) = pack(system%plpl_encounter%r2(i, 1:nenc_old), lmask(1:nenc_old))
                  system%plpl_encounter%v1(i, 1:nencmin) = pack(system%plpl_encounter%v1(i, 1:nenc_old), lmask(1:nenc_old))
                  system%plpl_encounter%v2(i, 1:nencmin) = pack(system%plpl_encounter%v2(i, 1:nenc_old), lmask(1:nenc_old))
               end do
            end if
         end if
      end associate

      return
   end subroutine swiftest_util_rearray_pl


   module subroutine swiftest_util_rescale_system(self, param, mscale, dscale, tscale)
      !! author: David A. Minton
      !!
      !! Rescales an nbody system to a new set of units. Inputs are the multipliers on the mass (mscale), distance (dscale), and time units (tscale). 
      !! Rescales all united quantities in the system, as well as the mass conversion factors, gravitational constant, and Einstein's constant in the parameter object.
      implicit none
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters. Returns with new values of the scale vactors and GU
      real(DP),                     intent(in)    :: mscale, dscale, tscale !! Scale factors for mass, distance, and time units, respectively. 
      ! Internals
      real(DP) :: vscale

      param%MU2KG = param%MU2KG * mscale
      param%DU2M = param%DU2M * dscale
      param%TU2S = param%TU2S * tscale

      ! Calculate the G for the system units
      param%GU = GC / (param%DU2M**3 / (param%MU2KG * param%TU2S**2))

      if (param%lgr) then
         ! Calculate the inverse speed of light in the system units
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


   module subroutine swiftest_util_resize_arr_char_string(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. nnew = 0 will deallocate.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                     intent(in)    :: nnew !! New size
      ! Internals
      character(len=STRMAX), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
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
      if (nold > 0) then
         if (nnew > nold) then
            tmp(1:nold) = arr(1:nold)
            tmp(nold+1:nnew) = ""
         else
            tmp(1:nnew) = arr(1:nnew)
         end if
      else
         tmp(1:nnew) = ""
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_resize_arr_char_string


   module subroutine swiftest_util_resize_arr_DP(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of double precision type. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                        intent(in)    :: nnew !! New size
      ! Internals
      real(DP), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size
      real(DP), parameter :: init_val = 0.0_DP

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
      if (nold > 0) then
         if (nnew > nold) then
            tmp(1:nold) = arr(1:nold)
            tmp(nold+1:nnew) = init_val
         else
            tmp(1:nnew) = arr(1:nnew)
         end if
      else
         tmp(1:nnew) = init_val
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_resize_arr_DP


   module subroutine swiftest_util_resize_arr_DPvec(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of double precision vectors of size (NDIM, n). Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                          intent(in)    :: nnew !! New size
      ! Internals
      real(DP), dimension(:,:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size
      real(DP), dimension(NDIM), parameter :: init_val = 0.0_DP
      integer(I4B) :: i

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr, dim=2)
      else
         nold = 0
      end if

      if (nnew == nold) return
      
      allocate(tmp(NDIM, nnew))
      if (nold > 0) then
         if (nnew > nold) then
            tmp(:,1:nold) = arr(:,1:nold)
            do i = nold+1, nnew
               tmp(:,i) = init_val(:)
            end do
         else
            tmp(:,1:nnew) = arr(:,1:nnew)
         end if
      else
         do i = 1, nnew
            tmp(:, i) = init_val(:)
         end do
      end if
      call move_alloc(tmp, arr)

      return

      return
   end subroutine swiftest_util_resize_arr_DPvec


   module subroutine swiftest_util_resize_arr_I4B(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of integer type. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                            intent(in)    :: nnew !! New size
      ! Internals
      integer(I4B), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size
      integer(I4B), parameter :: init_val = -1

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
      if (nold > 0) then
         if (nnew > nold) then
            tmp(1:nold) = arr(1:nold)
            tmp(nold+1:nnew) = init_val
         else
            tmp(1:nnew) = arr(1:nnew)
         end if
      else
         tmp(1:nnew) = init_val
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_resize_arr_I4B


   module subroutine swiftest_util_resize_arr_info(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                         intent(in)    :: nnew !! New size
      ! Internals
      type(swiftest_particle_info), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
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
         call util_copy_particle_info_arr(arr(1:nold), tmp(1:nold))
      else
         call util_copy_particle_info_arr(arr(1:nnew), tmp(1:nnew))
      end if

      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_resize_arr_info


   module subroutine swiftest_util_resize_arr_logical(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of logical type. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                       intent(in)    :: nnew !! New size
      ! Internals
      logical, dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size
      logical, parameter :: init_val = .false.

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
      if (nold > 0) then
         if (nnew > nold) then
            tmp(1:nold) = arr(1:nold)
            tmp(nold+1:nnew) = init_val
         else
            tmp(1:nnew) = arr(1:nnew)
         end if
      else
         tmp(1:nnew) = init_val
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_resize_arr_logical


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

      call util_resize_body(self, nnew)

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

      call util_resize_body(self, nnew)

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
      !! Sets the value of msys and the vector mass quantities based on the total mass of the system
      implicit none
      ! Arguments
      class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest nobdy system object

      self%Gmtot = self%cb%Gmass + sum(self%pl%Gmass(1:self%pl%nbody), self%pl%status(1:self%pl%nbody) /= INACTIVE)

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


   module subroutine swiftest_util_set_particle_info(self, name, particle_type, status, origin_type, origin_time, collision_id, origin_rh,&
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

      call self%xv2el(cb) 
      self%rhill(1:self%nbody) = self%a(1:self%nbody) * (self%Gmass(1:self%nbody) / cb%Gmass / 3)**THIRD 

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


   module subroutine swiftest_util_snapshot_system(self, param, system, t, arg)
      !! author: David A. Minton
      !!
      !! Takes a snapshot of the system for later file storage
      implicit none
      ! Arguments
      class(swiftest_storage(*)),   intent(inout) :: self   !! Swiftest storage object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object to store
      real(DP),                     intent(in), optional :: t      !! Time of snapshot if different from system time
      character(*),                 intent(in), optional :: arg    !! Optional argument (needed for extended storage type used in collision snapshots)

      self%iframe = self%iframe + 1
      self%nt = self%iframe
      self%frame(self%iframe) = system ! Store a snapshot of the system for posterity
      self%nid = self%nid + 1 ! Central body
      if (allocated(system%pl)) self%nid = self%nid + system%pl%nbody
      if (allocated(system%tp)) self%nid = self%nid + system%tp%nbody

      return
   end subroutine swiftest_util_snapshot_system


   module function swiftest_util_solve_linear_system_d(A,b,n,lerr) result(x)
      !! Author: David A. Minton
      !!
      !! Solves the linear equation of the form A*x = b for x. 
      !!   A is an (n,n) arrays
      !!   x and b are (n) arrays
      !! Uses Gaussian elimination, so will have issues if system is ill-conditioned.
      !! Uses quad precision intermidiate values, so works best on small arrays.
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      integer(I4B),             intent(in)  :: n
      real(DP), dimension(:,:), intent(in)  :: A
      real(DP), dimension(:),   intent(in)  :: b
      logical,                  intent(out) :: lerr
      ! Result
      real(DP), dimension(n)                :: x
      ! Internals
      real(QP), dimension(:), allocatable :: qx
      type(ieee_status_type) :: original_fpe_status
      logical, dimension(:), allocatable :: fpe_flag 

      call ieee_get_status(original_fpe_status) ! Save the original floating point exception status
      call ieee_set_flag(ieee_all, .false.) ! Set all flags to quiet
      allocate(fpe_flag(size(ieee_usual)))

      qx = solve_wbs(ge_wpp(real(A, kind=QP), real(b, kind=QP)))

      call ieee_get_flag(ieee_usual, fpe_flag)
      lerr = any(fpe_flag) 
      if (lerr .or. (any(abs(qx) > huge(x))) .or. (any(abs(qx) < tiny(x)))) then
         x = 0.0_DP
      else
         x = real(qx, kind=DP)
      end if
      call ieee_set_status(original_fpe_status)

      return
   end function util_solve_linear_system_d


   module function swiftest_util_solve_linear_system_q(A,b,n,lerr) result(x)
      !! Author: David A. Minton
      !!
      !! Solves the linear equation of the form A*x = b for x. 
      !!   A is an (n,n) arrays
      !!   x and b are (n) arrays
      !! Uses Gaussian elimination, so will have issues if system is ill-conditioned.
      !! Uses quad precision intermidiate values, so works best on small arrays.
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      integer(I4B),             intent(in) :: n
      real(QP), dimension(:,:), intent(in) :: A
      real(QP), dimension(:),   intent(in) :: b
      logical,                  intent(out) :: lerr
      ! Result
      real(QP), dimension(n)  :: x
      ! Internals
      type(ieee_status_type) :: original_fpe_status
      logical, dimension(:), allocatable :: fpe_flag 

      call ieee_get_status(original_fpe_status) ! Save the original floating point exception status
      call ieee_set_flag(ieee_all, .false.) ! Set all flags to quiet
      allocate(fpe_flag(size(ieee_usual)))

      x = solve_wbs(ge_wpp(A, b))

      call ieee_get_flag(ieee_usual, fpe_flag)
      lerr = any(fpe_flag) 
      if (lerr) x = 0.0_DP
      call ieee_set_status(original_fpe_status) 

      return
   end function util_solve_linear_system_q

   function solve_wbs(u) result(x) ! solve with backward substitution
      !! Based on code available on Rosetta Code: https://rosettacode.org/wiki/Gaussian_elimination#Fortran
      use, intrinsic :: ieee_exceptions
      use swiftest
      implicit none
      ! Arguments
      real(QP), intent(in), dimension(:,:), allocatable  :: u
      ! Result
      real(QP), dimension(:), allocatable :: x
      ! Internals
      integer(I4B)             :: i,n

      n = size(u, 1)
      if (allocated(x)) deallocate(x)
      if (.not.allocated(x)) allocate(x(n))
      if (any(abs(u) < tiny(1._DP)) .or. any(abs(u) > huge(1._DP))) then 
         x(:) = 0._DP
         return
      end if
      call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      do i = n, 1, -1 
         x(i) = (u(i, n + 1) - sum(u(i, i + 1:n) * x(i + 1:n))) / u(i, i)
      end do   
      return
   end function solve_wbs


   function ge_wpp(A, b) result(u) ! gaussian eliminate with partial pivoting
      !! Solve  Ax=b  using Gaussian elimination then backwards substitution.
      !!   A being an n by n matrix.
      !!   x and b are n by 1 vectors. 
      !! Based on code available on Rosetta Code: https://rosettacode.org/wiki/Gaussian_elimination#Fortran
      use, intrinsic :: ieee_exceptions
      use swiftest
      implicit none
      ! Arguments
      real(QP), dimension(:,:), intent(in) :: A
      real(QP), dimension(:),   intent(in) :: b
      ! Result
      real(QP), dimension(:,:), allocatable :: u
      ! Internals
      integer(I4B) :: i,j,n,p
      real(QP)     ::  upi

      n = size(a, 1)
      allocate(u(n, (n + 1)))
      u = reshape([A, b], [n, n + 1])
      call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      do j = 1, n
         p = maxloc(abs(u(j:n, j)), 1) + j - 1 ! maxloc returns indices between (1, n - j + 1)
         if (p /= j) u([p, j], j) = u([j, p], j)
         u(j + 1:, j) = u(j + 1:, j) / u(j, j)
         do i = j + 1, n + 1
            upi = u(p, i)
            if (p /= j) u([p, j], i) = u([j, p], i)
            u(j + 1:n, i) = u(j + 1:n, i) - upi * u(j + 1:n, j)
         end do
      end do
      return
   end function ge_wpp


   module function swiftest_util_solve_rkf45(f, y0in, t1, dt0, tol) result(y1)
      !! author: David A. Minton
      !!
      !! Implements the 4th order Runge-Kutta-Fehlberg ODE solver for initial value problems of the form f=dy/dt, y0 = y(t=0), solving for y1 = y(t=t1). Uses a 5th order adaptive step size control.
      !! Uses a lambda function object as defined in the lambda_function module
      implicit none
      ! Arguments   
      class(lambda_obj),      intent(inout) :: f    !! lambda function object that has been initialized to be a function of derivatives. The object will return with components lastarg and lasteval set
      real(DP), dimension(:), intent(in)    :: y0in !! Initial value at t=0
      real(DP),               intent(in)    :: t1   !! Final time
      real(DP),               intent(in)    :: dt0  !! Initial step size guess
      real(DP),               intent(in)    :: tol  !! Tolerance on solution
      ! Result
      real(DP), dimension(:), allocatable   :: y1  !! Final result
      ! Internals
      integer(I4B),                          parameter :: MAXREDUX = 1000 !! Maximum number of times step size can be reduced
      real(DP),                              parameter :: DTFAC = 0.95_DP !! Step size reduction safety factor (Value just under 1.0 to prevent adaptive step size control from discarding steps too aggressively)
      integer(I4B),                          parameter :: RKS = 6         !! Number of RK stages
      real(DP),     dimension(RKS, RKS - 1), parameter :: rkf45_btab = reshape( & !! Butcher tableau for Runge-Kutta-Fehlberg method
         (/       1./4.,       1./4.,          0.,            0.,           0.,           0.,&
                  3./8.,      3./32.,      9./32.,            0.,           0.,           0.,&
                12./13., 1932./2197., -7200./2197.,  7296./2197.,           0.,           0.,&
                     1.,   439./216.,          -8.,   3680./513.,   -845./4104.,          0.,&
                  1./2.,     -8./27.,           2., -3544./2565.,   1859./4104.,    -11./40./), shape(rkf45_btab))
      real(DP), dimension(RKS),  parameter   :: rkf4_coeff =  (/ 25./216., 0., 1408./2565. ,  2197./4104. , -1./5.,      0. /)
      real(DP), dimension(RKS),  parameter   :: rkf5_coeff =  (/ 16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55. /)
      real(DP), dimension(:, :), allocatable :: k                !! Runge-Kutta coefficient vector
      real(DP), dimension(:),   allocatable  :: ynorm            !! Normalized y value used for adaptive step size control
      real(DP), dimension(:),   allocatable  :: y0           !! Value of y at the beginning of each substep
      integer(I4B)                           :: Nvar             !! Number of variables in problem
      integer(I4B)                           :: rkn              !! Runge-Kutta loop index
      real(DP)                               :: t, x1, dt, trem      !! Current time, step size and total time remaining
      real(DP)                               :: s, yerr, yscale  !!  Step size reduction factor, error in dependent variable, and error scale factor
      integer(I4B)                           :: i

      allocate(y0, source=y0in)
      allocate(y1, mold=y0)
      allocate(ynorm, mold=y0)
      Nvar = size(y0)
      allocate(k(Nvar, RKS))

      dt = dt0

      trem = t1
      t = 0._DP
      do
         yscale = norm2(y0(:))
         do i = 1, MAXREDUX
            select type(f)
            class is (lambda_obj_tvar)
               do rkn = 1, RKS
                  y1(:) = y0(:) + matmul(k(:, 1:rkn - 1), rkf45_btab(2:rkn, rkn - 1))
                  if (rkn == 1) then
                     x1 = t
                  else
                     x1 = t + rkf45_btab(1,rkn-1)
                  end if
                  k(:, rkn) = dt * f%evalt(y1(:), t)
               end do
            class is (lambda_obj)
               do rkn = 1, RKS
                  y1(:) = y0(:) + matmul(k(:, 1:rkn - 1), rkf45_btab(2:rkn, rkn - 1))
                  k(:, rkn) = dt * f%eval(y1(:))
               end do
            end select
            ! Now determine if the step size needs adjusting
            ynorm(:) = matmul(k(:,:), (rkf5_coeff(:) - rkf4_coeff(:))) / yscale
            yerr = norm2(ynorm(:)) 
            s = (tol / (2 * yerr))**(0.25_DP)
            dt = min(s * DTFAC * dt, trem) ! Alter step size either up or down, but never bigger than the remaining time
            if (s >= 1.0_DP) exit ! Good step!
            if (i == MAXREDUX) then
               write(*,*) "Something has gone wrong in util_solve_rkf45!! Step size reduction has gone too far this time!"
               call util_exit(FAILURE)
            end if
         end do
      
         ! Compute new value then step ahead in time
         y1(:) = y0(:) + matmul(k(:, :), rkf4_coeff(:))
         trem = trem - dt
         t = t + dt
         if (trem <= 0._DP) exit
         y0(:) = y1(:)
      end do

      return
   end function util_solve_rkf45



   module subroutine swiftest_util_sort_body(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a Swiftest body structure in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self      !! Swiftest body object
      character(*),         intent(in)    :: sortby    !! Sorting attribute
      logical,              intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
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
         case("lfirst", "nbody", "ldiscard", "rh", "vh", "rb", "vb", "ah", "aobl", "atide", "agr")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not found!'
            return
         end select

         call body%rearrange(ind)

      end associate

      return
   end subroutine swiftest_util_sort_body


   pure module subroutine swiftest_util_sort_dp(arr)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array in place into ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(inout) :: arr

      call qsort_DP(arr)

      return
   end subroutine swiftest_util_sort_dp


   pure module subroutine swiftest_util_sort_index_dp(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array by index in ascending numerical order using quick sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here).
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine swiftest_allocates it.
      !!
      implicit none
      ! Arguments
      real(DP),     dimension(:),              intent(in)    :: arr
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      integer(I4B) :: n, i
      real(DP), dimension(:), allocatable :: tmparr

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1, n)]
      end if
      allocate(tmparr, mold=arr)
      tmparr(:) = arr(ind(:))
      call qsort_DP(tmparr, ind)
   
      return
   end subroutine swiftest_util_sort_index_dp


   recursive pure subroutine swiftest_qsort_DP(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array by index in ascending numerical order using quicksort sort.
      !!
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(inout)           :: arr
      integer(I4B),dimension(:),intent(out), optional :: ind
      !! Internals
      integer :: iq

      if (size(arr) > 1) then
         if (present(ind)) then
            call partition_DP(arr, iq, ind)
            call qsort_DP(arr(:iq-1),ind(:iq-1))
            call qsort_DP(arr(iq:),  ind(iq:))
         else
            call partition_DP(arr, iq)
            call qsort_DP(arr(:iq-1))
            call qsort_DP(arr(iq:))
         end if
      end if

      return
   end subroutine swiftest_qsort_DP

 
   pure subroutine swiftest_partition_DP(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on DP type
      !!
      implicit none
      ! Arguments
      real(DP),     intent(inout), dimension(:)           :: arr
      integer(I4B), intent(inout), dimension(:), optional :: ind
      integer(I4B), intent(out)                           :: marker
      ! Internals
      integer(I4B) :: i, j, itmp, narr, ipiv
      real(DP) :: temp
      real(DP) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2
      x = arr(ipiv)
      i = 0
      j = narr + 1
   
      do
         j = j - 1
         do
            if (arr(j) <= x) exit
            j = j - 1
         end do
         i = i + 1
         do
            if (arr(i) >= x) exit
            i = i + 1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine swiftest_partition_DP
 

   pure module subroutine swiftest_util_sort_i4b(arr)
      !! author: David A. Minton
      !!
      !! Sort input integer array in place into ascending numerical order using quick sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here)
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:), intent(inout) :: arr

      call qsort_I4B(arr)

      return
   end subroutine swiftest_util_sort_i4b


   pure module subroutine swiftest_util_sort_index_I4B(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input integer array by index in ascending numerical order using quicksort.
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine swiftest_allocates it.
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:),              intent(in)  :: arr
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      integer(I4B) :: n, i
      integer(I4B), dimension(:), allocatable :: tmparr

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1, n)]
      end if
      allocate(tmparr, mold=arr)
      tmparr(:) = arr(ind(:))
      call qsort_I4B(tmparr, ind)

      return
   end subroutine swiftest_util_sort_index_I4B


   pure module subroutine swiftest_util_sort_index_I4B_I8Bind(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input integer array by index in ascending numerical order using quicksort.
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine swiftest_allocates it.
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:),              intent(in)  :: arr
      integer(I8B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      integer(I8B) :: n, i
      integer(I4B), dimension(:), allocatable :: tmparr

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1_I8B, n)]
      end if
      allocate(tmparr, mold=arr)
      tmparr(:) = arr(ind(:))
      call qsort_I4B_I8Bind(tmparr, ind)

      return
   end subroutine swiftest_util_sort_index_I4B_I8Bind


   recursive pure subroutine swiftest_qsort_I4B(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input I4B array by index in ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:), intent(inout)          :: arr
      integer(I4B), dimension(:), intent(out),  optional :: ind
      ! Internals
      integer(I4B) :: iq

      if (size(arr) > 1) then
         if (present(ind)) then
            call partition_I4B(arr, iq, ind)
            call qsort_I4B(arr(:iq-1),ind(:iq-1))
            call qsort_I4B(arr(iq:),  ind(iq:))
         else
            call partition_I4B(arr, iq)
            call qsort_I4B(arr(:iq-1))
            call qsort_I4B(arr(iq:))
         end if
      end if

      return
   end subroutine swiftest_qsort_I4B

   recursive pure subroutine swiftest_qsort_I4B_I8Bind(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input I4B array by index in ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:), intent(inout)          :: arr
      integer(I8B), dimension(:), intent(out),  optional :: ind
      ! Internals
      integer(I8B) :: iq

      if (size(arr) > 1_I8B) then
         if (present(ind)) then
            call partition_I4B_I8Bind(arr, iq, ind)
            call qsort_I4B_I8Bind(arr(:iq-1_I8B),ind(:iq-1_I8B))
            call qsort_I4B_I8Bind(arr(iq:),  ind(iq:))
         else
            call partition_I4B_I8Bind(arr, iq)
            call qsort_I4B_I8Bind(arr(:iq-1_I8B))
            call qsort_I4B_I8Bind(arr(iq:))
         end if
      end if

      return
   end subroutine swiftest_qsort_I4B_I8Bind


   recursive pure subroutine swiftest_qsort_I8B_I8Bind(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input I8B array by index in ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      integer(I8B), dimension(:), intent(inout)          :: arr
      integer(I8B), dimension(:), intent(out),  optional :: ind
      ! Internals
      integer(I8B) :: iq

      if (size(arr) > 1_I8B) then
         if (present(ind)) then
            call partition_I8B_I8Bind(arr, iq, ind)
            call qsort_I8B_I8Bind(arr(:iq-1_I8B),ind(:iq-1_I8B))
            call qsort_I8B_I8Bind(arr(iq:),  ind(iq:))
         else
            call partition_I8B_I8Bind(arr, iq)
            call qsort_I8B_I8Bind(arr(:iq-1_I8B))
            call qsort_I8B_I8Bind(arr(iq:))
         end if
      end if

      return
   end subroutine swiftest_qsort_I8B_I8Bind

 
   pure subroutine swiftest_partition_I4B(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on I4B type
      !!
      implicit none
      ! Arguments
      integer(I4B), intent(inout), dimension(:)           :: arr
      integer(I4B), intent(inout), dimension(:), optional :: ind
      integer(I4B), intent(out)                           :: marker
      ! Internals
      integer(I4B) :: i, j, itmp, narr, ipiv
      integer(I4B) :: temp
      integer(I4B) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2
      x = arr(ipiv)
      i = 0
      j = narr + 1
   
      do
         j = j - 1
         do
            if (arr(j) <= x) exit
            j = j - 1
         end do
         i = i + 1
         do
            if (arr(i) >= x) exit
            i = i + 1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine swiftest_partition_I4B

   pure subroutine swiftest_partition_I4B_I8Bind(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on I4B type
      !!
      implicit none
      ! Arguments
      integer(I4B), intent(inout), dimension(:)           :: arr
      integer(I8B), intent(inout), dimension(:), optional :: ind
      integer(I8B), intent(out)                           :: marker
      ! Internals
      integer(I8B) :: i, j, itmp, narr, ipiv
      integer(I4B) :: temp
      integer(I8B) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2_I8B
      x = arr(ipiv)
      i = 0_I8B
      j = narr + 1_I8B
   
      do
         j = j - 1_I8B
         do
            if (arr(j) <= x) exit
            j = j - 1_I8B
         end do
         i = i + 1_I8B
         do
            if (arr(i) >= x) exit
            i = i + 1_I8B
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1_I8B
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine swiftest_partition_I4B_I8Bind

   pure subroutine swiftest_partition_I8B_I8Bind(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on I8B type with I8B index
      !!
      implicit none
      ! Arguments
      integer(I8B), intent(inout), dimension(:)           :: arr
      integer(I8B), intent(inout), dimension(:), optional :: ind
      integer(I8B), intent(out)                           :: marker
      ! Internals
      integer(I8B) :: i, j, itmp, narr, ipiv
      integer(I8B) :: temp
      integer(I8B) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2_I8B
      x = arr(ipiv)
      i = 0_I8B
      j = narr + 1_I8B
   
      do
         j = j - 1_I8B
         do
            if (arr(j) <= x) exit
            j = j - 1_I8B
         end do
         i = i + 1_I8B
         do
            if (arr(i) >= x) exit
            i = i + 1_I8B
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1_I8B
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine swiftest_partition_I8B_I8Bind


   pure module subroutine swiftest_util_sort_sp(arr)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array in place into ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      real(SP), dimension(:), intent(inout) :: arr

      call qsort_SP(arr)

      return
   end subroutine swiftest_util_sort_sp


   pure module subroutine swiftest_util_sort_index_sp(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array by index in ascending numerical order using quicksort.
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine swiftest_allocates it.
      !!
      implicit none
      ! Arguments
      real(SP),     dimension(:),              intent(in)    :: arr
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      integer(I4B) :: n, i
      real(SP), dimension(:), allocatable :: tmparr

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1, n)]
      end if
      allocate(tmparr, mold=arr)
      tmparr(:) = arr(ind(:))
      call qsort_SP(tmparr, ind)
   
      return
   end subroutine swiftest_util_sort_index_sp


   recursive pure subroutine swiftest_qsort_SP(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array by index in ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      real(SP), dimension(:), intent(inout)           :: arr
      integer(I4B),dimension(:),intent(out), optional :: ind
      !! Internals
      integer :: iq

      if (size(arr) > 1) then
         if (present(ind)) then
            call partition_SP(arr, iq, ind)
            call qsort_SP(arr(:iq-1),ind(:iq-1))
            call qsort_SP(arr(iq:),  ind(iq:))
         else
            call partition_SP(arr, iq)
            call qsort_SP(arr(:iq-1))
            call qsort_SP(arr(iq:))
         end if
      end if

      return
   end subroutine swiftest_qsort_SP


   pure subroutine swiftest_partition_SP(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on SP type
      !!
      implicit none
      ! Arguments
      real(SP),     intent(inout), dimension(:)           :: arr
      integer(I4B), intent(inout), dimension(:), optional :: ind
      integer(I4B), intent(out)                           :: marker
      ! Internals
      integer(I4B) :: i, j, itmp, narr, ipiv
      real(SP) :: temp
      real(SP) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2
      x = arr(ipiv)
      i = 0
      j = narr + 1
   
      do
         j = j - 1
         do
            if (arr(j) <= x) exit
            j = j - 1
         end do
         i = i + 1
         do
            if (arr(i) >= x) exit
            i = i + 1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine swiftest_partition_SP


   module subroutine swiftest_util_sort_pl(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a Swiftest massive body object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self      !! Swiftest massive body object
      character(*),       intent(in)    :: sortby    !! Sorting attribute
      logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
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
         case("rbeg", "rend", "vbeg", "Ip", "rot", "k_plpl", "nplpl")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class
            call util_sort_body(pl, sortby, ascending)
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
      logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
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
         case("peri")
            call util_sort(direction * tp%peri(1:ntp), ind)
         case("atp")
            call util_sort(direction * tp%atp(1:ntp), ind)
         case("isperi")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class
            call util_sort_body(tp, sortby, ascending)
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
      integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      associate(n => self%nbody)
         call util_sort_rearrange(self%id,       ind, n)
         call util_sort_rearrange(self%lmask,    ind, n)
         call util_sort_rearrange(self%info,     ind, n)
         call util_sort_rearrange(self%status,   ind, n)
         call util_sort_rearrange(self%ldiscard, ind, n)
         call util_sort_rearrange(pl%lcollision, ind, n)
         call util_sort_rearrange(pl%lencounter, ind, n)
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


   pure module subroutine swiftest_util_sort_rearrange_arr_char_string(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of character string in-place from an index list.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B),          dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                                     intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      character(len=STRMAX), dimension(:), allocatable                :: tmp !! Temporary copy of arry used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_sort_rearrange_arr_char_string


   pure module subroutine swiftest_util_sort_rearrange_arr_DP(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of DP type in-place from an index list.
      implicit none
      ! Arguments
      real(DP),     dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B), dimension(:),              intent(in)  :: ind !! Index to rearrange against
      integer(I4B),                            intent(in)  :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      real(DP), dimension(:), allocatable :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_sort_rearrange_arr_DP


   pure module subroutine swiftest_util_sort_rearrange_arr_DPvec(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of (NDIM,n) DP-type vectors in-place from an index list.
      implicit none
      ! Arguments
      real(DP),     dimension(:,:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B), dimension(:),                intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                              intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      real(DP), dimension(:,:), allocatable :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(:,1:n) = arr(:, ind)
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_sort_rearrange_arr_DPvec


   pure module subroutine swiftest_util_sort_rearrange_arr_I4B(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of integers in-place from an index list.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                             intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      integer(I4B), dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_sort_rearrange_arr_I4B

   pure module subroutine swiftest_util_sort_rearrange_arr_I4B_I8Bind(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of integers in-place from an index list.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I8B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I8B),                            intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      integer(I4B), dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0_I8B) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_sort_rearrange_arr_I4B_I8Bind


   pure module subroutine swiftest_util_sort_rearrange_arr_logical(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of logicals in-place from an index list.
      implicit none
      ! Arguments
      logical,      dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                            intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      logical, dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_sort_rearrange_arr_logical


   pure module subroutine swiftest_util_sort_rearrange_arr_logical_I8Bind(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of logicals in-place from an index list.
      implicit none
      ! Arguments
      logical,      dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I8B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I8B),                            intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      logical, dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_sort_rearrange_arr_logical_I8Bind


   module subroutine swiftest_util_sort_rearrange_arr_info(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of particle information type in-place from an index list.
      implicit none
      ! Arguments
      type(swiftest_particle_info),  dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B),                  dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                                             intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      type(swiftest_particle_info),  dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)

      call util_copy_particle_info_arr(arr, tmp, ind)
      call move_alloc(tmp, arr)

      return
   end subroutine swiftest_util_sort_rearrange_arr_info


   module subroutine swiftest_util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange Swiftest massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      class(swiftest_pl),               intent(inout) :: self !! Swiftest massive body object
      integer(I4B),       dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

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

         call util_sort_rearrange_body(pl, ind)
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
      integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      associate(tp => self, ntp => self%nbody)
         call util_sort_rearrange(tp%nplenc,  ind, ntp)

         if (allocated(tp%k_pltp)) deallocate(tp%k_pltp)

         call util_sort_rearrange_body(tp, ind)
      end associate

      return
   end subroutine swiftest_util_sort_rearrange_tp


   module subroutine swiftest_util_spill_arr_char_string(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of type character strings
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,               dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                                          intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist
      character(len=STRMAX), dimension(:), allocatable                :: tmp          !! Array of values to keep 

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
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
   end subroutine swiftest_util_spill_arr_char_string
   

   module subroutine swiftest_util_spill_arr_DP(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of type DP
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      real(DP), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,  dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
      logical,                             intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist
      real(DP), dimension(:), allocatable                :: tmp          !! Array of values to keep 

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
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
   end subroutine swiftest_util_spill_arr_DP


   module subroutine swiftest_util_spill_arr_DPvec(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of DP vectors with shape (NDIM, n)
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      real(DP), dimension(:,:), allocatable, intent(inout) :: discards     !! Array discards
      logical,  dimension(:),                intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: i, nspill, nkeep, nlist
      real(DP), dimension(:,:), allocatable                :: tmp          !! Array of values to keep 

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (.not.allocated(discards)) then
         allocate(discards(NDIM, nspill))
      else if (size(discards, dim=2) /= nspill) then
         deallocate(discards)
         allocate(discards(NDIM, nspill))
      end if

      do i = 1, NDIM
         discards(i,:) = pack(keeps(i,1:nlist), lspill_list(1:nlist))
      end do
      if (ldestructive) then
         if (nkeep > 0) then
            allocate(tmp(NDIM, nkeep))
            do i = 1, NDIM
               tmp(i, :) = pack(keeps(i, 1:nlist), .not. lspill_list(1:nlist))
            end do
            call move_alloc(tmp, keeps)
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine swiftest_util_spill_arr_DPvec


   module subroutine swiftest_util_spill_arr_I4B(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of type I4B
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      integer(I4B), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,      dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                                 intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist
      integer(I4B), dimension(:), allocatable                :: tmp          !! Array of values to keep 

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
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
   end subroutine swiftest_util_spill_arr_I4B


   module subroutine swiftest_util_spill_arr_I8B(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of type I4B
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      integer(I8B), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      integer(I8B), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,      dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                                 intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist
      integer(I8B), dimension(:), allocatable                :: tmp          !! Array of values to keep 

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
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
   end subroutine swiftest_util_spill_arr_I8B


   module subroutine swiftest_util_spill_arr_info(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of particle origin information types
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,                       dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
      logical,                                                  intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: i, nspill, nkeep, nlist
      integer(I4B), dimension(:), allocatable :: idx
      type(swiftest_particle_info), dimension(:), allocatable :: tmp

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
      if (.not.allocated(discards)) then
         allocate(discards(nspill))
      else if (size(discards) /= nspill) then
         deallocate(discards)
         allocate(discards(nspill))
      end if

      allocate(idx(nspill))
      idx(:) = pack([(i, i = 1, nlist)], lspill_list)
      call util_copy_particle_info_arr(keeps, discards, idx)
      if (ldestructive) then
         if (nkeep > 0) then
            deallocate(idx)
            allocate(idx(nkeep))
            allocate(tmp(nkeep))
            idx(:) = pack([(i, i = 1, nlist)], .not. lspill_list)
            call util_copy_particle_info_arr(keeps, tmp, idx)
            call move_alloc(tmp, keeps)
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine swiftest_util_spill_arr_info


   module subroutine swiftest_util_spill_arr_logical(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of logicals
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      logical, dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical, dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                            intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or no
      ! Internals
      integer(I4B) :: nspill, nkeep, nlist
      logical, dimension(:), allocatable                :: tmp          !! Array of values to keep 

      nkeep = count(.not.lspill_list(:))
      nspill = count(lspill_list(:))
      nlist = size(lspill_list(:))

      if (.not.allocated(keeps) .or. nspill == 0) return
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
   end subroutine swiftest_util_spill_arr_logical


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
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      ! Internals
      integer(I4B) :: nbody_old

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
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
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list

      associate(keeps => self)
         select type (discards) ! The standard requires us to select the type of both arguments in order to access all the components
         class is (swiftest_pl)
            !> Spill components specific to the massive body class
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

            call util_spill_body(keeps, discards, lspill_list, ldestructive)
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
      class(swiftest_tp),    intent(inout) :: self        !! Swiftest test particle object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discardse
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list

      associate(keeps => self, ntp => self%nbody)
         select type(discards)
         class is (swiftest_tp)
            !> Spill components specific to the test particle class
            call util_spill(keeps%nplenc,  discards%nplenc,  lspill_list, ldestructive)
            call util_spill_body(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on swiftest_tp'
         end select
      end associate

      return
   end subroutine swiftest_util_spill_tp


   module subroutine swiftest_util_unique_DP(input_array, output_array, index_map)
      !! author: David A. Minton
      !!
      !! Takes an input unsorted integer array and returns a new array of sorted, unique values (DP version)
      implicit none
      ! Arguments
      real(DP),     dimension(:),              intent(in)  :: input_array  !! Unsorted input array 
      real(DP),     dimension(:), allocatable, intent(out) :: output_array !! Sorted array of unique values 
      integer(I4B), dimension(:), allocatable, intent(out) :: index_map    !! An array of the same size as input_array that such that any for any index i, output_array(index_map(i)) = input_array(i)       
      ! Internals
      real(DP), dimension(:), allocatable :: unique_array
      integer(I4B) :: n
      real(DP) :: lo, hi

      allocate(unique_array, mold=input_array)
      allocate(index_map(size(input_array)))
      lo = minval(input_array) - 1
      hi = maxval(input_array)

      n = 0
      do 
         n = n + 1
         lo = minval(input_array(:), mask=input_array(:) > lo)
         unique_array(n) = lo
         where(input_array(:) == lo) index_map(:) = n
         if (lo >= hi) exit
      enddo
      allocate(output_array(n), source=unique_array(1:n)) 

      return
   end subroutine swiftest_util_unique_DP


   module subroutine swiftest_util_unique_I4B(input_array, output_array, index_map)
      !! author: David A. Minton
      !!
      !! Takes an input unsorted integer array and returns a new array of sorted, unique values (I4B version)
      implicit none
      ! Arguments
      integer(I4B), dimension(:),              intent(in)  :: input_array  !! Unsorted input array 
      integer(I4B), dimension(:), allocatable, intent(out) :: output_array !! Sorted array of unique values
      integer(I4B), dimension(:), allocatable, intent(out) :: index_map    !! An array of the same size as input_array that such that any for any index i, output_array(index_map(i)) = input_array(i)     
      ! Internals
      integer(I4B), dimension(:), allocatable :: unique_array
      integer(I4B) :: n, lo, hi

      allocate(unique_array, mold=input_array)
      allocate(index_map, mold=input_array)
      lo = minval(input_array) - 1
      hi = maxval(input_array)

      n = 0
      do 
         n = n + 1
         lo = minval(input_array(:), mask=input_array(:) > lo)
         unique_array(n) = lo
         where(input_array(:) == lo) index_map(:) = n
         if (lo >= hi) exit
      enddo
      allocate(output_array(n), source=unique_array(1:n)) 

      return
   end subroutine swiftest_util_unique_I4B


   module subroutine swiftest_util_valid_id_system(self, param)
      !! author: David A. Minton
      !!
      !! Validate massive body and test particle ids
      !! subroutine swiftest_causes program to exit with error if any ids are not unique
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_valid.f90
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                  :: i
      integer(I4B), dimension(:), allocatable :: idarr

      associate(cb => self%cb, pl => self%pl, npl => self%pl%nbody, tp => self%tp, ntp => self%tp%nbody)
         allocate(idarr(1+npl+ntp))
         idarr(1) = cb%id
         do i = 1, npl
            idarr(1+i) = pl%id(i)
         end do
         do i = 1, ntp
            idarr(1+npl+i) = tp%id(i)
         end do
         call util_sort(idarr)
         do i = 1, npl + ntp 
            if (idarr(i) == idarr(i+1)) then
               write(*, *) "Swiftest error:"
               write(*, *) "   more than one body/particle has id = ", idarr(i)
               call util_exit(FAILURE)
            end if
         end do
         param%maxid = max(param%maxid, maxval(idarr))
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
      write(*, 200) VERSION_NUMBER
      200 format(/, "************* Swiftest: Version ", f3.1, " *************", //, &
            "Based off of Swifter:", //,                                         &
            "Authors:", //,                                                      &
            "    The Purdue University Swiftest Development team ", /,           &
            "    Lead by David A. Minton ", /,                                   &
            "    Single loop blocking by Jacob R. Elliott", /,                   &
            "    Fragmentation by Carlisle A. Wishard and", //,                  &
            "    Jennifer L. L. Poutplin                 ", //,                  &
            "Please address comments and questions to:", //,                     &
            "    David A. Minton", /,                                            &
            "    Department Earth, Atmospheric, & Planetary Sciences ",/,        &
            "    Purdue University", /,                                          &
            "    550 Stadium Mall Drive", /,                                     &
            "    West Lafayette, Indiana 47907", /,                              &
            "    765-250-8034 ", /,                                              &
            "    daminton@purdue.edu", /,                                        &
            "Special thanks to Hal Levison and Martin Duncan for the original",/,&
            "SWIFTER and SWIFT codes that made this possible.", //,              &
            "************************************************", /)


      100 FORMAT(/,  "************* SWIFTER: Version ", F3.1, " *************", //, &
                  "Authors:", //,                                                &
                  "    Martin Duncan: Queen's University", /,                    &
                  "    Hal Levison  : Southwest Research Institute", //,         &
                  "Please address comments and questions to:", //,               &
                  "    Hal Levison or David Kaufmann", /,                        &
                  "    Department of Space Studies", /,                          &
                  "    Southwest Research Institute", /,                         &
                  "    1050 Walnut Street, Suite 400", /,                        &
                  "    Boulder, Colorado  80302", /,                             &
                  "    303-546-0290 (HFL), 720-240-0119 (DEK)", /,               &
                  "    303-546-9687 (fax)", /,                                   &
                  "    hal@gort.boulder.swri.edu (HFL)", /,                      &
                  "    kaufmann@boulder.swri.edu (DEK)", //,                     &
                  "************************************************", /)

      return
   end subroutine swiftest_util_version

end submodule s_util