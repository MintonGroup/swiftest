!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (encounter) s_encounter_check
   use swiftest
contains

   module subroutine encounter_check_all_plpl(param, npl, r, v, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(base_parameters),                  intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies
      real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter 
      real(DP),     dimension(:,:),            intent(in)    :: r      !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: v      !! Velocity vectors of massive bodies
      real(DP),                                intent(in)    :: dt     !! Step size
      integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x

      if (param%lencounter_sas_plpl) then
         call encounter_check_all_sort_and_sweep_plpl(npl, r, v, renc, dt, nenc, index1, index2, lvdotr)
      else
         call encounter_check_all_triangular_plpl(npl, r, v, renc, dt, nenc, index1, index2, lvdotr) 
      end if

      return
   end subroutine encounter_check_all_plpl


   module subroutine encounter_check_all_plplm(param, nplm, nplt, rplm, vplm, rplt, vplt, rencm, renct, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between fully interacting massive bodies partially interacting massive bodies. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(base_parameters),                  intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: nplm   !! Total number of fully interacting massive bodies 
      integer(I4B),                            intent(in)    :: nplt   !! Total number of partially interacting masive bodies (GM < GMTINY) 
      real(DP),     dimension(:,:),            intent(in)    :: rplm   !! Position vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vplm   !! Velocity vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: rplt   !! Position vectors of partially interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vplt   !! Velocity vectors of partially interacting massive bodies
      real(DP),     dimension(:),              intent(in)    :: rencm  !! Critical radii of fully interacting massive bodies that defines an encounter
      real(DP),     dimension(:),              intent(in)    :: renct  !! Critical radii of partially interacting massive bodies that defines an encounter
      real(DP),                                intent(in)    :: dt     !! Step size
      integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      logical,      dimension(:), allocatable :: plmplt_lvdotr !! Logical flag indicating the sign of v .dot. x in the plm-plt group
      integer(I4B), dimension(:), allocatable :: plmplt_index1 !! List of indices for body 1 in each encounter in the plm-plt group
      integer(I4B), dimension(:), allocatable :: plmplt_index2 !! List of indices for body 2 in each encounter in the plm-lt group
      integer(I8B)                            :: plmplt_nenc   !! Number of encounters of the plm-plt group
      class(base_parameters), allocatable :: tmp_param     !! Temporary parameter structure to turn off adaptive timer for the pl-pl phase if necessary
      integer(I8B), dimension(:), allocatable :: ind
      integer(I4B), dimension(:), allocatable :: itmp
      logical, dimension(:), allocatable :: ltmp

      allocate(tmp_param, source=param)

      ! Start with the pl-pl group
      call encounter_check_all_plpl(tmp_param, nplm, rplm, vplm, rencm, dt, nenc, index1, index2, lvdotr)

      if (param%lencounter_sas_plpl) then
         call encounter_check_all_sort_and_sweep_plplm(nplm, nplt, rplm, vplm, rplt, vplt, rencm, renct, dt, &
                                                       plmplt_nenc, plmplt_index1, plmplt_index2, plmplt_lvdotr)
      else
         call encounter_check_all_triangular_plplm(nplm, nplt, rplm, vplm, rplt, vplt, rencm, renct, dt, &
                                                   plmplt_nenc, plmplt_index1, plmplt_index2, plmplt_lvdotr) 
      end if

      if (plmplt_nenc > 0) then ! Consolidate the two lists
         allocate(itmp(nenc+plmplt_nenc))
         if (nenc > 0) itmp(1:nenc) = index1(1:nenc)
         itmp(nenc+1:nenc+plmplt_nenc) = plmplt_index1(1:plmplt_nenc)
         call move_alloc(itmp, index1)
         allocate(itmp(nenc+plmplt_nenc))
         if (nenc > 0) itmp(1:nenc) = index2(1:nenc)
         itmp(nenc+1:nenc+plmplt_nenc) = plmplt_index2(1:plmplt_nenc) + nplm ! Be sure to shift these indices back to their natural range
         call move_alloc(itmp, index2)
         allocate(ltmp(nenc+plmplt_nenc))
         if (nenc > 0) ltmp(1:nenc) = lvdotr(1:nenc)
         ltmp(nenc+1:nenc+plmplt_nenc) = plmplt_lvdotr(1:plmplt_nenc)
         call move_alloc(ltmp, lvdotr)
         nenc = nenc + plmplt_nenc

         call util_sort(index1, ind)
         call util_sort_rearrange(index1, ind, nenc)
         call util_sort_rearrange(index2, ind, nenc)
         call util_sort_rearrange(lvdotr, ind, nenc)

      end if

      return
   end subroutine encounter_check_all_plplm


   module subroutine encounter_check_all_pltp(param, npl, ntp, rpl, vpl, rtp, vtp, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies and test particles. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(base_parameters),                  intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies 
      integer(I4B),                            intent(in)    :: ntp    !! Total number of test particles 
      real(DP),     dimension(:,:),            intent(in)    :: rpl    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vpl    !! Velocity vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: rtp    !! Position vectors of test particlse
      real(DP),     dimension(:,:),            intent(in)    :: vtp    !! Velocity vectors of test particles
      real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter
      real(DP),                                intent(in)    :: dt     !! Step size
      integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x

      if (param%lencounter_sas_pltp) then
         call encounter_check_all_sort_and_sweep_pltp(npl, ntp, rpl, vpl, rtp, vtp, renc, dt, nenc, index1, index2, lvdotr)
      else
         call encounter_check_all_triangular_pltp(npl, ntp, rpl, vpl, rtp, vtp, renc, dt, nenc, index1, index2, lvdotr) 
      end if

      return
   end subroutine encounter_check_all_pltp


   subroutine encounter_check_all_sort_and_sweep_plpl(npl, r, v, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. 
      !! This is the sort and sweep version
      !! References: Adapted from Christer Ericson's _Real Time Collision Detection_
      !!
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: npl    !! Total number of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: r      !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: v      !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)  :: renc   !! Critical radii of massive bodies that defines an encounter 
      real(DP),                                intent(in)  :: dt     !! Step size
      integer(I8B),                            intent(out) :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      integer(I4B) :: i, n
      integer(I4B), save :: npl_last = 0
      type(encounter_bounding_box), save :: boundingbox
      real(DP), dimension(npl) :: rmin,rmax
      real(DP) :: rmag

      if (npl == 0) return

      ! If this is the first time through, build the index lists
      n = 2 * npl
      if (npl_last /= npl) then
         call boundingbox%setup(npl, npl_last)
         npl_last = npl
      end if

#ifdef DOCONLOC
      do concurrent (i = 1:npl) shared(r,renc,rmin,rmax) local(rmag)
#else
      do concurrent (i = 1:npl)
#endif
         rmag = norm2(r(:,i))
         rmax(i) = rmag + RSWEEP_FACTOR * renc(i)
         rmin(i) = rmag - RSWEEP_FACTOR * renc(i)
      end do

      call boundingbox%aabb%sort(npl, [rmin,rmax]) 

      call boundingbox%sweep(npl, r, v, renc, dt, nenc, index1, index2, lvdotr) 

      return
   end subroutine encounter_check_all_sort_and_sweep_plpl


   subroutine encounter_check_all_sort_and_sweep_plplm(nplm, nplt, rplm, vplm, rplt, vplt, rencm, renct, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies and test particles. 
      !! This is the sort and sweep version
      !! References: Adapted from Christer Ericson's _Real Time Collision Detection_
      !!
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: nplm   !! Total number of fully interacting massive bodies 
      integer(I4B),                            intent(in)  :: nplt   !! Total number of partially interacting masive bodies (GM < GMTINY) 
      real(DP),     dimension(:,:),            intent(in)  :: rplm   !! Position vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vplm   !! Velocity vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: rplt   !! Position vectors of partially interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vplt   !! Velocity vectors of partially interacting massive bodies
      real(DP),     dimension(:),              intent(in)  :: rencm  !! Critical radii of fully interacting massive bodies that defines an encounter
      real(DP),     dimension(:),              intent(in)  :: renct  !! Critical radii of partially interacting massive bodies that defines an encounter
      real(DP),                                intent(in)  :: dt     !! Step size
      integer(I8B),                            intent(out) :: nenc   !! Total number of encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      type(encounter_bounding_box), save :: boundingbox
      integer(I4B) :: i, n, ntot
      integer(I4B), save :: ntot_last = 0
      real(DP), dimension(nplm+nplt) :: rmin,rmax
      real(DP) :: rmag

      ! If this is the first time through, build the index lists
      if ((nplm == 0) .or. (nplt == 0)) return

      ntot = nplm + nplt
      n = 2 * ntot
      if (ntot /= ntot_last) then
         call boundingbox%setup(ntot, ntot_last)
         ntot_last = ntot
      end if

#ifdef DOCONLOC
      do concurrent (i = 1:nplm) shared(rmin,rmax,rplm,rencm) local(rmag)
#else
      do concurrent (i = 1:nplm)
#endif
         rmag = norm2(rplm(:,i))
         rmax(i) = rmag + RSWEEP_FACTOR * rencm(i)
         rmin(i) = rmag - RSWEEP_FACTOR * rencm(i)
      end do
#ifdef DOCONLOC
      do concurrent (i = 1:nplt) shared(rmin,rmax,rplt,renct) local(rmag)
#else
      do concurrent (i = 1:nplt)
#endif
         rmag = norm2(rplt(:,i))
         rmax(nplm+i) = rmag + RSWEEP_FACTOR * renct(i)
         rmin(nplm+i) = rmag - RSWEEP_FACTOR * renct(i)
      end do

      call boundingbox%aabb%sort(ntot, [rmin, rmax])

      call boundingbox%sweep(nplm, nplt, rplm, vplm, rplt, vplt, rencm, renct, dt, nenc, index1, index2, lvdotr) 

      return
   end subroutine encounter_check_all_sort_and_sweep_plplm


   subroutine encounter_check_all_sort_and_sweep_pltp(npl, ntp, rpl, vpl, rtp, vtp, rencpl, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies and test particles. 
      !! This is the sort and sweep version
      !! References: Adapted from Christer Ericson's _Real Time Collision Detection_
      !!
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: npl    !! Total number of massive bodies 
      integer(I4B),                            intent(in)  :: ntp    !! Total number of test particles 
      real(DP),     dimension(:,:),            intent(in)  :: rpl    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vpl    !! Velocity vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: rtp    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vtp    !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)  :: rencpl !! Critical radii of massive bodies that defines an encounter
      real(DP),                                intent(in)  :: dt     !! Step size
      integer(I8B),                            intent(out) :: nenc   !! Total number of encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      type(encounter_bounding_box), save :: boundingbox
      integer(I4B) :: i, n, ntot
      integer(I4B), save :: ntot_last = 0
      real(DP), dimension(npl+ntp) :: rmin,rmax
      real(DP), dimension(ntp) :: renctp
      real(DP) :: rmag

      ! If this is the first time through, build the index lists
      if ((ntp == 0) .or. (npl == 0)) return

      ntot = npl + ntp
      n = 2 * ntot
      if (ntot /= ntot_last) then
         call boundingbox%setup(ntot, ntot_last)
         ntot_last = ntot
      end if

      renctp(:) = 0.0_DP

#ifdef DOCONLOC
      do concurrent (i = 1:npl) shared(rmin,rmax,rpl,rencpl) local(rmag)
#else
      do concurrent (i = 1:npl)
#endif
         rmag = norm2(rpl(:,i))
         rmax(i) = rmag + RSWEEP_FACTOR * rencpl(i)
         rmin(i) = rmag - RSWEEP_FACTOR * rencpl(i)
      end do
#ifdef DOCONLOC
      do concurrent (i = 1:ntp) shared(rmin,rmax,rtp,renctp) local(rmag)
#else
      do concurrent (i = 1:ntp) 
#endif
         rmag = norm2(rtp(:,i))
         rmax(npl+i) = rmag + RSWEEP_FACTOR * renctp(i)
         rmin(npl+i) = rmag - RSWEEP_FACTOR * renctp(i)
      end do

      call boundingbox%aabb%sort(ntot, [rmin, rmax])

      call boundingbox%sweep(npl, ntp, rpl, vpl, rtp, vtp, rencpl, renctp, dt, nenc, index1, index2, lvdotr) 

      return
   end subroutine encounter_check_all_sort_and_sweep_pltp


   pure subroutine encounter_check_all_sweep_one(i, n, xi, yi, zi, vxi, vyi, vzi, x, y, z, vx, vy, vz, renci, renc, dt, &
                                                 ind_arr, lgood, nenci, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between the ith body and all other bodies.
      !! This is used in the narrow phase of the sort & sweep algorithm
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)    :: i             !! Index of the ith body that is being checked
      integer(I4B),                            intent(in)    :: n             !! Total number of bodies being checked
      real(DP),                                intent(in)    :: xi, yi, zi    !! Position vector components of the ith body
      real(DP),                                intent(in)    :: vxi, vyi, vzi !! Velocity vector components of the ith body
      real(DP),     dimension(:),              intent(in)    :: x, y, z       !! Arrays of position vector components of all bodies
      real(DP),     dimension(:),              intent(in)    :: vx, vy, vz    !! Arrays of velocity vector components of all bodies
      real(DP),                                intent(in)    :: renci         !! Encounter radius of the ith body
      real(DP),     dimension(:),              intent(in)    :: renc          !! Array of encounter radii of all bodies
      real(DP),                                intent(in)    :: dt            !! Step size
      integer(I4B), dimension(:),              intent(in)    :: ind_arr       !! Index array [1, 2, ..., n]
      logical,      dimension(:),              intent(in)    :: lgood         !! Logical array mask where true values correspond to bodies selected in the broad phase
      integer(I8B),                            intent(out)   :: nenci         !! Total number of encountering bodies
      integer(I4B), dimension(:), allocatable, intent(inout) :: index1        !! Array of indices of the ith body of size nenci [i, i, ..., i]
      integer(I4B), dimension(:), allocatable, intent(inout) :: index2        !! Array of indices of the encountering bodies of size nenci 
      logical,      dimension(:), allocatable, intent(inout) :: lvdotr        !! v.dot.r direction array
      ! Internals
      integer(I4B) :: j
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc12
      logical, dimension(n) :: lencounteri, lvdotri

      lencounteri(:) = .false.
#ifdef DOCONLOC
      do concurrent(j = 1:n, lgood(j)) shared(lgood,lencounteri,lvdotri,x,y,z,vx,vy,vz,renci,renc) local(xr,yr,zr,vxr,vyr,vzr,renc12)
#else
      do concurrent(j = 1:n, lgood(j))
#endif
         xr = x(j) - xi
         yr = y(j) - yi
         zr = z(j) - zi
         vxr = vx(j) - vxi
         vyr = vy(j) - vyi
         vzr = vz(j) - vzi
         renc12 = renci + renc(j)
         call encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc12, dt, lencounteri(j), lvdotri(j))
      end do
      if (any(lencounteri(:))) then
         nenci = count(lencounteri(:))
         allocate(lvdotr(nenci), index1(nenci), index2(nenci))
         index1(:) = i
         index2(:) = pack(ind_arr(1:n), lencounteri(1:n)) 
         lvdotr(:) = pack(lvdotri(1:n), lencounteri(1:n)) 
      end if

      return
   end subroutine encounter_check_all_sweep_one


   subroutine encounter_check_all_triangular_one(i, n, xi, yi, zi, vxi, vyi, vzi, x, y, z, vx, vy, vz, renci, renc, &
                                                      dt, ind_arr, lenci)
      !! author: David A. Minton
      !!
      !! Check for encounters between the ith body and all other bodies.
      !! This is the upper triangular (double loop) version.
      implicit none
      ! Arguments
      integer(I4B),                       intent(in)  :: i             !! Index of the ith body that is being checked
      integer(I4B),                       intent(in)  :: n             !! Total number of bodies being checked
      real(DP),                           intent(in)  :: xi, yi, zi    !! Position vector components of the ith body
      real(DP),                           intent(in)  :: vxi, vyi, vzi !! Velocity vector components of the ith body
      real(DP),             dimension(:), intent(in)  :: x, y, z       !! Arrays of position vector components of all bodies
      real(DP),             dimension(:), intent(in)  :: vx, vy, vz    !! Arrays of velocity vector components of all bodies
      real(DP),                           intent(in)  :: renci         !! Encounter radius of the ith body
      real(DP),             dimension(:), intent(in)  :: renc          !! Array of encounter radii of all bodies
      real(DP),                           intent(in)  :: dt            !! Step size
      integer(I4B),         dimension(:), intent(in)  :: ind_arr       !! Index array [1, 2, ..., n]
      class(encounter_list),              intent(out) :: lenci         !! Output encounter lists containing number of encounters, the v.dot.r direction array, and the index list of encountering bodies 
      ! Internals
      integer(I4B) :: j
      integer(I8B) :: nenci
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc12
      logical, dimension(n) :: lencounteri, lvdotri

#ifdef DOCONLOC
      do concurrent(j = i+1:n) shared(lencounteri, lvdotri, renci, renc) local(xr,yr,zr,vxr,vyr,vzr,renc12)
#else
      do concurrent(j = i+1:n)
#endif
         xr = x(j) - xi
         yr = y(j) - yi
         zr = z(j) - zi
         vxr = vx(j) - vxi
         vyr = vy(j) - vyi
         vzr = vz(j) - vzi
         renc12 = renci + renc(j)
         call encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc12, dt, lencounteri(j), lvdotri(j))
      end do
      nenci = count(lencounteri(i+1:n))
      if (nenci > 0_I8B) then
         allocate(lenci%lvdotr(nenci), lenci%index1(nenci), lenci%index2(nenci))
         lenci%nenc = nenci
         lenci%index1(:) = i
         lenci%index2(:) = pack(ind_arr(i+1:n), lencounteri(i+1:n)) 
         lenci%lvdotr(:) = pack(lvdotri(i+1:n), lencounteri(i+1:n)) 
      end if

      return
   end subroutine encounter_check_all_triangular_one


   subroutine encounter_check_all_triangular_plpl(npl, r, v, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Split off from the main subroutine for performance
      !! This is the upper triangular (double loop) version.
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: npl    !! Total number of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: r      !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: v      !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)  :: renc   !! Critical radii of massive bodies that defines an encounter 
      real(DP),                                intent(in)  :: dt     !! Step size
      integer(I8B),                            intent(out) :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      integer(I4B) :: i
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      type(collision_list_plpl), dimension(npl) :: lenc

      call swiftest_util_index_array(ind_arr, npl) 

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(r, v, renc, lenc, ind_arr) &
      !$omp firstprivate(npl, dt)
      do i = 1,npl
         call encounter_check_all_triangular_one(i, npl, r(1,i), r(2,i), r(3,i), &
                                                         v(1,i), v(2,i), v(3,i), &
                                                         r(1,:), r(2,:), r(3,:), &
                                                         v(1,:), v(2,:), v(3,:), &
                                                         renc(i), renc(:), dt, ind_arr(:), lenc(i))
         if (lenc(i)%nenc > 0) lenc(i)%index1(:) = i
      end do
      !$omp end parallel do

      call encounter_check_collapse_ragged_list(lenc, npl, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_all_triangular_plpl


   subroutine encounter_check_all_triangular_plplm(nplm, nplt, rplm, vplm, rplt, vplt, rencm, renct, dt, &
                                                   nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies and test particles. Split off from the main subroutine for performance
      !! This is the upper triangular (double loop) version.
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: nplm   !! Total number of fully interacting massive bodies 
      integer(I4B),                            intent(in)  :: nplt   !! Total number of partially interacting masive bodies (GM < GMTINY) 
      real(DP),     dimension(:,:),            intent(in)  :: rplm   !! Position vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vplm   !! Velocity vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: rplt   !! Position vectors of partially interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vplt   !! Velocity vectors of partially interacting massive bodies
      real(DP),     dimension(:),              intent(in)  :: rencm  !! Critical radii of fully interacting massive bodies that defines an encounter
      real(DP),     dimension(:),              intent(in)  :: renct  !! Critical radii of partially interacting massive bodies that defines an encounter
      real(DP),                                intent(in)  :: dt     !! Step size
      integer(I8B),                            intent(out) :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      integer(I4B) :: i
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      type(collision_list_plpl), dimension(nplm) :: lenc

      call swiftest_util_index_array(ind_arr, nplt)

      !$omp parallel do default(private) schedule(dynamic)&
      !$omp shared(rplm, vplm, rplt, vplt, rencm, renct, lenc, ind_arr) &
      !$omp firstprivate(nplm, nplt, dt)
      do i = 1, nplm
         call encounter_check_all_triangular_one(0, nplt, rplm(1,i), rplm(2,i), rplm(3,i), &
                                                          vplm(1,i), vplm(2,i), vplm(3,i), &
                                                          rplt(1,:), rplt(2,:), rplt(3,:), &
                                                          vplt(1,:), vplt(2,:), vplt(3,:), &
                                                          rencm(i), renct(:), dt, ind_arr(:), lenc(i))
         if (lenc(i)%nenc > 0) lenc(i)%index1(:) = i
      end do
      !$omp end parallel do

      call encounter_check_collapse_ragged_list(lenc, nplm, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_all_triangular_plplm


   subroutine encounter_check_all_triangular_pltp(npl, ntp, rpl, vpl, rtp, vtp, renc, dt, &
                                                  nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies and test particles. Split off from the main subroutine for performance
      !! This is the upper triangular (double loop) version.
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: npl    !! Total number of massive bodies 
      integer(I4B),                            intent(in)  :: ntp    !! Total number of test particles 
      real(DP),     dimension(:,:),            intent(in)  :: rpl    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vpl    !! Velocity vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: rtp    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vtp    !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)  :: renc   !! Critical radii of massive bodies that defines an encounter
      real(DP),                                intent(in)  :: dt     !! Step size
      integer(I8B),                            intent(out) :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      integer(I4B) :: i
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      type(collision_list_pltp), dimension(npl) :: lenc
      real(DP), dimension(ntp) :: renct

      call swiftest_util_index_array(ind_arr, ntp)
      renct(:) = 0.0_DP

      !$omp parallel do default(private) schedule(dynamic)&
      !$omp shared(rpl, vpl, rtp, vtp, renc, renct, lenc, ind_arr) &
      !$omp firstprivate(npl, ntp, dt)
      do i = 1, npl
         call encounter_check_all_triangular_one(0, ntp, rpl(1,i), rpl(2,i), rpl(3,i), &
                                                         vpl(1,i), vpl(2,i), vpl(3,i), &
                                                         rtp(1,:), rtp(2,:), rtp(3,:), &
                                                         vtp(1,:), vtp(2,:), vtp(3,:), &
                                                         renc(i), renct(:), dt, ind_arr(:), lenc(i))
         if (lenc(i)%nenc > 0) lenc(i)%index1(:) = i
      end do
      !$omp end parallel do

      call encounter_check_collapse_ragged_list(lenc, npl, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_all_triangular_pltp


   elemental module subroutine encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc, dt, lencounter, lvdotr)
      !! author: David A. Minton
      !!
      !! Determine whether a test particle and planet are having or will have an encounter within the next time step
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk_ind.f90
      !! Adapted from Hal Levison's Swift routine rmvs_chk_ind.f
      implicit none
      ! Arguments
      real(DP), intent(in)  :: xr, yr, zr    !! Relative distance vector components
      real(DP), intent(in)  :: vxr, vyr, vzr !! Relative velocity vector components
      real(DP), intent(in)  :: renc        !! Square of the critical encounter distance
      real(DP), intent(in)  :: dt            !! Step size
      logical,  intent(out) :: lencounter    !! Flag indicating that an encounter has occurred
      logical,  intent(out) :: lvdotr        !! Logical flag indicating the direction of the v .dot. r vector
      ! Internals
      real(DP) :: r2crit, r2min, r2, v2, vdotr, tmin

      r2 = xr**2 + yr**2 + zr**2
      r2crit = renc**2

      if (r2 > r2crit) then 
         vdotr = vxr * xr + vyr * yr + vzr * zr
         if (vdotr > 0.0_DP) then
            r2min = r2
         else
            v2 = vxr**2 + vyr**2 + vzr**2
            if (v2 <= VSMALL) then
               r2min = r2
            else
               tmin = -vdotr / v2
         
               if (tmin < dt) then
                  r2min = r2 - vdotr**2 / v2
               else
                  r2min = r2 + 2 * vdotr * dt + v2 * dt**2
               end if
            end if
         end if
      else
         vdotr = -1.0_DP
         r2min = r2
      end if

      lvdotr = (vdotr < 0.0_DP)
      lencounter = lvdotr .and. (r2min <= r2crit) 

      return
   end subroutine encounter_check_one


   module subroutine encounter_check_collapse_ragged_list(ragged_list, n1, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!    
      !! Collapses a ragged index list (one encounter list per body) into a pair of index arrays and a vdotr logical array (optional)
      implicit none
      ! Arguments
      class(encounter_list), dimension(:),              intent(in)            :: ragged_list !! The ragged encounter list
      integer(I4B),                                     intent(in)            :: n1          !! Number of bodies 1
      integer(I8B),                                     intent(out)           :: nenc        !! Total number of encountersj 
      integer(I4B),          dimension(:), allocatable, intent(out)           :: index1      !! Array of indices for body 1
      integer(I4B),          dimension(:), allocatable, intent(out)           :: index2      !! Array of indices for body 1
      logical,               dimension(:), allocatable, intent(out), optional :: lvdotr      !! Array indicating which bodies are approaching
      ! Internals
      integer(I4B) :: i
      integer(I8B) :: j1, j0, nenci
      integer(I8B), dimension(n1) :: ibeg

      associate(nenc_arr => ragged_list(:)%nenc)
         nenc = sum(nenc_arr(:))
      end associate
      if (nenc == 0) return

      allocate(index1(nenc))
      allocate(index2(nenc))
      if (present(lvdotr)) allocate(lvdotr(nenc))

      j0 = 1
      do i = 1, n1
         nenci = ragged_list(i)%nenc
         if (nenci == 0) cycle
         ibeg(i) = j0
         j0 = j0 + nenci
      end do

      !$omp parallel do simd default(private) schedule(simd: static)&
      !$omp shared(ragged_list, index1, index2, ibeg, lvdotr) &
      !$omp firstprivate(n1)
      do i = 1,n1
         if (ragged_list(i)%nenc == 0_I8B) cycle
         nenci = ragged_list(i)%nenc
         j0 = ibeg(i)
         j1 = j0 + nenci - 1_I8B
         index1(j0:j1) = ragged_list(i)%index1(:)
         index2(j0:j1) = ragged_list(i)%index2(:)
         if (present(lvdotr)) lvdotr(j0:j1) = ragged_list(i)%lvdotr(:)
      end do
      !$omp end parallel do simd

      return
   end subroutine encounter_check_collapse_ragged_list


   subroutine encounter_check_remove_duplicates(n, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Takes the candidate encounter lists that came out of the sort & sweep method and remove any duplicates.
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)    :: n          !! Number of bodies 
      integer(I8B),                            intent(inout) :: nenc       !! Number of encountering bodies (input is the broad phase value, output is the final narrow phase value)
      integer(I4B), dimension(:), allocatable, intent(inout) :: index1     !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(inout) :: index2     !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(inout) :: lvdotr     !! Logical flag indicating the sign of v .dot. x
      ! Internals
      integer(I4B) :: i, i0
      integer(I8B) :: j, k, klo, khi, nenci
      integer(I4B), dimension(:), allocatable :: itmp
      integer(I8B), dimension(:), allocatable :: ind
      integer(I8B), dimension(n) :: ibeg, iend
      logical, dimension(:), allocatable :: ltmp
      logical, dimension(nenc) :: lencounter

      if (nenc == 0) then
         if (allocated(index1)) deallocate(index1)
         if (allocated(index2)) deallocate(index2)
         if (allocated(lvdotr)) deallocate(lvdotr)
         return
      end if

      call util_sort(index1, ind)
      call util_sort_rearrange(index1, ind, nenc)
      call util_sort_rearrange(index2, ind, nenc)
      call util_sort_rearrange(lvdotr, ind, nenc)

      ! Get the bounds on the bodies in the first index
      ibeg(:) = n
      iend(:) = 1_I8B
      i0 = index1(1) 
      ibeg(i0) = 1_I8B
      do k = 2_I8B, nenc 
         i = index1(k)
         if (i /= i0) then
            iend(i0) = k - 1_I8B
            ibeg(i) = k
            i0 = i
         end if
         if (k == nenc) iend(i) = k
      end do

      lencounter(:) = .true.
      ! Sort on the second index and remove duplicates 
      if (allocated(itmp)) deallocate(itmp)
      allocate(itmp, source=index2)
#ifdef DOCONLOC
      do concurrent(i = 1:n, iend(i) - ibeg(i) > 0_I8B) shared(iend,ibeg,index2,lencounter,itmp) local(klo,khi,nenci,j)
#else
      do concurrent(i = 1:n, iend(i) - ibeg(i) > 0_I8B)
#endif
         klo = ibeg(i)
         khi = iend(i)
         nenci = khi - klo + 1_I8B
         if (allocated(ind)) deallocate(ind)
         call util_sort(index2(klo:khi), ind)
         index2(klo:khi) = itmp(klo - 1_I8B + ind(:))
         do j = klo + 1_I8B, khi
            if (index2(j) == index2(j - 1_I8B)) lencounter(j) = .false. 
         end do
      end do

      if (count(lencounter(:)) == nenc) return
      nenc = count(lencounter(:)) ! Count the true number of encounters

      if (allocated(itmp)) deallocate(itmp)
      allocate(itmp(nenc))
      itmp(:) = pack(index1(:), lencounter(:))
      call move_alloc(itmp, index1)

      allocate(itmp(nenc))
      itmp(:) = pack(index2(:), lencounter(:))
      call move_alloc(itmp, index2)

      allocate(ltmp(nenc))
      ltmp(:) = pack(lvdotr(:), lencounter(:))
      call move_alloc(ltmp, lvdotr)

      return
   end subroutine encounter_check_remove_duplicates


   pure module subroutine encounter_check_sort_aabb_1D(self, n, extent_arr)
      !! author: David A. Minton
      !!
      !! Sorts the bounding box extents along a single dimension prior to the sweep phase. 
      !! This subroutine sets the sorted index array (ind) and the beginning/ending index list (beg & end)
      implicit none
      ! Arguments
      class(encounter_bounding_box_1D), intent(inout) :: self       !! Bounding box structure along a single dimension
      integer(I4B),                     intent(in)    :: n          !! Number of bodies with extents
      real(DP), dimension(:),           intent(in)    :: extent_arr !! Array of extents of size 2*n
      ! Internals
      integer(I8B) :: i, k

      call util_sort(extent_arr, self%ind)

#ifdef DOCONLOC
      do concurrent(k = 1_I8B:2_I8B * n) shared(self,n) local(i)
#else
      do concurrent(k = 1_I8B:2_I8B * n)
#endif
         i = self%ind(k)
         if (i <= n) then
            self%ibeg(i) = k
         else
            self%iend(i - n) = k
         end if
      end do

      return
   end subroutine encounter_check_sort_aabb_1D 


   module subroutine encounter_check_sweep_aabb_double_list(self, n1, n2, r1, v1, r2, v2, renc1, renc2, dt, &
                                                            nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Sweeps the sorted bounding box extents and returns the true encounters (combines broad and narrow phases)
      !! Double list version (e.g. pl-tp or plm-plt)
      implicit none
      ! Arguments
      class(encounter_bounding_box),           intent(inout) :: self       !! Multi-dimensional bounding box structure
      integer(I4B),                            intent(in)    :: n1         !! Number of bodies 1
      integer(I4B),                            intent(in)    :: n2         !! Number of bodies 2
      real(DP),     dimension(:,:),            intent(in)    :: r1, v1     !! Array of position and velocity vectorrs for bodies 1
      real(DP),     dimension(:,:),            intent(in)    :: r2, v2     !! Array of position and velocity vectorrs for bodies 2
      real(DP),     dimension(:),              intent(in)    :: renc1      !! Radius of encounter regions of bodies 1
      real(DP),     dimension(:),              intent(in)    :: renc2      !! Radius of encounter regions of bodies 2
      real(DP),                                intent(in)    :: dt         !! Step size
      integer(I8B),                            intent(out)   :: nenc       !! Total number of encounter candidates
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1     !! List of indices for body 1 in each encounter candidate pair
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2     !! List of indices for body 2 in each encounter candidate pair
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr     !! Logical array indicating which pairs are approaching
      ! Internals
      integer(I4B) :: ii, i, ntot, nbox, dim
      logical, dimension(n1+n2) :: loverlap
      logical, dimension(2*(n1+n2)) :: llist1
      integer(I4B), dimension(2*(n1+n2)) :: ext_ind
      type(collision_list_pltp), dimension(n1+n2) :: lenc         !! Array of encounter lists (one encounter list per body)
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      integer(I8B) :: ibeg, iend
      real(DP), dimension(2*(n1+n2)) :: xind, yind, zind, vxind, vyind, vzind, rencind

      ntot = n1 + n2
      call swiftest_util_index_array(ind_arr, ntot)

      loverlap(:) = (self%aabb%ibeg(:) + 1_I8B) < (self%aabb%iend(:) - 1_I8B)
      where(self%aabb%ind(:) > ntot)
         ext_ind(:) = self%aabb%ind(:) - ntot
      elsewhere
         ext_ind(:) = self%aabb%ind(:)
      endwhere
      llist1(:) = ext_ind(:) <= n1
      where(.not.llist1(:)) ext_ind(:) = ext_ind(:) - n1


      dim = 1
      where(llist1(:))
         xind(:) = r1(1,ext_ind(:))
         yind(:) = r1(2,ext_ind(:))
         zind(:) = r1(3,ext_ind(:))
         vxind(:) = v1(1,ext_ind(:))
         vyind(:) = v1(2,ext_ind(:))
         vzind(:) = v1(3,ext_ind(:))
         rencind(:) = renc1(ext_ind(:))
      elsewhere
         xind(:) = r2(1,ext_ind(:))
         yind(:) = r2(2,ext_ind(:))
         zind(:) = r2(3,ext_ind(:))
         vxind(:) = v2(1,ext_ind(:))
         vyind(:) = v2(2,ext_ind(:))
         vzind(:) = v2(3,ext_ind(:))
         rencind(:) = renc2(ext_ind(:))
      endwhere

      where(.not.loverlap(:)) lenc(:)%nenc = 0
      !$omp parallel default(private) &
      !$omp shared(self, ext_ind, lenc, loverlap, r1, v1, r2, v2, renc1, renc2, xind, yind, zind, vxind, vyind, vzind, rencind, llist1) &
      !$omp firstprivate(ntot, n1, n2, dt, dim) 
     
      ! Do the first group of bodies (i is in list 1, all the others are from list 2)
      !$omp do schedule(static)
      do i = 1, n1
         if (loverlap(i)) then
            ibeg =  self%aabb%ibeg(i) + 1_I8B
            iend =  self%aabb%iend(i) - 1_I8B
            nbox = int(iend - ibeg, kind=I4B) + 1
            call encounter_check_all_sweep_one(i, nbox, r1(1,i), r1(2,i), r1(3,i), v1(1,i), v1(2,i), v1(3,i), &
                                                         xind(ibeg:iend), yind(ibeg:iend), zind(ibeg:iend),&
                                                         vxind(ibeg:iend), vyind(ibeg:iend), vzind(ibeg:iend), &
                                                         renc1(i), rencind(ibeg:iend), dt, ext_ind(ibeg:iend), &
                                                         .not.llist1(ibeg:iend), lenc(i)%nenc, lenc(i)%index1, lenc(i)%index2, lenc(i)%lvdotr)
         end if
      end do
      !$omp end do nowait

      ! Do the second group of bodies (i is in list 2, all the others are in list 1)
      !$omp do schedule(static)
      do i = n1+1, ntot
         if (loverlap(i)) then
            ibeg =  self%aabb%ibeg(i) + 1_I8B
            iend =  self%aabb%iend(i) - 1_I8B
            nbox = int(iend - ibeg, kind=I4B) + 1
            ii = i - n1
            call encounter_check_all_sweep_one(ii, nbox, r2(1,ii), r2(2,ii), r2(3,ii), v2(1,ii), v2(2,ii), v2(3,ii), &
                                                          xind(ibeg:iend), yind(ibeg:iend), zind(ibeg:iend),&
                                                          vxind(ibeg:iend), vyind(ibeg:iend), vzind(ibeg:iend), &
                                                          renc2(ii), rencind(ibeg:iend), dt, ext_ind(ibeg:iend), &
                                                          llist1(ibeg:iend), lenc(i)%nenc, lenc(i)%index2, lenc(i)%index1, lenc(i)%lvdotr)
         end if
      end do
      !$omp end do nowait

      !$omp end parallel

      call encounter_check_collapse_ragged_list(lenc, ntot, nenc, index1, index2, lvdotr)

      call encounter_check_remove_duplicates(ntot, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_sweep_aabb_double_list


   module subroutine encounter_check_sweep_aabb_single_list(self, n, r, v, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Sweeps the sorted bounding box extents and returns the true encounters (combines broad and narrow phases)
      !! Single list version (e.g. pl-pl)
      implicit none
      ! Arguments
      class(encounter_bounding_box),           intent(inout) :: self       !! Multi-dimensional bounding box structure
      integer(I4B),                            intent(in)    :: n          !! Number of bodies
      real(DP),     dimension(:,:),            intent(in)    :: r, v       !! Array of position and velocity vectors 
      real(DP),     dimension(:),              intent(in)    :: renc       !! Radius of encounter regions of bodies 1
      real(DP),                                intent(in)    :: dt         !! Step size
      integer(I8B),                            intent(out)   :: nenc       !! Total number of encounter candidates
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1     !! List of indices for one body in each encounter candidate pair
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2     !! List of indices for the other body in each encounter candidate pair
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr     !! Logical array indicating which pairs are approaching
      ! Internals
      integer(I4B) :: i, dim, itmp, nbox
      integer(I8B) :: k
      logical, dimension(n) :: loverlap
      logical, dimension(2*n) :: lencounteri
      real(DP), dimension(2*n) :: xind, yind, zind, vxind, vyind, vzind, rencind
      integer(I4B), dimension(2*n) :: ext_ind
      type(collision_list_plpl), dimension(n) :: lenc         !! Array of encounter lists (one encounter list per body)
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      integer(I8B) :: ibeg, iend

      call swiftest_util_index_array(ind_arr, n)
      dim = 1

      ! Sweep the intervals for each of the massive bodies along one dimension
      ! This will build a ragged pair of index lists inside of the lenc data structure
      where(self%aabb%ind(:) > n)
         ext_ind(:) = self%aabb%ind(:) - n
      elsewhere
         ext_ind(:) = self%aabb%ind(:)
      endwhere

      xind(:) = r(1,ext_ind(:))
      yind(:) = r(2,ext_ind(:))
      zind(:) = r(3,ext_ind(:))
      vxind(:) = v(1,ext_ind(:))
      vyind(:) = v(2,ext_ind(:))
      vzind(:) = v(3,ext_ind(:))
      rencind(:) = renc(ext_ind(:))

      loverlap(:) = (self%aabb%ibeg(:) + 1_I8B) < (self%aabb%iend(:) - 1_I8B)
      where(.not.loverlap(:)) lenc(:)%nenc = 0

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(self, ext_ind, lenc, loverlap, r, v, renc, xind, yind, zind, vxind, vyind, vzind, rencind) &
      !$omp firstprivate(n, dt, dim) 
      do i = 1, n
         if (loverlap(i)) then
            ibeg =  self%aabb%ibeg(i) + 1_I8B
            iend =  self%aabb%iend(i) - 1_I8B
            nbox = int(iend - ibeg, kind=I4B) + 1
            lencounteri(ibeg:iend) = .true.
            call encounter_check_all_sweep_one(i, nbox, r(1,i), r(2,i), r(3,i), v(1,i), v(2,i), v(3,i), &
                                                      xind(ibeg:iend), yind(ibeg:iend), zind(ibeg:iend),&
                                                      vxind(ibeg:iend), vyind(ibeg:iend), vzind(ibeg:iend), &
                                                      renc(i), rencind(ibeg:iend), dt, ext_ind(ibeg:iend), &
                                                      lencounteri(ibeg:iend), lenc(i)%nenc, lenc(i)%index1, lenc(i)%index2, lenc(i)%lvdotr)
            end if
      end do
      !$omp end parallel do

      call encounter_check_collapse_ragged_list(lenc, n, nenc, index1, index2, lvdotr)

      ! By convention, we always assume that index1 < index2, and so we must swap any that are out of order
#ifdef DOCONLOC
      do concurrent(k = 1_I8B:nenc, index1(k) > index2(k)) shared(index1,index2) local(itmp)
#else
      do concurrent(k = 1_I8B:nenc, index1(k) > index2(k))
#endif
         itmp = index1(k)
         index1(k) = index2(k)
         index2(k) = itmp
      end do

      call encounter_check_remove_duplicates(n, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_sweep_aabb_single_list

end submodule s_encounter_check
