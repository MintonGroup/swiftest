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

   module subroutine encounter_check_all_plpl(param, npl, x, v, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(base_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: x      !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: v      !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter 
      real(DP),                                intent(in)    :: dt     !! Step size
      integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      !type(interaction_timer), save :: itimer
      logical, save :: lfirst = .true.
      logical, save :: skipit = .false. ! This will be used to ensure that the sort & sweep subroutine gets called at least once before timing it so that the extent array is nearly sorted when it is timed
      integer(I8B) :: nplpl = 0_I8B

      ! if (param%ladaptive_encounters_plpl .and. (.not. skipit)) then
      !    nplpl = (npl * (npl - 1) / 2) 
      !    if (nplpl > 0) then
      !       if (lfirst) then  
      !          write(itimer%loopname, *) "encounter_check_all_plpl"
      !          write(itimer%looptype, *) "ENCOUNTER_PLPL"
      !          lfirst = .false.
      !          itimer%step_counter = INTERACTION_TIMER_CADENCE
      !       else 
      !          if (itimer%netcdf_io_check(param, nplpl)) call itimer%time_this_loop(param, nplpl)
      !       end if
      !    else
      !       param%lencounter_sas_plpl = .false.
      !    end if
      ! end if
      
      ! if (param%lencounter_sas_plpl) then
      !    call encounter_check_all_sort_and_sweep_plpl(npl, x, v, renc, dt, nenc, index1, index2, lvdotr) 
      ! else
         call encounter_check_all_triangular_plpl(npl, x, v, renc, dt, nenc, index1, index2, lvdotr) 
      ! end if

      ! if (skipit) then
      !    skipit = .false.
      ! else
      !    if (param%ladaptive_encounters_plpl .and. nplpl > 0) then 
      !       if (itimer%is_on) then
      !          call itimer%adapt(param, nplpl)
      !          skipit = .true.
      !       end if
      !    end if
      ! end if

      return
   end subroutine encounter_check_all_plpl


   module subroutine encounter_check_all_plplm(param, nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between fully interacting massive bodies partially interacting massive bodies. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(base_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: nplm   !! Total number of fully interacting massive bodies 
      integer(I4B),                            intent(in)    :: nplt   !! Total number of partially interacting masive bodies (GM < GMTINY) 
      real(DP),     dimension(:,:),            intent(in)    :: xplm   !! Position vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vplm   !! Velocity vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: xplt   !! Position vectors of partially interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vplt   !! Velocity vectors of partially interacting massive bodies
      real(DP),     dimension(:),              intent(in)    :: rencm  !! Critical radii of fully interacting massive bodies that defines an encounter
      real(DP),     dimension(:),              intent(in)    :: renct  !! Critical radii of partially interacting massive bodies that defines an encounter
      real(DP),                                intent(in)    :: dt     !! Step size
      integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      ! type(interaction_timer), save :: itimer
      logical, save :: lfirst = .true.
      logical, save :: skipit = .false.
      integer(I8B) :: nplplm = 0_I8B
      integer(I4B) :: npl
      logical,      dimension(:), allocatable :: plmplt_lvdotr !! Logical flag indicating the sign of v .dot. x in the plm-plt group
      integer(I4B), dimension(:), allocatable :: plmplt_index1 !! List of indices for body 1 in each encounter in the plm-plt group
      integer(I4B), dimension(:), allocatable :: plmplt_index2 !! List of indices for body 2 in each encounter in the plm-lt group
      integer(I8B)                            :: plmplt_nenc   !! Number of encounters of the plm-plt group
      class(base_parameters), allocatable :: tmp_param     !! Temporary parameter structure to turn off adaptive timer for the pl-pl phase if necessary
      integer(I8B), dimension(:), allocatable :: ind
      integer(I4B), dimension(:), allocatable :: itmp
      logical, dimension(:), allocatable :: ltmp

      ! if (param%ladaptive_encounters_plpl .and. (.not. skipit)) then
      !    npl = nplm + nplt
      !    nplplm = nplm * npl - nplm * (nplm + 1) / 2 
      !    if (nplplm > 0) then
      !       if (lfirst) then  
      !          write(itimer%loopname, *) "encounter_check_all_plpl"
      !          write(itimer%looptype, *) "ENCOUNTER_PLPL"
      !          lfirst = .false.
      !          itimer%step_counter = INTERACTION_TIMER_CADENCE
      !       else 
      !          if (itimer%netcdf_io_check(param, nplplm)) call itimer%time_this_loop(param, nplplm)
      !       end if
      !    else
      !       param%lencounter_sas_plpl = .false.
      !    end if
      ! end if

      allocate(tmp_param, source=param)
      select type(tmp_param)
      class is (swiftest_parameters)
         tmp_param%system_history%nc%lfile_is_open = .false.
      end select

      ! Turn off adaptive encounter checks for the pl-pl group
      tmp_param%ladaptive_encounters_plpl = .false.

      ! Start with the pl-pl group
      call encounter_check_all_plpl(tmp_param, nplm, xplm, vplm, rencm, dt, nenc, index1, index2, lvdotr)

      ! if (param%lencounter_sas_plpl) then
         ! call encounter_check_all_sort_and_sweep_plplm(nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, &
          !                                             plmplt_nenc, plmplt_index1, plmplt_index2, plmplt_lvdotr)
      ! else
         call encounter_check_all_triangular_plplm(nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, plmplt_nenc, plmplt_index1, plmplt_index2, plmplt_lvdotr) 
      ! end if

      ! if (skipit) then
      !    skipit = .false.
      ! else
      !    if (param%ladaptive_encounters_plpl .and. nplplm > 0) then 
      !       if (itimer%is_on) then
      !          call itimer%adapt(param, nplplm)
      !          skipit = .true.
      !       end if
      !    end if
      ! end if

      if (plmplt_nenc > 0) then ! Consolidate the two lists
         allocate(itmp(nenc+plmplt_nenc))
         itmp(1:nenc) = index1(1:nenc)
         itmp(nenc+1:nenc+plmplt_nenc) = plmplt_index1(1:plmplt_nenc)
         call move_alloc(itmp, index1)
         allocate(itmp(nenc+plmplt_nenc))
         itmp(1:nenc) = index2(1:nenc)
         itmp(nenc+1:nenc+plmplt_nenc) = plmplt_index2(1:plmplt_nenc) + nplm ! Be sure to shift these indices back to their natural range
         call move_alloc(itmp, index2)
         allocate(ltmp(nenc+plmplt_nenc))
         ltmp(1:nenc) = lvdotr(1:nenc)
         ltmp(nenc+1:nenc+plmplt_nenc) = plmplt_lvdotr(1:plmplt_nenc)
         call move_alloc(ltmp, lvdotr)
         nenc = nenc + plmplt_nenc

         call swiftest_util_sort(index1, ind)
         call swiftest_util_sort_rearrange(index1, ind, nenc)
         call swiftest_util_sort_rearrange(index2, ind, nenc)
         call swiftest_util_sort_rearrange(lvdotr, ind, nenc)

      end if

      return
   end subroutine encounter_check_all_plplm


   module subroutine encounter_check_all_pltp(param, npl, ntp, xpl, vpl, xtp, vtp, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies and test particles. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(base_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies 
      integer(I4B),                            intent(in)    :: ntp    !! Total number of test particles 
      real(DP),     dimension(:,:),            intent(in)    :: xpl    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vpl    !! Velocity vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: xtp    !! Position vectors of test particlse
      real(DP),     dimension(:,:),            intent(in)    :: vtp    !! Velocity vectors of test particles
      real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter
      real(DP),                                intent(in)    :: dt     !! Step size
      integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      ! type(interaction_timer), save :: itimer
      logical, save :: lfirst = .true.
      logical, save :: lsecond = .false.
      integer(I8B) :: npltp = 0_I8B

      ! if (param%ladaptive_encounters_pltp) then
      !    npltp = npl * ntp
      !    if (npltp > 0) then
      !       if (lfirst) then  
      !          write(itimer%loopname, *) "encounter_check_all_pltp"
      !          write(itimer%looptype, *) "ENCOUNTER_PLTP"
      !          lfirst = .false.
      !          lsecond = .true.
      !       else
      !          if (lsecond) then ! This ensures that the encounter check methods are run at least once prior to timing. Sort and sweep improves on the second pass due to the bounding box extents needing to be nearly sorted 
      !             call itimer%time_this_loop(param, npltp)
      !             lsecond = .false.
      !          else if (itimer%netcdf_io_check(param, npltp)) then
      !             lsecond = .true.
      !             itimer%is_on = .false.
      !          end if
      !       end if
      !    else
      !       param%lencounter_sas_pltp = .false.
      !    end if
      ! end if

      ! if (param%lencounter_sas_pltp) then
      !    call encounter_check_all_sort_and_sweep_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, nenc, index1, index2, lvdotr)
      ! else
         call encounter_check_all_triangular_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, nenc, index1, index2, lvdotr) 
      ! end if

      ! if (.not.lfirst .and. param%ladaptive_encounters_pltp .and. npltp > 0) then 
      !    if (itimer%is_on) call itimer%adapt(param, npltp)
      ! end if

      return
   end subroutine encounter_check_all_pltp


   subroutine encounter_check_all_sort_and_sweep_plpl(npl, x, v, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. 
      !! This is the sort and sweep version
      !! References: Adapted from Christer Ericson's _Real Time Collision Detection_
      !!
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: npl    !! Total number of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: x      !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: v      !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)  :: renc   !! Critical radii of massive bodies that defines an encounter 
      real(DP),                                intent(in)  :: dt     !! Step size
      integer(I8B),                            intent(out) :: nenc   !! Total number of encounters
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      integer(I4B) :: dim, n
      integer(I4B), save :: npl_last = 0
      type(encounter_bounding_box), save :: boundingbox
      integer(I2B), dimension(npl) :: vshift_min, vshift_max

      if (npl == 0) return

      ! If this is the first time through, build the index lists
      n = 2 * npl
      if (npl_last /= npl) then
         call boundingbox%setup(npl, npl_last)
         npl_last = npl
      end if

      !$omp parallel do default(private) schedule(static) &
      !$omp shared(x, v, renc, boundingbox) &
      !$omp firstprivate(dt, npl, n)
      do dim = 1, SWEEPDIM
         where(v(dim,1:npl) < 0.0_DP)
            vshift_min(1:npl) = 1
            vshift_max(1:npl) = 0
         elsewhere
            vshift_min(1:npl) = 0
            vshift_max(1:npl) = 1
         end where
         call boundingbox%aabb(dim)%sort(npl, [x(dim,1:npl) - renc(1:npl) + vshift_min(1:npl) * v(dim,1:npl) * dt, &
                                               x(dim,1:npl) + renc(1:npl) + vshift_max(1:npl) * v(dim,1:npl) * dt])
      end do
      !$omp end parallel do 

      call boundingbox%sweep(npl, x, v, renc, dt, nenc, index1, index2, lvdotr) 

      return
   end subroutine encounter_check_all_sort_and_sweep_plpl


   subroutine encounter_check_all_sort_and_sweep_plplm(nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, nenc, index1, index2, lvdotr)
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
      real(DP),     dimension(:,:),            intent(in)  :: xplm   !! Position vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vplm   !! Velocity vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: xplt   !! Position vectors of partially interacting massive bodies
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
      integer(I4B) :: dim, n, ntot
      integer(I4B), save :: ntot_last = 0
      integer(I2B), dimension(nplm) :: vplmshift_min, vplmshift_max
      integer(I2B), dimension(nplt) :: vpltshift_min, vpltshift_max

      ! If this is the first time through, build the index lists
      if ((nplm == 0) .or. (nplt == 0)) return

      ntot = nplm + nplt
      n = 2 * ntot
      if (ntot /= ntot_last) then
         call boundingbox%setup(ntot, ntot_last)
         ntot_last = ntot
      end if

      !$omp parallel do default(private) schedule(static) &
      !$omp shared(xplm, xplt, vplm, vplt, rencm, renct, boundingbox) &
      !$omp firstprivate(dt, nplm, nplt, ntot)
      do dim = 1, SWEEPDIM
         where(vplm(dim,1:nplm) < 0.0_DP)
            vplmshift_min(1:nplm) = 1
            vplmshift_max(1:nplm) = 0
         elsewhere
            vplmshift_min(1:nplm) = 0
            vplmshift_max(1:nplm) = 1
         end where

         where(vplt(dim,1:nplt) < 0.0_DP)
            vpltshift_min(1:nplt) = 1
            vpltshift_max(1:nplt) = 0
         elsewhere
            vpltshift_min(1:nplt) = 0
            vpltshift_max(1:nplt) = 1
         end where

         call boundingbox%aabb(dim)%sort(ntot, [xplm(dim,1:nplm) - rencm(1:nplm) + vplmshift_min(1:nplm) * vplm(dim,1:nplm) * dt, &
                                                xplt(dim,1:nplt) - renct(1:nplt) + vpltshift_min(1:nplt) * vplt(dim,1:nplt) * dt, &
                                                xplm(dim,1:nplm) + rencm(1:nplm) + vplmshift_max(1:nplm) * vplm(dim,1:nplm) * dt, &
                                                xplt(dim,1:nplt) + renct(1:nplt) + vpltshift_max(1:nplt) * vplt(dim,1:nplt) * dt])
      end do
      !$omp end parallel do

      call boundingbox%sweep(nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, nenc, index1, index2, lvdotr) 

      return
   end subroutine encounter_check_all_sort_and_sweep_plplm


   subroutine encounter_check_all_sort_and_sweep_pltp(npl, ntp, xpl, vpl, xtp, vtp, rencpl, dt, nenc, index1, index2, lvdotr)
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
      real(DP),     dimension(:,:),            intent(in)  :: xpl    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vpl    !! Velocity vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: xtp    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vtp    !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)  :: rencpl !! Critical radii of massive bodies that defines an encounter
      real(DP),                                intent(in)  :: dt     !! Step size
      integer(I8B),                            intent(out) :: nenc   !! Total number of encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      ! Internals
      type(encounter_bounding_box), save :: boundingbox
      integer(I4B) :: dim, n, ntot
      integer(I4B), save :: ntot_last = 0
      integer(I2B), dimension(npl) :: vplshift_min, vplshift_max
      integer(I2B), dimension(ntp) :: vtpshift_min, vtpshift_max
      real(DP), dimension(ntp) :: renctp 

      ! If this is the first time through, build the index lists
      if ((ntp == 0) .or. (npl == 0)) return

      ntot = npl + ntp
      n = 2 * ntot
      if (ntot /= ntot_last) then
         call boundingbox%setup(ntot, ntot_last)
         ntot_last = ntot
      end if

      renctp(:) = 0.0_DP
      
      !$omp parallel do default(private) schedule(static) &
      !$omp shared(xpl, xtp, vpl, vtp, rencpl, renctp, boundingbox) &
      !$omp firstprivate(dt, npl, ntp, ntot)
      do dim = 1, SWEEPDIM
         where(vpl(dim,1:npl) < 0.0_DP)
            vplshift_min(1:npl) = 1
            vplshift_max(1:npl) = 0
         elsewhere
            vplshift_min(1:npl) = 0
            vplshift_max(1:npl) = 1
         end where

         where(vtp(dim,1:ntp) < 0.0_DP)
            vtpshift_min(1:ntp) = 1
            vtpshift_max(1:ntp) = 0
         elsewhere
            vtpshift_min(1:ntp) = 0
            vtpshift_max(1:ntp) = 1
         end where

         call boundingbox%aabb(dim)%sort(ntot, [xpl(dim,1:npl) - rencpl(1:npl) + vplshift_min(1:npl) * vpl(dim,1:npl) * dt, &
                                                xtp(dim,1:ntp) - renctp(1:ntp) + vtpshift_min(1:ntp) * vtp(dim,1:ntp) * dt, &
                                                xpl(dim,1:npl) + rencpl(1:npl) + vplshift_max(1:npl) * vpl(dim,1:npl) * dt, &
                                                xtp(dim,1:ntp) + renctp(1:ntp) + vtpshift_max(1:ntp) * vtp(dim,1:ntp) * dt])
      end do
      !$omp end parallel do

      call boundingbox%sweep(npl, ntp, xpl, vpl, xtp, vtp, rencpl, renctp, dt, nenc, index1, index2, lvdotr) 

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
      do concurrent(j = 1:n, lgood(j))
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


   pure subroutine encounter_check_all_triangular_one(i, n, xi, yi, zi, vxi, vyi, vzi, x, y, z, vx, vy, vz, renci, renc, &
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

      do concurrent(j = i+1:n)
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


   subroutine encounter_check_all_triangular_plpl(npl, x, v, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Split off from the main subroutine for performance
      !! This is the upper triangular (double loop) version.
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: npl    !! Total number of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: x      !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: v      !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)  :: renc  !! Critical radii of massive bodies that defines an encounter 
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
      !$omp shared(x, v, renc, lenc, ind_arr) &
      !$omp firstprivate(npl, dt)
      do i = 1,npl
         call encounter_check_all_triangular_one(i, npl, x(1,i), x(2,i), x(3,i), &
                                                         v(1,i), v(2,i), v(3,i), &
                                                         x(1,:), x(2,:), x(3,:), &
                                                         v(1,:), v(2,:), v(3,:), &
                                                         renc(i), renc(:), dt, ind_arr(:), lenc(i))
         if (lenc(i)%nenc > 0) lenc(i)%index1(:) = i
      end do
      !$omp end parallel do

      call encounter_check_collapse_ragged_list(lenc, npl, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_all_triangular_plpl


   subroutine encounter_check_all_triangular_plplm(nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, &
                                                   nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies and test particles. Split off from the main subroutine for performance
      !! This is the upper triangular (double loop) version.
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: nplm   !! Total number of fully interacting massive bodies 
      integer(I4B),                            intent(in)  :: nplt   !! Total number of partially interacting masive bodies (GM < GMTINY) 
      real(DP),     dimension(:,:),            intent(in)  :: xplm   !! Position vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vplm   !! Velocity vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: xplt   !! Position vectors of partially interacting massive bodies
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
      !$omp shared(xplm, vplm, xplt, vplt, rencm, renct, lenc, ind_arr) &
      !$omp firstprivate(nplm, nplt, dt)
      do i = 1, nplm
         call encounter_check_all_triangular_one(0, nplt, xplm(1,i), xplm(2,i), xplm(3,i), &
                                                          vplm(1,i), vplm(2,i), vplm(3,i), &
                                                          xplt(1,:), xplt(2,:), xplt(3,:), &
                                                          vplt(1,:), vplt(2,:), vplt(3,:), &
                                                          rencm(i), renct(:), dt, ind_arr(:), lenc(i))
         if (lenc(i)%nenc > 0) lenc(i)%index1(:) = i
      end do
      !$omp end parallel do

      call encounter_check_collapse_ragged_list(lenc, nplm, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_all_triangular_plplm


   subroutine encounter_check_all_triangular_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, &
                                                  nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies and test particles. Split off from the main subroutine for performance
      !! This is the upper triangular (double loop) version.
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: npl    !! Total number of massive bodies 
      integer(I4B),                            intent(in)  :: ntp    !! Total number of test particles 
      real(DP),     dimension(:,:),            intent(in)  :: xpl    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: vpl    !! Velocity vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: xtp    !! Position vectors of massive bodies
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
      !$omp shared(xpl, vpl, xtp, vtp, renc, renct, lenc, ind_arr) &
      !$omp firstprivate(npl, ntp, dt)
      do i = 1, npl
         call encounter_check_all_triangular_one(0, ntp, xpl(1,i), xpl(2,i), xpl(3,i), &
                                                         vpl(1,i), vpl(2,i), vpl(3,i), &
                                                         xtp(1,:), xtp(2,:), xtp(3,:), &
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
            tmin = -vdotr / v2
         
            if (tmin < dt) then
               r2min = r2 - vdotr**2 / v2
            else
               r2min = r2 + 2 * vdotr * dt + v2 * dt**2
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
      integer(I4B),                                    intent(in)            :: n1          !! Number of bodies 1
      integer(I8B),                                    intent(out)           :: nenc        !! Total number of encountersj 
      integer(I4B),         dimension(:), allocatable, intent(out)           :: index1      !! Array of indices for body 1
      integer(I4B),         dimension(:), allocatable, intent(out)           :: index2      !! Array of indices for body 1
      logical,              dimension(:), allocatable, intent(out), optional :: lvdotr      !! Array indicating which bodies are approaching
      ! Internals
      integer(I4B) :: i
      integer(I8B) :: j1, j0, nenci
      integer(I4B), dimension(n1) :: ibeg

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

      call swiftest_util_sort(index1, ind)
      call swiftest_util_sort_rearrange(index1, ind, nenc)
      call swiftest_util_sort_rearrange(index2, ind, nenc)
      call swiftest_util_sort_rearrange(lvdotr, ind, nenc)

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
      do concurrent(i = 1:n, iend(i) - ibeg(i) > 0_I8B)
         klo = ibeg(i)
         khi = iend(i)
         nenci = khi - klo + 1_I8B
         if (allocated(ind)) deallocate(ind)
         call swiftest_util_sort(index2(klo:khi), ind)
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

      call swiftest_util_sort(extent_arr, self%ind)

      do concurrent(k = 1_I8B:2_I8B * n)
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
      logical, dimension(SWEEPDIM,n1+n2) :: loverlap_by_dimension
      logical, dimension(SWEEPDIM,2*(n1+n2)) :: llist1
      integer(I4B), dimension(SWEEPDIM,2*(n1+n2)) :: ext_ind
      type(collision_list_pltp), dimension(n1+n2) :: lenc         !! Array of encounter lists (one encounter list per body)
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      integer(I8B) :: ibeg, iend
      real(DP), dimension(2*(n1+n2)) :: xind, yind, zind, vxind, vyind, vzind, rencind

      ntot = n1 + n2
      call swiftest_util_index_array(ind_arr, ntot)

      do concurrent(dim = 1:SWEEPDIM)
         loverlap_by_dimension(dim,:) = (self%aabb(dim)%ibeg(:) + 1_I8B) < (self%aabb(dim)%iend(:) - 1_I8B)
         where(self%aabb(dim)%ind(:) > ntot)
            ext_ind(dim,:) = self%aabb(dim)%ind(:) - ntot
         elsewhere
            ext_ind(dim,:) = self%aabb(dim)%ind(:)
         endwhere
      end do
      llist1(:,:) = ext_ind(:,:) <= n1
      where(.not.llist1(:,:)) ext_ind(:,:) = ext_ind(:,:) - n1

      loverlap(:) = loverlap_by_dimension(1,:)
      do dim = 2, SWEEPDIM
         loverlap(:) = loverlap(:) .and. loverlap_by_dimension(dim,:)
      end do

      dim = 1
      where(llist1(dim,:))
         xind(:) = r1(1,ext_ind(dim,:))
         yind(:) = r1(2,ext_ind(dim,:))
         zind(:) = r1(3,ext_ind(dim,:))
         vxind(:) = v1(1,ext_ind(dim,:))
         vyind(:) = v1(2,ext_ind(dim,:))
         vzind(:) = v1(3,ext_ind(dim,:))
         rencind(:) = renc1(ext_ind(dim,:))
      elsewhere
         xind(:) = r2(1,ext_ind(dim,:))
         yind(:) = r2(2,ext_ind(dim,:))
         zind(:) = r2(3,ext_ind(dim,:))
         vxind(:) = v2(1,ext_ind(dim,:))
         vyind(:) = v2(2,ext_ind(dim,:))
         vzind(:) = v2(3,ext_ind(dim,:))
         rencind(:) = renc2(ext_ind(dim,:))
      endwhere

      where(.not.loverlap(:)) lenc(:)%nenc = 0
      !$omp parallel default(private) &
      !$omp shared(self, ext_ind, lenc, loverlap, r1, v1, r2, v2, renc1, renc2, xind, yind, zind, vxind, vyind, vzind, rencind, llist1) &
      !$omp firstprivate(ntot, n1, n2, dt, dim) 
     
      ! Do the first group of bodies (i is in list 1, all the others are from list 2)
      !$omp do schedule(static)
      do i = 1, n1
         if (loverlap(i)) then
            ibeg =  self%aabb(dim)%ibeg(i) + 1_I8B
            iend =  self%aabb(dim)%iend(i) - 1_I8B
            nbox = iend - ibeg + 1
            call encounter_check_all_sweep_one(i, nbox, r1(1,i), r1(2,i), r1(3,i), v1(1,i), v1(2,i), v1(3,i), &
                                                         xind(ibeg:iend), yind(ibeg:iend), zind(ibeg:iend),&
                                                         vxind(ibeg:iend), vyind(ibeg:iend), vzind(ibeg:iend), &
                                                         renc1(i), rencind(ibeg:iend), dt, ext_ind(dim,ibeg:iend), &
                                                         .not.llist1(dim,ibeg:iend), lenc(i)%nenc, lenc(i)%index1, lenc(i)%index2, lenc(i)%lvdotr)
         end if
      end do
      !$omp end do nowait

      ! Do the second group of bodies (i is in list 2, all the others are in list 1)
      !$omp do schedule(static)
      do i = n1+1, ntot
         if (loverlap(i)) then
            ibeg =  self%aabb(dim)%ibeg(i) + 1_I8B
            iend =  self%aabb(dim)%iend(i) - 1_I8B
            nbox = iend - ibeg + 1
            ii = i - n1
            call encounter_check_all_sweep_one(ii, nbox, r1(1,ii), r1(2,ii), r1(3,ii), v2(1,ii), v2(2,ii), v2(3,ii), &
                                                          xind(ibeg:iend), yind(ibeg:iend), zind(ibeg:iend),&
                                                          vxind(ibeg:iend), vyind(ibeg:iend), vzind(ibeg:iend), &
                                                          renc2(ii), rencind(ibeg:iend), dt, ext_ind(dim,ibeg:iend), &
                                                          llist1(dim,ibeg:iend), lenc(i)%nenc, lenc(i)%index2, lenc(i)%index1, lenc(i)%lvdotr)
         end if
      end do
      !$omp end do nowait

      !$omp end parallel

      call encounter_check_collapse_ragged_list(lenc, ntot, nenc, index1, index2, lvdotr)

      call encounter_check_remove_duplicates(ntot, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_sweep_aabb_double_list


   module subroutine encounter_check_sweep_aabb_single_list(self, n, x, v, renc, dt, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!
      !! Sweeps the sorted bounding box extents and returns the true encounters (combines broad and narrow phases)
      !! Single list version (e.g. pl-pl)
      implicit none
      ! Arguments
      class(encounter_bounding_box),           intent(inout) :: self       !! Multi-dimensional bounding box structure
      integer(I4B),                            intent(in)    :: n          !! Number of bodies
      real(DP),     dimension(:,:),            intent(in)    :: x, v       !! Array of position and velocity vectors 
      real(DP),     dimension(:),              intent(in)    :: renc       !! Radius of encounter regions of bodies 1
      real(DP),                                intent(in)    :: dt         !! Step size
      integer(I8B),                            intent(out)   :: nenc       !! Total number of encounter candidates
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1     !! List of indices for one body in each encounter candidate pair
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2     !! List of indices for the other body in each encounter candidate pair
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr     !! Logical array indicating which pairs are approaching
      ! Internals
      integer(I4B) :: i, nbox, dim
      integer(I8B) :: k, itmp
      logical, dimension(n) :: loverlap
      logical, dimension(2*n) :: lencounteri
      real(DP), dimension(2*n) :: xind, yind, zind, vxind, vyind, vzind, rencind
      integer(I4B), dimension(SWEEPDIM,2*n) :: ext_ind
      type(collision_list_plpl), dimension(n) :: lenc         !! Array of encounter lists (one encounter list per body)
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      integer(I8B) :: ibeg, iend

      call swiftest_util_index_array(ind_arr, n)
      dim = 1

      ! Sweep the intervals for each of the massive bodies along one dimension
      ! This will build a ragged pair of index lists inside of the lenc data structure
      where(self%aabb(dim)%ind(:) > n)
         ext_ind(dim,:) = self%aabb(1)%ind(:) - n
      elsewhere
         ext_ind(dim,:) = self%aabb(1)%ind(:)
      endwhere

      xind(:) = x(1,ext_ind(dim,:))
      yind(:) = x(2,ext_ind(dim,:))
      zind(:) = x(3,ext_ind(dim,:))
      vxind(:) = v(1,ext_ind(dim,:))
      vyind(:) = v(2,ext_ind(dim,:))
      vzind(:) = v(3,ext_ind(dim,:))
      rencind(:) = renc(ext_ind(dim,:))

      loverlap(:) = (self%aabb(dim)%ibeg(:) + 1_I8B) < (self%aabb(dim)%iend(:) - 1_I8B)
      where(.not.loverlap(:)) lenc(:)%nenc = 0

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(self, ext_ind, lenc, loverlap, x, v, renc, xind, yind, zind, vxind, vyind, vzind, rencind) &
      !$omp firstprivate(n, dt, dim) 
      do i = 1, n
         if (loverlap(i)) then
            ibeg =  self%aabb(dim)%ibeg(i) + 1_I8B
            iend =  self%aabb(dim)%iend(i) - 1_I8B
            nbox = iend - ibeg + 1
            lencounteri(ibeg:iend) = .true.
            call encounter_check_all_sweep_one(i, nbox, x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i), &
                                                      xind(ibeg:iend), yind(ibeg:iend), zind(ibeg:iend),&
                                                      vxind(ibeg:iend), vyind(ibeg:iend), vzind(ibeg:iend), &
                                                      renc(i), rencind(ibeg:iend), dt, ext_ind(dim,ibeg:iend), &
                                                      lencounteri(ibeg:iend), lenc(i)%nenc, lenc(i)%index1, lenc(i)%index2, lenc(i)%lvdotr)
            end if
      end do
      !$omp end parallel do

      call encounter_check_collapse_ragged_list(lenc, n, nenc, index1, index2, lvdotr)

      ! By convention, we always assume that index1 < index2, and so we must swap any that are out of order
      do concurrent(k = 1_I8B:nenc, index1(k) > index2(k))
         itmp = index1(k)
         index1(k) = index2(k)
         index2(k) = itmp
      end do

      call encounter_check_remove_duplicates(n, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_sweep_aabb_single_list

end submodule s_encounter_check
