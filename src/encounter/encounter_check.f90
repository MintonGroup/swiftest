submodule (encounter_classes) s_encounter_check
   use swiftest
contains

   module subroutine encounter_check_all(nenc, index1, index2, x1, v1, x2, v2, renc1, renc2, dt, lencounter, lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between n pairs of bodies. 
      !! This implementation is general for any type of body. For instance, for massive bodies, you would pass x1 = x2, for test particles renc2 is an array of zeros, etc.
      !!
      implicit none
      ! Arguments
      integer(I4B),                 intent(in)  :: nenc       !! Number of encounters in the encounter lists
      integer(I4B), dimension(:),   intent(in)  :: index1     !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:),   intent(in)  :: index2     !! List of indices for body 2 in each encounter1
      real(DP),     dimension(:,:), intent(in)  :: x1, v1     !! Array of indices of bodies 1
      real(DP),     dimension(:,:), intent(in)  :: x2, v2     !! Array of indices of bodies 2
      real(DP),     dimension(:),   intent(in)  :: renc1      !! Radius of encounter regions of bodies 1
      real(DP),     dimension(:),   intent(in)  :: renc2      !! Radius of encounter regions of bodies 2
      real(DP),                     intent(in)  :: dt         !! Step size
      logical,      dimension(:),   intent(out) :: lencounter !! Logical array indicating which pairs are in a close encounter state
      logical,      dimension(:),   intent(out) :: lvdotr     !! Logical array indicating which pairs are approaching
      ! Internals
      integer(I4B) :: i, j, k
      real(DP)     :: xr, yr, zr, vxr, vyr, vzr, renc12

      !$omp parallel do simd default(firstprivate) schedule(static)&
      !$omp shared(lencounter, lvdotr, index1, index2, x1, v1, x2, v2) &
      !$omp lastprivate(i, j, xr, yr, zr, vxr, vyr, vzr, renc12)
      do k = 1, nenc
         i = index1(k)
         j = index2(k)
         xr = x2(1, j) - x1(1, i)
         yr = x2(2, j) - x1(2, i)
         zr = x2(3, j) - x1(3, i)
         vxr = v2(1, j) - v1(1, i)
         vyr = v2(2, j) - v1(2, i)
         vzr = v2(3, j) - v1(3, i)
         renc12 = renc1(i) + renc2(j)
         call encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc12, dt, lencounter(k), lvdotr(k)) 
      end do
      !$omp end parallel do simd

      return
   end subroutine encounter_check_all


   module subroutine encounter_check_all_plpl(param, npl, x, v, renc, dt, lvdotr, index1, index2, nenc)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(swiftest_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: x      !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: v      !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter 
      real(DP),                                intent(in)    :: dt     !! Step size
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out)   :: nenc   !! Total number of encounters
      ! Internals
      type(interaction_timer), save :: itimer
      logical, save :: lfirst = .true.
      logical, save :: skipit = .false. ! This will be used to ensure that the sort & sweep subroutine gets called at least once before timing it so that the extent array is nearly sorted when it is timed
      integer(I8B) :: nplpl = 0_I8B

      if (param%ladaptive_encounters_plpl .and. (.not. skipit)) then
         nplpl = (npl * (npl - 1) / 2) 
         if (nplpl > 0) then
            if (lfirst) then  
               write(itimer%loopname, *) "encounter_check_all_plpl"
               write(itimer%looptype, *) "ENCOUNTER_PLPL"
               lfirst = .false.
               itimer%step_counter = INTERACTION_TIMER_CADENCE
            else 
               if (itimer%check(param, nplpl)) call itimer%time_this_loop(param, nplpl)
            end if
         else
            param%lencounter_sas_plpl = .false.
         end if
      end if

      if (param%lencounter_sas_plpl) then
         call encounter_check_all_sort_and_sweep_plpl(npl, x, v, renc, dt, lvdotr, index1, index2, nenc) 
      else
         call encounter_check_all_triangular_plpl(npl, x, v, renc, dt, lvdotr, index1, index2, nenc) 
      end if

      if (skipit) then
         skipit = .false.
      else
         if (param%ladaptive_encounters_plpl .and. nplpl > 0) then 
            if (itimer%is_on) then
               call itimer%adapt(param, nplpl)
               skipit = .true.
            end if
         end if
      end if

      return
   end subroutine encounter_check_all_plpl


   module subroutine encounter_check_all_plplm(param, nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, lvdotr, index1, index2, nenc)
      !! author: David A. Minton
      !!
      !! Check for encounters between fully interacting massive bodies partially interacting massive bodies. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(swiftest_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: nplm   !! Total number of fully interacting massive bodies 
      integer(I4B),                            intent(in)    :: nplt   !! Total number of partially interacting masive bodies (GM < GMTINY) 
      real(DP),     dimension(:,:),            intent(in)    :: xplm   !! Position vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vplm   !! Velocity vectors of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: xplt   !! Position vectors of partially interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vplt   !! Velocity vectors of partially interacting massive bodies
      real(DP),     dimension(:),              intent(in)    :: rencm  !! Critical radii of fully interacting massive bodies that defines an encounter
      real(DP),     dimension(:),              intent(in)    :: renct  !! Critical radii of partially interacting massive bodies that defines an encounter
      real(DP),                                intent(in)    :: dt     !! Step size
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out)   :: nenc   !! Total number of encounters
      ! Internals
      type(interaction_timer), save :: itimer
      logical, save :: lfirst = .true.
      logical, save :: skipit = .false.
      integer(I8B) :: nplplm = 0_I8B
      integer(I4B) :: npl, i
      logical,      dimension(:), allocatable :: plmplt_lvdotr !! Logical flag indicating the sign of v .dot. x in the plm-plt group
      integer(I4B), dimension(:), allocatable :: plmplt_index1 !! List of indices for body 1 in each encounter in the plm-plt group
      integer(I4B), dimension(:), allocatable :: plmplt_index2 !! List of indices for body 2 in each encounter in the plm-lt group
      integer(I4B)                            :: plmplt_nenc   !! Number of encounters of the plm-plt group
      class(swiftest_parameters), allocatable :: tmp_param     !! Temporary parameter structure to turn off adaptive timer for the pl-pl phase if necessary
      integer(I4B), dimension(:), allocatable :: itmp, ind
      logical, dimension(:), allocatable :: ltmp

      if (param%ladaptive_encounters_plpl .and. (.not. skipit)) then
         npl = nplm + nplt
         nplplm = nplm * npl - nplm * (nplm + 1) / 2 
         if (nplplm > 0) then
            if (lfirst) then  
               write(itimer%loopname, *) "encounter_check_all_plpl"
               write(itimer%looptype, *) "ENCOUNTER_PLPL"
               lfirst = .false.
               itimer%step_counter = INTERACTION_TIMER_CADENCE
            else 
               if (itimer%check(param, nplplm)) call itimer%time_this_loop(param, nplplm)
            end if
         else
            param%lencounter_sas_plpl = .false.
         end if
      end if

      allocate(tmp_param, source=param)

      ! Turn off adaptive encounter checks for the pl-pl group
      tmp_param%ladaptive_encounters_plpl = .false.

      ! Start with the pl-pl group
      call encounter_check_all_plpl(tmp_param, nplm, xplm, vplm, rencm, dt, lvdotr, index1, index2, nenc)

      if (param%lencounter_sas_plpl) then
         call encounter_check_all_sort_and_sweep_plplm(nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, plmplt_lvdotr, plmplt_index1, plmplt_index2, plmplt_nenc)
      else
         call encounter_check_all_triangular_plplm(nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, plmplt_lvdotr, plmplt_index1, plmplt_index2, plmplt_nenc) 
      end if

      if (skipit) then
         skipit = .false.
      else
         if (param%ladaptive_encounters_plpl .and. nplplm > 0) then 
            if (itimer%is_on) then
               call itimer%adapt(param, nplplm)
               skipit = .true.
            end if
         end if
      end if

      if (plmplt_nenc > 0) then ! Consolidate the two lists
         allocate(itmp(nenc+plmplt_nenc))
         itmp(1:nenc) = index1(1:nenc)
         itmp(nenc+1:nenc+plmplt_nenc) = plmplt_index1(1:plmplt_nenc)
         call move_alloc(itmp, index1)
         allocate(itmp(nenc+plmplt_nenc))
         itmp(1:nenc) = index2(1:nenc)
         itmp(nenc+1:nenc+plmplt_nenc) = plmplt_index2(1:plmplt_nenc)
         call move_alloc(itmp, index2)
         allocate(ltmp(nenc+plmplt_nenc))
         ltmp(1:nenc) = lvdotr(1:nenc)
         ltmp(nenc+1:nenc+plmplt_nenc) = plmplt_lvdotr(1:plmplt_nenc)
         call move_alloc(ltmp, lvdotr)
         nenc = nenc + plmplt_nenc

         call util_sort(index1, itmp)
         call util_sort_rearrange(index1, itmp, nenc)
         call util_sort_rearrange(index2, itmp, nenc)
         call util_sort_rearrange(lvdotr, itmp, nenc)
      end if

      return
   end subroutine encounter_check_all_plplm


   module subroutine encounter_check_all_pltp(param, npl, ntp, xpl, vpl, xtp, vtp, renc, dt, lvdotr, index1, index2, nenc)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies and test particles. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(swiftest_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies 
      integer(I4B),                            intent(in)    :: ntp    !! Total number of test particles 
      real(DP),     dimension(:,:),            intent(in)    :: xpl    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vpl    !! Velocity vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: xtp    !! Position vectors of test particlse
      real(DP),     dimension(:,:),            intent(in)    :: vtp    !! Velocity vectors of test particles
      real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter
      real(DP),                                intent(in)    :: dt     !! Step size
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out)   :: nenc   !! Total number of encounters
      ! Internals
      type(interaction_timer), save :: itimer
      logical, save :: lfirst = .true.
      logical, save :: lsecond = .false.
      integer(I8B) :: npltp = 0_I8B

      if (param%ladaptive_encounters_pltp) then
         npltp = npl * ntp
         if (npltp > 0) then
            if (lfirst) then  
               write(itimer%loopname, *) "encounter_check_all_pltp"
               write(itimer%looptype, *) "ENCOUNTER_PLTP"
               lfirst = .false.
               lsecond = .true.
            else
               if (lsecond) then ! This ensures that the encounter check methods are run at least once prior to timing. Sort and sweep improves on the second pass due to the bounding box extents needing to be nearly sorted 
                  call itimer%time_this_loop(param, npltp)
                  lsecond = .false.
               else if (itimer%check(param, npltp)) then
                  lsecond = .true.
                  itimer%is_on = .false.
               end if
            end if
         else
            param%lencounter_sas_pltp = .false.
         end if
      end if

      if (param%lencounter_sas_pltp) then
         call encounter_check_all_sort_and_sweep_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, lvdotr, index1, index2, nenc)
      else
         call encounter_check_all_triangular_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, lvdotr, index1, index2, nenc) 
      end if

      if (.not.lfirst .and. param%ladaptive_encounters_pltp .and. npltp > 0) then 
         if (itimer%is_on) call itimer%adapt(param, npltp)
      end if

      return
   end subroutine encounter_check_all_pltp


   subroutine encounter_check_reduce_broadphase(n, nenc, index1, index2, lencounter, lvdotr)
      !! author: David A. Minton
      !!
      !! Takes the candidate encounter lists that came out of the broad phase and narrow it down to the true encounters
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)    :: n          !! Number of bodies 
      integer(I4B),                            intent(inout) :: nenc       !! Number of encountering bodies (input is the broad phase value, output is the final narrow phase value)
      integer(I4B), dimension(:), allocatable, intent(inout) :: index1     !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(inout) :: index2     !! List of indices for body 2 in each encounter
      logical,      dimension(:), allocatable, intent(inout) :: lencounter !! Logical flag indicating which of the pairs identified in the broad phase was selected in the narrow phase
      logical,      dimension(:), allocatable, intent(inout) :: lvdotr     !! Logical flag indicating the sign of v .dot. x
      ! Internals
      integer(I4B) :: i, i0, j, k
      integer(I4B), dimension(:), allocatable :: itmp
      logical, dimension(n) :: lfresh
      integer(I4B), dimension(:), allocatable :: ind
      logical, dimension(:), allocatable :: ltmp

      nenc = count(lencounter(:)) ! Count the true number of encounters

      allocate(itmp(nenc))
      itmp(:) = pack(index1(:), lencounter(:))
      call move_alloc(itmp, index1)

      allocate(itmp(nenc))
      itmp(:) = pack(index2(:), lencounter(:))
      call move_alloc(itmp, index2)

      allocate(ltmp(nenc))
      ltmp(:) = pack(lvdotr(:), lencounter(:))
      call move_alloc(ltmp, lvdotr)

      if (allocated(lencounter)) deallocate(lencounter)
      allocate(lencounter(nenc))
      lencounter(:) = .true.

      call util_sort(index1, ind)
      lfresh(:) = .true.

      ! Remove duplicates
      do k = 1, nenc 
         i = index1(ind(k))
         j = index2(ind(k))
         if (k == 1) then
            lfresh(j) = .false.
         else
            i0 = index1(ind(k-1))
            if (i /= i0) lfresh(:) = .true.
            if (.not.lfresh(j)) lencounter(ind(k)) = .false.
            lfresh(j) = .false.
         end if
      end do

      if (count(lencounter(:)) == nenc) return
      nenc = count(lencounter(:)) ! Count the true number of encounters

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
   end subroutine encounter_check_reduce_broadphase


   subroutine encounter_check_all_sort_and_sweep_plpl(npl, x, v, renc, dt, lvdotr, index1, index2, nenc)
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
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out) :: nenc   !! Total number of encounters
      ! Internals
      integer(I4B) :: i, dim, n
      integer(I4B), save :: npl_last = 0
      type(encounter_bounding_box), save :: boundingbox
      logical, dimension(:), allocatable :: lencounter
      integer(I2B), dimension(npl) :: vshift_min, vshift_max
      type(walltimer) :: timer 

      if (npl == 0) return

      ! If this is the first time through, build the index lists
      n = 2 * npl
      if (npl_last /= npl) then
         call boundingbox%setup(npl, npl_last)
         npl_last = npl
      end if

      call timer%reset()
      call timer%start()
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
      call timer%stop()
      write(*,*) "plpl  sort  : ",timer%count_stop_step - timer%count_start_step

      call timer%reset()
      call timer%start()
      call boundingbox%sweep(npl, nenc, index1, index2)
      call timer%stop()
      write(*,*) "plpl  sweep : ",timer%count_stop_step - timer%count_start_step

      if (nenc > 0) then
         ! Now that we have identified potential pairs, use the narrow-phase process to get the final values
         allocate(lencounter(nenc))
         allocate(lvdotr(nenc))

         call timer%reset()
         call timer%start()
         call encounter_check_all(nenc, index1, index2, x, v, x, v, renc, renc, dt, lencounter, lvdotr)
         call timer%stop()
         write(*,*) "plpl  check : ",timer%count_stop_step - timer%count_start_step

         call timer%reset()
         call timer%start()
         call encounter_check_reduce_broadphase(npl, nenc, index1, index2, lencounter, lvdotr)
         call timer%stop()
         write(*,*) "plpl  reduce: ",timer%count_stop_step - timer%count_start_step
         deallocate(lencounter)
      end if

      return
   end subroutine encounter_check_all_sort_and_sweep_plpl


   subroutine encounter_check_all_sort_and_sweep_plplm(nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, lvdotr, index1, index2, nenc)
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
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out) :: nenc   !! Total number of encounter
      ! Internals
      type(encounter_bounding_box), save :: boundingbox
      integer(I4B) :: i, dim, n, ntot
      integer(I4B), save :: ntot_last = 0
      logical, dimension(:), allocatable :: lencounter
      integer(I2B), dimension(nplm) :: vplmshift_min, vplmshift_max
      integer(I2B), dimension(nplt) :: vpltshift_min, vpltshift_max
      type(walltimer) :: timer 

      ! If this is the first time through, build the index lists
      if ((nplm == 0) .or. (nplt == 0)) return

      ntot = nplm + nplt
      n = 2 * ntot
      if (ntot /= ntot_last) then
         call boundingbox%setup(ntot, ntot_last)
         ntot_last = ntot
      end if
    
      call timer%reset()
      call timer%start()
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
      call timer%stop()
      write(*,*) "plplm sort  : ",timer%count_stop_step - timer%count_start_step

      call timer%reset()
      call timer%start()
      call boundingbox%sweep(nplm, nplt, nenc, index1, index2)
      call timer%stop()
      write(*,*) "plplm sweep : ",timer%count_stop_step - timer%count_start_step

      if (nenc > 0) then
         ! Shift tiny body indices back into the range of the input position and velocity arrays
         index2(:) = index2(:) - nplm

         ! Now that we have identified potential pairs, use the narrow-phase process to get the final values
         allocate(lencounter(nenc))
         allocate(lvdotr(nenc))

         call timer%reset()
         call timer%start()
         call encounter_check_all(nenc, index1, index2, xplm, vplm, xplt, vplt, rencm, renct, dt, lencounter, lvdotr)
         call timer%stop()
         write(*,*) "plplm check : ",timer%count_stop_step - timer%count_start_step

         ! Shift the tiny body indices back to their natural range
         index2(:) = index2(:) + nplm

         call timer%reset()
         call timer%start()
         call encounter_check_reduce_broadphase(ntot, nenc, index1, index2, lencounter, lvdotr)
         call timer%stop()
         write(*,*) "plplm reduce: ",timer%count_stop_step - timer%count_start_step
      end if
      return
   end subroutine encounter_check_all_sort_and_sweep_plplm


   subroutine encounter_check_all_sort_and_sweep_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, lvdotr, index1, index2, nenc)
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
      real(DP),     dimension(:),              intent(in)  :: renc   !! Critical radii of massive bodies that defines an encounter
      real(DP),                                intent(in)  :: dt     !! Step size
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out) :: nenc   !! Total number of encounter
      ! Internals
      type(encounter_bounding_box), save :: boundingbox
      integer(I4B) :: i, dim, n, ntot
      integer(I4B), save :: ntot_last = 0
      logical, dimension(:), allocatable :: lencounter
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
      
      !$omp parallel do default(private) schedule(static) &
      !$omp shared(xpl, xtp, vpl, vtp, renc, boundingbox) &
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

         call boundingbox%aabb(dim)%sort(ntot, [xpl(dim,1:npl) - renc(1:npl) + vplshift_min(1:npl) * vpl(dim,1:npl) * dt, &
                                                xtp(dim,1:ntp)               + vtpshift_min(1:ntp) * vtp(dim,1:ntp) * dt, &
                                                xpl(dim,1:npl) + renc(1:npl) + vplshift_max(1:npl) * vpl(dim,1:npl) * dt, &
                                                xtp(dim,1:ntp)               + vtpshift_max(1:ntp) * vtp(dim,1:ntp) * dt])
      end do
      !$omp end parallel do

      call boundingbox%sweep(npl, ntp, nenc, index1, index2)

      if (nenc > 0) then
         ! Shift test particle indices back into the proper range
         index2(:) = index2(:) - npl

         ! Now that we have identified potential pairs, use the narrow-phase process to get the final values
         allocate(lencounter(nenc))
         allocate(lvdotr(nenc))
         renctp(:) = 0.0_DP

         call encounter_check_all(nenc, index1, index2, xpl, vpl, xtp, vtp, renc, renctp, dt, lencounter, lvdotr)

         call encounter_check_reduce_broadphase(ntot, nenc, index1, index2, lencounter, lvdotr)

      end if

      return
   end subroutine encounter_check_all_sort_and_sweep_pltp


   pure subroutine encounter_check_all_triangular_one(i, n, xi, yi, zi, vxi, vyi, vzi, x, y, z, vx, vy, vz, renci, renc, dt, ind_arr, lenci)
      implicit none
      ! Arguments
      integer(I4B),                       intent(in)  :: i
      integer(I4B),                       intent(in)  :: n
      real(DP),                           intent(in)  :: xi, yi, zi
      real(DP),                           intent(in)  :: vxi, vyi, vzi
      real(DP),             dimension(:), intent(in)  :: x, y, z
      real(DP),             dimension(:), intent(in)  :: vx, vy, vz
      real(DP),                           intent(in)  :: renci
      real(DP),             dimension(:), intent(in)  :: renc
      real(DP),                           intent(in)  :: dt
      integer(I4B),         dimension(:), intent(in)  :: ind_arr
      type(encounter_list),               intent(out) :: lenci
      ! Internals
      integer(I4B) :: j, nenci
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
      if (nenci > 0) then
         allocate(lenci%lvdotr(nenci), lenci%index2(nenci))
         lenci%nenc = nenci
         lenci%lvdotr(:) = pack(lvdotri(i+1:n), lencounteri(i+1:n)) 
         lenci%index2(:) = pack(ind_arr(i+1:n), lencounteri(i+1:n)) 
      end if

      return
   end subroutine encounter_check_all_triangular_one


   subroutine encounter_check_all_triangular_plpl(npl, x, v, renc, dt, lvdotr, index1, index2, nenc)
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
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out) :: nenc   !! Total number of encounters
      ! Internals
      integer(I4B) :: i, j, k, nenci, j0, j1
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc12
      logical, dimension(npl) :: lencounteri, lvdotri
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      type(encounter_list), dimension(npl) :: lenc
      type(walltimer) :: timer 

      call util_index_array(ind_arr, npl) 

      call timer%reset()
      call timer%start()
      !$omp parallel do default(private) schedule(static)&
      !$omp shared(x, v, renc, lenc, ind_arr) &
      !$omp firstprivate(npl, dt)
      do i = 1,npl
         call encounter_check_all_triangular_one(i, npl, x(1,i), x(2,i), x(3,i), &
                                                         v(1,i), v(2,i), v(3,i), &
                                                         x(1,:), x(2,:), x(3,:), &
                                                         v(1,:), v(2,:), v(3,:), &
                                                         renc(i), renc(:), dt, ind_arr(:), lenc(i))
      end do
      !$omp end parallel do
      call timer%stop()
      write(*,*) "plpl  triang: ",timer%count_stop_step - timer%count_start_step

      call encounter_check_collapse_ragged_list(lenc, npl, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_all_triangular_plpl


   subroutine encounter_check_all_triangular_plplm(nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, lvdotr, index1, index2, nenc)
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
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out) :: nenc   !! Total number of encounters
      ! Internals
      integer(I4B) :: i, j, nenci, j0, j1
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc12
      logical, dimension(nplt) :: lencounteri, lvdotri
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      type(encounter_list), dimension(nplm) :: lenc
      type(walltimer) :: timer 

      call util_index_array(ind_arr, nplt)

      call timer%reset()
      call timer%start()
      !$omp parallel do default(private) schedule(static)&
      !$omp shared(xplm, vplm, xplt, vplt, rencm, renct, lenc, ind_arr) &
      !$omp firstprivate(nplm, nplt, dt)
      do i = 1, nplm
         call encounter_check_all_triangular_one(0, nplt, xplm(1,i), xplm(2,i), xplm(3,i), &
                                                          vplm(1,i), vplm(2,i), vplm(3,i), &
                                                          xplt(1,:), xplt(2,:), xplt(3,:), &
                                                          vplt(1,:), vplt(2,:), vplt(3,:), &
                                                          rencm(i), renct(:), dt, ind_arr(:), lenc(i))
      end do
      !$omp end parallel do
      call timer%stop()
      write(*,*) "plplm triang: ",timer%count_stop_step - timer%count_start_step

      call encounter_check_collapse_ragged_list(lenc, nplm, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_all_triangular_plplm


   subroutine encounter_check_all_triangular_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, lvdotr, index1, index2, nenc)
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
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out) :: nenc   !! Total number of encounters
      ! Internals
      integer(I4B) :: i, j, nenci, j0, j1
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc1, renc2
      logical, dimension(ntp) :: lencounteri, lvdotri
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      type(encounter_list), dimension(npl) :: lenc
      real(DP), dimension(ntp) :: renct

      call util_index_array(ind_arr, ntp)
      renct(:) = 0.0_DP

      !$omp parallel do default(private) schedule(static)&
      !$omp shared(xpl, vpl, xtp, vtp, renc, renct, lenc, ind_arr) &
      !$omp firstprivate(npl, ntp, dt)
      do i = 1, npl
         call encounter_check_all_triangular_one(0, ntp, xpl(1,i), xpl(2,i), xpl(3,i), &
                                                         vpl(1,i), vpl(2,i), vpl(3,i), &
                                                         xtp(1,:), xtp(2,:), xtp(3,:), &
                                                         vtp(1,:), vtp(2,:), vtp(3,:), &
                                                         renc(i), renct(:), dt, ind_arr(:), lenc(i))
      end do
      !$omp end parallel do

      call encounter_check_collapse_ragged_list(lenc, npl, nenc, index1, index2, lvdotr)

      return
   end subroutine encounter_check_all_triangular_pltp


   module pure subroutine encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc, dt, lencounter, lvdotr)
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
      type(encounter_list), dimension(:),              intent(in)            :: ragged_list !! The ragged encounter list
      integer(I4B),                                    intent(in)            :: n1          !! Number of bodies 1
      integer(I4B),                                    intent(out)           :: nenc        !! Total number of encountersj 
      integer(I4B),         dimension(:), allocatable, intent(out)           :: index1      !! Array of indices for body 1
      integer(I4B),         dimension(:), allocatable, intent(out)           :: index2      !! Array of indices for body 1
      logical,              dimension(:), allocatable, intent(out), optional :: lvdotr      !! Array indicating which bodies are approaching
      ! Internals
      integer(I4B) :: i, j0, j1, nenci
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

      !$omp parallel do default(private) &
      !$omp shared(ragged_list, index1, index2, ibeg, lvdotr) &
      !$omp firstprivate(n1)
      do i = 1,n1
         if (ragged_list(i)%nenc == 0) cycle
         nenci = ragged_list(i)%nenc
         j0 = ibeg(i)
         j1 = j0 + nenci - 1
         index1(j0:j1) = i
         index2(j0:j1) = ragged_list(i)%index2(:)
         if (present(lvdotr)) lvdotr(j0:j1) = ragged_list(i)%lvdotr(:)
      end do
      !$omp end parallel do

      return
   end subroutine encounter_check_collapse_ragged_list


   pure subroutine encounter_check_make_ragged_list(lencounteri, ind_arr, nenc,index2)
      implicit none
      ! Arguments
      logical, dimension(:), intent(in) :: lencounteri
      integer(I4B), dimension(:), intent(in) :: ind_arr
      integer(I4B),                            intent(out) :: nenc
      integer(I4B), dimension(:), allocatable, intent(out) :: index2

      nenc = count(lencounteri(:))
      if (nenc > 0) then
         allocate(index2(nenc))
         index2(:) = pack(ind_arr(:), lencounteri(:)) 
      end if

      return
   end subroutine encounter_check_make_ragged_list

  
   module pure subroutine encounter_check_sort_aabb_1D(self, n, extent_arr)
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
      integer(I4B) :: i, j, k, ibox, jbox
      integer(I4B), dimension(2) :: extent
      logical, dimension(:), allocatable :: lfresh

      call util_sort(extent_arr, self%ind)

      do concurrent(k = 1:2*n)
         i = self%ind(k)
         if (i <= n) then
            self%ibeg(i) = k
         else
            self%iend(i - n) = k
         end if
      end do

      return
   end subroutine encounter_check_sort_aabb_1D 

   
   module subroutine encounter_check_sweep_aabb_double_list(self, n1, n2, nenc, index1, index2)
      !! author: David A. Minton
      !!
      !! Sweeps the sorted bounding box extents and returns the encounter candidates
      implicit none
      ! Arguments
      class(encounter_bounding_box),           intent(inout) :: self     !! Multi-dimensional bounding box structure
      integer(I4B),                            intent(in)    :: n1       !! Number of bodies 1
      integer(I4B),                            intent(in)    :: n2       !! Number of bodies 2
      integer(I4B),                            intent(out)   :: nenc     !! Total number of encounter candidates
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1   !! List of indices for body 1 in each encounter candidate pair
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2   !! List of indices for body 2 in each encounter candidate pair
      !Internals
      integer(I4B) :: i, k, ntot, dim
      type(encounter_list), dimension(n1+n2) :: lenc         !! Array of encounter lists (one encounter list per body)
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      integer(I4B), dimension(:), allocatable :: ibeg, iend
      type(walltimer) :: timer 

      ntot = n1 + n2
      call util_index_array(ind_arr, ntot)
      allocate(ibeg(SWEEPDIM * ntot))
      allocate(iend(SWEEPDIM * ntot))
      do dim = 1, SWEEPDIM
         ibeg((dim - 1) * ntot + 1:dim * ntot) = self%aabb(dim)%ibeg(:)
         iend((dim - 1) * ntot + 1:dim * ntot) = self%aabb(dim)%iend(:)
      end do
      ! Sweep the intervals for each of the massive bodies along one dimension
      ! This will build a ragged pair of index lists inside of the lenc data structure
      call timer%reset()
      call timer%start()
      call encounter_check_sweep_aabb_all_double_list(n1, n2, self%aabb(1)%ind(:), reshape(ibeg(:), [SWEEPDIM, ntot]), reshape(iend(:), [SWEEPDIM, ntot]), ind_arr(:), lenc(:))
      call timer%stop()
      write(*,*) "sweep double: ",timer%count_stop_step - timer%count_start_step

      call encounter_check_collapse_ragged_list(lenc, ntot, nenc, index1, index2)

      ! Reorder the pairs and sort the first index in order to remove any duplicates
      do concurrent(k = 1:nenc, index2(k) < index1(k))
         i = index2(k)
         index2(k) = index1(k)
         index1(k) = i
      end do

      return
   end subroutine encounter_check_sweep_aabb_double_list


   module subroutine encounter_check_sweep_aabb_single_list(self, n, nenc, index1, index2)
      !! author: David A. Minton
      !!
      !! Sweeps the sorted bounding box extents and returns the encounter candidates. Mutual encounters
      !! allowed. That is, all bodies are from the same list
      implicit none
      ! Arguments
      class(encounter_bounding_box),           intent(inout) :: self    !! Multi-dimensional bounding box structure
      integer(I4B),                            intent(in)    :: n       !! Number of bodies 1
      integer(I4B),                            intent(out)   :: nenc    !! Total number of encounter candidates
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1  !! List of indices for body 1 in each encounter candidate pair
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2  !! List of indices for body 2 in each encounter candidate pair
      !Internals
      Integer(I4B) :: i, k, dim
      type(encounter_list), dimension(n) :: lenc         !! Array of encounter lists (one encounter list per body)
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      type(walltimer) :: timer 
      integer(I4B), dimension(:), allocatable :: ibeg, iend

      call util_index_array(ind_arr, n)
      allocate(ibeg(SWEEPDIM * n))
      allocate(iend(SWEEPDIM * n))
      do i = 1, n
         do dim = 1, SWEEPDIM
            ibeg((i - 1) * SWEEPDIM + dim) = self%aabb(dim)%ibeg(i)
            iend((i - 1) * SWEEPDIM + dim) = self%aabb(dim)%iend(i)
         end do
      end do

      ! Sweep the intervals for each of the massive bodies along one dimension
      ! This will build a ragged pair of index lists inside of the lenc data structure
      call encounter_check_sweep_aabb_all_single_list(n, self%aabb(1)%ind(:), reshape(ibeg(:), [SWEEPDIM, n]), reshape(iend(:), [SWEEPDIM, n]), ind_arr(:), lenc(:))

      call encounter_check_collapse_ragged_list(lenc, n, nenc, index1, index2)

      ! Reorder the pairs and sort the first index in order to remove any duplicates
      do concurrent(k = 1:nenc, index2(k) < index1(k))
         i = index2(k)
         index2(k) = index1(k)
         index1(k) = i
      end do

      return
   end subroutine encounter_check_sweep_aabb_single_list


   subroutine encounter_check_sweep_aabb_all_double_list(n1, n2, ext_ind, ibeg, iend, ind_arr, lenc)
      !! author: David A. Minton
      !!
      !! Performs the loop part of the sweep operation. Double list version (e.g. pl-tp or plm-plt)
      implicit none
      ! Arguments
      integer(I4B),                         intent(in)    :: n1, n2     !! Number of bodies
      integer(I4B),         dimension(:),   intent(in)    :: ext_ind    !! Sorted index array of extents
      integer(I4B),         dimension(:,:), intent(in)    :: ibeg, iend !! Beginning and ending index lists in the n-dimensions
      integer(I4B),         dimension(:),   intent(in)    :: ind_arr    !! index array for mapping body 2 indexes
      type(encounter_list), dimension(:),   intent(inout) :: lenc       !! Encounter list for the ith body   
      ! Internals
      integer(I4B) :: i, ntot
      logical, dimension(n1+n2) :: lencounteri
      integer(I4B), dimension(SWEEPDIM) :: ibegi, iendi

      ntot = n1 + n2
      !$omp parallel do default(private) schedule(dynamic)&
      !$omp shared(ext_ind, ibeg, iend, ind_arr, lenc) &
      !$omp firstprivate(ntot, n1, n2) 
      do i = 1, ntot
         ibegi(1) = ibeg(1,i) + 1
         iendi(1) = iend(1,i) - 1
         if (iendi(1) >= ibegi(1)) then
            ibegi(2:SWEEPDIM) = ibeg(2:SWEEPDIM,i) 
            iendi(2:SWEEPDIM) = iend(2:SWEEPDIM,i) 
            call encounter_check_sweep_aabb_one_double_list(i, n1, n2, ntot, ext_ind(:), ibegi(:), iendi(:), ibeg(:,:), iend(:,:), lencounteri(:))
            call encounter_check_make_ragged_list(lencounteri(:), ind_arr(:), lenc(i)%nenc, lenc(i)%index2)
         else
            lenc(i)%nenc = 0
         end if
      end do
      !$omp end parallel do

      return
   end subroutine encounter_check_sweep_aabb_all_double_list


   subroutine encounter_check_sweep_aabb_all_single_list(n, ext_ind, ibeg, iend, ind_arr, lenc)
      !! author: David A. Minton
      !!
      !! Performs the loop part of the sweep operation. Single list version (e.g. pl-pl)
      implicit none
      ! Arguments
      integer(I4B),                         intent(in)    :: n          !! Number of bodies
      integer(I4B),         dimension(:),   intent(in)    :: ext_ind    !! Sorted index array of extents
      integer(I4B),         dimension(:,:), intent(in)    :: ibeg, iend !! Beginning and ending index lists in the n-dimension
      integer(I4B),         dimension(:),   intent(in)    :: ind_arr    !! index array for mapping body 2 indexes
      type(encounter_list), dimension(:),   intent(inout) :: lenc       !! Encounter list for the ith body   
      ! Internals
      integer(I4B) :: i
      logical, dimension(n) :: lencounteri
      integer(I4B), dimension(SWEEPDIM) :: ibegi, iendi

      !$omp parallel do default(private) schedule(dynamic)&
      !$omp shared(ext_ind, ibeg, iend, ind_arr, lenc) &
      !$omp firstprivate(n) &
      !$omp lastprivate(ibegi, iendi, lencounteri) 
      do i = 1, n
         ibegi(1) = ibeg(1,i) + 1
         iendi(1) = iend(1,i) - 1
         if (iendi(1) >= ibegi(1)) then
            ibegi(2:SWEEPDIM) = ibeg(2:SWEEPDIM,i) 
            iendi(2:SWEEPDIM) = iend(2:SWEEPDIM,i) 
            call encounter_check_sweep_aabb_one_single_list(n, ext_ind(:), ibegi(:), iendi(:), ibeg(:,:), iend(:,:), lencounteri(:))
            call encounter_check_make_ragged_list(lencounteri(:), ind_arr(:), lenc(i)%nenc, lenc(i)%index2)
         else
            lenc(i)%nenc = 0
         end if
      end do
      !$omp end parallel do

      return
   end subroutine encounter_check_sweep_aabb_all_single_list


   pure subroutine encounter_check_sweep_aabb_one_double_list(i, n1, n2, ntot, ext_ind, ibegi, iendi, ibeg, iend, lencounteri)
      !! author: David A. Minton
      !!
      !! Performs a sweep operation on a single body. Encounters from the same lists not allowed (e.g. pl-tp encounters only)
      implicit none
      ! Arguments
      integer(I4B),                 intent(in)  :: i            !! The current index of the ith body
      integer(I4B),                 intent(in)  :: n1           !! Number of bodies 1
      integer(I4B),                 intent(in)  :: n2           !! Number of bodies 2
      integer(I4B),                 intent(in)  :: ntot         !! n1 + n2
      integer(I4B), dimension(:),   intent(in)  :: ext_ind      !! Sorted index array of extents
      integer(I4B), dimension(:),   intent(in)  :: ibegi, iendi !! The beginning and ending indices of the ith bounding box in the n-dimensions
      integer(I4B), dimension(:,:), intent(in)  :: ibeg, iend   !! Beginning and ending index lists in the n-dimensions
      logical, dimension(:),        intent(out) :: lencounteri    !! Encounter list for the ith body
      ! Internals
      integer(I4B) :: j, jbox, dim, jlo, jhi
      integer(I4B), dimension(:), allocatable :: box

      lencounteri(:) = .false.
      jlo = ibegi(1)
      jhi = iendi(1)
      allocate(box(jlo:jhi))
      box(:) = ext_ind(jlo:jhi)
      where(box(:) > ntot)
         box(:) = box(:) - ntot
      endwhere
      
      do concurrent (jbox = jlo:jhi)
         j = box(jbox)
         if (((i <= n1) .and. (j <= n1)) .or. ((i > n1) .and. (j > n1))) cycle  ! only pairs from the two different lists allowed
         ! Check the other dimensions
         lencounteri(j) = all((iend(2:SWEEPDIM,j) > ibegi(2:SWEEPDIM)) .and. (ibeg(2:SWEEPDIM,j) < iendi(2:SWEEPDIM)))
      end do

      return
   end subroutine encounter_check_sweep_aabb_one_double_list
 

   pure subroutine encounter_check_sweep_aabb_one_single_list(n, ext_ind, ibegi, iendi, ibeg, iend, lencounteri)
      !! author: David A. Minton
      !!
      !! Performs a sweep operation on a single body. Mutual encounters allowed (e.g. pl-pl)
      implicit none
      ! Arguments
      integer(I4B),                 intent(in)  :: n            !! Number of bodies
      integer(I4B), dimension(:),   intent(in)  :: ext_ind      !! Sorted index array of extents
      integer(I4B), dimension(:),   intent(in)  :: ibegi, iendi !! The beginning and ending indices of the ith bounding box in the n-dimensions
      integer(I4B), dimension(:,:), intent(in)  :: ibeg, iend   !! Beginning and ending index lists in the x-dimension
      logical, dimension(:),        intent(out) :: lencounteri  !! Encounter list for the ith body
      ! Internals
      integer(I4B) :: j, jbox, dim, jlo, jhi
      integer(I4B), dimension(:), allocatable :: box

      lencounteri(:) = .false.

      jlo = ibegi(1)
      jhi = iendi(1)

      allocate(box(jlo:jhi))
      box(:) = ext_ind(jlo:jhi)

      where(box(:) > n)
         box(:) = box(:) - n
      endwhere

      do concurrent(jbox = jlo:jhi) ! Sweep forward until the end of the interval
         j = box(jbox)
         ! Check the other dimensions
         lencounteri(j) = all((iend(2:SWEEPDIM,j) > ibegi(2:SWEEPDIM)) .and. (ibeg(2:SWEEPDIM,j) < iendi(2:SWEEPDIM)))
      end do

      return
   end subroutine encounter_check_sweep_aabb_one_single_list

end submodule s_encounter_check
