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
         renc12 = renc1(i) + renc1(j)
         call encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc12, dt, lencounter(k), lvdotr(k)) 
      end do
      !$omp end parallel do simd

      return
   end subroutine encounter_check_all


   module subroutine encounter_check_all_plpl(param, npl, nplm, x, v, renc, dt, lvdotr, index1, index2, nenc)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Choose between the standard triangular or the Sort & Sweep method based on user inputs
      !!
      implicit none
      ! Arguments
      class(swiftest_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
      integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies
      integer(I4B),                            intent(in)    :: nplm   !! Number of fully interacting massive bodies
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
      logical, save :: lsecond = .false.
      integer(I8B) :: nplplm = 0_I8B

      if (param%ladaptive_encounters) then
         nplplm = nplm * npl - nplm * (nplm + 1) / 2 
         if (nplplm > 0) then
            if (lfirst) then  
               write(itimer%loopname, *) "encounter_check_all_plpl"
               write(itimer%looptype, *) "ENCOUNTER"
               lfirst = .false.
               lsecond = .true.
            else
               if (lsecond) then ! This ensures that the encounter check methods are run at least once prior to timing. Sort and sweep improves on the second pass due to the bounding box extents needing to be nearly sorted 
                  call itimer%time_this_loop(param, nplplm)
                  lsecond = .false.
               else if (itimer%check(param, nplplm)) then
                  lsecond = .true.
                  itimer%is_on = .false.
               end if
            end if
         else
            param%lencounter_sas = .false.
         end if
      end if

      if (param%lencounter_sas) then
         call encounter_check_all_sort_and_sweep_plpl(npl, nplm, x, v, renc, dt, lvdotr, index1, index2, nenc) 
      else
         call encounter_check_all_triangular_plpl(npl, nplm, x, v, renc, dt, lvdotr, index1, index2, nenc) 
      end if

      if (.not.lfirst .and. param%ladaptive_encounters .and. nplplm > 0) then 
         if (itimer%is_on) call itimer%adapt(param, nplplm)
      end if

      return
   end subroutine encounter_check_all_plpl


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
      real(DP),     dimension(:,:),            intent(in)    :: xtp    !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)    :: vtp    !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter
      real(DP),                                intent(in)    :: dt     !! Step size
      logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out)   :: nenc   !! Total number of encounters


      call encounter_check_all_triangular_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, lvdotr, index1, index2, nenc) 

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

      ! Reorder the pairs and sort the first index in order to remove any duplicates
      do concurrent(k = 1:nenc, index2(k) < index1(k))
         i = index2(k)
         index2(k) = index1(k)
         index1(k) = i
      end do

      call util_sort(index1, ind)
      lfresh(:) = .true.

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


   subroutine encounter_check_all_sort_and_sweep_plpl(npl, nplm, x, v, renc, dt, lvdotr, index1, index2, nenc)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. 
      !! This is the sort and sweep version
      !! References: Adapted from Christer Ericson's _Real Time Collision Detection_
      !!
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: npl    !! Total number of massive bodies
      integer(I4B),                            intent(in)  :: nplm   !! Number of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: x      !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: v      !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)  :: renc   !! Critical radii of massive bodies that defines an encounter 
      real(DP),                                intent(in)  :: dt     !! Step size
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out) :: nenc   !! Total number of encounters
      ! Internals
      integer(I4B) :: i, j, k, m, dim, n, i0
      logical, dimension(npl) :: lfresh
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      type(encounter_list), dimension(npl) :: lenc
      integer(I4B), dimension(:), allocatable :: tmp, ind
      integer(I4B), save :: npl_last = 0
      type(encounter_bounding_box), save :: boundingbox
      logical, dimension(:), allocatable :: lencounter
      integer(I2B), dimension(npl) :: vshift_min, vshift_max
      integer(I4B) :: ybegi, yendi

      if (npl <= 1) return

      ! If this is the first time through, build the index lists
      n = 2 * npl
      if (npl_last /= npl) then
         if (allocated(ind_arr)) deallocate(ind_arr)
         allocate(ind_arr(npl))
         ind_arr(:) = [(i, i = 1, npl)]

         call boundingbox%setup(npl, npl_last)
      end if

      !$omp taskloop default(private) num_tasks(SWEEPDIM) &
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
         call boundingbox%aabb(dim)%sort(npl, [x(dim,1:npl)-renc(1:npl)+vshift_min(1:npl)*v(dim,1:npl)*dt, &
                                               x(dim,1:npl)+renc(1:npl)+vshift_max(1:npl)*v(dim,1:npl)*dt])
      end do
      !$omp end taskloop

      call boundingbox%sweep(npl, ind_arr, nenc, index1, index2)

      if (nenc > 0) then
         ! Now that we have identified potential pairs, use the narrow-phase process to get the final values
         allocate(lencounter(nenc))
         allocate(lvdotr(nenc))

         call encounter_check_all(nenc, index1, index2, x, v, x, v, renc, renc, dt, lencounter, lvdotr)

         call encounter_check_reduce_broadphase(npl, nenc, index1, index2, lencounter, lvdotr)
      end if

      return
   end subroutine encounter_check_all_sort_and_sweep_plpl


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
      type(encounter_list), dimension(npl) :: lenc
      type(encounter_bounding_box),  save :: boundingbox
      integer(I4B) :: i, j, k, nenci, j0, j1, dim, ibox, jbox, n, ntot
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc12
      logical, dimension(ntp) :: lencounteri
      logical, dimension(npl) :: lfresh
      integer(I4B), dimension(:), allocatable, save :: tpind_arr
      integer(I4B), dimension(:), allocatable :: tmp
      integer(I4B), save :: ntot_last = 0
      integer(I4B), save :: n_last = 0
      integer(I4B), save :: npl_last = 0
      integer(I4B), save :: ntp_last = 0
      logical, dimension(:), allocatable :: lencounter, lvdotr_final
      integer(I2B), dimension(npl) :: vplshift_min, vplshift_max
      integer(I2B), dimension(ntp) :: vtpshift_min, vtpshift_max
      integer(I4B) :: ybegi, yendi
      real(DP), dimension(ntp) :: renctp

      ! If this is the first time through, build the index lists
      if ((ntp == 0) .or. (npl == 0)) return

      ntot = npl + ntp
      n = 2 * ntot
      if ((ntp_last /= ntp) .or. (npl_last /= npl)) then
         if (ntp_last /= ntp) then
            if (allocated(tpind_arr)) deallocate(tpind_arr)
            allocate(tpind_arr(ntp))
            tpind_arr(1:ntp) = [(i, i = 1, ntp)]
         end if

         call boundingbox%setup(ntot, ntot_last)

         npl_last = npl
         ntp_last = ntp
         n_last = n
      end if
      
      !$omp taskloop default(private) num_tasks(SWEEPDIM) &
      !$omp shared(xpl, xtp, vpl, vtp, renc, boundingbox) &
      !$omp firstprivate(dt, npl, n)
      do dim = 1, SWEEPDIM
         where(vpl(dim,1:npl) < 0.0_DP)
            vplshift_min(1:npl) = 1
            vplshift_max(1:npl) = 0
         elsewhere
            vplshift_min(1:npl) = 0
            vplshift_max(1:npl) = 1
         end where

         where(vtp(i,1:ntp) < 0.0_DP)
            vtpshift_min(1:ntp) = 1
            vtpshift_max(1:ntp) = 0
         elsewhere
            vtpshift_min(1:ntp) = 0
            vtpshift_max(1:ntp) = 1
         end where

         call boundingbox%aabb(dim)%sort(npl+ntp, [xpl(dim,1:npl)-renc(1:npl)+vplshift_min(1:npl)*vpl(i,1:npl)*dt, &
                                                   xtp(dim,1:ntp)              +vtpshift_min(1:ntp)*vtp(i,1:ntp)*dt, &
                                                   xpl(dim,1:npl)+renc(1:npl)+vplshift_max(1:npl)*vpl(i,1:npl)*dt, &
                                                   xtp(dim,1:ntp)              +vtpshift_max(1:ntp)*vtp(i,1:ntp)*dt])
      end do
      !$omp end taskloop

      call boundingbox%sweep(npl, ntp, tpind_arr, nenc, index1, index2)

      if (nenc > 0) then
         ! Now that we have identified potential pairs, use the narrow-phase process to get the final values
         allocate(lencounter(nenc))
         allocate(lvdotr(nenc))
         renctp(:) = 0.0_DP

         call encounter_check_all(nenc, index1, index2, xpl, vpl, xtp, vtp, renc, renctp, dt, lencounter, lvdotr)

         call encounter_check_reduce_broadphase(npl, nenc, index1, index2, lencounter, lvdotr)
      end if

      return
   end subroutine encounter_check_all_sort_and_sweep_pltp


   subroutine encounter_check_all_triangular_plpl(npl, nplm, x, v, renc, dt, lvdotr, index1, index2, nenc)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Split off from the main subroutine for performance
      !! This is the upper triangular (double loop) version.
      implicit none
      ! Arguments
      integer(I4B),                            intent(in)  :: npl    !! Total number of massive bodies
      integer(I4B),                            intent(in)  :: nplm   !! Number of fully interacting massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: x      !! Position vectors of massive bodies
      real(DP),     dimension(:,:),            intent(in)  :: v      !! Velocity vectors of massive bodies
      real(DP),     dimension(:),              intent(in)  :: renc  !! Critical radii of massive bodies that defines an encounter 
      real(DP),                                intent(in)  :: dt     !! Step size
      logical,      dimension(:), allocatable, intent(out) :: lvdotr !! Logical flag indicating the sign of v .dot. x
      integer(I4B), dimension(:), allocatable, intent(out) :: index1 !! List of indices for body 1 in each encounter
      integer(I4B), dimension(:), allocatable, intent(out) :: index2 !! List of indices for body 2 in each encounter
      integer(I4B),                            intent(out) :: nenc   !! Total number of encounters
      ! Internals
      integer(I4B) :: i, j, nenci, j0, j1
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc12
      logical, dimension(npl) :: lencounteri, lvdotri
      integer(I4B), dimension(npl) :: ind_arr
      type(encounter_list), dimension(nplm) :: lenc
   
      ind_arr(:) = [(i, i = 1, npl)]
      !$omp parallel do default(private) schedule(static)&
      !$omp shared(npl, nplm, x, v, renc, dt, lenc, ind_arr) &
      !$omp lastprivate(xr, yr, zr, vxr, vyr, vzr, renc12, lencounteri, lvdotri)
      do i = 1, nplm
         do concurrent(j = i+1:npl)
            xr = x(1, j) - x(1, i)
            yr = x(2, j) - x(2, i)
            zr = x(3, j) - x(3, i)
            vxr = v(1, j) - v(1, i)
            vyr = v(2, j) - v(2, i)
            vzr = v(3, j) - v(3, i)
            renc12 = renc(i) + renc(j)
            call encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc12, dt, lencounteri(j), lvdotri(j))
         end do
         nenci = count(lencounteri(i+1:npl))
         if (nenci > 0) then
            allocate(lenc(i)%lvdotr(nenci), lenc(i)%index2(nenci))
            lenc(i)%nenc = nenci
            lenc(i)%lvdotr(:) = pack(lvdotri(i+1:npl), lencounteri(i+1:npl)) 
            lenc(i)%index2(:) = pack(ind_arr(i+1:npl), lencounteri(i+1:npl)) 
         end if
      end do
      !$omp end parallel do

      associate(nenc_arr => lenc(:)%nenc)
         nenc = sum(nenc_arr(1:nplm))
      end associate
      if (nenc > 0) then
         allocate(lvdotr(nenc))
         allocate(index1(nenc))
         allocate(index2(nenc))
         j0 = 1
         do i = 1, nplm
            if (lenc(i)%nenc > 0) then
               j1 = j0 + lenc(i)%nenc - 1
               lvdotr(j0:j1) = lenc(i)%lvdotr(:)
               index1(j0:j1) = i
               index2(j0:j1) = lenc(i)%index2(:)
               j0 = j1 + 1
            end if
         end do
      end if

      return
   end subroutine encounter_check_all_triangular_plpl


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
      integer(I4B), dimension(ntp) :: ind_arr
      type(encounter_list), dimension(npl) :: lenc

      ind_arr(:) = [(i, i = 1, ntp)]
      !$omp parallel do default(private) schedule(static)&
      !$omp shared(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, lenc, ind_arr) &
      !$omp lastprivate(xr, yr, zr, vxr, vyr, vzr, lencounteri, lvdotri)
      do i = 1, npl
         do concurrent(j = 1:ntp)
            xr = xtp(1, j) - xpl(1, i)
            yr = xtp(2, j) - xpl(2, i)
            zr = xtp(3, j) - xpl(3, i)
            vxr = vtp(1, j) - vpl(1, i)
            vyr = vtp(2, j) - vpl(2, i)
            vzr = vtp(3, j) - vpl(3, i)
            call encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc(i), dt, lencounteri(j), lvdotri(j))
         end do
         nenci = count(lencounteri(1:ntp))
         if (nenci > 0) then
            allocate(lenc(i)%lvdotr(nenci), lenc(i)%index2(nenci))
            lenc(i)%nenc = nenci
            lenc(i)%lvdotr(:) = pack(lvdotri(1:ntp), lencounteri(1:ntp)) 
            lenc(i)%index2(:) = pack(ind_arr(1:ntp), lencounteri(1:ntp)) 
         end if
      end do
      !$omp end parallel do

      associate(nenc_arr => lenc(:)%nenc)
         nenc = sum(nenc_arr(1:npl))
      end associate
      if (nenc > 0) then
         allocate(lvdotr(nenc))
         allocate(index1(nenc))
         allocate(index2(nenc))
         j0 = 1
         do i = 1, npl
            if (lenc(i)%nenc > 0) then
               j1 = j0 + lenc(i)%nenc - 1
               lvdotr(j0:j1) = lenc(i)%lvdotr(:)
               index1(j0:j1) = i
               index2(j0:j1) = lenc(i)%index2(:)
               j0 = j1 + 1
            end if
         end do
      end if

      return
   end subroutine encounter_check_all_triangular_pltp


   module pure subroutine encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc, dt, lencounter, lvdotr)
      !$omp declare simd(encounter_check_one)
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
      lencounter = (r2 < r2crit) 
      if (lencounter) return 

      vdotr = vxr * xr + vyr * yr + vzr * zr
      lvdotr = (vdotr < 0.0_DP)
      if (.not.lvdotr) return
     
      v2 = vxr**2 + vyr**2 + vzr**2
      tmin = -vdotr / v2

      if (tmin < dt) then
         r2min = r2 - vdotr**2 / v2
      else
         r2min = r2 + 2 * vdotr * dt + v2 * dt**2
      end if
      lencounter = (r2min <= r2crit)

      return
   end subroutine encounter_check_one

end submodule s_encounter_check
