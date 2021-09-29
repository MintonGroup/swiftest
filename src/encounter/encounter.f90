submodule (swiftest_classes) s_encounter
   use swiftest
contains

   module subroutine encounter_check_all_sort_and_sweep_plpl(npl, nplm, x, v, renc, dt, lvdotr, index1, index2, nenc)
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
      integer(I4B) :: i, j, k, nenci, j0, j1, dim, ibox, jbox, n, n_last
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc12
      logical, dimension(npl) :: lencounteri, lfresh
      integer(I4B), dimension(:), allocatable, save :: ind_arr
      type lenctype
         logical, dimension(:), allocatable :: lvdotr
         integer(I4B), dimension(:), allocatable :: index2
         integer(I4B) :: nenc
      end type
      type(lenctype), dimension(nplm) :: lenc
      integer(I4B), dimension(:), allocatable :: tmp
      integer(I4B), save :: npl_last = 0
      type boundingBox
         integer(I4B), dimension(:), allocatable :: ind !! Sorted minimum/maximum extent indices
      end type
      type(boundingBox), dimension(NDIM), save :: aabb
      logical, dimension(:), allocatable :: lenc_final, lvdotr_final
      integer(I2B), dimension(npl) :: vshift_min, vshift_max

      if (npl <= 1) return

      ! If this is the first time through, build the index lists
      n = 2 * npl
      n_last = 2 * npl_last
      if (npl_last /= npl) then
         if (allocated(ind_arr)) deallocate(ind_arr)
         allocate(ind_arr(npl))
         ind_arr(:) = [(i, i = 1, npl)]
         if (npl > npl_last) then ! The number of bodies has grown. Resize and append the new bodies
            do i = 1, NDIM
               allocate(tmp(n))
               if (npl_last > 0) tmp(1:n_last) = aabb(i)%ind(1:n_last)
               call move_alloc(tmp, aabb(i)%ind)
               aabb(i)%ind(n_last+1:n) = [(k, k = n_last+1, n)]
            end do
         else ! The number of bodies has gone down. Resize and chop of the old indices
            do i = 1, NDIM
               allocate(tmp(n))
               tmp(1:n) = pack(aabb(i)%ind(1:n_last), aabb(i)%ind(1:n_last) <= n)
               call move_alloc(tmp, aabb(i)%ind)
            end do
         end if
         npl_last = npl
      end if

      do i = 1, NDIM
         where(v(i,1:npl) < 0.0_DP)
            vshift_min(1:npl) = 1
            vshift_max(1:npl) = 0
         elsewhere
            vshift_min(1:npl) = 0
            vshift_max(1:npl) = 1
         end where
         call util_sort([x(i,1:npl)-renc(1:npl)+vshift_min(1:npl)*v(i,1:npl)*dt, &
                         x(i,1:npl)+renc(1:npl)+vshift_max(1:npl)*v(i,1:npl)*dt], aabb(i)%ind)
      end do

      ! Sweep the intervals
      dim = 1
      lfresh(:) = .true. ! This will prevent double-counting of pairs
      do ibox = 1, n
         i = aabb(dim)%ind(ibox)
         if (i > npl) i = i - npl ! If this is an endpoint index, shift it to the correct range
         if (i > nplm) cycle ! Not fully interacting, so move on
         if (.not.lfresh(i)) cycle ! This body has already been evaluated, so move on
         lencounteri(:) = .false.
         do jbox = ibox + 1, n ! Sweep forward until the end of the interval
            j = aabb(dim)%ind(jbox)
            if (j > npl) j = j - npl ! If this is an endpoint index, shift it to the correct range
            if (j == i) exit ! We've reached the end of this interval
            if (lfresh(j)) lencounteri(j) = .true. ! An overlapping box is found that has not previously been tallied
         end do
         lfresh(i) = .false. ! This body has now been processed, so it should no longer show up in future encounter lists
         nenci = count(lencounteri(:))
         if (nenci > 0) then
            allocate(lenc(i)%index2(nenci))
            lenc(i)%nenc = nenci
            lenc(i)%index2(:) = pack(ind_arr(:), lencounteri(:)) 
         end if 
      end do
      
      associate(nenc_arr => lenc(:)%nenc)
         nenc = sum(nenc_arr(1:nplm))
      end associate

      if (nenc > 0) then
         allocate(index1(nenc))
         allocate(index2(nenc))
         allocate(lenc_final(nenc))
         allocate(lvdotr_final(nenc))
         j0 = 1
         do i = 1, nplm
            if (lenc(i)%nenc > 0) then
               j1 = j0 + lenc(i)%nenc - 1
               index1(j0:j1) = i
               index2(j0:j1) = lenc(i)%index2(:)
               j0 = j1 + 1
            end if
         end do
         ! Now that we have identified potential pairs, use the narrow-phase process to get the final values
         lenc_final(:) = .true.

         do k = 1, nenc
            i = index1(k)
            j = index2(k)
            xr = x(1, j) - x(1, i)
            yr = x(2, j) - x(2, i)
            zr = x(3, j) - x(3, i)
            vxr = v(1, j) - v(1, i)
            vyr = v(2, j) - v(2, i)
            vzr = v(3, j) - v(3, i)
            renc12 = renc(i) + renc(j)
            call encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc12, dt, lenc_final(k), lvdotr_final(k)) 
         end do

         nenc = count(lenc_final(:)) ! Count the true number of encounters
         allocate(tmp(nenc))
         tmp(:) = pack(index1(:), lenc_final(:))
         call move_alloc(tmp, index1)
         allocate(tmp(nenc))
         tmp(:) = pack(index2(:), lenc_final(:))
         call move_alloc(tmp, index2)
         allocate(lvdotr(nenc))
         lvdotr(:) = pack(lvdotr_final(:), lenc_final(:))
      end if

      return
   end subroutine encounter_check_all_sort_and_sweep_plpl


   module subroutine encounter_check_all_sort_and_sweep_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, lvdotr, index1, index2, nenc)
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
      integer(I4B) :: i, j, k, nenci, j0, j1, dim, ibox, jbox, n, ntot
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc12
      logical, dimension(ntp) :: lencounteri
      logical, dimension(npl) :: lfresh
      integer(I4B), dimension(:), allocatable, save :: tpind_arr
      type lenctype
         logical, dimension(:), allocatable :: lvdotr
         integer(I4B), dimension(:), allocatable :: index2
         integer(I4B) :: nenc
      end type
      type(lenctype), dimension(npl) :: lenc
      integer(I4B), dimension(:), allocatable :: tmp
      integer(I4B), save :: ntot_last = 0
      integer(I4B), save :: n_last = 0
      integer(I4B), save :: npl_last = 0
      integer(I4B), save :: ntp_last = 0
      type boundingBox
         integer(I4B), dimension(:), allocatable :: ind !! Sorted minimum/maximum extent indices
      end type
      type(boundingBox), dimension(NDIM), save :: aabb
      logical, dimension(:), allocatable :: lenc_final, lvdotr_final
      integer(I2B), dimension(npl) :: vplshift_min, vplshift_max
      integer(I2B), dimension(ntp) :: vtpshift_min, vtpshift_max

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

         if (n > n_last) then ! The number of bodies has grown. Resize and append the new bodies
            do i = 1, NDIM
               allocate(tmp(n))
               if (npl_last > 0) tmp(1:n_last) = aabb(i)%ind(1:n_last)
               call move_alloc(tmp, aabb(i)%ind)
               aabb(i)%ind(n_last+1:n) = [(k, k = n_last+1, n)]
            end do
         else if (n < n_last) then ! The number of bodies has gone down. Resize and chop of the old indices
            do i = 1, NDIM
               allocate(tmp(n))
               tmp(1:n) = pack(aabb(i)%ind(1:n_last), aabb(i)%ind(1:n_last) <= n)
               call move_alloc(tmp, aabb(i)%ind)
            end do
         end if

         npl_last = npl
         ntp_last = ntp
         n_last = n
      end if

      do i = 1, NDIM
         where(vpl(i,1:npl) < 0.0_DP)
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

         call util_sort([xpl(i,1:npl)-renc(1:npl)+vplshift_min(1:npl)*vpl(i,1:npl)*dt, &
                         xtp(i,1:ntp)            +vtpshift_min(1:ntp)*vtp(i,1:ntp)*dt, &
                         xpl(i,1:npl)+renc(1:npl)+vplshift_max(1:npl)*vpl(i,1:npl)*dt, &
                         xtp(i,1:ntp)            +vtpshift_max(1:ntp)*vtp(i,1:ntp)*dt], aabb(i)%ind)
      end do

      ! Sweep the intervals
      dim = 1
      lfresh(:) = .true. ! This will allow us to skip end-points of processed massive bodies
      do ibox = 1, n
         i = aabb(dim)%ind(ibox)
         if (i > ntot) i = i - ntot ! If this is an endpoint index, shift it to the correct range
         if (i > npl) cycle ! This is a test particle, so move on
         if (.not.lfresh(i)) cycle ! This body has already been evaluated, so move on
         lencounteri(:) = .false.
         do jbox = ibox + 1, n ! Sweep forward until the end of the interval
            j = aabb(dim)%ind(jbox)
            if (j > ntot) j = j - ntot ! If this is an endpoint index, shift it to the correct range
            if (j == i) exit ! We've reached the end of this interval
            if (j > npl) then ! this is an unprocessed test particle
               j = j - npl ! Shift into the range of the test particles
               lencounteri(j) = .true. ! An overlapping tp/pl is found that has not previously been tallied
            end if
         end do
         lfresh(i) = .false. ! This body has now been processed, so it should no longer show up in future encounter lists
         nenci = count(lencounteri(:))
         if (nenci > 0) then
            allocate(lenc(i)%index2(nenci))
            lenc(i)%nenc = nenci
            lenc(i)%index2(:) = pack(tpind_arr(:), lencounteri(:)) 
         end if 
      end do
      
      associate(nenc_arr => lenc(:)%nenc)
         nenc = sum(nenc_arr(1:npl))
      end associate

      if (nenc > 0) then
         allocate(index1(nenc))
         allocate(index2(nenc))
         allocate(lenc_final(nenc))
         allocate(lvdotr_final(nenc))
         j0 = 1
         do i = 1, npl
            if (lenc(i)%nenc > 0) then
               j1 = j0 + lenc(i)%nenc - 1
               index1(j0:j1) = i
               index2(j0:j1) = lenc(i)%index2(:)
               j0 = j1 + 1
            end if
         end do
         ! Now that we have identified potential pairs, use the narrow-phase process to get the final values
         lenc_final(:) = .true.

         do k = 1, nenc
            i = index1(k)
            j = index2(k)
            xr = xtp(1, j) - xpl(1, i)
            yr = xtp(2, j) - xpl(2, i)
            zr = xtp(3, j) - xpl(3, i)
            vxr = vtp(1, j) - vpl(1, i)
            vyr = vtp(2, j) - vpl(2, i)
            vzr = vtp(3, j) - vpl(3, i)
            call encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc(i), dt, lenc_final(k), lvdotr_final(k)) 
         end do

         nenc = count(lenc_final(:)) ! Count the true number of encounters
         allocate(tmp(nenc))
         tmp(:) = pack(index1(:), lenc_final(:))
         call move_alloc(tmp, index1)
         allocate(tmp(nenc))
         tmp(:) = pack(index2(:), lenc_final(:))
         call move_alloc(tmp, index2)
         allocate(lvdotr(nenc))
         lvdotr(:) = pack(lvdotr_final(:), lenc_final(:))
      end if

      return
   end subroutine encounter_check_all_sort_and_sweep_pltp


   module subroutine encounter_check_all_triangular_plpl(npl, nplm, x, v, renc, dt, lvdotr, index1, index2, nenc)
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
      type lenctype
         logical, dimension(:), allocatable :: lvdotr
         integer(I4B), dimension(:), allocatable :: index2
         integer(I4B) :: nenc
      end type
      type(lenctype), dimension(nplm) :: lenc
   
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


   module subroutine encounter_check_all_triangular_pltp(npl, ntp, xpl, vpl, xtp, vtp, renc, dt, lvdotr, index1, index2, nenc)
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
      type lenctype
         logical, dimension(:), allocatable :: lvdotr
         integer(I4B), dimension(:), allocatable :: index2
         integer(I4B) :: nenc
      end type
      type(lenctype), dimension(npl) :: lenc


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
      !! Adapted from David E. Kaufmann's Swifter routine: encounter_check_one.f90
      !! Adapted from Hal Levison's Swift routine encounter_check_one.f
      implicit none
      ! Arguments
      real(DP), intent(in)  :: xr, yr, zr    !! Relative distance vector components
      real(DP), intent(in)  :: vxr, vyr, vzr !! Relative velocity vector components
      real(DP), intent(in)  :: renc        !! Square of the critical encounter distance
      real(DP), intent(in)  :: dt            !! Step size
      logical,  intent(out) :: lencounter    !! Flag indicating that an encounter has occurred
      logical,  intent(out) :: lvdotr        !! Logical flag indicating the direction of the v .dot. r vector
      ! Internals
      real(DP) :: r2crit, r2min, r2, v2, vdotr

      r2 = xr**2 + yr**2 + zr**2
      r2crit = renc**2
      lencounter = (r2 < r2crit) 
      if (lencounter) return 

      vdotr = vxr * xr + vyr * yr + vzr * zr
      lvdotr = (vdotr < 0.0_DP)
      if (.not.lvdotr) return
     
      v2 = vxr**2 + vyr**2 + vzr**2

      if (-vdotr < v2 * dt) then
         r2min = r2 - vdotr**2 / v2
      else
         r2min = r2 + 2 * vdotr * dt + v2 * dt**2
      end if
      lencounter = (r2min <= r2crit)

      return
   end subroutine encounter_check_one

end submodule s_encounter
