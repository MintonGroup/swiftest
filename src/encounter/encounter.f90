submodule (swiftest_classes) s_encounter
   use swiftest
contains


   module subroutine encounter_check_all_flat_plpl(nplplm, k_plpl, x, v, renc, dt, lencounter, loc_lvdotr)
      !! author: David A. Minton
      !!
      !! Check for encounters between massive bodies. Split off from the main subroutine for performance.
      !! This is the flat (single loop) version.
      implicit none
      ! Arguments
      integer(I8B),                 intent(in)  :: nplplm     !! Total number of plm-pl encounters to test
      integer(I4B), dimension(:,:), intent(in)  :: k_plpl     !! List of all pl-pl encounters
      real(DP),     dimension(:,:), intent(in)  :: x          !! Position vectors of massive bodies
      real(DP),     dimension(:,:), intent(in)  :: v          !! Velocity vectors of massive bodies
      real(DP),     dimension(:),   intent(in)  :: renc      !! Hill's radii of massive bodies
      real(DP),                     intent(in)  :: dt         !! Step size
      logical,      dimension(:),   intent(out) :: lencounter !! Logical array indicating which pair is in an encounter state
      logical,      dimension(:),   intent(out) :: loc_lvdotr !! Logical array indicating the sign of v .dot. x for each encounter
      ! Internals
      integer(I8B) :: k
      integer(I4B) :: i, j
      real(DP) :: xr, yr, zr, vxr, vyr, vzr, renc12
      
      !$omp parallel do simd default(private) schedule(static)&
      !$omp shared(nplplm, k_plpl, x, v, renc, dt, lencounter, loc_lvdotr) &
      !$omp lastprivate(xr, yr, zr, vxr, vyr, vzr, renc12)
      do k = 1_I8B, nplplm
         i = k_plpl(1, k)
         j = k_plpl(2, k)
         xr = x(1, j) - x(1, i)
         yr = x(2, j) - x(2, i)
         zr = x(3, j) - x(3, i)
         vxr = v(1, j) - v(1, i)
         vyr = v(2, j) - v(2, i)
         vzr = v(3, j) - v(3, i)
         renc12 = renc(i) + renc(j)
         call encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc12, dt, lencounter(k), loc_lvdotr(k))
      end do
      !$omp end parallel do simd

      return
   end subroutine encounter_check_all_flat_plpl


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
      real(DP),     dimension(:),              intent(in)  :: renc  !! Critical radii of massive bodies that defines an encounter
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
