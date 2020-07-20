submodule(whm_classes) whm_drift
contains

   module procedure whm_drift_pl
      !! author: David A. Minton
      !!
      !! Loop through planets and call Danby drift routine
      !!
      !! Adapted from Hal Levison's Swift routine drift.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift.f90
      use swiftest
      implicit none
      integer(I4B)          :: i
      real(DP)     :: energy, vmag2, rmag  !! Variables used in GR calculation
      integer(I4B), dimension(:), allocatable  :: iflag
      real(DP), dimension(:), allocatable :: dtp
      real(DP) :: px,py,pz,vx,vy,vz

      associate(npl => self%nbody, xj => self%xj, vj => self%vj, status => self%status, mu => self%muj, &
               lgr => config%lgr, c2 => config%inv_c2)
         if (npl == 0) return

         allocate(iflag(npl))
         iflag = 0
         allocate(dtp(npl))
         iflag(:) = 0
         
         do i = 1, npl
            if (lgr) then
               rmag = norm2(xj(:, i))
               vmag2 = dot_product(vj(:, i),  vj(:, i))
               energy = 0.5_DP * vmag2 - mu(i) / rmag
               dtp(i) = dt * (1.0_DP + 3 * c2 * energy)
            else
               dtp(i) = dt
            end if
         end do
         !dir$ parallel always
         do i = 1, npl
            call drift_one(mu(i), xj(1, i), xj(2, i), xj(3, i), vj(1, i), vj(2, i), vj(3, i), dtp(i), iflag(i),&
                           px, py, pz, vx, vy, vz)
            xj(1, i) = px
            xj(2, i) = py
            xj(3, i) = pz
            vj(1, i) = vx
            vj(2, i) = vy
            vj(3, i) = vz
         end do 
         if (any(iflag(1:npl) /= 0)) then
            do i = 1, npl
               if (iflag(i) /= 0) then
                  write(*, *) " Planet ", self%name(i), " is lost!!!!!!!!!!"
                  write(*, *) xj
                  write(*, *) vj
                  write(*, *) " stopping "
               end if
            end do
            call util_exit(FAILURE)
         end if
      end associate

      return

   end procedure whm_drift_pl

   module procedure whm_drift_tp
      !! author: David A. Minton
      !!
      !! Loop through test particles and call Danby drift routine
      !!
      !! Adapted from Hal Levison's Swift routine drift_tp.f 
      !! Includes 
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift_tp.f90
      use swiftest
      implicit none
      integer(I4B)                            :: i   
      real(DP)     ::  energy, vmag2, rmag  !! Variables used in GR calculation
      integer(I4B), dimension(:), allocatable  :: iflag
      real(DP), dimension(:), allocatable :: dtp
      real(DP), dimension(:,:), allocatable :: xnew, vnew
      real(DP) :: px,py,pz,vx,vy,vz

      associate(ntp => self%nbody, xh => self%xh, vh => self%vh, status => self%status, mu => self%mu, &
               lgr => config%lgr, c2 => config%inv_c2)
         if (ntp == 0) return
         allocate(iflag(ntp))
         iflag = 0
         allocate(dtp(ntp))
         allocate(xnew, mold = self%xh)
         allocate(vnew, mold = self%vh)
         iflag(:) = 0
         do i = 1, ntp
            if (lgr) then
               rmag = norm2(xh(:, i))
               vmag2 = dot_product(vh(:, i), vh(:, i))
               energy = 0.5_DP * vmag2 - mu(i) / rmag
               dtp(i) = dt * (1.0_DP + 3 * c2 * energy)
            else
               dtp(i) = dt
            end if
         end do
         call annotate_site_begin( "drift_tp_loop" ) 
         !$omp simd safelen(1200)
         do i = 1, ntp
            call annotate_iteration_task( "drift_tp_loop" )
            call drift_one(mu(i), xh(1, i), xh(2, i), xh(3, i), vh(1, i), vh(2, i), vh(3, i), dtp(i), iflag(i),&
                           px, py, pz, vx, vy, vz)
            xh(1, i) = px
            xh(2, i) = py
            xh(3, i) = pz
            vh(1, i) = vx
            vh(2, i) = vy
            vh(3, i) = vz
         end do
         call annotate_site_end 
         if (any(iflag(1:ntp) /= 0)) then
            do i = 1, ntp
               if (iflag(i) /= 0) then
                  write(*, *) "Particle ", self%name(i), " lost due to error in Danby drift"
                  status(i) = DISCARDED_DRIFTERR
               end if
            end do
         end if
      end associate

      return

   end procedure whm_drift_tp
end submodule whm_drift
