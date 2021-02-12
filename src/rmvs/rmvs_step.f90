submodule(rmvs_classes) s_rmvs_step
contains
   module procedure rmvs_step_system
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step.f90
      use swiftest
      implicit none
      logical, save :: lfirst = .true.
      logical :: lencounter
      real(DP) :: rts
      real(DP), dimension(:,:), allocatable :: xbeg, xend, vbeg
 
      associate(ntp => tp%nbody, npl => pl%nbody, t => config%t, dt => config%dt, &
         xh => pl%xh, vh => pl%vh, xj => pl%xj, vj => pl%vj, ah => pl%ah,  eta => pl%eta, & ! These two lines of associations aid in debugging with gdb
         xht => tp%xh, vht => tp%vh, aht => tp%ah, irij3 => tp%irij3) 
         allocate(xbeg, source=pl%xh)
         allocate(vbeg, source=pl%vh)
         if (allocated(tp%xbeg)) deallocate(tp%xbeg)
         allocate(tp%xbeg, source=xbeg)
         if (allocated(tp%vbeg)) deallocate(tp%vbeg) 
         allocate(tp%vbeg, source=vbeg) 
         ! ****** Check for close encounters ***** !
         rts = RHSCALE
         lencounter = tp%encounter_check(cb, pl, dt, rts)
         if (lencounter) then
            call pl%setup_encounter(cb, tp)
            pl%xout(:,:,0) = tp%xbeg(:,:)
            pl%vout(:,:,0) = tp%vbeg(:,:)
            call pl%step(cb, config, t) 
            pl%xout(:,:,NTENC) = pl%xh(:,:)
            pl%vout(:,:,NTENC) = pl%vh(:,:)
            allocate(xend, source=pl%xh)
            call pl%interp_out(cb, dt)
            call pl%step_out(cb, tp, dt, config)
            where (tp%status(:) == ACTIVE)
               tp%status(:) = INACTIVE
            elsewhere (tp%status(:) == INACTIVE)
               tp%status(:) = ACTIVE
            end where
            if (allocated(tp%xbeg)) deallocate(tp%xbeg)
            allocate(tp%xbeg, source=xbeg)
            if (allocated(tp%xend)) deallocate(tp%xend)
            allocate(tp%xend, source=xend)
            call tp%step(cb, pl, config, t)
            deallocate(tp%xbeg, tp%xend)
            where (tp%status(:) == INACTIVE) 
               tp%status(:) = ACTIVE
            end where
            call pl%destruct_encounter()
         else
            call whm_step_system(cb, pl, tp, config)
         end if
      end associate

   end procedure rmvs_step_system 

   module procedure rmvs_step_in_pl
      !! author: David A. Minton
      !!
      !! Step active test particles ahead in the inner encounter region
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step_in.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step_in.f90
      use swiftest
      implicit none

      logical                                      :: lfirsttp
      integer(I4B)                                 :: i, j, k, nenc, link
      real(DP)                                     :: mu, rhill, dti, time
      real(DP), dimension(:), allocatable          :: irh
      real(DP), dimension(:, :), allocatable       :: xh_original 

   ! executable code
      dti = dt / NTPHENC
      associate(npl => self%nbody)
         if (cb%j2rp2 /= 0.0_DP) then
            allocate(xh_original, source=self%xh)
            allocate(irh(npl))
            do i = 0, NTPHENC
               self%xh(:,:) = self%xin(:,:,i) ! Temporarily replace heliocentric position with inner substep values to calculate the oblateness terms
               do j = 1, npl
                  irh(j) = 1.0_DP / norm2(self%xh(:, j))
               end do
               call self%obl_acc(cb, irh)
               self%aoblin(:,:,i) = self%aobl(:,:) ! Save the oblateness acceleration on the planet for this substep
            end do
            ! Put back the original heliocentric position for the planets
            self%xh(:,:) = xh_original(:,:)
            deallocate(xh_original, irh)
         end if
         do i = 1, npl
            nenc = self%nenc(i) 
            if (nenc > 0) then
            ! There are inner encounters with this planet...switch to planetocentric coordinates to proceed
            ! Determine initial planetocentric positions and velocities for those test particles encountering this planet
               link = self%tpenc1P(i)
               do j = 1, nenc
                  self%tpenc(i)%xh(:, j) = tp%xh(:, link) - self%xin(:, i, 0)
                  self%tpenc(i)%vh(:, j) = tp%vh(:, link) - self%vin(:, i, 0)
                  link = tp%tpencP(link)
               end do
               ! Determine planetocentric positions for planets and sun at all interpolated points in inner encounter
               do k = 0, NTPHENC
                  do j = 1, npl
                     if (j == i) then ! We will substitute the Sun in the array location occupied by the encountering planet
                        self%plenc(i, k)%xh(:, j) = cb%xin(:) - self%xin(:, i, k)
                     else
                        self%plenc(i, k)%xh(:, j) = self%xin(:, i, k) - self%xin(:, i, k)
                     end if
                  end do
               end do
               time = config%t
               mu = self%Gmass(i)
               rhill = self%rhill(i)
               call self%tpenc(i)%peri_pass(cb, self, time, dti, .true., 0, nenc, i, config) 
      ! now step the encountering test particles fully through the inner encounter
               lfirsttp = .true.
               do j = 1, NTPHENC ! Integrate over the encounter region, using the "substitute" planetocentric systems at each level
                  self%tpenc(i)%lfirst = .true.
                  allocate(self%tpenc(i)%xbeg, source=self%plenc(i,j-1)%xh)
                  allocate(self%tpenc(i)%xend, source=self%plenc(i,j)%xh)
                  call self%tpenc(i)%step(self%cbenc(i), self%plenc(i,j), config, time)
                  time = config%t + j * dti
                  call self%tpenc(i)%peri_pass(cb, self, time, dti, .false., j, nenc, i, config) 
                  deallocate(self%tpenc(i)%xbeg)
                  deallocate(self%tpenc(i)%xend)
               end do
               link = self%tpenc1P(i)
               do j = 1, nenc
                  ! Copy the results of the integration back over
                  tp%xh(:, link) = self%xin(:, i, NTPHENC) + self%tpenc(i)%xh(:,j)
                  tp%vh(:, link) = self%vin(:, i, NTPHENC) + self%tpenc(i)%vh(:,j)
                  if (tp%status(link) == ACTIVE) tp%status(link) = INACTIVE
                  link = tp%tpencP(link)

               end do
            end if
         end do
      end associate

      return

   end procedure rmvs_step_in_pl

   module procedure rmvs_step_out
      !! author: David A. Minton
      !!
      !! Step ACTIVE test particles ahead in the outer encounter region
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step_out.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step_out.f90
      use swiftest
      implicit none
      integer(I4B)       :: i, j, k, nenc, itpp
      real(DP)           :: dto, time
      type(rmvs_pl)      :: rmvs_plep

   ! executable code
      dto = dt / NTENC
      associate(npl => self%nbody, ntp => tp%nbody, t => config%t)
         where(tp%plencp(:) == 0)
            tp%status(:) = INACTIVE
         elsewhere  
            tp%lperi(:) = .false.
         end where
         do i = 1, NTENC
            time = t + (i - 1) * dto
            call rmvs_step_out2(cb, self, tp, time, dto, i, config)
            do j = 1, npl
               nenc = self%nenc(j) 
               if (nenc > 0) then
                  itpp = self%tpenc1P(j)
                  do k = 1, nenc
                     if (tp%status(itpp) == INACTIVE) tp%status(itpp) = ACTIVE
                     itpp = tp%tpencp(itpp)
                  end do
               end if
            end do
         end do
      end associate
      return

      end procedure rmvs_step_out   

      module procedure rmvs_step_out2
         !! author: David A. Minton
         !!
         !! Step active test particles ahead in the outer encounter region, setting up and calling the inner region
         !!    integration if necessary
         !! 
         !! Adapted from Hal Levison's Swift routine rmvs3_step_out.f
         !! Adapted from David E. Kaufmann's Swifter routine rmvs_step_out.f90 
   
         use swiftest
         implicit none
   
         logical, save   :: lmalloc = .true.
         integer(I4B)    :: i
         real(DP)        :: rts
         logical         :: lencounter

         if (allocated(tp%xbeg)) deallocate(tp%xbeg)
         if (allocated(tp%vbeg)) deallocate(tp%vbeg)
         if (allocated(tp%xend)) deallocate(tp%xend)
         allocate(tp%xbeg, source=pl%xout(:, :, index-1))
         allocate(tp%vbeg, source=pl%vout(:, :, index-1))
         allocate(tp%xend, source=pl%xout(:, :, index))
         rts = RHPSCALE
         lencounter = tp%encounter_check(cb, pl, dt, rts) 
         if (lencounter) then
            ! Interpolate planets in inner encounter region
            pl%xin(:,:,0) = tp%xbeg(:,:)
            pl%vin(:,:,0) = tp%vbeg(:,:)
            pl%xin(:,:,NTPHENC) = tp%xend(:, :)
            pl%vin(:,:,NTPHENC) = pl%vout(:, :, index)
            call pl%interp_in(cb, dt)
            call pl%step_in(cb, tp, config, dt)
            tp%lfirst = .true.
            call tp%step(cb, pl, config, t)
         else
            call tp%step(cb, pl, config, t)
         end if
         deallocate(tp%xbeg, tp%vbeg, tp%xend)
         return
   
      end procedure rmvs_step_out2

end submodule s_rmvs_step
