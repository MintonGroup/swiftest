submodule(rmvs_classes) rmvs_step
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
      !real(DP), dimension(:,:), allocatable :: xbeg, vbeg, xend
 
      associate(ntp => tp%nbody, npl => pl%nbody, t => config%t, dt => config%dt, &
         xh => pl%xh, vh => pl%vh, xj => pl%xj, vj => pl%vj, ah => pl%ah,  eta => pl%eta, & ! These two lines of associations aid in debugging with gdb
         xht => tp%xh, vht => tp%vh, aht => tp%ah, irij3 => tp%irij3) 
         allocate(tp%xbeg, source=pl%xh)
         allocate(tp%vbeg, source=pl%vh)
         ! ****** Check for close encounters ***** !
         rts = RHSCALE
         lencounter = tp%encounter_check(cb, pl, dt, rts)
         if (lencounter) then
            pl%xout(:,:,0) = tp%xbeg(:,:)
            pl%vout(:,:,0) = tp%vbeg(:,:)
            call pl%step(cb, config, t) 
            pl%xout(:,:,NTENC) = pl%xh(:,:)
            pl%vout(:,:,NTENC) = pl%vh(:,:)
            call pl%interp_out(cb, dt)
            call pl%step_out(cb, tp, dt, config)
            where (tp%status(:) == ACTIVE)
               tp%status(:) = INACTIVE
            elsewhere (tp%status(:) == INACTIVE)
               tp%status(:) = ACTIVE
            end where
            call tp%step(cb, pl, config, t)
            where (tp%status(:) == INACTIVE) 
               tp%status(:) = ACTIVE
            end where
         else
            deallocate(tp%xbeg)
            deallocate(tp%vbeg)
            call whm_step_system(cb, pl, tp, config)
         end if
      end associate

   end procedure rmvs_step_system 

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
      dto = dt/NTENC
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
  
         allocate(tp%xbeg, source=pl%xout(:, :, index-1))
         allocate(tp%vbeg, source=pl%vout(:, :, index-1))
         allocate(tp%xend, source=pl%xout(:, :, index-1))
         rts = RHPSCALE
         lencounter = tp%encounter_check(cb, pl, dt, rts) 
         if (lencounter) then
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

end submodule rmvs_step
