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
         call tp%set_beg_end(xbeg = xbeg, vbeg = vbeg)
         ! ****** Check for close encounters ***** !
         rts = RHSCALE
         lencounter = tp%encounter_check(cb, pl, dt, rts)
         if (lencounter) then
            pl%xout(:,:,0) = tp%xbeg(:,:)
            pl%vout(:,:,0) = tp%vbeg(:,:)
            call pl%step(cb, config, t, dt) 
            pl%xout(:,:,NTENC) = pl%xh(:,:)
            pl%vout(:,:,NTENC) = pl%vh(:,:)
            call tp%set_beg_end(xend = pl%xh)
            call pl%interp_out(cb, dt)
            call pl%step_out(cb, tp, dt, config)
            call tp%reverse_status()
            call tp%set_beg_end(xbeg = xbeg, xend = xend)
            tp%lfirst = .true.
            call tp%step(cb, pl, config, t, dt)
            call tp%reverse_status()
         else
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
      dto = dt / NTENC
      associate(pl => self, npl => self%nbody, ntp => tp%nbody, t => config%t, xht => tp%xh, vht => tp%vh)
         where(tp%plencp(:) == 0)
            tp%status(:) = INACTIVE
         elsewhere  
            tp%lperi(:) = .false.
         end where
         do i = 1, NTENC
            time = t + (i - 1) * dto
            call rmvs_step_out2(cb, pl, tp, time, dto, i, config)
            do j = 1, npl
               nenc = pl%nenc(j) 
               if (nenc > 0) then
                  itpp = pl%tpenc1P(j)
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

         associate(xht => tp%xh, vht => tp%vh, xbeg => tp%xbeg, xend => tp%xend, aht => tp%ah)

            call tp%set_beg_end(xbeg = pl%xout(:, :, index-1), &
                                vbeg = pl%vout(:, :, index-1), &
                                xend = pl%xout(:, :, index))
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
               call tp%step(cb, pl, config, t, dt)
            else
               call tp%step(cb, pl, config, t, dt)
            end if
         end associate
         return
   
      end procedure rmvs_step_out2

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

   
         dti = dt / NTPHENC
         associate(pl => self, npl => self%nbody, xht => tp%xh, vht => tp%vh)
            if (cb%j2rp2 /= 0.0_DP) call pl%obl_acc_in(cb)
            call pl%make_planetocentric(cb, tp)
            do i = 1, npl
               nenc = pl%nenc(i) 
               if (nenc > 0) then
               ! There are inner encounters with this planet...switch to planetocentric coordinates to proceed
                  time = config%t
                  mu = pl%Gmass(i)
                  rhill = pl%rhill(i)
                  call pl%tpenc(i)%peri_pass(cb, pl, time, dti, .true., 0, nenc, i, config) 
         ! now step the encountering test particles fully through the inner encounter
                  lfirsttp = .true.
                 
                  associate(xpc => self%tpenc(i)%xh, vpc => self%tpenc(i)%vh, apc => self%tpenc(i)%ah)
                  do j = 1, NTPHENC ! Integrate over the encounter region, using the "substitute" planetocentric systems at each level
                     pl%tpenc(i)%lfirst = .true.
                     call pl%tpenc(i)%set_beg_end(xbeg = pl%plenc(i,j-1)%xh, xend = pl%plenc(i,j)%xh )
                     call pl%tpenc(i)%step(pl%cbenc(i), pl%plenc(i,j), config, time, dti)
                     time = config%t + j * dti
                     call pl%tpenc(i)%peri_pass(cb, pl, time, dti, .false., j, nenc, i, config) 
                  end do
                  end associate

               end if
            end do
            call pl%end_planetocentric(tp)

         end associate
   
         return
   
      end procedure rmvs_step_in_pl

      module procedure rmvs_step_in_obl_acc
         !! author: David A. Minton
         !!
         !! Compute the oblateness acceleration in the inner encounter region with planets 
         !! 
   
         use swiftest
         implicit none
         integer(I4B) :: i
         real(DP), dimension(:, :), allocatable       :: xh_original 

         associate(pl => self)
            allocate(xh_original, source=pl%xh)
            do i = 0, NTPHENC
               pl%xh(:,:) = pl%xin(:,:,i) ! Temporarily replace heliocentric position with inner substep values to calculate the oblateness terms
               pl%aoblin(:,:,i) = pl%aobl(:,:) ! Save the oblateness acceleration on the planet for this substep
            end do
            call self%obl_acc(cb)
         ! Put back the original heliocentric position for the planets
            pl%xh(:,:) = xh_original(:,:)
            deallocate(xh_original)
         end associate
   
         return
      end procedure rmvs_step_in_obl_acc

      module procedure rmvs_step_make_planetocentric
         !! author: David A. Minton
         !!
         !! When encounters are detected, this method will call the interpolation methods for the planets and 
         !! creates a Swiftest test particle structure for each planet's encountering test particles to simplify the 
         !! planetocentric calculations. This subroutine is not based on an existing one from Swift and Swifter
         !!
         implicit none
         integer(I4B) :: i, j, k, link, nenc
   
         associate(pl => self, npl => self%nbody, nenc => self%nenc)
            allocate(pl%tpenc(npl))
            allocate(pl%plenc(npl, 0:NTPHENC))
            allocate(pl%cbenc(npl))
            do i = 1, npl
               if (nenc(i) > 0) then
                  ! There are inner encounters with this planet...first make the planet a central body
                  pl%cbenc(i)%Gmass = pl%Gmass(i)
   
                  ! Next create an encountering test particle structure
                  call pl%tpenc(i)%setup(nenc(i))  
                  link = pl%tpenc1P(i)
                  do j = 1, nenc(i)
                     pl%tpenc(i)%name(j) = tp%name(link)
                     pl%tpenc(i)%status(j) = tp%status(link)
                     pl%tpenc(i)%xh(:, j) = tp%xh(:, link) - pl%xin(:, i, 0)
                     pl%tpenc(i)%vh(:, j) = tp%vh(:, link) - pl%vin(:, i, 0)
                     link = tp%tpencP(link)
                  end do
   
                  call pl%tpenc(i)%set_mu(pl%cbenc(i)) ! Make sure that the test particles get the proper value of mu 
                  
                  ! Now create a planetocentric "planet" structure containing the *other* planets (plus the Sun) in it
                  do j = 0, NTPHENC
                     call pl%plenc(i, j)%setup(npl)
                     pl%plenc(i, j)%status(:) = pl%status(:)
                  end do
                  do k = 0, NTPHENC
                     do j = 1, npl
                        if (j == i) then ! We will substitute the Sun in the array location occupied by the encountering planet
                           pl%plenc(i, k)%name(j) = 0
                           pl%plenc(i, k)%Gmass(j) = cb%Gmass
                           pl%plenc(i, k)%mass(j) = cb%mass
                           pl%plenc(i, k)%radius(j) = cb%radius
                           pl%plenc(i, k)%xh(:, j) = cb%xin(:) - self%xin(:, i, k)
                        else
                           pl%plenc(i, k)%name(j) = self%name(i)
                           pl%plenc(i, k)%Gmass(j) = pl%Gmass(i)
                           pl%plenc(i, k)%mass(j) = pl%mass(i)
                           pl%plenc(i, k)%xh(:, j) = pl%xin(:, i, k) - pl%xin(:, i, k)
                        end if
                     end do
                  end do
   
               end if
   
            end do
         end associate
   
   
      end procedure rmvs_step_make_planetocentric
   
      module procedure rmvs_step_end_planetocentric
         !! author: David A. Minton
         !!
         !! Deallocates all of the encountering particle data structures for next time
         !!
         use swiftest
         implicit none

         integer(I4B) :: link, i, j

         associate(pl => self, nenc => self%nenc, npl => self%nbody)

            do i = 1, npl
               link = pl%tpenc1P(i)
               do j = 1, nenc(i)
                  ! Copy the results of the integration back over
                  tp%xh(:, link) = pl%xin(:, i, NTPHENC) + pl%tpenc(i)%xh(:,j)
                  tp%vh(:, link) = pl%vin(:, i, NTPHENC) + pl%tpenc(i)%vh(:,j)
                  if (tp%status(link) == ACTIVE) tp%status(link) = INACTIVE
                  link = tp%tpencP(link)
               end do
            end do
      
      
            deallocate(pl%tpenc)
            deallocate(pl%cbenc)
            deallocate(pl%plenc)

         end associate
   
         return
   
      end procedure rmvs_step_end_planetocentric
   
   
 
end submodule s_rmvs_step
