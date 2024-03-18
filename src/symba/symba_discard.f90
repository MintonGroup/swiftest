! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (symba) s_symba_discard
   use swiftest
contains

   subroutine symba_discard_cb_pl(pl, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Check to see if planets should be discarded based on their positions relative to the central body.
      !! If a body gets flagged here when it has also been previously flagged for a collision with another massive body,
      !! its collisional status will be revoked. Discards due to colliding with or escaping the central body take precedence 
      !! over pl-pl collisions
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_sun.f90
      !! Adapted from Hal Levison's Swift routine discard_massive5.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: pl     !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i
      real(DP)     :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2
      character(len=STRMAX) :: idstr, timestr, message
   
      associate(npl => pl%nbody, cb => nbody_system%cb)
         call nbody_system%set_msys()
         rmin2 = param%rmin**2
         rmax2 = param%rmax**2
         rmaxu2 = param%rmaxu**2
         do i = 1, npl
            if (pl%status(i) == ACTIVE) then
               rh2 = dot_product(pl%rh(:,i), pl%rh(:,i))
               if ((param%rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
                  pl%ldiscard(i) = .true.
                  pl%lcollision(i) = .false. 
                  pl%status(i) = DISCARDED_RMAX
                  write(idstr, *) pl%id(i)
                  write(timestr, *) nbody_system%t
                  write(message, *) trim(adjustl(pl%info(i)%name)) // " (" // trim(adjustl(idstr)) // ")" // &
                                    " too far from the central body at t = " // trim(adjustl(timestr))
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "")
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, &
                                                   "***********************************************************" // &
                                                   "***********************************************************")
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, & 
                                                   "***********************************************************" // &
                                                   "***********************************************************")
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "")
                  call pl%info(i)%set_value(status="DISCARDED_RMAX", discard_time=nbody_system%t, discard_rh=pl%rh(:,i), &
                                            discard_vh=pl%vh(:,i))
               else if ((param%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
                  pl%ldiscard(i) = .true.
                  pl%lcollision(i) = .false. 
                  pl%status(i) = DISCARDED_RMIN
                  write(idstr, *) pl%id(i)
                  write(timestr, *) nbody_system%t
                  write(message, *) trim(adjustl(pl%info(i)%name)) // " ("  // trim(adjustl(idstr)) // ")" // &
                                    " too close to the central body at t = " // trim(adjustl(timestr))
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "")
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, &
                                                    "************************************************************" // &
                                                    "************************************************************")
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, &
                                                   "************************************************************" // &
                                                   "************************************************************")
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "")
                  call pl%info(i)%set_value(status="DISCARDED_RMIN", discard_time=nbody_system%t, discard_rh=pl%rh(:,i), &
                                            discard_vh=pl%vh(:,i), discard_body_id=cb%id)
               else if (param%rmaxu >= 0.0_DP) then
                  rb2 = dot_product(pl%rb(:,i), pl%rb(:,i))
                  vb2 = dot_product(pl%vb(:,i), pl%vb(:,i))
                  energy = 0.5_DP * vb2 - nbody_system%Gmtot / sqrt(rb2)
                  if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
                     pl%ldiscard(i) = .true.
                     pl%lcollision(i) = .false. 
                     pl%status(i) = DISCARDED_RMAXU
                     write(idstr, *) pl%id(i)
                     write(timestr, *) nbody_system%t
                     write(message, *) trim(adjustl(pl%info(i)%name)) // " (" // trim(adjustl(idstr)) // ")" // &
                                       " is unbound and too far from barycenter at t = " // trim(adjustl(timestr))
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, "")
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, &
                                                      "************************************************************" // &
                                                      "************************************************************")
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, &
                                                      "************************************************************" // &
                                                      "************************************************************")
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, "")
                     call pl%info(i)%set_value(status="DISCARDED_RMAXU", discard_time=nbody_system%t, discard_rh=pl%rh(:,i), &
                                               discard_vh=pl%vh(:,i))
                  end if
               end if
            end if
         end do
      end associate
   
      return
   end subroutine symba_discard_cb_pl


   subroutine symba_discard_conserve_energy_and_momentum(pl, nbody_system, param, ipl, lescape_body)
      !! author: David A. Minton
      !! 
      !! Conserves nbody_system momentum when a body is lost from the nbody_system or collides with central body
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout) :: pl
      class(symba_nbody_system), intent(inout) :: nbody_system
      class(swiftest_parameters),   intent(inout) :: param
      integer(I4B),              intent(in)    :: ipl
      logical,                   intent(in)         :: lescape_body
      ! Internals
      real(DP), dimension(NDIM) :: Lpl, L_total, Lcb, xcom, vcom, drot0, drot1
      real(DP)                  :: pe, be, ke_orbit, ke_spin, becb0, becb1
      integer(I4B)              :: i, oldstat
   
      select type(cb => nbody_system%cb)
      class is (symba_cb)
   
         ! Add the potential, binding, and kinetic energy of the lost body to the records
         pe = -cb%Gmass * pl%mass(ipl) / norm2(pl%rb(:, ipl) - cb%rb(:))
         be = -3*pl%Gmass(ipl) * pl%mass(ipl) / (5 * pl%radius(ipl))
         ke_orbit = 0.5_DP * pl%mass(ipl) * dot_product(pl%vb(:, ipl), pl%vb(:, ipl)) 
         if (param%lrotation) then
            ke_spin  = 0.5_DP * pl%mass(ipl) * pl%radius(ipl)**2 * pl%Ip(3, ipl) * dot_product(pl%rot(:, ipl), pl%rot(:, ipl))
         else
            ke_spin = 0.0_DP
         end if
   
         if (lescape_body) then
            nbody_system%GMescape = nbody_system%GMescape + pl%Gmass(ipl)
            do i = 1, pl%nbody
               if (i == ipl) cycle
               pe = pe - pl%Gmass(i) * pl%mass(ipl) / norm2(pl%rb(:, ipl) - pl%rb(:, i))
            end do

            nbody_system%E_collisions  = nbody_system%E_collisions + ke_orbit + ke_spin + pe + be
            nbody_system%E_untracked  = nbody_system%E_untracked - (ke_orbit + ke_spin + pe + be)
   
            L_total(:) = 0.0_DP
            do i = 1, pl%nbody
               Lpl(:) = pL%mass(i) * (pl%rb(:,i) .cross. pl%vb(:, i))
               L_total(:) = L_total(:) + Lpl(:)
            end do
            L_total(:) = L_total(:) + cb%mass * (cb%rb(:) .cross. cb%vb(:))
            call pl%b2h(cb)
            oldstat = pl%status(ipl)
            pl%status(ipl) = INACTIVE
            call pl%h2b(cb)
            pl%status(ipl) = oldstat
            do i = 1, pl%nbody
               if (i == ipl) cycle
               Lpl(:) = pl%mass(i) * (pl%rb(:,i) .cross. pl%vb(:, i))
               L_total(:) = L_total(:) - Lpl(:) 
            end do 
            L_total(:) = L_total(:) - cb%mass * (cb%rb(:) .cross. cb%vb(:))
            nbody_system%L_escape(:) = nbody_system%L_escape(:) + L_total(:)
            if (param%lrotation) nbody_system%L_escape(:) = nbody_system%L_escape + pl%mass(ipl) * pl%radius(ipl)**2 &
                                                                    * pl%Ip(3, ipl) * pl%rot(:, ipl)
   
         else
            xcom(:) = (pl%mass(ipl) * pl%rb(:, ipl) + cb%mass * cb%rb(:)) / (cb%mass + pl%mass(ipl))
            vcom(:) = (pl%mass(ipl) * pl%vb(:, ipl) + cb%mass * cb%vb(:)) / (cb%mass + pl%mass(ipl))
            Lpl(:) = (pl%rb(:,ipl) - xcom(:)) .cross. (pL%vb(:,ipl) - vcom(:))
            if (param%lrotation) Lpl(:) = Lpl(:) + pl%radius(ipl)**2 * pl%Ip(3,ipl) * pl%rot(:, ipl)
            Lpl(:) = pl%mass(ipl) * Lpl(:)
     
            Lcb(:) = cb%mass * ((cb%rb(:) - xcom(:)) .cross. (cb%vb(:) - vcom(:)))
   
            ke_orbit = ke_orbit + 0.5_DP * cb%mass * dot_product(cb%vb(:), cb%vb(:)) 
            if (param%lrotation) ke_spin = ke_spin + 0.5_DP * cb%mass * cb%radius**2 * cb%Ip(3) * dot_product(cb%rot(:), cb%rot(:))
            ! Update mass of central body to be consistent with its total mass
            becb0 = -(3 * cb%Gmass * cb%mass) / (5 * cb%radius)
            cb%dGM = cb%dGM + pl%Gmass(ipl)
            cb%dR = cb%dR + 1.0_DP / 3.0_DP * (pl%radius(ipl) / cb%radius)**3 - 2.0_DP / 9.0_DP * (pl%radius(ipl) / cb%radius)**6
            cb%Gmass = cb%GM0 + cb%dGM
            cb%mass = cb%Gmass / param%GU
            cb%radius = cb%R0 + cb%dR
            param%rmin = cb%radius
            becb1 = -(3 * cb%Gmass * cb%mass) / (5 * cb%radius)

            ! Add planet angular momentum to central body accumulator
            cb%dL(:) = Lpl(:) + cb%dL(:) + Lcb(:)
            ! Update rotation of central body to by consistent with its angular momentum 
            if (param%lrotation) then
               drot0(:) = cb%L0(:)/ (cb%Ip(3) * cb%mass * cb%radius**2)  
               drot1(:) = cb%dL(:) / (cb%Ip(3) * cb%mass * cb%radius**2)
               cb%rot(:) = drot0(:) + drot1(:)
               ke_spin  = ke_spin - 0.5_DP * cb%mass * cb%radius**2 * cb%Ip(3) * dot_product(cb%rot(:), cb%rot(:)) 
            end if
            cb%rb(:) = xcom(:)
            cb%vb(:) = vcom(:)
            ke_orbit = ke_orbit - 0.5_DP * cb%mass * dot_product(cb%vb(:), cb%vb(:)) 

            ! Add the change in central body binding energy to the collision energy tracker
            nbody_system%E_collisions  = nbody_system%E_collisions + (becb1 - becb0)
         end if
         call pl%b2h(cb)
   
      end select
      return
   end subroutine symba_discard_conserve_energy_and_momentum


   subroutine symba_discard_nonplpl(pl, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Check to see if planets should be discarded based on their positions or because they are unbound
      !!
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_massive5.f 
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: pl     !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      ! Internals
      logical, dimension(pl%nbody) ::  ldiscard
      integer(I4B) :: i, nstart, nend, nsub
      class(symba_pl), allocatable            :: plsub
    
      ! First check for collisions with the central body
      associate(npl => pl%nbody, cb => nbody_system%cb, pl_discards => nbody_system%pl_discards)
         if (npl == 0) return 
         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or.  (param%rmaxu >= 0.0_DP)) then
            call symba_discard_cb_pl(pl, nbody_system, param)
         end if
         if (param%qmin >= 0.0_DP) call symba_discard_peri_pl(pl, nbody_system, param)
      end associate

      return
   end subroutine symba_discard_nonplpl


   subroutine symba_discard_nonplpl_conservation(pl, nbody_system, param)
      !! author: David A. Minton
      !!
      !! If there are any bodies that are removed due to either colliding with the central body or escaping the systme,
      !! we need to track the conserved quantities with the nbody_system bookkeeping terms.
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout) :: pl     !! SyMBA test particle object
      class(symba_nbody_system), intent(inout) :: nbody_system !! SyMBA nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B)                            :: i, ndiscard, dstat
      logical                                 :: lescape
      logical, dimension(pl%nbody)            :: discard_l_pl
      integer(I4B), dimension(:), allocatable :: discard_index_list

      associate(npl => pl%nbody)
         discard_l_pl(1:npl) = pl%ldiscard(1:npl) .and. .not. pl%lcollision(1:npl) ! These are bodies that are discarded but not 
                                                                                   ! flagged as pl-pl collision
         ndiscard = count(discard_l_pl(:)) 
         allocate(discard_index_list(ndiscard))
         discard_index_list(:) = pack([(i, i = 1, npl)], discard_l_pl(1:npl))
         do i = 1, ndiscard
            dstat = pl%status(discard_index_list(i)) 
            if ((dstat == DISCARDED_RMIN) .or. (dstat == DISCARDED_PERI)) then
               lescape = .false.
            else if ((dstat == DISCARDED_RMAX) .or. (dstat == DISCARDED_RMAXU)) then
               lescape = .true.
            else 
               cycle
            end if
            ! Conserve all the quantities
            call symba_discard_conserve_energy_and_momentum(pl, nbody_system, param, discard_index_list(i), lescape)
         end do
      end associate

      return
   end subroutine symba_discard_nonplpl_conservation


   subroutine symba_discard_peri_pl(pl, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_peri_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_mass_peri.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: pl     !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      logical               :: lfirst_orig
      integer(I4B)          :: i
      character(len=STRMAX) :: timestr, idstr, message


      lfirst_orig = pl%lfirst
      pl%lfirst = nbody_system%lfirst_peri
      if (nbody_system%lfirst_peri) then
         call pl%get_peri(nbody_system, param)          
         nbody_system%lfirst_peri = .false.
      else
         call pl%get_peri(nbody_system, param)          
         do i = 1, pl%nbody
            if (pl%status(i) == ACTIVE) then
               if ((pl%isperi(i) == 0) .and. (pl%nplenc(i)== 0)) then
                  if ((pl%atp(i) >= param%qmin_alo) .and. (pl%atp(i) <= param%qmin_ahi) .and. (pl%peri(i) <= param%qmin)) then
                     pl%ldiscard(i) = .true.
                     pl%lcollision(i) = .false.
                     pl%status(i) = DISCARDED_PERI
                     write(timestr, *) nbody_system%t
                     write(idstr, *) pl%id(i)
                     write(message, *) trim(adjustl(pl%info(i)%name)) // " (" // trim(adjustl(idstr)) // &
                                 ") perihelion distance too small at t = " // trim(adjustl(timestr)) 
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                     call pl%info(i)%set_value(status="DISCARDED_PERI", discard_time=nbody_system%t, &
                                               discard_rh=pl%rh(:,i), discard_vh=pl%vh(:,i), discard_body_id=nbody_system%cb%id)
                  end if
               end if
            end if
         end do
      end if
      pl%lfirst = lfirst_orig
   
      return
   end subroutine symba_discard_peri_pl


   module subroutine symba_discard_pl(self, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Call the various flavors of discards for massive bodies in SyMBA runs, including discards due to colliding with the central
      !! body or escaping the nbody_system
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      real(DP) :: E_orbit_before, E_orbit_after
   
      select type(nbody_system)
      class is (symba_nbody_system)
         select type(param)
         class is (swiftest_parameters)
            associate(pl => self, plpl_encounter => nbody_system%plpl_encounter, plpl_collision => nbody_system%plpl_collision)
               call pl%vb2vh(nbody_system%cb) 
               call pl%rh2rb(nbody_system%cb)
               !call plpl_encounter%write(pl, pl, param) TODO: write the encounter list writer for NetCDF

               call symba_discard_nonplpl(self, nbody_system, param)

               if (.not.any(pl%ldiscard(:))) return

               if (param%lenergy) then
                  call nbody_system%get_energy_and_momentum(param)
                  E_orbit_before = nbody_system%te
               end if

               call symba_discard_nonplpl_conservation(self, nbody_system, param)

               call pl%rearray(nbody_system, param)

               if (param%lenergy) then
                  call nbody_system%get_energy_and_momentum(param)
                  E_orbit_after = nbody_system%te
                  nbody_system%E_collisions = nbody_system%E_collisions + (E_orbit_after - E_orbit_before)
               end if

            end associate
         end select 
      end select

      return
   end subroutine symba_discard_pl

end submodule s_symba_discard