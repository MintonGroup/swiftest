!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle) s_fraggle_generate
   use swiftest
   use symba

contains

   module subroutine fraggle_generate(self, nbody_system, param, t)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic disruption collision
      !! 
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      real(DP),                 intent(in)    :: t            !! Time of collision
      ! Internals
      integer(I4B)          :: i, ibiggest, nfrag
      character(len=STRMAX) :: message 

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(pl => nbody_system%pl)
      class is (swiftest_pl)
         associate(impactors => self%impactors, status => self%status)
            select case (impactors%regime) 
            case (COLLRESOLVE_REGIME_HIT_AND_RUN)
               call self%hitandrun(nbody_system, param, t)
               return
            case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
               call self%merge(nbody_system, param, t) ! Use the default collision model, which is merge
               return
            case(COLLRESOLVE_REGIME_DISRUPTION)
               message = "Disruption between"
            case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
               message = "Supercatastrophic disruption between"
            case default 
               write(*,*) "Error in swiftest_collision, unrecognized collision regime"
               call base_util_exit(FAILURE)
            end select
            call self%set_mass_dist(param) 
            call self%disrupt(nbody_system, param, t)

            associate (fragments => self%fragments)
               ! Populate the list of new bodies
               nfrag = fragments%nbody
               write(message, *) nfrag
               select case(impactors%regime)
               case(COLLRESOLVE_REGIME_DISRUPTION)
                  status = DISRUPTED
                  ibiggest = impactors%id(maxloc(pl%Gmass(impactors%id(:)), dim=1))
                  fragments%id(1) = pl%id(ibiggest)
                  fragments%id(2:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag - 1)]
                  param%maxid = fragments%id(nfrag)
               case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
                  status = SUPERCATASTROPHIC
                  fragments%id(1:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag)]
                  param%maxid = fragments%id(nfrag)
               end select

               call collision_resolve_mergeaddsub(nbody_system, param, t, status)
 
            end associate
         end associate
      end select
      end select
      return
   end subroutine fraggle_generate


   module subroutine fraggle_generate_disrupt(self, nbody_system, param, t, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Generates a nbody_system of fragments in barycentric coordinates that conserves energy and momentum.
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self         !! Fraggle system object the outputs will be the fragmentation 
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      real(DP),                 intent(in)    :: t            !! Time of collision 
      logical, optional,        intent(out)   :: lfailure     !! Answers the question: Should this have been a merger instead?
       ! Internals
      logical                              :: lk_plpl, lfailure_local
      logical, dimension(size(IEEE_ALL))   :: fpe_halting_modes, fpe_quiet_modes
      real(DP)                             :: dE, dL
      character(len=STRMAX)                :: message
      real(DP), parameter                  :: fail_scale_initial = 1.001_DP

      ! The minimization and linear solvers can sometimes lead to floating point exceptions. Rather than halting the code entirely if this occurs, we
      ! can simply fail the attempt and try again. So we need to turn off any floating point exception halting modes temporarily 
      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      fpe_quiet_modes(:) = .false.
      call ieee_set_halting_mode(IEEE_ALL,fpe_quiet_modes)

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
      associate(impactors => self%impactors, fragments => self%fragments, nfrag => self%fragments%nbody, pl => nbody_system%pl)

         write(message,*) nfrag
         call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Fraggle generating " // trim(adjustl(message)) // " fragments.")

         if (param%lflatten_interactions) then
            lk_plpl = allocated(pl%k_plpl)
            if (lk_plpl) deallocate(pl%k_plpl)
         else 
            lk_plpl = .false.
         end if
         call ieee_set_flag(ieee_all, .false.) ! Set all fpe flags to quiet

         call self%set_natural_scale()
         lfailure_local = .false.
         call self%get_energy_and_momentum(nbody_system, param, phase="before")
         self%fail_scale = fail_scale_initial
         call fraggle_generate_pos_vec(self)
         call fraggle_generate_rot_vec(self, nbody_system, param)
         call fraggle_generate_vel_vec(self, nbody_system, param, lfailure_local)
         call self%get_energy_and_momentum(nbody_system, param, phase="after")
         call self%set_original_scale()

         dE = self%te(2) - self%te(1)
         dL = .mag.(self%L_total(:,2) - self%L_total(:,1))

         write(message,*) dE
         call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Estimated energy change: " // trim(adjustl(message)))
         write(message,*) dL
         call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Estimated angular momentum change: " // trim(adjustl(message)))


         ! Restore the big array
         if (lk_plpl) call pl%flatten(param)
         if (present(lfailure)) lfailure = lfailure_local
      end associate
      end select
      end select
      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Restore the original halting modes

      return 
   end subroutine fraggle_generate_disrupt


   module subroutine fraggle_generate_hitandrun(self, nbody_system, param, t) 
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic hit-and-run collision
      !! 
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self         !! Collision system object
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      real(DP),                 intent(in)    :: t            !! Time of collision
      ! Result
      integer(I4B)                            :: status       !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                            :: i, ibiggest, jtarg, jproj, nfrag
      logical                                 :: lpure 
      character(len=STRMAX) :: message

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(pl => nbody_system%pl)
      class is (swiftest_pl)
         associate(impactors => self%impactors)
            call collision_io_collider_message(nbody_system%pl, impactors%id, message)
            if (impactors%mass(1) > impactors%mass(2)) then
               jtarg = 1
               jproj = 2
            else
               jtarg = 2
               jproj = 1
            end if

            ! The Fraggle disruption model (and its extended types allow for non-pure hit and run. 
            if (impactors%mass_dist(2) > 0.9_DP * impactors%mass(jproj)) then ! Pure hit and run, so we'll just keep the two bodies untouched
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Pure hit and run. No new fragments generated.")
               nfrag = 0
               call self%collision_basic%hitandrun(nbody_system, param, t)
               lpure = .true.
               return
            end if
            lpure = .false.
            call self%set_mass_dist(param)
            message = "Hit and run between"
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, trim(adjustl(message)))

            ! Generate the position and velocity distributions of the fragments
            call self%disrupt(nbody_system, param, t, lpure)
            nfrag = self%fragments%nbody

            write(message, *) nfrag
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Generating " // trim(adjustl(message)) // " fragments")

            ibiggest = impactors%id(maxloc(pl%Gmass(impactors%id(:)), dim=1))
            self%fragments%id(1) = pl%id(ibiggest)
            self%fragments%id(2:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag - 1)]
            param%maxid = self%fragments%id(nfrag)
            status = HIT_AND_RUN_DISRUPT
            call collision_resolve_mergeaddsub(nbody_system, param, t, status)
         end associate
      end select
      end select

      return
   end subroutine fraggle_generate_hitandrun


   module subroutine fraggle_generate_pos_vec(collider)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Initializes the position vectors of the fragments around the center of mass based on the collision style.
      !! For hit and run with disruption, the fragments are generated in a random cloud around the smallest of the two colliders (body 2)
      !! For disruptive collisions, the fragments are generated in a random cloud around the impact point. Bodies are checked for overlap and
      !! regenerated if they overlap.
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: collider !! Fraggle collision system object
      ! Internals
      real(DP)  :: dis, direction, rdistance
      real(DP), dimension(NDIM,2) :: fragment_cloud_center
      real(DP), dimension(2) :: fragment_cloud_radius
      logical, dimension(collider%fragments%nbody) :: loverlap
      real(DP), dimension(collider%fragments%nbody) :: mass_rscale
      integer(I4B) :: i, j, loop, istart
      logical :: lsupercat, lhitandrun
      integer(I4B), parameter :: MAXLOOP = 20000
      real(DP), parameter :: rdistance_scale_factor = 1.0_DP ! Scale factor to apply to distance scaling of cloud centers in the event of overlap
                                                             ! The distance is chosen to be close to the original locations of the impactors
                                                             ! but far enough apart to prevent a collisional cascade between fragments 
      real(DP), parameter :: cloud_size_scale_factor = 3.0_DP ! Scale factor to apply to the size of the cloud relative to the distance from the impact point. 
                                                               ! A larger value puts more space between fragments initially

      associate(fragments => collider%fragments, impactors => collider%impactors, nfrag => collider%fragments%nbody)
         lsupercat = (impactors%regime == COLLRESOLVE_REGIME_SUPERCATASTROPHIC) 
         lhitandrun = (impactors%regime == COLLRESOLVE_REGIME_HIT_AND_RUN) 

         ! We will treat the first two fragments of the list as special cases. 
         ! Place the first two bodies at the centers of the two fragment clouds, but be sure they are sufficiently far apart to avoid overlap
         loverlap(:) = .true.
         if (lhitandrun) then
            rdistance = impactors%radius(2)
         else if (lsupercat) then
            rdistance = 0.5_DP * sum(impactors%radius(:))
         else
            rdistance = 10 * impactors%radius(2)
         end if
         ! Give the fragment positions a random value that is scaled with fragment mass so that the more massive bodies tend to be closer to the impact point
         ! Later, velocities will be scaled such that the farther away a fragment is placed from the impact point, the higher will its velocity be.
         call random_number(mass_rscale)
         mass_rscale(:) = (mass_rscale(:) + 1.0_DP) / 2
         mass_rscale(:) = mass_rscale(:) * (fragments%mtot / fragments%mass(:))**(0.125_DP) ! The power is arbitrary. It just gives the velocity a small mass dependence
         mass_rscale(:) = mass_rscale(:) / maxval(mass_rscale(:))

         do loop = 1, MAXLOOP
            if (.not.any(loverlap(:))) exit
            if (lhitandrun) then
               fragment_cloud_radius(:) = impactors%radius(:)
               fragment_cloud_center(:,1) = impactors%rc(:,1) 
               fragment_cloud_center(:,2) = impactors%rc(:,2) + rdistance * impactors%bounce_unit(:)
            else if (lsupercat) then
               fragment_cloud_center(:,1) = impactors%rc(:,1) 
               fragment_cloud_center(:,2) = impactors%rc(:,2) 
               fragment_cloud_radius(:) = cloud_size_scale_factor * rdistance 
            else
               fragment_cloud_center(:,1) = impactors%rbimp(:) - impactors%radius(1) * impactors%y_unit(:)
               fragment_cloud_center(:,2) = impactors%rbimp(:) + impactors%radius(2) * impactors%y_unit(:)
               fragment_cloud_radius(:) = cloud_size_scale_factor * rdistance
            end if
            if (lsupercat) then
               istart = 1
            else
               fragments%rc(:,1) = fragment_cloud_center(:,1)
               istart = 2
            end if

            do concurrent(i = istart:nfrag, loverlap(i))
               ! Make a random cloud
               call random_number(fragments%rc(:,i))

               ! Make the fragment cloud symmertic about 0
               fragments%rc(:,i) = 2 *(fragments%rc(:,i) - 0.5_DP)

               j = fragments%origin_body(i)

               ! Scale the cloud size
               fragments%rc(:,i) = fragment_cloud_radius(j) * mass_rscale(:) * .unit.fragments%rc(:,i)

               ! Shift to the cloud center coordinates
               fragments%rc(:,i) = fragments%rc(:,i) + fragment_cloud_center(:,j)

               ! Make sure that the fragments are positioned away from the impact point
               direction = dot_product(fragments%rc(:,i) - impactors%rbimp(:), fragment_cloud_center(:,j) - impactors%rbimp(:))
               if (direction < 0.0_DP) then
                  fragments%rc(:,i) = fragments%rc(:,i) - fragment_cloud_center(:,j)
                  fragments%rc(:,i) = -fragments%rc(:,i) + fragment_cloud_center(:,j)
               end if

            end do

            ! Check for any overlapping bodies.
            loverlap(:) = .false.
            do j = 1, nfrag
               do i = j + 1, nfrag
                  dis = .mag.(fragments%rc(:,j) - fragments%rc(:,i))
                  loverlap(i) = loverlap(i) .or. (dis <= (fragments%radius(i) + fragments%radius(j))) 
               end do
            end do
            rdistance = rdistance * collider%fail_scale
         end do

         call collision_util_shift_vector_to_origin(fragments%mass, fragments%rc)
         call collider%set_coordinate_system()

         do concurrent(i = 1:nfrag)
            fragments%rb(:,i) = fragments%rc(:,i) + impactors%rbcom(:)
         end do

      end associate

      return
   end subroutine fraggle_generate_pos_vec


   module subroutine fraggle_generate_rot_vec(collider, nbody_system, param)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Calculates the spins of a collection of fragments such that they conserve angular momentum 
      implicit none
      ! Arguments
      class(collision_fraggle),     intent(inout) :: collider     !! Fraggle collision system object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      ! Internals
      real(DP), dimension(NDIM) :: Lbefore, Lafter, L_spin, rotdir
      real(DP) :: v_init, v_final, mass_init, mass_final, rotmag
      real(DP), parameter :: random_scale_factor = 0.01_DP !! The relative scale factor to apply to the random component of the rotation vector
      integer(I4B) :: i

      associate(fragments => collider%fragments, impactors => collider%impactors, nfrag => collider%fragments%nbody)

         ! Torque the first body based on the change in angular momentum betwen the pre- and post-impact system assuming a simple bounce
         mass_init = impactors%mass(2)
         mass_final = sum(fragments%mass(2:nfrag))
         v_init = .mag.(impactors%vb(:,2) - impactors%vb(:,1))
         v_final = sqrt(mass_init / mass_final * v_init**2 - 2 * impactors%Qloss / mass_final)

         Lbefore(:) = mass_init * (impactors%rb(:,2) - impactors%rb(:,1)) .cross. (impactors%vb(:,2) - impactors%vb(:,1))
          
         Lafter(:) = mass_final * (impactors%rb(:,2) - impactors%rb(:,1)) .cross. (v_final * impactors%bounce_unit(:))
         L_spin(:) = impactors%L_spin(:,1) + (Lbefore(:) - Lafter(:))
         fragments%rot(:,1) = L_spin(:) / (fragments%mass(1) * fragments%radius(1)**2 * fragments%Ip(3,1))

         ! Add in some random spin noise. The magnitude will be scaled by the before-after amount and the direction will be random
         do concurrent(i = 2:nfrag)
            call random_number(rotdir)
            rotdir = rotdir - 0.5_DP
            rotdir = .unit. rotdir
            call random_number(rotmag)
            rotmag = random_scale_factor * rotmag * .mag. (Lbefore(:) - Lafter(:)) / ((nfrag - 1) * fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i))
            fragments%rot(:,i) = rotmag * rotdir
         end do

      end associate

      return
   end subroutine fraggle_generate_rot_vec


   module subroutine fraggle_generate_vel_vec(collider, nbody_system, param, lfailure)
      !! Author:  David A. Minton
      !!
      !! Generates an initial velocity distribution. For disruptions, the velocity magnitude is set to be
      !! 2x the escape velocity of the colliding pair. For hit and runs the velocity magnitude is set to be
      !! 2x the escape velocity of the smallest of the two bodies.
      implicit none
      ! Arguments
      class(collision_fraggle),     intent(inout) :: collider     !! Fraggle collision system object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      logical,                      intent(out)   :: lfailure     !! Did the velocity computation fail?
      ! Internals
      integer(I4B) :: i, j, loop, try, istart, nfrag
      logical :: lhitandrun, lsupercat
      real(DP), dimension(NDIM) :: vimp_unit, rimp, vrot, L_residual
      real(DP) :: vmag, vesc,  E_residual, ke_avail, ke_remove, E_residual_min, fscale, f_spin, f_orbit
      integer(I4B), dimension(collider%fragments%nbody) :: vsign
      real(DP), dimension(collider%fragments%nbody) :: vscale, ke_rot_remove
      real(DP), parameter     :: vmin_initial_factor = 1.5_DP ! For the initial "guess" of fragment velocities, this is the maximum velocity relative to escape velocity that the fragments will have
      real(DP), parameter     :: vmax_initial_factor = 5.0_DP ! For the initial "guess" of fragment velocities, this is the maximum velocity relative to escape velocity that the fragments will have
      integer(I4B), parameter :: MAXLOOP = 100
      integer(I4B), parameter :: MAXTRY  = 20
      real(DP), parameter :: TOL = 1.0_DP
      class(collision_fraggle), allocatable :: collider_local
      character(len=STRMAX) :: message

      associate(impactors => collider%impactors)
         nfrag = collider%fragments%nbody
         lhitandrun = (impactors%regime == COLLRESOLVE_REGIME_HIT_AND_RUN) 
         lsupercat = (impactors%regime == COLLRESOLVE_REGIME_SUPERCATASTROPHIC) 

         allocate(collider_local, source=collider)
         associate(fragments => collider_local%fragments)

            ! The fragments will be divided into two "clouds" based on identified origin body. 
            ! These clouds will collectively travel like two impactors bouncing off of each other. 
            where(fragments%origin_body(:) == 1)
               vsign(:) = -1
            elsewhere
               vsign(:) = 1
            end where

            ! The minimum fragment velocity will be set by the escape velocity
            if (lhitandrun) then
               vesc = sqrt(2 * impactors%Gmass(2) / impactors%radius(2))
            else
               vesc = sqrt(2 * sum(impactors%Gmass(:)) / sum(impactors%radius(:)))
            end if

            ! Scale the magnitude of the velocity by the distance from the impact point
            ! This will reduce the chances of fragments colliding with each other immediately, and is more physically correct  
            do concurrent(i = 1:nfrag)
               rimp(:) = fragments%rc(:,i) - impactors%rbimp(:) 
               vscale(i) = .mag. rimp(:) / (.mag. (impactors%rb(:,2) - impactors%rb(:,1)))
            end do

            ! Set the velocity scale factor to span from vmin/vesc to vmax/vesc
            vscale(:) = vscale(:)/minval(vscale(:))
            fscale = log(vmax_initial_factor - vmin_initial_factor + 1.0_DP) / log(maxval(vscale(:)))
            vscale(:) = vscale(:)**fscale + vmin_initial_factor - 1.0_DP

            ! Set the velocities of all fragments using all of the scale factors determined above
            do concurrent(i = 1:nfrag)
               j = fragments%origin_body(i)
               vrot(:) = impactors%rot(:,j) .cross. (fragments%rc(:,i) - impactors%rc(:,j))
               if (lhitandrun) then
                  if (i == 1) then
                     fragments%vc(:,1) = impactors%vc(:,1)
                  else
                     vmag = .mag.impactors%vc(:,2) / maxval(vscale(:))
                     fragments%vc(:,i) = vmag * vscale(i) * impactors%bounce_unit(:) * vsign(i) + vrot(:)
                  end if
               else
                  ! Add more velocity dispersion to disruptions vs hit and runs.
                  vmag = vesc * vscale(i) 
                  rimp(:) = fragments%rc(:,i) - impactors%rbimp(:)
                  vimp_unit(:) = .unit. rimp(:)
                  fragments%vc(:,i) = vmag * vscale(i) * vimp_unit(:) + vrot(:) 
               end if
            end do

            ! Hit and run collisions should only affect the runner
            if (lhitandrun) then
               istart = 2
            else 
               istart = 1
            end if

            ! Every time the collision-frame velocities are altered, we need to be sure to shift everything back to the center-of-mass frame
            call collision_util_shift_vector_to_origin(fragments%mass, fragments%vc)            
            call fragments%set_coordinate_system()
            E_residual_min = huge(1.0_DP)
            lfailure = .false.
            outer: do try = 1, MAXTRY
               do loop = 1, MAXLOOP
                  call collider_local%get_energy_and_momentum(nbody_system, param, phase="after")
                  ! Check for any residual angular momentum, and if there is any, put it into spin of the largest body
                  E_residual = collider_local%te(2) - collider_local%te(1) 
                  E_residual = E_residual + impactors%Qloss
                  L_residual(:) = collider_local%L_total(:,2) - collider_local%L_total(:,1)
                  fragments%L_spin(:,1) = fragments%L_spin(:,1) - L_residual(:) 
                  fragments%rot(:,1) = fragments%L_spin(:,1) / (fragments%mass(1) * fragments%radius(1)**2 * fragments%Ip(:,1)) 

                  call collider_local%get_energy_and_momentum(nbody_system, param, phase="after")
                  E_residual = collider_local%te(2) + impactors%Qloss - collider_local%te(1) 
                  L_residual(:) = (collider_local%L_total(:,2) - collider_local%L_total(:,1)) 
                  if ((abs(E_residual) < abs(E_residual_min)) .or. ((E_residual <= 0.0_DP) .and. (E_residual_min >= 0.0_DP))) then ! This is our best case so far. Save it for posterity
                     if (allocated(collider%fragments)) deallocate(collider%fragments)
                     allocate(collider%fragments, source=fragments)
                     E_residual_min = E_residual
                     if (abs(E_residual) <= tiny(0.0_DP)) exit outer
                     if ((E_residual < 0.0_DP) .and. (impactors%Qloss / abs(E_residual)) > TOL) exit outer
                  end if

                  ke_avail = fragments%ke_orbit_tot - 0.5_DP * fragments%mtot * vesc**2

                  ! Remove a constant amount of velocity from the bodies so we don't shift the center of mass and screw up the momentum 
                  f_spin = (fragments%ke_spin_tot )/ (fragments%ke_spin_tot + fragments%ke_orbit_tot)
                  f_orbit = 1.0_DP - f_spin

                  ke_remove = min(f_orbit * E_residual, ke_avail)
                  fscale = sqrt((fragments%ke_orbit_tot - f_orbit * ke_remove)/fragments%ke_orbit_tot)
                  fragments%vc(:,:) = fscale * fragments%vc(:,:)

                  ke_remove = min(E_residual - ke_remove, fragments%ke_spin_tot)
                  ke_rot_remove(:) = ke_remove * (fragments%ke_spin(:) / fragments%ke_spin_tot)
                  where(ke_rot_remove(:) > fragments%ke_spin(:)) ke_rot_remove(:) = fragments%ke_spin(:) 
                  do concurrent(i = 1:nfrag, fragments%ke_spin(i) > epsilon(1.0_DP))
                     fscale = sqrt((fragments%ke_spin(i) - ke_rot_remove(i))/fragments%ke_spin(i))
                     fragments%rot(:,i) = fscale * fragments%rot(:,i)
                  end do

                  ! Update the unit vectors and magnitudes for the fragments based on their new orbits and rotations
                  call collision_util_shift_vector_to_origin(fragments%mass, fragments%vc)            
                  call fragments%set_coordinate_system()

               end do
               ! We didn't converge. Try another configuration and see if we get a better result
               call fraggle_generate_pos_vec(collider_local)
               call fraggle_generate_rot_vec(collider_local, nbody_system, param)
               collider_local%fail_scale = collider_local%fail_scale*1.01_DP
            end do outer
            lfailure = E_residual < 0.0_DP
            if (lfailure) then
               write(message,*) "Fraggle failed to converge after ",try*loop," steps. This collision may have added energy."
            else 
               write(message,*) "Fraggle succeeded after ",try*loop," steps."
            end if
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)

         end associate
      end associate
      return
   end subroutine fraggle_generate_vel_vec


end submodule s_fraggle_generate
