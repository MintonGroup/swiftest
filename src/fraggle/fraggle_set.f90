submodule(fraggle_classes) s_fraggle_set
   use swiftest
contains

   module subroutine fraggle_set_budgets_fragments(self, colliders)
      !! author: David A. Minton
      !!
      !! Sets the energy and momentum budgets of the fragments based on the collider values and the before/after values of energy and momentum
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: self      !! Fraggle fragment system object
      class(fraggle_colliders), intent(inout) :: colliders !! Fraggle collider system object
      ! Internals
      real(DP) :: dEtot
      real(DP), dimension(NDIM) :: dL

      associate(frag => self)

         dEtot = frag%Etot_after - frag%Etot_before 
         dL(:) = frag%Ltot_after(:) - frag%Ltot_before(:)

         frag%L_budget(:) = -dL(:)
         frag%ke_budget = -(dEtot - 0.5_DP * frag%mtot * dot_product(frag%vbcom(:), frag%vbcom(:))) - frag%Qloss 

      end associate
      return
   end subroutine fraggle_set_budgets_fragments


   module subroutine fraggle_set_mass_dist_fragments(self, colliders, param)
      !! author: David A. Minton
      !!
      !! Sets the mass of fragments based on the mass distribution returned by the regime calculation.
      !! This subroutine must be run after the the setup rourtine has been run on the fragments
      implicit none
      ! Arguments
      class(fraggle_fragments),     intent(inout) :: self      !! Fraggle fragment system object
      class(fraggle_colliders),     intent(inout) :: colliders !! Fraggle collider system object
      class(swiftest_parameters),   intent(in)    :: param     !! Current Swiftest run configuration parameters
      ! Internals
      integer(I4B)              :: i, jproj, jtarg, nfrag, istart
      real(DP), dimension(2)    :: volume
      real(DP), dimension(NDIM) :: Ip_avg
      real(DP) :: mfrag, mremaining, min_mfrag
      real(DP), parameter :: BETA = 2.85_DP
      integer(I4B), parameter :: NFRAGMAX = 100  !! Maximum number of fragments that can be generated
      integer(I4B), parameter :: NFRAGMIN = 7 !! Minimum number of fragments that can be generated (set by the fraggle_generate algorithm for constraining momentum and energy)
      integer(I4B), parameter :: NFRAG_SIZE_MULTIPLIER = 3 !! Log-space scale factor that scales the number of fragments by the collisional system mass
     
      associate(frag => self)
         ! Get mass weighted mean of Ip and density
         volume(1:2) = 4._DP / 3._DP * PI * colliders%radius(1:2)**3
         Ip_avg(:) = (colliders%mass(1) * colliders%Ip(:,1) + colliders%mass(2) * colliders%Ip(:,2)) / frag%mtot
         if (colliders%mass(1) > colliders%mass(2)) then
            jtarg = 1
            jproj = 2
         else
            jtarg = 2
            jproj = 1
         end if

         select type(param)
         class is (symba_parameters)
            min_mfrag = (param%min_GMfrag / param%GU) 
            nfrag = ceiling(NFRAG_SIZE_MULTIPLIER  * log(frag%mtot / min_mfrag))
            nfrag = max(min(nfrag, NFRAGMAX), NFRAGMIN)
         class default
            min_mfrag = 0.0_DP
            nfrag = NFRAGMAX
         end select
  
         select case(frag%regime)
         case(COLLRESOLVE_REGIME_DISRUPTION, COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
            istart = 2
         case(COLLRESOLVE_REGIME_HIT_AND_RUN)
            istart = 1
         case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE) 
            call frag%setup(1, param)
            frag%mass(1) = frag%mass_dist(1)
            frag%radius(1) = colliders%radius(jtarg)
            frag%density(1) = frag%mass_dist(1) / volume(jtarg)
            frag%Ip(:, 1) = colliders%Ip(:,1)
            return
         case default
            write(*,*) "fraggle_set_mass_dist_fragments error: Unrecognized regime code",frag%regime
         end select

         i = istart 
         mremaining = max(frag%mtot - sum(frag%mass_dist(1:istart)), 0.0_DP)
         do while (i <= NFRAGMAX)
            mfrag = (i - istart + 1)**(-3._DP / BETA) * frag%mass_dist(istart)
            if (mremaining - mfrag < 0.0_DP) exit
            mremaining = mremaining - mfrag
            i = i + 1
         end do

         call frag%setup(nfrag, param)
         frag%mass(1:istart) = frag%mass_dist(1:istart)
         mremaining = max(frag%mtot - sum(frag%mass_dist(1:istart)), 0.0_DP)
         do i = istart + 1, nfrag
            mfrag = (i - istart + 1)**(-3._DP / BETA) * frag%mass_dist(istart)
            mfrag = min(mfrag, mremaining)
            frag%mass(i) = mfrag
            mremaining = mremaining - mfrag
         end do
         select case(frag%regime)
         case(COLLRESOLVE_REGIME_HIT_AND_RUN) 
            frag%mass(1) = frag%mass_dist(1)
            frag%radius(1) = colliders%radius(jtarg)
            frag%density(1) = frag%mass_dist(1) / volume(jtarg)
            frag%Ip(:, 1) = colliders%Ip(:,1)
            istart = 2
         case default
            istart = 1
         end select
         if (mremaining > 0.0_DP) then
            ! Distribute remaining mass among the fragments
            mfrag = 1._DP + mremaining / sum(frag%mass(istart:nfrag))
            do i = istart, nfrag
               frag%mass(i) = frag%mass(i) * mfrag
            end do
            mremaining = frag%mtot - sum(frag%mass(1:nfrag))
            frag%mass(nfrag) = frag%mass(nfrag) + mremaining
         end if
         frag%density(istart:nfrag) = frag%mtot / sum(volume(:))
         frag%radius(istart:nfrag) = (3 * frag%mass(istart:nfrag) / (4 * PI * frag%density(istart:nfrag)))**(1.0_DP / 3.0_DP)
         do i = istart, nfrag
            frag%Ip(:, i) = Ip_avg(:)
         end do

      end associate

      return
   end subroutine fraggle_set_mass_dist_fragments


   module subroutine fraggle_set_coordinate_system(self, colliders)
      !! author: David A. Minton
      !!
      !! Defines the collisional coordinate system, including the unit vectors of both the system and individual fragments.
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: self      !! Fraggle fragment system object
      class(fraggle_colliders), intent(inout) :: colliders !! Fraggle collider system object
      ! Internals
      integer(I4B) :: i
      real(DP), dimension(NDIM) ::  x_cross_v, delta_r, delta_v, Ltot
      real(DP)   :: r_col_norm, v_col_norm
      real(DP), dimension(NDIM, self%nbody) :: L_sigma

      associate(frag => self, nfrag => self%nbody)
         delta_v(:) = colliders%vb(:, 2) - colliders%vb(:, 1)
         v_col_norm = .mag. delta_v(:)
         delta_r(:) = colliders%xb(:, 2) - colliders%xb(:, 1)
         r_col_norm = .mag. delta_r(:)
   
         ! We will initialize fragments on a plane defined by the pre-impact system, with the z-axis aligned with the angular momentum vector
         ! and the y-axis aligned with the pre-impact distance vector.
         Ltot = colliders%L_orbit(:,1) + colliders%L_orbit(:,2) + colliders%L_spin(:,1) + colliders%L_spin(:,2)
         frag%y_coll_unit(:) = delta_r(:) / r_col_norm 
         frag%z_coll_unit(:) = Ltot(:) / (.mag. Ltot(:))
         ! The cross product of the y- by z-axis will give us the x-axis
         frag%x_coll_unit(:) = frag%y_coll_unit(:) .cross. frag%z_coll_unit(:)
   
         if (.not.any(frag%x_coll(:,:) > 0.0_DP)) return
         frag%rmag(:) = .mag. frag%x_coll(:,:)
   
         call random_number(L_sigma(:,:)) ! Randomize the tangential velocity direction. This helps to ensure that the tangential velocity doesn't completely line up with the angular momentum vector,
                                          ! otherwise we can get an ill-conditioned system
         do concurrent(i = 1:nfrag, frag%rmag(i) > 0.0_DP)
            frag%v_r_unit(:, i) = frag%x_coll(:, i) / frag%rmag(i)
            frag%v_n_unit(:, i) = frag%z_coll_unit(:) + 2e-1_DP * (L_sigma(:,i) - 0.5_DP)
            frag%v_n_unit(:, i) = frag%v_n_unit(:, i) / (.mag. frag%v_n_unit(:, i))
            frag%v_t_unit(:, i) = frag%v_n_unit(:, i) .cross. frag%v_r_unit(:, i)
            frag%v_t_unit(:, i) = frag%v_t_unit(:, i) / (.mag. frag%v_t_unit(:, i))
         end do
      end associate

      return
   end subroutine fraggle_set_coordinate_system


   ! module subroutine symba_set_collresolve_colliders(self, cb, pl, idx)
   !    !! author: David A. Minton
   !    !!
   !    !! Calculate the two-body equivalent values given a set of input collider indices
   !    use swiftest_classes, only : swiftest_nbody_system
   !    implicit none
   !    ! Arguments
   !    class(fraggle_colliders),               intent(inout) :: self !! Fraggle collider object
   !    class(symba_cb),                     intent(in)    :: cb   !! Swiftest central body object system object
   !    class(symba_pl),                     intent(in)    :: pl   !! Swiftest central body object system object
   !    integer(I4B),             dimension(:), intent(in)    :: idx  !! Index array of bodies from the pl object to use to calculate a "two-body equivalent" collisional pair
   !    ! Internals
   !    real(DP), dimension(NDIM, 2)  :: mxc, vc
   !    real(DP), dimension(NDIM) :: vcom, xcom

   !    associate(colliders => self)

   !       ! Compute orbital angular momentum of pre-impact system
   !       xcom(:) = (colliders%mass(1) * colliders%xb(:, 1) + colliders%mass(2) * colliders%xb(:, 2)) / sum(colliders%mass(:))
   !       vcom(:) = (colliders%mass(1) * colliders%vb(:, 1) + colliders%mass(2) * colliders%vb(:, 2)) / sum(colliders%mass(:))
   !       mxc(:, 1) = colliders%mass(1) * (colliders%xb(:, 1) - xcom(:))
   !       mxc(:, 2) = colliders%mass(2) * (colliders%xb(:, 2) - xcom(:))
   !       vc(:, 1) = colliders%vb(:, 1) - vcom(:)
   !       vc(:, 2) = colliders%vb(:, 2) - vcom(:)

   !       colliders%L_orbit(:,:) = mxc(:,:) .cross. vc(:,:)

   !    end associate

   !    return
   ! end subroutine symbe_set_collresolve_colliders


   module subroutine fraggle_set_natural_scale_factors(self, colliders)
      !! author: David A. Minton
      !!
      !! Scales dimenional quantities to ~O(1) with respect to the collisional system. 
      !! This scaling makes it easier for the non-linear minimization to converge on a solution
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: self      !! Fraggle fragment system object
      class(fraggle_colliders), intent(inout) :: colliders !! Fraggle collider system object
      ! Internals
      integer(I4B) :: i

      associate(frag => self)
         ! Find the center of mass of the collisional system	
         frag%xbcom(:) = (colliders%mass(1) * colliders%xb(:,1) + colliders%mass(2) * colliders%xb(:,2)) / frag%mtot
         frag%vbcom(:) = (colliders%mass(1) * colliders%vb(:,1) + colliders%mass(2) * colliders%vb(:,2)) / frag%mtot

         ! Set scale factors
         frag%Escale = 0.5_DP * (colliders%mass(1) * dot_product(colliders%vb(:,1), colliders%vb(:,1)) + colliders%mass(2) * dot_product(colliders%vb(:,2), colliders%vb(:,2)))
         frag%dscale = sum(colliders%radius(:))
         frag%mscale = frag%mtot 
         frag%vscale = sqrt(frag%Escale / frag%mscale) 
         frag%tscale = frag%dscale / frag%vscale 
         frag%Lscale = frag%mscale * frag%dscale * frag%vscale

         ! Scale all dimensioned quantities of colliders and fragments
         frag%xbcom(:) = frag%xbcom(:) / frag%dscale
         frag%vbcom(:) = frag%vbcom(:) / frag%vscale
         colliders%xb(:,:) = colliders%xb(:,:) / frag%dscale
         colliders%vb(:,:) = colliders%vb(:,:) / frag%vscale
         colliders%mass(:) = colliders%mass(:) / frag%mscale
         colliders%radius(:) = colliders%radius(:) / frag%dscale
         colliders%L_spin(:,:) = colliders%L_spin(:,:) / frag%Lscale

         do i = 1, 2
            colliders%rot(:,i) = colliders%L_spin(:,i) / (colliders%mass(i) * colliders%radius(i)**2 * colliders%Ip(3, i))
         end do

         frag%mtot = frag%mtot / frag%mscale
         frag%mass = frag%mass / frag%mscale
         frag%radius = frag%radius / frag%dscale
         frag%Qloss = frag%Qloss / frag%Escale
      end associate

      return
   end subroutine fraggle_set_natural_scale_factors


   module subroutine fraggle_set_original_scale_factors(self, colliders)
      !! author: David A. Minton
      !!
      !! Restores dimenional quantities back to the system units
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: self      !! Fraggle fragment system object
      class(fraggle_colliders), intent(inout) :: colliders !! Fraggle collider system object
      ! Internals
      integer(I4B) :: i
      logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      call ieee_set_halting_mode(IEEE_ALL,.false.)

      associate(frag => self)

         ! Restore scale factors
         frag%xbcom(:) = frag%xbcom(:) * frag%dscale
         frag%vbcom(:) = frag%vbcom(:) * frag%vscale
   
         colliders%mass = colliders%mass * frag%mscale
         colliders%radius = colliders%radius * frag%dscale
         colliders%xb = colliders%xb * frag%dscale
         colliders%vb = colliders%vb * frag%vscale
         colliders%L_spin = colliders%L_spin * frag%Lscale
         do i = 1, 2
            colliders%rot(:,i) = colliders%L_spin(:,i) * (colliders%mass(i) * colliders%radius(i)**2 * colliders%Ip(3, i))
         end do
         frag%Qloss = frag%Qloss * frag%Escale
   
         frag%mtot = frag%mtot * frag%mscale
         frag%mass = frag%mass * frag%mscale
         frag%radius = frag%radius * frag%dscale
         frag%rot = frag%rot / frag%tscale
         frag%x_coll = frag%x_coll * frag%dscale
         frag%v_coll = frag%v_coll * frag%vscale
   
         do i = 1, frag%nbody
            frag%xb(:, i) = frag%x_coll(:, i) + frag%xbcom(:)
            frag%vb(:, i) = frag%v_coll(:, i) + frag%vbcom(:)
         end do
   
         frag%mscale = 1.0_DP
         frag%dscale = 1.0_DP
         frag%vscale = 1.0_DP
         frag%tscale = 1.0_DP
         frag%Lscale = 1.0_DP
         frag%Escale = 1.0_DP
      end associate
      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)
   
      return
   end subroutine fraggle_set_original_scale_factors


end submodule s_fraggle_set