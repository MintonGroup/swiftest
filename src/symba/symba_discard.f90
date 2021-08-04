submodule (symba_classes) s_symba_discard
   use swiftest
contains

   subroutine symba_discard_cb_pl(pl, system, param)
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
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i, j
      real(DP)     :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2
   
      associate(npl => pl%nbody, cb => system%cb)
         call system%set_msys()
         rmin2 = param%rmin**2
         rmax2 = param%rmax*2
         rmaxu2 = param%rmaxu**2
         do i = 1, npl
            if (pl%status(i) == ACTIVE) then
               rh2 = dot_product(pl%xh(:,i), pl%xh(:,i))
               if ((param%rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
                  pl%ldiscard(i) = .true.
                  pl%lcollision(i) = .false. 
                  pl%status(i) = DISCARDED_RMAX
                  write(*, *) "Massive body ",  pl%id(i), " too far from the central body at t = ", param%t
               else if ((param%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
                  pl%ldiscard(i) = .true.
                  pl%lcollision(i) = .false. 
                  pl%status(i) = DISCARDED_RMIN
                  write(*, *) "Massive body ", pl%id(i), " too close to the central body at t = ", param%t
               else if (param%rmaxu >= 0.0_DP) then
                  rb2 = dot_product(pl%xb(:,i), pl%xb(:,i))
                  vb2 = dot_product(pl%vb(:,i), pl%vb(:,i))
                  energy = 0.5_DP * vb2 - system%Gmtot / sqrt(rb2)
                  if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
                     pl%ldiscard(i) = .true.
                     pl%lcollision(i) = .false. 
                     pl%status(i) = DISCARDED_RMAXU
                     write(*, *) "Massive body ", pl%id(i), " is unbound and too far from barycenter at t = ", param%t
                  end if
               end if
            end if
         end do
      end associate
   
      return
   end subroutine symba_discard_cb_pl


   subroutine symba_discard_nonplpl(pl, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if planets should be discarded based on their positions or because they are unbound
      !s
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_massive5.f 
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: pl     !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
    
      ! First check for collisions with the central body
      associate(npl => pl%nbody, cb => system%cb)
         if (npl == 0) return 
         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or. &
             (param%rmaxu >= 0.0_DP) .or. ((param%qmin >= 0.0_DP) .and. (param%qmin_coord == "BARY"))) then
            call pl%h2b(cb) 
         end if
         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or.  (param%rmaxu >= 0.0_DP)) then
            call symba_discard_cb_pl(pl, system, param)
         end if
         if (param%qmin >= 0.0_DP .and. npl > 0) call symba_discard_peri_pl(pl, system, param)
      end associate

      return
   end subroutine symba_discard_nonplpl


   subroutine symba_discard_peri_pl(pl, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_peri_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_mass_peri.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: pl     !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      logical, save      :: lfirst = .true.
      logical            :: lfirst_orig
      integer(I4B)       :: i


      lfirst_orig = pl%lfirst
      pl%lfirst = lfirst
      if (lfirst) then
         call pl%get_peri(system, param)          
         lfirst = .false.
      else
         call pl%get_peri(system, param)          
         do i = 1, pl%nbody
            if (pl%status(i) == ACTIVE) then
               if ((pl%isperi(i) == 0) .and. (pl%nplenc(i)== 0)) then
                  if ((pl%atp(i) >= param%qmin_alo) .and. (pl%atp(i) <= param%qmin_ahi) .and. (pl%peri(i) <= param%qmin)) then
                     pl%ldiscard(i) = .true.
                     pl%lcollision(i) = .false.
                     pl%status(i) = DISCARDED_PERI
                     write(*, *) "Particle ", pl%id(i), " perihelion distance too small at t = ", param%t
                  end if
               end if
            end if
         end do
      end if
      pl%lfirst = lfirst_orig
   
      return
   
   end subroutine symba_discard_peri_pl


   module subroutine symba_discard_pl(self, system, param)
      !! author: David A. Minton
      !!
      !! Call the various flavors of discards for massive bodies in SyMBA runs, including discards due to colling with the central body, 
      !! escaping the system, or colliding with each other.
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
   
      select type(system)
      class is (symba_nbody_system)
         call symba_discard_nonplpl(self, system, param)
         call system%plplenc_list%scrub_non_collision(system, param)
      end select

      return
   end subroutine symba_discard_pl

end submodule s_symba_discard