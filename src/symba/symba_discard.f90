submodule (symba_classes) s_symba_discard
   use swiftest
contains

   module subroutine symba_discard_pl(self, system, param)
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
    
      ! First check for collisions with the central body
      associate(pl => self, npl => self%nbody, cb => system%cb)
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

      select type(param)
      class is (symba_parameters)
         if (param%lfragmentation) then

         end if

      end select

      return
   end subroutine symba_discard_pl


   subroutine symba_discard_cb_pl(pl, system, param)
      !! author: David A. Minton
      !!
      !!  Check to see if planets should be discarded based on their positions relative to the central body 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_cb.f90
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
                  pl%status(i) = DISCARDED_RMAX
                  write(*, *) "Massive body ",  pl%id(i), " too far from the central body at t = ", param%t
               else if ((param%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
                  pl%ldiscard(i) = .true.
                  pl%status(i) = DISCARDED_RMIN
                  write(*, *) "Massive body ", pl%id(i), " too close to the central body at t = ", param%t
               else if (param%rmaxu >= 0.0_DP) then
                  rb2 = dot_product(pl%xb(:,i), pl%xb(:,i))
                  vb2 = dot_product(pl%vb(:,i), pl%vb(:,i))
                  energy = 0.5_DP * vb2 - system%msys / sqrt(rb2)
                  if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
                     pl%ldiscard(i) = .true.
                     pl%status(i) = DISCARDED_RMAXU
                     write(*, *) "Massive body ", pl%id(i), " is unbound and too far from barycenter at t = ", param%t
                  end if
               end if
            end if
         end do
      end associate
   
      return
   end subroutine symba_discard_cb_pl


   subroutine symba_discard_peri_pl(pl, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_peri.f90
      !! Adapted from Hal Levison's Swift routine discard_peri.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: pl     !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      return
   end subroutine symba_discard_peri_pl

end submodule s_symba_discard