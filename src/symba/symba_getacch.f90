submodule (symba_classes) s_symba_getacch
   use swiftest
contains
   module subroutine symba_getacch_pl(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine symba_getacch.f90
      !! Adapted from Hal Levison's Swift routine symba5_getacch.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)              :: k
      real(DP)                  :: rji2, rlim2, faci, facj
      real(DP), dimension(NDIM) :: dx

      select type(system)
      class is (symba_nbody_system)
         associate(pl => self, cb => system%cb, plplenc_list => system%plplenc_list, nplplenc => system%plplenc_list%nenc)
            call helio_getacch_pl(pl, system, param, t, lbeg)
            ! Remove accelerations from encountering pairs
            do k = 1, nplplenc
               associate(i => plplenc_list%index1(k), j => plplenc_list%index2(k))
                  dx(:) = pl%xh(:, j) - pl%xh(:, i)
                  rji2 = dot_product(dx(:), dx(:))
                  rlim2 = (pl%radius(i) + pl%radius(j))**2
                  if (rji2 > rlim2) then
                     irij3 = 1.0_DP / (rji2 * sqrt(rji2))
                     faci = pl%Gmass(i) * irij3
                     facj = pl%Gmass(j) * irij3
                     pl%ah(:, i) = pl%ah(:, i) - facj * dx(:)
                     pl%ah(:, j) = pl%ah(:, j) + faci * dx(:)
                  end if
               end associate
            end do
         end associate
      end select

      return
      end subroutine symba_getacch_pl

      module subroutine symba_getacch_tp(self, system, param, t, lbeg)
         !! author: David A. Minton
         !!
         !! Compute heliocentric accelerations of test particles
         !!
         !! Adapted from David E. Kaufmann's Swifter routine symba_getacch_tp.f90
         !! Adapted from Hal Levison's Swift routine symba5_getacch.f
         implicit none
         ! Arguments
         class(symba_tp),              intent(inout) :: self   !! SyMBA test particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
         ! Internals
         integer(I4B)              :: k
         real(DP)                  :: rji2, fac, rlim2
         real(DP), dimension(NDIM) :: dx
     
         select type(system)
         class is (symba_nbody_system)
            associate(tp => self, cb => system%cb, pl => system%pl, pltpenc_list => system%pltpenc_list, npltpenc => system%pltpenc_list%nenc)
               call helio_getacch_tp(tp, system, param, t, lbeg)
               ! Remove accelerations from encountering pairs
               do k = 1, npltpenc
                  associate(i => pltpenc_list%index1(k), j => pltpenc_list%index2(k))
                     if (tp%status(j) == ACTIVE) THEN
                        dx(:) = tp%xh(:,j) - pl%xh(:,i)
                        rji2 = dot_product(dx(:), dx(:))
                        rlim2 = (pl%radius(i))**2
                        if (rji2 > rlim2) then
                           fac = pl%Gmass(i) / (rji2 * sqrt(rji2))
                           tp%ah(:,j) = tp%ah(:,j) + fac * dx(:)
                        end if
                     end IF
                  end associate
               end DO
            end associate
         end select
         return
      end subroutine symba_getacch_tp

end submodule s_symba_getacch
