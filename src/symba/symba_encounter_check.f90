submodule (symba_classes) s_symba_encounter_check
   use swiftest
contains
   module function symba_encounter_check_pl(self, system, dt, irec) result(lany_encounter)
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout)  :: self           !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout)  :: system         !! SyMBA nbody system object
      real(DP),                  intent(in)     :: dt             !! step size
      integer(I4B),              intent(in)     :: irec           !! Current recursion level
      ! Result
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      ! Internals
      real(DP)                                  :: r2crit, vdotr, r2, v2, tmin, r2min, term2
      integer(I4B)                              :: j, nenc_old
      integer(I8B)                              :: k
      real(DP),     dimension(NDIM)             :: xr, vr
      integer(I4B), dimension(:,:), allocatable :: ind
      logical,      dimension(:),   allocatable :: lencounter, loc_lvdotr
   
      associate(pl => self, npl => self%nbody, nplpl => self%nplpl)
         allocate(lencounter(nplpl), loc_lvdotr(nplpl))
         lencounter(:) = .false.
   
         term2 = RHSCALE * (RSHELL**irec)
   
         do k = 1, nplpl
            associate(i => pl%k_plpl(1, k), j => pl%k_plpl(2, k))
               xr(:) = pl%xh(:, j) - pl%xh(:, i)
               r2 = dot_product(xr(:), xr(:)) 
               r2crit = ((pl%rhill(i) + pl%rhill(i)) * term2)**2
               vr(:) = pl%vh(:, j) - pl%vh(:, i)
               vdotr = dot_product(vr(:), xr(:))
               if (r2 < r2crit) then
                  lencounter(k) = .true.
                  loc_lvdotr(k) = (vdotr < 0.0_DP)
               else
                  if (vdotr < 0.0_DP) then
                     v2 = dot_product(vr(:), vr(:))
                     tmin = -vdotr /  v2
                     if (tmin < dt) then
                        r2min = r2 - vdotr * vdotr / v2
                     else
                        r2min = r2 + 2 * vdotr * dt + v2 * dt * dt
                     end if
                     r2min = min(r2min, r2)
                     if (r2min <= r2crit) then
                        lencounter(k) = .true.
                        loc_lvdotr(k) = (vdotr < 0.0_DP)
                     end if
                  end if
               end if
            end associate
         end do

         lany_encounter = any(lencounter(:))
         if (lany_encounter) then 
            associate(plplenc_list => system%plplenc_list, nenc => system%plplenc_list%nenc)
               nenc_old = nenc
               call plplenc_list%resize(nenc_old + count(lencounter(:)))
               plplenc_list%status(nenc_old+1:nenc) = ACTIVE
               plplenc_list%level(nenc_old+1:nenc) = irec
               plplenc_list%lvdotr(nenc_old+1:nenc) = pack(loc_lvdotr(:), lencounter(:))
               plplenc_list%index1(nenc_old+1:nenc) = pack(pl%k_plpl(1,:), lencounter(:))
               plplenc_list%index2(nenc_old+1:nenc) = pack(pl%k_plpl(2,:), lencounter(:))
               pl%lencounter(plplenc_list%index1(nenc_old+1:nenc)) = .true.
               pl%lencounter(plplenc_list%index2(nenc_old+1:nenc)) = .true.
            end associate
         end if
      end associate
      return
   end function symba_encounter_check_pl

   module function symba_encounter_check_tp(self, system, dt) result(lencounter)
      implicit none
      ! Arguments
      class(symba_tp),           intent(inout) :: self       !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      ! Result
      logical                                  :: lencounter !! Returns true if there is at least one close encounter      

      lencounter = .false.
      return
   end function symba_encounter_check_tp

end submodule s_symba_encounter_check