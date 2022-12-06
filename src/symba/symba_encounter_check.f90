!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (symba_classes) s_symba_encounter_check
   use swiftest
contains

   module function symba_encounter_check_pl(self, param, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between massive bodies.
      !!
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout)  :: self   !! SyMBA test particle object  
      class(swiftest_parameters), intent(inout)  :: param  !! Current swiftest run configuration parameters
      class(symba_nbody_system),  intent(inout)  :: system !! SyMBA nbody system object
      real(DP),                   intent(in)     :: dt     !! step size
      integer(I4B),               intent(in)     :: irec   !! Current recursion level
      ! Result
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      ! Internals
      integer(I8B) :: k, nenc
      integer(I4B) :: i, j, npl, nplm, nplt
      logical, dimension(:), allocatable :: lvdotr
      integer(I4B), dimension(:), allocatable :: index1, index2
 
      lany_encounter = .false.
      if (self%nbody == 0) return

      associate(pl => self, plplenc_list => system%plplenc_list, cb => system%cb)

         npl = pl%nbody
         nplm = pl%nplm
         nplt = npl - nplm

         call pl%set_renc(irec)

         if (nplt == 0) then
            call encounter_check_all_plpl(param, npl, pl%rh, pl%vb, pl%renc, dt, nenc, index1, index2, lvdotr)
         else
            call encounter_check_all_plplm(param, nplm, nplt, pl%rh(:,1:nplm), pl%vb(:,1:nplm), pl%rh(:,nplm+1:npl), &
                  pl%vb(:,nplm+1:npl), pl%renc(1:nplm), pl%renc(nplm+1:npl), dt, nenc, index1, index2, lvdotr)
         end if
         
         lany_encounter = nenc > 0_I8B
         if (lany_encounter) then
            call plplenc_list%resize(nenc)
            call move_alloc(lvdotr, plplenc_list%lvdotr)
            call move_alloc(index1, plplenc_list%index1)
            call move_alloc(index2, plplenc_list%index2)
         end if

         if (lany_encounter) then
            do k = 1_I8B, nenc
               i = plplenc_list%index1(k)
               j = plplenc_list%index2(k)
               plplenc_list%id1(k) = pl%id(i)
               plplenc_list%id2(k) = pl%id(j)
               plplenc_list%status(k) = ACTIVE
               plplenc_list%level(k) = irec
               plplenc_list%x1(:,k) = pl%rh(:,i)
               plplenc_list%x2(:,k) = pl%rh(:,j)
               plplenc_list%v1(:,k) = pl%vb(:,i) - cb%vb(:)
               plplenc_list%v2(:,k) = pl%vb(:,j) - cb%vb(:)
               pl%lencounter(i) = .true.
               pl%lencounter(j) = .true.
               pl%levelg(i) = irec
               pl%levelm(i) = irec
               pl%levelg(j) = irec
               pl%levelm(j) = irec
               pl%nplenc(i) = pl%nplenc(i) + 1
               pl%nplenc(j) = pl%nplenc(j) + 1
            end do
         end if

      end associate

      return
   end function symba_encounter_check_pl


   module function symba_encounter_check(self, param, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between test particles and massive bodies in the plplenc and pltpenc list.
      !! Note: This method works for the polymorphic symba_pltpenc and symba_plplenc types.
      !!
      !! Adapted from portions of David E. Kaufmann's Swifter routine: symba_step_recur.f90
      implicit none
      ! Arguments
      class(symba_encounter),     intent(inout) :: self           !! SyMBA pl-pl encounter list object
      class(swiftest_parameters), intent(inout) :: param          !! Current swiftest run configuration parameters
      class(symba_nbody_system),  intent(inout) :: system         !! SyMBA nbody system object
      real(DP),                   intent(in)    :: dt             !! step size
      integer(I4B),               intent(in)    :: irec           !! Current recursion level 
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter  
      ! Internals
      integer(I4B)              :: i, j, k, lidx, nenc_enc
      real(DP), dimension(NDIM) :: xr, vr
      logical                   :: isplpl
      real(DP)                  :: rlim2, rji2, rcrit12
      logical, dimension(:), allocatable :: lencmask, lencounter
      integer(I4B), dimension(:), allocatable :: eidx

      lany_encounter = .false.
      if (self%nenc == 0) return

      select type(self)
      class is (symba_plplenc)
         isplpl = .true.
      class is (symba_pltpenc)
         isplpl = .false.
      end select

      select type(pl => system%pl)
      class is (symba_pl)
         select type(tp => system%tp)
         class is (symba_tp)
            allocate(lencmask(self%nenc))
            lencmask(:) = (self%status(1:self%nenc) == ACTIVE) .and. (self%level(1:self%nenc) == irec - 1)
            nenc_enc = count(lencmask(:))
            if (nenc_enc == 0) return

            call pl%set_renc(irec)

            allocate(eidx(nenc_enc))
            allocate(lencounter(nenc_enc))
            eidx(:) = pack([(k, k = 1, self%nenc)], lencmask(:))
            lencounter(:) = .false.
            if (isplpl) then
               do concurrent(lidx = 1:nenc_enc)
                  k = eidx(lidx)
                  i = self%index1(k)
                  j = self%index2(k)
                  xr(:) = pl%rh(:,j) - pl%rh(:,i)
                  vr(:) = pl%vb(:,j) - pl%vb(:,i)
                  rcrit12 = pl%renc(i) + pl%renc(j)
                  call encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), rcrit12, dt, lencounter(lidx), self%lvdotr(k))
                  if (lencounter(lidx)) then
                     rlim2 = (pl%radius(i) + pl%radius(j))**2
                     rji2 = dot_product(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                     lencounter(lidx) = rji2 > rlim2
                  end if
               end do
            else
               do concurrent(lidx = 1:nenc_enc)
                  k = eidx(lidx)
                  i = self%index1(k)
                  j = self%index2(k)
                  xr(:) = tp%rh(:,j) - pl%rh(:,i)
                  vr(:) = tp%vb(:,j) - pl%vb(:,i)
                  call encounter_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%renc(i), dt, &
                                           lencounter(lidx), self%lvdotr(k))
                  if (lencounter(lidx)) then
                     rlim2 = (pl%radius(i))**2
                     rji2 = dot_product(xr(:), xr(:))! Check to see if these are physically overlapping bodies first, which we should ignore
                     lencounter(lidx) = rji2 > rlim2
                  end if
               end do
            end if
            lany_encounter = any(lencounter(:))
            if (lany_encounter) then
               nenc_enc = count(lencounter(:))
               eidx(1:nenc_enc) = pack(eidx(:), lencounter(:))
               do lidx = 1, nenc_enc
                  k = eidx(lidx)
                  i = self%index1(k)
                  j = self%index2(k)
                  pl%levelg(i) = irec
                  pl%levelm(i) = MAX(irec, pl%levelm(i))
                  if (isplpl) then
                     pl%levelg(j) = irec
                     pl%levelm(j) = MAX(irec, pl%levelm(j))
                  else
                     tp%levelg(j) = irec
                     tp%levelm(j) = MAX(irec, tp%levelm(j))
                  end if
                  self%level(k) = irec
               end do
            end if   
         end select
      end select

      return
   end function symba_encounter_check


   module function symba_encounter_check_tp(self, param, system, dt, irec) result(lany_encounter)
      !! author: David A. Minton
      !!
      !! Check for an encounter between test particles and massive bodies.
      !!
      implicit none
      ! Arguments
      class(symba_tp),            intent(inout) :: self   !! SyMBA test particle object  
      class(swiftest_parameters), intent(inout) :: param  !! Current swiftest run configuration parameters
      class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
      real(DP),                   intent(in)    :: dt     !! step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
      ! Result
      logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      ! Internals
      integer(I4B)                              :: plind, tpind
      integer(I8B)                              :: k, nenc
      logical,      dimension(:),   allocatable :: lvdotr
      integer(I4B), dimension(:),   allocatable :: index1, index2
 
      lany_encounter = .false.
      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, pl => system%pl, npl => system%pl%nbody)
         call pl%set_renc(irec)
         call encounter_check_all_pltp(param, npl, ntp, pl%rh, pl%vb, tp%rh, tp%vb, pl%renc, dt, nenc, index1, index2, lvdotr) 
   
         lany_encounter = nenc > 0
         if (lany_encounter) then 
            associate(pltpenc_list => system%pltpenc_list)
               call pltpenc_list%resize(nenc)
               pltpenc_list%status(1:nenc) = ACTIVE
               pltpenc_list%level(1:nenc) = irec
               call move_alloc(index1, pltpenc_list%index1)
               call move_alloc(index2, pltpenc_list%index2)
               call move_alloc(lvdotr, pltpenc_list%lvdotr)
               pltpenc_list%id1(1:nenc) = pl%id(pltpenc_list%index1(1:nenc))
               pltpenc_list%id2(1:nenc) = tp%id(pltpenc_list%index2(1:nenc))
               select type(pl)
               class is (symba_pl)
                  pl%lencounter(1:npl) = .false.
                  do k = 1_I8B, nenc
                     plind = pltpenc_list%index1(k)
                     tpind = pltpenc_list%index2(k)
                     pl%lencounter(plind) = .true.
                     pl%levelg(plind) = irec
                     pl%levelm(plind) = irec
                     tp%levelg(tpind) = irec
                     tp%levelm(tpind) = irec
                     pl%ntpenc(plind) = pl%ntpenc(plind) + 1
                     tp%nplenc(tpind) = tp%nplenc(tpind) + 1
                  end do
               end select
            end associate
         end if
      end associate

      return
   end function symba_encounter_check_tp

end submodule s_symba_encounter_check
