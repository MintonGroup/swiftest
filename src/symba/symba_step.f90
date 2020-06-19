submodule (symba) s_symba_step
contains
   module procedure symba_step
   !! author: David A. Minton
   !!
   !! Step planets and active test particles ahead in democratic heliocentric coordinates, descending the recursive
   !!   branch if necessary to handle possible close encounters
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_step.f90
   !! Adapted from Hal Levison's Swift routine symba5_step_pl.f
use swiftest
implicit none
   logical           :: lencounter, lvdotr
   integer(I4B)          :: i, j, irec, nplm
   real(DP), dimension(NDIM) :: xr, vr
   logical, save         :: lfirst = .true.
   
! executable code

    do i = 1,npl
      symba_plA%nplenc(i) = 0
      symba_plA%ntpenc(i) = 0
      symba_plA%levelg(i) = -1
      symba_plA%levelm(i) = -1
      symba_plA%index_parent(i) = i
      symba_plA%index_child(:,i) = 0
   end do
   do i =1,ntp
      symba_tpA%nplenc(i) = 0 
      symba_tpA%levelg(i) = -1
      symba_tpA%levelm(i) = -1
   end do 


   !there should be some parallel bits in here

   nplplenc = 0
   npltpenc = 0
   if (symba_plA%mass(1) < config%mtiny) then
      nplm = 0
   else
      nplm = 1
   end if
   irec = 0

! all this needs to be changed to the tree search function for encounters

   do i = 2, npl
      if (symba_plA%mass(i) < config%mtiny) exit
      nplm = nplm + 1
      do j = i + 1, npl
         xr(:) = symba_plA%xh(:,j) - symba_plA%xh(:,i)
         vr(:) = symba_plA%vh(:,j) - symba_plA%vh(:,i)
         call symba_chk(xr(:), vr(:), symba_plA%rhill(i), &
            symba_plA%rhill(j), dt, irec, lencounter, lvdotr)
         if (lencounter) then
            nplplenc = nplplenc + 1
            if (nplplenc > nenmax) then
               write(*, *) "Swiftest Error:"
               write(*, *) "   pl-pl encounter list is full."
               write(*, *) "   stopping..."
               call util_exit(failure)
            end if
            plplenc_list%status(nplplenc) = ACTIVE
            plplenc_list%lvdotr(nplplenc) = lvdotr
            plplenc_list%level(nplplenc) = irec
            plplenc_list%index1(nplplenc) = i
            plplenc_list%index2(nplplenc) = j
            symba_plA%lmerged(i) = .false.
            symba_plA%nplenc(i) = symba_plA%nplenc(i) + 1
            symba_plA%levelg(i) = irec
            symba_plA%levelm(i) = irec
            symba_plA%nchild(i) = 0 
            symba_plA%lmerged(j) = .false.
            symba_plA%nplenc(j) = symba_plA%nplenc(j) + 1
            symba_plA%levelg(j) = irec
            symba_plA%levelm(j) = irec
            symba_plA%nchild(j) = 0
         end if
      end do
      do j = 1, ntp
         xr(:) = symba_tpA%xh(:,j) - symba_plA%xh(:,i)
         vr(:) = symba_tpA%vh(:,j) - symba_plA%vh(:,i)
         call symba_chk(xr(:), vr(:), symba_plA%rhill(i), 0.0_DP, dt, irec, lencounter, lvdotr)
         if (lencounter) then
            npltpenc = npltpenc + 1
            symba_plA%ntpenc(i) = symba_plA%ntpenc(i) + 1
            symba_plA%levelg(i) = irec
            symba_plA%levelm(i) = irec
            symba_tpA%nplenc(j) = symba_tpA%nplenc(j) + 1
            symba_tpA%levelg(j) = irec
            symba_tpA%levelm(j) = irec
            if (npltpenc > nenmax) then
               write(*, *) "Swiftest Error:"
               write(*, *) "   pl-tp encounter list is full."
               write(*, *) "   stopping..."
               call util_exit(failure)
            end if
            pltpenc_list%status(npltpenc) = ACTIVE
            pltpenc_list%lvdotr(npltpenc) = lvdotr
            pltpenc_list%level(npltpenc) = irec
            pltpenc_list%indexpl(npltpenc) = i
            pltpenc_list%indextp(npltpenc) = j
         end if
      end do
   end do

! end of things that need to be changed in the tree

   lencounter = ((nplplenc > 0) .or. (npltpenc > 0))
   if (lencounter) then
      call symba_step_interp(config%lextra_force, config%lclose, t, npl, nplm, config%nplmax, &
         ntp, config%ntpmax, symba_plA, symba_tpA, config%j2rp2, config%j4rp4,   &
         dt, eoffset, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, &
         nmergesub, mergeadd_list, mergesub_list, config%encounter_file, config%out_type, &
         fragmax, param)
      lfirst = .true.
   else
      call symba_step_helio(lfirst, config%lextra_force, t, npl, nplm, config%nplmax, ntp,&
         config%ntpmax, symba_plA, symba_tpA, &
         config%j2rp2, config%j4rp4, dt)
   end if

   return

   end procedure symba_step
end submodule s_symba_step
