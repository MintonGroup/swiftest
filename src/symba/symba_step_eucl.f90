submodule (symba) s_symba_step_eucl
contains
   module procedure symba_step_eucl
   !! author: Jacob R. Elliott
   !!
   !! Step planets and active test particles ahead in democratic heliocentric coordinates, descending the recursive 
   !!   branch if necessary to handle possible close encounters. 
   !! Uses the single-loop blocking to evaluate the Euclidean distance matri
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_step.f90
   !! Adapted from Hal Levison's Swift routine symba5_step_pl.f
use swiftest
implicit none
   logical(lgt)          :: lencounter, lvdotr
   integer(I4B)          :: i, j, irec, nplm, k, counter
   integer(I4B), allocatable :: plpl_encounters_indices(:), pltp_encounters_indices(:)
   real(DP), dimension(ndim) :: xr, vr
   
   integer(I4B), allocatable, dimension(:) :: pltp_encounters, pltp_lvdotr
   integer(I4B), allocatable, dimension(:) :: plpl_encounters, plpl_lvdotr
   
! executable code

   ! initialize massive bodies
   symba_plA%nplenc(1:npl) = 0 ! number of massive body encounters this particular massive body has
   symba_plA%ntpenc(1:npl) = 0 ! number of test particle encounters this particle massive body has
   symba_plA%levelg(1:npl) = -1 ! 
   symba_plA%levelm(1:npl) = -1 ! 
   symba_plA%index_parent(1:npl) = (/ (i, i=1,npl)/)
   symba_plA%index_child(:,1:npl) = 0

   ! initialize test particles
   symba_tpA%nplenc(1:ntp) = 0 
   symba_tpA%levelg(1:ntp) = -1
   symba_tpA%levelm(1:ntp) = -1

   nplplenc = 0 ! number of encounters in the entire run 
   npltpenc = 0

! all this needs to be changed to the tree search function for encounters
   allocate(plpl_encounters(num_plpl_comparisons))
   allocate(plpl_lvdotr(num_plpl_comparisons))
   plpl_encounters = 0
   plpl_lvdotr = 0

   ! call util_dist_eucl_plpl(npl,symba_plA%xh, num_plpl_comparisons, k_plpl, dist_plpl_array) 
   ! call util_dist_eucl_plpl(npl,symba_plA%vh, num_plpl_comparisons, k_plpl, vel_plpl_array) 
   call symba_chk_eucl(num_plpl_comparisons, k_plpl, symba_plA, dt, plpl_encounters, plpl_lvdotr, nplplenc)

   ! here i'll order the encounters
   ! nplplenc = count(plpl_encounters > 0)
   ! print *,'step nplplenc: ',nplplenc
   if(nplplenc>0)then

      allocate(plpl_encounters_indices(nplplenc))

      ! plpl_encounters_indices = pack(plpl_encounters,plpl_encounters > 0)
      ! so it turns out this is significantly faster than the pack command
      counter = 1
      do k = 1,num_plpl_comparisons
         if(plpl_encounters(k).gt.0)then
            plpl_encounters_indices(counter) = k
            counter = counter + 1
         endif
      enddo

      symba_plA%lmerged(k_plpl(1,plpl_encounters_indices)) = .false. ! they have not merged yet
      symba_plA%nplenc(k_plpl(1,plpl_encounters_indices)) = symba_plA%nplenc(k_plpl(1,plpl_encounters_indices)) + 1 ! number of particles that massive body "i" has close encountered
      symba_plA%levelg(k_plpl(1,plpl_encounters_indices)) = 0 ! recursion level
      symba_plA%levelm(k_plpl(1,plpl_encounters_indices)) = 0 ! recursion level
      symba_plA%nchild(k_plpl(1,plpl_encounters_indices)) = 0 
      ! for the j particle
      symba_plA%lmerged(k_plpl(2,plpl_encounters_indices)) = .false.
      symba_plA%nplenc(k_plpl(2,plpl_encounters_indices)) = symba_plA%nplenc(k_plpl(2,plpl_encounters_indices)) + 1
      symba_plA%levelg(k_plpl(2,plpl_encounters_indices)) = 0
      symba_plA%levelm(k_plpl(2,plpl_encounters_indices)) = 0
      symba_plA%nchild(k_plpl(2,plpl_encounters_indices)) = 0

      plplenc_list%status(1:nplplenc) = active ! you are in an encounter
      plplenc_list%lvdotr(1:nplplenc) = plpl_lvdotr(plpl_encounters_indices)! flag of relative accelerations to say if there will be a close encounter in next timestep 
      plplenc_list%level(1:nplplenc)  = 0 ! recursion level
      plplenc_list%index1(1:nplplenc) = k_plpl(1,plpl_encounters_indices) ! index of first massive body in encounter
      plplenc_list%index2(1:nplplenc) = k_plpl(2,plpl_encounters_indices) ! index of second massive body in encounter
      deallocate(plpl_encounters_indices)
   endif
   
   deallocate(plpl_encounters, plpl_lvdotr)

   if(ntp>0)then
       allocate(pltp_encounters(num_pltp_comparisons))
       allocate(pltp_lvdotr(num_pltp_comparisons))


       pltp_encounters = 0
       pltp_lvdotr = 0

      ! call util_dist_eucl_pltp(npl, ntp, symba_plA%xh, symba_tpA%xh, &
      !    num_pltp_comparisons, k_pltp, dist_pltp_array)
      ! call util_dist_eucl_pltp(npl, ntp, symba_plA%vh, symba_tpA%vh, &
      !    num_pltp_comparisons, k_pltp, vel_pltp_array)
      call symba_chk_eucl_pltp(num_pltp_comparisons, k_pltp, symba_plA, symba_tpA, dt, pltp_encounters, pltp_lvdotr, npltpenc)
   
      ! npltpenc = count(pltp_encounters > 0)
      ! print *,'step npltpenc: ',npltpenc
      if(npltpenc>0)then

         allocate(pltp_encounters_indices(npltpenc))

         counter = 1
         do k = 1,num_pltp_comparisons
            if(pltp_encounters(k).gt.0)then
               pltp_encounters_indices(counter) = k
               counter = counter + 1
            endif
         enddo

         symba_plA%ntpenc(k_pltp(1,pltp_encounters_indices)) = symba_plA%ntpenc(k_pltp(1,pltp_encounters_indices)) + 1
         symba_plA%levelg(k_pltp(1,pltp_encounters_indices)) = 0
         symba_plA%levelm(k_pltp(1,pltp_encounters_indices)) = 0

         symba_tpA%nplenc(k_pltp(2,pltp_encounters_indices)) = symba_tpA%nplenc(k_pltp(2,pltp_encounters_indices)) + 1
         symba_tpA%levelg(k_pltp(2,pltp_encounters_indices)) = 0
         symba_tpA%levelm(k_pltp(2,pltp_encounters_indices)) = 0

         pltpenc_list%status(1:npltpenc) = active
         pltpenc_list%lvdotr(1:npltpenc) = pltp_lvdotr(pltp_encounters_indices)
         pltpenc_list%level(1:npltpenc) = 0
         pltpenc_list%indexpl(1:npltpenc) = k_pltp(1,pltp_encounters_indices)
         pltpenc_list%indextp(1:npltpenc) = k_pltp(1,pltp_encounters_indices)

         deallocate(pltp_encounters_indices)
      endif

      deallocate(pltp_encounters, pltp_lvdotr)
   endif
   
! end of things that need to be changed in the tree

   nplm = count(symba_plA%mass > mtiny)
   ! flag to see if there was an encounter
   lencounter = ((nplplenc > 0) .or. (npltpenc > 0))

   if (lencounter) then ! if there was an encounter, we need to enter symba_step_interp to see if we need recursion
      call symba_step_interp_eucl(lextra_force, lclose, t, npl, nplm, config%nplmax, ntp, config%ntpmax, symba_plA, symba_tpA, config%j2rp2, config%j4rp4,&
         dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,&
         mergesub_list, encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)
      lfirst = .true.
   else ! otherwise we can just advance the particles
      call symba_step_helio(lfirst, lextra_force, t, npl, nplm, config%nplmax, ntp, config%ntpmax, symba_plA, symba_tpA, &
         config%j2rp2, config%j4rp4, dt)
   end if

   return

   end procedure symba_step_eucl
end submodule s_symba_step_eucl
