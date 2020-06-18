submodule (symba) s_symba_merge_pl
contains
   module procedure symba_merge_pl
   !! author: David A. Minton
   !!
   !! Check whether or not bodies are colliding or on collision path
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_merge_pl.f90
   !! Adapted from Hal Levison's Swift routine symba5_merge.f
   use swiftest
   implicit none
   logical(lgt)           :: lmerge
   integer(I4B)           :: i, j, k, stat1, stat2, index1, index2, index_keep, index_rm, indexchild
   integer(I4B)           :: index1_child, index2_child, index1_parent, index2_parent, index_big1, index_big2
   integer(I4B)           :: name1, name2
   real(DP)             :: r2, rlim, rlim2, vdotr, tcr2, dt2, mtot, a, e, q, m1, m2, mtmp, mmax 
   real(DP)             :: eold, enew, rad1, rad2, mass1, mass2
   real(DP), dimension(ndim)    :: xr, vr, x1, v1, x2, v2, xnew, vnew
   integer(I4B), dimension(npl) :: array_index1_child, array_index2_child, array_keep_child, array_rm_child

! executable code
   lmerge = .false.

   index1 = plplenc_list%index1(index_enc)
   index2 = plplenc_list%index2(index_enc)

   rlim = symba_plA%radius(index1) + symba_plA%radius(index2)
   xr(:) = symba_plA%xh(:,index2) - symba_plA%xh(:,index1)
   r2 = dot_product(xr(:), xr(:))
   rlim2 = rlim*rlim
   ! checks if bodies are actively colliding in this time step
   if (rlim2 >= r2) then 
      lmerge = .true.
   ! if they are not actively colliding in  this time step, 
   !checks if they are going to collide next time step based on velocities and q
   else 
      vr(:) = symba_plA%vb(:,index2) - symba_plA%vb(:,index1)
      vdotr = dot_product(xr(:), vr(:))
      if (plplenc_list%lvdotr(index_enc) .and. (vdotr > 0.0_DP)) then 
         tcr2 = r2/dot_product(vr(:), vr(:))
         dt2 = dt*dt
         if (tcr2 <= dt2) then
            mtot = symba_plA%mass(index1) + symba_plA%mass(index2)
            call orbel_xv2aeq(xr(:), vr(:), mtot, a, e, q)
            if (q < rlim) lmerge = .true.
         end if
         ! if no collision is going to happen, write as close encounter, not  merger
         if (.not. lmerge) then
            if (encounter_file /= "") then
               name1 = symba_plA%name(index1)
               m1 = symba_plA%mass(index1)
               rad1 = symba_plA%radius(index1)
               x1(:) = symba_plA%xh(:,index1)
               v1(:) = symba_plA%vb(:,index1) - vbs(:)
               name2 = symba_plA%name(index2)
               m2 = symba_plA%mass(index2)
               rad2 = symba_plA%radius(index2)
               x2(:) = symba_plA%xh(:,index2)
               v2(:) = symba_plA%vb(:,index2) - vbs(:)

               call io_write_encounter(t, name1, name2, m1, m2, rad1, rad2, x1(:), x2(:), &
                  v1(:), v2(:), encounter_file, out_type)
            end if
         end if
      end if
   end if
   !set up the merger for symba_discard_merge_pl 
   if (lmerge) then
      symba_plA%lmerged(index1) = .true.
      symba_plA%lmerged(index2) = .true.
      index1_parent = symba_plA%index_parent(index1)
      m1 = symba_plA%mass(index1_parent)
      mass1 = m1 
      rad1 = symba_plA%radius(index1_parent)
      x1(:) = m1*symba_plA%xh(:,index1_parent)
      v1(:) = m1*symba_plA%vb(:,index1_parent)
      mmax = m1
      name1 = symba_plA%name(index1_parent)
      index_big1 = index1_parent
      stat1 = symba_plA%status(index1_parent)
      array_index1_child(1:npl) = symba_plA%index_child(1:npl,index1_parent)
      do i = 1, symba_plA%nchild(index1_parent) ! initialize an array of children
         index1_child = array_index1_child(i)
         mtmp = symba_plA%mass(index1_child)
         if (mtmp > mmax) then
            mmax = mtmp
            name1 = symba_plA%name(index1_child)
            index_big1 = index1_child
            stat1 = symba_plA%status(index1_child)
         end if
         m1 = m1 + mtmp
         x1(:) = x1(:) + mtmp*symba_plA%xh(:,index1_child)
         v1(:) = v1(:) + mtmp*symba_plA%vb(:,index1_child)
      end do
      x1(:) = x1(:)/m1
      v1(:) = v1(:)/m1
      index2_parent = symba_plA%index_parent(index2)
      m2 = symba_plA%mass(index2_parent)
      mass2 = m2
      rad2 = symba_plA%radius(index2_parent)
      x2(:) = m2*symba_plA%xh(:,index2_parent)
      v2(:) = m2*symba_plA%vb(:,index2_parent)
      mmax = m2
      name2 = symba_plA%name(index2_parent)
      index_big2 = index2_parent
      stat2 = symba_plA%status(index2_parent)
      array_index2_child(1:npl) = symba_plA%index_child(1:npl,index2_parent)
      do i = 1, symba_plA%nchild(index2_parent)
         index2_child = array_index2_child(i)
         mtmp = symba_plA%mass(index2_child)
         if (mtmp > mmax) then
            mmax = mtmp
            name2 = symba_plA%name(index2_child)
            index_big2 = index2_child
            stat2 = symba_plA%status(index2_child)
         end if
         m2 = m2 + mtmp
         x2(:) = x2(:) + mtmp*symba_plA%xh(:,index2_child)
         v2(:) = v2(:) + mtmp*symba_plA%vb(:,index2_child)
      end do
      x2(:) = x2(:)/m2
      v2(:) = v2(:)/m2
      mtot = m1 + m2
      xnew(:) = (m1*x1(:) + m2*x2(:))/mtot
      vnew(:) = (m1*v1(:) + m2*v2(:))/mtot
      write(*, *) "merging particles ", name1, " and ", name2, " at time t = ",t
      nmergesub = nmergesub + 1
      mergesub_list%name(nmergesub) = name1
      mergesub_list%status(nmergesub) = merged
      mergesub_list%xh(:,nmergesub) = x1(:)
      mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
      mergesub_list%mass(nmergesub) = mass1
      mergesub_list%radius(nmergesub) = rad1
      nmergesub = nmergesub + 1
      mergesub_list%name(nmergesub) = name2
      mergesub_list%status(nmergesub) = merged
      mergesub_list%xh(:,nmergesub) = x2(:)
      mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
      mergesub_list%mass(nmergesub) = mass2
      mergesub_list%radius(nmergesub) = rad2
      nmergeadd = nmergeadd + 1
      if (m2 > m1) then
         index_keep = index_big2
         index_rm = index_big1
         mergeadd_list%name(nmergeadd) = name2
         mergeadd_list%status(nmergeadd) = stat2

      else
         index_keep = index_big1
         index_rm = index_big2
         mergeadd_list%name(nmergeadd) = name1
         mergeadd_list%status(nmergeadd) = stat1

      end if
      mergeadd_list%ncomp(nmergeadd) = 2
      mergeadd_list%xh(:,nmergeadd) = xnew(:)
      mergeadd_list%vh(:,nmergeadd) = vnew(:) - vbs(:)
      eold = 0.5_DP*(m1*dot_product(v1(:), v1(:)) + m2*dot_product(v2(:), v2(:)))
      xr(:) = x2(:) - x1(:)
      eold = eold - m1*m2/sqrt(dot_product(xr(:), xr(:)))
      enew = 0.5_DP*mtot*dot_product(vnew(:), vnew(:))
      eoffset = eoffset + eold - enew

      !write(*,*) "symba_merge_pl.f90 name", mergeadd_list%name(nmergeadd)
      !write(*,*) "symba_merge_pl.f90 xh", mergeadd_list%xh(:,nmergeadd)
      !write(*,*) "symba_merge_pl.f90 vh", mergeadd_list%vh(:,nmergeadd)
      !write(*,*) "symba_merge_pl.f90 eoffset", eoffset
      do k = 1, nplplenc                          !go through the encounter list and for particles actively encoutering, get their children
         if (plplenc_list%status(k) == active) then
            do i = 0, symba_plA%nchild(index1_parent)
               if (i == 0) then 
                  index1_child = index1_parent
               else
                  index1_child = array_index1_child(i)
               end if 
               do j = 0, symba_plA%nchild(index2_parent)
                  if (j == 0) then
                     index2_child = index2_parent
                  else
                     index2_child = array_index2_child(j)
                  end if
                  if ((index1_child == plplenc_list%index1(k)) .and. (index2_child == plplenc_list%index2(k))) then
                     plplenc_list%status(k) = merged
                  else if ((index1_child == plplenc_list%index2(k)) .and. (index2_child == plplenc_list%index1(k))) then
                     plplenc_list%status(k) = merged
                  end if
               end do
            end do
         end if
      end do
      symba_plA%xh(:,index1_parent) = xnew(:)
      symba_plA%vb(:,index1_parent) = vnew(:)
      symba_plA%xh(:,index2_parent) = xnew(:) 
      symba_plA%vb(:,index2_parent) = vnew(:) 
      array_keep_child(1:npl) = symba_plA%index_child(1:npl,index1_parent)
      do i = 1, symba_plA%nchild(index1_parent)
         indexchild = array_keep_child(i)
         symba_plA%xh(:,indexchild) = xnew(:)
         symba_plA%vb(:,indexchild) = vnew(:)
      end do

      symba_plA%index_child((symba_plA%nchild(index1_parent)+1),index1_parent) = index2_parent
      array_rm_child(1:npl) = symba_plA%index_child(1:npl,index2_parent)
      symba_plA%index_parent(index2) = index1_parent

      do i = 1, symba_plA%nchild(index2_parent)
         symba_plA%index_parent(array_rm_child(i)) = index1_parent
         indexchild = array_rm_child(i)
         symba_plA%xh(:,indexchild) = xnew(:)
         symba_plA%vb(:,indexchild) = vnew(:)
      end do
      do i = 1, symba_plA%nchild(index2_parent)
         symba_plA%index_child(symba_plA%nchild(index1_parent)+i+1,index1_parent)= array_rm_child(i)
      end do 
      symba_plA%nchild(index1_parent) = symba_plA%nchild(index1_parent) + symba_plA%nchild(index2_parent) + 1
   end if

   return

   end procedure symba_merge_pl
end submodule s_symba_merge_pl
