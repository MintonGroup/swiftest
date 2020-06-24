submodule (symba) s_symba_casemerge 
contains
   module procedure symba_casemerge 
   !! author: Jennifer L. L. Pouplin and Carlisle A. Wishard
   !!
   !! Merge planetss
   !! Adapted from David E. Kaufmann's Swifter routine: symba_merge_pl .f90
   !! Adapted from Hal Levison's Swift routine symba5_merge.f
use swiftest
implicit none
 
   integer(I4B)           :: i, j, k, stat1, stat2, index1, index2, indexchild
   integer(I4B)           :: index1_child, index2_child, index1_parent, index2_parent
   integer(I4B)           :: name1, name2
   real(DP)             :: mtot
   real(DP)             :: eold, enew, mass1, mass2
   real(DP), dimension(NDIM)    :: xr, xnew, vnew
   integer(I4B), dimension(npl) :: array_keep_child, array_rm_child

! executable code
         index1 = plplenc_list%index1(index_enc)
         index2 = plplenc_list%index2(index_enc)
         index1_parent = symba_plA%index_parent(index1)
         index2_parent = symba_plA%index_parent(index2)
         mtot = m1 + m2
         xnew(:) = (m1*x1(:) + m2*x2(:))/mtot
         vnew(:) = (m1*v1(:) + m2*v2(:))/mtot
         name1 = symba_plA%name(index1)
         name2 = symba_plA%name(index2)
         mass1 = symba_plA%mass(index1)
         mass2 = symba_plA%mass(index2)
         stat1 = symba_plA%status(index1)
         stat2 = symba_plA%status(index2)
         write(*, *) "merging particles ", name1, " and ", name2, " at time t = ",t
         nmergesub = nmergesub + 1
         mergesub_list%name(nmergesub) = name1
         mergesub_list%status(nmergesub) = MERGED
         mergesub_list%xh(:,nmergesub) = x1(:)
         mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
         mergesub_list%mass(nmergesub) = mass1
         mergesub_list%radius(nmergesub) = rad1
         nmergesub = nmergesub + 1
         mergesub_list%name(nmergesub) = name2
         mergesub_list%status(nmergesub) = MERGED
         mergesub_list%xh(:,nmergesub) = x2(:)
         mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
         mergesub_list%mass(nmergesub) = mass2
         mergesub_list%radius(nmergesub) = rad2
         nmergeadd = nmergeadd + 1
         if (m2 > m1) then
            mergeadd_list%name(nmergeadd) = name2
            mergeadd_list%status(nmergeadd) = stat2

         else
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

         do k = 1, nplplenc
            if (plplenc_list%status(k) == ACTIVE) then
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
                        plplenc_list%status(k) = MERGED
                     else if ((index1_child == plplenc_list%index2(k)) .and. &
                        (index2_child == plplenc_list%index1(k))) then
                        plplenc_list%status(k) = MERGED
                     end if
                  end do
               end do
            end if
         end do

         symba_plA%xh(:,index1_parent) = xnew(:)
         symba_plA%vb(:,index1_parent) = vnew(:)
         symba_plA%xh(:,index2_parent) = xnew(:) 
         symba_plA%vb(:,index2_parent) = vnew(:)

         ! the children of parent one are the children we are keeping
         array_keep_child(1:npl) = symba_plA%index_child(1:npl,index1_parent)
         ! go through the children of the kept parent and add those children to the array of kept children
         do i = 1, symba_plA%nchild(index1_parent)
            indexchild = array_keep_child(i)
            symba_plA%xh(:,indexchild) = xnew(:)
            symba_plA%vb(:,indexchild) = vnew(:)
         end do
         ! the removed parent is assigned as a new child to the list of children of the kept parent
         ! gives kept parent a new child 
         symba_plA%index_child((symba_plA%nchild(index1_parent)+1),index1_parent) = index2_parent
         array_rm_child(1:npl) = symba_plA%index_child(1:npl,index2_parent)
         ! the parent of the removed parent is assigned to be the kept parent 
         ! gives removed parent a new parent
         symba_plA%index_parent(index2) = index1_parent
         ! go through the children of the removed parent and add those children to the array of removed children 
         do i = 1, symba_plA%nchild(index2_parent)
            symba_plA%index_parent(array_rm_child(i)) = index1_parent
            indexchild = array_rm_child(i)
            symba_plA%xh(:,indexchild) = xnew(:)
            symba_plA%vb(:,indexchild) = vnew(:)
         end do
         ! go through the children of the removed parent and add those children to the list of children of the kept parent
         do i = 1, symba_plA%nchild(index2_parent)
            symba_plA%index_child(symba_plA%nchild(index1_parent)+i+1,index1_parent)= array_rm_child(i)
         end do 
         ! updates the number of children of the kept parent
         symba_plA%nchild(index1_parent) = symba_plA%nchild(index1_parent) + symba_plA%nchild(index2_parent) + 1

   return 
   end procedure symba_casemerge
end submodule s_symba_casemerge
