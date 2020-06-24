submodule (symba) s_symba_getacch
contains
   module procedure symba_getacch
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of planets
   !!      Accelerations in an encounter are not included here
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: symba_getacch.f90
   !! Adapted from Hal Levison's Swift routine symba5_getacch.f
use swiftest
implicit none
   integer(I4B)                     :: i, j, index_i, index_j
   real(DP)                       :: rji2, irij3, faci, facj, r2
   real(DP), dimension(NDIM)            :: dx
   real(DP), dimension(npl)             :: irh
   real(DP), dimension(NDIM, npl)         :: aobl

! executable code

   do i = 2, npl
      symba_plA%ah(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do

   do i = 2, nplm
      do j = i + 1, npl
         if ((.not. symba_plA%lmerged(i)) .or. (.not. symba_plA%lmerged(j)) .or. &
             (symba_plA%index_parent(i) /= symba_plA%index_parent(j))) then
            dx(:) = symba_plA%xh(:,j) - symba_plA%xh(:,i)
            rji2 = dot_product(dx(:), dx(:))
            irij3 = 1.0_DP/(rji2*sqrt(rji2))
            if (irij3  .ne. irij3 ) then 
               write(*,*) "dx==0 for pl: ", i, "name:", symba_plA%name(i), &
            "and pl:", j, "name:", symba_plA%name(j)
            write(*,*) "dx==0 for pl: ", i, "xh:", symba_plA%xh(1,i), &
            "and pl:", j, "xh:", symba_plA%xh(1,j)
               write(*,*) "parent pl 1:", symba_plA%name(symba_plA%index_parent(i))
               write(*,*) "parent pl 2:", symba_plA%name(symba_plA%index_parent(j))
               stop
            end if
            faci = symba_plA%mass(i)*irij3
            facj = symba_plA%mass(j)*irij3
            symba_plA%ah(:,i) = symba_plA%ah(:,i) + facj*dx(:)
            symba_plA%ah(:,j) = symba_plA%ah(:,j) - faci*dx(:)
         end if
      end do
   end do

   do i = 1, nplplenc
      index_i = plplenc_list%index1(i)
      index_j = plplenc_list%index2(i)
      if ((.not. symba_plA%lmerged(index_i)) .or. (.not. symba_plA%lmerged(index_j))  &
          .or. (symba_plA%index_parent(index_i) /= symba_plA%index_parent(index_j))) then !need to update parent/children
         dx(:) = symba_plA%xh(:,index_j) - symba_plA%xh(:,index_i)
         rji2 = dot_product(dx(:), dx(:))
         irij3 = 1.0_DP/(rji2*sqrt(rji2))
         if (irij3  .ne. irij3 ) then 
               write(*,*) "dx==0 for pl: ", i, "name:", symba_plA%name(i), &
            "and pl:", j, "name:", symba_plA%name(j)
            write(*,*) "dx==0 for pl: ", i, "xh:", symba_plA%xh(1,i), &
            "and pl:", j, "xh:", symba_plA%xh(1,j)
               write(*,*) "parent pl 1:", symba_plA%name(symba_plA%index_parent(i))
               write(*,*) "parent pl 2:", symba_plA%name(symba_plA%index_parent(j))
               stop
         end if
         faci = symba_plA%mass(index_i)*irij3
         facj = symba_plA%mass(index_j)*irij3
         symba_plA%ah(:,index_i) = symba_plA%ah(:,index_i) - facj*dx(:)
         symba_plA%ah(:,index_j) = symba_plA%ah(:,index_j) + faci*dx(:)
      end if
   end do
   if (config%j2rp2 /= 0.0_DP) then
      !if (lmalloc) then
          !allocate(xh(NDIM, npl),aobl(NDIM, npl), irh(npl))
         !lmalloc = .false.
      !end if
      do i = 2, npl
         r2 = dot_product(symba_plA%xh(:,i), symba_plA%xh(:,i))
         irh(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc(symba_plA, config%j2rp2, config%j4rp4, symba_plA%xh(:,:), irh, aobl)
      do i = 2, npl
         symba_plA%ah(:,i) = symba_plA%ah(:,i) + aobl(:, i) - aobl(:, 1)
      end do
   end if
   if (lextra_force) call symba_user_getacch(t, npl, symba_plA)

   return

   end procedure symba_getacch
end submodule s_symba_getacch
