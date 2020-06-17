submodule (symba) s_symba_getacch
contains
   module procedure symba_getacch
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of planets
   !!      Accelerations in an encounter are not included here
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_getacch.f90
   !! Adapted from Hal Levison's Swift routine symba5_getacch.f
use swiftest
implicit none
   integer(I4B)                     :: i, j, index_i, index_j
   real(DP)                       :: rji2, irij3, faci, facj, r2
   real(DP), dimension(ndim)            :: dx
   real(DP), dimension(npl)             :: irh
   real(DP), dimension(ndim, npl)         :: aobl

! executable code

   do i = 2, npl
      symba_pla%helio%ah(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do

   do i = 2, nplm
      do j = i + 1, npl
         if ((.not. symba_pla%lmerged(i)) .or. (.not. symba_pla%lmerged(j)) .or. &
             (symba_pla%index_parent(i) /= symba_pla%index_parent(j))) then
            dx(:) = symba_pla%helio%swiftest%xh(:,j) - symba_pla%helio%swiftest%xh(:,i)
            rji2 = dot_product(dx(:), dx(:))
            irij3 = 1.0_DP/(rji2*sqrt(rji2))
            if (irij3  .ne. irij3 ) then 
               write(*,*) "dx==0 for pl: ", i, "name:", symba_pla%helio%swiftest%name(i), &
            "and pl:", j, "name:", symba_pla%helio%swiftest%name(j)
            write(*,*) "dx==0 for pl: ", i, "xh:", symba_pla%helio%swiftest%xh(1,i), &
            "and pl:", j, "xh:", symba_pla%helio%swiftest%xh(1,j)
               write(*,*) "parent pl 1:", symba_pla%helio%swiftest%name(symba_pla%index_parent(i))
               write(*,*) "parent pl 2:", symba_pla%helio%swiftest%name(symba_pla%index_parent(j))
               stop
            end if
            faci = symba_pla%helio%swiftest%mass(i)*irij3
            facj = symba_pla%helio%swiftest%mass(j)*irij3
            symba_pla%helio%ah(:,i) = symba_pla%helio%ah(:,i) + facj*dx(:)
            symba_pla%helio%ah(:,j) = symba_pla%helio%ah(:,j) - faci*dx(:)
         end if
      end do
   end do

   do i = 1, nplplenc
      index_i = plplenc_list%index1(i)
      index_j = plplenc_list%index2(i)
      if ((.not. symba_pla%lmerged(index_i)) .or. (.not. symba_pla%lmerged(index_j))  &
          .or. (symba_pla%index_parent(index_i) /= symba_pla%index_parent(index_j))) then !need to update parent/children
         dx(:) = symba_pla%helio%swiftest%xh(:,index_j) - symba_pla%helio%swiftest%xh(:,index_i)
         rji2 = dot_product(dx(:), dx(:))
         irij3 = 1.0_DP/(rji2*sqrt(rji2))
         if (irij3  .ne. irij3 ) then 
               write(*,*) "dx==0 for pl: ", i, "name:", symba_pla%helio%swiftest%name(i), &
            "and pl:", j, "name:", symba_pla%helio%swiftest%name(j)
            write(*,*) "dx==0 for pl: ", i, "xh:", symba_pla%helio%swiftest%xh(1,i), &
            "and pl:", j, "xh:", symba_pla%helio%swiftest%xh(1,j)
               write(*,*) "parent pl 1:", symba_pla%helio%swiftest%name(symba_pla%index_parent(i))
               write(*,*) "parent pl 2:", symba_pla%helio%swiftest%name(symba_pla%index_parent(j))
               stop
         end if
         faci = symba_pla%helio%swiftest%mass(index_i)*irij3
         facj = symba_pla%helio%swiftest%mass(index_j)*irij3
         symba_pla%helio%ah(:,index_i) = symba_pla%helio%ah(:,index_i) - facj*dx(:)
         symba_pla%helio%ah(:,index_j) = symba_pla%helio%ah(:,index_j) + faci*dx(:)
      end if
   end do
   if (j2rp2 /= 0.0_DP) then
      !if (lmalloc) then
          !allocate(xh(ndim, npl),aobl(ndim, npl), irh(npl))
         !lmalloc = .false.
      !end if
      do i = 2, npl
         r2 = dot_product(symba_pla%helio%swiftest%xh(:,i), symba_pla%helio%swiftest%xh(:,i))
         irh(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc(symba_pla%helio%swiftest, j2rp2, j4rp4, symba_pla%helio%swiftest%xh(:,:), irh, aobl)
      do i = 2, npl
         symba_pla%helio%ah(:,i) = symba_pla%helio%ah(:,i) + aobl(:, i) - aobl(:, 1)
      end do
   end if
   if (lextra_force) call symba_user_getacch(t, npl, symba_pla)

   return

   end procedure symba_getacch
end submodule s_symba_getacch
