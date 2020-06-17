submodule (symba) s_symba_discard_merge_pl
contains
   module procedure symba_discard_merge_pl
   !! author: David A. Minton
   !!
   !! Merge planets
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_discard_merge_pl.f90
   !! Adapted from Hal Levison's Swift routine discard_mass_merge.f
use swiftest
implicit none
   integer(I4B)            :: i, j, nchild, indexchild, enc_big, index1, index2, indexk 
   real(DP)              :: m, mmax, mtot, r, r3, mu, energy, ap, v2, msun
   real(DP), dimension(ndim)   :: x, v, vbs
   integer(I4B), dimension(npl)  :: array_child

! executable code
   msun = symba_pla%helio%swiftest%mass(1)
   vbs(:) = symba_pla%helio%swiftest%vb(:,1)
   do i = 1, nplplenc
      if (plplenc_list%status(i) == merged) then
         index1 = plplenc_list%index1(i)
         index2 = plplenc_list%index2(i)
         ! this if statement is for if lfragmentation = false
         if ((symba_pla%helio%swiftest%status(index1) == active) .and.                                &
             (symba_pla%helio%swiftest%status(index2) == active)) then

            enc_big = plplenc_list%index1(i)

            m = symba_pla%helio%swiftest%mass(enc_big)
            r = symba_pla%helio%swiftest%radius(enc_big)
            r3 = r**3
            mmax = m
            mtot = m
            x(:) = m*symba_pla%helio%swiftest%xh(:,enc_big)
            v(:) = m*symba_pla%helio%swiftest%vb(:,enc_big)
            indexk = enc_big

            nchild = symba_pla%nchild(enc_big)
            array_child(1:npl) = symba_pla%index_child(1:npl,enc_big)

            do j = 1, nchild
               indexchild = array_child(j)
               m = symba_pla%helio%swiftest%mass(indexchild)
               r = symba_pla%helio%swiftest%radius(indexchild)
               r3 = r3 + r**3
               mtot = mtot + m
               x(:) = x(:) + m*symba_pla%helio%swiftest%xh(:,indexchild)
               v(:) = v(:) + m*symba_pla%helio%swiftest%vb(:,indexchild)
               if (m > mmax) then
                  mmax = m
                  indexk = indexchild
               end if
            end do
            x(:) = x(:)/mtot
            v(:) = v(:)/mtot
            r = r3**(1.0_DP/3.0_DP)
            symba_pla%helio%swiftest%mass(indexk) = mtot
            symba_pla%helio%swiftest%radius(indexk) = r
            symba_pla%helio%swiftest%xh(:,indexk) = x(:)
            symba_pla%helio%swiftest%vb(:,indexk) = v(:)
            symba_pla%helio%swiftest%vh(:,indexk) = v(:) - vbs(:)
            mu = msun*mtot/(msun + mtot)
            r = sqrt(dot_product(x(:), x(:)))
            v(:) = symba_pla%helio%swiftest%vh(:,indexk)
            v2 = dot_product(v(:), v(:))
            energy = -1.0_DP*msun*mtot/r + 0.5_DP*mu*v2
            ap = -1.0_DP*msun*mtot/(2.0_DP*energy)
            symba_pla%helio%swiftest%rhill(indexk) = ap*(((mu/msun)/3.0_DP)**(1.0_DP/3.0_DP))
            array_child(1:npl) = symba_pla%index_child(1:npl,enc_big)
            indexchild = enc_big
            ldiscard = .true.
            do j = 0, nchild
               if (indexchild /= indexk) then
                  symba_pla%helio%swiftest%status(indexchild) = merged
               end if
               indexchild = array_child(j+1)
            end do

         else if ((symba_pla%helio%swiftest%status(index1) == disruption) .and.    &                              
             (symba_pla%helio%swiftest%status(index2) == disruption)) then 

            enc_big = plplenc_list%index1(i)
            nchild = symba_pla%nchild(enc_big)
            array_child(1:npl) = symba_pla%index_child(1:npl,enc_big)
            do j = 1, nchild
               symba_pla%helio%swiftest%status(array_child(j)) = inactive
            end do
            ldiscard = .true.
         else if ((symba_pla%helio%swiftest%status(index1) == supercatastrophic) .and.   &                               
             (symba_pla%helio%swiftest%status(index2) == supercatastrophic)) then 
            
            enc_big = plplenc_list%index1(i)
            nchild = symba_pla%nchild(enc_big)
            array_child(1:npl) = symba_pla%index_child(1:npl,enc_big)
            do j = 1, nchild
               symba_pla%helio%swiftest%status(array_child(j)) = inactive
            end do
            ldiscard = .true.
         else if ((symba_pla%helio%swiftest%status(index1) == hit_and_run) .and.    &                            
             (symba_pla%helio%swiftest%status(index2) == hit_and_run)) then 

            enc_big = plplenc_list%index1(i)
            nchild = symba_pla%nchild(enc_big)
            array_child(1:npl) = symba_pla%index_child(1:npl,enc_big)
            do j = 1, nchild
               symba_pla%helio%swiftest%status(array_child(j)) = inactive
            end do
            ldiscard = .true.
         else if ((symba_pla%helio%swiftest%status(index1) == graze_and_merge) .and.  &                              
             (symba_pla%helio%swiftest%status(index2) == graze_and_merge)) then 

            enc_big = plplenc_list%index1(i)
            nchild = symba_pla%nchild(enc_big)
            array_child(1:npl) = symba_pla%index_child(1:npl,enc_big)
            do j = 1, nchild
               symba_pla%helio%swiftest%status(array_child(j)) = inactive
            end do
            ldiscard = .true.
         end if
      end if
   end do

   return
   
   end procedure symba_discard_merge_pl
end submodule s_symba_discard_merge_pl
