submodule (symba) s_symba_discard_merge_pl
contains
   module procedure symba_discard_merge_pl
   !! author: David A. Minton
   !!
   !! Merge planets
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_merge_pl.f90
   !! Adapted from Hal Levison's Swift routine discard_mass_merge.f
use swiftest
implicit none
   integer(I4B)            :: i, j, nchild, indexchild, enc_big, index1, index2, indexk 
   real(DP)              :: m, mmax, mtot, r, r3, mu, energy, ap, v2, msun
   real(DP), dimension(NDIM)   :: x, v, vbs
   integer(I4B), dimension(npl)  :: array_child

! executable code
   msun = symba_plA%mass(1)
   vbs(:) = symba_plA%vb(:,1)
   do i = 1, nplplenc
      if (plplenc_list%status(i) == MERGED) then
         index1 = plplenc_list%index1(i)
         index2 = plplenc_list%index2(i)
         ! this if statement is for if lfragmentation = false
         if ((symba_plA%status(index1) == ACTIVE) .and.                                &
             (symba_plA%status(index2) == ACTIVE)) then

            enc_big = plplenc_list%index1(i)

            m = symba_plA%mass(enc_big)
            r = symba_plA%radius(enc_big)
            r3 = r**3
            mmax = m
            mtot = m
            x(:) = m*symba_plA%xh(:,enc_big)
            v(:) = m*symba_plA%vb(:,enc_big)
            indexk = enc_big

            nchild = symba_plA%nchild(enc_big)
            array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)

            do j = 1, nchild
               indexchild = array_child(j)
               m = symba_plA%mass(indexchild)
               r = symba_plA%radius(indexchild)
               r3 = r3 + r**3
               mtot = mtot + m
               x(:) = x(:) + m*symba_plA%xh(:,indexchild)
               v(:) = v(:) + m*symba_plA%vb(:,indexchild)
               if (m > mmax) then
                  mmax = m
                  indexk = indexchild
               end if
            end do
            x(:) = x(:)/mtot
            v(:) = v(:)/mtot
            r = r3**(1.0_DP/3.0_DP)
            symba_plA%mass(indexk) = mtot
            symba_plA%radius(indexk) = r
            symba_plA%xh(:,indexk) = x(:)
            symba_plA%vb(:,indexk) = v(:)
            symba_plA%vh(:,indexk) = v(:) - vbs(:)
            mu = msun*mtot/(msun + mtot)
            r = sqrt(dot_product(x(:), x(:)))
            v(:) = symba_plA%vh(:,indexk)
            v2 = dot_product(v(:), v(:))
            energy = -1.0_DP*msun*mtot/r + 0.5_DP*mu*v2
            ap = -1.0_DP*msun*mtot/(2.0_DP*energy)
            symba_plA%rhill(indexk) = ap*(((mu/msun)/3.0_DP)**(1.0_DP/3.0_DP))
            array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
            indexchild = enc_big
            ldiscard = .true.
            do j = 0, nchild
               if (indexchild /= indexk) then
                  symba_plA%status(indexchild) = MERGED
               end if
               indexchild = array_child(j+1)
            end do

         else if ((symba_plA%status(index1) == DISRUPTION) .and.    &                              
             (symba_plA%status(index2) == DISRUPTION)) then 

            enc_big = plplenc_list%index1(i)
            nchild = symba_plA%nchild(enc_big)
            array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
            do j = 1, nchild
               symba_plA%status(array_child(j)) = INACTIVE
            end do
            ldiscard = .true.
         else if ((symba_plA%status(index1) == SUPERCATASTROPHIC) .and.   &                               
             (symba_plA%status(index2) == SUPERCATASTROPHIC)) then 
            
            enc_big = plplenc_list%index1(i)
            nchild = symba_plA%nchild(enc_big)
            array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
            do j = 1, nchild
               symba_plA%status(array_child(j)) = INACTIVE
            end do
            ldiscard = .true.
         else if ((symba_plA%status(index1) == HIT_AND_RUN) .and.    &                            
             (symba_plA%status(index2) == HIT_AND_RUN)) then 

            enc_big = plplenc_list%index1(i)
            nchild = symba_plA%nchild(enc_big)
            array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
            do j = 1, nchild
               symba_plA%status(array_child(j)) = INACTIVE
            end do
            ldiscard = .true.
         else if ((symba_plA%status(index1) == GRAZE_AND_MERGE) .and.  &                              
             (symba_plA%status(index2) == GRAZE_AND_MERGE)) then 

            enc_big = plplenc_list%index1(i)
            nchild = symba_plA%nchild(enc_big)
            array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
            do j = 1, nchild
               symba_plA%status(array_child(j)) = INACTIVE
            end do
            ldiscard = .true.
         end if
      end if
   end do

   return
   
   end procedure symba_discard_merge_pl
end submodule s_symba_discard_merge_pl
