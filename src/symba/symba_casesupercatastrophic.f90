submodule (symba) s_symba_casesupercatastrophic 
contains
   module procedure symba_casesupercatastrophic 
   !! author: Jennifer L. L. Pouplin and Carlisle A. Wishard
   !!
   !! Supercatastrophic disruption event
   use swiftest
implicit none
   integer(I4B)                       :: nfrag, i, k, index1, index2, frags_added
   integer(I4B)                       :: index1_parent, index2_parent
   integer(I4B)                       :: name1, name2, nstart
   real(DP)                         :: mtot, msun, avg_d, d_p1, d_p2, semimajor_encounter, e, q, semimajor_inward
   real(DP)                         :: rhill_p1, rhill_p2, r_circle, theta, radius1, radius2, r_smallestcircle
   real(DP)                         :: m_rem, m_test, mass1, mass2, enew, eold, a, b, v_col
   real(DP)                         :: x_com, y_com, z_com, vx_com, vy_com, vz_com
   real(DP)                         :: x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag, m1m2_10
   real(DP), dimension(ndim)                :: vnew, xr, mv, l, kk, p

   !temporary
   interface 
    function cross_product_supercatastrophic(ar1,ar2) result(ans)
       use swiftest
       implicit none
       real(DP),dimension(3),intent(in) :: ar1,ar2
       real(DP),dimension(3)         :: ans
    end function cross_product_supercatastrophic
   end interface

! executable code
   
   write(*,*) "entering casesupercatastrophic"
   ! set the maximum number of fragments to be added in a supercatastrophic disruption collision (nfrag)
   nfrag = 10
   ! pull in the information about the two particles involved in the collision  
   index1 = plplenc_list%index1(index_enc)
   index2 = plplenc_list%index2(index_enc)
   index1_parent = symba_pla%index_parent(index1)
   index2_parent = symba_pla%index_parent(index2)
   name1 = symba_pla%helio%swiftest%name(index1)
   name2 = symba_pla%helio%swiftest%name(index2)
   mass1 = symba_pla%helio%swiftest%mass(index1) ! the mass of the first particle in the collision not including all it's children
   mass2 = symba_pla%helio%swiftest%mass(index2)
   radius1 = symba_pla%helio%swiftest%radius(index1)
   radius2 = symba_pla%helio%swiftest%radius(index2)
   msun = symba_pla%helio%swiftest%mass(1)

   ! find com
   x_com = ((x1(1) * m1) + (x2(1) * m2)) / (m1 + m2)
   y_com = ((x1(2) * m1) + (x2(2) * m2)) / (m1 + m2)
   z_com = ((x1(3) * m1) + (x2(3) * m2)) / (m1 + m2)

   vx_com = ((v1(1) * m1) + (v2(1) * m2)) / (m1 + m2)
   vy_com = ((v1(2) * m1) + (v2(2) * m2)) / (m1 + m2)
   vz_com = ((v1(3) * m1) + (v2(3) * m2)) / (m1 + m2)

   ! find collision velocity
   v_col = norm2(v2(:) - v1(:))

   ! find energy pre-frag
   eold = 0.5_DP*(m1*dot_product(v1(:), v1(:)) + m2*dot_product(v2(:), v2(:)))
   xr(:) = x2(:) - x1(:)
   eold = eold - (m1*m2/(sqrt(dot_product(xr(:), xr(:)))))

   write(*, *) "supercatastrophic disruption between particles ", name1, " and ", name2, " at time t = ",t
   
   ! add both particles involved in the collision to mergesub_list
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = supercatastrophic
   mergesub_list%xh(:,nmergesub) = x1(:)
   mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = radius1
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = supercatastrophic
   mergesub_list%xh(:,nmergesub) = x2(:)
   mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass2
   mergesub_list%radius(nmergesub) = radius2

   ! go through the encounter list and look for particles actively encoutering in this timestep
   ! prevent them from having further encounters in this timestep by setting status in plplenc_list to merged
   do k = 1, nplplenc 
    if ((plplenc_list%status(k) == active) .and. &
       ((index1 == plplenc_list%index1(k) .or. index2 == plplenc_list%index2(k)) .or. &
       (index2 == plplenc_list%index1(k) .or. index1 == plplenc_list%index2(k)))) then
        plplenc_list%status(k) = merged
    end if
   end do

   ! set the status of the particles in symba_pla to disruption
   symba_pla%helio%swiftest%status(index1) = supercatastrophic
   symba_pla%helio%swiftest%status(index2) = supercatastrophic

   l(:) = (v2(:) - v1(:)) / norm2(v2(:)-v1(:))
   p(:) = cross_product_supercatastrophic(xr(:) / norm2(xr(:)), l(:))
   kk(:) = cross_product_supercatastrophic(l(:),p(:))

   ! calculate the positions of the new fragments in a circle with a radius large enough to space
   ! all fragments apart by a distance of rhill_p1 + rhill_p2
   rhill_p1 = symba_pla%helio%swiftest%rhill(index1_parent)
   rhill_p2 = symba_pla%helio%swiftest%rhill(index2_parent)
   r_smallestcircle = (rhscale * rhill_p1 + rhscale * rhill_p2) / (2.0_DP * sin(pi /2.0_DP))

   ! check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   call orbel_xv2aeq(x1, v1, msun, semimajor_encounter, e, q)
   ! if they are going to be added interior to this orbit, give a warning
   if (semimajor_inward > (semimajor_encounter - r_smallestcircle)) then
    write(*,*) "warning in symba_casesupercatastrophic: timestep is too large to resolve fragments."
   end if
   ! if not, continue through all possible fragments to be added
    mtot = 0.0_DP ! running total mass of new fragments
    mv = 0.0_DP   ! running sum of m*v of new fragments to be used in com calculation
    frags_added = 0 ! running total number of new fragments
    m1m2_10 = 0.1_DP * (m1 + m2) ! one tenth the total initial mass of the system used to check the size of the fragments
    nstart = nmergeadd

    d_p1 = (3.0_DP * m1) / (4.0_DP * pi * (rad1 ** 3.0_DP))
    d_p2 = (3.0_DP * m2) / (4.0_DP * pi * (rad2 ** 3.0_DP))
    avg_d = ((m1 * d_p1) + (m2 * d_p2)) / (m1 + m2)

       ! if we are adding the first and largest fragment (lr), check to see if its mass is smaller than one tenth the total
       ! mass of the system aka if it is too small to resolve. if so, add a fragment with a mass of one tenth the total mass 
       ! of the system and calculate its radius.
       if ((mres(1) < m1m2_10)) then
        do i = 1, nfrag
         frags_added = frags_added + 1
         nmergeadd = nmergeadd + 1
         mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
         mergeadd_list%status(nmergeadd) = supercatastrophic
         mergeadd_list%ncomp(nmergeadd) = 2
         mergeadd_list%mass(nmergeadd) = m1m2_10
         mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * pi * avg_d))  & 
            ** (1.0_DP / 3.0_DP)
         mtot = mtot + mergeadd_list%mass(nmergeadd) 
        end do 
       end  if
       ! if we are adding the first and largest fragment (lr), check to see if its mass is larger than one tenth the total 
       ! mass of the system aka if it is large enough to resolve. if so, its mass and radius should be taken from 
       ! util_regime.
       if ((mres(1) > m1m2_10)) then
        frags_added = frags_added + 1
        nmergeadd = nmergeadd + 1
        mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
        mergeadd_list%status(nmergeadd) = supercatastrophic
        mergeadd_list%ncomp(nmergeadd) = 2
        mergeadd_list%mass(nmergeadd) = mres(1)
        mergeadd_list%radius(nmergeadd) = rres(1)
        mtot = mtot + mergeadd_list%mass(nmergeadd) 
        do i = 2, nfrag
         frags_added = frags_added + 1
         nmergeadd = nmergeadd + 1
         mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
         mergeadd_list%status(nmergeadd) = supercatastrophic
         mergeadd_list%ncomp(nmergeadd) = 2
         mergeadd_list%mass(nmergeadd) = (m1 + m2 - mres(1)) / (nfrag - 1.0_DP)
         mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * pi * avg_d))  & 
            ** (1.0_DP / 3.0_DP)  
         mtot = mtot + mergeadd_list%mass(nmergeadd)
        end do
       end if 


       !end if
       ! if we are adding more than one fragment
       !if ((i > 1) .and. (mres(1) > m1m2_10)) then
        ! m_rem is the mass needed to be "made up for" in fragments, mres(1) is the mass of the largest fragments 
        ! that has already been added, and m1 and m2 are the masses of the original particles involved in the collision.
       !   m_rem = (m1 + m2) - (mergeadd_list%mass(nmergeadd))
        ! check if these fragments will be large enough to be resolved
       !   if (m_rem > (1.0_DP / 10.0_DP)*mres(1))) then

         ! if yes, add a fragment using durda et al 2007 figure 2 supercatastrophic: n = (1.5e5)e(-1.3*d) for the mass
       !    frags_added = frags_added + 1
       !    nmergeadd = nmergeadd + 1
       !    mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
       !    mergeadd_list%status(nmergeadd) = supercatastrophic
       !    mergeadd_list%ncomp(nmergeadd) = 2
       !    mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 1) 

         ! create a "test mass" using durda et al 2007 figure 2 supercatastrophic: n = (1.5e5)e(-1.3*d)
         !m_test = (((- 1.0_DP / 2.6_DP) * log(i / (1.5_DP * 10.0_DP ** 5))) ** 3.0_DP) * ((4.0_DP / 3.0_DP) &
         !   * pi * avg_d)
         ! if the test mass is smaller than the mass that needs to be "made up for", add it.
         !if (m_test < m_rem) then
         !   mergeadd_list%mass(nmergeadd) = m_test
         ! if not, aka if the test mass is too large, then add a fragment with a mass equal to the difference between
         ! the sum of the mass of the parents and the total mass already added. 
         !else
         !   mergeadd_list%mass(nmergeadd) = (m1 + m2) - mtot 
         !end if

         ! calculate the radius of the fragments using the weighted average density of the parents. 
       !    mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * pi * avg_d))  & 
       !       ** (1.0_DP / 3.0_DP) 
       !    mtot = mtot + mergeadd_list%mass(nmergeadd)
       !   else 
       !    mergeadd_list%mass(nmergeadd) = mergeadd_list%mass(nmergeadd) + m_rem
       !    mergeadd_list%radius(nmergeadd) = (((3.0_DP/4.0_DP) * pi) * (mergeadd_list%mass(nmergeadd) / avg_d)) &
       !       ** (1.0_DP / 3.0_DP)
       !   end if
       !end if

   r_circle = (rhscale * rhill_p1 + rhscale * rhill_p2) / (2.0_DP * sin(pi / frags_added))
   theta = (2.0_DP * pi) / frags_added

   do i=1, frags_added

    ! increment around the circle for positions of fragments
    x_frag = (r_circle * cos(theta * i))*l(1) + (r_circle * sin(theta * i))*p(1) + x_com
    y_frag = (r_circle * cos(theta * i))*l(2) + (r_circle * sin(theta * i))*p(2) + y_com
    z_frag = (r_circle * cos(theta * i))*l(3) + (r_circle * sin(theta * i))*p(3) + z_com

    a = v_col * (m1 + m2) * (1.0_DP / mergeadd_list%mass(nstart + i))

    vx_frag = ((a * cos(theta * i))*l(1)) + ((a * sin(theta * i))*p(1)) + vx_com
    vy_frag = ((a * cos(theta * i))*l(2)) + ((a * sin(theta * i))*p(2)) + vy_com
    vz_frag = ((a * cos(theta * i))*l(3)) + ((a * sin(theta * i))*p(3)) + vz_com

       !conservation of angular momentum for velocities of fragments
       !a = (((x1(2) * v1(3) * m1) - (x1(3) * v1(2) * m1)) + ((x2(2) * v2(3) * m2) - (x2(3) * v2(2) * m2))) &
       !   / mergeadd_list%mass(nmergeadd)
       !b = (((x1(3) * v1(1) * m1) - (x1(1) * v1(3) * m1)) + ((x2(3) * v2(1) * m2) - (x2(1) * v2(3) * m2))) &
       !   / mergeadd_list%mass(nmergeadd)
       !vx_frag = ((1.0_DP / frags_added) * (b / z_frag)) - vbs(1)
       !vy_frag = ((1.0_DP / frags_added) * (-a / z_frag)) - vbs(2)
       !vz_frag = vz_com - vbs(3)

    mergeadd_list%xh(1,nstart + i) = x_frag
    mergeadd_list%xh(2,nstart + i) = y_frag 
    mergeadd_list%xh(3,nstart + i) = z_frag                                
    mergeadd_list%vh(1,nstart + i) = vx_frag
    mergeadd_list%vh(2,nstart + i) = vy_frag
    mergeadd_list%vh(3,nstart + i) = vz_frag
    ! tracking linear momentum.
    mv = mv + (mergeadd_list%mass(nstart + i) * mergeadd_list%vh(:,nstart + i))
   end do

   write(*, *) "number of fragments added: ", frags_added
   ! calculate energy after frag                                             
   vnew(:) = mv / mtot    ! com of new fragments                   
   enew = 0.5_DP*mtot*dot_product(vnew(:), vnew(:))
   eoffset = eoffset + eold - enew

   ! update fragmax to account for new fragments
   fragmax = fragmax + frags_added
   write(*,*) "leaving casesupercatastrophic"

   return 
 end procedure symba_casesupercatastrophic
end submodule s_symba_casesupercatastrophic
