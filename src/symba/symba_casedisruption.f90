submodule (symba) s_symba_casedisruption 
contains
   module procedure symba_casedisruption 
   !! author: Jennifer L. L. Pouplin and Carlisle A. Wishard
   !!
   !! Disruptive Collision
   use swiftest
implicit none
   integer(I4B)                       :: nfrag, i, k, index1, index2, frags_added
   integer(I4B)                       :: index1_parent, index2_parent
   integer(I4B)                       :: name1, name2, nstart
   real(DP)                         :: mtot, msun, avg_d, d_p1, d_p2, semimajor_encounter, e, q, semimajor_inward
   real(DP)                         :: rhill_p1, rhill_p2, r_circle, theta, radius1, radius2, r_smallestcircle
   real(DP)                         :: m_rem, m_test, mass1, mass2, enew, eold, a, b, v_col
   real(DP)                         :: x_com, y_com, z_com, vx_com, vy_com, vz_com
   real(DP)                         :: x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag
   real(DP), dimension(ndim)                :: vnew, xr, mv, l, kk, p

   !temporary
   interface 
    function cross_product_disruption(ar1,ar2) result(ans)
       use swiftest
       implicit none
       real(DP),dimension(3),intent(in) :: ar1,ar2
       real(DP),dimension(3)         :: ans
    end function cross_product_disruption
   end interface

! executable code

   write(*,*) "entering casedisruption"
   ! set the maximum number of fragments to be added in a disruption collision (nfrag)
   nfrag = 5 
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
   
   write(*, *) "disruption between particles ", name1, " and ", name2, " at time t = ",t
   
   ! add both particles involved in the collision to mergesub_list
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = disruption
   mergesub_list%xh(:,nmergesub) = x1(:)
   mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = radius1
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = disruption
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
   symba_pla%helio%swiftest%status(index1) = disruption
   symba_pla%helio%swiftest%status(index2) = disruption

   l(:) = (v2(:) - v1(:)) / norm2(v2(:)-v1(:))
   p(:) = cross_product_disruption(xr(:) / norm2(xr(:)), l(:))
   kk(:) = cross_product_disruption(l(:),p(:))

   rhill_p1 = symba_pla%helio%swiftest%rhill(index1_parent)
   rhill_p2 = symba_pla%helio%swiftest%rhill(index2_parent)
   r_smallestcircle = (rhscale * rhill_p1 + rhscale * rhill_p2) / (2.0_DP * sin(pi / 2.0_DP))

   ! check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   call orbel_xv2aeq(x1, v1, msun, semimajor_encounter, e, q)
   ! if they are going to be added interior to this orbit, give a warning
   if (semimajor_inward > (semimajor_encounter - r_smallestcircle)) then
    write(*,*) "warning in symba_casedisruption: timestep is too large to resolve fragments."
   end if
   ! if not, continue through all possible fragments to be added
   mtot = 0.0_DP ! running total mass of new fragments
   mv = 0.0_DP   ! running sum of m*v of new fragments to be used in com calculation
   frags_added = 0 ! running total number of new fragments
   nstart = nmergeadd ! start of new fragments in mergeadd_list

   d_p1 = (3.0_DP * m1) / (4.0_DP * pi * (rad1 ** 3.0_DP))
   d_p2 = (3.0_DP * m2) / (4.0_DP * pi * (rad2 ** 3.0_DP))
   avg_d = ((m1 * d_p1) + (m2 * d_p2)) / (m1 + m2)

   !do i = 1, nfrag
    ! if we are adding the first and largest fragment (lr), it's mass and radius should be taken from 
    ! util_regime while it's position and velocity should be calculated on the circle of radius 
    ! r_circle as described above.
    !if (i == 1) then
       frags_added = frags_added + 1
       nmergeadd = nmergeadd + 1
       mergeadd_list%status(nmergeadd) = disruption
       mergeadd_list%ncomp(nmergeadd) = 2
       mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
       mergeadd_list%mass(nmergeadd) = mres(1)
       mergeadd_list%radius(nmergeadd) = rres(1)
       mtot = mtot + mergeadd_list%mass(nmergeadd)                   
    !end if
    ! if we are adding the second fragment (slr), it's mass and radius should be taken from 
    ! util_regime while it's position and velocity should be calculated on the circle of 
    ! radius r_circle as described above.
    if ((mres(2) > (1.0_DP / 3.0_DP)*mres(1))) then
       write(*,*) "casedisruption 1st if"
       ! frags_added is the actual number of fragments added to the simulation vs nfrag which is the total possible
       frags_added = frags_added + 1
       nmergeadd = nmergeadd + 1
       mergeadd_list%status(nmergeadd) = disruption
       mergeadd_list%ncomp(nmergeadd) = 2
       mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
       mergeadd_list%mass(nmergeadd) = mres(2)
       mergeadd_list%radius(nmergeadd) = rres(2)
       mtot = mtot + mergeadd_list%mass(nmergeadd)
       do i = 3, nfrag
        write(*,*) "casedisruption 1st do"
        frags_added = frags_added + 1
        nmergeadd = nmergeadd + 1
        mergeadd_list%status(nmergeadd) = disruption
        mergeadd_list%ncomp(nmergeadd) = 2
        mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
        m_rem = (m1 + m2) - (mres(1) + mres(2))
        mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 1) 
        mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * pi * avg_d))  & 
         ** (1.0_DP / 3.0_DP) 
        mtot = mtot + mergeadd_list%mass(nmergeadd) 
       end do                  
    end if

    if ((mres(2) < (1.0_DP / 3.0_DP)*mres(1))) then
       write(*,*) "casedisruption 2nd if"   
       do i = 2, nfrag
        write(*,*) "casedisruption 2nd do"
        m_rem = (m1 + m2) - mres(1)
        frags_added = frags_added + 1
        nmergeadd = nmergeadd + 1
        mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
        mergeadd_list%status(nmergeadd) = disruption
        mergeadd_list%ncomp(nmergeadd) = 2
        mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 1) 
        mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * pi * avg_d))  & 
         ** (1.0_DP / 3.0_DP)  
        mtot = mtot + mergeadd_list%mass(nmergeadd)
       end do
    end if

    ! if we are doing more fragments
    !if ((i > 2) .and. (mres(2) > (1.0_DP / 3.0_DP)*mres(1))) then
       ! m_rem is the mass needed to be "made up for" in fragments, mres(1) and mres(2) are the mass of the largest 
       ! and second largest fragments that have already been added, and m1 and m2 are the masses of the original 
       ! particles involved in the collision.
       !frags_added = frags_added + 1
       !nmergeadd = nmergeadd + 1
       !mergeadd_list%status(nmergeadd) = disruption
       !mergeadd_list%ncomp(nmergeadd) = 2
       !mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
       !m_rem = (m1 + m2) - (mres(1) + mres(2))
       !mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 1) 
       !mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * pi * avg_d))  & 
       !   ** (1.0_DP / 3.0_DP) 
       !mtot = mtot + mergeadd_list%mass(nmergeadd) 

        !write(*,*) "casedisruption mres(1) and mres(2)", mres(1), mres(2)

        ! check if these fragments will be large enough to be resolved    
        !if (m_rem > (m2) / 100.0_DP) then
         ! if yes, add a fragment using durda et al 2007 figure 2 supercatastrophic: n = (1.5e5)e(-1.3*d) for the mass
        
        ! create a "test mass" using durda et al 2007 figure 2 supercatastrophic: n = (1.5e5)e(-1.3*d)
        !m_test = (((- 1.0_DP / 2.6_DP) * log(i / (1.5_DP * 10.0_DP ** 5.0_DP))) ** 3.0_DP) * ((4.0_DP / 3.0_DP) &
        !   * pi * avg_d)
        ! if the test mass is smaller than the mass that needs to be "made up for", add it. 
        !if (m_test < m_rem) then
        !    mergeadd_list%mass(nmergeadd) = m_test
         ! if not, aka if the test mass is too large, then add a fragment with a mass equal to the difference between
         ! the sum of the mass of the parents and the total mass already added.   
        !else
        !   mergeadd_list%mass(nmergeadd) = m_rem !(m1 + m2) - mtot 
        !end if 

        ! calculate the radius of the fragments using the weighted average density of the parents. 
        !mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * pi * avg_d))  & 
        !   ** (1.0_DP / 3.0_DP) 
        !mtot = mtot + mergeadd_list%mass(nmergeadd) 

        ! if these fragments will not be large enough to be resolved, add the remaining mass that we need to 
        ! "make up for" to the mass of the most recent fragment and recalculate the radius. 
        !else 
        !   mergeadd_list%mass(nmergeadd) = mergeadd_list%mass(nmergeadd) + m_rem
        !   mergeadd_list%radius(nmergeadd) = (((3.0_DP/4.0_DP) * pi) * (mergeadd_list%mass(nmergeadd) / avg_d)) &
        !    ** (1.0_DP / 3.0_DP)                                    
        !end if
    !end if
   !end do

    ! calculate the positions of the new fragments in a circle with a radius large enough to space
    ! all fragments apart by a distance of rhill_p1 + rhill_p2
   r_circle = (2.0_DP * rhill_p1 + 2.0_DP * rhill_p2) / (2.0_DP * sin(pi / frags_added))
   theta = (2.0_DP * pi) / frags_added

   do i=1, frags_added
    write(*,*) "casedisruption 3rd do"

       !write(*,*) "casedisruption mfrag/mtot", mergeadd_list%mass(nstart + i) / (m1 + m2)

       ! increment around the circle for positions of fragments
    x_frag = (r_circle * cos(theta * i))*l(1) + (r_circle * sin(theta * i))*p(1) + x_com
    y_frag = (r_circle * cos(theta * i))*l(2) + (r_circle * sin(theta * i))*p(2) + y_com
    z_frag = (r_circle * cos(theta * i))*l(3) + (r_circle * sin(theta * i))*p(3) + z_com

    a = v_col * (m1 + m2) * (1.0_DP / mergeadd_list%mass(nstart + i))

    vx_frag = ((a * cos(theta * i))*l(1)) + ((a * sin(theta * i))*p(1)) + vx_com
    vy_frag = ((a * cos(theta * i))*l(2)) + ((a * sin(theta * i))*p(2)) + vy_com
    vz_frag = ((a * cos(theta * i))*l(3)) + ((a * sin(theta * i))*p(3)) + vz_com

    write(*,*) "casedisruption vx_frag", vx_frag
    write(*,*) "casedisruption vy_frag", vy_frag
    write(*,*) "casedisruption vz_frag", vz_frag

       !vx_frag = ((1.0_DP / frags_added) * (1.0_DP / mergeadd_list%mass(nstart + i)) * ((m1 * v1(1)) + (m2 * v2(1)))) + vx_com !- vbs(1)
       !vy_frag = ((1.0_DP / frags_added) * (1.0_DP / mergeadd_list%mass(nstart + i)) * ((m1 * v1(2)) + (m2 * v2(2)))) + vy_com !- vbs(2)
       !vz_frag = ((1.0_DP / frags_added) * (1.0_DP / mergeadd_list%mass(nstart + i)) * ((m1 + v1(3)) + (m2 * v2(3)))) + vz_com !- vbs(3)

       !write(*,*) "casedisruption vx_frag", vx_frag
       !write(*,*) "casedisruption vy_frag", vy_frag
       !write(*,*) "casedisruption vz_frag", vz_frag

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
   write(*,*) "leaving casedisruption"
   return 
 end procedure symba_casedisruption
end submodule s_symba_casedisruption
