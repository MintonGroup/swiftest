submodule (symba) s_symba_fragmentation 
contains
   module procedure symba_fragmentation 
   !! author: Jennifer L. L. Pouplin and Carlisle A. Wishard
   !!
   !! Generate new particles resulting from collisions
use swiftest
implicit none
 
   integer(I4B)             :: model, nres, i, itarg, iproj
   real(DP), dimension(3)       :: mres, rres
   real(DP), dimension(ndim, 3)   :: pres, vres
   integer(I4B)             :: regime 
   integer(I4B)             :: index1, index2, index1_child, index2_child, index1_parent, index2_parent
   integer(I4B)             :: name1, name2, index_big1, index_big2, stat1, stat2
   real(DP)               :: r2, rlim, rlim2, vdotr, tcr2, dt2, a, e, q
   real(DP)               :: rad1, rad2, m1, m2, den1, den2, vol1, vol2, vchild, dentarg, denproj, dentot, mcenter
   real(DP)               :: mass1, mass2, mmax, mtmp, mtot, m1_si, m2_si
   real(DP), dimension(ndim)    :: xr, vr, x1, v1, x2, v2, x1_si, x2_si, v1_si, v2_si, xproj, xtarg, vproj, vtarg
   real(DP)               :: den1_si, den2_si, rad1_si, rad2_si, rproj, rtarg
   logical(lgt)             :: lfrag_add, lmerge
   integer(I4B), dimension(npl)   :: array_index1_child, array_index2_child
   real(DP)               :: mlr, mslr, mtarg, mproj
   !real(DP)               :: k2 = 2.959122082855911e-4 ! in si units
   !real(DP)               :: msun = 1.98847e30 ! in si units
   !real(DP)               :: au = 1.495978707e11 ! in si units
   !real(DP)               :: year = 3.154e7 ! in si units


! executable code

   lmerge = .false.
   lfrag_add = .false.
   ! model 2 is the model for collresolve_resolve (ls12)
   model = 2

   index1 = plplenc_list%index1(index_enc)
   index2 = plplenc_list%index2(index_enc)

   rlim = symba_pla%helio%swiftest%radius(index1) + symba_pla%helio%swiftest%radius(index2)
   xr(:) = symba_pla%helio%swiftest%xh(:,index2) - symba_pla%helio%swiftest%xh(:,index1)
   r2 = dot_product(xr(:), xr(:))
   rlim2 = rlim*rlim
   ! checks if bodies are actively colliding in this time step
   if (rlim2 >= r2) then 
      lfrag_add = .true.
   ! if they are not actively colliding in  this time step, 
   !checks if they are going to collide next time step based on velocities and q
   else 
      vr(:) = symba_pla%helio%swiftest%vb(:,index2) - symba_pla%helio%swiftest%vb(:,index1)
      vdotr = dot_product(xr(:), vr(:))
      if (plplenc_list%lvdotr(index_enc) .and. (vdotr > 0.0_DP)) then 
         tcr2 = r2/dot_product(vr(:), vr(:))
         dt2 = dt*dt
         if (tcr2 <= dt2) then
            mtot = symba_pla%helio%swiftest%mass(index1) + symba_pla%helio%swiftest%mass(index2)
            call orbel_xv2aeq(xr(:), vr(:), mtot, a, e, q)
            if (q < rlim) lfrag_add = .true.
         end if
         ! if no collision is going to happen, write as close encounter, not  merger
         if (.not. lfrag_add) then
            if (encounter_file /= "") then
               name1 = symba_pla%helio%swiftest%name(index1)
               m1 = symba_pla%helio%swiftest%mass(index1)
               rad1 = symba_pla%helio%swiftest%radius(index1)
               x1(:) = symba_pla%helio%swiftest%xh(:,index1)
               v1(:) = symba_pla%helio%swiftest%vb(:,index1) - vbs(:)
               name2 = symba_pla%helio%swiftest%name(index2)
               m2 = symba_pla%helio%swiftest%mass(index2)
               rad2 = symba_pla%helio%swiftest%radius(index2)
               x2(:) = symba_pla%helio%swiftest%xh(:,index2)
               v2(:) = symba_pla%helio%swiftest%vb(:,index2) - vbs(:)

               call io_write_encounter(t, name1, name2, m1, m2, rad1, rad2, x1(:), x2(:), &
                  v1(:), v2(:), encounter_file, out_type)
            end if
         end if
      end if
   end if

   nres = 2
   if (lfrag_add) then 
      symba_pla%lmerged(index1) = .true.
      symba_pla%lmerged(index2) = .true.
      index1_parent = symba_pla%index_parent(index1)
      m1 = symba_pla%helio%swiftest%mass(index1_parent)
      mass1 = m1 
      x1(:) = m1*symba_pla%helio%swiftest%xh(:,index1_parent)
      v1(:) = m1*symba_pla%helio%swiftest%vb(:,index1_parent)
      mmax = m1
      name1 = symba_pla%helio%swiftest%name(index1_parent)
      index_big1 = index1_parent
      stat1 = symba_pla%helio%swiftest%status(index1_parent)
      array_index1_child(1:npl) = symba_pla%index_child(1:npl,index1_parent)
      
      vol1 =  ((4.0_DP / 3.0_DP) * pi * symba_pla%helio%swiftest%radius(index1_parent)**3.0_DP)
      do i = 1, symba_pla%nchild(index1_parent) ! initialize an array of children
         index1_child = array_index1_child(i)
         mtmp = symba_pla%helio%swiftest%mass(index1_child)
         vchild = ((4.0_DP / 3.0_DP) * pi * symba_pla%helio%swiftest%radius(index1_child)**3.0_DP)
         vol1 = vol1 + vchild
         if (mtmp > mmax) then
            mmax = mtmp
            name1 = symba_pla%helio%swiftest%name(index1_child)
            index_big1 = index1_child
            stat1 = symba_pla%helio%swiftest%status(index1_child)
         end if
         m1 = m1 + mtmp
         x1(:) = x1(:) + mtmp*symba_pla%helio%swiftest%xh(:,index1_child)
         v1(:) = v1(:) + mtmp*symba_pla%helio%swiftest%vb(:,index1_child)
      end do
      den1 =  m1 / vol1
      rad1 = ((3.0_DP * m1) / (den1 * 4.0_DP * pi)) ** (1.0_DP / 3.0_DP)
      x1(:) = x1(:)/m1
      v1(:) = v1(:)/m1

      index2_parent = symba_pla%index_parent(index2)
      m2 = symba_pla%helio%swiftest%mass(index2_parent)
      mass2 = m2
      rad2 = symba_pla%helio%swiftest%radius(index2_parent)
      x2(:) = m2*symba_pla%helio%swiftest%xh(:,index2_parent)
      v2(:) = m2*symba_pla%helio%swiftest%vb(:,index2_parent)
      mmax = m2
      name2 = symba_pla%helio%swiftest%name(index2_parent)
      index_big2 = index2_parent
      stat2 = symba_pla%helio%swiftest%status(index2_parent)
      array_index2_child(1:npl) = symba_pla%index_child(1:npl,index2_parent)

      vol2 = ((4.0_DP / 3.0_DP) * pi * symba_pla%helio%swiftest%radius(index2_parent)**3.0_DP)
      do i = 1, symba_pla%nchild(index2_parent)
         index2_child = array_index2_child(i)
         mtmp = symba_pla%helio%swiftest%mass(index2_child)
         vchild =  ((4.0_DP / 3.0_DP) * pi * symba_pla%helio%swiftest%radius(index2_child)**3.0_DP)
         vol2 = vol2 + vchild
         if (mtmp > mmax) then
            mmax = mtmp
            name2 = symba_pla%helio%swiftest%name(index2_child)
            index_big2 = index2_child
            stat2 = symba_pla%helio%swiftest%status(index2_child)
         end if
         m2 = m2 + mtmp
         x2(:) = x2(:) + mtmp*symba_pla%helio%swiftest%xh(:,index2_child)
         v2(:) = v2(:) + mtmp*symba_pla%helio%swiftest%vb(:,index2_child)
      end do
      gu = gc / (du2m**3 / (mu2kg * tu2s**2))
      den2 =  m2 / vol2
      rad2 = ((3.0_DP * m2) / (den2 * 4.0_DP * pi)) ** (1.0_DP / 3.0_DP)
      x2(:) = x2(:)/m2
      v2(:) = v2(:)/m2

      m1_si = (m1 / gu) * mu2kg 
      m2_si = (m2 / gu) * mu2kg
      rad1_si = rad1 * du2m
      rad2_si = rad2 * du2m
      x1_si(:) = x1(:) * du2m
      x2_si(:) = x2(:) * du2m
      v1_si(:) = v1(:) * du2m / tu2s
      v2_si(:) = v2(:) * du2m / tu2s
      den1_si = (den1 / gu) * mu2kg / (du2m ** 3.0_DP)
      den2_si = (den2 / gu) * mu2kg / (du2m ** 3.0_DP)

      mres(:) = 0.0_DP
      rres(:) = 0.0_DP
      pres(:,:) = 0.0_DP
      vres(:,:) = 0.0_DP

      if (m1_si > m2_si) then 
         itarg = index1
         iproj = index2
         dentarg = den1_si
         denproj = den2_si
         mtarg = m1_si
         mproj = m2_si
         rtarg = rad1_si
         rproj = rad2_si
         xtarg(:) = x1_si(:)
         xproj(:) = x2_si(:)
         vtarg(:) = v1_si(:)
         vproj(:) = v2_si(:)
      else
         itarg = index2
         iproj = index1
         dentarg = den2_si
         denproj = den1_si
         mtarg = m2_si
         mproj = m1_si
         rtarg = rad2_si
         rproj = rad1_si
         xtarg(:) = x2_si(:)
         xproj(:) = x1_si(:)
         vtarg(:) = v2_si(:)
         vproj(:) = v1_si(:)
      end if
      mtot = m1_si + m2_si
      dentot = (m1_si *den1_si +m2_si*den2_si )/ mtot
      mcenter = symba_pla%helio%swiftest%mass(1) * mu2kg / gu

      !regime = collresolve_resolve(model,mtarg,mproj,rtarg,rproj,xtarg,xproj, vtarg,vproj, nres, mres, rres, pres, vres)

      call util_regime(mcenter, mtarg, mproj, rtarg, rproj, xtarg, xproj, vtarg, vproj, dentarg, denproj, regime, mlr, mslr)

      mres(1) = mlr
      mres(2) = mslr 
      mres(3) = mtot - mlr - mslr
      rres(1) = (3.0_DP * mres(1)  / (4.0_DP * pi * dentarg)) *(1.0_DP/3.0_DP)
      rres(1) = (3.0_DP * mres(1)  / (4.0_DP * pi * dentarg)) ** (1.0_DP/3.0_DP)
      rres(2) = (3.0_DP * mres(2)  / (4.0_DP * pi * denproj)) ** (1.0_DP/3.0_DP)
      rres(3) = (3.0_DP * mres(2)  / (4.0_DP * pi * dentot)) ** (1.0_DP/3.0_DP)

      mres(:) = (mres(:) / mu2kg) * gu
      rres(:) = rres(:) / du2m

      call symba_caseresolve(t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
      npl, symba_pla, nplplenc, plplenc_list, regime, nplmax, ntpmax, fragmax, mres, rres, array_index1_child, &
      array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)

   end if 
   return

   end procedure symba_fragmentation
end submodule s_symba_fragmentation
