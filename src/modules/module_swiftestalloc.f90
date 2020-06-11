module module_swiftestalloc
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Module containing subroutines that allocate and initialize the Swiftest data structures
   !!
   use swiftest_globals
   use swiftest_data_structures

   contains 




   subroutine helio_pl_allocate(helio_plA, npl)
      use module_helio
      implicit none

      integer(I4B), intent(in)            :: npl
      type(helio_pl), intent(inout)        :: helio_plA

      if (npl <= 0) return
      allocate(helio_plA%ah(NDIM,npl))
      allocate(helio_plA%ahi(NDIM,npl))
      helio_plA%ah = 0.0_DP
      helio_plA%ahi = 0.0_DP
      call helio_plA%swiftest%alloc(npl)
      return
   end subroutine helio_pl_allocate


   subroutine symba_pl_allocate(symba_plA, npl)
      use module_symba
      implicit none

      integer(I4B), intent(in)            :: npl
      type(symba_pl), intent(inout)        :: symba_plA

      if (npl <= 0) return
      allocate(symba_plA%lmerged(npl))
      allocate(symba_plA%nplenc(npl))
      allocate(symba_plA%ntpenc(npl))
      allocate(symba_plA%levelg(npl))
      allocate(symba_plA%levelm(npl))
      allocate(symba_plA%nchild(npl))
      allocate(symba_plA%isperi(npl))
      allocate(symba_plA%peri(npl))
      allocate(symba_plA%atp(npl))
      allocate(symba_plA%index_parent(npl))
      allocate(symba_plA%index_child(npl,npl))

      symba_plA%lmerged = .false.
      symba_plA%nplenc = 0
      symba_plA%ntpenc = 0
      symba_plA%levelg = 0
      symba_plA%levelm = 0
      symba_plA%nchild = 0
      symba_plA%isperi = 0
      symba_plA%peri = 0.0_DP
      symba_plA%atp = 0.0_DP
      symba_plA%index_parent = 1
      symba_plA%index_child = 1
      call helio_pl_allocate(symba_plA%helio,npl)
      return
   end subroutine symba_pl_allocate

   subroutine symba_plplenc_allocate(plplenc_list, nplplenc)
      use module_symba
      implicit none

      integer(I4B), intent(in)                :: nplplenc
      type(symba_plplenc), intent(inout)        :: plplenc_list

      if (nplplenc <= 0) return
      allocate(plplenc_list%lvdotr(nplplenc))
      allocate(plplenc_list%status(nplplenc))
      allocate(plplenc_list%level(nplplenc))
      allocate(plplenc_list%index1(nplplenc))
      allocate(plplenc_list%index2(nplplenc))
      allocate(plplenc_list%enc_child(nplplenc))
      allocate(plplenc_list%enc_parent(nplplenc))

      plplenc_list%lvdotr = .false.
      plplenc_list%status = 0
      plplenc_list%level = 1
      plplenc_list%index1 = 1
      plplenc_list%index2 = 1
      plplenc_list%enc_child = 1
      plplenc_list%enc_parent = 1
      return
   end subroutine symba_plplenc_allocate

   subroutine symba_merger_allocate(mergeadd_list, nmergeadd)
      use module_symba
      implicit none

      integer(I4B), intent(in)                :: nmergeadd
      type(symba_merger), intent(inout)        :: mergeadd_list

      if (nmergeadd <= 0) return
      allocate(mergeadd_list%name(nmergeadd))
      allocate(mergeadd_list%index_ps(nmergeadd))
      allocate(mergeadd_list%status(nmergeadd))
      allocate(mergeadd_list%ncomp(nmergeadd))
      allocate(mergeadd_list%xh(NDIM,nmergeadd))
      allocate(mergeadd_list%vh(NDIM,nmergeadd))
      allocate(mergeadd_list%mass(nmergeadd))
      allocate(mergeadd_list%radius(nmergeadd))

      mergeadd_list%name = 0
      mergeadd_list%index_ps = 1
      mergeadd_list%status = 0
      mergeadd_list%ncomp = 0
      mergeadd_list%xh = 0.0_DP
      mergeadd_list%vh = 0.0_DP
      mergeadd_list%mass = 0.0_DP
      mergeadd_list%radius = 0.0_DP

      return
   end subroutine symba_merger_allocate

   subroutine helio_tp_allocate(helio_tpA, ntp)
      use module_helio
      implicit none

      integer(I4B), intent(in)            :: ntp
      type(helio_tp), intent(inout)       :: helio_tpA

      if (ntp <= 0) return
      allocate(helio_tpA%ah(NDIM,ntp))
      allocate(helio_tpA%ahi(NDIM,ntp))

      helio_tpA%ah = 0.0_DP
      helio_tpA%ahi = 0.0_DP
      call helio_tpA%swiftest%alloc(ntp)

      return
   end subroutine helio_tp_allocate


   subroutine symba_tp_allocate(symba_tpA, ntp)
      use module_symba
      implicit none

      integer(I4B), intent(in)            :: ntp
      type(symba_tp), intent(inout)       :: symba_tpA

      if (ntp <= 0) return

      allocate(symba_tpA%nplenc(ntp))
      allocate(symba_tpA%levelg(ntp))
      allocate(symba_tpA%levelm(ntp))

      symba_tpA%nplenc = 0
      symba_tpA%levelg = 0
      symba_tpA%levelm = 0
      call helio_tp_allocate(symba_tpA%helio,ntp)
      return
   end subroutine symba_tp_allocate

   subroutine symba_pltpenc_allocate(pltpenc_list, npltpenc)
      use module_symba
      implicit none

      integer(I4B), intent(in)                :: npltpenc
      type(symba_pltpenc), intent(inout)        :: pltpenc_list

      if (npltpenc <= 0) return

      allocate(pltpenc_list%lvdotr(npltpenc))
      allocate(pltpenc_list%status(npltpenc))
      allocate(pltpenc_list%level(npltpenc))
      allocate(pltpenc_list%indexpl(npltpenc))
      allocate(pltpenc_list%indextp(npltpenc))

      pltpenc_list%lvdotr = .false.
      pltpenc_list%status = 0
      pltpenc_list%level = 0
      pltpenc_list%indexpl = 1
      pltpenc_list%indextp = 1
      return
   end subroutine symba_pltpenc_allocate

!___________________________





   subroutine helio_pl_deallocate(helio_plA)
      use module_helio
      implicit none

      type(helio_pl), intent(inout)        :: helio_plA

      deallocate(helio_plA%ah)
      deallocate(helio_plA%ahi)
      call helio_plA%swiftest%dealloc()
      return
   end subroutine helio_pl_deallocate


   subroutine symba_pl_deallocate(symba_plA)
      use module_symba
      implicit none

      type(symba_pl), intent(inout)        :: symba_plA

      deallocate(symba_plA%lmerged)
      deallocate(symba_plA%nplenc)
      deallocate(symba_plA%ntpenc)
      deallocate(symba_plA%levelg)
      deallocate(symba_plA%levelm)
      deallocate(symba_plA%nchild)
      deallocate(symba_plA%isperi)
      deallocate(symba_plA%peri)
      deallocate(symba_plA%atp)
      deallocate(symba_plA%index_parent)
      deallocate(symba_plA%index_child)
      call helio_pl_deallocate(symba_plA%helio)
      return
   end subroutine symba_pl_deallocate

   subroutine symba_plplenc_deallocate(plplenc_list)
      use module_symba
      implicit none

      type(symba_plplenc), intent(inout)        :: plplenc_list

      deallocate(plplenc_list%lvdotr)
      deallocate(plplenc_list%status)
      deallocate(plplenc_list%level)
      deallocate(plplenc_list%index1)
      deallocate(plplenc_list%index2)
      deallocate(plplenc_list%enc_child)
      deallocate(plplenc_list%enc_parent)
      return
   end subroutine symba_plplenc_deallocate

   subroutine symba_merger_deallocate(mergeadd_list)
      use module_symba
      implicit none

      type(symba_merger), intent(inout)        :: mergeadd_list

      deallocate(mergeadd_list%name)
      deallocate(mergeadd_list%index_ps)
      deallocate(mergeadd_list%status)
      deallocate(mergeadd_list%ncomp)
      deallocate(mergeadd_list%xh)
      deallocate(mergeadd_list%vh)
      deallocate(mergeadd_list%mass)
      deallocate(mergeadd_list%radius)
      return
   end subroutine symba_merger_deallocate

   subroutine helio_tp_deallocate(helio_tpA)
      use module_helio
      implicit none

      type(helio_tp), intent(inout)        :: helio_tpA

      deallocate(helio_tpA%ah)
      deallocate(helio_tpA%ahi)
      call helio_tpA%swiftest%dealloc()
      return
   end subroutine helio_tp_deallocate


   subroutine symba_tp_deallocate(symba_tpA)
      use module_symba
      implicit none

      type(symba_tp), intent(inout)        :: symba_tpA

      deallocate(symba_tpA%nplenc)
      deallocate(symba_tpA%levelg)
      deallocate(symba_tpA%levelm)

      return
   end subroutine symba_tp_deallocate

   subroutine symba_pltpenc_deallocate(pltpenc_list)
      use module_symba
      implicit none

      type(symba_pltpenc), intent(inout)        :: pltpenc_list

      deallocate(pltpenc_list%lvdotr)
      deallocate(pltpenc_list%status)
      deallocate(pltpenc_list%level)
      deallocate(pltpenc_list%indexpl)
      deallocate(pltpenc_list%indextp)
      return
   end subroutine symba_pltpenc_deallocate


end module module_swiftestalloc




