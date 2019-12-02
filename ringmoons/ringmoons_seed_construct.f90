!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_construct
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Constructs the satellite seeds. Given the number of seeds, this will find the mass of each seed needed to generate
!                seeds spaced by spacing_factor times the Hill's radius between the Fluid Roche Limit (FRL) and the outer bin of the
!                ring
!
!  Input
!    Arguments : 
!                
!    Teringinal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : CALL ringmoons_seed_construct(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_seed_construct(swifter_pl1P,ring,seeds)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_construct
      implicit none

! Arguments
      type(swifter_pl),pointer :: swifter_pl1P
      type(ringmoons_ring), intent(inout) :: ring
      type(ringmoons_seeds), intent(inout) :: seeds

! Internals
      integer(I4B)                        :: i,j,seed_bin,inner_outer_sign,nbin,Nactive
      real(DP)                            :: a, dGm, Gmleft
      logical(LGT)                        :: open_space,destructo,spawnbin
      real(DP), parameter                 :: dzone_width = 0.01_DP ! Width of the destruction zone as a fraction of the RRL distance
      integer(I4B)                        :: dzone_inner,dzone_outer ! inner and outer destruction zone bins
      real(DP)                            :: ndz,deltaL,Lorig,c,b,dr,Lring_orig
      real(DP)                            :: Gm_min,R_min
      

! Executable code
  
      
      ! First convert any recently destroyed satellites into ring material
      destructo = .false.
      do i = 1,seeds%N
         if ((.not.seeds%active(i)).and.(seeds%Gm(i) > 0.0_DP)) then
            write(*,*) 'Destruction activated!',i,seeds%a(i),seeds%Gm(i)
            destructo = .true.
            Lring_orig = sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
            Gmleft = seeds%Gm(i)
            Lorig = Gmleft * sqrt((swifter_pl1P%mass + Gmleft) * seeds%a(i)) 
            deltaL = 0.0_DP
            c = dzone_width * ring%RRL ! Create an approximately Gaussian distribution of mass
            a = Gmleft / (sqrt(2 * PI) * c)
            do j = 0,(ring%N - ring%iRRL)
               do inner_outer_sign = -1,1,2
                  nbin = ring%iRRL + inner_outer_sign * j
                  if ((nbin > 0).and.(nbin < ring%N).and.(Gmleft > 0.0_DP)) then
                     dr = 0.5_DP * ring%X(nbin) * ring%deltaX
                     dGm = min(Gmleft,a * dr *exp(-(ring%r(nbin) - ring%RRL)**2 / (2 * c**2)))
                     ring%Gm(nbin) = ring%Gm(nbin) + dGm
                     Gmleft = Gmleft - dGm
                     deltaL = deltaL + (dGm * ring%Iz(nbin) * ring%w(nbin))
                  end if
                  if (j == 0) exit
               end do
               if (Gmleft == 0.0_DP) exit
            end do 
            j = ring%iRRL
            deltaL = Lorig - deltaL
            do nbin = 1,2
               dGm = deltaL / (ring%Iz(j) * ring%w(j) - ring%Iz(j + 1) * ring%w(j + 1))
               ring%Gm(j) = ring%Gm(j) + dGm
               ring%Gm(j + 1) = ring%Gm(j + 1) - dGm
               deltaL = Lring_orig + Lorig - sum(ring%Gm(:) * ring%Iz(:) * ring%w(:)) 
            end do
            ring%Gsigma(:) = ring%Gm(:) / ring%deltaA(:)
            seeds%Gm(i) = 0.0_DP
         end if
      end do
      !if (destructo) then 
         Nactive = count(seeds%active(:))
         seeds%a(1:Nactive) = pack(seeds%a(:),seeds%active(:))
         seeds%Gm(1:Nactive) = pack(seeds%Gm(:),seeds%active(:))
         seeds%active(1:Nactive) = .true.
         if (size(seeds%active) > Nactive) seeds%active(Nactive+1:size(seeds%active)) = .false.
         seeds%N = Nactive
         call ringmoons_update_seeds(swifter_pl1P,ring,seeds)
         call ringmoons_update_ring(swifter_pl1P,ring)
      !end if
      
      ! Make seeds small enough to fit into each bin 
      do i = ring%iRRL,ring%N
         if (ring%Gm(i) < N_DISK_FACTOR * ring%Gm_pdisk(i)) cycle ! don't consider bins that don't have enough mass
         if (ring%Q(i) > 1._DP) cycle  ! don't consider bins that are gravitationally stable
         if (any(seeds%rbin(:) == i .and. seeds%active(:))) cycle ! don't consider bins that already have a seed

         spawnbin = (i >= ring%iFRL)  ! Always spawn seeds at the FRL 

         ! See Tajeddine et al. (2017) section 2.3. Spawn seed if aggregates are 1% the gap opening mass
         R_min = 0.01_DP * (3.3e5 / DU2CM) *  (ring%nu(i) / (100 * TU2S / DU2CM**2))
         if (ring%r_pdisk(i) > R_min) spawnbin = .true.

         if (spawnbin) then
            a = ring%r(i)
            dGm = ring%Gm_pdisk(i)
            call ringmoons_seed_spawn(swifter_pl1P,ring,seeds,a,dGm)
         end if
      end do     
      !read(*,*)
       
      call ringmoons_update_seeds(swifter_pl1P,ring,seeds)
          

end subroutine ringmoons_seed_construct
