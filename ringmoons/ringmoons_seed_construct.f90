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
      integer(I4B)                        :: i,j,seed_bin,inner_outer_sign,nbin,Nactive,rbin
      real(DP)                            :: a, dGm, Gmleft
      logical(LGT)                        :: open_space,destructo,spawnbin
      real(DP), parameter                 :: dzone_width = 0.025_DP ! Width of the destruction zone as a fraction of the RRL distance
      integer(I4B)                        :: dzone_inner,dzone_outer ! inner and outer destruction zone bins
      real(DP)                            :: ndz,Lseed_orig,c,b,dr,Lring,deltaL
      real(DP)                            :: Gm_min,R_min,dt,dtleft
      logical(lgt)                        :: stepfail
      type(ringmoons_ring)                :: seedring
      type(ringmoons_seeds)               :: tmpseeds
      

! Executable code
  
      
      ! First convert any recently destroyed satellites into ring material
      destructo = .false.

      do i = 1, seeds%N
         if (seeds%active(i)) then
            if (seeds%a(i) <= ring%RRL) then   ! Destroy the satellite!
               write(*,*) 'We are on our way to destruction!'
               DESTRUCTION_EVENT = .true.
               DESTRUCTION_COUNTER = 0
               seeds%active(i) = .false.
            end if
         end if
      end do

      do i = 1,seeds%N
         if ((.not.seeds%active(i)).and.(seeds%Gm(i) > 0.0_DP)) then
            write(*,*) 'Destruction activated!',i,seeds%a(i),seeds%Gm(i)
            destructo = .true.
            Lseed_orig = seeds%Gm(i) * sqrt((swifter_pl1P%mass + seeds%Gm(i)) * seeds%a(i)) 

            seedring%N = ring%N
            tmpseeds%N = 0
            call ringmoons_allocate(seedring,tmpseeds)
            seedring = ring
            seedring%Gm(:) = 0.0_DP

            Gmleft = seeds%Gm(i)
            c = dzone_width * seeds%a(i) ! Create an approximately Gaussian distribution of mass
            rbin = ringmoons_ring_bin_finder(ring,seeds%a(i))
            a = Gmleft / (sqrt(2 * PI) * c)
            do j = 0,(ring%N - rbin)
               do inner_outer_sign = -1,1,2
                  nbin = rbin + inner_outer_sign * j
                  if ((nbin > seedring%inside).and.(nbin < seedring%N).and.(Gmleft > 0.0_DP)) then
                     dr = 0.5_DP * seedring%X(nbin) * seedring%deltaX
                     dGm = min(Gmleft,a * dr * exp(-(seedring%r(nbin) - seeds%a(i))**2 / (2 * c**2)))
                     seedring%Gm(nbin) = seedring%Gm(nbin) + dGm
                     Gmleft = Gmleft - dGm
                  end if
                  if (j == 0) exit
               end do
               if (Gmleft == 0.0_DP) exit
            end do 
            !j = seedring%iRRL
            ! Offset in angular momentum
            Lring = sum(seedring%Gm(:) * seedring%Iz(:) * seedring%w(:))
            deltaL = Lseed_orig - Lring

            ! Apply a torque to the temporary ring to bring it back to the seed's original angular momentum
            dt = 1._DP
            seedring%nu(1:ring%N) =  1._DP / (16 * 12 * dt / (seedring%deltaX)**2) / seedring%X2(1:ring%N)
            seedring%Gsigma(:) = seedring%Gm(:) / seedring%deltaA(:)
            where (seedring%Gm(:) > 0.0_DP)
               seedring%Torque(:) = deltaL * seedring%Gm(:) / sum(seedring%Gm(:))
            elsewhere
               seedring%Torque(:) = 0.0_DP
            end where
            dtleft = dt
            dt =  ringmoons_ring_timestep(swifter_pl1P,seedring,dt) 
            do 
               call ringmoons_sigma_solver(seedring,swifter_pl1P%mass,dt,stepfail)
               dtleft = dtleft - dt
               if (dtleft <= 0.0_DP) exit
               dt = min(dtleft,dt)
            end do
      
            ring%Gm(:) = ring%Gm(:) + seedring%Gm(:) 
            ring%Gsigma(:) = ring%Gm(:) / ring%deltaA(:)
            seeds%Gm(i) = 0.0_DP
            call ringmoons_deallocate(seedring,tmpseeds)
         end if
      end do
      Nactive = count(seeds%active(:))
      seeds%a(1:Nactive) = pack(seeds%a(:),seeds%active(:))
      seeds%Gm(1:Nactive) = pack(seeds%Gm(:),seeds%active(:))
      seeds%active(1:Nactive) = .true.
      if (size(seeds%active) > Nactive) seeds%active(Nactive+1:size(seeds%active)) = .false.
      seeds%N = Nactive
      seeds%rbin(1:seeds%N) = ringmoons_ring_bin_finder(ring,seeds%a(1:seeds%N))
      
      ! Make seeds small enough to fit into each bin 
      do i = ring%iFRL,ring%N
         spawnbin = .true.
         if (ring%Gsigma(i)*MU2GM/DU2CM**2/GU < 1e-2_DP) spawnbin = .false. ! don't consider bins that don't have enough mass
         if (any(seeds%rbin(:) == i .and. seeds%active(:))) spawnbin = .false. ! don't consider bins that already have a seed
         !! See Tajeddine et al. (2017) section 2.3. Spawn seed if aggregates are 1% the gap opening mass
         !R_min = 0.01_DP * (3.3e5 / DU2CM) *  (ring%nu(i) / (100 * TU2S / DU2CM**2))
         !if (ring%r_pdisk(i) > R_min) spawnbin = .true.
         if (spawnbin) then
            a = ring%r(i)
            dGm = ring%Gm_pdisk(i) * 1
            call ringmoons_seed_spawn(swifter_pl1P,ring,seeds,a,dGm)
         end if
      end do     
      where (seeds%active(:))
         seeds%rbin(:) = ringmoons_ring_bin_finder(ring,seeds%a(:))
      elsewhere
         seeds%rbin(:) = 0
      end where

          

end subroutine ringmoons_seed_construct
