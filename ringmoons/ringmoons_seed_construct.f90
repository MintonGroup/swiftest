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
      integer(I4B)                        :: i,j,seed_bin,inner_outer_sign,nbin
      real(DP)                            :: a, dGm, Gmleft
      logical(LGT)                        :: open_space
      real(DP), parameter                 :: dzone_width = 0.01_DP ! Width of the destruction zone as a fraction of the RRL distance
      integer(I4B)                        :: dzone_inner,dzone_outer ! inner and outer destruction zone bins
      real(DP)                            :: ndz,deltaL,Lorig,c,b,dr,Lring_orig
      

! Executable code
  
      
      ! First convert any recently destroyed satellites into ring material
      do i = 1,seeds%N
         if ((.not.seeds%active(i)).and.(seeds%Gm(i) > 0.0_DP)) then
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
                     dr = ring%router(nbin) - ring%rinner(nbin)
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
            call ringmoons_viscosity(ring)
            seeds%Gm(i) = 0.0_DP
         end if
      end do

      seeds%N = count(seeds%active(:))
      seeds%a(:) = pack(seeds%a(:),seeds%active(:))
      seeds%Gm(:) = pack(seeds%Gm(:),seeds%active(:))
      seeds%Rhill(:) = pack(seeds%Rhill(:),seeds%active(:))
      seeds%rbin(:) = pack(seeds%rbin(:),seeds%active(:))
      seeds%fz_bin_inner(:) = pack(seeds%fz_bin_inner(:),seeds%active(:))
      seeds%fz_bin_outer(:) = pack(seeds%fz_bin_outer(:),seeds%active(:))
      seeds%Torque(:) = pack(seeds%Torque(:),seeds%active(:))
      seeds%Ttide(:) = pack(seeds%Ttide(:),seeds%active(:))
      seeds%active(1:seeds%N) = .true. 
      


      ! Make seeds small enough to fit into each bin 
      do i = ring%iFrl,ring%N
         if (ring%Gm(i) > INITIAL_MASS_FACTOR * ring%Gm_pdisk) then 
            open_space = .true.
            do j = 1, seeds%N
               if ((i >= seeds%fz_bin_inner(j)) .and. (i <= seeds%fz_bin_outer(j))) then
                  open_space = .false. ! There is already a seed with a feeding zone here
                  exit
               end if
            end do
            if (open_space) then
               a = ring%r(i)
               dGm = INITIAL_MASS_FACTOR * ring%Gm_pdisk !3 * swifter_pl1P%mass * ((ring%router(i) - ring%rinner(i)) / (1.25_DP * FEEDING_ZONE_FACTOR * a))**3
               call ringmoons_seed_spawn(swifter_pl1P,ring,seeds,a,dGm)
            end if
         end if
      end do      

          

end subroutine ringmoons_seed_construct
