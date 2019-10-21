!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step the fluid ring in time
!
!  Input
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                t              : time
!                npl            : number of planets
!                nplmax         : maximum allowed number of planets
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                j2rp2          : J2 * R**2 for the Sun
!                j4rp4          : J4 * R**4 for the Sun
!                dt             : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                eoffset        : energy offset (net energy lost to ring)
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL ringmoons_step(lfirst, t, npl, nplmax, ntp, ntpmax, symba_pl1P, j2rp2, j4rp4, eoffset, dt)
!
!  Notes       : Adapted from Andy Hesselbrock's RING-MOONS Python scripts
!
!**********************************************************************************************************************************
SUBROUTINE ringmoons_step(lfirst, t, rmin,npl, nplmax, symba_pl1P, j2rp2, j4rp4, eoffset, dt,ring)

! Modules
     USE module_parameters
     USE module_symba
     USE module_ringmoons
     USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_step
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(INOUT)                      :: lfirst
     INTEGER(I4B), INTENT(IN)                         :: npl, nplmax
     REAL(DP), INTENT(IN)                             :: t, rmin,j2rp2, j4rp4, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset
     TYPE(symba_pl), POINTER                          :: symba_pl1P
     TYPE(ringmoons_ring),INTENT(INOUT) :: ring
! Internals
      integer(I4B) :: i
      real(DP),parameter :: s2y = 365.25_DP * 24 * 60 * 60

! Executable code
      if (lfirst) then
         call ringmoons_viscosity(symba_pl1P%helio%swifter%mass,rmin,ring)
         open(unit=22,file='test.ic',status='replace')
         do i = 1,ring%N
            write(22,*) ring%r(i), ring%sigma(i), ring%nu(i)
         end do
         close(22)
      end if
      call ringmoons_pde_solver(symba_pl1P%helio%swifter%mass,rmin,dt,ring)

      if ((t >= 0.999*1e3 * s2y)  .and. (t < 1.001*1e3 * s2y)) then
         open(unit=22,file='test.1e3y',status='replace')
         do i = 1,ring%N
            write(22,*) ring%r(i), ring%sigma(i), ring%nu(i) 
         end do
         close(22)
      else if ((t >= 0.999*1e4 * s2y)  .and. (t < 1.001*1e4 * s2y)) then
         open(unit=22,file='test.1e4y',status='replace')
         do i = 1,ring%N
            write(22,*) ring%r(i), ring%sigma(i), ring%nu(i)
         end do
         close(22)
      else if ((t >= 0.999*1e5 * s2y)  .and. (t < 1.001*1e5 * s2y)) then
         open(unit=22,file='test.1e5y',status='replace')
         do i = 1,ring%N
            write(22,*) ring%r(i), ring%sigma(i), ring%nu(i)
         end do
         close(22)
      end if
       


     RETURN

END SUBROUTINE ringmoons_step
!**********************************************************************************************************************************
!
!  Author(s)   : David A. Minton
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
