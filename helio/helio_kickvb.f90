!**********************************************************************************************************************************
!
!  Unit Name   : helio_kickvb
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Kick barycentric velocities of planets
!
!  Input
!    Arguments : npl        : number of planets
!                helio_pl1P : pointer to head of helio planet structure linked-list
!                dt         : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_kickvb(npl, helio_pl1P, dt)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_kickvb(npl, helio_plA, dt)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_kickvb
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: dt
     TYPE(helio_pl), DIMENSION(:), INTENT(INOUT) :: helio_plA

! Internals
     INTEGER(I4B)              :: i

! Executable code
     !Removed by D. Minton
     !helio_plP => helio_pl1P
     !^^^^^^^^^^^^^^^^^^^^
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
     !$OMP PRIVATE(i,helio_plP,swifter_plP) &
     !$OMP SHARED(npl,helio_pl1P,dt) 
     DO i = 2, npl
          !Removed by D. Minton
          !helio_plP => helio_plP%nextP
          !^^^^^^^^^^^^^^^^^^^^
          !Added by D. Minton
          !helio_plP => helio_pl1P%helio_plPA(i)%thisP

          helio_plA%swiftest%vb(:,i) = helio_plA%swiftest%vb(:,i) + helio_plA%ah(:,i)*dt
     END DO
     !$OMP END PARALLEL DO

     RETURN

END SUBROUTINE helio_kickvb
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
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
