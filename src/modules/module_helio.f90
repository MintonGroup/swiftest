!**********************************************************************************************************************************
!
!  Unit Name   : helio
!  Unit Type   : module
!  Project     : SWIFTEST
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the Democratic Heliocentric Method
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Invocation  : N/A
!
!  Notes       : 
!
!**********************************************************************************************************************************
module helio

   use swiftest_globals
   use nbody_data_structures
   implicit none

   type helio_pl
       real(DP),     dimension(:,:),   allocatable :: ah     ! total heliocentric acceleration
       real(DP),     dimension(:,:),   allocatable :: ahi    ! heliocentric acceleration due to interactions
       type(swiftest_pl)                           :: swiftest  ! swifter massive body structure
    end type helio_pl

    type helio_tp
       real(DP),     dimension(:,:),   allocatable :: ah       ! total heliocentric acceleration
       real(DP),     dimension(:,:),   allocatable :: ahi      ! heliocentric acceleration due to interactions
       type(swiftest_tp)                           :: swiftest  ! swifter test particle structure
    end type helio_tp
end module helio
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
