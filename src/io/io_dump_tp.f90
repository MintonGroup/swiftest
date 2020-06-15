!
!  Unit Name   : io_dump_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Dump test particle data to file
!
!  Input
!    Arguments : ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : to dump file
!                ntp          : number of active test particles
!                id           : test particle identifier        (from test particle structure, for each test particle)
!                xh           : heliocentric position           (from test particle structure, for each test particle)
!                vh           : heliocentric velocity           (from test particle structure, for each test particle)
!
!  Invocation  : CALL io_dump_tp(ntp, swifter_tp1P)
!
!  Notes       : Adapted from Martin Duncan's Swift routine io_dump_tp.f
!
!**********************************************************************************************************************************
subroutine io_dump_tp(ntp, swiftest_tpa)
   use swiftest, except_this_one => io_dump_tp
   implicit none

   integer(I4B), intent(in)            :: ntp
   type(swiftest_tp), intent(inout)    :: swiftest_tpa

   integer(I4B)                        :: i, iu, ierr
   integer(I4B), save                  :: idx = 1
   integer(I4B), parameter             :: lun = 7

   open(unit = LUN, file = DUMP_TP_FILE(idx), form = "UNFORMATTED", status = 'REPLACE', iostat = ierr)

   if (ierr /= 0) then
        write(*, *) "swiftest error:"
        write(*, *) "   unable to open binary dump file ", trim(dump_tp_file(idx))
        call util_exit(failure)
   end if
   write(LUN) ntp
   write(LUN) swiftest_tpA%name(:)
   write(LUN) swiftest_tpA%xh(:,:)
   write(LUN) swiftest_tpA%vh(:,:)
   close(LUN)
   idx = idx + 1
   if (idx > 2) idx = 1

   return

END subroutine io_dump_tp
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
