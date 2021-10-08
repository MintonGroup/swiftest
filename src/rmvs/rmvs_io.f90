submodule (rmvs_classes) s_rmvs_io
   use swiftest
contains

   module subroutine rmvs_io_write_encounter(t, id1, id2, Gmass1, Gmass2, radius1, radius2, &
                                 xh1, xh2, vh1, vh2, enc_out)
      !! author: David A. Minton
      !!
      !! Write close encounter data from RMVS to output binary files
      !!  There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_write_encounter.f90
      !! Adapted from Hal Levison's Swift routine io_write_encounter.f
      implicit none
      ! Arguments
      integer(I4B),           intent(in) :: id1, id2
      real(DP),               intent(in) :: t, Gmass1, Gmass2, radius1, radius2
      real(DP), dimension(:), intent(in) :: xh1, xh2, vh1, vh2
      character(*),           intent(in) :: enc_out
      ! Internals
      logical , save    :: lfirst = .true.
      integer(I4B)        :: ierr

      if (enc_out == "") return

      open(unit = LUN, file = enc_out, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr)
      if ((ierr /= 0) .and. lfirst) then
         open(unit = LUN, file = enc_out, status = 'NEW', form = 'UNFORMATTED', iostat = ierr)
      end if
      if (ierr /= 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to open binary encounter file"
         call util_exit(FAILURE)
      end if
      lfirst = .false.
      call encounter_io_write_frame(LUN, t, id1, id2, Gmass1, Gmass2, radius1, radius2, xh1, xh2, vh1, vh2)
      close(unit = LUN, iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to close binary encounter file"
         call util_exit(FAILURE)
      end if

      return
   end subroutine rmvs_io_write_encounter

end submodule s_rmvs_io