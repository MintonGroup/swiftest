module io
   !! author: David A. Minton
   !! todo: Replace XDR with HDF5 
   !!
   !! Module containing all input/output subroutine interface blocks 
   use swiftest_globals
   use swiftest_data_structures

   interface

      module function io_get_token(buffer, ifirst, ilast, ierr) result(token)
         character(len=*), intent(in)     :: buffer         !! Input string buffer
         integer(I4B), intent(inout)      :: ifirst         !! Index of the buffer at which to start the search for a token
         integer(I4B), intent(out)        :: ilast          !! Index of the buffer at the end of the returned token
         integer(I4B), intent(out)        :: ierr           !! Error code
         character(len=:),allocatable     :: token          !! Returned token string
      end function io_get_token

      module subroutine io_write_frame(t, swiftest_plA, swiftest_tpA, outfile, out_type, out_form, out_stat)
         real(DP), intent(in)             :: t              !! Current time of simulation
         type(swiftest_pl), intent(inout) :: swiftest_plA   !! Swiftest massive body structure
         type(swiftest_tp), intent(inout) :: swiftest_tpA   !! Swiftest test particle structure
         character(*), intent(in)         :: outfile        !! Name of output file
         character(*), intent(in)         :: out_type       !! Output file format type (REAL4, REAL8 - see swiftest module for 
                                                            !!    symbolic name definitions)
         character(*), intent(in)         :: out_form       !! Output format type (EL, XV,- see swiftest module for symbolic name definitions)
         character(*), intent(in)         :: out_stat       !! Output status code (NEW, APPEND)
      end subroutine io_write_frame

      module subroutine io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
         integer(I4B), intent(in) :: iu            !! Output file unit number
         real(DP), intent(in)     :: t             !! Current time of simulation
         integer(I4B), intent(in) :: npl           !! Number of massive bodies
         integer(I4B), intent(in) :: ntp           !! Number of test particles
         integer(I4B), intent(in) :: iout_form     !! Output format type (EL, XV,- see swiftest module for symbolic name definitions)
         character(*), intent(in) :: out_type      !! Output file format type (REAL4, REAL8 - see swiftest module for symbolic name definitions)
      end subroutine io_write_hdr

      subroutine io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, nmergesub, symba_plA, & 
         discard_plA, discard_tpA, mergeadd_list, mergesub_list, fname, lbig_discard)
         use swiftest_globals
         use swiftest_data_structures
         use symba
         implicit none
         logical(lgt), intent(in)       :: lbig_discard
         integer(I4B), intent(in)       :: npl, ntp, nsppl, nsptp, nmergeadd, nmergesub
         real(DP), intent(in)         :: t, mtiny
         character(*), intent(in)       :: fname
         type(symba_pl), intent(inout)      :: symba_plA
         type(swiftest_tp), intent(inout)      :: discard_tpA
         type(swiftest_pl), intent(inout)      :: discard_plA
         type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
      end subroutine io_discard_write_symba

      subroutine io_dump_pl(npl, swiftest_plA, lclose, lrhill_present)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         logical(lgt), intent(in)   :: lclose, lrhill_present
         integer(I4B), intent(in)   :: npl
         type(swiftest_pl), intent(inout):: swiftest_plA
      end subroutine io_dump_pl

      subroutine io_dump_tp(ntp, swiftest_tpA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)   :: ntp
         type(swiftest_tp), intent(inout):: swiftest_tpA
      end subroutine io_dump_tp

      function io_read_encounter(t, name1, name2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file, out_type)
         use swiftest_globals
         implicit none
         integer(I4B)         :: io_read_encounter
         integer(I4B), intent(out)     :: name1, name2
         real(DP), intent(out)      :: t, mass1, mass2
         real(DP), dimension(ndim), intent(out) :: xh1, xh2, vh1, vh2
         character(*), intent(in)      :: encounter_file,out_type
      end function io_read_encounter

      function io_read_hdr(iu, t, npl, ntp, iout_form, out_type)
         use swiftest_globals
         implicit none
         integer(I4B)      :: io_read_hdr
         integer(I4B), intent(in)   :: iu
         integer(I4B), intent(out)  :: npl, ntp, iout_form
         real(DP), intent(out)   :: t
         character(*), intent(in)   :: out_type
      end function io_read_hdr

      function io_read_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, mass, radius)
         use swiftest_globals
         implicit none
         integer(I4B)      :: io_read_line
         integer(I4B), intent(in)   :: iu
         integer(I4B), intent(out)   :: name
         real(DP), intent(out)    :: d1, d2, d3, d4, d5, d6
         real(DP), optional, intent(out) :: mass, radius
         character(*), intent(in)   :: out_type
      end function io_read_line

      subroutine io_write_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
         xh1, xh2, vh1, vh2, encounter_file, out_type)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in)     :: name1, name2
         real(DP), intent(in)      :: t, mass1, mass2, radius1, radius2
         real(DP), dimension(ndim), intent(in) :: xh1, xh2, vh1, vh2
         character(*), intent(in)     :: encounter_file, out_type
      end subroutine io_write_encounter

      subroutine io_write_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, mass, radius)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in)   :: iu, name
         real(DP), intent(in)    :: d1, d2, d3, d4, d5, d6
         real(DP), optional, intent(in) :: mass, radius
         character(*), intent(in)   :: out_type
      end subroutine io_write_line    


         
   end interface

end module io


