module io
   !! author: David A. Minton
   !! todo: Replace XDR with HDF5 
   !!
   !! Module containing all input/output subroutine interface blocks 
   use module_globals
   use module_swiftest
   use user

   interface

      module subroutine io_getn(param,swiftest_plA,swiftest_tpA)
         type(user_input_parameters),intent(inout) :: param      !! Input collection of user-defined parameters
         type(swiftest_pl), intent(inout)  :: swiftest_plA  !! Swiftest data structure to store number of massive bodies
         type(swiftest_tp), intent(inout)  :: swiftest_tpA  !! Swiftest data structure to store number of test particles
      end subroutine io_getn

      module subroutine io_read_pl_in(param,swiftest_plA) 
         type(user_input_parameters),intent(in) :: param         !! Input collection of user-defined parameters
         type(swiftest_pl), intent(inout)  :: swiftest_plA  !! Swiftest data structure to store massive body initial conditions
      end subroutine io_read_pl_in

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
         
   end interface

end module io


