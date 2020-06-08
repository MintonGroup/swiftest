module io
   !! author: David A. Minton
   !! todo: Replace XDR with HDF5 
   !!
   !! Module containing all input/output subroutine interface blocks 
   use swiftest
   use module_swiftest

   interface

      module function io_get_token(buffer, ifirst, ilast, ierr) result(token)
         character(len=*), intent(in)     :: buffer         !! Input string buffer
         integer(I4B), intent(inout)      :: ifirst         !! Index of the buffer at which to start the search for a token
         integer(I4B), intent(out)        :: ilast          !! Index of the buffer at the end of the returned token
         integer(I4B), intent(out)        :: ierr           !! Error code
         character(len=:),allocatable     :: token          !! Returned token string
      end function io_get_token

      module subroutine io_getn(param,swiftest_plA,swiftest_tpA)
         type(input_parameters),intent(inout) :: param      !! Input collection of user-defined parameters
         type(swiftest_pl), intent(inout)  :: swiftest_plA  !! Swiftest data structure to store number of massive bodies
         type(swiftest_tp), intent(inout)  :: swiftest_tpA  !! Swiftest data structure to store number of test partifles
      end subroutine io_getn

      module function io_read_param_in(inparfile) result(param)
         character(*), intent(in)         :: inparfile     !! Parameter input file name (typically param.in)
         type(input_parameters)           :: param         !! Output collection of user-defined parameters
      end function io_read_param_in

      module subroutine io_read_pl_in(param,swiftest_plA) 
         type(input_parameters),intent(in) :: param         !! Input collection of user-defined parameters
         type(swiftest_pl), intent(inout)  :: swiftest_plA  !! Swiftest data structure to store massive body initial conditions
      end subroutine io_read_pl_in

      module subroutine io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
         integer(I4B), intent(in) :: iu            !! Output file unit number
         integer(I4B), intent(in) :: npl           !! Number of massive bodies
         integer(I4B), intent(in) :: ntp           !! Number of test particles
         integer(I4B), intent(in) :: iout_form     !! Output format type (EL, XV, FILT - see swiftest module for symbolic name definitions)
         real(DP), intent(in)     :: t             !! Current time of simulation
         character(*), intent(in) :: out_type      !! Output file format type (REAL4, REAL8 - see swiftest module for symbolic name definitions)
      end subroutine io_write_hdr         
         
   end interface

end module io


