module io
   !! Module containing all input/output subroutine interface blocks 
   use module_parameters
   implicit none

   interface
      module subroutine io_read_param_in(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type,&
            istep_out, outfile, &
            out_type, out_form, out_stat, istep_dump, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,&
            qmin_ahi, encounter_file, mtiny, feature, ring_outfile)
         implicit none
         integer(I4B), intent(out) :: nplmax, ntpmax, istep_out, istep_dump
         real(DP), intent(out)     :: t0, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
         character(*), intent(in)  :: inparfile
         character(*), intent(out) :: qmin_coord, encounter_file, inplfile, intpfile, in_type, outfile, out_type, out_form,out_stat
         real(DP), intent(out), optional :: mtiny 
         type(feature_list), intent(out) :: feature
         character(*), intent(out), optional :: ring_outfile
      end subroutine io_read_param_in
   end interface

end module io


