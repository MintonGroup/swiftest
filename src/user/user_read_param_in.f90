submodule (user) s_user_read_param_in
contains
   module procedure user_read_param_in
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Read in parameters for the integration
   !!
   !! Adapted from David E. Kaufmann's Swifter routine io_init_param.f90
   !! Adapted from Martin Duncan's Swift routine io_init_param.f
   !$ use omp_lib
   !use util, only: util_exit ! IMPLEMENTATION TBD
   implicit none

   integer(I4B), parameter :: LUN = 7                 !! Unit number of input file
   integer(I4B)            :: ierr = 0                !! Input error code
   character(STRMAX)       :: error_message           !! Error message in UDIO procedure

   ! Read in name of parameter file
   write(*, *) 'Parameter data file is ', trim(adjustl(inparfile))
   write(*, *) ' '
   100 format(A)
   open(unit = LUN, file = inparfile, status = 'old', iostat = ierr)
   if (ierr /= 0) then
      write(*,*) 'Swiftest error: ', ierr
      write(*,*) '   Unable to open file ',trim(adjustl(inparfile))
      call util_exit(FAILURE)
   end if

   read(LUN,'(DT)', iostat= ierr, iomsg = error_message) param
   if (ierr /= 0) then
      write(*,*) 'Swiftest error reading ', trim(adjustl(inparfile))
      write(*,*) ierr,trim(adjustl(error_message))
      call util_exit(FAILURE)
   end if

   close(LUN)

   write(*,*) "NPLMAX         = ",param%nplmax
   write(*,*) "NTPMAX         = ",param%ntpmax
   write(*,*) "T0             = ",param%t0
   write(*,*) "TSTOP          = ",param%tstop
   write(*,*) "DT             = ",param%dt
   write(*,*) "PL_IN          = ",trim(adjustl(param%inplfile))
   write(*,*) "TP_IN          = ",trim(adjustl(param%intpfile))
   write(*,*) "IN_TYPE        = ",trim(adjustl(param%in_type))
   write(*,*) "ISTEP_OUT      = ",param%istep_out
   write(*,*) "BIN_OUT        = ",trim(adjustl(param%outfile))
   write(*,*) "OUT_TYPE       = ",trim(adjustl(param%out_type))
   write(*,*) "OUT_FORM       = ",trim(adjustl(param%out_form))
   write(*,*) "OUT_STAT       = ",trim(adjustl(param%out_stat))
   write(*,*) "ISTEP_DUMP     = ",param%istep_dump
   write(*,*) "J2             = ",param%j2rp2
   write(*,*) "J4             = ",param%j4rp4
   write(*,*) "CHK_CLOSE      = ",param%feature%lclose
   write(*,*) "CHK_RMIN       = ",param%rmin
   write(*,*) "CHK_RMAX       = ",param%rmax
   write(*,*) "CHK_EJECT      = ",param%rmaxu
   write(*,*) "CHK_QMIN       = ",param%qmin
   write(*,*) "CHK_QMIN_COORD = ",trim(adjustl(param%qmin_coord))
   write(*,*) "CHK_QMIN_RANGE = ",param%qmin_alo, param%qmin_ahi
   write(*,*) "ENC_OUT        = ",trim(adjustl(param%encounter_file))
   write(*,*) "EXTRA_FORCE    = ",param%feature%lextra_force
   write(*,*) "BIG_DISCARD    = ",param%feature%lbig_discard
   write(*,*) "RHILL_PRESENT  = ",param%feature%lrhill_present
   ierr = 0

   ! Added by D. Minton
   MU2KG = param%MU2KG
   TU2S  = param%TU2S 
   DU2M = param%DU2M
   ! The fragmentation model requires the user to set the unit system explicitly.
   write(*,*) "FRAGMENTATION  = ",param%feature%lfragmentation
   if (param%feature%lfragmentation) then
      write(*,*) "MU2KG          = ",MU2KG
      write(*,*) "TU2S           = ",TU2S 
      write(*,*) "DU2M          = ",DU2M
      if ((MU2KG < 0.0_DP) .or. (TU2S < 0.0_DP) .or. (DU2M < 0.0_DP)) then
         write(*,*) 'Invalid unit conversion factor'
         write(*,*) 'MU2KG: ',MU2KG
         write(*,*) 'TU2S: ',TU2S
         write(*,*) 'DU2M: ',DU2M
         ierr = -1
      end if
   end if 
   !Added mtiny to the argument list rather than from the terminal
   if (param%mtiny < 0.0_DP) then
      write(*,*) "Invalid MTINY: ",param%mtiny
      ierr = -1
   else
      write(*,*) "MTINY          = ",param%mtiny   
   end if
   if (param%feature%lpython) write(*,*) "PYTHON         = ",param%feature%lpython
   if (param%feature%lenergy) write(*,*) "ENERGY         = ",param%feature%lenergy
   if (param%feature%lringmoons) write(*,*) "RINGMOONS      = ",param%feature%lringmoons

   if (ierr < 0) then
      write(*, 100) "Input parameter(s) failed check"
      call util_exit(FAILURE)
   end if

   !> Define the maximum number of threads
   nthreads = 1            ! In the *serial* case
   !$ nthreads = omp_get_max_threads() ! In the *parallel* case
   !$ write(*,'(a)')   ' OpenMP parameters:'
   !$ write(*,'(a)')   ' ------------------'
   !$ write(*,'(a,i3,/)') ' Number of threads  = ', nthreads   

   return 

   end procedure user_read_param_in

end submodule s_user_read_param_in
