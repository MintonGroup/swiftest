submodule(fraggle_classes) s_fraggle_io
   use swiftest

contains

   module subroutine fraggle_io_log_generate(frag)
      !! author: David A. Minton
      !!
      !! Writes a log of the results of the fragment generation
      implicit none
      ! Arguments
      class(fraggle_fragments),   intent(in) :: frag
      ! Internals
      integer(I4B) :: i
      character(STRMAX) :: errmsg
      character(len=*), parameter :: fmtlabel = "(A14,10(ES11.4,1X,:))"

      open(unit=LUN, file=FRAGGLE_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *, err = 667, iomsg = errmsg)
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "           Fraggle fragment generation results"
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN,    "(' dL_tot should be very small' )")
      write(LUN,fmtlabel) ' dL_tot      |', (.mag.(frag%Ltot_after(:) - frag%Ltot_before(:))) / (.mag.frag%Ltot_before(:))
      write(LUN,        "(' dE_tot should be negative and equal to Qloss' )")
      write(LUN,fmtlabel) ' dE_tot      |', (frag%Etot_after - frag%Etot_before) / abs(frag%Etot_before)
      write(LUN,fmtlabel) ' Qloss       |', -frag%Qloss / abs(frag%Etot_before)
      write(LUN,fmtlabel) ' dE - Qloss  |', (frag%Etot_after - frag%Etot_before + frag%Qloss) / abs(frag%Etot_before)
      write(LUN,        "(' -------------------------------------------------------------------------------------')")
      write(LUN, *) "Individual fragment values (collisional system natural units)"
      write(LUN, *) "mass"
      do i = 1, frag%nbody
         write(LUN, *) i, frag%mass(i)
      end do
      write(LUN, *) "x_coll"
      do i = 1, frag%nbody
         write(LUN, *) i, frag%x_coll(:,i)
      end do
      write(LUN, *) "v_coll"
      do i = 1, frag%nbody
         write(LUN, *) i, frag%v_coll(:,i)
      end do
      write(LUN, *) "xb"
      do i = 1, frag%nbody
         write(LUN, *) i, frag%xb(:,i)
      end do
      write(LUN, *) "vb"
      do i = 1, frag%nbody
         write(LUN, *) i, frag%vb(:,i)
      end do
      write(LUN, *) "rot"
      do i = 1, frag%nbody
         write(LUN, *) i, frag%rot(:,i)
      end do

      close(LUN)

      return
      667 continue
      write(*,*) "Error writing Fraggle message to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_generate


   ! module subroutine io_log_one_message(FRAGGLE_LOG_OUT, message)
   !    !! author: David A. Minton
   !    !!
   !    !! Writes a single message to the fraggle log file
   !    implicit none
   !    ! Arguments
   !    character(len=*), intent(in) :: message
   !    ! Internals
   !    character(STRMAX) :: errmsg

   !    open(unit=LUN, file=FRAGGLE_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
   !    write(LUN, *) trim(adjustl(message)) 
   !    close(LUN)

   !    return
   !    667 continue
   !    write(*,*) "Error writing Fraggle message to log file: " // trim(adjustl(errmsg))
   ! end subroutine fraggle_io_log_one_message


   module subroutine fraggle_io_log_pl(pl, param)
      !! author: David A. Minton
      !!
      !! Writes a single message to the fraggle log file
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(in) :: pl    !! Swiftest massive body object (only the new bodies generated in a collision)
      class(swiftest_parameters), intent(in) :: param !! Current swiftest run configuration parameters
      ! Internals
      integer(I4B) :: i
      character(STRMAX) :: errmsg

      open(unit=LUN, file=FRAGGLE_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *, err = 667, iomsg = errmsg)

      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "           Fraggle fragment final body properties"
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "id, name"
      do i = 1, pl%nbody
         write(LUN, *) i, pl%id(i), pl%info(i)%name
      end do
      write(LUN, *) "mass, Gmass"
      do i = 1, pl%nbody
         write(LUN, *) i, pl%mass(i), pl%Gmass(i)
      end do
      write(LUN, *) "radius"
      do i = 1, pl%nbody
         write(LUN, *) i, pl%radius(i)
      end do
      write(LUN, *) "xb"
      do i = 1, pl%nbody
         write(LUN, *) i, pl%xb(:,i)
      end do
      write(LUN, *) "vb"
      do i = 1, pl%nbody
         write(LUN, *) i, pl%vb(:,i)
      end do
      write(LUN, *) "xh"
      do i = 1, pl%nbody
         write(LUN, *) i, pl%xh(:,i)
      end do
      write(LUN, *) "vh"
      do i = 1, pl%nbody
         write(LUN, *) i, pl%vh(:,i)
      end do

      if (param%lrotation) then
         write(LUN, *) "rot"
         do i = 1, pl%nbody
            write(LUN, *) i, pl%rot(:,i)
         end do
         write(LUN, *) "Ip"
         do i = 1, pl%nbody
            write(LUN, *) i, pl%Ip(:,i)
         end do
      end if

      if (param%ltides) then
         write(LUN, *) "Q"
         do i = 1, pl%nbody
            write(LUN, *) i, pl%Q(i)
         end do
         write(LUN, *) "k2"
         do i = 1, pl%nbody
            write(LUN, *) i, pl%k2(i)
         end do
         write(LUN, *) "tlag"
         do i = 1, pl%nbody
            write(LUN, *) i, pl%tlag(i)
         end do
      end if

      close(LUN)

      return
      667 continue
      write(*,*) "Error writing Fraggle message to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_pl


   module subroutine fraggle_io_log_regime(colliders, frag)
      !! author: David A. Minton
      !!
      !! Writes a log of the results of the collisional regime determination
      implicit none
      ! Arguments
      class(fraggle_colliders),   intent(in) :: colliders !! Fraggle collider system object
      class(fraggle_fragments),   intent(in) :: frag      !! Fraggle fragment object
      ! Internals
      character(STRMAX) :: errmsg

      open(unit=LUN, file=FRAGGLE_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *, err = 667, iomsg = errmsg)
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "           Fraggle collisional regime determination results"
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "----------------------- Collider information -----------------------" 
      write(LUN, *) "True number of colliders : ",colliders%ncoll
      write(LUN, *) "Index list of true colliders  : ",colliders%idx(1:colliders%ncoll)
      write(LUN, *) "-------------------- Two-body equialent values ---------------------"
      write(LUN, *) "mass1    : ",colliders%mass(1)
      write(LUN, *) "radius1  : ",colliders%radius(1)
      write(LUN, *) "xb1      : ",colliders%xb(:,1)
      write(LUN, *) "vb1      : ",colliders%vb(:,1)
      write(LUN, *) "rot1     : ",colliders%rot(:,1)
      write(LUN, *) "Ip1      : ",colliders%Ip(:,1)
      write(LUN, *) "L_spin1  : ",colliders%L_spin(:,1)
      write(LUN, *) "L_orbit1 : ",colliders%L_orbit(:,1)
      write(LUN, *) "mass2    : ",colliders%mass(2)
      write(LUN, *) "radius2  : ",colliders%radius(2)
      write(LUN, *) "xb2      : ",colliders%xb(:,2)
      write(LUN, *) "vb2      : ",colliders%vb(:,2)
      write(LUN, *) "rot2     : ",colliders%rot(:,2)
      write(LUN, *) "Ip2      : ",colliders%Ip(:,2)
      write(LUN, *) "L_spin2  : ",colliders%L_spin(:,2)
      write(LUN, *) "L_orbit2 : ",colliders%L_orbit(:,2)
      write(LUN, *) "------------------------------ Regime -----------------------------"
      select case(frag%regime) 
      case(COLLRESOLVE_REGIME_MERGE)
         write(LUN, *) "Merge"
      case(COLLRESOLVE_REGIME_DISRUPTION)
         write(LUN, *) "Disruption"
      case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
         write(LUN, *) "Supercatastrophic disruption"
      case(COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
         write(LUN, *) "Graze and merge"
      case(COLLRESOLVE_REGIME_HIT_AND_RUN)
         write(LUN, *) "Hit and run"
      end select
      write(LUN, *) "----------------------- Fragment information ----------------------"
      write(LUN, *) "Total mass of fragments      : ", frag%mtot
      write(LUN, *) "Largest fragment mass        : ", frag%mass_dist(1)
      write(LUN, *) "Second-largest fragment mass : ", frag%mass_dist(2)
      write(LUN, *) "Remaining fragment mass      : ", frag%mass_dist(3)
      write(LUN, *) "Center of mass position      : ", frag%xbcom(:)
      write(LUN, *) "Center of mass velocity      : ", frag%vbcom(:)
      write(LUN, *) "Energy loss                  : ", frag%Qloss
      write(LUN, *) "--------------------------------------------------------------------"
      close(LUN)

      return
      667 continue
      write(*,*) "Error writing Fraggle regime information to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_regime


   ! module subroutine fraggle_io_log_start(param)
   !    !! author: David A. Minton
   !    !!
   !    !! Checks to see if the Fraggle log file needs to be replaced if this is a new run, or appended if this is a restarted run
   !    implicit none
   !    ! Arguments
   !    class(swiftest_parameters), intent(in) :: param
   !    ! Internals
   !    character(STRMAX) :: errmsg
   !    logical           :: fileExists

   !    inquire(file=FRAGGLE_LOG_OUT, exist=fileExists)
   !    if (.not.param%lrestart .or. .not.fileExists) then
   !       open(unit=LUN, file=FRAGGLE_LOG_OUT, status="REPLACE", err = 667, iomsg = errmsg)
   !       write(LUN, *, err = 667, iomsg = errmsg) "Fraggle logfile"
   !    end if
   !    close(LUN)

   !    return

   !    667 continue
   !    write(*,*) "Error writing Fraggle log file: " // trim(adjustl(errmsg))
   ! end subroutine fraggle_io_log_start

end submodule s_fraggle_io