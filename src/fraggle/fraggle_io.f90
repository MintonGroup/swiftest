submodule(fraggle_classes) s_fraggle_io
   use swiftest

contains

   module subroutine fraggle_io_log_regime(param, colliders, frag)
      !! author: David A. Minton
      !!
      !! Writes a log of the results of the collisional regime determination
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(in) :: param
      class(fraggle_colliders),   intent(in) :: colliders
      class(fraggle_fragments),   intent(in) :: frag
      ! Internals
      character(STRMAX) :: errmsg

      open(unit=FRAGGLE_LOG_UNIT, file=FRAGGLE_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(FRAGGLE_LOG_UNIT, *) 
      write(FRAGGLE_LOG_UNIT, *) "--------------------------------------------------------------------"
      write(FRAGGLE_LOG_UNIT, *) "           Fraggle collisional regime determination results"
      write(FRAGGLE_LOG_UNIT, *) "--------------------------------------------------------------------"
      write(FRAGGLE_LOG_UNIT, *) "----------------------- Collider information -----------------------" 
      write(FRAGGLE_LOG_UNIT, *) "True number of colliders : ",colliders%ncoll
      write(FRAGGLE_LOG_UNIT, *) "Index list of true colliders  : ",colliders%idx(1:colliders%ncoll)
      write(FRAGGLE_LOG_UNIT, *) "-------------------- Two-body equialent values ---------------------"
      write(FRAGGLE_LOG_UNIT, *) "mass1    : ",colliders%mass(1)
      write(FRAGGLE_LOG_UNIT, *) "radius1  : ",colliders%radius(1)
      write(FRAGGLE_LOG_UNIT, *) "xb1      : ",colliders%xb(:,1)
      write(FRAGGLE_LOG_UNIT, *) "vb1      : ",colliders%vb(:,1)
      write(FRAGGLE_LOG_UNIT, *) "rot1     : ",colliders%rot(:,1)
      write(FRAGGLE_LOG_UNIT, *) "Ip1      : ",colliders%Ip(:,1)
      write(FRAGGLE_LOG_UNIT, *) "L_spin1  : ",colliders%L_spin(:,1)
      write(FRAGGLE_LOG_UNIT, *) "L_orbit1 : ",colliders%L_orbit(:,1)
      write(FRAGGLE_LOG_UNIT, *) "mass2    : ",colliders%mass(2)
      write(FRAGGLE_LOG_UNIT, *) "radius2  : ",colliders%radius(2)
      write(FRAGGLE_LOG_UNIT, *) "xb2      : ",colliders%xb(:,2)
      write(FRAGGLE_LOG_UNIT, *) "vb2      : ",colliders%vb(:,2)
      write(FRAGGLE_LOG_UNIT, *) "rot2     : ",colliders%rot(:,2)
      write(FRAGGLE_LOG_UNIT, *) "Ip2      : ",colliders%Ip(:,2)
      write(FRAGGLE_LOG_UNIT, *) "L_spin2  : ",colliders%L_spin(:,2)
      write(FRAGGLE_LOG_UNIT, *) "L_orbit2 : ",colliders%L_orbit(:,2)
      write(FRAGGLE_LOG_UNIT, *) "------------------------------ Regime -----------------------------"
      select case(frag%regime) 
      case(COLLRESOLVE_REGIME_MERGE)
         write(FRAGGLE_LOG_UNIT, *) "Merge"
      case(COLLRESOLVE_REGIME_DISRUPTION)
         write(FRAGGLE_LOG_UNIT, *) "Disruption"
      case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
         write(FRAGGLE_LOG_UNIT, *) "Supercatastrophic disruption"
      case(COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
         write(FRAGGLE_LOG_UNIT, *) "Graze and merge"
      case(COLLRESOLVE_REGIME_HIT_AND_RUN)
         write(FRAGGLE_LOG_UNIT, *) "Hit and run"
      end select
      write(FRAGGLE_LOG_UNIT, *) "----------------------- Fragment information ----------------------"
      write(FRAGGLE_LOG_UNIT, *) "Total mass of fragments      : ", frag%mtot
      write(FRAGGLE_LOG_UNIT, *) "Largest fragment mass        : ", frag%mass_dist(1)
      write(FRAGGLE_LOG_UNIT, *) "Second-largest fragment mass : ", frag%mass_dist(2)
      write(FRAGGLE_LOG_UNIT, *) "Remaining fragment mass      : ", frag%mass_dist(3)
      write(FRAGGLE_LOG_UNIT, *) "Remaining fragment mass      : ", frag%mass_dist(3)
      write(FRAGGLE_LOG_UNIT, *) "Center of mass position      : ", frag%xbcom(:)
      write(FRAGGLE_LOG_UNIT, *) "Center of mass velocity      : ", frag%vbcom(:)
      write(FRAGGLE_LOG_UNIT, *) "Energy loss                  : ",frag%Qloss
      write(FRAGGLE_LOG_UNIT, *) "--------------------------------------------------------------------"

      return
      667 continue
      write(*,*) "Error writing Fraggle regime information to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_regime


   module subroutine fraggle_io_log_start(param)
      !! author: David A. Minton
      !!
      !! Checks to see if the Fraggle log file needs to be replaced if this is a new run, or appended if this is a restarted run
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(in) :: param
      ! Internals
      character(STRMAX) :: errmsg
      logical           :: fileExists

      inquire(file=FRAGGLE_LOG_OUT, exist=fileExists)
      if (.not.param%lrestart .or. .not.fileExists) then
         open(unit=FRAGGLE_LOG_UNIT, file=FRAGGLE_LOG_OUT, status="REPLACE", err = 667, iomsg = errmsg)
         write(FRAGGLE_LOG_UNIT, *, err = 667, iomsg = errmsg) "Fraggle logfile"
      end if

      return

      667 continue
      write(*,*) "Error writing Fraggle log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_start

end submodule s_fraggle_io