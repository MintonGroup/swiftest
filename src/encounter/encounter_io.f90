submodule (encounter_classes) s_encounter_io
   use swiftest
contains


   module subroutine encounter_io_write_frame(iu, t, id1, id2, Gmass1, Gmass2, radius1, radius2, xh1, xh2, vh1, vh2)
      !! author: David A. Minton
      !!
      !! Write a single frame of close encounter data to output binary files
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_write_encounter.f90
      !! Adapted from Hal Levison's Swift routine io_write_encounter.f
      implicit none
      ! Arguments
      integer(I4B),           intent(in) :: iu               !! Open file unit number
      real(DP),               intent(in) :: t                !! Time of encounter
      integer(I4B),           intent(in) :: id1, id2         !! ids of the two encountering bodies
      real(DP),               intent(in) :: Gmass1, Gmass2   !! G*mass of the two encountering bodies
      real(DP),               intent(in) :: radius1, radius2 !! Radii of the two encountering bodies
      real(DP), dimension(:), intent(in) :: xh1, xh2         !! Heliocentric position vectors of the two encountering bodies 
      real(DP), dimension(:), intent(in) :: vh1, vh2         !! Heliocentric velocity vectors of the two encountering bodies  
      ! Internals
      character(len=STRMAX)   :: errmsg

      write(iu, err=667, iomsg=errmsg) t
      write(iu, err=667, iomsg=errmsg) id1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), Gmass1, radius1
      write(iu, err=667, iomsg=errmsg) id2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), Gmass2, radius2

      return
      667 continue
      write(*,*) "Error writing encounter file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine

   module subroutine encounter_io_write_list(self, pl, encbody, param)
      implicit none
      ! Arguments
      class(encounter_list),  intent(in) :: self    !! Swiftest encounter list object
      class(swiftest_pl),         intent(in) :: pl      !! Swiftest massive body object
      class(swiftest_body),       intent(in) :: encbody !! Encountering body - Swiftest generic body object (pl or tp) 
      class(swiftest_parameters), intent(in) :: param   !! Current run configuration parameters 
      ! Internals
      logical , save          :: lfirst = .true.
      integer(I4B)            :: k, ierr
      character(len=STRMAX)   :: errmsg

      if (param%enc_out == "" .or. self%nenc == 0) return

      open(unit=LUN, file=param%enc_out, status='OLD', position='APPEND', form='UNFORMATTED', iostat=ierr, iomsg=errmsg)
      if (ierr /= 0) then
         if (lfirst) then
            open(unit=LUN, file=param%enc_out, status='NEW', form='UNFORMATTED', err=667, iomsg=errmsg)
         else
            goto 667
         end if
      end if
      lfirst = .false.

      associate(ind1 => self%index1, ind2 => self%index2)
         select type(encbody)
         class is (swiftest_pl)
            do k = 1, self%nenc
               call encounter_io_write_frame(LUN, self%t(k), &
                                             pl%id(ind1(k)),     encbody%id(ind2(k)), &
                                             pl%Gmass(ind1(k)),  encbody%Gmass(ind2(k)), &
                                             pl%radius(ind1(k)), encbody%radius(ind2(k)), &
                                             self%x1(:,k),       self%x2(:,k), &
                                             self%v1(:,k),       self%v2(:,k))
            end do
         class is (swiftest_tp)
            do k = 1, self%nenc
               call encounter_io_write_frame(LUN, self%t(k), &
                                             pl%id(ind1(k)),     encbody%id(ind2(k)), &
                                             pl%Gmass(ind1(k)),  0.0_DP, &
                                             pl%radius(ind1(k)), 0.0_DP, &
                                             self%x1(:,k),       self%x2(:,k), &
                                             self%v1(:,k),       self%v2(:,k))
            end do 
         end select
      end associate

      close(unit = LUN, err = 667, iomsg = errmsg)

      return
      667 continue
      write(*,*) "Error writing encounter file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine encounter_io_write_list

end submodule s_encounter_io