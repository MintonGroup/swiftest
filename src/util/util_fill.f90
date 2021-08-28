submodule (swiftest_classes) s_util_fill
   use swiftest
contains

   module subroutine util_fill_arr_char_string(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of type character strings
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      character(len=STRMAX), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,               dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))

      return
   end subroutine util_fill_arr_char_string

   module subroutine util_fill_arr_DP(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of type DP
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      real(DP), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,  dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))

      return
   end subroutine util_fill_arr_DP

   module subroutine util_fill_arr_DPvec(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of DP vectors with shape (NDIM, n)
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      real(DP), dimension(:,:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,  dimension(:),                intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! Internals
      integer(I4B) :: i

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      do i = 1, NDIM
         keeps(i,:) = unpack(keeps(i,:),   .not.lfill_list(:), keeps(i,:))
         keeps(i,:) = unpack(inserts(i,:),      lfill_list(:), keeps(i,:))
      end do

      return
   end subroutine util_fill_arr_DPvec

   module subroutine util_fill_arr_I4B(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of type I4B
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      integer(I4B), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,      dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))

      return
   end subroutine util_fill_arr_I4B


   module subroutine util_fill_arr_info(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of particle origin information types
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      type(swiftest_particle_info), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,                      dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! Internals
      integer(I4B), dimension(:), allocatable  :: insert_idx
      integer(I4B) :: i, nkeep, ninsert

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      nkeep = size(keeps)
      ninsert = count(lfill_list)

      allocate(insert_idx(ninsert))

      insert_idx(:) = pack([(i, i = 1, nkeep)], lfill_list)
      call util_copy_particle_info_arr(inserts, keeps, insert_idx)

      return
   end subroutine util_fill_arr_info


   module subroutine util_fill_arr_logical(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of logicals
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      logical, dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical, dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))

      return
   end subroutine util_fill_arr_logical


   module subroutine util_fill_body(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest generic particle structure into an old one. 
      !! This is the inverse of a spill operation.
      implicit none
      ! Arguments
      class(swiftest_body),  intent(inout) :: self       !! Swiftest generic body object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! internals
      integer(I4B) :: i

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Fill all the common components
      associate(keeps => self)
         call util_fill(keeps%id, inserts%id, lfill_list)
         call util_fill(keeps%info, inserts%info, lfill_list)
         call util_fill(keeps%status, inserts%status, lfill_list)
         call util_fill(keeps%ldiscard, inserts%ldiscard, lfill_list)
         call util_fill(keeps%lmask, inserts%lmask, lfill_list)
         call util_fill(keeps%mu, inserts%mu, lfill_list)
         call util_fill(keeps%xh, inserts%xh, lfill_list)
         call util_fill(keeps%vh, inserts%vh, lfill_list)
         call util_fill(keeps%xb, inserts%xb, lfill_list)
         call util_fill(keeps%vb, inserts%vb, lfill_list)
         call util_fill(keeps%ah, inserts%ah, lfill_list)
         call util_fill(keeps%aobl, inserts%aobl, lfill_list)
         call util_fill(keeps%agr, inserts%agr, lfill_list)
         call util_fill(keeps%atide, inserts%atide, lfill_list)
         call util_fill(keeps%a, inserts%a, lfill_list)
         call util_fill(keeps%e, inserts%e, lfill_list)
         call util_fill(keeps%inc, inserts%inc, lfill_list)
         call util_fill(keeps%capom, inserts%capom, lfill_list)
         call util_fill(keeps%omega, inserts%omega, lfill_list)
         call util_fill(keeps%capm, inserts%capm, lfill_list)
           
         ! This is the base class, so will be the last to be called in the cascade. 
         keeps%nbody = size(keeps%id(:))
      end associate
     
      return
   end subroutine util_fill_body


   module subroutine util_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest massive body structure into an old one. 
      !! This is the inverse of a spill operation.
      implicit none
      ! Arguments
      class(swiftest_pl),    intent(inout) :: self       !! Swiftest massive body object
      class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! Internals
      integer(I4B) :: i

      associate(keeps => self)

      select type (inserts) ! The standard requires us to select the type of both arguments in order to access all the components
         class is (swiftest_pl)
            !> Fill components specific to the massive body class
            call util_fill(keeps%mass, inserts%mass, lfill_list)
            call util_fill(keeps%Gmass, inserts%Gmass, lfill_list)
            call util_fill(keeps%rhill, inserts%rhill, lfill_list)
            call util_fill(keeps%radius, inserts%radius, lfill_list)
            call util_fill(keeps%density, inserts%density, lfill_list)
            call util_fill(keeps%k2, inserts%k2, lfill_list)
            call util_fill(keeps%Q, inserts%Q, lfill_list)
            call util_fill(keeps%tlag, inserts%tlag, lfill_list)
            call util_fill(keeps%xbeg, inserts%xbeg, lfill_list)
            call util_fill(keeps%vbeg, inserts%vbeg, lfill_list)
            call util_fill(keeps%Ip, inserts%Ip, lfill_list)
            call util_fill(keeps%rot, inserts%rot, lfill_list)
            
            call util_fill_body(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on swiftest_pl'
         end select
      end associate

      return
   end subroutine util_fill_pl


   module subroutine util_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new Swiftest test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      implicit none
      ! Arguments
      class(swiftest_tp),    intent(inout) :: self       !! Swiftest test particle object
      class(swiftest_body),  intent(in)    :: inserts    !! Swiftest body object to be inserted
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (swiftest_tp)
            !> Spill components specific to the test particle class
            call util_fill(keeps%isperi, inserts%isperi, lfill_list)
            call util_fill(keeps%peri, inserts%peri, lfill_list)
            call util_fill(keeps%atp, inserts%atp, lfill_list)

            call util_fill_body(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on swiftest_tp'
         end select
      end associate

      return
   end subroutine util_fill_tp

end submodule s_util_fill