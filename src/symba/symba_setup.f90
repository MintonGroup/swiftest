submodule(symba_classes) s_symba_setup
   use swiftest
contains
   module subroutine symba_setup_pl(self,n)
      !! author: David A. Minton
      !!
      !! Allocate SyMBA test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine symba_setup.f90
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self !! SyMBA test particle object
      integer(I4B),    intent(in)    :: n    !! Number of massive bodies to allocate
      ! Internals
      integer(I4B)                   :: i,j

      !> Call allocation method for parent class
      associate(pl => self)
         call helio_setup_pl(pl, n) 
         if (n <= 0) return
      end associate
      return
   end subroutine symba_setup_pl 

   module subroutine symba_setup_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize an SyMBA nbody system from files and sets up the planetocentric structures.
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self    !! SyMBA system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i, j

      ! Call parent method
      call helio_setup_system(self, param)

      ! Set up the pl-tp planetocentric encounter structures for pl and cb. The planetocentric tp structures are 
      ! generated as necessary during close encounter steps.
      select type(pl => self%pl)
      class is(symba_pl)
         select type(cb => self%cb)
         class is (symba_cb)
            select type (tp => self%tp)
            class is (symba_tp)


            end select
         end select
      end select
   
   end subroutine symba_setup_system

   module subroutine symba_setup_tp(self,n)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(symba_tp), intent(inout) :: self !! SyMBA test particle object
      integer,         intent(in)    :: n    !! Number of test particles to allocate

      !> Call allocation method for parent class
      call helio_setup_tp(self, n) 
      if (n <= 0) return
      return
   end subroutine symba_setup_tp

end submodule s_symba_setup
