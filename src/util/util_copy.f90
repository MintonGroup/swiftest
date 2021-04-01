submodule (swiftest_classes) s_util_copy
   use swiftest
contains
   module procedure util_copy_cb
      !! author: David A. Minton
      !!
      !! Copies elements of one Swiftest central body object to another. The mask is ignored for this case
      !! 
      implicit none
      select type(src)
      class is(swiftest_cb)
         self%mass = src%mass 
         self%Gmass = src%Gmass
         self%radius = src%radius
         self%density = src%density
         self%j2rp2 = src%j2rp2
         self%j4rp4 = src%j4rp4
         self%aobl = src%aobl
         self%xb = src%xb 
         self%vb = src%vb
         self%Ip = src%Ip 
         self%rot = src%rot
         self%k2 = src%k2
         self%Q = src%Q
      class default
         write(*,*) 'Error copying swiftest_cb object. Incompatible type.'
         call util_exit(FAILURE)
      end select
      return
   end procedure util_copy_cb

   module procedure util_copy_body
      !! author: David A. Minton
      !!
      !! Copies elements of one Swiftest generic body object to another. The elements are determined by a mask.
      !! 
      implicit none
      
      select type(src)
      class is(swiftest_body)
         where (mask(:))
            self%name(:) = src%name(:)
            self%status(:) = src%status(:)
            self%ldiscard(:) = src%ldiscard(:)
            self%xh(1,:) = src%xh(1,:)
            self%xh(2,:) = src%xh(2,:)
            self%xh(3,:) = src%xh(3,:)
            self%vh(1,:) = src%vh(1,:)
            self%vh(2,:) = src%vh(2,:)
            self%vh(3,:) = src%vh(3,:)
            self%xb(1,:) = src%xb(1,:)
            self%xb(2,:) = src%xb(2,:)
            self%xb(3,:) = src%xb(3,:)
            self%vb(1,:) = src%vb(1,:)
            self%vb(2,:) = src%vb(2,:)
            self%vb(3,:) = src%vb(3,:)
            self%ah(1,:) = src%ah(1,:)
            self%ah(2,:) = src%ah(2,:)
            self%ah(3,:) = src%ah(3,:)
            self%aobl(1,:) = src%aobl(1,:)
            self%aobl(3,:) = src%aobl(2,:)
            self%aobl(2,:) = src%aobl(3,:)
            self%ir3h(:) = src%ir3h(:)
            self%a(:) = src%a(:)
            self%e(:) = src%e(:)
            self%inc(:) = src%inc(:)
            self%capom(:) = src%capom(:)
            self%omega(:) = src%omega(:)
            self%capm(:) = src%capm(:)
            self%mu(:) = src%mu(:)
         end where
      class default
         write(*,*) 'Error copying swiftest_body object. Incompatible type.'
         call util_exit(FAILURE)
      end select

   end procedure util_copy_body

   module procedure util_copy_pl
      !! author: David A. Minton
      !!
      !! Copies elements of one Swiftest massive body object to another. The elements are determined by a mask.
      !! 
      implicit none
      select type(src)
      class is(swiftest_pl)
         where (mask(:))
            self%mass(:) = src%mass(:)
            self%Gmass(:) = src%Gmass(:)
            self%rhill(:) = src%rhill(:)
            self%radius(:) = src%radius(:)
            self%density(:) = src%density(:)
            self%Ip(1,:) = src%Ip(1,:)
            self%Ip(2,:) = src%Ip(2,:)
            self%Ip(3,:) = src%Ip(3,:)
            self%rot(1,:) = src%rot(1,:)
            self%rot(2,:) = src%rot(2,:)
            self%rot(3,:) = src%rot(3,:)
            self%k2(:) = src%k2(:)
            self%Q(:) = src%Q(:)
         end where
         call util_copy_body(self, src, mask)
      class default
         write(*,*) 'Error copying swiftest_pl object. Incompatible type.'
         call util_exit(FAILURE)
      end select

   end procedure util_copy_pl

   module procedure util_copy_tp
       !! author: David A. Minton
      !!
      !! Copies elements of one swiftest test particle object to another. The elements are determined by a mask.
      implicit none

      select type(src)
      class is(swiftest_tp)
         where (mask(:))
            self%isperi(:) = src%isperi(:)
            self%peri(:)   = src%peri(:)
            self%atp(:)   = src%atp(:)
         end where
         call util_copy_body(self, src, mask)
      class default
         write(*,*) 'Error copying swiftest_tp object. Incompatible type.'
         call util_exit(FAILURE)
      end select

      return
   end procedure util_copy_tp

   module procedure util_copy_system
      !! author: David A. Minton
      !!
      !! Copies elements of one Swiftest nbody system object to another. The elements are determined by a mask. 
      !! In this case the mask contains the pl elements followed by the tp elements.
      !! 
      implicit none

      associate(pl => self%pl, tp => self%tp, npl => self%pl%nbody, ntp => self%tp%nbody)

         select type(src)
         class is(swiftest_nbody_system)
            call pl%copy(src%pl, mask(1:npl))
            call tp%copy(src%tp, mask(npl+1:npl+ntp))
         class default
            write(*,*) 'Error copying swiftest_nbody_system object. Incompatible type.'
            call util_exit(FAILURE)
         end select


      end associate
   end procedure util_copy_system
end submodule s_util_copy