!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_laplace_coefficient
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description :  Evaluates the Laplace coefficient for alpha<1.0  using F(x) series given in Eq. 6.68 from Murray & Dermott (1999) or
!                 G(y) series (solution to problem 6.2) depending on how close alpha is to 1
!                 optional argument ver lets you specify the method:
!                 ver = 1 : F(x) series
!                 ver = 2 : G(y) series
!                 ver = 3 : Direct integration using Simpson's rule
!
!  Input
!    Arguments : 
!                
!    Teringinal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : kappa = ringmoons_laplace_coefficient(y)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
recursive function ringmoons_laplace_coefficient(alpha,j,s,n) result(ans)

! Modules
      use module_parameters
      IMPLICIT NONE

! Arguments
   real(DP),intent(in) :: alpha,s
   integer, intent(in) :: j,n
   real(DP) :: ans

! Internals
   real(DP) :: T1,T2,T3,T4
   
   if (n==0) then
      ans=laplace(alpha,j,s)
      return
   else if (n==1) then
      T1=laplace(alpha,j-1,s+1._DP)
      T2=-2*alpha*laplace(alpha,j,s+1._DP)
      T3=laplace(alpha,j+1,s+1._DP)
      ans=s*(T1+T2+T3)
      return
   else 
      T1=ringmoons_laplace_coefficient(alpha,j-1,s+1._DP,n-1)
      T2=-2*alpha*ringmoons_laplace_coefficient(alpha,j,s+1._DP,n-1)
      T3=ringmoons_laplace_coefficient(alpha,j+1,s+1._DP,n-1)
      T4=-2*(n-1)*ringmoons_laplace_coefficient(alpha,j,s+1._DP,n-2)
      ans=s*(T1+T2+T3+T4)
   end if


   return
contains

   function laplace(alpha,jp,s,ver) result(ans)

   implicit none
   real(DP),intent(in) :: alpha,s
   integer, intent(in) :: jp
   integer,intent(in),optional :: ver
   real(DP)           :: ans
   real(DP) :: num,denom,T1,tmp,F,F0
   real(DP),parameter :: tol=1e-25_DP
   real(DP),parameter :: alphaG=0.9_DP ! Switchover alpha to go to G series instead of F series
   integer :: i,j,k,v
   real(DP),parameter :: x1=0.4_DP ! For G series near alpha=1. Intermediate value that converges fast using F series
   real(DP),parameter :: x2=0.8_DP ! For G series near alpha=1. Intermediate value that converges fast using F series
   real(DP),dimension(2,2) :: G,Ginv
   real(DP) :: A0,C2,F1,F2,det

   j=abs(jp)
   if (.not.present(ver)) then
      if (alpha<alphaG) then 
         v=1
      else
         v=2
      end if
   else 
      v=ver
   end if
   select case(v)
   case(1) ! F(x) series
      T1=1._DP
      do i=0,j-1
         num=s+i
         denom=(j-i)*1._DP
         T1=T1*num/denom
      end do
      T1=T1*alpha**j
      F = Fseries(alpha**2,s,j,tol)
      ans=2._DP*T1*F
   case(2) ! G(y) series: based on solution to problem 6.2 in Murray & Dermott (1999)

      ! Use an "easy" problem for the F(x) series to bootstrap a solution to the G(y) series coefficients A0 and C2
      G(1,1) = Gseries(1._DP-x1,1._DP,0._DP,s,jp,tol)
      G(1,2) = Gseries(1._DP-x1,0._DP,1._DP,s,jp,tol)
      G(2,1) = Gseries(1._DP-x2,1._DP,0._DP,s,jp,tol)
      G(2,2) = Gseries(1._DP-x2,0._DP,1._DP,s,jp,tol)
      det=(G(1,1)*G(2,2)-G(1,2)*G(2,1))
      Ginv(1,1) = G(2,2)
      Ginv(1,2) = -G(1,2)
      Ginv(2,1) = -G(2,1)
      Ginv(2,2) = G(1,1)
      Ginv = Ginv / det
      F1=Fseries(x1,s,jp,tol)
      F2=Fseries(x2,s,jp,tol)

      A0 = Ginv(1,1)*F1+Ginv(1,2)*F2
      C2 = Ginv(2,1)*F1+Ginv(2,2)*F2
      F = Gseries(1._DP-alpha**2,A0,C2,s,j,tol)

      T1=1._DP
      do i=0,j-1
         num=s+i
         denom=(j-i)*1._DP
         T1=T1*num/denom
      end do
      T1=T1*alpha**j
      ans=2._DP*T1*F

   case(3) ! Simpson's rule (this is very slow, but is included here for testing purposes
      ans=1._DP/pi*simp(bfunci,0._DP,2._DP*pi,10000000,alpha,jp,s)
   end select

   return
   end function laplace

   function Fseries(x,s,jp,tol) result(ans)
   implicit none
   real(DP),intent(in) :: x,s,tol
   integer,intent(in) :: jp
   real(DP) :: ans
   real(DP) :: F0,F,tmp,num,denom
   integer :: i,k

   F0=0._DP
   F=1._DP
   k=1
   do 
      if (abs(1._DP-F/F0)<tol) exit 
      tmp=1.0
      do i=1,k
         num=(s+(i-1))*(s+(jp+i-1))
         denom=i*(jp+i)
         tmp=tmp*num/denom*x
      end do
      F0=F
      F=F+tmp
      k=k+1
   end do
   ans = F
   return
   end function Fseries

   function Gseries(y,A0,C2,s,jp,tol) result(ans)
   ! implements G(y) series (based on answer to problem 6.2 in Murray &  Dermott (1999))
   ! C2 is either B0 (for s=1/2) or A_{2s-1} (for s>1/2)
   implicit none
   real(DP),intent(in) :: y,A0,C2,s,tol
   integer,intent(in) :: jp
   real(DP) :: ans
   integer :: s2,l,tic,toc
   integer,parameter :: lmax=10000000
   real(DP) :: G0,Bl2s1,Bl2s2,Bl,Al
   real(DP),dimension(:),allocatable :: Bold

   s2 = nint(2*s)
   allocate(Bold(s2))
   ans = 0._DP
   G0=-1._DP
   Al=A0
   l = 0
   if (s2==1) then
      Bl=C2
      do l=0,lmax
         G0 = ans
         ans = ans + Al*y**(l-s2+1) + log(y)*Bl*y**l
         if (abs(1._DP-ans/G0)<tol) exit 
         Bold(1)=Bl
         Bl=Blp1(Bl,s,jp,l)
         Al=Alp1(Al,Bold(1),Bl,s,jp,l)
      end do
   else
      Bl=B0func(A0,s,jp)
      Bl2s1=0._DP
      Bl2s2=0._DP
      tic=1
      Bold(tic) = Bl
      do l=0,lmax
         G0 = ans
         ans = ans + Al*y**(l-s2+1) + log(y)*Bl*y**l
         if ((l>s2-2) .and. (abs(1._DP-ans/G0)<tol)) exit 
         if (l/=s2-2) then 
            if (l-s2+1>=0) then
               toc = tic - s2 + 1
               if (toc <= 0) toc = size(Bold) + toc
               Bl2s1 = Bold(toc)
            end if
            if (l-s2+2>=0) then
               toc = tic - s2 + 2
               if (toc <= 0) toc = size(Bold) + toc
               Bl2s2 = Bold(toc)
            end if
            Al=Alp1(Al,Bl2s1,Bl2s2,s,jp,l)
         else
            Al = C2
         end if
         Bl=Blp1(Bl,s,jp,l)
         tic = tic + 1
         if (tic > size(Bold)) tic = 1
         Bold(tic)=Bl ! tic keeps a looping array of old Bl values
      end do

   end if


   deallocate(Bold)

   return
   end function Gseries

   function Alp1(Al,Bl2s1,Bl2s2,s,jp,l) result(ans)
   ! evaluates B_l coefficients (based on answer to problem 6.2 in Murray &  Dermott (1999))
   implicit none
   real(DP),intent(in) :: Al,Bl2s1,Bl2s2,s
   integer,intent(in) :: jp,l
   real(DP) :: ans
   integer :: s2,n,start

   s2=nint(2*s)
   ans = Al*((l-s2+1)*(l+jp+1._DP)+s*(s+jp))+Bl2s1*(2*l-s2+jp+2)-Bl2s2*(2*l-s2+3)
   ans = ans / ((l+1)*(l-s2+2._DP))

   return
   end function Alp1


   function Blp1(Bl,s,jp,l) result(ans)
   ! evaluates B_l coefficients (based on answer to problem 6.2 in Murray & Dermott (1999) )
   implicit none
   real(DP),intent(in) :: Bl,s
   integer,intent(in) :: jp,l
   real(DP) :: ans
   integer :: n

   ans = Bl*(l*(2*s+jp+l)+s*(s+jp))/((l+1)*(2*s+l))

   return
   end function Blp1

   function B0func(A0,s,jp) result(ans)
   implicit none
   real(DP),intent(in) :: A0,s
   integer,intent(in) :: jp
   real(DP) :: ans
   integer :: n,s2

   s2 = nint(2*s)
   ans = A0
   do n=1,s2-2
      ans = ans*((n-s2)*(n+jp)+s*(s+jp))/(n*(n-s2+1._DP))
   end do
   ans = ans * (1._DP-s2-jp+s*(s+jp))/(s2-1._DP)
  
   return
   end function B0func



   real(DP) function simp(f,a,b,m,alpha,jp,s)
   implicit none
   real(DP), external :: f
   real(DP), intent(in) :: a,b,alpha,s
   integer, intent(in) :: m,jp
   real(DP) :: f1,f2,f3,f4,ans,h
   integer :: k

   h=(b-a)/(3.0*m)

   ans=0.0
   do k=1,m
      f1=f(a+(3._DP*k-3._DP)*h,alpha,jp,s)
      f2=f(a+(3._DP*k-2._DP)*h,alpha,jp,s)
      f3=f(a+(3._DP*k-1._DP)*h,alpha,jp,s)
      f4=f(a+3._DP*k*h,alpha,jp,s)
      ans = ans +  (f1+3._DP*f2+3._DP*f3+f4)
   end do
   simp = ans * 3._DP/8._DP * h
   return 
   end function

   real(DP) function bfunci(psi,alpha,jp,s)
   implicit none
   integer,intent(in) :: jp
   real(DP), intent(in) :: s,psi,alpha
   
   bfunci=cos(real(jp,8)*psi)/(1._DP-2._DP*alpha*cos(psi)+alpha**2)**(s)
   
   return
   end function
   



end function ringmoons_laplace_coefficient
