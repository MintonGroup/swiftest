module symba_classes
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Adapted from David E. Kaufmann's Swifter routine: helio.f90
   use swiftest_globals
   use helio_classes, only : helio_cb, helio_pl, helio_tp, helio_nbody_system
   implicit none


   !********************************************************************************************************************************
   !  symba_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   type, public, extends(helio_nbody_system) :: symba_nbody_system
   contains
   end type symba_nbody_system

   !********************************************************************************************************************************
   ! symba_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Helio central body particle class
   type, public, extends(helio_cb) :: symba_cb
   contains
   end type symba_cb

   !********************************************************************************************************************************
   !                                    symba_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio massive body particle class
   type, public, extends(helio_pl) :: symba_pl
   contains
   end type symba_pl

   !********************************************************************************************************************************
   !                                    symba_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio test particle class
   type, public, extends(helio_tp) :: symba_tp
   contains
   end type symba_tp
end module symba_classes