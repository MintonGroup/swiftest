submodule (swiftest_classes) s_step

contains
   module procedure step_system
      !! author: David A. Minton
      !!
      !! Wrapper function that selects the correct step function for the integrator type.
      !! The nested select statements serve to make this method a "bridge" between the polymorphic swiftest_nbody_system class
      !! in which the cb, pl, and tp components are allocatable abstract classes and the actual integrator-specific methods that are
      !! called internally. Before this point, the actual types of cb, pl, and tp are ambiguous. The select type constructs remove the 
      !! ambiguity.
      !! 
      use swiftest
      implicit none

      select type(system => self)
      class is (whm_nbody_system) 
         select type(cb => self%cb)
         class is (whm_cb)
            select type(pl => self%pl)
            class is (whm_pl)
               select type(tp => self%tp)
               class is (whm_tp)
                  call whm_step_system(cb, pl, tp, config)
               end select
            end select
         end select
      class is (rmvs_nbody_system) 
         select type(cb => self%cb)
         class is (rmvs_cb)
            select type(pl => self%pl)
            class is (rmvs_pl)
               select type(tp => self%tp)
               class is (rmvs_tp)
                  call rmvs_step_system(cb, pl, tp, config)
               end select
            end select
         end select
      class is (helio_nbody_system) 
         select type(cb => self%cb)
         class is (helio_cb)
            select type(pl => self%pl)
            class is (helio_pl)
               select type(tp => self%tp)
               class is (helio_tp)
                  call helio_step_system(cb, pl, tp, config)
               end select
            end select
         end select
      end select
      return
   end procedure step_system

end submodule s_step
