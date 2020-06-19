submodule (nbody_data_structures) s_nbody_read_initial_conditions_file
contains
   module procedure nbody_read_initial_conditions_file
      !! author: David A. Minton
      !!
      !! Wrapper to select which file input method to use for either tp or pl 
      use swiftest
      implicit none

      select type(self)
      type is (helio_tp)
         call io_read_tp_in(self,config)
      type is (symba_tp)
         call io_read_tp_in(self,config)
      type is (helio_pl)
         call io_read_pl_in(self,config)
      type is (symba_pl)
         call io_read_pl_in(self,config)
      end select

      return
   end procedure nbody_read_initial_conditions_file

end submodule s_nbody_read_initial_conditions_file

