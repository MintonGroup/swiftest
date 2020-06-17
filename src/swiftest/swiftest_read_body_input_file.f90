submodule (swiftest_data_structures) s_swiftest_read_body_input_file
contains
   module procedure swiftest_read_body_input_file
   !! author: David A. Minton
   !!
   !! Reads the input file, selecting whether to use the tp or pl method depending on type
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
   end procedure swiftest_read_body_input_file

end submodule s_swiftest_read_body_input_file

