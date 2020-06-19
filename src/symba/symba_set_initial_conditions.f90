submodule (symba) s_symba_set_initial_conditions
contains
   module procedure symba_set_initial_conditions
      !! author: David A. Minton
      !!
      !! Sets up initial conditions for a run. Currently, it reads in all the ICs from input files, but future versions could also
      !!    initialize them in some other way
      use swiftest
      implicit none

      ! read in the total number of bodies from the input files
      call symba_plA%read_from_file(config)
      call symba_tpA%read_from_file(config)

      ! Save central body mass in vector form so that elemental functions can be evaluated with it
      call symba_tpA%set_vec(symba_plA%mass(1),config%dt)
      call symba_plA%set_vec(symba_plA%mass(1),config%dt)

      ! Save system mass to both objects
      call symba_plA%set_msys(symba_plA)
      call symba_tpA%set_msys(symba_plA)

      ! create arrays of data structures big enough to store the number of bodies we are adding
      call mergeadd_list%alloc(10*npl)!DM: Why 10*npl?
      call mergesub_list%alloc(npl)
      call plplenc_list%alloc(10*npl)!DM: See ^
      call pltpenc_list%alloc(ntp)!DM: See ^

      ! reads in initial conditions of all massive bodies from input file
      ! reorder by mass 
      call symba_plA%reorder()
      call util_valid(symba_plA, symba_tpA)
   
      return
   end procedure symba_set_initial_conditions
end submodule s_symba_set_initial_conditions
