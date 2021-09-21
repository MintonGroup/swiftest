#******************************************************************************
#
#  Unit Name   : Makefile
#  Unit Type   : makefile
#  Project     : SWIFTEST
#  Package     : N/A
#  Language    : GNU makefile syntax
#
#  Description : Controls, via the make program, the building of the Swifter
#                modules, library, drivers, and tools, as well as initiating
#                the build of the FXDR library by means of its own makefile
#
#  Input
#    Arguments : Zero or more of the following targets:
#                (1) all     : builds modules, entire Swifter library, FXDR
#                              library, Swifter drivers and tools
#                (2) mod     : builds modules
#                (3) lib     : builds entire Swifter library
#                (4) libdir  : compiles local directory source and adds the
#                              resulting objects to the Swifter library
#                (5) drivers : builds Swifter drivers
#                (6) tools   : builds Swifter tools
#                (7) bin     : compiles local directory source and installs
#                              resulting executables to $(SWIFTEST_HOME)/bin
#                (8) clean   : removes all soft links to Makefile and
#                              Makefile.Defines from subdirectories of
#                              $(SWIFTEST_HOME), removes the entire contents
#                              of $(SWIFTEST_HOME)/lib and $(SWIFTEST_HOME)/bin,
#                              and removes the include file installed by the
#                              FXDR makefile
#    Terminal  : none
#    File      : Makefile.Defines
#
#  Output
#    Arguments : none
#    Terminal  : status messages
#    File      : none
#
#  Invocation  : make [all|mod|lib|libdir|drivers|tools|bin|clean]
#
#  Notes       : The use of the above arguments as phony targets inside the
#                makefile precludes their use as base names of Swifter drivers
#                or tools
#
#******************************************************************************

SWIFTEST_MODULES =   swiftest_globals.f90 \
                     swiftest_operators.f90 \
							lambda_function.f90\
                     swiftest_classes.f90 \
                     fraggle_classes.f90 \
                     whm_classes.f90 \
                     rmvs_classes.f90 \
                     helio_classes.f90 \
                     symba_classes.f90 \
							walltime_classes.f90 \
                     swiftest.f90 


include Makefile.Defines

MKL_ROOT = /apps/spack/bell/apps/intel-parallel-studio/cluster.2019.5-intel-19.0.5-4brgqlf/mkl/lib
IMKL = -I$(MKLROOT)/include
LMKL = -L$(MKLROOT)/lib/intel64 -qopt-matmul 

MODULES         = $(SWIFTEST_MODULES) $(USER_MODULES) 

.PHONY : all mod lib libdir fast drivers bin clean force 

% : %.f90 force
	$(FORTRAN) $(FFLAGS) -I$(SWIFTEST_HOME)/include -I$(NETCDF_FORTRAN_HOME)/include $(IMKL) $< -o $@ \
	  -L$(SWIFTEST_HOME)/lib -lswiftest -L$(NETCDF_FORTRAN_HOME)/lib -lnetcdf -lnetcdff $(LMKL)
	$(INSTALL_PROGRAM) $@ $(SWIFTEST_HOME)/bin
	rm -f $@

all:
	cd $(SWIFTEST_HOME); \
	  make mod; \
	  make lib; \
	  make fast; \
	  make drivers; \

mod:
	cd $(SWIFTEST_HOME)/src/modules/; \
	  $(FORTRAN) $(FFLAGS) -I$(SWIFTEST_HOME)/include -I$(NETCDF_FORTRAN_HOME)/include $(IMKL) -c $(MODULES); \
	  $(AR) rv $(SWIFTEST_HOME)/lib/libswiftest.a *.o; \
	  $(INSTALL_DATA) *.mod *.smod $(SWIFTEST_HOME)/include; \
	  rm -f *.o *.mod  *.smod

lib:
	cd $(SWIFTEST_HOME)/src/discard; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/gr; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/io; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/kick; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/netcdf; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/obl; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/operators; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/setup; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/tides; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/whm; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/rmvs; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/helio; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/symba; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/user; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTEST_HOME)/src/walltime; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make libdir

fast:
	cd $(SWIFTEST_HOME)/src/fraggle; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make fastdir

	cd $(SWIFTEST_HOME)/src/util; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make fastdir

	cd $(SWIFTEST_HOME)/src/orbel; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make fastdir

	cd $(SWIFTEST_HOME)/src/drift; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make fastdir

	cd $(SWIFTEST_HOME)/src/helio; \
		$(FORTRAN) $(FFASTFLAGS) -I$(SWIFTEST_HOME)/include -I$(NETCDF_FORTRAN_HOME)/include $(IMKL) -c helio_drift.f90; \
		$(AR) rv $(SWIFTEST_HOME)/lib/libswiftest.a *.o *.smod; \
		$(INSTALL_DATA) *.smod $(SWIFTEST_HOME)/include; \
		rm -f *.o *.smod	

	cd $(SWIFTEST_HOME)/src/rmvs; \
		$(FORTRAN) $(FFASTFLAGS) -I$(SWIFTEST_HOME)/include -I$(NETCDF_FORTRAN_HOME)/include $(IMKL) -c rmvs_encounter_check.f90; \
		$(AR) rv $(SWIFTEST_HOME)/lib/libswiftest.a *.o *.smod; \
		$(INSTALL_DATA) *.smod $(SWIFTEST_HOME)/include; \
		rm -f *.o *.smod	

	cd $(SWIFTEST_HOME)/src/symba; \
		$(FORTRAN) $(FFASTFLAGS) -I$(SWIFTEST_HOME)/include -I$(NETCDF_FORTRAN_HOME)/include $(IMKL) -c symba_encounter_check.f90; \
		$(AR) rv $(SWIFTEST_HOME)/lib/libswiftest.a *.o *.smod; \
		$(INSTALL_DATA) *.smod $(SWIFTEST_HOME)/include; \
		rm -f *.o *.smod

libdir:
	$(FORTRAN) $(FFLAGS) -I$(SWIFTEST_HOME)/include -I$(NETCDF_FORTRAN_HOME)/include $(IMKL) -c *.f90; \
	$(AR) rv $(SWIFTEST_HOME)/lib/libswiftest.a *.o *.smod; \
	$(INSTALL_DATA) *.smod $(SWIFTEST_HOME)/include; \
	rm -f *.o *.smod

fastdir:
	$(FORTRAN) $(FFASTFLAGS) -I$(SWIFTEST_HOME)/include -I$(NETCDF_FORTRAN_HOME)/include $(IMKL) -c *.f90; \
	$(AR) rv $(SWIFTEST_HOME)/lib/libswiftest.a *.o *.smod; \
	$(INSTALL_DATA) *.smod $(SWIFTEST_HOME)/include; \
	rm -f *.o *.smod

drivers:
	cd $(SWIFTEST_HOME)/src/main; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTEST_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTEST_HOME)/Makefile .; \
	  make bin

bin: *.f90
	make $(basename $^)

clean:
	cd $(SWIFTEST_HOME)/src/modules;  rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/discard; rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/drift;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/fraggle;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/gr;          rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/helio;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/io;      rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/kick;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/main;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/netcdf;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/obl;     rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/operators;      rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/orbel;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/rmvs;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/setup;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/symba;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/tides;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/user;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/util;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/walltime;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/src/whm;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTEST_HOME)/bin;     rm -f swiftest_*
	cd $(SWIFTEST_HOME)/bin;     rm -f tool_*
	cd $(SWIFTEST_HOME)/lib;     rm -f lib*
	cd $(SWIFTEST_HOME)/include; rm -f *.mod *.smod 
	cd $(COLLRESOLVE_HOME); rm -rf autom4te.cache aux Makefile stamp-h1 configure config.status config.h config.log aclocal.m4 lib* *.in *.o *.lo cambioni2019/*.o cambioni2019/*.lo


force:

#******************************************************************************
#
#  Author(s)   : David E. Kaufmann
#
#  Revision Control System (RCS) Information
#
#  Source File : $RCSfile: Makefile,v $
#  Full Path   : $Source: /d1/kaufmann/development/RCS/Makefile,v $
#  Revision    : $Revision: 0.1 $
#  Date        : $Date: 2003/04/15 22:56:34 $
#  Programmer  : $Author: kaufmann $
#  Locked By   : $Locker: kaufmann $
#  State       : $State: Exp $
#
#  Modification History:
#
#  $Log: Makefile,v $
#  Revision 0.1  2003/04/15 22:56:34  kaufmann
#  Initial implementation
#
#
#******************************************************************************


