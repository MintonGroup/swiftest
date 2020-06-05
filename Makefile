#******************************************************************************
#
#  Unit Name   : Makefile
#  Unit Type   : makefile
#  Project     : SWIFTER
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
#                (5) fxdr    : builds FXDR library by invoking its makefile
#                (6) drivers : builds Swifter drivers
#                (7) tools   : builds Swifter tools
#                (8) bin     : compiles local directory source and installs
#                              resulting executables to $(SWIFTER_HOME)/bin
#                (9) clean   : removes all soft links to Makefile and
#                              Makefile.Defines from subdirectories of
#                              $(SWIFTER_HOME), removes the entire contents
#                              of $(SWIFTER_HOME)/lib and $(SWIFTER_HOME)/bin,
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
#  Invocation  : make [all|mod|lib|libdir|fxdr|drivers|tools|bin|clean]
#
#  Notes       : The use of the above arguments as phony targets inside the
#                makefile precludes their use as base names of Swifter drivers
#                or tools
#
#******************************************************************************

SWIFTER_MODULES = module_parameters.f90 \
		  module_swiftest.f90 \
		  module_swifter.f90 \
		  module_helio.f90 \
        module_nrutil.f90 \
		  module_symba.f90 \
		  module_fxdr.f90 \
		  module_swiftestalloc.f90 \
        module_interfaces.f90

include Makefile.Defines

MODULES         = $(SWIFTER_MODULES) $(USER_MODULES)

.PHONY : all mod lib libdir fxdr collresolve drivers tools bin clean force 

% : %.f90 force
	$(FORTRAN) $(FFLAGS) -I$(SWIFTER_HOME)/include $< -o $@ \
	  -L$(SWIFTER_HOME)/lib -lswifter -lfxdr -lcollresolve
	$(INSTALL_PROGRAM) $@ $(SWIFTER_HOME)/bin
	rm -f $@

all:
	cd $(SWIFTER_HOME); \
	  make mod; \
	  make lib; \
	  make fxdr; \
	  make collresolve; \
	  make drivers; \
	  make tools

mod:
	cd $(SWIFTER_HOME)/src/module; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  $(FORTRAN) $(FFLAGS) -I$(SWIFTER_HOME)/include -c $(MODULES); \
	  $(AR) rv $(SWIFTER_HOME)/lib/libswifter.a *.o; \
	  $(INSTALL_DATA) *.mod $(SWIFTER_HOME)/include; \
	  rm -f *.o *.mod

lib:
	cd $(SWIFTER_HOME)/src/coord; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTER_HOME)/src/discard; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTER_HOME)/src/drift; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTER_HOME)/src/helio; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTER_HOME)/src/io; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTER_HOME)/src/obl; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTER_HOME)/src/orbel; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTER_HOME)/src/rmvs; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTER_HOME)/src/symba; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir
	cd $(SWIFTER_HOME)/src/util; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make libdir

libdir:
	$(FORTRAN) $(FFLAGS) -I$(SWIFTER_HOME)/include -c *.f90
	$(AR) rv $(SWIFTER_HOME)/lib/libswifter.a *.o
	rm -f *.o

fxdr:
	cd $(SWIFTER_HOME)/src/fxdr; \
	  rm -f Makefile.Defines; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
     chmod -w test_read_only.xdr ; \
	  make -f Makefile.fxdr; \
	  make -f Makefile.fxdr test; \
	  make -f Makefile.fxdr install; \
	  make -f Makefile.fxdr clean

collresolve:
	cd $(COLLRESOLVE_HOME); \
	  autoreconf --install;\
	  ./configure --prefix=$(SWIFTER_HOME);\
	  make; \
	  make install


drivers:
	cd $(SWIFTER_HOME)/src/main; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make bin

tools:
	cd $(SWIFTER_HOME)/src/tool; \
	  rm -f Makefile.Defines Makefile; \
	  ln -s $(SWIFTER_HOME)/Makefile.Defines .; \
	  ln -s $(SWIFTER_HOME)/Makefile .; \
	  make bin

bin: *.f90
	make $(basename $^)

clean:
	cd $(SWIFTER_HOME)/src/module;  rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/coord;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/discard; rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/drift;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/helio;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/io;      rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/obl;     rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/orbel;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/rmvs;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/symba;   rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/util;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/fxdr;    rm -f Makefile.Defines  *.gc*
	cd $(SWIFTER_HOME)/src/main;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/src/tool;    rm -f Makefile.Defines Makefile *.gc*
	cd $(SWIFTER_HOME)/bin;     rm -f swifter_*
	cd $(SWIFTER_HOME)/bin;     rm -f tool_*
	cd $(SWIFTER_HOME)/lib;     rm -f lib*
	cd $(SWIFTER_HOME)/include; rm -f *.mod fxdr.inc collresolve.h
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
