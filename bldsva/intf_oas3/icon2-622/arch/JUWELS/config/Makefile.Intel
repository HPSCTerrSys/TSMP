export

SHELL = /bin/bash

ARCH = x86_64
OS   = linux-gnu

HOST = juwels01
SITE = fzj.de

srcdir = .

oasisdir = __oasisdir__
include __oasismakefile__
LIBPSMILE=$(ARCHDIR)/lib/libpsmile.${CHAN}.a $(ARCHDIR)/lib/libmct.a $(ARCHDIR)/lib/libmpeu.a $(ARCHDIR)/lib/libscrip.a

prefix = $(shell pwd)
exec_prefix = build/x86_64-unknown-linux-gnu

bindir = ${exec_prefix}/bin
libdir = ${exec_prefix}/lib

NETCDFFROOT    = $(EBROOTNETCDFMINFORTRAN)
NETCDFROOT     = $(EBROOTNETCDF)
NETCDFLIBPATH  = $(NETCDFROOT)/lib64
NETCDF_LIB     = -L$(NETCDFFROOT)/lib -lnetcdff -L$(NETCDFLIBPATH)/ -lnetcdf
NETCDF_INCLUDE = -I$(NETCDFROOT)/include -I$(NETCDFFROOT)/include

CDIROOT     = 
CDILIBPATH  = 
CDI_LIB     = 
CDI_INCLUDE = 

OASIS_LIB      = __cpllib__
OASIS_INCLUDE  = __cplinc__

HDF5ROOT       = $(EBROOTHDF5)
HDF5_LIB       = -L$(HDF5ROOT)/lib -lhdf5_hl -lhdf5
HDF5_INCLUDE   = -I$(HDF5ROOT)/include

SZIPROOT       = $(EBROOTSZIP)
SZIP_LIB       = -L$(SZIPROOT)/lib -lsz
SZIP_INCLUDE   = -I$(SZIPROOT)/include

ZLIBROOT       = $(EBROOTZLIB)
ZLIB_LIB       = -L$(ZLIBROOT)/lib -lz
ZLIB_INCLUDE   = -I$(ZLIBROOT)/include

GRIBAPIROOT    =$(EBROOTGRIB_API)
GRIBAPI_LIB    = -L$(GRIBAPIROOT)/lib -lgrib_api
GRIBAPI_INCLUDE= -I$(GRIBAPIROOT)/include

SCTROOT        = 
SCT_LIB        = 
SCT_INCLUDE    = 

XML2ROOT       = 
XML2_LIB       = 
XML2_INCLUDE   = 

MPIROOT        = $(mpiPath)
MPI_LIB        = -L$(MPIROOT)/lib -L$(MPIROOT)/lib -lmpi
MPI_INCLUDE    = -I$(MPIROOT)/include

LAPACKROOT     = 
LAPACK_LIB_PATH= 
LAPACK_LIB     = -L$(prefix)/$(exec_prefix)/lib -llapack -lblas

PROFILE_LIB     = 
PROFILE_INCLUDE = 

OTHER_LIBS     = 

LIBS           = -L$(prefix)/$(exec_prefix)/lib -lsupport -lself -lmtime $(LAPACK_LIB) $(NETCDF_LIB) $(CDI_LIB) $(HDF5_LIB) $(SZIP_LIB) $(ZLIB_LIB) $(GRIBAPI_LIB) $(XML2_LIB) $(MPI_LIB) $(PROFILE_LIB) $(SCT_LIB) $(OTHER_LIBS) $(OASIS_LIB)
INCLUDE        = -I$(prefix)/$(exec_prefix)/include -I$(prefix)/src/include $(MPI_INCLUDE) $(NETCDF_INCLUDE) $(CDI_INCLUDE) $(HDF5_INCLUDE) $(SZIP_INCLUDE) $(ZLIB_INCLUDE) $(GRIBAPI_INCLUDE) $(XML2_INCLUDE) $(PROFILE_INCLUDE) $(SCT_INCLUDE) $(OASIS_INCLUDE)
INCLUDES       = $(INCLUDE)

AS             = as

CC             = __comCC__
CFLAGS         = $(INCLUDE) -std=gnu99 -march=native -O2 -DHAVE_LIBNETCDF -DHAVE_CF_INTERFACE -DHAVE_LIBSZ -DHAVE_NETCDF4 -DHAVE_LIBHDF5 -DHAVE_NETCDF2 -DpgiFortran -D__ICON__
FC             = __comF90__
FFLAGS         = $(INCLUDES) -J$(prefix)/$(exec_prefix)/module -I$(prefix)/$(exec_prefix)/src -I$(prefix)/$(exec_prefix)/support -I$(prefix)/$(exec_prefix)/module __cplflg__ -g -O2 -msse2 -cpp -std=f2008 -fimplicit-none -fmax-identifier-length=63 -ffree-line-length-512 -Wall -Wcharacter-truncation -Wconversion -Wunderflow -Wunused-parameter -D__LOOP_EXCHANGE -D__ICON__ -D__LOOP_EXCHANGE -D__withoas__ -D__NO_ICON_OCEAN__ -D__NO_JSBACH__ -D__NO_ICON_TESTBED__
FlibFLAGS      = $(INCLUDES) -J$(prefix)/$(exec_prefix)/module 
F77            = __comF77__
F77FLAGS       = __cplflg__ -march=native -g -O2 -ffast-math -mpc64 -cpp -D__LOOP_EXCHANGE -fbacktrace -fstack-protector-all -finit-real=nan -finit-integer=-2147483648 -finit-character=127 -ffree-line-length-512 -D__LOOP_EXCHANGE

AR             = ar
ARFLAGS        = crv

LDFLAGS        = -I$(prefix)/$(exec_prefix)/src -I$(prefix)/$(exec_prefix)/support -I$(prefix)/$(exec_prefix)/module -g -march=native -O3 -ffast-math -D__LOOP_EXCHANGE -D__ICON__

SRCDIRS        = externals/self/src externals/mtime/src blas lapack support src
OBJDIRS        = build/x86_64-unknown-linux-gnu/externals/self/src build/x86_64-unknown-linux-gnu/externals/mtime/src build/x86_64-unknown-linux-gnu/blas build/x86_64-unknown-linux-gnu/lapack build/x86_64-unknown-linux-gnu/support build/x86_64-unknown-linux-gnu/src

ECHO_N = -n

.PHONY: doc

all:
	@for DIR in $(OBJDIRS) ;\
	  do \
	    back=`pwd`; \
	    cd $$DIR && $(MAKE) ; status=$$? ; \
	    if [ $$status != 0 ] ; then \
	      echo "Exit status from make was $$status" ; exit $$status ; \
	    fi ; \
	    cd $$back ; \
	  done 
model:
	@for DIR in $(OBJDIRS) ;\
	  do LASTDIR=$$DIR ;\
	done ;\
	back=`pwd` ;\
	cd $$LASTDIR && $(MAKE) icon  ;\
	cd $$back

one:
	@for DIR in $(OBJDIRS) ;\
	  do LASTDIR=$$DIR ;\
	done ;\
	back=`pwd` ;\
	cd $$LASTDIR && $(MAKE) $(name)  ;\
	cd $$back


install:
	@for DIR in $(OBJDIRS) ;\
	  do \
	  (cd $$DIR && $(MAKE) install ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done

clean:
	@for DIR in $(OBJDIRS) ;\
	  do \
	  (cd $$DIR && $(MAKE) clean ; if [ $$? != 0 ] ; then \
	        echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	  done
	-rm -f ${exec_prefix}/bin/*
	-rm -f ${exec_prefix}/lib/*.a
	-rm -f ${exec_prefix}/module/*.mod  
	-rm -f ${exec_prefix}/src/*.o   
	-rm -rf html/[a-z]*

distclean:
	-rm -rf build
	-rm Makefile
	-rm build_command
	-rm config.log
	-rm config.status
	-rm config/config.h
	-rm config/mh-config-use
	-rm config/set-up.info
	-rm -rf doc/html
	-rm -rf doc/latex
	-rm -rf html/[a-z]*
	-rm -f src/messy

doc:
	doxygen doc/resources/doxyfile_icon_html
	@echo 
	@echo "Start of HTML documentation: doc/html/index.html"
	@echo 

pdf: 
	doxygen doc/resources/doxyfile_icon_pdf

index:
	-rm -rf html/[a-z]*
	scripts/f2html_scripts/f2html.pl -f scripts/f2html_scripts/fgenrc -d html $(SRCDIRS)

checkstyle:
	scripts/code_processing/process_src -v


check:
	@. config/set-up.info ; \
	mpi_total_procs=4 ; \
	for testprog in $(exec_prefix)/bin/test_* ; do \
	echo $(ECHO_N) "$$testprog: " | sed -e 's/^.*\///'; \
	case "$$testprog" in  *_mp|*_mpi) \
	eval $$use_mpi_startrun \"\$$testprog\" && echo "SUCCESS" \
	|| echo "ERROR" ;; \
	*) "$$testprog" && echo "SUCCESS" || echo "ERROR" ;; \
	esac ; done
