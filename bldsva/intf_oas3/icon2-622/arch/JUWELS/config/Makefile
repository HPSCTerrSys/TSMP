export

SHELL = /bin/bash

ARCH = x86_64
OS   = linux-gnu

HOST = $(shell hostname)
SITE = fz-juelich.de

srcdir = .


prefix = $(shell pwd)
exec_prefix = build/x86_64-unknown-linux-gnu

bindir = ${exec_prefix}/bin
libdir = ${exec_prefix}/lib

NETCDFFROOT    = $(EBROOTNETCDFMINFORTRAN)
NETCDFROOT     = $(EBROOTNETCDF)
NETCDFLIBPATH  = $(EBROOTNETCDF)/lib64
NETCDF_LIB     = -L$(NETCDFFROOT)/lib -lnetcdff -L$(NETCDFROOT)/lib64 -lnetcdf
NETCDF_INCLUDE = -I$(EBROOTNETCDF)/include -I$(EBROOTNETCDFMINFORTRAN)/include

CDIROOT     = 
CDILIBPATH  = 
CDI_LIB     = 
CDI_INCLUDE = 

HDF5ROOT       = $(EBROOTHDF5)
HDF5_LIB       = -L$(HDF5ROOT)/lib -lhdf5_hl -lhdf5
HDF5_INCLUDE   = -I$(HDF5ROOT)/include

SZIPROOT       = $(EBROOTSZIP)
SZIP_LIB       = -L$(SZIPROOT)/lib -lsz
SZIP_INCLUDE   = -I$(SZIPROOT)/include

ZLIBROOT       = $(EBROOTZLIB)
ZLIB_LIB       = -L$(ZLIBROOT)/lib -lz
ZLIB_INCLUDE   = -I$(ZLIBROOT)/include

GRIBAPIROOT    = $(EBROOTGRIB_API)
GRIBAPI_LIB    = -L$(GRIBAPIROOT)/lib -lgrib_api
GRIBAPI_INCLUDE= -I$(GRIBAPIROOT)/include

SCTROOT        = 
SCT_LIB        = 
SCT_INCLUDE    = 

XML2ROOT       = 
XML2_LIB       = 
XML2_INCLUDE   = 

YAXTROOT        = 
YAXT_LIB        = 
YAXT_INCLUDE    = 

MPIROOT        = $(EBROOTPSMPI)
MPI_LIB        = -L$(MPIROOT)/lib -L/lib -lmpi -lmpifort
MPI_INCLUDE    = -I$(MPIROOT)/intel64/include -I$(MPIROOT)/include/

LAPACKROOT     = 
LAPACK_LIB_PATH= 
LAPACK_LIB     = -L../lib -llapack -lblas

PROFILE_LIB     = 
PROFILE_INCLUDE = 

SERIALBOX2ROOT        = 
SERIALBOX2_LIB        = 
SERIALBOX2_INCLUDE    = 

OTHER_LIBS     = 

LIBS           = -L../lib -lsupport   -lmtime -lself $(LAPACK_LIB) $(NETCDF_LIB) $(CDI_LIB) $(HDF5_LIB) $(SZIP_LIB) $(ZLIB_LIB) $(GRIBAPI_LIB) $(XML2_LIB) $(YAXT_LIB) $(MPI_LIB) $(PROFILE_LIB) $(SCT_LIB) $(OTHER_LIBS) $(SERIALBOX2_LIB)
INCLUDE        = -I../include -I../../../src/include  $(YAXT_INCLUDE) $(MPI_INCLUDE) $(NETCDF_INCLUDE) $(CDI_INCLUDE) $(HDF5_INCLUDE) $(SZIP_INCLUDE) $(ZLIB_INCLUDE) $(GRIBAPI_INCLUDE) $(XML2_INCLUDE) $(PROFILE_INCLUDE) $(SCT_INCLUDE) $(SERIALBOX2_INCLUDE)
INCLUDES       = $(INCLUDE)

AS             = as

CC             = __comCC__
CFLAGS         = $(INCLUDE) -std=gnu99 -xCORE-AVX2 -O3 -no-prec-div -DHAVE_LIBNETCDF -DHAVE_NETCDF4 -DUSE_MPI -DHAVE_C_NETCDF -DgFortran -D__ICON_ -I/p/scratch/cslts/ghasemi1/tsmp_icon205/tsmp_icon/terrsysmp/icon-nwp_JUWELS_1.4.1MCT_icon/build/x86_64-unknown-linux-gnu/include
FC             = __comF90__
FFLAGS         = $(INCLUDES) -module  ../module -I../src -I../support -I../module -g -xCORE-AVX2 -O3 -no-prec-div -pc64 -fpp -D__LOOP_EXCHANGE -assume realloc_lhs -D__ICON__ -D__NO_JSBACH__ -DUSE_MPI -D__NO_ICON_OCEAN__ -D__NO_ICON_TESTBED__  -fltconsistency -traceback  -D__LOOP_EXCHANGE -DHAVE_FC_ATTRIBUTE_CONTIGUOUS


FlibFLAGS      = $(INCLUDES) -module  ../module 
F77            = __comF77__
F77FLAGS       = -D__LOOP_EXCHANGE -assume realloc_lhs -O3 -no-prec-div -xCORE-AVX2 -fltconsistency -g -pc64 -fpp -D__LOOP_EXCHANGE -D__LOOP_EXCHANGE

AR             = xiar
ARFLAGS        = crv

LDFLAGS        = 

SRCDIRS        = externals/self/src externals/mtime/src  blas lapack support src
OBJDIRS        =  build/x86_64-unknown-linux-gnu/externals/self/src build/x86_64-unknown-linux-gnu/externals/mtime/src build/x86_64-unknown-linux-gnu/blas build/x86_64-unknown-linux-gnu/lapack build/x86_64-unknown-linux-gnu/support build/x86_64-unknown-linux-gnu/src

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
	-rm -f ${exec_prefix}/src/*.a
	-rm -f ${exec_prefix}/src/*.f90
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
	-rm -f src/messy_interface
	-rm -f src/lnd_phy_jsbach 

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
