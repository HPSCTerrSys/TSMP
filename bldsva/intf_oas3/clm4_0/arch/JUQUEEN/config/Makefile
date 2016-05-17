.SUFFIXES:
.SUFFIXES: .F90 .o .c .f90
# name of macros file - but default this is generic

VPATH := $(OAS3_DIR)

# Determine whether to compile threaded or not
ifeq ($(strip $(BUILD_THREADED)),TRUE)
   compile_threaded = true
endif
ifeq ($(strip $(SMP)),TRUE)
   compile_threaded = true
endif


ifndef MOD_SUFFIX
   MOD_SUFFIX := mod
endif


ifeq ($(strip $(MACFILE)),)
   MACFILE  := Macros
endif

# Machine specific macros file 
# This must be included before any settings are overwritten 
# But must be AFTER any definitions it uses are defined.
# So be careful if moving this either earlier or later in the makefile!!!
include $(MACFILE)

INCLDIR += __comflg__
INCLDIR += -I$(EXEROOT)/csm_share -I$(EXEROOT)/mct/mct -I$(EXEROOT)/mct/mpeu 
INCLDIR += -I$(EXEROOT)/gptl -I$(EXEROOT)/pio
INCLDIR += -I$(CCSMROOT)/models/lnd/clm/src/util_share/
LDFLAGS += __ldflg__


ifeq ($(strip $(MPILIB)), mpi-serial)
  CC      := $(SCC)
  FC      := $(SFC)
  MPIFC   := $(SFC)
  MPICC   := $(SCC)
  INCLDIR += -I$(EXEROOT)/mct/mpi-serial
else
  CC := $(MPICC)
  FC := $(MPIFC)
  CPPDEFS += -WF,-DHAVE_MPI
endif
ifdef CPRE
  FPPDEFS := $(patsubst -D%,$(CPRE)%,$(CPPDEFS)) 
else
  FPPDEFS := $(CPPDEFS)
endif



OBJS = oas_clm_vardef.o oas_clm_init.o oas_clm_snd.o  \
       oas_clm_rcv.o oas_clm_define.o oas_clm_finalize.o

liboas3.a: $(OBJS)
	$(AR) ruv $@ $(OBJS)



.c.o:
	$(CC) -c $(INCLDIR) $(INCS) $(CFLAGS) $(CPPDEFS) $<
.F.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FPPDEFS) $(FIXEDFLAGS) $<
.f90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS)  $<
.F90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FPPDEFS) $(FREEFLAGS)  $<

mostlyclean:
	$(RM) -f *.f *.f90 

clean:
	$(RM) -f *.f *.f90 *.d *.$(MOD_SUFFIX) $(OBJS)


install: liboas3.a
	cp -p *.$(MOD_SUFFIX) $(INCROOT) 
	cp -p liboas3.a $(LIBROOT)

