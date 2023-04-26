####################################################
#
#    Makefile for the Fortran 90 program
#
####################################################
#
#
####################################################
#
#    Declaration of Variables
#
####################################################
#
.SILENT:
#
STDROOT      = __cosmoroot__
SRCDIR       = $(STDROOT)/src
OASDIR       = $(STDROOT)/src/oas3
OBJDIR       = $(STDROOT)/obj
WORKDIR      = $(STDROOT)/work
WITHOAS      = __withoas__
#
#
####################################################
#  include (or build) an appropriate compiler file
#     some compiler files are available in LOCAL
#     and can be copied to Fopts
####################################################

include $(STDROOT)/Fopts

####################################################
#  Declaration of the Object Files
####################################################

include $(STDROOT)/ObjFiles
include $(STDROOT)/ObjFiles.twomom
####################################################
#  Actions that can be performed
####################################################
#
exe:    purpar
#
purpar:  check_is_not_mpifh / $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ) $(OASOBJ)
	echo linking parallel program
	( $(LDPAR) $(LDFLG) -o $(PROGRAM)_pur $(DATAOBJ) $(UTILOBJ) $(EXTOBJ) \
                                    $(PUBLICUTILS) $(SRCOBJ) $(OASOBJ) $(LIB) $(LIBPATH) )
#
allpar:  check_is_not_mpifh / $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ) $(ASSOBJ) $(SATOBJ)
	echo linking parallel program
	( $(LDPAR) $(LDFLG) -o $(PROGRAM)_all $(DATAOBJ) $(UTILOBJ) $(EXTOBJ) \
                                    $(PUBLICUTILS) $(SRCOBJ) $(ASSOBJ) $(SATOBJ) $(LIB) )
#
nudpar:  check_is_not_mpifh / $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ) $(ASSOBJ)
	echo linking parallel program
	( $(LDPAR) $(LDFLG) -o $(PROGRAM)_nud $(DATAOBJ) $(UTILOBJ) $(EXTOBJ) \
                                    $(PUBLICUTILS) $(SRCOBJ) $(ASSOBJ) $(LIB) )
#
satpar:  check_is_not_mpifh / $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ) $(SATOBJ)
	echo linking parallel program
	( $(LDPAR) $(LDFLG) -o $(PROGRAM)_sat $(DATAOBJ) $(UTILOBJ) $(EXTOBJ) \
                                    $(PUBLICUTILS) $(SRCOBJ) $(SATOBJ) $(LIB) )
#
###############################################################################
#
purseq:  check_is_mpifh / $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ)  \
         $(OBJDIR)/dummy_mpi.o
	echo linking sequential program
	( $(LDSEQ) $(LDFLG) -o $(PROGRAM)_purseq $(DATAOBJ) $(UTILOBJ)           \
              $(OBJDIR)/dummy_mpi.o $(EXTOBJ) $(PUBLICUTILS) $(SRCOBJ) $(LIB) )
#
allseq:  check_is_mpifh / $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ)  \
         $(ASSOBJ) $(SATOBJ) $(OBJDIR)/dummy_mpi.o
	( $(LDSEQ) $(LDFLG) -o $(PROGRAM)_allseq $(DATAOBJ) $(UTILOBJ)           \
              $(OBJDIR)/dummy_mpi.o $(EXTOBJ) $(PUBLICUTILS) $(SRCOBJ) $(ASSOBJ) $(SATOBJ) $(LIB) )
#
nudseq:  check_is_mpifh / $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ)  \
         $(ASSOBJ) $(OBJDIR)/dummy_mpi.o
	echo linking sequential program
	( $(LDSEQ) $(LDFLG) -o $(PROGRAM)_nudseq $(DATAOBJ) $(UTILOBJ)           \
              $(OBJDIR)/dummy_mpi.o $(EXTOBJ) $(PUBLICUTILS) $(SRCOBJ) $(ASSOBJ) $(LIB) )
#
satseq:  check_is_mpifh / $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ)  \
         $(SATOBJ) $(OBJDIR)/dummy_mpi.o
	echo linking sequential program
	( $(LDSEQ) $(LDFLG) -o $(PROGRAM)_satseq $(DATAOBJ) $(UTILOBJ)           \
              $(OBJDIR)/dummy_mpi.o $(EXTOBJ) $(PUBLICUTILS) $(SRCOBJ) $(SATOBJ) $(LIB) )
#
###############################################################################
#
$(OBJDIR):
	mkdir -p $(OBJDIR)
#
check_is_mpifh:
	if [ ! -f $(SRCDIR)/mpif.h ]; then  \
          echo "ERROR: mpif.h is not contained in src/ !";  \
          exit 1;  \
        fi
#
check_is_not_mpifh:
	if [ -f $(SRCDIR)/mpif.h ]; then  \
          echo "ERROR: there is an mpif.h in src/ !";  \
          exit 1;  \
        fi
#
addmods:
	for DATEI in `ls -x $(WORKDIR)` ; do \
	   cp -p $(WORKDIR)/$${DATEI}  $(SRCDIR) ; \
	done
#
clean:
	echo cleaning up
	rm -f $(PROGRAM)*
	rm -f $(OBJDIR)/*
#
####################################################
#
#     Dependencies of the Data Modules
#
####################################################

include ObjDependencies.twomom
include ObjDependencies
