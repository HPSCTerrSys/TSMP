####################################################
#
#    Makefile for the Fortran 90 program
#
####################################################
#
#    Architecture + configuration
#
#    cluma       : HRZ cluster
#    cluma_debug : HRZ cluster with debug options
#    sx9         : NEC SX9 at DWD
#    sx9_debug   : NEC SX9 at DWD with debug options
#    blgq        : Ibm Blue Gene/Q at JSC
#
####################################################
#
MACHINE_CFG  = blgq
#
AGB_MODS = -DHYMACS 
#
####################################################
#
#    Changes to original version
#
####################################################
# additional modules:
#
# HYMACS, VK 2012/03/15
# ---------------------
# precomp. flag: -DHYMACS
# global vars  : hymacs_...
# modules+dep. :
#
# src_conv_hymacs
#    organize_physics
# data_hymacs
#    src_allocation
#    src_conv_hymacs
#    src_setup_vartab
#    src_runge_kutta
#
# TRACER, AK 2012/04/11
# ---------------------
# global vars : tracer...
# modules+dep. :
#
# data_tracer
#   src_allocation
#   src_setup_vartab
#   lmorg
#   src_tracer_supply
#   src_advection_rk
#   src_conv_shallow
#   organize_physics
#   src_conv_hymacs
#   src_conv_tiedtke
#   src_runge_kutta
#   src_slow_tendencies
#
# src_tracer_supply
#   lmorg
#   src_advection_rk
#   src_conv_shallow
#   src_conv_hymacs
#
####################################################
#
#    Declaration of variables and configurations
#
####################################################
#
.SILENT:
#
STDROOT      = ${PWD}
SRCDIR       = $(STDROOT)/src
OASDIR       = $(STDROOT)/src/oas3
OBJDIR       = $(STDROOT)/obj
WORKDIR      = $(SVAROOT)/run 
PROGRAM      = lmparbin
#
####################################################
#  HRZ cluster production
####################################################
ifeq ($(MACHINE_CFG), cluma)
INCDIR       = ${NETCDF}/include 
# compiler
F90          = ${MPIDIR}/bin/mpif90
COMFLG1      = -c -O0 -I$(INCDIR) -cpp -DGRIBDWD $(AGB_MODS) ${CPPDEF2}
COMFLG2      = $(COMFLG1)
COMFLG3      = $(COMFLG1)
COMFLG4      = $(COMFLG1)
COMFLG5      = $(COMFLG1) ${PSMILE_INCDIR} ${CPPDEF1} ${CPPDEF2}
#
# loader
LDPAR        = ${MPIDIR}/bin/mpif90 
LDSEQ        = 
LDFLG        =
#
# libs
LIB          =  ${LIBGRIB1_COS}/libgrib1.a \
                ${LIBOASIS} ${LIBPSMILE} -L${NETCDF}/lib -lnetcdf -lnetcdff
endif
####################################################
#  HRZ cluster debug
####################################################
ifeq ($(MACHINE_CFG), cluma_debug)
INCDIR       = ${NETCDF}/include 
# compiler
F90          = ${MPIDIR}/bin/mpif90 
COMFLG1      = -c -O0 -I$(INCDIR) -g -C -debug all -traceback -warn all -warn nounused -cpp -DGRIBDWD $(AGB_MODS) ${CPPDEF2}
COMFLG2      = $(COMFLG1)
COMFLG3      = $(COMFLG1)
COMFLG4      = $(COMFLG1)
COMFLG5      = $(COMFLG1) ${PSMILE_INCDIR} ${CPPDEF1} ${CPPDEF2}
#
# loader
LDPAR        = ${MPIDIR}/bin/mpif90
LDSEQ        = 
LDFLG        =
#
# libs
#
LIB          = ${LIBGRIB1_COS}/libgrib1.a \
                ${LIBOASIS} ${LIBPSMILE} -L${NETCDF}/lib -lnetcdf -lnetcdff 
              
endif
####################################################
#  BlueGeneQ JUQUEEN production
####################################################
ifeq ($(MACHINE_CFG), blgq)
INCDIR       = ${NETCDF}/include

# compiler
F90          = ${PREP_C} ${MPIDIR}/bin/mpixlf90_r
F902         = ${PREP_C} ${MPIDIR}/bin/mpixlf90_r
# 
#
# linker
LDPAR        = ${PREP_C} ${MPIDIR}/bin/mpixlf90_r
LDSEQ        = ${PREP_C} ${MPIDIR}/bin/mpixlf90_r
#
COMFLG1      = -c $(OPT_C) -qextname=flush -qsuffix=cpp=f90 -WF,-DGRIBDWD -WF,$(AGB_MODS) -I${INCDIR} ${CPPDEF2} 
COMFLG2      = $(COMFLG1)
COMFLG3      = $(COMFLG1)
COMFLG4      = $(COMFLG1) 
COMFLG5      = $(COMFLG1) ${PSMILE_INCDIR} ${CPPDEF1}

# link-flags
LDFLG        =                                     # VK: Linux
#
LIB          =  ${LIBGRIB1_COS}/libgrib1.a \
                ${LIBOASIS} ${LIBPSMILE} -L${NETCDF}/lib -lnetcdf  
endif
####################################################
#  NEC SX production
####################################################
ifeq ($(MACHINE_CFG), sx9)
# compiler
F90          = sxmpif90 -clear -c -EP -w                                                                      \
               -I $(SRCDIR) -I/usr/local/pkg-sx9/include  -I/e/uhome/uschaett/include                         \
               -Wf"-pvctl neighbors vr256 noassume loopcnt=2000000 vwork=stack                                \
                          vworksz=100M novlchk noverrchk chgpwr nomsg"                                        \
               -Wf,"-K a -M setall flovf inv zdiv"                                                            \
               -Wf"-L fmtlist transform summary mrgmsg"                                                       \
               -DNECSX -DGRIBDWD
#              -DNUDGING -DNETCDF -DRTTOV9
#              -Wf"-s -msg o"                                                                                 \
#              -Wf"-L fmtlist transform summary mrgmsg"                                                       \
#              -Wf,"-O nomsg -init stack=nan heap=nan"                                                        \
#              -Wf,"-O nomsg"                                                                                 \
#
LDPAR        = sxmpif90
LDSEQ        =
#
COMFLG1      = -Cvopt -Pstack                                                                                           \
               -Wf,"-i noerrchk -O move" -Wf,"-pvctl cond_mem_opt ifopt loopfusion outerstrip"                          \
               -pi auto exp=udsdx,udsdx_up5_xy,coe_th,cu_cond,flake_radflux,flake_driver,flake_buoypar                  \
                   -pi  exp=flake_snowdensity,flake_snowheatconduct,SfcFlx_rhoair                                       \
                 expin=../src/numeric_utilities_rk.f90                                                                  \
                 expin=../src/src_flake.f90                                                                             \
                 expin=../src/numeric_utilities.f90                                                                     \
                       noexp=mean_over_box,mean_cov_over_box,vert_avg,interpol_sl_trilin,lap_2,lap_4am,lap_4aml         \
                 expin=../src/src_radiation.f90                                                                         \
                 expin=../src/utilities.f90                                                                             \
                 expin=../src/src_flake.f90                                                                             \
                 expin=../src/src_flake_sfcflx.f90
COMFLG2      = -Cvopt -Pstack                                                                                           \
               -pi auto exp=udsdx,udsdx_up5_xy,coe_th,cu_cond                                                           \
                 expin=../src/numeric_utilities_rk.f90                                                                  \
                 expin=../src/numeric_utilities.f90 noexp=mean_over_box,mean_cov_over_box,vert_avg,interpol_sl_trilin   \
                 expin=../src/src_radiation.f90                                                                         \
                 expin=../src/utilities.f90
COMFLG3      = -Cvsafe -Pstack -Wf,"-pvctl noneighbors"
COMFLG4      = $(COMFLG3)
#
LDFLG       =                
AR          = sxar
LIB          = /usr/local/pkg-sx9/libdwd/lib/libgrib1.a                    \
               /usr/local/pkg-sx9/libdwd/lib/libmisc.a                     \
               /usr/local/pkg-sx9/for0adm/lib/librttov_ifc.a               \
               /usr/local/pkg-sx9/for0adm/lib/librttov9.3.a                \
               /usr/local/pkg-sx9/for0adm/lib/librttov9.3_parallel.a       \
               /usr/local/pkg-sx9/lib/libnetcdf.a
#              /usr/local/pkg-sx9/for0adm/lib/libRTTOV7_synsat_vector.a    
endif
####################################################
#  NEC SX debug
####################################################
ifeq ($(MACHINE_CFG), sx9_debug)
# compiler
#
F90          = sxmpif90 -clear -c -EP -w                                                                      \
               -I $(SRCDIR) -I/usr/local/pkg-sx9/include  -I/e/uhome/uschaett/include                         \
               -Wf"-pvctl neighbors vr256 noassume loopcnt=2000000 vwork=stack                                \
                          vworksz=100M novlchk noverrchk chgpwr nomsg"                                        \
               -Wf,"-K a -M setall flovf inv zdiv"                                                            \
               -Wf"-L fmtlist transform summary mrgmsg"                                                       \
               -DNECSX -DGRIBDWD
#              -DNUDGING -DNETCDF -DRTTOV9
#              -Wf"-s -msg o"                                                                                 \
#              -Wf"-L fmtlist transform summary mrgmsg"                                                       \
#              -Wf,"-O nomsg -init stack=nan heap=nan"                                                        \
#              -Wf,"-O nomsg"                                                                                 \
#
LDPAR        = sxmpif90
LDSEQ        =
#
COMFLG1      = -Cvopt -Pstack                                                                                           \
               -Wf,"-i noerrchk -O move" -Wf,"-pvctl cond_mem_opt ifopt loopfusion outerstrip"                          \
               -pi auto exp=udsdx,udsdx_up5_xy,coe_th,cu_cond,flake_radflux,flake_driver,flake_buoypar                  \
                   -pi  exp=flake_snowdensity,flake_snowheatconduct,SfcFlx_rhoair                                       \
                 expin=../src/numeric_utilities_rk.f90                                                                  \
                 expin=../src/src_flake.f90                                                                             \
                 expin=../src/numeric_utilities.f90                                                                     \
                       noexp=mean_over_box,mean_cov_over_box,vert_avg,interpol_sl_trilin,lap_2,lap_4am,lap_4aml         \
                 expin=../src/src_radiation.f90                                                                         \
                 expin=../src/utilities.f90                                                                             \
                 expin=../src/src_flake.f90                                                                             \
                 expin=../src/src_flake_sfcflx.f90
COMFLG2     = -g -eC -Wf"-init stack=nan heap=nan"
COMFLG3      = -Cvsafe -Pstack -Wf,"-pvctl noneighbors"
COMFLG4      = $(COMFLG3)
#
LDFLG       = -p -ftrace -Wl,"-d -m -sum" -Wl,"-f nan"
AR          = sxar
LIB          = /usr/local/pkg-sx9/libdwd/lib/libgrib1.a                    \
               /usr/local/pkg-sx9/libdwd/lib/libmisc.a                     \
               /usr/local/pkg-sx9/for0adm/lib/librttov_ifc.a               \
               /usr/local/pkg-sx9/for0adm/lib/librttov9.3.a                \
               /usr/local/pkg-sx9/for0adm/lib/librttov9.3_parallel.a       \
               /usr/local/pkg-sx9/lib/libnetcdf.a
#              /usr/local/pkg-sx9/for0adm/lib/libRTTOV7_synsat_vector.a  
endif
####################################################
#
#     Declaration of the Object Files
#
####################################################
#
DATAOBJ      =   $(OBJDIR)/data_1d_global.o          \
                 $(OBJDIR)/data_constants.o          \
                 $(OBJDIR)/data_convection.o         \
                 $(OBJDIR)/data_fields.o             \
                 $(OBJDIR)/data_flake.o              \
                 $(OBJDIR)/data_gscp.o               \
                 $(OBJDIR)/data_io.o                 \
                 $(OBJDIR)/data_modelconfig.o        \
                 $(OBJDIR)/data_parallel.o           \
                 $(OBJDIR)/data_parameters.o         \
                 $(OBJDIR)/data_radiation.o          \
                 $(OBJDIR)/data_runcontrol.o         \
                 $(OBJDIR)/data_soil.o               \
                 $(OBJDIR)/data_turbulence.o         \
                 $(OBJDIR)/data_hymacs.o             \
                 $(OBJDIR)/data_tracer.o
#
UTILOBJ      =   $(OBJDIR)/environment.o             \
                 $(OBJDIR)/info_lm_f90.o             \
                 $(OBJDIR)/io_utilities.o            \
                 $(OBJDIR)/meteo_utilities.o         \
                 $(OBJDIR)/numeric_utilities.o       \
                 $(OBJDIR)/numeric_utilities_rk.o    \
                 $(OBJDIR)/parallel_utilities.o      \
                 $(OBJDIR)/pp_utilities.o            \
                 $(OBJDIR)/time_utilities.o          \
                 $(OBJDIR)/turbulence_utilities.o    \
                 $(OBJDIR)/utilities.o
#
EXTOBJ        =  $(OBJDIR)/dfi_initialization.o      \
                 $(OBJDIR)/hori_diffusion.o          \
                 $(OBJDIR)/fast_waves.o              \
                 $(OBJDIR)/fast_waves_rk.o           \
                 $(OBJDIR)/near_surface.o            \
                 $(OBJDIR)/organize_data.o           \
                 $(OBJDIR)/organize_diagnosis.o      \
                 $(OBJDIR)/organize_dynamics.o       \
                 $(OBJDIR)/organize_eps.o            \
                 $(OBJDIR)/organize_physics.o        \
                 $(OBJDIR)/semi_implicit.o           \
                 $(OBJDIR)/slow_tendencies.o         \
                 $(OBJDIR)/lmorg.o
#
SRCOBJ       =   $(OBJDIR)/src_2timelevel.o          \
                 $(OBJDIR)/src_advection_rk.o        \
                 $(OBJDIR)/src_allocation.o          \
                 $(OBJDIR)/src_artifdata.o           \
                 $(OBJDIR)/src_conv_tiedtke.o        \
                 $(OBJDIR)/src_conv_kainfri.o        \
                 $(OBJDIR)/src_conv_shallow.o        \
                 $(OBJDIR)/src_conv_hymacs.o         \
                 $(OBJDIR)/src_diagbudget.o          \
                 $(OBJDIR)/src_flake.o               \
                 $(OBJDIR)/src_gridpoints.o          \
                 $(OBJDIR)/src_gscp.o                \
                 $(OBJDIR)/src_input.o               \
                 $(OBJDIR)/src_integrals.o           \
                 $(OBJDIR)/src_flake_sfcflx.o        \
                 $(OBJDIR)/src_leapfrog.o            \
                 $(OBJDIR)/src_meanvalues.o          \
                 $(OBJDIR)/src_output.o              \
                 $(OBJDIR)/src_radiation.o           \
                 $(OBJDIR)/src_relaxation.o          \
                 $(OBJDIR)/src_runge_kutta.o         \
                 $(OBJDIR)/src_seaice.o              \
                 $(OBJDIR)/src_setup.o               \
                 $(OBJDIR)/src_setup_vartab.o        \
                 $(OBJDIR)/src_slow_tendencies_rk.o  \
                 $(OBJDIR)/src_soil.o                \
                 $(OBJDIR)/src_soil_multlay.o        \
                 $(OBJDIR)/src_spectral_nudging.o    \
                 $(OBJDIR)/src_sso.o                 \
                 $(OBJDIR)/src_tracer_supply.o       \
                 $(OBJDIR)/src_turbulence.o          \
                 $(OBJDIR)/turbulence_interface.o    \
                 $(OBJDIR)/turbulence_tran.o         \
                 $(OBJDIR)/turbulence_diff.o
#
ASSOBJ =         $(OBJDIR)/data_nudge_all.o          \
                 $(OBJDIR)/data_lheat_nudge.o        \
                 $(OBJDIR)/data_lhn_diag.o           \
                 $(OBJDIR)/data_nudge_gather.o       \
                 $(OBJDIR)/data_nudge_local.o        \
                 $(OBJDIR)/data_nudge_spread.o       \
                 $(OBJDIR)/data_obs_process.o        \
                 $(OBJDIR)/data_obs_record.o         \
                 $(OBJDIR)/organize_assimilation.o   \
                 $(OBJDIR)/src_correl_cutoff.o       \
                 $(OBJDIR)/src_gather_info.o         \
                 $(OBJDIR)/src_lheat_nudge.o         \
                 $(OBJDIR)/src_lheating.o            \
                 $(OBJDIR)/src_mult_local.o          \
                 $(OBJDIR)/src_mult_spread.o         \
                 $(OBJDIR)/src_nudging.o             \
                 $(OBJDIR)/src_obs_1dvar_org.o       \
                 $(OBJDIR)/src_obs_proc_air.o        \
                 $(OBJDIR)/src_obs_proc_aof.o        \
                 $(OBJDIR)/src_obs_proc_cdf.o        \
                 $(OBJDIR)/src_obs_proc_util.o       \
                 $(OBJDIR)/src_obs_processing.o      \
                 $(OBJDIR)/src_sfcana.o              \
                 $(OBJDIR)/src_sing_local.o          \
                 $(OBJDIR)/src_sing_spread.o
#
SATOBJ =         $(OBJDIR)/data_satellites.o         \
                 $(OBJDIR)/src_sat_tbs.o             \
                 $(OBJDIR)/src_sat_rttov.o           \
                 $(OBJDIR)/organize_satellites.o
#
OASOBJ       =   $(OBJDIR)/oas_cos_vardef.o          \
                 $(OBJDIR)/oas_cos_define.o          \
                 $(OBJDIR)/oas_cos_init.o            \
                 $(OBJDIR)/oas_cos_finalize.o        \
                 $(OBJDIR)/oas_cos_rcv.o             \
                 $(OBJDIR)/oas_cos_snd.o             \
                 $(OBJDIR)/receive_fld_2clm.o        \
                 $(OBJDIR)/send_fld_2clm.o
#
PUBLICUTILS =    $(OBJDIR)/mpe_io.o                  \
                 $(OBJDIR)/dummy_db.o
#
####################################################
#
#     Actions that can be performed
#
####################################################
#
exe:    purexe
#
purexe:  $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ) $(OASOBJ)
	echo linking parallel program
	if [ -f $(SRCDIR)/mpif.h ]; then \
	  echo "ERROR: there is an mpif.h in src/ !"; \
	  exit 1; \
	fi
	( $(LDPAR) $(LDFLG) -o $(PROGRAM)_pur $(DATAOBJ) $(UTILOBJ) $(EXTOBJ) \
                                    $(PUBLICUTILS) $(SRCOBJ) $(OASOBJ) $(LIB) )
#
allexe:  $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ) $(ASSOBJ) $(SATOBJ)
	echo linking parallel program
	if [ -f $(SRCDIR)/mpif.h ]; then \
	  echo "ERROR: there is an mpif.h in src/ !"; \
	  exit 1; \
	fi
	( $(LDPAR) $(LDFLG) -o $(PROGRAM)_all $(DATAOBJ) $(UTILOBJ) $(EXTOBJ) \
                                    $(PUBLICUTILS) $(SRCOBJ) $(ASSOBJ) $(SATOBJ) $(LIB) )
#
nudexe:  $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ) $(ASSOBJ)
	echo linking parallel program
	if [ -f $(SRCDIR)/mpif.h ]; then \
	  echo "ERROR: there is an mpif.h in src/ !"; \
	  exit 1; \
	fi
	( $(LDPAR) $(LDFLG) -o $(PROGRAM)_nud $(DATAOBJ) $(UTILOBJ) $(EXTOBJ) \
                                    $(PUBLICUTILS) $(SRCOBJ) $(ASSOBJ) $(LIB) )
#
satexe:  $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ) $(SATOBJ)
	echo linking parallel program
	if [ -f $(SRCDIR)/mpif.h ]; then \
	  echo "ERROR: there is an mpif.h in src/ !"; \
	  exit 1; \
	fi
	( $(LDPAR) $(LDFLG) -o $(PROGRAM)_sat $(DATAOBJ) $(UTILOBJ) $(EXTOBJ) \
                                    $(PUBLICUTILS) $(SRCOBJ) $(SATOBJ) $(LIB) )
#
seqexe:  $(OBJDIR) $(DATAOBJ) $(PUBLICUTILS) $(UTILOBJ) $(EXTOBJ) $(SRCOBJ) \
         $(OBJDIR)/dummy_mpi.o
	echo linking sequential program
	if [ ! -f $(SRCDIR)/mpif.h ]; then  \
	  echo "ERROR: mpif.h is not contained in src/ !";  \
	  exit 1;  \
	fi
	( $(LDSEQ) $(LDFLG) -o $(PROGRAM)_seq $(DATAOBJ) $(PUBLICUTILS)         \
         $(UTILOBJ) $(OBJDIR)/dummy_mpi.o $(EXTOBJ) $(SRCOBJ) $(LIB) )
#
clean:
	echo cleaning up
	rm -f $(PROGRAM)
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJDIR)/i.*
	rm -f $(OBJDIR)/*.T
	rm -f $(OBJDIR)/*.L
	rm -f $(OBJDIR)/*.mod
	rm -f $(OBJDIR)/*.lst
#
####################################################
#
#     Dependencies of the Data Modules
#
####################################################
#
#
$(OBJDIR)/data_1d_global.o:            $(SRCDIR)/data_1d_global.f90           \
           $(OBJDIR)/data_parameters.o
	echo compiling data_1d_global.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_1d_global.f90 )
#
$(OBJDIR)/data_constants.o:            $(SRCDIR)/data_constants.f90           \
           $(OBJDIR)/data_parameters.o
	echo compiling data_constants.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_constants.f90 )
#
$(OBJDIR)/data_convection.o:           $(SRCDIR)/data_convection.f90         \
           $(OBJDIR)/data_parameters.o
	echo compiling data_convection.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_convection.f90 )
#
$(OBJDIR)/data_fields.o:               $(SRCDIR)/data_fields.f90              \
           $(OBJDIR)/data_parameters.o
	echo compiling data_fields.f90      
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_fields.f90 )
#
$(OBJDIR)/data_flake.o:                $(SRCDIR)/data_flake.f90               \
           $(OBJDIR)/data_parameters.o
	echo compiling data_flake.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_flake.f90 )
#
$(OBJDIR)/data_gscp.o:                 $(SRCDIR)/data_gscp.f90               \
           $(OBJDIR)/data_parameters.o
	echo compiling data_gscp.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_gscp.f90 )
#
$(OBJDIR)/data_io.o:                   $(SRCDIR)/data_io.f90                  \
           $(OBJDIR)/data_parameters.o
	echo compiling data_io.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_io.f90 )
#
$(OBJDIR)/data_lheat_nudge.o:          $(SRCDIR)/data_lheat_nudge.f90         \
           $(OBJDIR)/data_parameters.o
	echo compiling data_lheat_nudge.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_lheat_nudge.f90 )
#
$(OBJDIR)/data_lhn_diag.o:             $(SRCDIR)/data_lhn_diag.f90            \
           $(OBJDIR)/data_parameters.o
	echo compiling data_lhn_diag.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_lhn_diag.f90 )
#
$(OBJDIR)/data_modelconfig.o:          $(SRCDIR)/data_modelconfig.f90         \
           $(OBJDIR)/data_parameters.o
	echo compiling data_modelconfig.f90 
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_modelconfig.f90 )
#
$(OBJDIR)/data_nudge_all.o:            $(SRCDIR)/data_nudge_all.f90           \
           $(OBJDIR)/data_parameters.o
	echo compiling data_nudge_all.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_nudge_all.f90 )
#
$(OBJDIR)/data_nudge_gather.o:         $(SRCDIR)/data_nudge_gather.f90        \
           $(OBJDIR)/data_parameters.o
	echo compiling data_nudge_gather.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_nudge_gather.f90 )
#
$(OBJDIR)/data_nudge_local.o:          $(SRCDIR)/data_nudge_local.f90         \
           $(OBJDIR)/data_parameters.o
	echo compiling data_nudge_local.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_nudge_local.f90 )
#
$(OBJDIR)/data_nudge_spread.o:         $(SRCDIR)/data_nudge_spread.f90        \
           $(OBJDIR)/data_parameters.o
	echo compiling data_nudge_spread.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_nudge_spread.f90 )
#
$(OBJDIR)/data_obs_process.o:          $(SRCDIR)/data_obs_process.f90         \
           $(OBJDIR)/data_parameters.o
	echo compiling data_obs_process.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_obs_process.f90 )
#
$(OBJDIR)/data_obs_record.o:           $(SRCDIR)/data_obs_record.f90          \
           $(OBJDIR)/data_parameters.o
	echo compiling data_obs_record.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_obs_record.f90 )
#
$(OBJDIR)/data_parallel.o:             $(SRCDIR)/data_parallel.f90            \
           $(OBJDIR)/data_parameters.o
	echo compiling data_parallel.f90  
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_parallel.f90 )
#
$(OBJDIR)/data_parameters.o:           $(SRCDIR)/data_parameters.f90
	echo compiling data_parameters.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_parameters.f90 )
#
$(OBJDIR)/data_radiation.o:            $(SRCDIR)/data_radiation.f90           \
           $(OBJDIR)/data_parameters.o
	echo compiling data_radiation.f90   
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_radiation.f90 )
#
$(OBJDIR)/data_runcontrol.o:           $(SRCDIR)/data_runcontrol.f90          \
           $(OBJDIR)/data_parameters.o
	echo compiling data_runcontrol.f90  
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_runcontrol.f90 )
#
$(OBJDIR)/data_satellites.o:           $(SRCDIR)/data_satellites.f90          \
           $(OBJDIR)/data_parameters.o
	echo compiling data_satellites.f90  
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_satellites.f90 )
#
$(OBJDIR)/data_soil.o:                 $(SRCDIR)/data_soil.f90                \
           $(OBJDIR)/data_parameters.o
	echo compiling data_soil.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_soil.f90 )
#
$(OBJDIR)/data_turbulence.o:           $(SRCDIR)/data_turbulence.f90          \
           $(OBJDIR)/data_parameters.o
	echo compiling data_turbulence.f90  
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_turbulence.f90 )
#
$(OBJDIR)/data_hymacs.o:       $(SRCDIR)/data_hymacs.f90      \
           $(OBJDIR)/data_parameters.o
	echo compiling data_hymacs.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_hymacs.f90 )
#
$(OBJDIR)/data_tracer.o:           $(SRCDIR)/data_tracer.f90      \
           $(OBJDIR)/data_parameters.o
	echo compiling data_tracer.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/data_tracer.f90 )
#
####################################################
#
#     Dependencies of OASIS utility 
#
####################################################
#
$(OBJDIR)/oas_cos_vardef.o:                  $(OASDIR)/oas_cos_vardef.F90     \
	$(OBJDIR)/data_parameters.o
	echo compiling oas_cos_vardef.F90
	( cd $(OBJDIR) && $(F90) $(COMFLG5) $(OASDIR)/oas_cos_vardef.F90 )
#
####################################################
#
#     Dependencies of the Utility Modules
#
####################################################
#
$(OBJDIR)/environment.o:		$(SRCDIR)/environment.f90        \
	$(OBJDIR)/data_parameters.o        $(OBJDIR)/oas_cos_vardef.o
	echo compiling environment.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG5) $(SRCDIR)/environment.f90 )
#
$(OBJDIR)/info_lm_f90.o:                     $(SRCDIR)/info_lm_f90.f90
	echo compiling info_lm_f90.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/info_lm_f90.f90 )
#
$(OBJDIR)/io_utilities.o:                    $(SRCDIR)/io_utilities.f90       \
          $(OBJDIR)/data_parameters.o                                         \
          $(OBJDIR)/mpe_io.o
	echo compiling io_utilities.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/io_utilities.f90 )
#
$(OBJDIR)/meteo_utilities.o:                 $(SRCDIR)/meteo_utilities.f90    \
          $(OBJDIR)/data_parameters.o
	echo compiling meteo_utilities.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/meteo_utilities.f90 )
#
$(OBJDIR)/numeric_utilities.o:               $(SRCDIR)/numeric_utilities.f90  \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/numeric_utilities_rk.o
	echo compiling numeric_utilities.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/numeric_utilities.f90 )
#
$(OBJDIR)/numeric_utilities_rk.o:           $(SRCDIR)/numeric_utilities_rk.f90 \
          $(OBJDIR)/data_parameters.o       $(OBJDIR)/data_runcontrol.o        \
          $(OBJDIR)/data_parallel.o         $(OBJDIR)/data_modelconfig.o       \
          $(OBJDIR)/data_fields.o           $(OBJDIR)/data_constants.o         \
          $(OBJDIR)/environment.o           $(OBJDIR)/parallel_utilities.o
	echo compiling numeric_utilities_rk.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/numeric_utilities_rk.f90 )
#
$(OBJDIR)/parallel_utilities.o:              $(SRCDIR)/parallel_utilities.f90 \
          $(OBJDIR)/data_parameters.o
	echo compiling parallel_utilities.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/parallel_utilities.f90 )
#
$(OBJDIR)/pp_utilities.o:                    $(SRCDIR)/pp_utilities.f90       \
          $(OBJDIR)/data_parameters.o                                         \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_constants.o       \
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/numeric_utilities.o    \
          $(OBJDIR)/environment.o
	echo compiling pp_utilities.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/pp_utilities.f90 )
#
$(OBJDIR)/turbulence_utilities.o:            $(SRCDIR)/turbulence_utilities.f90 \
          $(OBJDIR)/data_turbulence.o        $(OBJDIR)/data_parameters.o
	echo compiling turbulence_utilities.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/turbulence_utilities.f90 )
#
$(OBJDIR)/time_utilities.o:                  $(SRCDIR)/time_utilities.f90     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o
	echo compiling time_utilities.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/time_utilities.f90 )
#
$(OBJDIR)/utilities.o:                       $(SRCDIR)/utilities.f90          \
          $(OBJDIR)/data_parameters.o
	echo compiling utilities.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/utilities.f90 )
#
####################################################
#
#     Dependencies of the Public Utilities
#
####################################################
#
$(OBJDIR)/mpe_io.o:                          $(SRCDIR)/mpe_io.f90
	echo compiling mpe_io.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/mpe_io.f90 )
#
$(OBJDIR)/dummy_db.o:                        $(SRCDIR)/dummy_db.f90
	echo compiling dummy_db.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/dummy_db.f90 )
#
$(OBJDIR)/dummy_rttov.o:                     $(SRCDIR)/dummy_rttov.f90
	echo compiling dummy_rttov.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/dummy_rttov.f90 )
#
$(OBJDIR)/dummy_netcdf.o:                    $(SRCDIR)/dummy_netcdf.f90
	echo compiling dummy_netcdf.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/dummy_netcdf.f90 )
#
# only on sequential platforms
#
$(OBJDIR)/dummy_mpi.o:                       $(SRCDIR)/dummy_mpi.f90
	echo compiling dummy_mpi.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/dummy_mpi.f90 )
#
####################################################
#
#     Dependencies of the Source Modules
#
####################################################
#
$(OBJDIR)/src_2timelevel.o:                  $(SRCDIR)/src_2timelevel.f90     \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_lheat_nudge.o       $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_tracer.o            $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/meteo_utilities.o      \
          $(OBJDIR)/numeric_utilities.o      $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/src_lheating.o           $(OBJDIR)/time_utilities.o       \
          $(OBJDIR)/hori_diffusion.o
	echo compiling src_2timelevel.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_2timelevel.f90 )
#
$(OBJDIR)/src_advection_rk.o:                $(SRCDIR)/src_advection_rk.f90   \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o                                         \
          $(OBJDIR)/environment.o            $(OBJDIR)/hori_diffusion.o       \
          $(OBJDIR)/meteo_utilities.o        $(OBJDIR)/numeric_utilities.o    \
          $(OBJDIR)/numeric_utilities_rk.o   $(OBJDIR)/parallel_utilities.o
	echo compiling src_advection_rk.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_advection_rk.f90 )
#
$(OBJDIR)/src_allocation.o:                  $(SRCDIR)/src_allocation.f90     \
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_lheat_nudge.o       $(OBJDIR)/data_lhn_diag.o        \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
                                             $(OBJDIR)/data_hymacs.o          \
          $(OBJDIR)/data_tracer.o
	echo compiling src_allocation.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_allocation.f90 )
#
$(OBJDIR)/src_artifdata.o:                   $(SRCDIR)/src_artifdata.f90      \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_soil.o            \
          $(OBJDIR)/environment.o            $(OBJDIR)/meteo_utilities.o      \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/utilities.o
	echo compiling src_artifdata.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_artifdata.f90 )
#
$(OBJDIR)/src_conv_kainfri.o:                $(SRCDIR)/src_conv_kainfri.f90   \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o
	echo compiling src_conv_kainfri.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_conv_kainfri.f90 )
#
$(OBJDIR)/src_conv_tiedtke.o:                $(SRCDIR)/src_conv_tiedtke.f90   \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_convection.o      \
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_lheat_nudge.o     \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/src_lheating.o           $(OBJDIR)/data_tracer.o
	echo compiling src_conv_tiedtke.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_conv_tiedtke.f90 )
#
$(OBJDIR)/src_conv_shallow.o:                $(SRCDIR)/src_conv_shallow.f90   \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_convection.o      \
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_tracer.o
	echo compiling src_conv_shallow.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_conv_shallow.f90 )
#
$(OBJDIR)/src_conv_hymacs.o:                 $(SRCDIR)/src_conv_hymacs.f90    \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_convection.o      \
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/environment.o            $(OBJDIR)/data_hymacs.o          \
          $(OBJDIR)/data_tracer.o
	echo compiling src_conv_hymacs.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_conv_hymacs.f90 )
#
$(OBJDIR)/src_correl_cutoff.o:               $(SRCDIR)/src_correl_cutoff.f90  \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_nudge_gather.o      $(OBJDIR)/data_nudge_spread.o    \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o
	echo compiling src_correl_cutoff.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_correl_cutoff.f90 )
#
$(OBJDIR)/src_diagbudget.o:                  $(SRCDIR)/src_diagbudget.f90     \
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o
	echo compiling src_diagbudget.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_diagbudget.f90 )
#
$(OBJDIR)/src_flake.o:                       $(SRCDIR)/src_flake.f90          \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_flake.o                                              \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o                                         \
          $(OBJDIR)/src_flake_sfcflx.o
	echo compiling src_flake.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_flake.f90 )
#
$(OBJDIR)/src_gather_info.o:                 $(SRCDIR)/src_gather_info.f90    \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_nudge_gather.o      $(OBJDIR)/data_nudge_spread.o    \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/environment.o          \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/time_utilities.o
	echo compiling src_gather_info.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_gather_info.f90 )
#
$(OBJDIR)/src_gridpoints.o:                  $(SRCDIR)/src_gridpoints.f90     \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/data_soil.o              $(OBJDIR)/environment.o          \
          $(OBJDIR)/utilities.o
	echo compiling src_gridpoints.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_gridpoints.f90 )
#
$(OBJDIR)/src_gscp.o:                        $(SRCDIR)/src_gscp.f90           \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_gscp.o              $(OBJDIR)/data_lheat_nudge.o     \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/meteo_utilities.o        $(OBJDIR)/pp_utilities.o         \
          $(OBJDIR)/environment.o            $(OBJDIR)/src_lheating.o
	echo compiling src_gscp.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_gscp.f90 )
#
$(OBJDIR)/src_input.o:                       $(SRCDIR)/src_input.f90          \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_soil.o            \
          $(OBJDIR)/environment.o            $(OBJDIR)/io_utilities.o         \
          $(OBJDIR)/meteo_utilities.o        $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/time_utilities.o         $(OBJDIR)/utilities.o
	echo compiling src_input.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_input.f90 )
#
$(OBJDIR)/src_integrals.o:                   $(SRCDIR)/src_integrals.f90      \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/environment.o
	echo compiling src_integrals.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_integrals.f90 )
#
$(OBJDIR)/src_flake_sfcflx.o:                $(SRCDIR)/src_flake_sfcflx.f90   \
          $(OBJDIR)/data_flake.o             $(OBJDIR)/data_parameters.o
	echo compiling src_flake_sfcflx.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_flake_sfcflx.f90 )
#
#
$(OBJDIR)/src_leapfrog.o:                    $(SRCDIR)/src_leapfrog.f90       \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_lheat_nudge.o       $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_tracer.o            $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/meteo_utilities.o      \
          $(OBJDIR)/numeric_utilities.o      $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/src_lheating.o           $(OBJDIR)/time_utilities.o       \
          $(OBJDIR)/hori_diffusion.o
	echo compiling src_leapfrog.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_leapfrog.f90 )
#
$(OBJDIR)/src_lheat_nudge.o:                 $(SRCDIR)/src_lheat_nudge.f90    \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_lheat_nudge.o     \
          $(OBJDIR)/data_lhn_diag.o          $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/environment.o          \
          $(OBJDIR)/meteo_utilities.o        $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/pp_utilities.o           $(OBJDIR)/utilities.o            \
          $(OBJDIR)/time_utilities.o
	echo compiling src_lheat_nudge.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_lheat_nudge.f90 )
#
$(OBJDIR)/src_lheating.o:                    $(SRCDIR)/src_lheating.f90       \
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_lheat_nudge.o     \
                                             $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o
	echo compiling src_lheating.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_lheating.f90 )
#
$(OBJDIR)/src_meanvalues.o:                  $(SRCDIR)/src_meanvalues.f90     \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/utilities.o
	echo compiling src_meanvalues.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_meanvalues.f90 )
#
$(OBJDIR)/src_mult_local.o:                  $(SRCDIR)/src_mult_local.f90     \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_nudge_gather.o      $(OBJDIR)/data_nudge_local.o     \
          $(OBJDIR)/data_obs_record.o        $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o
	echo compiling src_mult_local.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_mult_local.f90 )
#
$(OBJDIR)/src_mult_spread.o:                 $(SRCDIR)/src_mult_spread.f90    \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_nudge_gather.o      $(OBJDIR)/data_nudge_spread.o    \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/src_correl_cutoff.o
	echo compiling src_mult_spread.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_mult_spread.f90 )
#
$(OBJDIR)/src_nudging.o:                     $(SRCDIR)/src_nudging.f90        \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_lheat_nudge.o     \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_nudge_spread.o      $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/meteo_utilities.o      \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/src_gather_info.o      \
          $(OBJDIR)/src_lheating.o           $(OBJDIR)/src_mult_local.o       \
          $(OBJDIR)/src_mult_spread.o        $(OBJDIR)/src_obs_1dvar_org.o    \
          $(OBJDIR)/src_obs_processing.o     $(OBJDIR)/src_obs_proc_cdf.o     \
          $(OBJDIR)/src_sing_local.o         $(OBJDIR)/src_sing_spread.o      \
          $(OBJDIR)/time_utilities.o
	echo compiling src_nudging.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_nudging.f90 )
#
$(OBJDIR)/src_obs_1dvar_org.o:               $(SRCDIR)/src_obs_1dvar_org.f90  \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_obs_process.o       $(OBJDIR)/data_obs_record.o      \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o                                         \
          $(OBJDIR)/environment.o            $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/src_obs_proc_util.o      $(OBJDIR)/utilities.o
	echo compiling src_obs_1dvar_org.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_obs_1dvar_org.f90 )
#
$(OBJDIR)/src_obs_proc_air.o:                $(SRCDIR)/src_obs_proc_air.f90   \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_obs_process.o       $(OBJDIR)/data_obs_record.o      \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o                                         \
          $(OBJDIR)/environment.o            $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/src_obs_proc_util.o      $(OBJDIR)/time_utilities.o
	echo compiling src_obs_proc_air.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_obs_proc_air.f90 )
#
$(OBJDIR)/src_obs_proc_aof.o:                $(SRCDIR)/src_obs_proc_aof.f90   \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_obs_process.o       $(OBJDIR)/data_obs_record.o      \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/environment.o          \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/src_obs_proc_util.o    \
          $(OBJDIR)/utilities.o
	echo compiling src_obs_proc_aof.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_obs_proc_aof.f90 )
#
$(OBJDIR)/src_obs_proc_cdf.o:                $(SRCDIR)/src_obs_proc_cdf.f90   \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_obs_process.o       $(OBJDIR)/data_obs_record.o      \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/environment.o          \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/src_obs_proc_util.o    \
          $(OBJDIR)/src_obs_proc_air.o       $(OBJDIR)/src_obs_processing.o   \
          $(OBJDIR)/utilities.o
	echo compiling src_obs_proc_cdf.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_obs_proc_cdf.f90 )
#
$(OBJDIR)/src_obs_proc_util.o:               $(SRCDIR)/src_obs_proc_util.f90  \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_obs_process.o       $(OBJDIR)/data_obs_record.o      \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o                                         \
          $(OBJDIR)/environment.o            $(OBJDIR)/parallel_utilities.o
	echo compiling src_obs_proc_util.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_obs_proc_util.f90 )
#
$(OBJDIR)/src_obs_processing.o:              $(SRCDIR)/src_obs_processing.f90 \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_obs_process.o       $(OBJDIR)/data_obs_record.o      \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/environment.o          \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/src_obs_proc_air.o     \
          $(OBJDIR)/src_obs_proc_aof.o       $(OBJDIR)/src_obs_proc_util.o    \
          $(OBJDIR)/utilities.o
	echo compiling src_obs_processing.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_obs_processing.f90 )
#
$(OBJDIR)/src_output.o:                      $(SRCDIR)/src_output.f90         \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_satellites.o      \
          $(OBJDIR)/src_artifdata.o                                           \
          $(OBJDIR)/environment.o            $(OBJDIR)/io_utilities.o         \
          $(OBJDIR)/numeric_utilities.o                                       \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/pp_utilities.o         \
          $(OBJDIR)/time_utilities.o         $(OBJDIR)/utilities.o
	echo compiling src_output.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_output.f90 )
#
$(OBJDIR)/src_radiation.o:                   $(SRCDIR)/src_radiation.f90      \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_flake.o             $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_radiation.o         $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/data_soil.o              $(OBJDIR)/data_turbulence.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/meteo_utilities.o      \
          $(OBJDIR)/utilities.o              $(OBJDIR)/data_tracer.o
	echo compiling src_radiation.f90 
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_radiation.f90 )
#
$(OBJDIR)/src_relaxation.o:                  $(SRCDIR)/src_relaxation.f90     \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_lheat_nudge.o     \
          $(OBJDIR)/data_modelconfig.o                                        \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/environment.o          \
          $(OBJDIR)/meteo_utilities.o        $(OBJDIR)/src_lheating.o         \
          $(OBJDIR)/utilities.o
	echo compiling src_relaxation.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_relaxation.f90 )
#
$(OBJDIR)/src_runge_kutta.o:                 $(SRCDIR)/src_runge_kutta.f90    \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_lheat_nudge.o       $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_hymacs.o            $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/meteo_utilities.o      \
          $(OBJDIR)/numeric_utilities_rk.o   $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/numeric_utilities.o      $(OBJDIR)/time_utilities.o       \
          $(OBJDIR)/hori_diffusion.o         $(OBJDIR)/fast_waves_rk.o        \
          $(OBJDIR)/src_advection_rk.o       $(OBJDIR)/src_lheating.o         \
          $(OBJDIR)/src_slow_tendencies_rk.o $(OBJDIR)/data_tracer.o
	echo compiling src_runge_kutta.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_runge_kutta.f90 )
#
$(OBJDIR)/src_sat_tbs.o:                     $(SRCDIR)/src_sat_tbs.f90        \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_satellites.o                                         \
          $(OBJDIR)/parallel_utilities.o
	echo compiling src_sat_tbs.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_sat_tbs.f90 )
#
$(OBJDIR)/src_sat_rttov.o:                   $(SRCDIR)/src_sat_rttov.f90      \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_satellites.o      \
          $(OBJDIR)/data_fields.o
	echo compiling src_sat_rttov.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_sat_rttov.f90 )
#
$(OBJDIR)/src_seaice.o:                      $(SRCDIR)/src_seaice.f90         \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_soil.o            \
          $(OBJDIR)/data_turbulence.o                                         \
          $(OBJDIR)/meteo_utilities.o
	echo compiling src_seaice.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_seaice.f90 )
#
$(OBJDIR)/src_setup.o:                       $(SRCDIR)/src_setup.f90          \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_convection.o      \
          $(OBJDIR)/data_gscp.o              $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_soil.o            \
          $(OBJDIR)/data_turbulence.o                                         \
          $(OBJDIR)/environment.o            $(OBJDIR)/io_utilities.o         \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/utilities.o            \
          $(OBJDIR)/src_sso.o
	echo compiling src_setup.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_setup.f90 )
#
$(OBJDIR)/src_setup_vartab.o:                $(SRCDIR)/src_setup_vartab.f90   \
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/environment.o                                             \
          $(OBJDIR)/data_tracer.o                                             \
          $(OBJDIR)/data_hymacs.o
	echo compiling src_setup_vartab.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG3) $(SRCDIR)/src_setup_vartab.f90 )
#
$(OBJDIR)/src_sfcana.o:                      $(SRCDIR)/src_sfcana.f90         \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_nudge_all.o         $(OBJDIR)/data_obs_record.o      \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/environment.o          \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/utilities.o
	echo compiling src_sfcana.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_sfcana.f90 )
#
$(OBJDIR)/src_sing_local.o:                  $(SRCDIR)/src_sing_local.f90     \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_nudge_gather.o      $(OBJDIR)/data_nudge_local.o     \
          $(OBJDIR)/data_obs_record.o        $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_turbulence.o
	echo compiling src_sing_local.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_sing_local.f90 )
#
$(OBJDIR)/src_sing_spread.o:                 $(SRCDIR)/src_sing_spread.f90    \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_nudge_gather.o      $(OBJDIR)/data_nudge_spread.o    \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/src_correl_cutoff.o
	echo compiling src_sing_spread.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_sing_spread.f90 )
#
$(OBJDIR)/src_slow_tendencies_rk.o:       $(SRCDIR)/src_slow_tendencies_rk.f90\
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/numeric_utilities_rk.o   $(OBJDIR)/numeric_utilities.o    \
          $(OBJDIR)/data_tracer.o
	echo compiling src_slow_tendencies_rk.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/src_slow_tendencies_rk.f90 )
#
$(OBJDIR)/src_soil.o:                        $(SRCDIR)/src_soil.f90           \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_soil.o            \
          $(OBJDIR)/meteo_utilities.o
	echo compiling src_soil.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_soil.f90 )
#
$(OBJDIR)/src_soil_multlay.o:                $(SRCDIR)/src_soil_multlay.f90   \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                                                 \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_soil.o            \
          $(OBJDIR)/meteo_utilities.o        $(OBJDIR)/environment.o
	echo compiling src_soil_multlay.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_soil_multlay.f90 )
#
$(OBJDIR)/src_spectral_nudging.o:          $(SRCDIR)/src_spectral_nudging.f90 \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_runcontrol.o                                         \
          $(OBJDIR)/environment.o            $(OBJDIR)/parallel_utilities.o
	echo compiling src_spectral_nudging.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_spectral_nudging.f90 )
#
$(OBJDIR)/src_sso.o:                         $(SRCDIR)/src_sso.f90            \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o
	echo compiling src_sso.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_sso.f90 )
#
$(OBJDIR)/src_turbulence.o:                  $(SRCDIR)/src_turbulence.f90     \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_flake.o             $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_turbulence.o
	echo compiling src_turbulence.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_turbulence.f90 )
#
$(OBJDIR)/src_tracer_supply.o:               $(SRCDIR)/src_tracer_supply.f90  \
          $(OBJDIR)/data_parameters.o                                         \
          $(OBJDIR)/data_tracer.o            $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/environment.o            $(OBJDIR)/meteo_utilities.o
	echo compiling src_tracer_supply.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/src_tracer_supply.f90 )
#
$(OBJDIR)/turbulence_interface.o:            $(SRCDIR)/turbulence_interface.f90\
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/data_turbulence.o                                         \
          $(OBJDIR)/src_artifdata.o                                           \
          $(OBJDIR)/turbulence_tran.o        $(OBJDIR)/turbulence_diff.o
	echo compiling turbulence_interface.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/turbulence_interface.f90 )
#
$(OBJDIR)/turbulence_tran.o:                 $(SRCDIR)/turbulence_tran.f90    \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_flake.o           \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/data_turbulence.o                                         \
          $(OBJDIR)/turbulence_utilities.o                                    \
          $(SRCDIR)/statement_functs.incf
	echo compiling turbulence_tran.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/turbulence_tran.f90 )
#
$(OBJDIR)/turbulence_diff.o:                 $(SRCDIR)/turbulence_diff.f90    \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_flake.o           \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/data_turbulence.o                                         \
          $(OBJDIR)/turbulence_utilities.o                                    \
          $(SRCDIR)/statement_functs.incf
	echo compiling turbulence_diff.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/turbulence_diff.f90 )
#
################################################################################
#
#  Dependencies of the external subroutines
#
################################################################################
#
$(OBJDIR)/dfi_initialization.o:              $(SRCDIR)/dfi_initialization.f90 \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                                                 \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/data_soil.o              $(OBJDIR)/data_satellites.o      \
          $(OBJDIR)/src_meanvalues.o         $(OBJDIR)/src_gridpoints.o       \
                                             $(OBJDIR)/environment.o          \
          $(OBJDIR)/meteo_utilities.o        $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/time_utilities.o         $(OBJDIR)/utilities.o
	echo compiling dfi_initialization.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/dfi_initialization.f90 )
#
$(OBJDIR)/fast_waves.o:                      $(SRCDIR)/fast_waves.f90         \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/time_utilities.o
	echo compiling fast_waves.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/fast_waves.f90 )
#
$(OBJDIR)/fast_waves_rk.o:                   $(SRCDIR)/fast_waves_rk.f90      \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/src_artifdata.o                                           \
          $(OBJDIR)/environment.o            $(OBJDIR)/numeric_utilities_rk.o \
          $(OBJDIR)/time_utilities.o
	echo compiling fast_waves_rk.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/fast_waves_rk.f90 )
#
$(OBJDIR)/hori_diffusion.o:                  $(SRCDIR)/hori_diffusion.f90     \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/numeric_utilities.o    \
          $(OBJDIR)/time_utilities.o         $(OBJDIR)/data_tracer.o
	echo compiling hori_diffusion.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG3) $(SRCDIR)/hori_diffusion.f90 )
#
$(OBJDIR)/near_surface.o:                    $(SRCDIR)/near_surface.f90       \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/data_turbulence.o      \
          $(OBJDIR)/meteo_utilities.o
	echo compiling near_surface.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/near_surface.f90 )
#
$(OBJDIR)/organize_assimilation.o:        $(SRCDIR)/organize_assimilation.f90 \
          $(OBJDIR)/data_lheat_nudge.o                                        \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/environment.o          \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/src_lheat_nudge.o      \
          $(OBJDIR)/src_nudging.o            $(OBJDIR)/src_sfcana.o
	echo compiling organize_assimilation.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/organize_assimilation.f90 )
#
$(OBJDIR)/organize_data.o:                   $(SRCDIR)/organize_data.f90      \
                                             $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_nudge_all.o       \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o                                         \
          $(OBJDIR)/environment.o            $(OBJDIR)/io_utilities.o         \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/src_artifdata.o        \
          $(OBJDIR)/src_artifdata.o          $(OBJDIR)/src_input.o            \
          $(OBJDIR)/src_output.o             $(OBJDIR)/src_setup_vartab.o     \
          $(OBJDIR)/mpe_io.o
	echo compiling organize_data.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/organize_data.f90 )
#
$(OBJDIR)/organize_diagnosis.o:              $(SRCDIR)/organize_diagnosis.f90 \
          $(OBJDIR)/data_io.o                $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o        $(OBJDIR)/environment.o          \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/utilities.o            \
          $(OBJDIR)/src_diagbudget.o         $(OBJDIR)/src_gridpoints.o       \
          $(OBJDIR)/src_integrals.o          $(OBJDIR)/src_meanvalues.o
	echo compiling organize_diagnosis.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/organize_diagnosis.f90 )
#
$(OBJDIR)/organize_dynamics.o:               $(SRCDIR)/organize_dynamics.f90  \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_tracer.o            $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parallel.o          $(OBJDIR)/data_parameters.o      \
          $(OBJDIR)/data_runcontrol.o                                         \
          $(OBJDIR)/environment.o            $(OBJDIR)/numeric_utilities.o    \
          $(OBJDIR)/parallel_utilities.o     $(OBJDIR)/utilities.o            \
          $(OBJDIR)/src_2timelevel.o         $(OBJDIR)/src_leapfrog.o         \
          $(OBJDIR)/src_relaxation.o         $(OBJDIR)/src_runge_kutta.o      \
          $(OBJDIR)/src_spectral_nudging.o   $(OBJDIR)/src_advection_rk.o     \
          $(OBJDIR)/fast_waves_rk.o
	echo compiling organize_dynamics.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/organize_dynamics.f90 )
#
$(OBJDIR)/organize_eps.o:                    $(SRCDIR)/organize_eps.f90       \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_runcontrol.o                                         \
          $(OBJDIR)/parallel_utilities.o
	echo compiling organize_eps.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/organize_eps.f90 )
#
$(OBJDIR)/organize_physics.o:                $(SRCDIR)/organize_physics.f90   \
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/src_artifdata.o          $(OBJDIR)/data_tracer.o          \
          $(OBJDIR)/src_conv_kainfri.o       $(OBJDIR)/src_conv_tiedtke.o     \
          $(OBJDIR)/src_conv_shallow.o       $(OBJDIR)/src_flake.o            \
          $(OBJDIR)/src_conv_hymacs.o        $(OBJDIR)/data_hymacs.o          \
          $(OBJDIR)/src_gscp.o               $(OBJDIR)/src_radiation.o        \
          $(OBJDIR)/src_soil_multlay.o       $(OBJDIR)/src_soil.o             \
          $(OBJDIR)/src_sso.o                $(OBJDIR)/src_seaice.o           \
          $(OBJDIR)/src_turbulence.o         $(OBJDIR)/turbulence_interface.o \
          $(OBJDIR)/time_utilities.o
	echo compiling organize_physics.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/organize_physics.f90 )
#
$(OBJDIR)/organize_satellites.o:             $(SRCDIR)/organize_satellites.f90\
          $(OBJDIR)/data_fields.o            $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/data_satellites.o        $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/src_sat_tbs.o            $(OBJDIR)/src_sat_rttov.o
	echo compiling organize_satellites.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/organize_satellites.f90 )
#
$(OBJDIR)/semi_implicit.o:                   $(SRCDIR)/semi_implicit.f90      \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/environment.o            $(OBJDIR)/parallel_utilities.o   \
          $(OBJDIR)/time_utilities.o
	echo compiling semi_implicit.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/semi_implicit.f90 )
#
$(OBJDIR)/slow_tendencies.o:                 $(SRCDIR)/slow_tendencies.f90    \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_flake.o             $(OBJDIR)/data_modelconfig.o     \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/data_turbulence.o
	echo compiling slow_tendencies.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG2) $(SRCDIR)/slow_tendencies.f90 )
#
################################################################################
#
#  Dependencies of OASIS utility continued 
#
################################################################################
#
$(OBJDIR)/oas_cos_define.o:                  $(OASDIR)/oas_cos_define.F90     \
	$(OBJDIR)/data_parameters.o        $(OBJDIR)/oas_cos_vardef.o       \
	$(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_fields.o          \
	$(OBJDIR)/data_parallel.o          $(OBJDIR)/parallel_utilities.o
	echo compiling oas_cos_define.F90
	( cd $(OBJDIR) && $(F90) $(COMFLG5)  $(OASDIR)/oas_cos_define.F90 )

$(OBJDIR)/oas_cos_init.o:                    $(OASDIR)/oas_cos_init.F90       \
	$(OBJDIR)/oas_cos_vardef.o
	echo compiling oas_cos_init.F90
	( cd $(OBJDIR) && $(F90) $(COMFLG5)  $(OASDIR)/oas_cos_init.F90 )

$(OBJDIR)/oas_cos_finalize.o:                $(OASDIR)/oas_cos_finalize.F90   \
	$(OBJDIR)/oas_cos_vardef.o
	echo compiling oas_cos_finalize.F90
	( cd $(OBJDIR) && $(F90) $(COMFLG5)  $(OASDIR)/oas_cos_finalize.F90 )

$(OBJDIR)/oas_cos_rcv.o:                    $(OASDIR)/oas_cos_rcv.F90         \
	$(OBJDIR)/data_parameters.o       $(OBJDIR)/oas_cos_vardef.o        \
	$(OBJDIR)/data_modelconfig.o
	echo compiling oas_cos_rcv.F90
	( cd $(OBJDIR) && $(F90) $(COMFLG5) $(OASDIR)/oas_cos_rcv.F90 )

$(OBJDIR)/oas_cos_snd.o:                    $(OASDIR)/oas_cos_snd.F90         \
	$(OBJDIR)/data_parameters.o       $(OBJDIR)/oas_cos_vardef.o        \
	$(OBJDIR)/data_modelconfig.o
	echo compiling oas_cos_snd.F90
	( cd $(OBJDIR) && $(F90) $(COMFLG5) $(OASDIR)/oas_cos_snd.F90 )

$(OBJDIR)/receive_fld_2clm.o:               $(OASDIR)/receive_fld_2clm.F90    \
	$(OBJDIR)/data_modelconfig.o      $(OBJDIR)/oas_cos_vardef.o        \
	$(OBJDIR)/data_parameters.o       $(OBJDIR)/data_parallel.o         \
	$(OBJDIR)/data_runcontrol.o       $(OBJDIR)/data_constants.o        \
	$(OBJDIR)/data_soil.o             $(OBJDIR)/data_fields.o
	echo compiling receive_fld_2clm.F90 
	( cd $(OBJDIR) && $(F90) $(COMFLG5) $(OASDIR)/receive_fld_2clm.F90 )

$(OBJDIR)/send_fld_2clm.o:                  $(OASDIR)/send_fld_2clm.F90       \
	$(OBJDIR)/data_modelconfig.o      $(OBJDIR)/oas_cos_vardef.o        \
	$(OBJDIR)/data_parameters.o       $(OBJDIR)/data_parallel.o         \
	$(OBJDIR)/data_runcontrol.o       $(OBJDIR)/data_fields.o
	echo compiling send_fld_2clm.F90 
	( cd $(OBJDIR) && $(F90) $(COMFLG5) $(OASDIR)/send_fld_2clm.F90 )
#
################################################################################
#
#  Dependencies of the Main Program
#
################################################################################
#
$(OBJDIR)/lmorg.o:                           $(SRCDIR)/lmorg.f90              \
          $(OBJDIR)/data_constants.o         $(OBJDIR)/data_fields.o          \
          $(OBJDIR)/data_flake.o             $(OBJDIR)/data_io.o              \
          $(OBJDIR)/data_modelconfig.o       $(OBJDIR)/data_parallel.o        \
          $(OBJDIR)/data_parameters.o        $(OBJDIR)/data_runcontrol.o      \
          $(OBJDIR)/data_soil.o              $(OBJDIR)/data_tracer.o          \
          $(OBJDIR)/environment.o            $(OBJDIR)/meteo_utilities.o      \
          $(OBJDIR)/src_allocation.o         $(OBJDIR)/src_artifdata.o        \
          $(OBJDIR)/src_setup.o              $(OBJDIR)/src_tracer_supply.o    \
          $(OBJDIR)/time_utilities.o         $(OBJDIR)/utilities.o            \
          $(OBJDIR)/mpe_io.o                 $(OBJDIR)/info_lm_f90.o
	echo compiling lmorg.f90
	( cd $(OBJDIR) && $(F90) $(COMFLG1) $(SRCDIR)/lmorg.f90 )
#
