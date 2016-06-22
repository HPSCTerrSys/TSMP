!+ Source module for the setup of the LM
!------------------------------------------------------------------------------

!!! IF PROBLEMS WITH NAMELISTS OCCUR: COMPILE ON PLATFORMS OTHER THAN NEC,
!!! THEN THE IOMSG-MECHANISM (FORTRAN 2003) WILL GIVE YOU MUCH MORE DETAILED
!!! ERROR MESSAGES TO FIND ERRORS IN NAMELIST FILES. FOR THIS, YOU
!!! NEED FAIRLY NEW VERSIONS OF COMPILERS LIKE INTEL (>= 11.0) 
!!! OR GFORTRAN (> 4.3 ??). 
!!! IF THIS IS THE CASE, JUST UNCOMMENT THE FOLLOWING LINE AND
!!! THE IOMSG-MECHANISM WILL BE ACTIVE FOR THE READING OF NAMELIST
!!! "INPUT_IDEAL".
!!!


!!! IN THE SUBROUTINE SEED_RANDOM_NUMBER() BELOW, THERE IS A CALL TO THE
!!! SYSTEM FUNCTION DATE_AND_TIME(). WITH THE PGI-COMPILER, THE DATA TYPES 
!!! OF THE ARGUMENTS OF THIS FUNCTION SEEM TO NEED THE SPECIFIC DATA TYPE 
!!! INTEGER*8, WHICH IS NOT THE CASE ON OTHER PLATFORMS. THEREFORE,
!!! THE DESTINCTION IS MADE BY HELP OF THE FOLLOWING PREPROCESSOR SWITCH
!!! (IF YOU WORK WITH PGI, UNCOMMENT THE FOLLOWING LINE OR ADD A
!!! "-D__PGI_FORTRAN__" TO YOUR COMPILER OPTIONS)
!!!
!#define __PGI_FORTRAN__
!!!
!!! HOWEVER, THIS STRANGE BEHAVIOR OF THE PGI-COMPILER IS NOT DOCUMENTED
!!! IN THE COMPILER DOCUMENTATION (IT SAYS "STANDARD INTEGER DATA TYPE" 
!!! INSTEAD)! SO YOU HAVE TO EXPERIMENT A LITTLE, IF YOU USE PGI!

!------------------------------------------------------------------------------

!!! NOTE FOR USING SOME PHYSICAL PARAMETERIZATIONS IN IDEALIZED RUNS:
!!!
!!! THIS IS THE LIST OF CURRENTLY USED EXTERNAL PARAMETERS IN THE COSMO MODEL:
!!! NOT ALL OF THEM ARE CURRENTLY IMPLEMENTED IN SRC_ARTIFDATA.F90!
!!! THOSE WITH A * ARE MISSING UP TO NOW, THE OTHERS MAY BE GIVEN A CONSTANT
!!! VALUE OR READ FROM AN ASCII FILE (2D).:
!!!
!!$  In any case needed ??? (may not hurt ...):
!!$
!!$      z0         ,    & ! surface roughness                             (  m  )
!!$      fr_land    ,    & ! fraction of land in a grid element            ( --  )
!!$      plcov      ,    & ! fraction of plant cover                         --
!!$      lai        ,    & ! leaf area index of plants                       --
!!$      soiltyp    ,    & ! type of the soil (keys 1-10)                  ( --  )
!!$
!!$    if (lsso)
!!$      sso_stdh   ,    & ! standard deviation of sub-grid scale orography( m   )
!!$      sso_gamma  ,    & ! anisotropy of sub-grid scale orography          --
!!$      sso_theta  ,    & ! angle betw. principal axis of orography and E ( rad )
!!$      sso_sigma  ,    & ! mean slope of sub-grid scale orography          --
!!$
!!$    if (lsoil)
!!$      rootdp     ,    & ! depth of the roots                            (  m  )
!!$
!!$    if (lstomata)
!!$ *    rsmin2d    ,    & ! minimum stomata resistance                    ( s/m )
!!$
!!$    if (lrad)
!!$      vio3       ,    & ! vertical integrated ozone contents            (pa O3)
!!$      hmo3       ,    & ! ozone maximum                                 ( pa  )
!!$
!!$    if (lrad .and. itype_aerosol == 2)
!!$ *    aer_su     ,    & ! monthly aerosol climatology sulfate drops     (0 - 1)
!!$ *    aer_du     ,    & ! monthly aerosol climatology total dust        (0 - 1)
!!$ *    aer_or     ,    & ! monthly aerosol climatology organic (water sol.)(0-1)
!!$ *    aer_bc     ,    & ! monthly aerosol climatology black carbon      (0 - 1)
!!$ *    aer_ss     ,    & ! monthly aerosol climatology sea salt          (0 - 1)
!!$
!!$    if (lrad .and. lemiss)
!!$ *    emis_rad   ,    & ! external thermal emissivity                   (0 - 1)
!!$
!!$    if (lrad .and. lradtopo)
!!$ *    skyview    ,    & ! sky view
!!$ *    slo_asp    ,    & ! slope aspect
!!$ *    slo_ang    ,    & ! slope angle
!!$ *    horizon    ,    & ! horizon
!!$
!!$    if (lforest)
!!$      for_e      ,    & ! ground fraction covered by evergreen forest     --
!!$      for_d      ,    & ! ground fraction covered by deciduous forest     --
!!$
!!$    if (llake)     !!! FLAKE OPTION NOT CORRECTLY INITIALIZED YET IN SRC_ARTIFDATA.F90! DO NOT USE IN IDEALIZED RUNS!
!!$ *    t_mnw_lk  ,     & ! mean temperature of the water column          (  K  )
!!$ *    t_wml_lk  ,     & ! mixed-layer temperature                       (  K  )
!!$ *    t_bot_lk  ,     & ! temperature at the water-bottom sediment
!!$                        ! interface                                     (  K  )
!!$ *    t_b1_lk   ,     & ! temperature at the bottom of the upper layer
!!$                        ! of the sediments                              (  K  )
!!$ *    c_t_lk    ,     & ! shape factor with respect to the
!!$                        ! temperature profile in lake thermocline       (  -  )
!!$ *    h_ml_lk   ,     & ! thickness of the mixed-layer                  (  m  )
!!$ *    h_b1_lk           ! thickness of the upper layer
!!$                        ! of bottom sediments                           (  m  )
!!$ *    fr_lake    ,    & ! lake fraction in a grid element [0,1]         (  -  )
!!$ *    depth_lk   ,    & ! lake depth (SET TO -1.0 AT NON-LAKE POINTS!)  (  m  )
!!$ *    fetch_lk   ,    & ! wind fetch over lake                          (  m  )
!!$ *    dp_bs_lk   ,    & ! thickness of the thermally active layer
!!$                        ! of bottom sediments                           (  m  )
!!$ *    t_bs_lk    ,    & ! climatological temperature at the bottom of
!!$                        ! the thermally active layer of sediments       (  K  )
!!$ *    gamso_lk          ! attenuation coefficient for
!!$                        ! solar radiation in lake water                 ( 1/m )
!!$
!!$    if (lseaice .or. llake)
!!$      h_ice     ,     & ! ice thickness                                 (  m  )
!!$      t_ice     ,     & ! temperature at the snow-ice or air-ice interface (  K  )
!!$
!!$    if (lseaice)    !!! NOT IMPLEMENTED YET !!! BOX EITHER ICE-FREE (H_ICE = 0) OR TOTALLY ICE COVERED (H_ICE > 0) !!!
!!$ *    fr_ice          & ! ice fraction for ocean/lake surfaces          (  -  )
!!$
!!$
!!$  IN ADDITION, THESE MORE VARYING THERMODYNAMIC SOIL AND SNOW VARIABLES 
!!$  ARE ALREADY INITIALIZED BELOW (EITHER AS CONSTANT VALUES OR READING
!!$  FROM A 2D ASCII FILE):
!!$  
!!$       t_snow    ,     & ! temperature of the snow-surface               (  k  )
!!$       t_snow_mult,    & ! temperature of the snow-surface               (  k  )
!!$       t_s       ,     & ! temperature of the soil, also surface temperature at water and ice points (  k  )
!!$       t_g       ,     & ! weighted surface temperature                  (  k  )
!!$       qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)
!!$       t_m       ,     & ! temperature between upper and medium 
!!$                         ! soil layer                                    (  k  )
!!$       t_cl      ,     & ! temperature beetween medium and lower
!!$                         ! soil layer                                    (  k  )
!!$       t_so      ,     & ! multi-layer soil temperature                  (  k  )
!!$       w_snow    ,     & ! snow water equivalent                         (m H2O)
!!$       w_i       ,     & ! water content of interception water           (m H2O)
!!$       w_g1      ,     & ! water content of the upper soil layer         (m H2O)
!!$       w_g2      ,     & ! water content of the medium soil layer        (m H2O)
!!$       w_g3      ,     & ! water content of the lower soil layer         (m H2O)
!!$                         ! (if NLWB=3, unused otherwise)
!!$       w_so      ,     & ! multi-layer soil moisture                  (m H2O)
!!$       w_so_ice  ,     & ! multi-layer soil ice                       (m H2O)
!!$       w_cl      ,     & ! climatological water content                  (m H2O) 
!!$       freshsnow ,     & ! weighting function indicating 'freshness' of snow
!!$


!------------------------------------------------------------------------------

!!! Debug mode: set below ldebug_artif = .true. and idbg_level to some value > 0.
!!!
!!! Currently implemented messages/mechanisms depending on idbg_level:
!!!   idbg_artif_level > 0 : print the subroutine name at the beginning of each subroutine
!!!   idbg_artif_level > 3 : additionally, write ASCII-files (or BIN-files on the NEC) containing
!!!                          the T- and QV- increment resp. heating rate for each artif. 
!!!                          temperature-, moisture or heating rate disturbance
!!!                          triggered in the simulation.
!!!                     *** On the NEC, you need the conversion program 
!!!                         "bin2ascii_convrates3d.f90" by Ulrich Blahak to
!!!                          generate ASCII-files from the BIN-files ***
!!!   idbg_artif_level > 4 : print additional checking output for the iterative hydrostatic
!!!                          pressure integration

MODULE src_artifdata


!------------------------------------------------------------------------------
!
! Description:
!
!   This module contains subroutines for the necessary "ingredients" to conduct
!   idealized simulations with the COSMO-Model:
!
!   - vertical coordinate specification
!   - specification of the reference atmosphere
!   - orography (read from ASCII-file or "analytical" hills)
!   - time-constant surface parameters (read from ASCII-file
!     or constant values for the whole domain)
!   - artificial inital and boundary data for the model:
!     vertical profiles from an ASCII-file or use analytical formulas
!   - enforce initial boundary layer for velocity near the ground,
!     based on an exponent wind profile
!   - hydrostatic pressure initialization
!   - free-slip lower BC for momentum and/or heat
!   - artificial convection triggers:
!     generate triggers in the initial data or within a time period later
!     during the model run.
!   - time-constant surface sensible and latent heat fluxes,
!     with possible added time-constant spatial white noise
!   - initial spatial noise in the lowest 100 hPa on w and T
!
!   The implementation is such that the different
!   ingredients for idealized simulations may be freely combined by
!   various namelist switches within the namelist INPUT_IDEAL.
!
!   There is an example runscript in the subdirectory *../run_ideal*
!   named *run_ideal*, which contains all namelist switches
!   for idealized runs in a documented fashion and which
!   may be used as a "cookbook" to set up your own idealized
!   runs. 
!
!   There is a separate PDF documentation *artif_docu.pdf*.
!
!   The user may also modify the below subroutines for his/her
!   own purpose if the idealized setup in mind cannot be
!   set up by the already implemented ingredients.
!
!
! Method:
!   See subroutines below
!
! Current Code Owner: DWD, Ulrich Blahak
!  phone:  +49  69  8062 2393
!  fax:    +49  69  8062 3721
!  email:  ulrich.blahak@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Guenter Doms
!  Initial release
! 1.4        1998/05/22 Guenther Doms
! Adaptions for the two-time level time integration scheme
! 1.5        1998/06/29 Guenther Doms
!  Printout of vertical coordinates now in routine 'print_vertcoord'
!  Adaptions to the call of routine 'reference_atmosphere' and 'gen_ini_data'
! 1.7        1998/07/16 Guenther Doms
!  Change in the calling list of routine 'calps'.
! 1.8        1998/08/03 Ulrich Schaettler
!  Correct ANSI violations.
! 1.10       1998/09/29 Ulrich Schaettler
!  Adapted call to routine reference_atmosphere.
! 1.20       1999/01/07 Guenhter Doms
!  Renaming of some global variables
! 1.29       1999/05/11 Ulrich Schaettler
!  Bug fixing
! 1.32       1999/08/24 Guenther Doms
!  Corrections for 2D runs and periodic boundary conditions.
! 1.34       1999/12/10 Ulrich Schaettler
!  Renaming of some boundary variables and changed calls to utility-routines
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed names for variables concerned to latitude or longitude and use 
!  of variables klvxxx now from data_modelconfig.
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists
! 2.14       2002/02/15 Ulrich Schaettler
!  Modifications to allow the use of the SLEVE coordinate
! 2.18       2002/07/16 Reinhold Schrodin
!  Eliminated variable rhde, use cf_snow instead
! 3.5        2003/09/02 Ulrich Schaettler
!  Adapted interface for routine exchg_boundaries
! 3.7        2004/02/18 Ulrich Schaettler
!  Adapted dimension of kzdims for routine exchg_boundaries
! 3.13       2004/12/03 Ulrich Schaettler
!  Adapted SR-name dc_topo => sleve_split_oro;
!  Adapted interfaces of gather_field and reference_atmosphere
! 3.14       2005/01/25 Ulrich Schaettler
!  Adjusted computation of llandmask: fr_land=0.5 is a land point now!!!
!  Editorial changes
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.18       2006/03/03 Jochen Foerstner
!  Changed treatment of w for lateral boundary relaxation
! 3.21       2006/12/04 Ulrich Schaettler
!  Included klv950, klv700 in interface to reference_atmosphere
! V3_24        2007/04/26 Michael Baldauf
!  Use lmulti_layer to determine which fields have to be initialized
! V4_1         2007/12/04 Ulrich Schaettler
!  Call to SR sleve_split_oro: introduced my_cart_id as argument
! V4_5         2008/09/10 Guenther Zaengl
!  Add new namelist for idealized runs (called only if lartif_data = .true.)
!  Include option for new reference atmosphere
!  More accurate initialization of perturbation pressure
! V4_8         2009/02/16 Guenther Zaengl
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation
!  Renamed IDEALCTL to ARTIFCTL (similar to GME)
! V4_11        2009/11/30 Ekaterina Machulskaya
!  Adaptations to use multi-layer snow model
! V4_12        2010/05/11 Guenther Zaengl, Ulrich Schaettler
!  Adaptations for reference atmosphere with constant BruntVaisala frequency
!  Exact initialization of perturbation pressure for itheta_adv=0
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_17        2011/02/24 Ulrich Blahak
!  Total reorganization of the file; many new namelist parameters;
!  documentation of the new namelist parameters in the exemplary namelist file *run_ideal*
!  and in a separate PDF-Document *artif_docu.pdf*;
!  warm bubbles and artificial heating rates are now imposed in lmorg.f90; 
!  namelist ARTIFCTL is now read in organize_data() after the gribout-namelist(s),
!  regardless of the setting of lartif_data; this is intended to
!  enable the specification of "warm bubbles", no-surface-flux-conditions and 
!  fixed turb. diffusion coefficients in real case simulations (DISABLED FOR NOW ---
!  IF INTENDED TO USE, CONTACT Ulrich Blahak). In this
!  case, not all the namelist parameters take effect, but only these:
!    - parameters for warm bubbles / temperature disturbances
!    - parameters for itype_turb = 100: tkvhfix, tkhhfix, tkvmfix, tkhmfix
!      (fixed diffusion coefficients)
!    - the switches "lnosurffluxes_m" and "lnosurffluxes_h" to enable runs without surface fluxes
!      (free-slip-condition and/or no surface heat/moisture fluxes).
!  Adapted interface of exchg_boundaries; corrected kzdims(1:20) -> kzdims(1:24);
!  eliminated my_peri_neigh.
! V4_18        2011/05/26 Ulrich Blahak
!  All comments are now in English.
!  Bugfix in the computation of the qvtens in SR artif_heatrate_dist()
!    (here, still the old approximate formula was used instead of the exact function
!    dqvdT_prh(), which is now called instead).
!  Introduced global registry strings for disturbance types (convection triggers),
!    to ease the implementation of own disturbance types by requiring their
!    "registration" at only one place in the code instead of 2 places before.
!  Introduced consistency checks of parameters for analytic T/QV profiles
!    to avoid cases where the "analytic" atmosphere has a finite height
!    (e.g., the polytrope atmosphere has a finite height!) and the
!    model top height is specified higher.
!  Corrected sign of shift_i, shift_j in read_ascii_field_2d().
!  Changed the meaning of the hill_width parameters from 1/2-width diameter
!    to 1/2-width radius (now half the value gives the same moutain
!    width as before). This is more consistent with the rationale for 
!    characteristic lengths elsewhere, e.g. for "warm bubble" radii.
!  Take into account the case N_nconst(k) == 0.0 when checking for
!    zz_top beeing smaller than the physical height of the N=const atmosphere
!    layer(s).
!  Added new namelist parameters t_soil_c, t_snow_c, wf_soil_c, w_i_c, w_snow_c
!    to initialize soil temperature (the same in all depths), 
!    snow temperature (the same in all depths), soil water 
!    saturation (the same in all depths),
!    interception water and snow water equivalent with constant values
!    different from t_s resp. 0.0, throughout the respective column. 
!    If these parameters are not specified, the "old" behaviour (t_s and 0.0)
!    is obtained, where t_s represents the atmospheric near ground temperature.
!    Also, if t_soil_c < 0.0 then t_soil=t_s, and if
!    t_snow_c < 0.0 then t_snow=t_s.
!  Renamed namelist parameters for constant soil- and atmosphere fields
!    from "*_ideal" to "*_c" for consistency.
!  Renamed "itype_soil" to "itype_soil_c" and introduced "itype_soil_tw".
!    If now "itype_soil_tw = 1", the fields t_soil=t_soil_c, 
!      wf_soil=wf_soil_c, w_i=w_i_c, w_snow=w_snow_c and t_snow=t_snow_c;
!    if "itype_soil_tw = 2", t_soil, wf_soil, w_snow, t_snow and w_i
!      are read from 2D ASCII files. 
!    The names of these ASCII files must be given in the new namelist
!      parameters tsoilfile, tsnowfile, wfsoilfile, wsnowfile, wifile.
!  Created the small F90-program "gen_soildata_ascii.f90" to produce
!    a sample ASCII data set suitable for "itype_soil_c = 2" and itype_soil_tw = 2".
!  Eliminated alb_rad and albradfile, since alb_rad is not an input parameter,
!    but calculated in the model depending on soil type, soil wetness, 
!    plant cover, snow, ...
!  Bugfix: ie, je --> ni, nj in SR calc_p_hydrostat_ana, calc_p_hydrostat_psts 
!    and calc_p_hydrostat_lf.
!  Deleted some unused variables and did some general code beautifying.
!  Default values changed: if no namelist parameters are specified (empty
!    namelist) the ICAO standard atmosphere with U=20 m/s everywhere results.
!  Bugfix: fixed error in add_orography_ana() with regard to the use
!    of zhillcutfact and a negative hillheight (valley):
!    variable "zhillcut" must be positive, i.e., "ABS(hillheight(ii))*zhillcutfact(ii)"
!    instead of "hillheight(ii)*zhillcutfact(ii)".
!  Implemented initialization of h_ice and t_ice for lseaice=.true. .or. llake=.true.
!  Implemented initialization of parameters for the SSO-scheme (lsso=.true.)
!  Implemented correct spherical metric into artificial hills/bubbles, if
!    lmetr=.true. In this case, the main axes of hills/bubbles are now
!    great circles, and x- and y- distances are measured on great circles
!    perpendicular to these main axes.
! V4_20        2011/08/31 Ulrich Blahak
!  On platforms other than the NEC SX, introduced IOMSG-identifier in OPEN and
!    READ of the namelist ARTIFCTL for much more detailed error messages. Sadly,
!    this feature is not available on the NEC, so it is capsuled by "#ifdef HAS_IOMSG".
!  Introduced possibility to specify pot. temp. profiles in radiosonde file
!    instead of ordinary temperature. Modified SR read_raso(), 
!    added new NL parameter "rasofile_t_is_theta".
!  Bugfix in SR read_ascii_field_2d(): wrong error code for runs with
!    > 1 PE leads to spurious model abort.
!  Bugfixes in the subroutines for heat/moisture disturbances: allocatable
!    "bub_tnoisebuffer" was used when it was not allocated in case of
!    "ladd_bubblenoise(ii)=.false."
! V4_21        2011/12/06 Ulrich Blahak
!  Removed debug output in inner loops of calc_p_hydrostat_xxxx()-routines
!    for better vectorization. Changed calls to distribute_values(...,imp_character,...)
!    for vector of character strings in input_artifctl().
!  Added validity checks for constant parameters gz0, fr_land, soiltyp,
!    plcov, lai, rootdp.
! V4_23        2012/05/10 Ulrich Blahak, Ulrich Schaettler
!  Eliminated unnecessary and potentially dangerous re-calculation
!   of rvd_m_o in the routines for hydrostatic pressure calculation
!   (calc_p_hydrostat_XXX).
!  Introduced call to init_grid_metrics (Michael Baldauf)
!  Adapted usage of obsolete features in WRITE Statements (Oli Fuhrer)
!  Adapted call to SR distribute_fields (added sender PE) (Uli Schaettler)
! V4_24        2012/06/22 Michael Baldauf
!  Adapted calls to SR reference_atmosphere_x
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!  Replaced qx-variables by using them from the tracer module
!  UB:
!  Added new namelist parameters "iseed_noise_t" and "iseed_noise_w":
!   these are the seeds for the random number generators for T- and W-noise
!   in the lowest 100 hPa, if ladd_noise_t/w = .true. and are simple
!   integer scalars. Before, these seeds were hardwired in the code (Uli Blahak)
!  Bugfix: Replaced istartpar,iendpar,... in tgcom() with 1,ie,..., otherwise
!    t_g is 0.0 in the interior boundary lines and model may crash in organize_radiation.
!  Bugfix: replaced istartpar,iendpar,... in init_w_followeta() with 1,ie,... so
!    that w is also set on the interior boundary lines. Additionally added
!    boundary exchange of w (necessary for periodic BCs).
!  Added support for new fast waves solver in calc_p_hydrostat_psts().
!  Check soiltyp with NINT-function
!  Removed dependency from ntstep for exchg_boundaries calls (Uli Schaettler)
! V4_26        2012/12/06 Anne Roches
!  Renaming of T_CLP_POSDEF to T_CLP_ON since only on and off are available
!  for the moment. (AR)
! V4_27        2013/03/19 Michael Baldauf, Astrid Kerkweg, Ulrich Schaettler
!  Introduced p0ref in argument list to SR reference_atmosphere_BV
!  MESSy interface introduced (AK)
! V4_28        2013/07/12 Ulrich Schaettler, Oliver Fuhrer
!  Use subroutines and variables for vertical grid and reference atmospheres 
!    from module vgrid_refatm_utils (US)
!  Clarified example given for usage of COSMO-Tracer construct (OF)
! V4_29        2013/10/04 Ulrich Blahak, Astrid Kerkweg, Ulrich Schaettler
!  Corrections for proper use of the new vertical coordinate and reference
!   atmospheres data types
!  Unification of MESSy interfaces and COSMO Tracer structure
!  Fixed computation of seed(i) in function seed_random_number(), because
!   the old method lead to a floating point exception with the gfortran compiler
!  Initialization of t_so, w_so, t_snow_mult, t_snow now only for land points
!   (SR init_tw_soil_snow_c()).
!  New namelist parameter t_surf_c for baseline initialization of surface temperature t_s.
!  Missing initialization of t_cl(nnow) for leapfrog and lbdclim
!  New subroutine set_idealized_surffluxes() for specifying time-constant surface fluxes;
!   corresponding new namelist parameters: lsensiflux_fix, llatentflux_fix, sensiflux_c,
!   latentflux_LzuS, H0_rel_noise, iseed_noise_H0
!  Namelist parameters for noise seeds (iseed_*) can now be set to -999, which
!   causes the system time to be used instead of a fixed seed.
!  Relaxed ASCII format of rasofiles and orography/soil initialization files:
!   Now arbitrary number of header lines possible, starting with '#' or '!'.
!   See format description in the example runscripts (folder RUNSCRIPTS/ in this
!   COSMO-model distribution!
! V4_30        2013/11/08 Ulrich Schaettler
!  Initialized izerror before using in several subroutines
! V5_1         2014-11-28 Ulrich Schaettler, Ulrich Blahak, Oliver Fuhrer
!  Adapted interfaces of reference_atmosphere-routines to the modifications
!  in vgrid_refatm_utils.
!  Dummy allocation of bub_glob_tnoisebuffer in gen_bubnoise() on nodes other than 0. Is
!   necessary at least for gfortran and OpenMPI. (UB)
!  Changed the format of some YUSPECIF entries for the CLM namelist tool.
!  Removed HAS_IOMSG, because this is implemented now throughout the model
!  Replaced ireals by wp (working precision) (OF)
!  Revised initialization of surface temperature for water and ice points
!   (set correct dummy values for t_snow / t_snow_multlay at water points; 
!    added t_water_c; tsurffile now contains the surface temperature for all surface/soil types
!    at all gridpoints, including water and ice; t_s, t_soil and t_ice is initialized correctly
!    from these variables / files, depending on soil type, water or ice surface).
!    Improved cross checks of the resulting surface- and soil parameters for internal model consistency.
!  Improved error handling during iteration for pressure initialization (calc_p_hydrostat_xxx)
!  Eliminated artificial convection triggers in the initial condition in case of restart runs
!  Initialization of hhl_prof
!  Return correct error status (ierrstat and not iz_err) in input_artifctl()
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

  USE data_parameters, ONLY :   &
       wp,        & ! KIND-type parameter for real variables
       irealgrib, & ! KIND-type parameter for standard real variables
       iintegers    ! KIND-type parameter for standard integer variables

!------------------------------------------------------------------------------


  USE data_modelconfig, ONLY :   &

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------
       ie_tot,       & ! number of grid points in zonal direction
       je_tot,       & ! number of grid points in meridional direction
       nlandpoints_tot,  & ! number of land points in the grid
       ie,           & ! number of grid points in zonal direction
       je,           & ! number of grid points in meridional direction
       ke,           & ! number of grid points in vertical direction
       ke_soil,      & ! number of layers in the multi-layer soil model
       ke_snow,      & ! number of layers in the multi-layer snow model
       ke1,          & ! KE+1
       nlandpoints,  & ! number of land points in the grid

! 2a. Variables for the new multi-layer soil model
! --------------------------------------------------------------------

       czmls,        & ! depth of the main soil layers in meters
       czhls,        & ! depth of the half-layers CPS
! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-B-grid.
!    
!   zonal direction
       istart,       & ! start index for the forecast of w, t, qv, qc and pp
       iend,         & ! end index for the forecast of w, t, qv, qc and pp

!   meridional direction
       jstart,       & ! start index for the forecast of w, t, qv, qc and pp
       jend,         & ! end index for the forecast of w, t, qv, qc and pp
       jstartpar,    & ! start index for computations in the parallel program
       jendpar         ! end index for computations in the parallel program

  USE data_modelconfig, ONLY :   &

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------
       pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
       pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
       polgam,       & ! angle between the north poles of the systems
       dlon,         & ! grid point distance in zonal direction (in degrees)
       dlat,         & ! grid point distance in meridional direction (in degrees)
       dt,           & ! large time step of the time integration
       startlon_tot, & ! transformed longitude of the lower left grid point
                                ! of the total domain (in degrees, E>0)
       startlat_tot, & ! transformed latitude of the lower left grid point
       startlon,     & ! transformed longitude of the lower left grid point
                                ! of this subdomain (in degrees, E>0)
       startlat,     & ! transformed latitude of the lower left grid point
                                ! of this subdomain (in degrees, N>0)
       degrad,       & ! factor for transforming degree to rad (pi/180)
       raddeg,       & ! factor for transforming rad to degree (180/pi)
       hhl_prof,     & ! a special hhl-profile


! 7. Layer index corresponding to a specified pressure (xxx HPa)
! ----------------------------------------------------
       klv950, klv850, klv800, klv700, klv500, klv400, klv300,  &
       lalloc_w_g3, lalloc_t_cl, lalloc_w_g3_bd, lalloc_t_s_bd, &

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------
       idt_qv, idt_qc,  idt_qi,  idt_qs,  idt_qg,  idt_qr,  idt_qh,    &
               idt_qnc, idt_qni, idt_qns, idt_qng, idt_qnr, idt_qnh

  ! end of data_modelconfig

  !------------------------------------------------------------------------------

  USE data_constants  , ONLY :   &

! 1. mathematical constants
! -------------------------
       pi,           & ! circle constant

! 2. physical constants and related variables
! -------------------------------------------
       t0_melt,      & ! 273.15 K
       r_d,          & ! gas constant for dry air
       rdv,          & ! r_d / r_v
       o_m_rdv,      & ! 1 - r_d/r_v
       rvd_m_o,      & ! r_v/r_d - 1
       lh_v,         & ! latent heat of vapourization
       cp_d,         & ! specific heat of dry air at constant pressure
       cpdr,         & ! 1 / cp_d
       rcpv,         & ! cp_v / cp_d - 1
       rcpl,         & ! cp_l / cp_d - 1
       rdocp,        & ! r_d / cp_d
       gamma,        & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
       g,            & ! acceleration due to gravity
       gq,           & ! g * g
       gh,           & ! g / 2
       r_earth,      & ! mean radius of the earth
       rho_w,        & ! density of water
       solc,         & ! solar constant

! 3. constants for parametrizations
! ---------------------------------
       p0ref,        & ! reference pressure for Exner-function (Pa)
       b1,           & ! variables for computing the saturation vapour pressure
       b2w,          & ! over water (w) and ice (i)
       b2i,          & ! over water (w) and ice (i)
       b3,           & !               -- " --
       b4w,          & !               -- " --
       b4i,          & !              -- " --
       b234w           ! b2w * (b3 - b4w)

! end of data_constants
!------------------------------------------------------------------------------

  USE data_runcontrol , ONLY :   &
       nuspecif,     & ! output unit for protocolling the task
       nstart,       & ! first time step of the forecast
       nstop,        & ! last time step of the forecast
       ntstep,       & ! actual time step
       nold,         & ! corresponds to ntstep - 1
       nnow,         & ! corresponds to ntstep
       nnew,         & ! corresponds to ntstep + 1
       nbd1,         & ! indices for permutation of the
       nbd2,         & ! two boundary time levels
       lartif_data,  & ! forecast with self-defined artificial data
       nbl_exchg,    & ! number of boundlines to exchange       
       lperi_x,        & ! if lartif_data=.TRUE.: periodic boundary conditions (.TRUE.) in x-dir.
       lperi_y,        & ! if lartif_data=.TRUE.: periodic boundary conditions (.TRUE.) in y-dir.
       lmetr,        & ! lartif_data=.TRUE.:  with metric terms
       l2dim,        & ! 2dimensional version
       l2tls,        & ! forecast with 2-TL integration scheme
       lgsp,         & ! forecast with gridscale precipitation
       itype_gscp,   & ! type of microphys. parametrization
       itype_turb,   & ! type of turbulent diffusion parametrization
       itype_tran,   & ! type of surface layer transport parameterization
       itype_fast_waves,&! type of fast waves solver
       lcond,        & ! forecast with condensation/evaporation
       lseaice,      & ! forecast with sea ice model
       llake,        & ! forecast with lake model
       lsso,         & ! forecast with sub-grid scale orography scheme
       lforest,      & ! if .true., run with forest (evergreen and deciduous)
       lrad,         & ! if .true., run with radiation scheme
       lsoil,        & ! if .true., run with soil model
       lw_freeslip,  & ! if .TRUE.: with free slip lateral boundary condition and
                                ! if .FALSE. specified lateral boundary values for w
       llm,          & ! if .TRUE., running with lowered upper boundary
       lmulti_layer, & ! run multi-layer soil model
       lmulti_snow,  & ! run multi-layer snow model
       rdheight,     & ! bottom height of Rayleigh damping layer
       lspubc          ! with Rayleigh damping in the upper levels

! end of data_runcontrol
!------------------------------------------------------------------------------

  USE data_soil       , ONLY :   &
       cporv,        & !  pore volume (fraction of volume)
       cf_snow,      & !  parameter for the calculation of the 
       cdzw12,       & !  thickness of upper soil water layer in 
                       !  two-layer model         
       cdzw22,       & !  thickness of lower soil water layer in 
                       !  two-layer model      
       cdzw13,       & !  thickness of upper soil water layer in 
                       !  three-layer model
       cdzw23,       & !  thickness of middle soil water layer in 
                       !  three-layer model 
       cdzw33          !  thickness of lower soil water layer in 
                       !  three-layer model

!------------------------------------------------------------------------------

  USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------
       rho0       ,    & ! reference density at the full model levels    (kg/m3)
       dp0        ,    & ! pressure thickness of model layers            ( Pa  )
       p0         ,    & ! reference pressure at full levels             ( Pa  )
       p0hl       ,    & ! reference pressure at half levels             ( Pa  )
       t0hl       ,    & ! reference temperature at half levels          ( K   )
       t0         ,    & ! reference temperature at full levels          ( K   )
       hhl        ,    & ! geometical height of half levels              ( m   )

! 2. external parameter fields                                        (unit)
! ----------------------------
       hsurf      ,    & ! geometical heigt of surface topography        (m)
       gz0        ,    & !           * g                                 (m2/s2)
       fr_land    ,    & ! fraction of land in a grid element              --
       soiltyp    ,    & ! type of the soil (keys 1-10)                     --
       vio3       ,    & ! vertical integrated ozone contents            (pa O3)
       hmo3       ,    & ! ozone maximum                                 ( pa  )
       plcov      ,    & ! fraction of plant cover                         --
       lai        ,    & ! leaf area index of plants                       --
       rootdp     ,    & ! depth of the roots                            ( m  )
       for_e      ,    & ! ground fraction covered by evergreen forest     --
       for_d      ,    & ! ground fraction covered by deciduous forest     --
       llandmask  ,    & ! landpoint mask
       sso_stdh   ,    & ! standard deviation of sub-grid scale orography ( m   )
       sso_gamma  ,    & ! anisotropy of sub-grid scale orography          --
       sso_theta  ,    & ! angle betw. principal axis of orography and E ( rad )
       sso_sigma  ,    & ! mean slope of sub-grid scale orography          --
       fr_lake    ,    & ! lake fraction in a grid element [0,1]         (  -- )
       depth_lk   ,    & ! lake depth (SET TO -1.0 AT NON-LAKE POINTS!)  (  m  )
       fetch_lk   ,    & ! wind fetch over lake                          (  m  )
       dp_bs_lk   ,    & ! depth of the thermally active layer
                         ! of bottom sediments                           (  m  )
       t_bs_lk    ,    & ! climatological temperature at the bottom of
                         ! the thermally active layer of sediments       (  K  )
       gamso_lk   ,    & ! attenuation coefficient for
                         ! solar radiation in lake water                 ( 1/m )
       rlat       ,    & ! geographical latitude                         ( rad )
       rlon       ,    & ! geographical longitude                        ( rad )
       crlat      ,    & ! cosine of transformed latitude
       acrlat     ,    & ! 1 / ( crlat * radius of the earth )           ( 1/m )
       tgrlat     ,    & ! tangens of transformed latitude                 --
       tch               ! transfer coefficient for heat and moisture      ( -- )

  USE data_fields     , ONLY :   &

! 3. prognostic variables                                             (unit)
! -----------------------
       u          ,    & ! zonal wind speed                              ( m/s )
       v          ,    & ! meridional wind speed                         ( m/s )
       w          ,    & ! vertical wind speed (defined on half levels)  ( m/s )
       t          ,    & ! temperature                                   (  k  )
       pp         ,    & ! deviation from the reference pressure         ( pa  )
       rho        ,    & ! total density of moist air                    (kg/m3)

! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    time derivative by diabatic and adiabatic processes 
!    without sound-wave terms
       ttens      ,    & ! t-tendency without sound-wave terms            ( K/s)

! 5. fields for surface values and soil model variables               (unit )
! -----------------------------------------------------
       ps        ,     & ! surface pressure                              ( pa  )
       t_snow    ,     & ! temperature of the snow-surface               (  k  )
       t_snow_mult,    & ! temperature of the snow-surface               (  k  )
       t_s       ,     & ! temperature of the ground surface             (  k  )
       t_g       ,     & ! weighted surface temperature                  (  k  )
       qv_s      ,     & ! specific water vapor content on the surface   (kg/kg)
       t_m       ,     & ! temperature between upper and medium 
                         ! soil layer                                    (  k  )
       t_cl      ,     & ! temperature between medium and lower
                         ! soil layer                                    (  k  )
       t_so      ,     & ! multi-layer soil temperature                  (  k  )
       w_snow    ,     & ! snow water equivalent                         (m H2O)
       w_i       ,     & ! water content of interception water           (m H2O)
       w_g1      ,     & ! water content of the upper soil layer         (m H2O)
       w_g2      ,     & ! water content of the medium soil layer        (m H2O)
       w_g3      ,     & ! water content of the lower soil layer         (m H2O)
                         ! (if NLWB=3, unused otherwise)
       w_so      ,     & ! multi-layer soil moisture                     (m H2O)
       w_so_ice  ,     & ! multi-layer soil ice                          (m H2O)
       w_cl      ,     & ! climatological water content                  (m H2O) 
       freshsnow ,     & ! weighting function indicating 'freshness' of snow
       t_ice     ,     & ! temperature at the snow-ice or air-ice interface (  K  )
       t_mnw_lk  ,     & ! mean temperature of the water column          (  K  )
       t_wml_lk  ,     & ! mixed-layer temperature                       (  K  )
       t_bot_lk  ,     & ! temperature at the water-bottom sediment
                         ! interface                                     (  K  )
       t_b1_lk   ,     & ! temperature at the bottom of the upper layer
                         ! of the sediments                              (  K  )
       c_t_lk    ,     & ! shape factor with respect to the
                         ! temperature profile in lake thermocline       (  -  )
       h_ice     ,     & ! ice thickness                                 (  m  )
       h_ml_lk   ,     & ! thickness of the mixed-layer                  (  m  )
       h_b1_lk   ,     & ! thickness of the upper layer
                         ! of bottom sediments                           (  m  )

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
       qrs               ! precipitation water (water loading)           (kg/kg)
       
  USE data_fields     , ONLY :   &

! 8. fields for the boundary values                                   (unit )
! ---------------------------------
       u_bd          , & ! boundary field for u                          ( m/s )
       v_bd          , & ! boundary field for v                          ( m/s )
       w_bd          , & ! boundary field for w                          ( m/s )
       t_bd          , & ! boundary field for t                          (  k  )
       pp_bd         , & ! boundary field for pp                         (  pa )
       qv_s_bd       , & ! boundary field for qv_s                       (kg/kg)
       t_snow_bd     , & ! boundary field for t_snow                     (  k  )
       t_s_bd        , & ! boundary field for t_s                        (  k  )
       t_m_bd        , & ! boundary field for t_m                        (  k  )
       t_cl_bd       , & ! boundary field for t_cl                       (  k  )
       w_snow_bd     , & ! boundary field for w_snow                     (m H2O)
       w_g1_bd       , & ! boundary field for w_g1                       (m H2O)
       w_g2_bd       , & ! boundary field for w_g2                       (m H2O)
       w_g3_bd       , & ! boundary field for w_g3                       (m H2O)
       hmo3_bd       , & ! boundary field for hmo3                       (m    )
       vio3_bd       , & ! boundary field for vio3                       (pa O3)
       w_cl_bd       , & ! boundary field for w_cl                       (m H2O)
       lai_bd        , & ! boundary field for lai                        ( --  )
       rootdp_bd     , & ! boundary field for rootdp                     (m    )
       plcov_bd          ! boundary field for plcov                      ( --  )

!------------------------------------------------------------------------------

  USE data_parallel,      ONLY :  &
       num_compute,     & ! number of compute PEs
       nproc,           & ! total number of processors: nprocx * nprocy
       nboundlines,     & ! number of boundary lines of the domain for which
                                ! no forecast is computed = overlapping boundary
                                ! lines of the subdomains
       ldatatypes,      & ! if .TRUE.: use MPI-Datatypes for some communications
       ltime_barrier,   & ! if .TRUE.: use additional barriers for determining the
                                ! load-imbalance
       ncomm_type,      & ! type of communication
       my_world_id,     & ! rank of this subdomain in the global communicator
       imp_integers,    & ! determines the correct INTEGER type used in the
                                ! model for MPI
       imp_reals,       & ! determines the correct REAL type used in the model
                                ! for MPI
       imp_character,   & ! determines the correct CHARACTER type used in the model
                                ! for MPI
       imp_logical,     & ! determines the correct LOGICAL type used in the model
                                ! for MPI
       icomm_world,     & ! communicator that belongs to igroup_world, i.e.
                                ! = MPI_COMM_WORLD
       my_cart_id,      & ! rank of this subdomain in the cartesian communicator
       my_cart_neigh,   & ! neighbors of this subdomain in the cartesian grid
       isubpos,         & ! positions of the subdomains in the total domain. Given
                                ! are the i- and the j-indices of the lower left and the
                                ! upper right grid point in the order
                                !                  i_ll, j_ll, i_ur, j_ur.
                                ! Only the interior of the domains are considered, not
                                ! the boundary lines.
       icomm_cart,      & ! communicator for the virtual cartesian topology
       iexch_req,       & ! stores the sends requests for the neighbor-exchange
                                ! that can be used by MPI_WAIT to identify the send
       intbuf, realbuf, charbuf, logbuf, & ! Buffers for distributing the namelists
       sendbuf,         & ! sending buffer for boundary exchange:
                                ! 1-4 are used for sending, 5-8 are used for receiving
       isendbuflen        ! length of one column of sendbuf

!------------------------------------------------------------------------------

  USE utilities,                ONLY: sleve_split_oro, uv2uvrot_vec, &
                                      linear_interpol, linear_interpol_vec, mittel_integral_vec

  USE parallel_utilities,       ONLY: gather_field, distribute_field, distribute_values, &
                                      global_values

  USE environment,              ONLY: exchg_boundaries, model_abort, &
                                      get_free_unit, release_unit

  USE meteo_utilities,          ONLY: tgcom, calps

  USE grid_metrics_utilities,   ONLY: init_grid_metrics, sqrtg_r_w, wgtfac

  USE vgrid_refatm_utils,       ONLY: reference_atmosphere, reference_atmosphere_2,      &
                                      reference_atmosphere_BVconst, lanalyt_calc_t0p0,   &
                                      vcoord, refatm, nfltvc, svc1, svc2,                &
                                      k_index_of_pressure_levels, khmax

  USE data_io, ONLY: root, ydate_ini, &
       lbdclim            ! boundary data in climate model     ! PIK  (D.Hauffe)
                          ! (in climate mode also some external parameters have
                          !  to be updated, which are held constant in forecast
                          !  mode; e.g. plant cover, root depth)

  USE data_turbulence, ONLY: &
       vel_min            ! minimal velocity scale [m/s]

!------------------------------------------------------------------------------

  USE src_tracer,              ONLY: trcr_get, trcr_new, trcr_errorstr, &
                                     trcr_get_ntrcr

  USE data_tracer,      ONLY:  T_ADV_OFF , T_ADV_ON  , T_DIFF_OFF, T_DIFF_ON, &
                               T_TURB_OFF, T_TURB_1D , T_CONV_OFF, T_CONV_ON, &
                               T_INI_ZERO, T_INI_FILE, T_LBC_ZERO, T_LBC_FILE,&
                               T_INI_USER, T_LBC_USER,                        &
                               T_LBC_ZEROGRAD        , T_LBC_CST             ,&
                               T_BBC_ZEROFLUX        , T_BBC_ZEROVAL         ,&
                               T_BBC_SURF_VAL        , T_RELAX_OFF           ,&
                               T_RELAX_FULL          , T_RELAX_INFLOW        ,&
                               T_DAMP_OFF            , T_DAMP_ON             ,&
                               T_CLP_OFF             , T_CLP_ON              ,&
                               T_ERR_NOTFOUND        , T_TURB_3D

!------------------------------------------------------------------------------

  IMPLICIT NONE

  !=================================================================================
  ! Global constants:
  !=================================================================================

  ! Base pressure for computing potential temperature [Pa]
  REAL (KIND=wp),     PARAMETER         ::   pt00 = 1.0E5_wp

  !=================================================================================
  ! "Registry" for artificial disturbances (e.g., warm bubbles, heating rate bubbles etc.)
  ! NOTE: only the disturbances, which are registered here, can be used in the
  ! respective subroutines in the "select case"-clauses !
  ! You may register your own "bubble(s)" here by adding the respective
  ! registry line(s)!
  !=================================================================================
  !
  ! 1) For subroutine set_tempdist()  
  !        (disturbances in the initial condidition of atmospheric fields)
  CHARACTER(len=20), PARAMETER  :: &
       reg_tempdist(5)        = (/ &
                                  'cos                 ', &
                                  'mcnider             ', &
                                  'hotspot             ', &
                                  'squall3D            ', &
                                  'SK94                '  &
                                /)
  !
  ! 2) For subroutine set_tempdist_tso()  
  !        (disturbances in the initial condidition of the soil for lsoil=.TRUE.)
  CHARACTER(len=20), PARAMETER  :: &
       reg_tempdist_tso(2)    = (/ &
                                  'cos-soil            ', &
                                  'hotspot-soil        '  &
                                /)
  !
  ! 3) For subroutine set_tempdist_bbc_ts()  
  !        (disturbances in the surface temperature T_S for lsoil=.FALSE.)
  CHARACTER(len=20), PARAMETER  :: &
       reg_tempdist_bbc_ts(1) = (/ &
                                  'hotspot-sfc         '  &
                                /)
  !
  ! 4) For subroutine artif_heatrate_dist()  
  !        (artif. heating rates in the atmosphere)
  CHARACTER(len=20), PARAMETER  :: &
       reg_heatrate_dist(4)   = (/ &
                                  'cos-hrd             ',  &
                                  'mcnider-hrd         ', &
                                  'hotspot-hrd         ', &
                                  'AS2005_hucmtexas-hrd'  &
                                /)
  !
  ! 5) For subroutine artif_heatrate_dist_tso()  
  !        (artif. heating rates in the soil for lsoil=.TRUE.)
  CHARACTER(len=20), PARAMETER  :: &
       reg_heatrate_dist_tso(2) = (/ &
                                  'cos-soil-hrd        ',  &
                                  'hotspot-soil-hrd    '  &
                                /)


  !=================================================================================
  ! Variables for IDEAL-namelist
  !=================================================================================

  ! Set debug flag and debug-level:
  LOGICAL                 :: ldebug_artif      ! Switch to turn on debug mode (providing more output)
  INTEGER(KIND=iintegers) :: idbg_artif_level  ! Debug level if ldebug_artif=.true.

  !.. Constant parameters:
  INTEGER(KIND=iintegers), PARAMETER :: lcbuf = 100   ! Has to be 100, do not change!

  ! Maximum number of mountains which can be defined in the namelist
  INTEGER(KIND=iintegers), PARAMETER :: nhill_max=50
  ! Maximum number of height levels allowed in namelist (NOTE: for grib output only <= 200 possible)
  INTEGER(KIND=iintegers), PARAMETER :: nvcoordvec_max = khmax
  ! Maximum number of bubbles which can be defined in the namelist
  INTEGER(KIND=iintegers), PARAMETER :: ntempdist_max=50   
  ! Maximum number of polytrope layers, which can be defined in the namelist
  INTEGER(KIND=iintegers), PARAMETER :: nlayers_poly_max=10   
  ! Maximum number of N-const. layers, which can be defined in the namelist
  INTEGER(KIND=iintegers), PARAMETER :: nlayers_nconst_max=10 
  ! Maximum number of linear wind profile layers, which can be defined in the namelist
  INTEGER(KIND=iintegers), PARAMETER :: nlayers_linwind_max=10  

  ! Global error variables
  INTEGER   (KIND=iintegers) :: ierror

  ! Generic tracer pointers:
  REAL(KIND=wp),     POINTER :: &
    ztrcr_now  (:,:,:)   => NULL(),       & ! tracer variable at nnow
    ztrcr_new  (:,:,:)   => NULL(),       & ! tracer variable at nnew
    ztrcr_bd   (:,:,:,:) => NULL()          ! tracer boundary field

  ! Microphysics tracers
  REAL(KIND=wp),     POINTER :: &
    qv       (:,:,:)   => NULL(),   & ! QV at nnew
    qc       (:,:,:)   => NULL(),   & ! QC at nnew
    qi       (:,:,:)   => NULL(),   & ! QI at nnew
    qr       (:,:,:)   => NULL(),   & ! QR at nnew
    qs       (:,:,:)   => NULL(),   & ! QS at nnew
    qg       (:,:,:)   => NULL(),   & ! QG at nnew
    qv_bd    (:,:,:,:) => NULL(),   & ! BC for QV
    qc_bd    (:,:,:,:) => NULL(),   & ! BC for QC
    qi_bd    (:,:,:,:) => NULL(),   & ! BC for QI
    qr_bd    (:,:,:,:) => NULL(),   & ! BC for QR
    qs_bd    (:,:,:,:) => NULL(),   & ! BC for QS
    qg_bd    (:,:,:,:) => NULL()      ! BC for QG

#ifdef TWOMOM_SB
  REAL(KIND=wp),     POINTER :: &
    qh       (:,:,:)   => NULL(),   & ! QH at nnew
    qnc      (:,:,:)   => NULL(),   & ! NCC at nnew
    qni      (:,:,:)   => NULL(),   & ! NCI at nnew
    qnr      (:,:,:)   => NULL(),   & ! NCR at nnew
    qns      (:,:,:)   => NULL(),   & ! NCS at nnew
    qng      (:,:,:)   => NULL(),   & ! NCG at nnew
    qnh      (:,:,:)   => NULL(),   & ! NCH at nnew
    qnc_bd   (:,:,:,:) => NULL(),   & ! BC for NCC
    qni_bd   (:,:,:,:) => NULL(),   & ! BC for NCi
    qnr_bd   (:,:,:,:) => NULL(),   & ! BC for NCR
    qns_bd   (:,:,:,:) => NULL(),   & ! BC for NCS
    qng_bd   (:,:,:,:) => NULL()      ! BC for NCG
#endif

  CHARACTER(len=lcbuf) :: zspacing_type   ! subtype of vertical coordinate system: 
  ! 'predefined', 'galchen', 'linear', 'vcoordvec'

  REAL(KIND=wp)     :: &
       exp_galchen,  & ! if ivctype = 2/3/4 and zspacing_type = 'galchen' this defines the
                       !   exponent in the Gal-Chen formula
       zz_top,       & ! Namelist-Parameter for input of model top (height of level 1)
       h_top,        & ! Actual model top (height of level 1)
       z0_c,         & ! roughness length for idealized simulation in m
       hcond_on,     & ! time to switch on Condensation/Evaporation and cloud microphysics
       vcflat,       & ! vertical coordinate parameter
       p0sl,         & ! reference pressure at sea level
       t0sl,         & ! reference temperature at sea level
       dt0lp,        & ! d (t0) / d (ln p0)
       delta_t,      & ! temp. difference between sea level and stratosphere (for irefatm = 2)
       h_scal,       & ! scale height for irefatm = 2
       bvref           ! Constant Brunt-Vaisala-Frequency for irefatm=3 [1/s]

  INTEGER(KIND=iintegers) ::  &
       ivctype,      & ! type of vertical coordinate
       irefatm         ! 1: old reference atm. based on dT/dlnp = const
                       ! 2: new reference atm. with asymptotically isothermal stratosphere


  REAL (KIND=wp)             ::        &
       vcoordvec(nvcoordvec_max) ! vector to specify vcoord explicitly in namelist

  INTEGER(KIND=iintegers) ::  &
       hill_combineaction(nhill_max)

  LOGICAL                    ::       &
       linitw_followeta, &
       lhill(nhill_max), lhill_2d(nhill_max), &
       linit_realoro

  ! i-Shift of the model domain center w.r.t. the center of the input orography data set from file
  INTEGER(KIND=iintegers) :: i_shift_realoro   
  ! j-Shift of the model domain center w.r.t. the center of the input orography data set from file
  INTEGER(KIND=iintegers) :: j_shift_realoro   

  CHARACTER(len=lcbuf) :: hill_type(nhill_max)

  CHARACTER(len=250) :: rasofile, orofile, z0file, frlandfile, &
       soiltypefile, plcovfile, laifile, rootdpfile, &
       forefile, fordfile, tsurffile

  REAL(KIND=wp)     :: hill_width_x(nhill_max), hill_width_y(nhill_max), &
       hillsideradius_y(nhill_max), hillheight(nhill_max), &
       hill_i(nhill_max), hill_j(nhill_max), &
       hill_rotangle(nhill_max), href_oro,  &
       zhillcutfact(nhill_max), &
       hillasym_x(nhill_max), hillasym_y(nhill_max)

  REAL(KIND=wp)     :: zo_boundary, exponent_windprof_boundary


  ! Defining the warm air bubble:

  INTEGER(KIND=iintegers) :: &
       ntstep_bubble(ntempdist_max) ,&     ! Time step of first occurence of the bubble
       ntstep_noise                 ,&     ! Time step of adding noise to the temperature field
       numbubs_tot                  ,&     ! Total number of temperature- and heatingrate disturbances
       numbubs_noise                ,&     ! Number of temperature- and heatingrate disturbances with noise
       numbubs_tbuf                 ,&     ! Number of temperature- and heatingrate disturbances needing a T-buffer
       bub_zeitzaehler(ntempdist_max) = 0  ! Timestep counter for every heating rate disturbance to determine 
                                           !    the time period for artificial heating

  REAL (KIND=wp),     ALLOCATABLE   ::       &
       bub_tnoisebuffer(:,:,:,:)   ! Field for artificial noise for warm bubbles.

  INTEGER (KIND=iintegers) :: &
       bub_index_noise(ntempdist_max)     ! Index variable to set up reduced bub_tnoisebuffer
                                          ! which only contains as much space as needed for the actually 
                                          !    defined bubbles with noise.

  REAL (KIND=wp),     ALLOCATABLE   ::       &
         theta_ini(:,:,:)            ! initital theta-profile

  LOGICAL                    ::       &
       lgsp_buffer, &
       lcond_buffer

  ! Defining the initial vertical profiles:
  INTEGER(KIND=iintegers) ::  &
       itype_artifprofiles       , &   ! switch to choose between file input or analytic profile or a combination of both
       itype_anaprof_tqv         , &   ! if analytic T/qv- profiles, the type of profile
       itype_anaprof_uv                ! if analytic u/v - profiles, the type of profile

  REAL(KIND=wp)     :: tkvhfix=0.0_wp   ! Constant vert. diffusion coeff. for heat [m**2/s]
  REAL(KIND=wp)     :: tkhhfix=0.0_wp   ! Constant horiz. diffusion coeff. for heat [m**2/s]
  REAL(KIND=wp)     :: tkvmfix=0.0_wp   ! Constant vert. diffusion coeff. for momentum [m**2/s]
  REAL(KIND=wp)     :: tkhmfix=0.0_wp   ! Constant horiz. diffusion coeff. for momentum [m**2/s]
  LOGICAL :: lnosurffluxes_m=.FALSE.     ! Switch to turn on free-slip lower BC by setting tcm to 0.0
  LOGICAL :: lnosurffluxes_h=.FALSE.     ! Switch to turn on no-surface-heat/moisture-flux lower BC by setting tch to 0.0
  
  
  LOGICAL :: ltempdist(ntempdist_max)   ! Switch(es) (up to 50) to release temperature disturbances
  CHARACTER(len=lcbuf) :: ctype_tempdist(ntempdist_max)   ! Type of temperature disturbance(s)
  REAL(KIND=wp)     :: htempdist(ntempdist_max)   ! Time for release of temperature disturbances [h]
  REAL(KIND=wp)     :: bub_centi(ntempdist_max)   ! Center (i) of temperature disturbances [grid points, real]
  REAL(KIND=wp)     :: bub_centj(ntempdist_max)   ! Center (j) of temperature disturbances [grid points, real]
  REAL(KIND=wp)     :: bub_centz(ntempdist_max)   ! Center (Z) of temperature disturbances [m]
  INTEGER(KIND=iintegers) :: bub_timespan(ntempdist_max)   ! Duration of release of temperature disturbance(s) [# time steps]
  REAL(KIND=wp)     :: bub_radx(ntempdist_max)   ! Length scale (X) of temperature disturbances [m]
  REAL(KIND=wp)     :: bub_rady(ntempdist_max)   ! Length scale (Y) of temperature disturbances [m]
  REAL(KIND=wp)     :: bub_radz(ntempdist_max)   ! Length scale (Z) of temperature disturbances [m]
  REAL(KIND=wp)     :: bub_rotangle(ntempdist_max)   ! Rotation angle of main axes of temperature disturbances [degrees]
  REAL(KIND=wp)     :: bub_dT(ntempdist_max)   ! Amplitude of the temperature disturbance [K]
  LOGICAL :: ladd_bubblenoise_t(ntempdist_max)   ! Switch(es) (up to 50) to overlay temperature disturbances with some noise
  REAL(KIND=wp)     :: bub_dT_bubblenoise(ntempdist_max) ! Rel. amplitude of the noise on the temperature disturbance(s) [K]
  REAL(KIND=wp)     :: bub_zi_mcnider(ntempdist_max)   ! Height parameter of the McNider-bubble [m]
  REAL(KIND=wp)     :: bub_zmax_mcnider(ntempdist_max)   ! Max. height of the McNider-bubble [m]
  REAL(KIND=wp)     :: bub_h0_mcnider(ntempdist_max)   ! Sens. heat flux for the McNider-bubble [W/m**2]
  REAL(KIND=wp)     :: bub_heatingrate(ntempdist_max)   ! Constant heating rate [K*s**-1]
  INTEGER(KIND=iintegers) :: bub_centk(ntempdist_max)   ! Center (k) of temperature disturbance(s= [grid point, integer]
  INTEGER(KIND=iintegers) :: nlayers_poly   ! Number of actually used polytrope layers
  REAL(KIND=wp)     :: h_poly(nlayers_poly_max)   ! Height boundaries of polytrope layers [m]
  REAL(KIND=wp)     :: t_poly(nlayers_poly_max)   ! temperatures at the lower boundaries of polytrope layers [K]
  !***************************************************************************************
  ! NOTE: THE GRADIENTS TGR_POLY RESP. RHGR_POLY ARE POSITIVE FOR *DECREASING* TEMPERATURE
  !       RESP. RELHUM WITH HEIGHT!
  !***************************************************************************************
  REAL(KIND=wp)     :: tgr_poly(nlayers_poly_max)   ! temperature gradients within the polytrope layers [K/m]
  REAL(KIND=wp)     :: rh_poly(nlayers_poly_max)   ! relhums at the lower boundaries of polytrope layers [-]
  REAL(KIND=wp)     :: rhgr_poly(nlayers_poly_max)   ! relhum gradients within the polytrope layers [1/m]
  REAL(KIND=wp)     :: u_infty   ! Scaling wind velocity for the artificial analytical wind profiles [m/s]
  REAL(KIND=wp)     :: href_wk   ! Scaling height (height of 70 % windspeed) for the Weisman-Klemp wind profile [m]
  REAL(KIND=wp)     :: hmin_wk   ! Base height for the  Weisman-Klemp wind profile [m]
  INTEGER(KIND=iintegers) :: nlayers_linwind   ! Number of actually used layers for the piecewise linear wind profile
  REAL(KIND=wp)     :: h_linwind(nlayers_linwind_max) ! Height boundaries of layers for the piecewise linear wind profile [m]
  REAL(KIND=wp)     :: u_linwind(nlayers_linwind_max) ! Windspeed at lower boundaries of piecewise linear wind layers [m/s]
  !***********************************************************************************
  ! NOTE: THE GRADIENT UGR_LINWIND IS POSITIVE FOR *INCREASING* WINDSPEED WITH HEIGHT!
  !***********************************************************************************
  REAL(KIND=wp)     :: ugr_linwind(nlayers_linwind_max) ! Windspeed gradients within the piecewise linear wind layers [1/s]
  REAL(KIND=wp)     :: h_tropo_wk   ! Tropopause height for the Weisman-Klemp thermodynamic profile [m]
  REAL(KIND=wp)     :: theta_0_wk   ! Pot. Temp. at z=hmin_wk for the Weisman-Klemp thermodynamic profile [K]
  REAL(KIND=wp)     :: theta_tropo_wk   ! Pot. Temp. at z=h_tropo_wk for the Weisman-Klemp thermodynamic profile [K]
  REAL(KIND=wp)     :: expo_theta_wk   ! Exponent of the pot. temp. profile for the Weisman-Klemp thermodynamic profile
  REAL(KIND=wp)     :: expo_relhum_wk   ! Exponent of the rel. hum. profile for the Weisman-Klemp thermodynamic profile
  REAL(KIND=wp)     :: rh_min_wk   ! Rel. hum. at z > h_tropo_wk for the Weisman-Klemp thermodynamic profile
  REAL(KIND=wp)     :: rh_max_wk   ! Max. Rel. hum. for the Weisman-Klemp thermodynamic profile
  REAL(KIND=wp)     :: qv_max_wk   ! Max. spec. hum. for the Weisman-Klemp thermodynamic profile
  INTEGER(KIND=iintegers) :: itype_topo   ! Type of analytical topography added to href or to previously read topo from file
  INTEGER(KIND=iintegers) :: itype_soil_c ! Type of specifying "constant" soil parameters
  REAL(KIND=wp)     :: fr_land_c   ! Constant land fraction for idealized soil parameters
  REAL(KIND=wp)     :: soiltyp_c   ! Constant soil type for idealized soil parameters
  REAL(KIND=wp)     :: plcov_c   ! Constant plant cover for idealized soil parameters
  REAL(KIND=wp)     :: lai_c   ! Constant LAI for idealized soil parameters
  REAL(KIND=wp)     :: rootdp_c   ! Constant root depth for idealized soil parameters [m]
  REAL(KIND=wp)     :: for_e_c   ! Constant ground fraction covered by evergreen forests for idealized soil parameters
  REAL(KIND=wp)     :: for_d_c   ! Constant ground fraction covered by deciduous forests for idealized soil parameters
  LOGICAL :: ladd_noise_t   ! Switch to add some noise to the T and w fields within the lowest 100 hPa with a vertical 
                            !     one-half-wave-cos-profile
  REAL(KIND=wp)     :: hadd_noise   ! Time to add the noise, if ladd_noise_t = .true.
  REAL(KIND=wp)     :: dT_noise   ! Max. Amplitude of the temperature noise, if ladd_noise_t = .true.
  INTEGER(KIND=iintegers) :: iseed_noise_t   ! (Optional) Seed for the random number generator for T-noise,
                                             ! if ladd_noise_t = .true.
  REAL(KIND=wp)     :: dW_noise   ! Max. Amplitude of the vertical velocity noise, if ladd_noise_t = .true.
  INTEGER(KIND=iintegers) :: iseed_noise_w   ! (Optional) Seed for the random number generator for W-noise,
                                             ! if ladd_noise_t = .true.
  LOGICAL :: lps_from_file ! Switch to indicate that surface pressure can be interpolated from pressure in radiosonde file,
                           !     because this pressure is reliable.
  LOGICAL :: rasofile_t_is_theta   ! Switch to indicate that radiosonde file contains theta instead of t.
  REAL(KIND=wp)     :: p_base_wk   ! Pressure at height hmin_wk for analytic Weisman-Klemp case [Pa]
  REAL(KIND=wp)     :: p_base_poly ! Pressure at height h_poly(1) for analytic polytrope layers case [Pa]
  REAL(KIND=wp)     :: t_tropo_wk ! Parameter in Weisman-Klemp theta-profile, specifying an estimate 
                                  ! of tropopause temperature
  LOGICAL :: lbub_rhconst(ntempdist_max)   ! For each bubble, specify const. rel. hum. during temperature disturbance
                                           ! heating  (.true.) or no moisture modification (.false.)
  REAL(KIND=wp)     :: p_base_nconst   ! Pressure at height h_nconst(1) for analytic N-const. layers case [Pa]
  REAL(KIND=wp)     :: theta0_base_nconst   ! Pot. temp. at height h_nconst(1) for analytic N-const. layers case [Pa]
  INTEGER(KIND=iintegers) :: nlayers_nconst   ! Number of actually used N-const. layers
  REAL(KIND=wp)     :: h_nconst(nlayers_nconst_max)   ! Height boundaries of N-const. layers [m]
  REAL(KIND=wp)     :: N_nconst(nlayers_nconst_max)   ! N (Brunt-Vaisala-frequencies) of N-const. layers [K]
  REAL(KIND=wp)     :: rh_nconst(nlayers_nconst_max)   ! relhums at the lower boundaries of N-const. layers [-]
  !***************************************************************************************
  ! NOTE: THE GRADIENTS RHGR_POLY ARE POSITIVE FOR *DECREASING* RELHUM WITH HEIGHT!
  !***************************************************************************************
  REAL(KIND=wp)     :: rhgr_nconst(nlayers_nconst_max)   ! relhum gradients within the N-const. layers [1/m]
  REAL(KIND=wp)     :: t_soil_c   ! Constant temperature of soil [K] (t_s if <0)
  REAL(KIND=wp)     :: t_water_c  ! Constant temperature of water [K] (t_s if <0)
  REAL(KIND=wp)     :: t_snow_c   ! Constant temperature of the snow [K] (t_s if <0)
  REAL(KIND=wp)     :: wf_soil_c  ! Constant water saturation of soil as fraction of pore volume (0.0-1.0)
  REAL(KIND=wp)     :: w_i_c      ! Constant water content of interception water [m H2O]
  REAL(KIND=wp)     :: w_snow_c   ! Constant snow water equivalent [m H2O]
  INTEGER(KIND=iintegers) :: itype_soil_tw   ! Type of specifying soil/snow temperature and water content
  CHARACTER(len=250) :: tsoilfile   ! Name of 2D ASCII file for reading t_soil
  CHARACTER(len=250) :: tsnowfile   ! Name of 2D ASCII file for reading t_snow
  CHARACTER(len=250) :: wfsoilfile   ! Name of 2D ASCII file for reading wf_soil
  CHARACTER(len=250) :: wsnowfile   ! Name of 2D ASCII file for reading w_snow
  CHARACTER(len=250) :: wifile   ! Name of 2D ASCII file for reading w_i
  CHARACTER(len=250) :: hicefile   ! Name of 2D ASCII file for reading h_ice
  CHARACTER(len=250) :: ssostdhfile   ! Name of 2D ASCII file for reading sso_stdh
  CHARACTER(len=250) :: ssogammafile   ! Name of 2D ASCII file for reading sso_gamma
  CHARACTER(len=250) :: ssothetafile   ! Name of 2D ASCII file for reading sso_theta
  CHARACTER(len=250) :: ssosigmafile   ! Name of 2D ASCII file for reading sso_sigma
  REAL(KIND=wp)     :: h_ice_c   ! ice thickness [m]
  REAL(KIND=wp)     :: t_ice_c   ! temperature at the snow-ice or air-ice interface [K]
  REAL(KIND=wp)     :: t_surf_c  ! Constant temperature of earth surface [K], t_s (if < 0, set to atmosph. T at the surface)
  LOGICAL :: lsensiflux_fix   ! Switch to turn on fixed sensible and latent heat fluxes at the surface (for lsoil=.false.)
  LOGICAL :: llatentflux_fix   ! In case of lsensiflux_fix=.true., turn on also a fixed latent heat flux at the surface 
                               ! (for lsoil=.false.)
  REAL(KIND=wp)     :: sensiflux_c   ! Fixed surface sensible heat flux [W m**-2]
  REAL(KIND=wp)     :: latentflux_LzuS   ! Ratio of latent to sensible surface heat flux in case of llatentflux_fix=.true.
  REAL(KIND=wp)     :: H0_rel_noise   ! Relative level of white noise on surface fluxes for lsensiflux_fix=.true.
  INTEGER(KIND=iintegers) :: iseed_noise_H0   ! (Optional) Seed for the random number generator for H0-noise, 
                                              ! if H0_rel_noise >= 1.0E-5
!!!namelisttag   used by the ed-script ed_addnl.in which automatically adds nlparams


  !=================================================================================

  !------------------------------------------------------------------------------

CONTAINS

  !------------------------------------------------------------------------------

  !============================================================================
  ! 
  ! Subroutine for getting  the microphysics tracers.
  ! 
  !============================================================================
  
  SUBROUTINE get_tracers

    !--------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine retrieves the globally needed microphysics tracers
    !   out of the tracer structure
    !
    !--------------------------------------------------------------------------

    ! Parameters
    CHARACTER(LEN=25)   :: yzroutine = 'get_tracers'
    CHARACTER (LEN=255) :: yzerrmsg

    ierror  = 0_iintegers
    yzerrmsg(:) = ' '
  
    ! Retrieve the microphysics tracers
    CALL trcr_get(ierror, idt_qv, ptr_tlev=nnew, ptr=qv, ptr_bd=qv_bd)
    IF ( ierror /= 0_iintegers ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qc, ptr_tlev=nnew, ptr=qc, ptr_bd=qc_bd)
    IF ( ierror /= 0_iintegers ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qi, ptr_tlev=nnew, ptr=qi, ptr_bd=qi_bd)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qr, ptr_tlev=nnew, ptr=qr, ptr_bd=qr_bd)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qs, ptr_tlev=nnew, ptr=qs, ptr_bd=qs_bd)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qg, ptr_tlev=nnew, ptr=qg, ptr_bd=qg_bd)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
#ifdef TWOMOM_SB
    CALL trcr_get(ierror, idt_qh, ptr_tlev=nnew, ptr=qh)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qnc, ptr_tlev=nnew, ptr=qnc, ptr_bd=qnc_bd)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror) 
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qnr, ptr_tlev=nnew, ptr=qnr, ptr_bd=qnr_bd)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qni, ptr_tlev=nnew, ptr=qni, ptr_bd=qni_bd)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror) 
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qns, ptr_tlev=nnew, ptr=qns, ptr_bd=qns_bd)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qng, ptr_tlev=nnew, ptr=qng, ptr_bd=qng_bd)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
    CALL trcr_get(ierror, idt_qnh, ptr_tlev=nnew, ptr=qnh)
    IF ( ierror /= 0_iintegers .AND. ierror /= T_ERR_NOTFOUND ) THEN
      yzerrmsg = trcr_errorstr(ierror)
      CALL model_abort(my_cart_id, ierror, yzerrmsg, yzroutine)
    ENDIF
#endif

  END SUBROUTINE get_tracers

  !=================================================================================

  !=================================================================================
  !
  ! Subroutine for reading the namelist ARTIFCTL from file 'INPUT_IDEAL'.
  !
  ! This function is called in organize_data() after reading the gribout-namelist(s)
  ! in any case, not just for lartif_data = .true. This is intended to, e.g.,
  ! enable, in real case simulations, the specification of "warm bubbles", or
  ! the specification of constant turbulent diffusion coefficients,
  ! or the use of free-slip BCs. 
  !
  ! If file 'INPUT_IDEAL' does not exist, nothing is read and all namelist
  ! parameters take on their default values, which are neutral to
  ! operational runs.
  !
  ! HOWEVER, TO ACTIVATE SUCH IDEALIZED FEATURES IN REAL CASE SIMULATIONS,
  ! THE USER HAS TO MANUALLY EDIT LMORG.F90 AND/OR ORGANIZE_PHYSICS.F90.
  ! 
  ! In this case, not all the namelist parameters take effect, but only these:
  ! 
  ! - parameters for warm bubbles
  ! - parameters for itype_turb = 100: tkvhfix, tkhhfix, tkvmfix, tkhmfix
  !   (fixed diffusion coefficients)
  ! - switches "lnosurffluxes_m" and "lnosurffluxes_h" to enable runs without surface fluxes
  !   (free-slip-condition and/or no surface heat and moisture fluxes)
  !   THESE ARE USED IN SRC_OUTPUT.F90, ORGANIZE_PHYSICS.F90 AND FAST_WAVES_RK.F90
  !
  !=================================================================================

  SUBROUTINE input_artifctl (nuspecif, ierrstat)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine organizes the input of the NAMELIST-group artifctl.
    !   The group contains variables for idealized simulations; 
    !   used only if lartif_data = .true.
    !
    !------------------------------------------------------------------------------

    ! Parameters
    INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
         nuspecif        ! Unit number for protocolling the task

    INTEGER   (KIND=iintegers),   INTENT (out)      ::        &
         ierrstat        ! error status

    INTEGER (KIND=iintegers)   :: i, k, ierr, iz_err, nrbuf, nibuf, ncbuf, nlbuf, nuin

    CHARACTER(len=250) :: nlartiffile, io_errmsg
    LOGICAL            :: nlfile_exists, bub_reg_check

    ! Variables for default values of namelist parameters:

    INTEGER (KIND=iintegers)   ::       &
         ivctype_d,     & ! type of vertical coordinate
         irefatm_d        ! 1: old reference atm. based on dT/dlnp = const
                          ! 2: new reference atm. with asymptotically isothermal stratosphere

    REAL (KIND=wp)             ::       &
         p0sl_d,         & ! reference pressure at sea level
         t0sl_d,         & ! reference temperature at sea level
         dt0lp_d,       & ! d (t0) / d (ln p0)
         delta_t_d,   & ! temp. difference between sea level and stratosphere (for irefatm = 2)
         h_scal_d        ! scale height for irefatm = 2

    CHARACTER(len=lcbuf) :: zspacing_type_d ! subtype of vertical coordinate system: 
    ! 'predefined', 'galchen', 'linear', 'vcoordvec'
    REAL(KIND=wp)     :: exp_galchen_d, & ! if ivctype = 2/3/4 and zspacing_type = 'galchen' this defines the
                                ! exponent in the Gal-Chen formula
         vcflat_d, &      ! height were model levels become flat
         zz_top_d, &      ! Model Top (height of level 1)
         z0_c_d, &    ! roughness length for idealized simulation in m
         hcond_on_d, &    ! time to switch on Condensation/Evaporation and cloud microphysics
         svc1_d, svc2_d

    REAL (KIND=wp)             ::        &
         vcoordvec_d(nvcoordvec_max)

    INTEGER (KIND=iintegers)   ::        &
         nfltvc_d


    ! Defaults for namelist variables (artificial orography):
    INTEGER(KIND=iintegers) ::  &
         hill_combineaction_d(nhill_max)

    LOGICAL                    ::       &
         linitw_followeta_d, &
         lhill_d(nhill_max), lhill_2d_d(nhill_max), &
         linit_realoro_d
    ! i-Shift of the model domain center w.r.t. the center of the input orography data set from file
    INTEGER(KIND=iintegers) :: i_shift_realoro_d   
    ! j-Shift of the model domain center w.r.t. the center of the input orography data set from file
    INTEGER(KIND=iintegers) :: j_shift_realoro_d   

    CHARACTER(len=lcbuf) :: hill_type_d(nhill_max)

    CHARACTER(len=250) :: rasofile_d, orofile_d, z0file_d, frlandfile_d, &
         soiltypefile_d, plcovfile_d, laifile_d, rootdpfile_d, &
         forefile_d, fordfile_d, tsurffile_d

    REAL(KIND=wp)     :: hill_width_x_d(nhill_max), hill_width_y_d(nhill_max), &
         hillsideradius_y_d(nhill_max), hillheight_d(nhill_max), &
         hill_i_d(nhill_max), hill_j_d(nhill_max), &
         hill_rotangle_d(nhill_max), href_oro_d,  &
         zhillcutfact_d(nhill_max), &
         hillasym_x_d(nhill_max), hillasym_y_d(nhill_max)

    REAL(KIND=wp)     :: zo_boundary_d, exponent_windprof_boundary_d, &
         bvref_tconst, zmax_nconst, zT0_nconst, zloghuge

    ! Defining the initial vertical profiles:
    INTEGER(KIND=iintegers) ::  &
         itype_artifprofiles_d       , & ! switch to choose between file input or analytic profile or a combination of both
         itype_anaprof_tqv_d         , & ! if analytic T/qv- profiles, the type of profile
         itype_anaprof_uv_d              ! if analytic u/v - profiles, the type of profile

    REAL(KIND=wp)     :: tkvhfix_d   ! Constant vert. diffusion coeff. for heat [m**2/s]
    REAL(KIND=wp)     :: tkhhfix_d   ! Constant horiz. diffusion coeff. for heat [m**2/s]
    REAL(KIND=wp)     :: tkvmfix_d   ! Constant vert. diffusion coeff. for momentum [m**2/s]
    REAL(KIND=wp)     :: tkhmfix_d   ! Constant horiz. diffusion coeff. for momentum [m**2/s]
    LOGICAL :: lnosurffluxes_m_d   ! Switch to turn on free-slip lower BC by setting tcm to 0.0
    LOGICAL :: lnosurffluxes_h_d   ! Switch to turn on no-surface-heat/moisture-flux lower BC by setting tch to 0.0
    LOGICAL :: ltempdist_d(ntempdist_max)   ! Switch(es) (up to 50) to release temperature disturbances
    CHARACTER(len=lcbuf) :: ctype_tempdist_d(ntempdist_max)   ! Type of temperature disturbance(s)
    REAL(KIND=wp)     :: htempdist_d(ntempdist_max)   ! Time for release of temperature disturbances [h]
    REAL(KIND=wp)     :: bub_centi_d(ntempdist_max)   ! Center (i) of temperature disturbances [grid points, real]
    REAL(KIND=wp)     :: bub_centj_d(ntempdist_max)   ! Center (j) of temperature disturbances [grid points, real]
    REAL(KIND=wp)     :: bub_centz_d(ntempdist_max)   ! Center (Z) of temperature disturbances [m]
    INTEGER(KIND=iintegers) :: bub_timespan_d(ntempdist_max)   ! Duration of release of temperature disturbance(s)
                                                               ! [# time steps]
    REAL(KIND=wp)     :: bub_radx_d(ntempdist_max)   ! Length scale (X) of temperature disturbances [m]
    REAL(KIND=wp)     :: bub_rady_d(ntempdist_max)   ! Length scale (Y) of temperature disturbances [m]
    REAL(KIND=wp)     :: bub_radz_d(ntempdist_max)   ! Length scale (Z) of temperature disturbances [m]
    REAL(KIND=wp)     :: bub_rotangle_d(ntempdist_max)   ! Rotation angle of main axes of temperature disturbances [degrees]
    REAL(KIND=wp)     :: bub_dT_d(ntempdist_max)   ! Amplitude of the temperature disturbance [K]
    LOGICAL :: ladd_bubblenoise_t_d(ntempdist_max) ! Switch(es) (up to 50) to overlay temperature disturbances
                                                   ! with some noise
    REAL(KIND=wp)     :: bub_dT_bubblenoise_d(ntempdist_max)   ! Rel. amplitude of the noise on the temperature
                                                               ! disturbance(s)[K]
    REAL(KIND=wp)     :: bub_zi_mcnider_d(ntempdist_max)   ! Height parameter of the McNider-bubble [m]
    REAL(KIND=wp)     :: bub_zmax_mcnider_d(ntempdist_max)   ! Max. height of the McNider-bubble [m]
    REAL(KIND=wp)     :: bub_h0_mcnider_d(ntempdist_max)   ! Sens. heat flux for the McNider-bubble [W/m**2]
    REAL(KIND=wp)     :: bub_heatingrate_d(ntempdist_max)   ! Constant heating rate [K*s**-1]
    INTEGER(KIND=iintegers) :: bub_centk_d(ntempdist_max)   ! Center (k) of temperature disturbance(s= [grid point, integer]
    INTEGER(KIND=iintegers) :: nlayers_poly_d   ! Number of actually used polytrope layers
    REAL(KIND=wp)     :: h_poly_d(nlayers_poly_max)   ! Height boundaries of polytrope layers [m]
    REAL(KIND=wp)     :: t_poly_d(nlayers_poly_max)   ! temperatures at the lower boundaries of polytrope layers [K]
    REAL(KIND=wp)     :: tgr_poly_d(nlayers_poly_max)   ! temperature gradients within the polytrope layers [K/m]
    REAL(KIND=wp)     :: rh_poly_d(nlayers_poly_max)   ! relhums at the lower boundaries of polytrope layers [-]
    REAL(KIND=wp)     :: rhgr_poly_d(nlayers_poly_max)   ! relhum gradients within the polytrope layers [1/m]
    REAL(KIND=wp)     :: u_infty_d   ! Scaling wind velocity for the artificial analytical wind profiles [m/s]
    REAL(KIND=wp)     :: href_wk_d   ! Scaling height (height of 70 % windspeed) for the Weisman-Klemp wind profile [m]
    REAL(KIND=wp)     :: hmin_wk_d   ! Base height for the  Weisman-Klemp wind profile [m]
    INTEGER(KIND=iintegers) :: nlayers_linwind_d   ! Number of actually used layers for the piecewise linear wind profile
    REAL(KIND=wp)     :: h_linwind_d(nlayers_linwind_max) ! Height boundaries of layers for the piecewise
                                                          ! linear wind profile [m]
    REAL(KIND=wp)     :: u_linwind_d(nlayers_linwind_max) ! Windspeed at lower boundaries of piecewise linear
                                                          ! wind layers [m/s]
    REAL(KIND=wp)     :: ugr_linwind_d(nlayers_linwind_max)   ! Windspeed gradients within the piecewise linear
                                                              ! wind layers [1/s]
    REAL(KIND=wp)     :: h_tropo_wk_d   ! Tropopause height for the Weisman-Klemp thermodynamic profile [m]
    REAL(KIND=wp)     :: theta_0_wk_d   ! Pot. Temp. at z=hmin_wk for the Weisman-Klemp thermodynamic profile [K]
    REAL(KIND=wp)     :: theta_tropo_wk_d   ! Pot. Temp. at z=h_tropo_wk for the Weisman-Klemp thermodynamic profile [K]
    REAL(KIND=wp)     :: expo_theta_wk_d   ! Exponent of the pot. temp. profile for the Weisman-Klemp thermodynamic profile
    REAL(KIND=wp)     :: expo_relhum_wk_d   ! Exponent of the rel. hum. profile for the Weisman-Klemp thermodynamic profile
    REAL(KIND=wp)     :: rh_min_wk_d   ! Rel. hum. at z > h_tropo_wk for the Weisman-Klemp thermodynamic profile
    REAL(KIND=wp)     :: rh_max_wk_d   ! Max. Rel. hum. for the Weisman-Klemp thermodynamic profile
    REAL(KIND=wp)     :: qv_max_wk_d   ! Max. spec. hum. for the Weisman-Klemp thermodynamic profile
    INTEGER(KIND=iintegers) :: itype_topo_d ! Type of analytical topography added to href or to a previously
                                            ! read topo from file
    INTEGER(KIND=iintegers) :: itype_soil_c_d ! Type of specifying "constant" soil parameters
    REAL(KIND=wp)     :: fr_land_c_d   ! Constant land fraction for idealized soil parameters
    REAL(KIND=wp)     :: soiltyp_c_d   ! Constant soil type for idealized soil parameters
    REAL(KIND=wp)     :: plcov_c_d   ! Constant plant cover for idealized soil parameters
    REAL(KIND=wp)     :: lai_c_d   ! Constant LAI for idealized soil parameters
    REAL(KIND=wp)     :: rootdp_c_d   ! Constant root depth for idealized soil parameters [m]
    REAL(KIND=wp)     :: for_e_c_d   ! Constant ground fraction covered by evergreen forests for idealized soil parameters
    REAL(KIND=wp)     :: for_d_c_d   ! Constant ground fraction covered by deciduous forests for idealized soil parameters
    LOGICAL :: ladd_noise_t_d   ! Switch to add some noise to the T and w fields within the lowest 100 hPa with a vertical 
                                !     one-half-wave-cos-profile
    REAL(KIND=wp)     :: hadd_noise_d   ! Time to add the noise, if ladd_noise_t = .true.
    REAL(KIND=wp)     :: dT_noise_d   ! Max. Amplitude of the temperature noise, if ladd_noise_t = .true.
    INTEGER(KIND=iintegers) :: iseed_noise_t_d   ! (Optional) Seed for the random number generator for T-noise, 
                                                 ! if ladd_noise_t = .true.
    REAL(KIND=wp)     :: dW_noise_d   ! Max. Amplitude of the vertical velocity noise, if ladd_noise_t = .true.
    INTEGER(KIND=iintegers) :: iseed_noise_w_d   ! (Optional) Seed for the random number generator for W-noise,
                                                 ! if ladd_noise_t = .true.
    LOGICAL :: lps_from_file_d  ! Switch to indicate that surface pressure can be interpolated from pressure in radiosonde 
                                !   file, because this pressure is reliable.
    LOGICAL :: rasofile_t_is_theta_d   ! Switch to indicate that radiosonde file contains theta instead of t.
    REAL(KIND=wp)     :: p_base_wk_d   ! Pressure at height hmin_wk for analytic Weisman-Klemp case [Pa]
    REAL(KIND=wp)     :: p_base_poly_d   ! Pressure at height h_poly(1) for analytic polytrope layers case [Pa]
    REAL(KIND=wp)     :: t_tropo_wk_d   ! Parameter in Weisman-Klemp theta-profile, specifying an estimate 
                                        ! of the tropopause temperature
    LOGICAL :: lbub_rhconst_d(ntempdist_max)   ! For each bubble, specify const. rel. hum. during temperature 
                                               ! disturbance heating (.true.) 
    ! or no moisture modification (.false.)
    REAL(KIND=wp)     :: p_base_nconst_d   ! Pressure at height h_nconst(1) for analytic N-const. layers case [Pa]
    REAL(KIND=wp)     :: theta0_base_nconst_d   ! Pot. temp. at height h_nconst(1) for analytic N-const. layers case [Pa]
    INTEGER(KIND=iintegers) :: nlayers_nconst_d   ! Number of actually used N-const. layers
    REAL(KIND=wp)     :: h_nconst_d(nlayers_nconst_max)   ! Height boundaries of N-const. layers [m]
    REAL(KIND=wp)     :: N_nconst_d(nlayers_nconst_max)   ! N (Brunt-Vaisala-frequencies) of N-const. layers [K]
    REAL(KIND=wp)     :: rh_nconst_d(nlayers_nconst_max)   ! relhums at the lower boundaries of N-const. layers [-]
    REAL(KIND=wp)     :: rhgr_nconst_d(nlayers_nconst_max)   ! relhum gradients within the N-const. layers [1/m]
    REAL(KIND=wp)     :: bvref_d   ! Constant Brunt-Vaisala-Frequency for irefatm=3 [1/s]
    LOGICAL :: ldebug_artif_d   ! Switch to turn on debug mode (providing more output)
    INTEGER(KIND=iintegers) :: idbg_artif_level_d   ! Debug level if ldebug_artif=.true.
    REAL(KIND=wp)     :: t_soil_c_d   ! Constant temperature of soil [K] (t_s if <0)
    REAL(KIND=wp)     :: t_water_c_d  ! Constant temperature of water [K] (t_s if <0)
    REAL(KIND=wp)     :: t_snow_c_d   ! Constant temperature of the snow [K] (t_s if <0)
    REAL(KIND=wp)     :: wf_soil_c_d  ! Constant water saturation of soil as fraction of pore volume (0.0-1.0)
    REAL(KIND=wp)     :: w_i_c_d      ! Constant water content of interception water [m H2O]
    REAL(KIND=wp)     :: w_snow_c_d   ! Constant snow water equivalent [m H2O]
    INTEGER(KIND=iintegers) :: itype_soil_tw_d   ! Type of specifying soil/snow temperature and water content
    CHARACTER(len=250) :: tsoilfile_d   ! Name of 2D ASCII file for reading t_soil
    CHARACTER(len=250) :: tsnowfile_d   ! Name of 2D ASCII file for reading t_snow
    CHARACTER(len=250) :: wfsoilfile_d   ! Name of 2D ASCII file for reading wf_soil
    CHARACTER(len=250) :: wsnowfile_d   ! Name of 2D ASCII file for reading w_snow
    CHARACTER(len=250) :: wifile_d   ! Name of 2D ASCII file for reading w_i
    CHARACTER(len=250) :: hicefile_d   ! Name of 2D ASCII file for reading h_ice
    CHARACTER(len=250) :: ssostdhfile_d   ! Name of 2D ASCII file for reading sso_stdh
    CHARACTER(len=250) :: ssogammafile_d   ! Name of 2D ASCII file for reading sso_gamma
    CHARACTER(len=250) :: ssothetafile_d   ! Name of 2D ASCII file for reading sso_theta
    CHARACTER(len=250) :: ssosigmafile_d   ! Name of 2D ASCII file for reading sso_sigma
    REAL(KIND=wp)     :: h_ice_c_d   ! ice thickness [m]
    REAL(KIND=wp)     :: t_ice_c_d   ! temperature at the snow-ice or air-ice interface [K]
    REAL(KIND=wp)     :: t_surf_c_d   ! Constant temperature of earth surface [K], t_s (if < 0, set to atmosph. T at the surface)
    LOGICAL :: lsensiflux_fix_d   ! Switch to turn on fixed sensible and latent heat fluxes at the surface (for lsoil=.false.)
    LOGICAL :: llatentflux_fix_d   ! In case of lsensiflux_fix=.true., turn on also a fixed latent heat flux at the surface 
                                   ! (for lsoil=.false.)
    REAL(KIND=wp)     :: sensiflux_c_d   ! Fixed surface sensible heat flux [W m**-2]
    REAL(KIND=wp)     :: latentflux_LzuS_d   ! Ratio of latent to sensible surface heat flux in case of llatentflux_fix=.true.
    REAL(KIND=wp)     :: H0_rel_noise_d   ! Relative level of white noise on surface fluxes for lsensiflux_fix=.true.
    INTEGER(KIND=iintegers) :: iseed_noise_H0_d   ! (Optional) Seed for the random number generator for H0-noise, 
                                                  ! if H0_rel_noise >= 1.0E-5
!!!defaulttag   used by the ed-script ed_addnl.in which automatically adds nlparams


    NAMELIST /ARTIFCTL/ ivctype, irefatm, p0sl, t0sl, dt0lp, delta_t, h_scal, &
         zspacing_type, exp_galchen, vcflat, zz_top, z0_c, vcoordvec, &
         nfltvc, svc1, svc2, &
         linitw_followeta, zo_boundary, exponent_windprof_boundary, &
         linit_realoro, orofile, i_shift_realoro , j_shift_realoro , &
         z0file, frlandfile, soiltypefile, plcovfile, laifile, rootdpfile, tsurffile, &
         forefile, fordfile, tsoilfile, tsnowfile, wfsoilfile, wsnowfile, wifile, hicefile, &
         ssostdhfile, ssogammafile, ssothetafile, ssosigmafile, &
         lhill, lhill_2d, hill_type, hill_width_x, hill_width_y, hillsideradius_y, &
         hillheight, hill_i, hill_j, hill_rotangle, href_oro, zhillcutfact, &
         hillasym_x, hillasym_y, hill_combineaction, &
         rasofile, itype_artifprofiles, itype_anaprof_tqv, itype_anaprof_uv ,&
         tkvhfix, tkhhfix, tkvmfix , tkhmfix, lnosurffluxes_m , lnosurffluxes_h , &
         ltempdist , ctype_tempdist , htempdist , &
         bub_centi , bub_centj , bub_centz , bub_timespan , &
         bub_radx , bub_rady , bub_radz , bub_rotangle , bub_dT , &
         ladd_bubblenoise_t , bub_dT_bubblenoise , &
         bub_zi_mcnider , bub_zmax_mcnider , bub_h0_mcnider , &
         bub_heatingrate , bub_centk , &
         nlayers_poly , h_poly , t_poly , tgr_poly , rh_poly , rhgr_poly , &
         u_infty , href_wk , hmin_wk , &
         nlayers_linwind , h_linwind , u_linwind , ugr_linwind , &
         h_tropo_wk , theta_0_wk , theta_tropo_wk , expo_theta_wk , expo_relhum_wk , &
         rh_min_wk , rh_max_wk , qv_max_wk , &
         itype_topo , itype_soil_c , fr_land_c , soiltyp_c , plcov_c , &
         lai_c , rootdp_c , for_e_c , for_d_c , &
         ladd_noise_t , hadd_noise , dT_noise , dW_noise , iseed_noise_t , iseed_noise_w , &
         lps_from_file , rasofile_t_is_theta , p_base_wk , p_base_poly , &
         t_tropo_wk , lbub_rhconst , p_base_nconst , theta0_base_nconst , &
         nlayers_nconst , h_nconst , N_nconst , rh_nconst , &
         rhgr_nconst , bvref , ldebug_artif , idbg_artif_level , &
         t_soil_c , t_water_c , t_snow_c , wf_soil_c , w_i_c , w_snow_c , itype_soil_tw, h_ice_c, t_ice_c , t_surf_c , &
         lsensiflux_fix , llatentflux_fix , sensiflux_c , latentflux_LzuS , H0_rel_noise , iseed_noise_H0 , &
!!!addlisttag   used by the ed-script ed_addnl.in which automatically adds nlparams
         hcond_on

    !===============================================================================

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. input_artifctl() ...'
    END IF

    ierrstat = 0

    !===============================================================================


    !-------------------------------------------------------------------------------
    !- Section 0: Preparations, e.g, initializing chars with ' ':
    !-------------------------------------------------------------------------------

    DO i=1, LEN(hill_type(1))
      hill_type(:)(i:i) = ' '
    END DO
    DO i=1, LEN(ctype_tempdist(1))
      ctype_tempdist(:)(i:i) = ' '
    END DO
!!!cinitag   used by the ed-script ed_addnl.in which automatically adds nlparams


    !-------------------------------------------------------------------------------
    !- Section 1: Initialize the default variables
    !-------------------------------------------------------------------------------


    ivctype_d     = 2_iintegers
    irefatm_d     = 2_iintegers
    p0sl_d        = 1.0E5_wp
    t0sl_d        = t0_melt + 15.0_wp
    dt0lp_d       = 42.0_wp
    delta_t_d     = 75.0_wp
    h_scal_d      = 10000.0_wp
    bvref_d       = 0.01_wp

    ! Parameters for the height-Koordinate:
    zspacing_type_d(:) = ' ';   zspacing_type_d = 'galchen'   ! 'predefined', 'galchen', or 'linear'
    vcoordvec_d(:)  = -999.99_wp
    vcflat_d        =   11000.0_wp
    zz_top_d        =   22000.0_wp
    exp_galchen_d   =   2.6_wp
    nfltvc_d = 100_iintegers
    svc1_d   = 8000.0_wp
    svc2_d   = 5000.0_wp

    ! Parameters for an artificial mountain:
    hill_i_d       = (ie_tot-1) * 0.5_wp + 1.0_wp  ! center of model domain
    hill_j_d       = (je_tot-1) * 0.5_wp + 1.0_wp  ! center of model domain
    href_oro_d   =   0.0_wp
    linit_realoro_d = .FALSE.
    zo_boundary_d = 0.0_wp                   ! [m] height of an artificial wind speed boundary layer, e.g. 400.0 m
    exponent_windprof_boundary_d = 0.333_wp  ! Exponent for the power law wind speed profile within the boundary layer

    ! i-Shift of the model domain center w.r.t. the center of the input orography data set from file
    i_shift_realoro_d  = 0_iintegers 
    ! j-Shift of the model domain center w.r.t. the center of the input orography data set from file
    j_shift_realoro_d  = 0_iintegers 
    orofile_d(:) = ' '; orofile_d = 'noname.dat'
    z0file_d(:) = ' '; z0file_d = 'noname.dat'
    frlandfile_d(:) = ' '; frlandfile_d = 'noname.dat'
    soiltypefile_d(:) = ' '; soiltypefile_d = 'noname.dat'
    plcovfile_d(:) = ' '; plcovfile_d = 'noname.dat'
    laifile_d(:) = ' '; laifile_d = 'noname.dat'
    rootdpfile_d(:) = ' '; rootdpfile_d = 'noname.dat'
    forefile_d(:) = ' '; forefile_d = 'noname.dat'
    fordfile_d(:) = ' '; fordfile_d = 'noname.dat'
    tsoilfile_d(:) = ' ';  tsoilfile_d  = 'noname.dat'
    tsnowfile_d(:) = ' ';  tsnowfile_d  = 'noname.dat'
    wfsoilfile_d(:) = ' ';  wfsoilfile_d  = 'noname.dat'
    wsnowfile_d(:) = ' ';  wsnowfile_d  = 'noname.dat'
    wifile_d(:) = ' ';   wifile_d  = 'noname.dat'
    hill_combineaction_d = 1  ! 1 = Sum of all hills, 2 = MAX of the single hills
    hill_width_x_d =   5000.0_wp
    hill_width_y_d =   5000.0_wp
    hillsideradius_y_d =   hill_width_x_d
    hillheight_d    =   1000.0_wp
    hill_rotangle_d =  0.0_wp
    hillasym_x_d    =   1.0_wp ! Asymmetrieparameter: hill_width_x gilt fuer Westseite,
                                   !                      hill_width_x/hillasym_x auf Ostseite
    hillasym_y_d    =   1.0_wp ! Asymmetrieparameter: hill_width_y gilt fuer Suedseite,
                                   !                      hill_width_y/hillasym_y auf Nordseite
    lhill_d        =   .FALSE.
    hill_type_d   =   'gauss'    ! 'bellshaped' or 'gauss'
    lhill_2d_d     =   .TRUE.
    zhillcutfact_d = 0.0_wp

    tkvhfix_d  = 300.0_wp
    tkhhfix_d  = 300.0_wp
    tkvmfix_d  = 300.0_wp
    tkhmfix_d  = 300.0_wp
    lnosurffluxes_m_d  = .FALSE.
    lnosurffluxes_h_d  = .FALSE.
    ltempdist_d  = .FALSE.
    ctype_tempdist_d  = 'cos'
    htempdist_d  = 0.0_wp
    bub_centi_d  = (ie_tot-1) * 0.5_wp + 1.0_wp  ! center of model domain
    bub_centj_d  = (je_tot-1) * 0.5_wp + 1.0_wp  ! center of model domain
    bub_centz_d  = 1400.0_wp
    bub_centk_d  = 5_iintegers
    bub_timespan_d  = 100_iintegers
    bub_radx_d  = 10000.0_wp
    bub_rady_d  = 10000.0_wp
    bub_radz_d  = 1400.0_wp
    bub_rotangle_d  = 0.0_wp
    bub_dT_d  = 3.0_wp
    bub_heatingrate_d  = 0.005_wp
    ladd_bubblenoise_t_d  = .FALSE.
    bub_dT_bubblenoise_d  = 0.1_wp
    bub_zi_mcnider_d  = 4000.0_wp
    bub_zmax_mcnider_d  = 6000.0_wp
    bub_h0_mcnider_d  = 600.0_wp
    lbub_rhconst_d  = .FALSE.

    ! Defaults for polytrope atmosphere (itype_anaprof_tqv = 2):
    !    3 Layers (PBL, free troposphere, tropopause)
    !      remember: one can have a maximum of nlayers_poly_max layers!
    nlayers_poly_d  = 3_iintegers
    h_poly_d(:)  = HUGE(1.0_wp)
      h_poly_d(1) = 0.0_wp
      h_poly_d(2) = 2000.0_wp
      h_poly_d(3) = 11000.0_wp
    t_poly_d(:)   = 288.16_wp
      t_poly_d(1) = 288.16_wp
      t_poly_d(2) = 275.16_wp
      t_poly_d(3) = 216.66_wp
    tgr_poly_d(:)   = 0.0_wp    ! positive for decreasing temperature with height
      tgr_poly_d(1) = 0.0095_wp ! - (t_poly_d(2)-t_poly_d(1))/(h_poly_d(2)-h_poly_d(1))
      tgr_poly_d(2) = 0.0065_wp ! - (t_poly_d(3)-t_poly_d(2))/(h_poly_d(3)-h_poly_d(2))
      tgr_poly_d(3) = 0.0_wp
    rh_poly_d(:)   = 0.5_wp
      rh_poly_d(1) = 0.0_wp     
      rh_poly_d(2) = 0.0_wp   
      rh_poly_d(3) = 0.0_wp   
    rhgr_poly_d(:)   = 0.0_wp                  ! positive for decreasing humidity with height
      rhgr_poly_d(1) = -0.0_wp / 2000.0_wp ! - (rh_poly_d(2)-rh_poly_d(1))/(h_poly_d(2)-h_poly_d(1))
      rhgr_poly_d(2) = 0.0_wp / 9000.0_wp ! - (rh_poly_d(3)-rh_poly_d(2))/(h_poly_d(3)-h_poly_d(2))
      rhgr_poly_d(3) = 0.0_wp
    u_infty_d  = 20.0_wp
    href_wk_d  = 3000.0_wp
    hmin_wk_d  = 0.0_wp
    nlayers_linwind_d  = 2_iintegers
    h_linwind_d(:)  = HUGE(1.0_wp)
      h_linwind_d(1) = 0.0_wp
      h_linwind_d(2) = 2500.0_wp
    u_linwind_d(:)  = 0.0_wp
      u_linwind_d(1) = 0.0_wp
      u_linwind_d(2) = 20.0_wp
    ugr_linwind_d(:)  = 0.0_wp  ! positive for increasing windspeed with height
      ugr_linwind_d(1) = 20.0_wp / 2500.0_wp  ! (u_linwind_d(2)-u_linwind_d(1))/(h_linwind_d(2)-h_linwind_d(1))
      ugr_linwind_d(2) = 0.0_wp
    h_tropo_wk_d  = 12000.0_wp
    theta_0_wk_d  = 300.0_wp
    theta_tropo_wk_d  = 343.0_wp
    expo_theta_wk_d  = 1.25_wp
    expo_relhum_wk_d  = 1.25_wp
    rh_min_wk_d  = 0.25_wp
    rh_max_wk_d  = 1.0_wp
    qv_max_wk_d  = 0.012_wp
    itype_topo_d  = 1_iintegers
    itype_soil_c_d  = 1_iintegers
    z0_c_d      =   0.1_wp
    fr_land_c_d  = 1.0_wp
    soiltyp_c_d  = 3.0_wp
    plcov_c_d  = 0.6_wp
    lai_c_d  = 3.0_wp
    rootdp_c_d  = 0.7_wp
    for_e_c_d  = 0.2_wp
    for_d_c_d  = 0.2_wp
    ladd_noise_t_d  = .FALSE.
    hadd_noise_d  = 0.0_wp
    dT_noise_d  = 0.02_wp
    dW_noise_d  = 0.02_wp
    iseed_noise_t_d  = 606_iintegers
    iseed_noise_w_d  = 607_iintegers
    lps_from_file_d  = .FALSE.
    rasofile_t_is_theta_d  = .FALSE.
    p_base_wk_d  = 100000.0_wp
    p_base_poly_d  = 100000.0_wp
    t_tropo_wk_d  = -999.99_wp
    p_base_nconst_d  = 100000.0_wp
    theta0_base_nconst_d  = 300.0_wp
    nlayers_nconst_d  = 3_iintegers
    h_nconst_d(:)  = HUGE(1.0_wp)
      h_nconst_d(1)  = 0.0_wp
      h_nconst_d(2)  = 1500.0_wp
      h_nconst_d(3)  = 12000.0_wp
    N_nconst_d(:)  = 0.01_wp
      N_nconst_d(1)  = 0.001_wp
      N_nconst_d(2)  = 0.01_wp
      N_nconst_d(3)  = 0.02_wp
    rh_nconst_d(:)  = 0.0_wp     ! positive for decreasing humidity with height
      rhgr_nconst_d(1)  = 0.0_wp ! - (rh_poly_d(2)-rh_poly_d(1))/(h_poly_d(2)-h_poly_d(1))
      rhgr_nconst_d(2)  = 0.0_wp ! - (rh_poly_d(3)-rh_poly_d(2))/(h_poly_d(3)-h_poly_d(2))
      rhgr_nconst_d(3)  = 0.0_wp
    ldebug_artif_d  = .FALSE.
    idbg_artif_level_d  = 0_iintegers
    t_soil_c_d  = -1.0_wp
    t_water_c_d  = -1.0_wp
    t_snow_c_d  = -1.0_wp
    wf_soil_c_d  = 0.0_wp
    w_i_c_d  = 0.0_wp
    w_snow_c_d  = 0.0_wp
    itype_soil_tw_d  = 1_iintegers
    tsurffile_d(:) = ' ';   tsurffile_d  = 'noname.dat'
    hicefile_d(:) = ' ';   hicefile_d  = 'noname.dat'
    ssostdhfile_d(:) = ' ';   ssostdhfile_d  = 'noname.dat'
    ssogammafile_d(:) = ' ';   ssogammafile_d  = 'noname.dat'
    ssothetafile_d(:) = ' ';   ssothetafile_d  = 'noname.dat'
    ssosigmafile_d(:) = ' ';   ssosigmafile_d  = 'noname.dat'
    h_ice_c_d  = 0.1_wp
    t_ice_c_d  = 270.0_wp
    t_surf_c_d  = -1.0_wp
    lsensiflux_fix_d  = .FALSE.
    llatentflux_fix_d  = .FALSE.
    sensiflux_c_d  = 0.0_wp
    latentflux_LzuS_d  = 0.4_wp
    H0_rel_noise_d  = 0.0_wp
    iseed_noise_H0_d  = 608_iintegers
!!!adddefaulttag   used by the ed-script ed_addnl.in which automatically adds nlparams

    itype_artifprofiles_d = 2    ! analytic profiles for t/qv and u/v
    itype_anaprof_tqv_d   = 2    ! polytrope layers
    itype_anaprof_uv_d    = 3    ! const. U


    ! If .true., initialize vertical wind speed such that the wind vector 
    ! follows the terrain following coordinte surfaces:
    linitw_followeta_d = .FALSE.

    rasofile_d(:) = ' '; rasofile_d = 'noname.dat'

    ! Parameter for warm bubble(s):
    hcond_on_d = 0.0_wp              ![hours]



    !-------------------------------------------------------------------------------
    !- Section 2: Initialize variables with defaults
    !-------------------------------------------------------------------------------

    ! Coordinate and reference atmosphere settings
    ivctype       = ivctype_d
    irefatm       = irefatm_d
    p0sl          = p0sl_d
    t0sl          = t0sl_d
    dt0lp         = dt0lp_d
    delta_t       = delta_t_d
    h_scal        = h_scal_d
    bvref         = bvref_d

    ! Other parameters
    zspacing_type =      zspacing_type_d
    vcflat =             vcflat_d
    zz_top =             zz_top_d
    vcoordvec =          vcoordvec_d
    exp_galchen =        exp_galchen_d
    nfltvc =             nfltvc_d
    svc1 =               svc1_d
    svc2 =               svc2_d

    hill_width_x =       hill_width_x_d
    hill_width_y =       hill_width_y_d
    hillheight =         hillheight_d
    hill_rotangle =      hill_rotangle_d
    hillsideradius_y =   hillsideradius_y_d
    hillasym_x =         hillasym_x_d
    hillasym_y =         hillasym_y_d
    href_oro =           href_oro_d
    hill_i =             hill_i_d
    hill_j =             hill_j_d
    lhill =              lhill_d
    lhill_2d =           lhill_2d_d
    hill_type =          hill_type_d
    linitw_followeta =   linitw_followeta_d
    zo_boundary   =      zo_boundary_d
    exponent_windprof_boundary =  exponent_windprof_boundary_d
    linit_realoro =      linit_realoro_d
    orofile =            orofile_d
    i_shift_realoro  =   i_shift_realoro_d
    j_shift_realoro  =   j_shift_realoro_d
    z0file =             z0file_d
    frlandfile =         frlandfile_d
    soiltypefile =       soiltypefile_d
    plcovfile =          plcovfile_d
    laifile =            laifile_d
    rootdpfile =         rootdpfile_d
    forefile =           forefile_d
    fordfile =           fordfile_d
    tsoilfile  = tsoilfile_d
    tsnowfile  = tsnowfile_d
    wfsoilfile  = wfsoilfile_d
    wsnowfile  = wsnowfile_d
    wifile  = wifile_d
    hill_combineaction = hill_combineaction_d
    zhillcutfact =       zhillcutfact_d
    rasofile =           rasofile_d
    hcond_on =           hcond_on_d
    itype_artifprofiles =  itype_artifprofiles_d
    itype_anaprof_tqv =    itype_anaprof_tqv_d
    itype_anaprof_uv =     itype_anaprof_uv_d
    tkvhfix  = tkvhfix_d
    tkhhfix  = tkhhfix_d
    tkvmfix  = tkvmfix_d
    tkhmfix  = tkhmfix_d
    lnosurffluxes_m  = lnosurffluxes_m_d
    lnosurffluxes_h  = lnosurffluxes_h_d
    ltempdist  = ltempdist_d
    ctype_tempdist  = ctype_tempdist_d
    htempdist  = htempdist_d
    bub_centi  = bub_centi_d
    bub_centj  = bub_centj_d
    bub_centz  = bub_centz_d
    bub_timespan  = bub_timespan_d
    bub_radx  = bub_radx_d
    bub_rady  = bub_rady_d
    bub_radz  = bub_radz_d
    bub_rotangle  = bub_rotangle_d
    bub_dT  = bub_dT_d
    ladd_bubblenoise_t  = ladd_bubblenoise_t_d
    bub_dT_bubblenoise  = bub_dT_bubblenoise_d
    bub_zi_mcnider  = bub_zi_mcnider_d
    bub_zmax_mcnider  = bub_zmax_mcnider_d
    bub_h0_mcnider  = bub_h0_mcnider_d
    bub_heatingrate  = bub_heatingrate_d
    bub_centk  = bub_centk_d
    nlayers_poly  = nlayers_poly_d
    h_poly  = h_poly_d
    t_poly  = t_poly_d
    tgr_poly  = tgr_poly_d
    rh_poly  = rh_poly_d
    rhgr_poly  = rhgr_poly_d
    u_infty  = u_infty_d
    href_wk  = href_wk_d
    hmin_wk  = hmin_wk_d
    nlayers_linwind  = nlayers_linwind_d
    h_linwind  = h_linwind_d
    u_linwind  = u_linwind_d
    ugr_linwind  = ugr_linwind_d
    h_tropo_wk  = h_tropo_wk_d
    theta_0_wk  = theta_0_wk_d
    theta_tropo_wk  = theta_tropo_wk_d
    expo_theta_wk  = expo_theta_wk_d
    expo_relhum_wk  = expo_relhum_wk_d
    rh_min_wk  = rh_min_wk_d
    rh_max_wk  = rh_max_wk_d
    qv_max_wk  = qv_max_wk_d
    itype_topo  = itype_topo_d
    itype_soil_c  = itype_soil_c_d
    z0_c =  z0_c_d
    fr_land_c  = fr_land_c_d
    soiltyp_c  = soiltyp_c_d
    plcov_c  = plcov_c_d
    lai_c  = lai_c_d
    rootdp_c  = rootdp_c_d
    for_e_c  = for_e_c_d
    for_d_c  = for_d_c_d
    ladd_noise_t  = ladd_noise_t_d
    hadd_noise  = hadd_noise_d
    dT_noise  = dT_noise_d
    dW_noise  = dW_noise_d
    iseed_noise_t  = iseed_noise_t_d
    iseed_noise_w  = iseed_noise_w_d
    lps_from_file  = lps_from_file_d
    rasofile_t_is_theta  = rasofile_t_is_theta_d
    p_base_wk  = p_base_wk_d
    p_base_poly  = p_base_poly_d
    t_tropo_wk  = t_tropo_wk_d
    lbub_rhconst  = lbub_rhconst_d
    p_base_nconst  = p_base_nconst_d
    theta0_base_nconst  = theta0_base_nconst_d
    nlayers_nconst  = nlayers_nconst_d
    h_nconst  = h_nconst_d
    N_nconst  = N_nconst_d
    rh_nconst  = rh_nconst_d
    rhgr_nconst  = rhgr_nconst_d
    ldebug_artif  = ldebug_artif_d
    idbg_artif_level  = idbg_artif_level_d
    t_soil_c  = t_soil_c_d
    t_water_c = t_water_c_d
    t_snow_c  = t_snow_c_d
    wf_soil_c  = wf_soil_c_d
    w_i_c  = w_i_c_d
    w_snow_c  = w_snow_c_d
    itype_soil_tw  = itype_soil_tw_d
    tsurffile  = tsurffile_d
    hicefile  = hicefile_d
    ssostdhfile  = ssostdhfile_d
    ssogammafile  = ssogammafile_d
    ssothetafile  = ssothetafile_d
    ssosigmafile  = ssosigmafile_d
    h_ice_c  = h_ice_c_d
    t_ice_c  = t_ice_c_d
    t_surf_c  = t_surf_c_d
    lsensiflux_fix  = lsensiflux_fix_d
    llatentflux_fix  = llatentflux_fix_d
    sensiflux_c  = sensiflux_c_d
    latentflux_LzuS  = latentflux_LzuS_d
    H0_rel_noise  = H0_rel_noise_d
    iseed_noise_H0  = iseed_noise_H0_d
!!!setdefaulttag   used by the ed-script ed_addnl.in which automatically adds nlparams

    !-------------------------------------------------------------------------------
    !- Section 3: Input of the namelist values
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !
    !- Check, if the file for namelist ARTIFCTL exists and if it is
    !  mandatory to read:
    !
    !-------------------------------------------------------------------------------

    nlartiffile(:) = ' '
    nlartiffile = 'INPUT_IDEAL'

    IF (my_world_id == 0) THEN

      INQUIRE ( file=TRIM(ADJUSTL(nlartiffile)) , exist=nlfile_exists )
      IF ( nlfile_exists ) THEN
        WRITE (*,*) '    INPUT OF THE NAMELIST FOR (SEMI-)IDEALIZED RUNS'
      ELSE
        IF ( lartif_data ) THEN
          WRITE (*,*) 'THE NAMELIST FILE "INPUT_IDEAL" FOR THIS IDEALIZED RUN DOES NOT EXIST!'
          ierrstat = 23456
          RETURN
        ELSE
          !.. In this case, the namelist is not mandatory, so the model run should
          !   continue as normal. To enable this, we cannot just RETURN here,
          !   because we are only on node 0, so we set a certain error status here
          !   and jump back later on all processors.
          WRITE (*,'(a)') '    NAMELIST FOR (SEMI-)IDEALIZED RUNS IS NOT READ '// &
               'BECAUSE FILE "INPUT_IDEAL" DOES NOT EXIST'
          iz_err = 8765
        END IF
      END IF

      !===============================================================================

      IF ( nlfile_exists ) THEN

        iz_err  = 0
        io_errmsg(:) = ' '
        CALL get_free_unit(nuin)
        OPEN(nuin , FILE=TRIM(ADJUSTL(nlartiffile)), FORM='FORMATTED', STATUS='OLD',  &
             IOSTAT=iz_err, IOMSG=io_errmsg)
        IF(iz_err /= 0) THEN
          WRITE(*,*) ' ERROR    *** Error while opening file INPUT_IDEAL *** '
          WRITE(*,*) TRIM(io_errmsg)
          iz_err = -3
        ELSE
          READ (nuin, artifctl, IOSTAT=iz_err, IOMSG=io_errmsg)
          IF (iz_err > 0) THEN
            WRITE(*,*) ' ERROR    *** Error while reading NAMELIST /ARTIFCTL/, probably wrong values ***'
            WRITE(*,*) TRIM(io_errmsg)
            iz_err = -2
          ELSE IF (iz_err < 0) THEN
            WRITE(*,*) ' ERROR    *** Premature EOF encountered while reading NAMELIST /ARTIFCTL/ ***'
            WRITE(*,*) TRIM(io_errmsg)
           iz_err = -1
          ENDIF
          CLOSE(nuin)
          CALL release_unit(nuin)
        ENDIF

      END IF

      !-------------------------------------------------------------------------------
      !- Section 4: Check values for errors and consistency
      !-------------------------------------------------------------------------------

      IF  (iz_err == 0) THEN

        IF ( ivctype < 1 .OR. ivctype > 4 ) THEN
          WRITE (*,'(a,/,a)') ' ERROR  *** Wrong value for ivctype: ',ivctype,' *** ', &
               '                ***   must be >= 1 and <= 4 !!!          *** '
          iz_err = 1011
        ENDIF

        IF ( irefatm < 1 .OR. irefatm > 3 ) THEN
          WRITE (*,'(a,/,a)') ' ERROR  *** Wrong value for irefatm: ',irefatm,' *** ', &
               '        ***   must be >= 1 and <= 3 !!!             *** '
          iz_err = 1013
        ENDIF

        SELECT CASE ( TRIM(ADJUSTL(zspacing_type)) )

        CASE ('vcoordvec', 'VCOORDVEC', 'Vcoordvec')

          DO k=1,nvcoordvec_max
            IF ( vcoordvec(k) < -900.0_wp ) EXIT
            IF ( k > 1 .AND. ivctype == 1 .AND. vcoordvec(k)-vcoordvec(MAX(1,k-1)) <= 0.0_wp ) EXIT
            IF ( k > 1 .AND. ivctype == 2 .AND. vcoordvec(k)-vcoordvec(MAX(1,k-1)) >= 0.0_wp ) EXIT
            IF ( k > 1 .AND. ivctype == 3 .AND. vcoordvec(k)-vcoordvec(MAX(1,k-1)) >= 0.0_wp ) EXIT
            IF ( k > 1 .AND. ivctype == 4 .AND. vcoordvec(k)-vcoordvec(MAX(1,k-1)) >= 0.0_wp ) EXIT
          ENDDO
          IF ( vcoordvec(k) < -900.0_wp ) k=k-1
          IF ( k /= ke+1 ) THEN
            WRITE (*,*) ' ERROR  *** Wrong values for vcoordvec: ',vcoordvec(1:MIN(nvcoordvec_max,k)),' *** '
            WRITE (*,'(a,i3,a)') '        *** must contain ke+1 valid levels for ivctype = ',ivctype,' !!! *** '
            iz_err = 1012
          ENDIF
          
          IF (ivctype == 1) THEN
            IF (vcoordvec(1) <= 0.0_wp .OR. vcoordvec(ke+1) /= 1.0_wp) THEN
              WRITE (*,*) ' ERROR  *** Wrong values for vcoordvec: ',vcoordvec(1:MIN(nvcoordvec_max,k)),' *** '
              WRITE (*,'(a,i3,a)') '        *** must contain ke+1 increasing values ', &
                   'between >0.0 and 1.0 for ivctype = ',ivctype,' !!! *** '
              iz_err = 1013
            END IF
          END IF

          IF ( ANY( ivctype == (/2,3,4/) ) ) THEN
            IF (vcoordvec(ke+1) /= 0.0_wp) THEN
              WRITE (*,*) ' ERROR  *** Wrong values for vcoordvec: ',vcoordvec(1:MIN(nvcoordvec_max,k)),' *** '
              WRITE (*,'(a,i3,a)') '        *** must contain ke+1 decreasing values ', &
                   'between >0.0 and =0.0 for ivctype = ',ivctype,' !!! *** '
              iz_err = 1014
            END IF
          END IF

        END SELECT

        !.. Check settings of rh_xxx and rhgr_xxx to avoid negative relative humidities:
        !-----------------------------------------------------------------------------
        IF (itype_artifprofiles == 2 .OR. itype_artifprofiles == 3) THEN
          SELECT CASE (itype_anaprof_tqv)

          CASE (2)
            IF (nlayers_poly > nlayers_poly_max) THEN
              WRITE (*,*) ' ERROR  *** nlayers_poly > nlayers_poly_max !!! *** '
              iz_err = 1018
            ENDIF
            IF (ANY(h_poly(2:nlayers_poly) < h_poly(1:nlayers_poly-1))) THEN
              WRITE (*,*) ' ERROR  *** h_poly (k+1) < h_poly(k) occured (h_poly must be increasing) !!! *** '
              iz_err = 1019
            ENDIF
            IF (ANY(rh_poly(1:nlayers_poly-1) < 0.0_wp) .OR. &
                 ANY(rh_poly(1:nlayers_poly-1) - rhgr_poly(1:nlayers_poly-1)*&
                 (h_poly(2:nlayers_poly)-h_poly(1:nlayers_poly-1)) < 0.0_wp)) THEN
              WRITE (*,*) ' ERROR  *** wrong values of rh_poly and/or rhgr_poly lead to relhum < 0 !!! *** '
              iz_err = 1020
            ENDIF
            IF (h_poly(nlayers_poly) <= zz_top) THEN
              IF (rh_poly(nlayers_poly) < 0.0_wp .OR. &
                   rh_poly(nlayers_poly) - rhgr_poly(nlayers_poly)*&
                   (zz_top-h_poly(nlayers_poly)) < 0.0_wp) THEN
                WRITE (*,'(2a)') ' ERROR  *** wrong values of rh_poly and/or ', &
                     'rhgr_poly at the top lead to relhum < 0 !!! *** '
                iz_err = 1021
              ENDIF
            ENDIF

          CASE (3)
            IF (nlayers_nconst > nlayers_nconst_max) THEN
              WRITE (*,*) ' ERROR  *** nlayers_nconst > nlayers_nconst_max !!! *** '
              iz_err = 1022
            ENDIF
            IF (ANY(h_nconst(2:nlayers_nconst) < h_nconst(1:nlayers_nconst-1))) THEN
              WRITE (*,*) ' ERROR  *** h_nconst (k+1) < h_nconst(k) occured (h_nconst must be increasing) !!! *** '
              iz_err = 1019
            ENDIF
            IF (ANY(rh_nconst(1:nlayers_nconst-1) < 0.0_wp) .OR. &
                 ANY(rh_nconst(1:nlayers_nconst-1) - rhgr_nconst(1:nlayers_nconst-1)*&
                 (h_nconst(2:nlayers_nconst)-h_nconst(1:nlayers_nconst-1)) < 0.0_wp)) THEN
              WRITE (*,*) ' ERROR  *** wrong values of rh_nconst and/or rhgr_nconst lead to relhum < 0 !!! *** '
              iz_err = 1020
            ENDIF

            IF (h_nconst(nlayers_nconst) <= zz_top) THEN
              IF (rh_nconst(nlayers_nconst) < 0.0_wp .OR. &
                   rh_nconst(nlayers_nconst) - rhgr_nconst(nlayers_nconst)*&
                   (zz_top-h_nconst(nlayers_nconst)) < 0.0_wp) THEN
                WRITE (*,'(2a)') ' ERROR  *** wrong values of rh_nconst and/or ', &
                     'rhgr_nconst at the top lead to relhum < 0 !!! *** '
                iz_err = 1025
              ENDIF
            ENDIF

          END SELECT
        END IF

        !.. Check if zz_top or other height specifications are
        !   above the maximum altitude of certain idealized atmospheres:
        !   (checking is done based on assuming a dry atmosphere ...)
        !--------------------------------------------------------------------------------
        IF (itype_artifprofiles == 2 .OR. itype_artifprofiles == 3) THEN

          SELECT CASE(itype_anaprof_tqv)
          CASE (1)
            ! Nothing to check because the const. tropopause temperature
            ! ensures that the pressure does not become 0 anywhere above
            ! the tropopause height!
            ! HOWEVER, THERE MIGHT BE PROBLEMS IF THE TROPOPAUSE HEIGHT
            ! IS VERY HIGH AND AT VERY LOW TEMPERATURE,
            ! OR CONVERSELY IF THE ATMOSPHERE IS VERY HIGH AND HAS
            ! A HIGH TEMPERATURE AND TOO HIGH RELATIVE HUMIDITY, SO THAT THE 
            ! RESULTING VAPOR PRESSURE WOULD BE HIGHER THAN THE TOTAL PRESSURE!

          CASE (2)
            ! Layered polytrope atmosphere has finite height, if the
            ! temperature decreases with height!
            ! Check the layers until the second last:
            DO k=1, nlayers_poly-1
              IF (tgr_poly(k) > 0.0_wp) THEN
                IF (h_poly(k+1) > h_poly(k)+t_poly(k)/tgr_poly(k)) THEN
                  WRITE (*,'(2a,i4,a)') ' INPUT_ARTIFCTL: ERROR  *** Combination of h_poly, t_poly and tgr_poly',&
                       ' will lead to p <= 0.0 in polytrope layer k = ',k,' !!! *** '
                  iz_err = 1026
                END IF
              END IF
            END DO
            ! and here the uppermost layer:
            IF (tgr_poly(nlayers_poly) > 0.0_wp) THEN
              IF (zz_top > h_poly(nlayers_poly) + &
                   t_poly(nlayers_poly)/tgr_poly(nlayers_poly)) THEN
                WRITE (*,'(2a,i4,a)') ' INPUT_ARTIFCTL: ERROR  *** Combination of h_poly, t_poly and tgr_poly',&
                     ' will lead to p <= 0.0 in polytrope layer k = ',nlayers_poly,' !!! *** '
                iz_err = 1027
              END IF
            END IF

          CASE (3)
            ! Layered N-const. atmosphere has finite height, 
            ! IF N < N_Tconst (bvref_tconst), otherwise no finite height!
            ! The check is a little bit complicated here, because
            ! the finite height depends on T_0 and not on theata_0,
            ! but we need p_0 to compute T_0 from theta_0.
            ! To make it not too complicated, we disregard humidity here.

            ! Check the layers until the second last:
            zT0_nconst = theta0_base_nconst * (p_base_nconst/pt00)**rdocp
            DO k=1, nlayers_nconst-1
              bvref_tconst = SQRT(g**2 / (zT0_nconst*cp_d))
              IF (N_nconst(k) < bvref_tconst) THEN
                IF (ABS(N_nconst(k)) > 1.0E-12_wp) THEN
                  zmax_nconst = h_nconst(k) - g/N_nconst(k)**2 * &
                       LOG(1.0_wp-(N_nconst(k)/bvref_tconst)**2)
                ELSE
                  ! Limiting case for N = 0:
                  zmax_nconst = h_nconst(k) + g/bvref_tconst**2
                END IF
                IF (h_nconst(k+1) > zmax_nconst) THEN
                  WRITE (*,'(3a,i4,a)') ' INPUT_ARTIFCTL: ERROR  *** Combination of h_nconst,',&
                       ' N_nconst and theta0_base_nconst',&
                       ' will lead to p <= 0.0 in N-const. layer k = ',k,' !!! *** '
                  iz_err = 1026
                END IF
              END IF
              ! Calculate T_0 for the next layer based on the current layer values:
              zT0_nconst = T_z_N_const_dry( h_nconst(k+1), h_nconst(k), N_nconst(k), &
                   zT0_nconst )
            END DO
            ! and here the uppermost layer:
            bvref_tconst = SQRT(g**2 / (zT0_nconst*cp_d))
            k = nlayers_nconst
            IF (N_nconst(k) < bvref_tconst) THEN
              IF (ABS(N_nconst(k)) > 1.0E-12_wp) THEN
                zmax_nconst = h_nconst(k) - g/N_nconst(k)**2 * &
                     LOG(1.0_wp-(N_nconst(k)/bvref_tconst)**2)
              ELSE
                ! Limiting case for N = 0:
                zmax_nconst = h_nconst(k) + g/bvref_tconst**2
              END IF
              IF (zz_top > zmax_nconst) THEN
                WRITE (*,'(3a,i4,a)') ' INPUT_ARTIFCTL: ERROR  *** Combination of h_nconst,',&
                     ' N_nconst and theta0_base_nconst',&
                     ' will lead to p <= 0.0 in N-const. layer k = ',k,' !!! *** '
                iz_err = 1028
              END IF
            END IF

            ! Error check for exp(-N^2/g*(z-z0)) if N^2/g*(z-z0) > log(huge(1.0_wp)-1.0)
            zloghuge = LOG(HUGE(1.0_wp)-1.0_wp)
            DO k=1, nlayers_nconst-1
              bvref_tconst = N_nconst(k)**2 / g * (h_nconst(k+1)-h_nconst(k))
              IF (bvref_tconst >= zloghuge) THEN
                WRITE (*,'(2a,i4,a)') ' INPUT_ARTIFCTL: ERROR  *** Combination of h_nconst, N_nconst',&
                     ' will lead to floating overflow in layer k = ',k,' !!! *** '
                iz_err = 1030
              END IF
            END DO

            bvref_tconst = N_nconst(nlayers_nconst)**2 / &
                 g * (zz_top - h_nconst(nlayers_nconst))
            IF (bvref_tconst >= zloghuge) THEN
              WRITE (*,'(2a,i4,a)') ' INPUT_ARTIFCTL: ERROR  *** Combination of h_nconst, N_nconst',&
                   ' will lead to floating overflow in layer k = ',nlayers_nconst,' !!! *** '
              iz_err = 1030
            END IF

          END SELECT
        END IF

        !.. Check whether there are warm bubbles released before condensation and microphysics
        !   are switched on, i.e., if any(htempdist < hcond_on). If that is the case, issue a
        !   warning. The model run is continued, however, because there might be rare occasions
        !   where such a setting is desired by the user.
        IF ( ANY(htempdist < hcond_on) ) THEN
          WRITE (*,'(a,/,a)') '  WARNING  *** one or several warm bubble(s) start before condensation *** ', &
               '           *** is switched on, because any(htempdist < hcond_on)       *** '
        END IF

        ! Check setting of ctype_tempdist:
        numbubs_tot = 0
        DO i=1, ntempdist_max
          IF (ltempdist(i)) THEN
            bub_reg_check = .FALSE.
            IF (ANY(reg_tempdist == TRIM(ctype_tempdist(i)))) THEN
              bub_reg_check = .TRUE.
            END IF
            IF (ANY(reg_tempdist_tso == TRIM(ctype_tempdist(i)))) THEN
              bub_reg_check = .TRUE.
            END IF
            IF (ANY(reg_tempdist_bbc_ts == TRIM(ctype_tempdist(i)))) THEN
              bub_reg_check = .TRUE.
            END IF
            IF (ANY(reg_heatrate_dist_tso == TRIM(ctype_tempdist(i)))) THEN
              bub_reg_check = .TRUE.
            END IF
            IF (ANY(reg_heatrate_dist == TRIM(ctype_tempdist(i)))) THEN
              bub_reg_check = .TRUE.
            END IF

            IF ( bub_reg_check ) THEN
              numbubs_tot = numbubs_tot + 1
            ELSE
              WRITE (*,*) 'ERROR: in input_artifctl(): ctype_tempdist = '// &
                   TRIM(ctype_tempdist(i))//' unknown!'
              CALL model_abort (my_world_id, 10002, 'ctype_tempdist = '// &
                   TRIM(ctype_tempdist(i))//' unknown! ' // &
                   'modify namelist parameter ctype_tempdist and start again!', &
                   'input_artifctl(), ctype_tempdist')
            END IF
          END IF
        END DO

        !.. Check whether orography is read from a file and at the
        !   same time periodic BCs are applied. This might not be intended,
        !   and may lead to catastrophic results if the real orography
        !   is not periodically continuable, so give a warning and break here.
        !
        ! IF THE USER IS SURE WHAT HE IS DOING, HE MIGHT COMMENT OUT THIS
        ! ERROR CHECK AND CONTINUE:
        !
        IF (linit_realoro .AND. (lperi_x .OR. lperi_y)) THEN
          WRITE (*,'(4(a,/),a)') &
               '   ERROR: linit_realoro=.true. combined with periodic BCs', &
               '          might lead to catastrophic results if the orography', &
               '          is not periodically continuable!', &
               '   WE BREAK HERE!  IF YOU KNOW WHAT YOU ARE DOING, YOU MIGHT', &
               '   COMMENT OUT THIS ERROR CHECK IN SRC_ARTIFDATA.F90 AND RECOMPILE'
          CALL model_abort (my_world_id, 10347, &
               ' linit_realoro=.true. combined with periodic BCs', &
               'input_artifctl(), linit_realoro and (lperi_x or lperi_y)')
        END IF

        !.. Constant surface fluxes can only be specified if the soil model is
        !   turned off:
        IF (lsoil .AND. lsensiflux_fix) THEN
          WRITE (*,*) 'ERROR: in input_artifctl(): both lsensiflux_fix=.true. and lsoil=.true. not allowed!'
          iz_err = 1033
        END IF

        !.. Constant surface fluxes can only be specified if .not. (itype_turb=3 .and. itype_tran=3):
        IF (itype_turb == 3 .AND. itype_tran == 3 .AND. lsensiflux_fix) THEN
          WRITE (*,*) 'ERROR: in input_artifctl: lsensiflux_fix=.true. with itype_turb=3 .and. itype_tran=3 not working!'
          iz_err = 1033
        END IF

      ENDIF   ! IF  (iz_err == 0) THEN

    END IF  ! IF (my_world_id == 0) THEN

    IF (nproc > 1) THEN
      ! distribute error status to all processors
      CALL distribute_values  (iz_err, 1, 0, imp_integers,  icomm_world, ierr)
    ENDIF

    !------------------------------------------------------------------------------
    !- Section 4a: Check error status:
    !------------------------------------------------------------------------------

    IF (iz_err == 8765 .OR. (iz_err /= 0 .AND. .NOT.lartif_data) ) THEN
      ! The namelist is not mandatory and could not be read, so just continue.
      WRITE (*,*) '    NAMELIST FOR (SEMI-)IDEALIZED RUNS IS NOT MANDATORY, SO WE JUST CONTINUE.'
      ierrstat = 0
      RETURN
    ELSE IF (iz_err /= 0) THEN
      ierrstat = iz_err
      RETURN
    ENDIF


    !------------------------------------------------------------------------------
    !- Section 5: Distribute variables to all nodes
    !------------------------------------------------------------------------------

    IF ( nproc > 1 ) THEN

      charbuf(:)(1:lcbuf) = ' '

      IF (my_world_id == 0) THEN

        ! distribute scalars by the usual buffer vectors:
        nibuf = 0_iintegers
        nibuf = nibuf +  1;  intbuf  (nibuf) = ivctype
        nibuf = nibuf +  1;  intbuf  (nibuf) = irefatm
        nibuf = nibuf +  1;  intbuf  (nibuf) = nfltvc
        nibuf = nibuf +  1;  intbuf  (nibuf) = itype_artifprofiles
        nibuf = nibuf +  1;  intbuf  (nibuf) = itype_anaprof_tqv
        nibuf = nibuf +  1;  intbuf  (nibuf) = itype_anaprof_uv
        nibuf = nibuf +  1;  intbuf  (nibuf) = nlayers_poly
        nibuf = nibuf +  1;  intbuf  (nibuf) = nlayers_linwind
        nibuf = nibuf +  1;  intbuf  (nibuf) = i_shift_realoro
        nibuf = nibuf +  1;  intbuf  (nibuf) = j_shift_realoro
        nibuf = nibuf +  1;  intbuf  (nibuf) = itype_topo
        nibuf = nibuf +  1;  intbuf  (nibuf) = itype_soil_c
        nibuf = nibuf +  1;  intbuf  (nibuf) = nlayers_nconst
        nibuf = nibuf +  1;  intbuf  (nibuf) = idbg_artif_level
        nibuf = nibuf +  1;  intbuf  (nibuf) = itype_soil_tw
        nibuf = nibuf +  1;  intbuf  (nibuf) = iseed_noise_t
        nibuf = nibuf +  1;  intbuf  (nibuf) = iseed_noise_w
        nibuf = nibuf +  1;  intbuf  (nibuf) = iseed_noise_H0
!!!ibuftag   used by the ed-script ed_addnl.in which automatically adds nlparams

        nrbuf = 0_iintegers
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = p0sl
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = t0sl
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = dt0lp
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = delta_t
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = h_scal
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = bvref
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = exp_galchen
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = zz_top
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = vcflat
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = href_oro
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = hcond_on
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = zo_boundary
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = exponent_windprof_boundary
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = svc1
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = svc2
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = tkvhfix
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = tkhhfix
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = tkvmfix
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = tkhmfix
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = u_infty
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = href_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = hmin_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = h_tropo_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = theta_0_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = theta_tropo_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = expo_theta_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = expo_relhum_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = rh_min_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = rh_max_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = qv_max_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = z0_c
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = fr_land_c
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = soiltyp_c
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = plcov_c
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = lai_c
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = rootdp_c
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = for_e_c
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = for_d_c
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = hadd_noise
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = dT_noise
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = dW_noise
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = p_base_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = p_base_poly
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = t_tropo_wk
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = p_base_nconst
        nrbuf = nrbuf +  1;  realbuf (nrbuf) = theta0_base_nconst
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = t_soil_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = t_water_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = t_snow_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = wf_soil_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = w_i_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = w_snow_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = h_ice_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = t_ice_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = t_surf_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = sensiflux_c
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = latentflux_LzuS
        nrbuf = nrbuf +  1;  realbuf  (nrbuf) = H0_rel_noise
!!!realbuftag   used by the ed-script ed_addnl.in which automatically adds nlparams

        ncbuf = 0_iintegers
        ncbuf = ncbuf +  1;  charbuf (ncbuf)(1:lcbuf) = zspacing_type
!!!charbuftag   used by the ed-script ed_addnl.in which automatically adds nlparams

        nlbuf = 0_iintegers
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = linitw_followeta
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = linit_realoro
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = lnosurffluxes_m
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = lnosurffluxes_h
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = ladd_noise_t
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = lps_from_file
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = rasofile_t_is_theta
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = ldebug_artif
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = lsensiflux_fix
        nlbuf = nlbuf +  1;  logbuf  (nlbuf) = llatentflux_fix
!!!logbuftag   used by the ed-script ed_addnl.in which automatically adds nlparams

      ENDIF

      ! First, distribute the buffer sizes:
      CALL distribute_values (nrbuf, 1, 0, imp_integers,     icomm_world, ierr)
      CALL distribute_values (nibuf, 1, 0, imp_integers,     icomm_world, ierr)
      CALL distribute_values (nlbuf, 1, 0, imp_integers,     icomm_world, ierr)
      CALL distribute_values (ncbuf, 1, 0, imp_integers,     icomm_world, ierr)

      ! Then distribute the actual buffers:
      CALL distribute_values (realbuf, nrbuf, 0, imp_reals,     icomm_world, ierr)
      CALL distribute_values (intbuf,  nibuf, 0, imp_integers,  icomm_world, ierr)
      CALL distribute_values (logbuf,  nlbuf, 0, imp_logical,   icomm_world, ierr)
      CALL distribute_values (charbuf, ncbuf, 0, imp_character, icomm_world, ierr)

      IF (my_world_id /= 0) THEN

        nibuf = 0_iintegers
        nibuf = nibuf +  1;  ivctype    = intbuf  (nibuf)
        nibuf = nibuf +  1;  irefatm    = intbuf  (nibuf)
        nibuf = nibuf +  1;  nfltvc     = intbuf  (nibuf)
        nibuf = nibuf +  1;  itype_artifprofiles = intbuf  (nibuf)
        nibuf = nibuf +  1;  itype_anaprof_tqv   = intbuf  (nibuf)
        nibuf = nibuf +  1;  itype_anaprof_uv    = intbuf  (nibuf)
        nibuf = nibuf +  1;  nlayers_poly    = intbuf  (nibuf)
        nibuf = nibuf +  1;  nlayers_linwind    = intbuf  (nibuf)
        nibuf = nibuf +  1;  i_shift_realoro    = intbuf  (nibuf)
        nibuf = nibuf +  1;  j_shift_realoro    = intbuf  (nibuf)
        nibuf = nibuf +  1;  itype_topo    = intbuf  (nibuf)
        nibuf = nibuf +  1;  itype_soil_c    = intbuf  (nibuf)
        nibuf = nibuf +  1;  nlayers_nconst    = intbuf  (nibuf)
        nibuf = nibuf +  1;  idbg_artif_level    = intbuf  (nibuf)
        nibuf = nibuf +  1;  itype_soil_tw    = intbuf  (nibuf)
        nibuf = nibuf +  1;  iseed_noise_t    = intbuf  (nibuf)
        nibuf = nibuf +  1;  iseed_noise_w    = intbuf  (nibuf)
        nibuf = nibuf +  1;  iseed_noise_H0    = intbuf  (nibuf)
!!!ibufbacktag   used by the ed-script ed_addnl.in which automatically adds nlparams

        nrbuf = 0_iintegers
        nrbuf = nrbuf +  1;  p0sl        = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  t0sl        = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  dt0lp       = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  delta_t     = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  h_scal      = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  bvref       = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  exp_galchen = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  zz_top      = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  vcflat      = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  href_oro    = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  hcond_on    = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  zo_boundary = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  exponent_windprof_boundary = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  svc1       = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  svc2       = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  tkvhfix    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  tkhhfix    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  tkvmfix    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  tkhmfix    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  u_infty    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  href_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  hmin_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  h_tropo_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  theta_0_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  theta_tropo_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  expo_theta_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  expo_relhum_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  rh_min_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  rh_max_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  qv_max_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  z0_c     = realbuf (nrbuf)
        nrbuf = nrbuf +  1;  fr_land_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  soiltyp_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  plcov_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  lai_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  rootdp_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  for_e_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  for_d_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  hadd_noise    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  dT_noise    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  dW_noise    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  p_base_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  p_base_poly    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  t_tropo_wk    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  p_base_nconst    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  theta0_base_nconst    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  t_soil_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  t_water_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  t_snow_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  wf_soil_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  w_i_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  w_snow_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  h_ice_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  t_ice_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  t_surf_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  sensiflux_c    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  latentflux_LzuS    = realbuf  (nrbuf)
        nrbuf = nrbuf +  1;  H0_rel_noise    = realbuf  (nrbuf)
!!!realbufbacktag   used by the ed-script ed_addnl.in which automatically adds nlparams

        ncbuf = 0_iintegers
        ncbuf = ncbuf +  1;  zspacing_type = charbuf (ncbuf)(1:lcbuf)
!!!charbufbacktag   used by the ed-script ed_addnl.in which automatically adds nlparams

        nlbuf = 0_iintegers
        nlbuf = nlbuf +  1;  linitw_followeta = logbuf (nlbuf)
        nlbuf = nlbuf +  1;  linit_realoro    = logbuf (nlbuf)
        nlbuf = nlbuf +  1;  lnosurffluxes_m  = logbuf  (nlbuf)
        nlbuf = nlbuf +  1;  lnosurffluxes_h  = logbuf  (nlbuf)
        nlbuf = nlbuf +  1;  ladd_noise_t     = logbuf  (nlbuf)
        nlbuf = nlbuf +  1;  lps_from_file    = logbuf  (nlbuf)
        nlbuf = nlbuf +  1;  rasofile_t_is_theta    = logbuf  (nlbuf)
        nlbuf = nlbuf +  1;  ldebug_artif    = logbuf  (nlbuf)
        nlbuf = nlbuf +  1;  lsensiflux_fix    = logbuf  (nlbuf)
        nlbuf = nlbuf +  1;  llatentflux_fix    = logbuf  (nlbuf)
!!!logbufbacktag   used by the ed-script ed_addnl.in which automatically adds nlparams

      ENDIF

      ! distribute fields of length nhill_max directly without using buffers:
      CALL distribute_values(zhillcutfact, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hill_i, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hill_j, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hill_width_x, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hill_width_y, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hillsideradius_y, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hillheight, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hill_rotangle, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hillasym_x, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hillasym_y, nhill_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(hill_combineaction, nhill_max, 0, imp_integers, icomm_world, ierr)
      CALL distribute_values(lhill, nhill_max, 0, imp_logical, icomm_world, ierr)
      CALL distribute_values(lhill_2d, nhill_max, 0, imp_logical, icomm_world, ierr)

      CALL distribute_values(hill_type, nhill_max, 0,imp_character, icomm_world, ierr)

      CALL distribute_path(orofile)
      CALL distribute_path(z0file)
      CALL distribute_path(frlandfile)
      CALL distribute_path(soiltypefile)
      CALL distribute_path(plcovfile)
      CALL distribute_path(laifile)
      CALL distribute_path(rootdpfile)
      CALL distribute_path(forefile)
      CALL distribute_path(fordfile)
      CALL distribute_path(tsoilfile)
      CALL distribute_path(tsnowfile)
      CALL distribute_path(wfsoilfile)
      CALL distribute_path(wsnowfile)
      CALL distribute_path(wifile)
      CALL distribute_path(rasofile)
      CALL distribute_values(vcoordvec, nvcoordvec_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(ltempdist, ntempdist_max, 0, imp_logical, icomm_world, ierr)

      CALL distribute_values(ctype_tempdist, ntempdist_max, 0, imp_character, icomm_world, ierr)

      CALL distribute_values(htempdist, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_centi, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_centj, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_centz, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_timespan, ntempdist_max, 0, imp_integers, icomm_world, ierr)
      CALL distribute_values(bub_radx, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_rady, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_radz, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_rotangle, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_dT, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(ladd_bubblenoise_t, ntempdist_max, 0, imp_logical, icomm_world, ierr)
      CALL distribute_values(bub_dT_bubblenoise, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_zi_mcnider, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_zmax_mcnider, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_h0_mcnider, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_heatingrate, ntempdist_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(bub_centk, ntempdist_max, 0, imp_integers, icomm_world, ierr)
      CALL distribute_values(h_poly, nlayers_poly_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(t_poly, nlayers_poly_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(tgr_poly, nlayers_poly_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(rh_poly, nlayers_poly_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(rhgr_poly, nlayers_poly_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(h_linwind, nlayers_linwind_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(u_linwind, nlayers_linwind_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(ugr_linwind, nlayers_linwind_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(lbub_rhconst, ntempdist_max, 0, imp_logical, icomm_world, ierr)
      CALL distribute_values(h_nconst, nlayers_nconst_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(N_nconst, nlayers_nconst_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(rh_nconst, nlayers_nconst_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_values(rhgr_nconst, nlayers_nconst_max, 0, imp_reals, icomm_world, ierr)
      CALL distribute_path(tsurffile)
      CALL distribute_path(hicefile)
      CALL distribute_path(ssostdhfile)
      CALL distribute_path(ssogammafile)
      CALL distribute_path(ssothetafile)
      CALL distribute_path(ssosigmafile)
!!!distribbacktag   used by the ed-script ed_addnl.in which automatically adds nlparams

    ENDIF  ! IF (. nproc > 1)

    !------------------------------------------------------------------------------
    ! Calculate some derived variables from namelist parameters:
    !------------------------------------------------------------------------------

    ! Time step of first occurence of temperature disturbance(s):
    ntstep_bubble(:) = INT(htempdist(:) * 3600 / dt)
    ! Time step of adding noise to the temperature field:
    ntstep_noise  = INT(hadd_noise * 3600 / dt)

    ! Total number of temperature- and heatingrate disturbances
    ! and ordering index of each defined disturbance within the 
    ! temperature- and noisebuffer fields
    numbubs_tot = 0_iintegers    ! computed again here on all PEs
    numbubs_noise = 0_iintegers
    bub_index_noise(:) = 0_iintegers
    DO i=1, ntempdist_max
      IF (ltempdist(i)) THEN
        numbubs_tot = numbubs_tot + 1
        IF (ladd_bubblenoise_t(i)) THEN
          numbubs_noise = numbubs_noise + 1
          bub_index_noise(i) = numbubs_noise
        END IF
      END IF
    END DO

    !------------------------------------------------------------------------------
    !- Section 6: Output of the namelist variables and their default values
    !------------------------------------------------------------------------------

    IF ( my_world_id == 0 ) THEN

      WRITE (nuspecif, '(A2)')  '  '
      WRITE (nuspecif, '(A23)') '0    NAMELIST: artifctl'
      WRITE (nuspecif, '(A23)') '     ------------------'
      WRITE (nuspecif, '(A2)')  '  '
      WRITE (nuspecif, '(A)')  'Variables for idealized simulations'
      WRITE (nuspecif, '(T8,A,T33,A,T52,A,T71,A)') 'Variable', 'Actual Value',   &
           'Default Value', 'Format'
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A)')                       &
           'ldebug_artif', ldebug_artif, ldebug_artif_d     , ' L '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'idbg_artif_level', idbg_artif_level, idbg_artif_level_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A)')                       &
           'ivctype'     ,ivctype  , ivctype_d       , ' I '
      WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A)')                       &
           'irefatm'     ,irefatm  , irefatm_d       , ' I '
      WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2,T71,A)')                       &
           'p0sl'   ,p0sl    , p0sl_d       , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           't0sl'   ,t0sl    , t0sl_d       , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'dt0lp'  ,dt0lp   , dt0lp_d      , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'delta_t',delta_t , delta_t_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.2,T52,F12.2,T71,A)')                       &
           'h_scal' ,h_scal  , h_scal_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'bvref', bvref, bvref_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'vcflat',vcflat , vcflat_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'zz_top',zz_top , zz_top_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'exp_galchen',exp_galchen , exp_galchen_d     , ' R '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'zspacing_type',TRIM(zspacing_type) , TRIM(zspacing_type_d), ' C '
      WRITE (nuspecif, '(T8,A,T33,I12  ,T52,I12  ,T71,A3)')                      &
           'nfltvc',nfltvc,nfltvc_d,' I '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
           'svc1',svc1,svc1_d,' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A3)')                      &
           'svc2',svc2,svc2_d,' R '
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A3)')                      &
           'linitw_followeta', linitw_followeta, linitw_followeta_d    ,' L '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'href_oro',href_oro , href_oro_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A3)')                      &
           'linit_realoro', linit_realoro, linit_realoro_d    ,' L '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'i_shift_realoro', i_shift_realoro, i_shift_realoro_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'j_shift_realoro', j_shift_realoro, j_shift_realoro_d     , ' I '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'orofile',     TRIM(orofile) , TRIM(orofile_d) , ' C '    
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'z0file',      TRIM(z0file) , TRIM(z0file_d) , ' C '    
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'frlandfile',  TRIM(frlandfile) , TRIM(frlandfile_d)  , ' C '   
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'soiltypefile',TRIM(soiltypefile) , TRIM(soiltypefile_d) , ' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'plcovfile',   TRIM(plcovfile) , TRIM(plcovfile_d)  , ' C '   
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'laifile',     TRIM(laifile) , TRIM(laifile_d)  , ' C '   
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'rootdpfile',  TRIM(rootdpfile) , TRIM(rootdpfile_d) , ' C '    
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'forefile',    TRIM(forefile) , TRIM(forefile_d) , ' C '    
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'fordfile',    TRIM(fordfile) , TRIM(fordfile_d) , ' C '    
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'tsurffile',   TRIM(tsurffile) ,TRIM(tsurffile_d), ' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'hicefile',    TRIM(hicefile) , TRIM(hicefile_d), ' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'ssostdhfile', TRIM(ssostdhfile) , TRIM(ssostdhfile_d), ' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'ssogammafile',TRIM(ssogammafile) , TRIM(ssogammafile_d), ' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'ssothetafile',TRIM(ssothetafile) , TRIM(ssothetafile_d), ' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'ssosigmafile',TRIM(ssosigmafile) , TRIM(ssosigmafile_d), ' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'rasofile',    TRIM(rasofile) , TRIM(rasofile_d)   , ' C '  
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'itype_artifprofiles', itype_artifprofiles, itype_artifprofiles_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'itype_anaprof_tqv', itype_anaprof_tqv, itype_anaprof_tqv_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'itype_anaprof_uv', itype_anaprof_uv, itype_anaprof_uv_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'tkvhfix', tkvhfix, tkvhfix_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'tkhhfix', tkhhfix, tkhhfix_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'tkvmfix', tkvmfix, tkvmfix_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'tkhmfix', tkhmfix, tkhmfix_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A)')                       &
           'lnosurffluxes_m', lnosurffluxes_m, lnosurffluxes_m_d     , ' L '
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A)')                       &
           'lnosurffluxes_h', lnosurffluxes_h, lnosurffluxes_h_d     , ' L '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'nlayers_poly', nlayers_poly, nlayers_poly_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'p_base_poly', p_base_poly, p_base_poly_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'u_infty', u_infty, u_infty_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'href_wk', href_wk, href_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'hmin_wk', hmin_wk, hmin_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'p_base_wk', p_base_wk, p_base_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'nlayers_linwind', nlayers_linwind, nlayers_linwind_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'h_tropo_wk', h_tropo_wk, h_tropo_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'theta_0_wk', theta_0_wk, theta_0_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'theta_tropo_wk', theta_tropo_wk, theta_tropo_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           't_tropo_wk', t_tropo_wk, t_tropo_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'expo_theta_wk', expo_theta_wk, expo_theta_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'expo_relhum_wk', expo_relhum_wk, expo_relhum_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'rh_min_wk', rh_min_wk, rh_min_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'rh_max_wk', rh_max_wk, rh_max_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'qv_max_wk', qv_max_wk, qv_max_wk_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'itype_topo', itype_topo, itype_topo_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'itype_soil_c', itype_soil_c, itype_soil_c_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'z0_c',z0_c , z0_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'fr_land_c', fr_land_c, fr_land_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'soiltyp_c', soiltyp_c, soiltyp_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'plcov_c', plcov_c, plcov_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'lai_c', lai_c, lai_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'rootdp_c', rootdp_c, rootdp_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'for_e_c', for_e_c, for_e_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'for_d_c', for_d_c, for_d_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           't_soil_c', t_soil_c, t_soil_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           't_water_c', t_water_c, t_water_c_d  , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           't_snow_c', t_snow_c, t_snow_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'wf_soil_c', wf_soil_c, wf_soil_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'w_i_c', w_i_c, w_i_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'w_snow_c', w_snow_c, w_snow_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'h_ice_c', h_ice_c, h_ice_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           't_ice_c', t_ice_c, t_ice_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A)')                       &
           'ladd_noise_t', ladd_noise_t, ladd_noise_t_d     , ' L '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'hadd_noise', hadd_noise, hadd_noise_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'dT_noise', dT_noise, dT_noise_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'dW_noise', dW_noise, dW_noise_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'iseed_noise_t', iseed_noise_t, iseed_noise_t_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'iseed_noise_w', iseed_noise_w, iseed_noise_w_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A)')                       &
           'lps_from_file', lps_from_file, lps_from_file_d     , ' L '
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A)')                       &
           'rasofile_t_is_theta', rasofile_t_is_theta, rasofile_t_is_theta_d     , ' L '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'p_base_nconst', p_base_nconst, p_base_nconst_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'theta0_base_nconst', theta0_base_nconst, theta0_base_nconst_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'nlayers_nconst', nlayers_nconst, nlayers_nconst_d     , ' I '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'itype_soil_tw', itype_soil_tw, itype_soil_tw_d     , ' I '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'tsoilfile', TRIM(tsoilfile) , TRIM(tsoilfile_d),' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'tsnowfile', TRIM(tsnowfile) , TRIM(tsnowfile_d),' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'wfsoilfile',TRIM(wfsoilfile) , TRIM(wfsoilfile_d),' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'wsnowfile', TRIM(wsnowfile) , TRIM(wsnowfile_d),' C '
      WRITE (nuspecif, '(T8,A,/,T33,A,/,T33,A,T71,A)')                       &
           'wifile',    TRIM(wifile) , TRIM(wifile_d),' C '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                   &
           't_surf_c', t_surf_c, t_surf_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A)')                       &
           'lsensiflux_fix', lsensiflux_fix, lsensiflux_fix_d     , ' L '
      WRITE (nuspecif, '(T8,A,T33,L12,T52,L12,T71,A)')                       &
           'llatentflux_fix', llatentflux_fix, llatentflux_fix_d     , ' L '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'sensiflux_c', sensiflux_c, sensiflux_c_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'latentflux_LzuS', latentflux_LzuS, latentflux_LzuS_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'H0_rel_noise', H0_rel_noise, H0_rel_noise_d     , ' R '
      WRITE (nuspecif, '(T8,A,T33,I12,T52,I12,T71,A)')                       &
           'iseed_noise_H0', iseed_noise_H0, iseed_noise_H0_d     , ' I '
!!!nuspeciftag   used by the ed-script ed_addnl.in which automatically adds nlparams

      DO k=1,ke+1
        IF (k == 1 .OR. vcoordvec(k) >= -900.0_wp) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")",T33,F12.4,T52,F12.4,T71,A3)')                    &
               'vcoordvec',k,vcoordvec(k),vcoordvec_d(k),' R '
        END IF
      ENDDO

      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")",T33,L12,T52,L12,T71,A3)')                      &
               'lhill', i, lhill(i), lhill_d(i)    ,' L '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")",T33,L12,T52,L12,T71,A3)')                      &
               'lhill_2d', i, lhill_2d(i), lhill_2d_d(i)    ,' L '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")",/,T33,A,/,T33,A,T71,A)')                       &
               'hill_type',i,TRIM(ADJUSTL(hill_type(i))) , TRIM(ADJUSTL(hill_type_d(i))),' C '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,f12.4,T52,f12.4,T71,A)')                       &
               'hill_i',i,hill_i(i) , hill_i_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,f12.4,T52,f12.4,T71,A)')                       &
               'hill_j',i,hill_j(i) , hill_j_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'hill_width_x',i,hill_width_x(i) , hill_width_x_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'hill_width_y',i,hill_width_y(i) , hill_width_y_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'hillsideradius_y',i,hillsideradius_y(i) , hillsideradius_y_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'hillheight',i,hillheight(i) , hillheight_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'hill_rotangle',i,hill_rotangle(i) , hill_rotangle_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'zhillcutfact',i,zhillcutfact(i) , zhillcutfact_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'hillasym_x',i,hillasym_x(i) , hillasym_x_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'hillasym_y',i,hillasym_y(i) , hillasym_y_d(i)     , ' R '
        END IF
      END DO
      DO i=1,nhill_max
        IF (lhill(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,i8,T52,i8,T71,A)')                       &
               'hill_combineaction',i,hill_combineaction(i) , hill_combineaction_d(i)     , ' I '
        END IF
      END DO
      WRITE (nuspecif, '(T8,A,T33,F12.4,T52,F12.4,T71,A)')                       &
           'hcond_on',hcond_on , hcond_on_d     , ' R '
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,L12,T52,L12,T71,A)')                       &
               'ltempdist',i,ltempdist(i) , ltempdist_d(i)     , ' L '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")",/,T33,A,/,T33,A,T71,A)')                       &
               'ctype_tempdist',i,TRIM(ADJUSTL(ctype_tempdist(i))) , TRIM(ADJUSTL(ctype_tempdist_d(i))),' C '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'htempdist',i,htempdist(i) , htempdist_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_centi',i,bub_centi(i) , bub_centi_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_centj',i,bub_centj(i) , bub_centj_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_centz',i,bub_centz(i) , bub_centz_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,i8,T52,i8,T71,A)')                       &
               'bub_timespan',i,bub_timespan(i) , bub_timespan_d(i)     , ' I '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_radx',i,bub_radx(i) , bub_radx_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_rady',i,bub_rady(i) , bub_rady_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_radz',i,bub_radz(i) , bub_radz_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_rotangle',i,bub_rotangle(i) , bub_rotangle_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_dT',i,bub_dT(i) , bub_dT_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,L12,T52,L12,T71,A)')                       &
               'ladd_bubblenoise_t',i,ladd_bubblenoise_t(i) , ladd_bubblenoise_t_d(i)     , ' L '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_dT_bubblenoise',i,bub_dT_bubblenoise(i) , bub_dT_bubblenoise_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_zi_mcnider',i,bub_zi_mcnider(i) , bub_zi_mcnider_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_zmax_mcnider',i,bub_zmax_mcnider(i) , bub_zmax_mcnider_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_h0_mcnider',i,bub_h0_mcnider(i) , bub_h0_mcnider_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
               'bub_heatingrate',i,bub_heatingrate(i) , bub_heatingrate_d(i)     , ' R '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,i8,T52,i8,T71,A)')                       &
               'bub_centk',i,bub_centk(i) , bub_centk_d(i)     , ' I '
        END IF
      END DO
      DO i=1,ntempdist_max
        IF (ltempdist(i)) THEN
          WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,L12,T52,L12,T71,A)')                       &
               'lbub_rhconst',i,lbub_rhconst(i) , lbub_rhconst_d(i)     , ' L '
        END IF
      END DO
      DO i=1,nlayers_poly
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
             'h_poly',i,h_poly(i) , h_poly_d(i)     , ' R '
      END DO
      DO i=1,nlayers_poly
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
             't_poly',i,t_poly(i) , t_poly_d(i)     , ' R '
      END DO
      DO i=1,nlayers_poly
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,ES12.5,T52,ES12.5,T71,A)')                       &
             'tgr_poly',i,tgr_poly(i) , tgr_poly_d(i)     , ' R '
      END DO
      DO i=1,nlayers_poly
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
             'rh_poly',i,rh_poly(i) , rh_poly_d(i)     , ' R '
      END DO
      DO i=1,nlayers_poly
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,ES12.5,T52,ES12.5,T71,A)')                       &
             'rhgr_poly',i,rhgr_poly(i) , rhgr_poly_d(i)     , ' R '
      END DO
      DO i=1,nlayers_linwind
        IF (h_linwind_d(i) == HUGE(1.0_wp)) h_linwind_d(i) = -999.99_wp
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
             'h_linwind',i,h_linwind(i) , h_linwind_d(i)     , ' R '
      END DO
      DO i=1,nlayers_linwind
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
             'u_linwind',i,u_linwind(i) , u_linwind_d(i)     , ' R '
      END DO
      DO i=1,nlayers_linwind
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,ES12.5,T52,ES12.5,T71,A)')                       &
             'ugr_linwind',i,ugr_linwind(i) , ugr_linwind_d(i)     , ' R '
      END DO
      DO i=1,nlayers_nconst
        IF (h_nconst_d(i) == HUGE(1.0_wp)) h_nconst_d(i) = -999.99_wp
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,F12.4,T52,F12.4,T71,A)')                       &
             'h_nconst',i,h_nconst(i) , h_nconst_d(i)     , ' R '
      END DO
      DO i=1,nlayers_nconst
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,ES12.5,T52,ES12.5,T71,A)')                       &
             'N_nconst',i,N_nconst(i) , N_nconst_d(i)     , ' R '
      END DO
      DO i=1,nlayers_nconst
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,ES12.5,T52,ES12.5,T71,A)')                       &
             'rh_nconst',i,rh_nconst(i) , rh_nconst_d(i)     , ' R '
      END DO
      DO i=1,nlayers_nconst
        WRITE (nuspecif, '(T8,A,"(",i3.3,")", T33,ES12.5,T52,ES12.5,T71,A)')                       &
             'rhgr_nconst',i,rhgr_nconst(i) , rhgr_nconst_d(i)     , ' R '
      END DO
!!!nuspecif2tag   used by the ed-script ed_addnl.in which automatically adds nlparams

    ENDIF

    RETURN

  CONTAINS

    ! This is a hack for distribute_values() of path variables with a length > 100 characters,
    !   because distribute_values() does only transfer the first 100 characters.
    SUBROUTINE distribute_path ( pathstring )

      CHARACTER(len=*), INTENT(inout):: pathstring
      INTEGER(kind=iintegers)        :: plen, i, badgelen, cs, ce

      badgelen = 100
      plen = LEN(pathstring) ! NOT len_trim(), because plen needs to be the same value on all processors!

      DO i=1, plen/badgelen+1
        cs = (i-1) * badgelen + 1
        ce = MIN(cs + badgelen - 1, plen)
        CALL distribute_values(pathstring(cs:ce), 100, 0, imp_character, icomm_world, ierr)
      END DO

    END SUBROUTINE distribute_path

  END SUBROUTINE input_artifctl


  !+ Module procedure in "src_artifdata" for generating artificial initial data

  SUBROUTINE gen_ini_data 

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   Artificial initial data are generated, if lartif_data = .TRUE.
    !   The data describe a hydrostatic layered atmosphere disturbed by an 
    !   artificial orography, the bell-shaped mountain, and / or a 
    !   non-balanced deviation of the temperature. The subroutine is written in a
    !   way that it serves for both environments, the parallel and the sequential.
    !
    !   The initial data are written on timelevel nnew (which becomes nnow in
    !   the beginning of the time loop). In case of leapfrog integration, the
    !   data are copied to limelevel nnow (which becomes nold in the beginning
    !   of the time loop).
    !
    ! Method:
    !   - Initialization of the vertical coordinate parameter
    !   - Computation of the artificial orography (bell-shaped mountain)
    !   - Initialization of constant fields (ozone content, soil description,...)
    !   - Computation of the constant reference atmosphere
    !   - Initialization of the wind fields
    !   - Initialization of thermodynamic variables
    !   - Initialization of soil variables
    !
    !------------------------------------------------------------------------------
    !
    ! Declarations:
    ! Local variables

    INTEGER (KIND=iintegers)   ::       &
         i, j, k, izerror, kzdims(24), &
         ierrstat, iztrcr

    CHARACTER (LEN= 80)        ::  yzroutine
    CHARACTER (LEN=250)        ::  yzerrmsg

    REAL (KIND=wp)             ::       &
         zvco_p20(21),      & ! pre-specified vertical coordinates of 20 p-layer LM
         zvco_p32(33),      & ! pre-specified vertical coordinates of 32 p-layer LM
         zvco_p35(36),      & ! pre-specified vertical coordinates of 35 p-layer LM
         zvco_p40(41),      & ! pre-specified vertical coordinates of 40 p-layer LM
         zvco_z20(21),      & ! pre-specified vertical coordinates of 20 z-layer LM
         zvco_z32(33),      & ! pre-specified vertical coordinates of 32 z-layer LM
         zvco_z35(36),      & ! pre-specified vertical coordinates of 35 z-layer LM
         zvco_z40(41),      & ! pre-specified vertical coordinates of 40 z-layer LM
         zvco_z50(51),      & ! CPS
         zrhf   (ie,je,ke), & !
         zrhs   (ie,je),    & !
         zml    (ie,je,ke), & ! height of main levels for mass points
         zml_u    (ie,je,ke), & ! height of main levels for u-points
         zml_v    (ie,je,ke), & ! height of main levels for v-points
         zak(ke1), zbk(ke1)

    REAL (KIND=wp)             ::       &
         zmaxqv, zrdm, zcpm

    LOGICAL                    ::       &
         ltheta_ini              ! has to be set to .true. if theta_ini is to be used instead of t

    REAL (KIND=wp),     ALLOCATABLE         ::       &
         hsurfs(:,:,:),    & ! Splitted topographies (SLEVE coordinate)
         hsurfs_tot(:,:,:),& ! Splitted topographies of full domain (SLEVE)
         hsurf_tot(:,:)      ! Full topography of full domain

    REAL (KIND=wp)             :: bvref_tconst,  & ! BV-freq. associated with const. T reference atmosphere
                                  hmax_refatm      ! Finite height of reference atmosphere (defined by p=0 and/or T=0, resp.)

    LOGICAL :: integral_averaging, calc_hhl

    INTEGER (KIND=iintegers) :: istata, corrcount

    !------------------------------------------------------------------------------
    !- End of header -
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE gen_ini_data
    !------------------------------------------------------------------------------

    yzroutine   = ' gen_ini_data'
    yzerrmsg(:) = ' '

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. '//TRIM(yzroutine)//'() ...'
    END IF

    kzdims(:) = 0_iintegers

    ALLOCATE( theta_ini( ie, je, ke) , STAT=istata)

    theta_ini = -999.99_wp

    ! update pointers to tracers (at nnew)
    CALL get_tracers()

    !-------------------------------------------------------
    ! 1. Initialization of the vertical coordinate parameter
    !-------------------------------------------------------
    !
    ! Data for some default vertical resolutions; other coordinate sets must
    ! be specified by the user

    zvco_p20 = (/ 0.02_wp , 0.060_wp, 0.110_wp, 0.160_wp, 0.220_wp, 0.280_wp, 0.340_wp, 0.400_wp, &
         0.465_wp, 0.535_wp, 0.600_wp, 0.660_wp, 0.720_wp, 0.780_wp, 0.840_wp, 0.890_wp, &
         0.925_wp, 0.950_wp, 0.975_wp, 0.992_wp, 1.000_wp /)

    zvco_p32 = (/ 0.020_wp, 0.046_wp, 0.077_wp, 0.113_wp, 0.153_wp, 0.193_wp, 0.233_wp, 0.273_wp, &
         0.313_wp, 0.353_wp, 0.393_wp, 0.433_wp, 0.473_wp, 0.513_wp, 0.553_wp, 0.593_wp, &
         0.633_wp, 0.673_wp, 0.712_wp, 0.749_wp, 0.783_wp, 0.813_wp, 0.839_wp, 0.862_wp, &
         0.883_wp, 0.903_wp, 0.922_wp, 0.940_wp, 0.956_wp, 0.970_wp, 0.982_wp, 0.992_wp, &
         1.000_wp /)

    zvco_p35 = (/ 0.0200_wp, 0.0400_wp, 0.0650_wp, 0.0930_wp, 0.1230_wp, 0.1540_wp, 0.1850_wp, &
         0.2160_wp, 0.2480_wp, 0.2810_wp, 0.3150_wp, 0.3510_wp, 0.3880_wp, 0.4250_wp, &
         0.4620_wp, 0.4990_wp, 0.5360_wp, 0.5730_wp, 0.6100_wp, 0.6470_wp, 0.6830_wp, &
         0.7180_wp, 0.7520_wp, 0.7840_wp, 0.8130_wp, 0.8390_wp, 0.8620_wp, 0.8830_wp, &
         0.9030_wp, 0.9220_wp, 0.9400_wp, 0.9560_wp, 0.9700_wp, 0.9820_wp, 0.9920_wp, &
         1.0000_wp  /)

    zvco_p40 = (/ 0.020_wp, 0.036_wp, 0.057_wp, 0.082_wp, 0.110_wp, 0.140_wp, 0.170_wp, 0.200_wp, &
         0.230_wp, 0.260_wp, 0.290_wp, 0.320_wp, 0.350_wp, 0.382_wp, 0.414_wp, 0.446_wp, &
         0.478_wp, 0.510_wp, 0.542_wp, 0.574_wp, 0.606_wp, 0.638_wp, 0.670_wp, 0.702_wp, &
         0.732_wp, 0.760_wp, 0.786_wp, 0.810_wp, 0.832_wp, 0.852_wp, 0.870_wp, 0.886_wp, &
         0.901_wp, 0.915_wp, 0.929_wp, 0.943_wp, 0.957_wp, 0.970_wp, 0.982_wp, 0.992_wp, &
         1.000_wp /)

    zvco_z20 = (/ 23580.0_wp, 18857.0_wp, 15616.0_wp, 13387.0_wp, 11357.0_wp, 9737.0_wp,   &
         8380.0_wp,  7209.0_wp,  6095.0_wp,  5033.0_wp,  4146.0_wp, 3397.0_wp,   &
         2703.0_wp,  2056.0_wp,  1451.0_wp,   974.0_wp,   653.0_wp,  430.0_wp,   &
         213.0_wp,    68.0_wp,     0.0_wp /)

    zvco_z32 = (/ 23580.0_wp, 20136.0_wp, 17579.0_wp, 15463.0_wp, 13663.0_wp, 12208.0_wp,   &
         10978.0_wp,  9911.0_wp,  8965.0_wp,  8113.0_wp,  7339.0_wp,  6627.0_wp,   &
         5968.0_wp,  5354.0_wp,  4779.0_wp,  4238.0_wp,  3727.0_wp,  3243.0_wp,   &
         2793.0_wp,  2386.0_wp,  2026.0_wp,  1719.0_wp,  1461.0_wp,  1239.0_wp,   &
         1040.0_wp,   854.0_wp,   681.0_wp,   520.0_wp,   378.0_wp,   256.0_wp,   &
         153.0_wp,    68.0_wp,     0.0_wp /)

    zvco_z35 = (/ 23580.4414_wp, 20773.3691_wp, 18455.5488_wp, 16559.7266_wp, 14970.4609_wp, &
         13623.1367_wp, 12477.7646_wp, 11478.0508_wp, 10561.6865_wp,  9712.8486_wp, &
         8919.9834_wp,  8154.0073_wp,  7431.7842_wp,  6764.6855_wp,  6144.3555_wp, &
         5564.2437_wp,  5019.1157_wp,  4504.7178_wp,  4017.5457_wp,  3554.6763_wp, &
         3125.2937_wp,  2725.8110_wp,  2353.2324_wp,  2015.3969_wp,  1719.1958_wp, &
         1461.1770_wp,  1238.5348_wp,  1039.6235_wp,   853.8964_wp,   680.6733_wp, &
         519.3526_wp,   378.1536_wp,   256.2476_wp,   152.9479_wp,    67.6839_wp, &
         0.0000_wp /)

    zvco_z40 = (/ 23580.0_wp, 21238.0_wp, 19111.0_wp, 17244.0_wp, 15617.0_wp, 14202.0_wp,   &
         13011.0_wp, 11978.0_wp, 11064.0_wp, 10243.0_wp,  9496.0_wp,  8809.0_wp,   &
         8174.0_wp,  7545.0_wp,  6958.0_wp,  6407.0_wp,  5889.0_wp,  5399.0_wp,   &
         4934.0_wp,  4491.0_wp,  4069.0_wp,  3665.0_wp,  3278.0_wp,  2906.0_wp,   &
         2571.0_wp,  2268.0_wp,  1995.0_wp,  1749.0_wp,  1530.0_wp,  1335.0_wp,   &
         1162.0_wp,  1012.0_wp,   872.0_wp,   744.0_wp,   618.0_wp,   493.0_wp,   &
         369.0_wp,   256.0_wp,   153.0_wp,    68.0_wp,     0.0_wp /)
!CPS
!    zvco_z50 =(/22000.0_wp, 21000.0_wp, 20028.0_wp, 19085.0_wp, 18170.0_wp, 17282.0_wp,      &
!                16421.0_wp, 15587.0_wp, 14780.0_wp, 13998.0_wp, 13242.0_wp, 12512.0_wp,      &
!                11807.0_wp, 11126.0_wp, 10470.0_wp,  9837.0_wp,  9228.0_wp,  8642.0_wp,      &
!                 8080.0_wp,  7539.0_wp,  7021.0_wp,  6525.0_wp,  6050.0_wp,  5596.0_wp,      &
!                 5162.0_wp,  4750.0_wp,  4357.0_wp,  3983.0_wp,  3630.0_wp,  3295.0_wp,      & 
!                 2978.0_wp,  2680.0_wp,  2400.0_wp,  2137.0_wp,  1891.0_wp,  1662.0_wp,      &
!                 1450.0_wp,  1253.0_wp,  1072.0_wp,   907.0_wp,   757.0_wp,   621.0_wp,      &
!                  500.0_wp,   392.0_wp,   298.0_wp,   217.0_wp,   150.0_wp,    94.0_wp,      &
!                   51.0_wp,    20.0_wp,     0.0_wp /)

   zvco_z50 =(/22000.0_wp, 21000.0_wp, 20028.570312_wp, 19085.359375_wp, 18170.0_wp, 17282.140625_wp,         &  
          16421.429688_wp, 15587.500000_wp, 14780.0_wp, 13998.570312_wp, 13242.859375_wp, 12512.5_wp,         &
          11807.136719_wp, 11126.429688_wp, 10470.0_wp, 9837.5_wp, 9228.570312_wp, 8642.859375_wp, 8080.0_wp, &
          7539.636719_wp , 7021.429688_wp , 6525.0_wp , 6050.0_wp, 5596.066406_wp, 5162.859375_wp, 4750.0_wp, &
          4357.136719_wp , 3983.929932_wp , 3630.0_wp , 3295.0_wp, 2978.570068_wp, 2680.360107_wp, 2400.0_wp, &
          2137.139893_wp , 1891.429932_wp , 1662.5_wp , 1450.0_wp, 1253.569824_wp, 1072.859863_wp, 907.5_wp,  &
          757.139893_wp  , 621.429932_wp  ,  500.0_wp , 392.5_wp, 298.569824_wp, 217.860001_wp, 150.0_wp,     &
          94.639999_wp   , 51.429993_wp   ,   20.0_wp , 0.0_wp /)
!CPS

    !-------------------------------------------------------
    !
    ! Store the parameters for the vertical grid and the reference atmosphere into
    ! the new structures vcoord, refatm

    vcoord%ivctype     = ivctype
    vcoord%ivcoord_id  = 0
    vcoord%nlevels     = ke+1
    vcoord%vcflat      = vcflat
    vcoord%kflat       = -1   ! will be set later by CALL reference_atmosphere_xxx()
    vcoord%vc_uuid(:)  = 'x'  ! ???
    vcoord%vert_coord(:) = -999.99_wp
    vcoord%sigm_coord(:) = -999.99_wp

    refatm%irefatm           =  irefatm
    refatm%irefatm_id        =  99
    refatm%p0sl              =  p0sl
    refatm%t0sl              =  t0sl
    refatm%dt0lp             =  dt0lp
    refatm%delta_t           =  delta_t
    refatm%h_scal            =  h_scal
    refatm%bvref             =  bvref

    !-------------------------------------------------------
    !
    ! Specification of the vertical model levels according
    ! to the namelist parameters read by SR input_artifctl():

    izerror  = 0
    yzerrmsg(:) = ' '

    IF ( vcoord%ivctype == 1 ) THEN ! Pressure based hybrid coordinate
      ! vcoord(k) cooresponds to pressure above sea-level
      ! normalized by POSL
      SELECT CASE (TRIM(ADJUSTL(zspacing_type)))
      CASE('predefined', 'PREDEFINED', 'Predefined')
        IF (ke == 20) THEN 
          vcoord%vcflat = 0.020_wp
          DO k = 1, ke1
            vcoord%sigm_coord(k) = zvco_p20(k)
          ENDDO
        ELSEIF (ke == 32) THEN 
          vcoord%vcflat = 0.220_wp
          DO k = 1, ke1
            vcoord%sigm_coord(k) = zvco_p32(k)
          ENDDO
        ELSEIF (ke == 35) THEN 
          vcoord%vcflat = 0.220_wp
          DO k = 1, ke1
            vcoord%sigm_coord(k) = zvco_p35(k)
          ENDDO
        ELSEIF (ke == 40) THEN 
          vcoord%vcflat = 0.220_wp
          DO k = 1, ke1
            vcoord%sigm_coord(k) = zvco_p40(k)
          ENDDO
        ELSE
          WRITE (*,*) '***   No vertical coordinates are pre-specified    ***'
          WRITE (*,*) '***  for ke=',ke,' layers and ivctype=',vcoord%ivctype,'  ***'
          WRITE (*,*) '*** Please edid subroutine '//TRIM(yzroutine)//' to proceed ***'
          STOP
        ENDIF
      CASE ('vcoordvec', 'VCOORDVEC', 'Vcoordvec')
        ! if vcoord is specified in namelist, copy it
        vcoord%sigm_coord(1:ke+1) = vcoordvec(1:ke+1)
      CASE default
        WRITE (*,*) '***   No vertical coordinates of type',TRIM(ADJUSTL(zspacing_type)),'    ***'
        WRITE (*,*) '***   defined for ivctype=',vcoord%ivctype,'  ***'
        WRITE (*,*) '***   Please edid subroutine '//TRIM(yzroutine)//' to proceed ***'
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      END SELECT

    ELSE IF ( vcoord%ivctype == 2 ) THEN ! Height-based hybrid Gal-Chen coordinate
      ! vcoord(k) cooresponds to nonnormalized height
      ! above sea-level (z=0)
      SELECT CASE (TRIM(ADJUSTL(zspacing_type)))
      CASE('predefined', 'PREDEFINED', 'Predefined')
        IF (ke == 20) THEN
          vcoord%vcflat = 11357.0_wp
          DO k = 1, ke1
            vcoord%vert_coord(k) = zvco_z20(k)
          ENDDO
        ELSEIF (ke == 32) THEN
          vcoord%vcflat = 11357.0_wp
          DO k = 1, ke1
            vcoord%vert_coord(k) = zvco_z32(k)
          ENDDO
        ELSEIF (ke == 35) THEN
          vcoord%vcflat = 11357.0_wp
          DO k = 1, ke1
            vcoord%vert_coord(k) = zvco_z35(k)
          ENDDO
        ELSEIF (ke == 40) THEN
          vcoord%vcflat = 11357.0_wp
          DO k = 1, ke1
            vcoord%vert_coord(k) = zvco_z40(k)
          ENDDO
!CPS
        ELSEIF (ke == 50) THEN
          vcoord%vcflat = 11000.0_wp
          DO k = 1, ke1
            vcoord%vert_coord(k) = zvco_z50(k)
          ENDDO
        ELSE
          WRITE (*,*) '***   No vertical coordinates are pre-specified    ***'
          WRITE (*,*) '***  for ke=',ke,' layers and ivctype=',vcoord%ivctype,'  ***'
          WRITE (*,*) '*** Please edit subroutine '//TRIM(yzroutine)//' to proceed ***'
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        END IF
      CASE ('linear', 'LINEAR', 'Linear')
        ! compute equidistant levels:
        DO k = 1, ke1
          vcoord%vert_coord(k) = zz_top * (1.0_wp - REAL(k-1,wp)/REAL(ke,wp))
        END DO
      CASE ('galchen', 'GALCHEN', 'Galchen')
        ! compute height of levels according to the formula
        ! z = zz_top * ( 2/pi*arccos( (k-1)/ke1 ) )**exp_galchen
        DO k = 1, ke1
          vcoord%vert_coord(k) = zz_top * &
               ( 2.0_wp/pi*ACOS( REAL(k-1,wp)/REAL(ke,wp) ) )**exp_galchen
        ENDDO
      CASE ('vcoordvec', 'VCOORDVEC', 'Vcoordvec')
        ! if vcoord is specified in namelist, copy it
        vcoord%vert_coord(1:ke+1) = vcoordvec(1:ke+1)
      CASE default
        WRITE (*,*) '***   No vertical coordinates of type',TRIM(ADJUSTL(zspacing_type)),'    ***'
        WRITE (*,*) '***   defined for ivctype=',vcoord%ivctype,'  ***'
        WRITE (*,*) '***   Please edid subroutine '//TRIM(yzroutine)//' to proceed ***'
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      END SELECT

    ELSE IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN ! SLEVE coordinates

      SELECT CASE (TRIM(ADJUSTL(zspacing_type)))
      CASE('predefined', 'PREDEFINED', 'Predefined')
        IF (ke == 20) THEN
          DO k = 1, ke1
            vcoord%vert_coord(k) = zvco_z20(k)
          ENDDO
          vcoord%vcflat = vcoord%vert_coord(1)
        ELSEIF (ke == 32) THEN
          DO k = 1, ke1
            vcoord%vert_coord(k) = zvco_z32(k)
          ENDDO
          vcoord%vcflat = vcoord%vert_coord(1)
        ELSEIF (ke == 35) THEN
          DO k = 1, ke1
            vcoord%vert_coord(k) = zvco_z35(k)
          ENDDO
          vcoord%vcflat = vcoord%vert_coord(1)
        ELSEIF (ke == 40) THEN
          DO k = 1, ke1
            vcoord%vert_coord(k) = zvco_z40(k)
          ENDDO
          vcoord%vcflat = vcoord%vert_coord(1)
        ELSE
          WRITE (*,*) '***   No vertical coordinates are pre-specified    ***'
          WRITE (*,*) '***  for ke=',ke,' layers and ivctype=',vcoord%ivctype,'  ***'
          WRITE (*,*) '*** Please edit subroutine '//TRIM(yzroutine)//' to proceed ***'
          CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
        END IF
      CASE ('linear', 'LINEAR', 'Linear')
        ! compute equidistant levels:
        DO k = 1, ke1
          vcoord%vert_coord(k) = zz_top * (1.0_wp - REAL(k-1,wp)/REAL(ke,wp))
        END DO
        vcoord%vcflat = vcoord%vert_coord(1)
      CASE ('galchen', 'GALCHEN', 'Galchen')
        ! compute height of levels according to the formula
        ! z = zz_top * ( 2/pi*arccos( (k-1)/ke1 ) )**exp_galchen
        DO k = 1, ke1
          vcoord%vert_coord(k) = zz_top * &
               ( 2.0_wp/pi*ACOS( REAL(k-1,wp)/REAL(ke,wp) ) )**exp_galchen
        ENDDO
        vcoord%vcflat = vcoord%vert_coord(1)
      CASE ('vcoordvec', 'VCOORDVEC', 'Vcoordvec')
        ! if vcoord is specified in namelist, copy it
        vcoord%vert_coord(1:ke+1) = vcoordvec(1:ke+1)
      CASE default
        WRITE (*,*) '***   No vertical coordinates of type', &
             TRIM(ADJUSTL(zspacing_type)),'    ***'
        WRITE (*,*) '***   defined for ivctype=',vcoord%ivctype,'  ***'
        WRITE (*,*) '***   Please edit subroutine '//TRIM(yzroutine)//' to proceed ***'
        CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
      END SELECT

    ENDIF

    ! .. Construct new (Version 5.1) vector hhl_prof, the reference height profile
    !     over the lowest orography height in the model domain which is equal to or above sea level,
    !     to avoid depressions. (Replaces vcoord%vert_coord in some of the other modules because of grib2)
    hhl_prof(0) = vcoord%vcflat
    hhl_prof(1:ke1) = vcoord%vert_coord(1:ke1)

    !======================================================================
    !  Catch further errors associated with coordinate specifications:
    !======================================================================

    !.. Check the values of vcoord against other namelist parameters:

    CALL check_vcoord()

    !----------------------------------------------------------------------
    ! 2. Computation of the artificial orography (hsurf) in form of a 
    !    bell-shaped  mountain  with the height zhmax (in meter) and half 
    !    extension zda (multiplication of grid length in x-direction)
    !----------------------------------------------------------------------

    !.. If namelist parameter linit_realoro=.true.:
    !
    !   Possibility to read an external orography field; needs the 
    !   function read_ascii_field_2d()
    !
    ! FORMAT OF THE FILE:
    !         - 1 header line (arbitrary content)
    !         - 1 line with field dimensions i,j
    !         - 1 long column with the data (first index varies first).
    !
    ! FOR EXAMPLE:
    !
    !>> BEGIN ASCII-FILE:
    ! # orography height [m]
    ! 461 421
    ! 0.001
    ! 0.002
    ! 0.001
    ! ...
    !<< END ASCII-FILE
    
    IF (linit_realoro) THEN

      CALL read_ascii_field_2d(TRIM(orofile), i_shift_realoro, j_shift_realoro, &
           href_oro, hsurf, ierrstat)
      IF (my_cart_id == 0) THEN
        IF (ierrstat /= 0 ) THEN
          WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
          CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
               'src_artifdata, reading hsurf from file '//TRIM(orofile)// &
               ' with read_ascii_field_2d()')
          RETURN
        END IF
      END IF

    ELSE

      ! Else set constant base height everywhere:
      hsurf(:,:) = href_oro

    END IF

    ! idealized hills/valleys might be specified or added to real orography:
    SELECT CASE (itype_topo)

    CASE (0)

      ! No (additional) hills.

    CASE (1)

      !========================================================================
      ! Add analytic hills/valleys to existing hsurf.
      !
      ! This is the mechanism which recognizes various namelist parameters
      ! associated with orography ('hill...').
      !
      ! What is actually done depends on the sign of the 'hillheight' and
      ! on the switch 'hill_combineaction' (1 = add/subtract, 
      !   2 = take max/min, dep. on hillsign):
      !
      ! You may extend it by adding your own 'hilltype' case to 
      !   the source code of add_orography_ana().
      !
      !========================================================================

      IF (ANY(lhill)) THEN

        CALL add_orography_ana(ierrstat)
        IF (ierrstat /= 0) THEN
          CALL model_abort (my_world_id, 10006, &
               'Error in add_orography_ana(), src_artifdata', &
               'add_orography_ana(), src_artifdata')
          RETURN
        END IF

      END IF

    CASE (2)

!!$========================================================================
!!$
!!$ Here is the possibility to implement own specialized orographies,
!!$ independent of the mechanism provided in add_orography_ana().
!!$
!!$ ADD YOUR OWN SUBROUTINE HERE!!!  
!!$
!!$     You have to set the field hsurf(1:ie,1:je) and you
!!$     have to do the same boundary exchange as in add_orography_ana(), see
!!$     the body of this subroutine below! 
!!$
!!$========================================================================

      WRITE (*,*) ' itype_topo = 2 not yet implemented ! '
      STOP        

    CASE default

      WRITE (*,'(a,i3)') 'ERROR: '//TRIM(yzroutine)//'(), src_artifdata, '// &
           'wrong itype_topo = ', itype_topo, &
           ' specified in namelist IDEAL'
      CALL model_abort (my_world_id, 10006, &
           'ERROR in '//TRIM(yzroutine)//'(), src_artifdata', &
           TRIM(yzroutine)//'(), src_artifdata')
      RETURN

    END SELECT


    !------------------------------------------------------------------------------
    ! 3. Initialization of constant fields (ozone content, soil description,...)
    !------------------------------------------------------------------------------


    SELECT CASE ( itype_soil_c )

    CASE (1)

      ! Initialize soil with constant values
      gz0    (:,:) = g * z0_c   ! z0 * g in m^2/s^2
      fr_land(:,:) = fr_land_c
      soiltyp(:,:) = soiltyp_c
      plcov  (:,:) = plcov_c
      lai    (:,:) = lai_c
      rootdp (:,:) = rootdp_c

      ! In case someone uses the radiation scheme:
      IF (lrad) THEN
        vio3    (:,:) = 0.06_wp   !vio3_c ! vertical integrated ozone content [Pa O3]
        hmo3    (:,:) = 4200.0_wp !hmo3_c ! ozone maximum                     [Pa]
      END IF

      ! For taking into account forests:
      IF (lforest) THEN
        for_e(:,:) = for_e_c
        for_d(:,:) = for_d_c
      ENDIF

      ! For sea ice modeling and lake modeling:
      IF (lseaice .OR. llake) THEN
        h_ice(:,:,nnew) = h_ice_c
      ENDIF

      ! For the subgrid scale orography scheme:
      IF (lsso) THEN
        sso_stdh  = 20.0_wp      ! standard deviation of sub-grid scale orography ( m   )
        sso_gamma = 0.0_wp       ! anisotropy of sub-grid scale orography          --
        sso_theta = 0.25_wp * pi ! angle betw. principal axis of orography and E  ( rad )
        sso_sigma = 0.05_wp      ! mean slope of sub-grid scale orography          --
      END IF

    CASE (2)

      !.. Read soil parameters from external ASCII files in the 
      !   same way as with the orography:
      !   (uses internal SR read_soilparams_from_files() )

      ! First: dummy initialization of some fields which,
      ! if not defined with a valid value, may lead to model crashes in src_gridpoints.f90:
      soiltyp(:,:) = soiltyp_c
      rootdp (:,:) = rootdp_c

      ! Then read the actual files:
      CALL read_soilparams_from_files()

      ! There is no GRIB or NETCDF input implemented. 
      ! The ASCII-files should be generated by the users themselves. This is
      ! achievable with some phantasy, the DWD web interface, int2lm, wgrib, matlab, idl, ...
      !
      ! Files for the following variables are necessary (probably more in the future):
      !
      ! IN ANY CASE:
      !   z0      : roughness length [m] (Give name in NL-parameter "z0file")
      !   fr_land : land fraction    [-] (                          "frlandfile")
      !   plcov   : plant cover      [-] (                          "plcovfile")
      !   lai     : leaf area index  [-] (                          "laifile")
      !   soiltyp : type of soil     [-] (                          "soiltypefile")
      ! if lsoil=.true.:
      !   rootdp  : root depth       [m] (                          "rootdpfile")
      ! if lforest=.true.:
      !   for_e   : aerea fraction of evergreen forest   [-] (      "forefile")
      !   for_d   : aerea fraction of deciduous forest   [-] (      "fordfile")
      ! if lseaice=.true. .or. llake=.true.:
      !   h_ice   : thickness of ice cover               [m] (      "hicefile")
      ! if llake=.true.:
      !   ... MISSING ...
      ! if lsso=.true.:
      !   sso_stdh : std. dev. of sub-grid scale orography [m] (    "ssostdhfile"  )
      !   sso_gamma: anisotr. of sub-grid scale orography [-] (     "ssogammafile" )
      !   sso_theta: angle betw. princ. axis of orogr. and E [rad] ("ssothetafile" )
      !   sso_sigma: mean slope of sub-grid scale orography [-] (   "ssosigmafile" )
      !
      ! For example, ASCII files may be extracted from the extpar-grib-files
      ! by wgrib or the script "grib_decode" from Ulrich Blahak.
      !
      !
      ! FORMAT OF THE FILES:
      !         - 1 header line (arbitrary content)
      !         - 1 line with field dimensions i,j
      !         - 1 long column with the data (first index varies first).
      !
      ! FOR EXAMPLE:
      !
      !>> BEGIN ASCII-FILE:
      ! # roughness length [m]
      ! 461 421
      ! 0.001
      ! 0.002
      ! 0.001
      ! ...
      !<< END ASCII-FILE

    CASE default

      CALL model_abort (my_world_id, 10006, &
           'Error in itype_soil_c, src_artifdata', &
           'src_artifdata, wrong value of itype_soil_c encounterd in '//TRIM(yzroutine)//'()')
      RETURN

    END SELECT

    ! Set depth_lk to -1.0, because the initialization of llake=.true. is currently not implemented in src_artifdata.f90:
    depth_lk = -1.0_wp

    ! Compute the landmask and the number of landpoints
    !   land point: fr_land >= 0.5, llandmask=.true.
    !   sea point:  fr_land <  0.5, depth_lk <  0.0, llandmask=.false.
    !   lake point: fr_land <  0.5, depth_lk >= 0.0, llandmask=.false.
    nlandpoints = 0
    DO j = 1,je
      DO i = 1,ie
        IF (fr_land(i,j) >= 0.5_wp) THEN
          llandmask(i,j) = .TRUE.
          nlandpoints    = nlandpoints + 1
        ELSE
          llandmask(i,j) = .FALSE.
        ENDIF
      ENDDO
    ENDDO


    IF (num_compute > 1) THEN
      nlandpoints_tot = COUNT( llandmask( istart:iend, jstart:jend ) )
      CALL global_values (nlandpoints_tot, 1, 'SUM', imp_integers, icomm_cart, -1, &
           yzerrmsg, izerror)
    ELSE
      nlandpoints_tot = COUNT( llandmask )
    ENDIF

    !==============================================================================
    !.. Consistency checks on surface parameters
    !     Some more should be included here in the future!
    !==============================================================================

    IF ( ANY(gz0 < 1.0E-20_wp ) ) THEN
      CALL model_abort (my_world_id, 13475, &
           'ERROR in constant parameter spec., src_artifdata', &
           'src_artifdata, unphysical values of gz0 < 1.0E-20 in '//TRIM(yzroutine)//'()')
    END IF

    IF ( ANY(NINT(soiltyp) < 1) .OR. ANY(NINT(soiltyp) > 10) ) THEN
      CALL model_abort (my_world_id, 13476, &
           'ERROR in constant parameter spec., src_artifdata', &
           'src_artifdata, wrong values for soiltyp in '//TRIM(yzroutine)//'(), should be >= 1.0 and <= 10.0')
    END IF

    IF ( ANY(fr_land < 0.0_wp) .OR. ANY(fr_land > 1.0_wp) ) THEN
      CALL model_abort (my_world_id, 13477, &
           'ERROR in constant parameter spec., src_artifdata', &
           'src_artifdata, wrong values for fr_land in '//TRIM(yzroutine)//'(), should be >= 0.0 and <= 1.0')
    END IF

    IF ( ANY(plcov < 0.0_wp) .OR. ANY(plcov > 1.0_wp) ) THEN
      CALL model_abort (my_world_id, 13478, &
           'ERROR in constant parameter spec., src_artifdata', &
           'src_artifdata, wrong values for plcov in '//TRIM(yzroutine)//'(), should be >= 0.0 and <= 1.0')
    END IF

    IF ( ANY(lai < 0.0_wp) ) THEN
      CALL model_abort (my_world_id, 13479, &
           'ERROR in constant parameter spec., src_artifdata', &
           'src_artifdata, wrong values for lai in '//TRIM(yzroutine)//'(), should be >= 0.0')
    END IF

    IF ( ANY(rootdp < 0.0_wp) ) THEN
      CALL model_abort (my_world_id, 13480, &
           'ERROR in constant parameter spec., src_artifdata', &
           'src_artifdata, wrong values for rootdp in '//TRIM(yzroutine)//'(), should be >= 0.0')
    END IF

    IF ( ANY( llandmask .AND. NINT(soiltyp,iintegers) > 8 ) ) THEN
      CALL model_abort (my_world_id, 13480, &
           'ERROR in constant parameter spec., src_artifdata', &
           'src_artifdata, wrong values for soiltyp at land points in '//TRIM(yzroutine)//'(), should be <= 8.0.'//&
           ' Hint: 9.0 (ocean) and 10.0 (sea ice) are only possible for fr_land < 0.5!')
    END IF

    IF (lseaice .OR. llake) THEN
      corrcount = 0
      DO j=1,je
        DO i=1,ie
          IF (.NOT.llake) THEN
            IF (llandmask(i,j) .AND. h_ice(i,j,nnew) > 0.0_wp) THEN
              corrcount = corrcount + 1
              h_ice(i,j,nnew) = 0.0_wp
            END IF
          ELSE
            IF (llandmask(i,j) .AND. &
                 depth_lk(i,j) < 0.0_wp .AND. h_ice(i,j,nnew) > 0.0_wp) THEN
              corrcount = corrcount + 1
              h_ice(i,j,nnew) = 0.0_wp
            END IF
          END IF
        END DO
      END DO
      IF (corrcount > 0) THEN
        WRITE (*,'(a,i5,a)') 'WARNING: '//TRIM(yzroutine)//': H_ICE corrected to 0.0 at ', &
             corrcount,' land points!'
      END IF
      ! Set very small values of h_ice to 0.0:
      WHERE (h_ice(:,:,nnew) < 1e-12) h_ice(:,:,nnew) = 0.0_wp
    END IF

    IF (lseaice .OR. llake) THEN
      corrcount = 0
      DO j=1,je
        DO i=1,ie
          IF (.NOT.llandmask(i,j)) THEN
            IF (h_ice(i,j,nnew) >= 1e-12_wp .AND. NINT(soiltyp(i,j),iintegers) /= 10) THEN
              corrcount = corrcount + 1
              soiltyp(i,j) = 10.0_wp    ! force soiltype to ice
            ELSE IF (h_ice(i,j,nnew) < 1e-12_wp  .AND. NINT(soiltyp(i,j),iintegers) /= 9) THEN
              corrcount = corrcount + 1
              soiltyp(i,j) = 9.0_wp     ! force soiltype to water
            END IF
          END IF
        END DO
      END DO
      IF (corrcount > 0) THEN
        WRITE (*,'(a,i5,a)') 'WARNING: '//TRIM(yzroutine)//': SOILTYP corrected to sea ice or water at ', &
             corrcount,' sea/lake points!'
      END IF
    ELSE
      corrcount = 0
      DO j=1,je
        DO i=1,ie
          IF (.NOT.llandmask(i,j)) THEN
            corrcount = corrcount + 1
            soiltyp(i,j) = 9.0_wp     ! force soiltype to water
          END IF
        END DO
      END DO
      IF (corrcount > 0) THEN
        WRITE (*,'(a,i5,a)') 'WARNING: '//TRIM(yzroutine)//': SOILTYP corrected to water at ', &
             corrcount,' sea points!'
      END IF
    END IF
    

    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    ! Only for SLEVE coordinate (ivctype = 3/4)
    ! 3.5 Splitting of the topography into a large-scale and a small-scale part
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------

    ! allocate initialize splitted topo parts
    ALLOCATE( hsurfs(ie,je,2) )
    hsurfs(:,:,:) = 0.0_wp

    IF ( ANY( vcoord%ivctype == (/3,4/) ) ) THEN

      ! allocate memory
      ALLOCATE( hsurfs_tot(ie_tot,je_tot,2), &
           hsurf_tot(ie_tot,je_tot) )

      ! initialize splitted topo parts
      hsurfs_tot(:,:,:) = 0.0_wp
      hsurf_tot(:,:)    = 0.0_wp

      ! collect full topo from all PE's
      IF (num_compute == 1 ) THEN     ! we are running on one PE
        hsurf_tot = hsurf
      ELSE                            ! we are running on multiple PE's
        CALL gather_field(hsurf,ie,je,hsurf_tot,ie_tot,je_tot,0,izerror)
      ENDIF

      ! split topo on PE 0
      IF (my_cart_id == 0) THEN
        CALL sleve_split_oro(hsurf_tot, hsurfs_tot, ie_tot, je_tot, nfltvc,  &
             0_iintegers, svc1, svc2, vcoord%vcflat, nuspecif,               &
             my_cart_id, izerror, yzerrmsg)
      ENDIF

      ! distribute splitted topo to all PE's
      IF (num_compute == 1) THEN      ! we are running on one PE
        hsurfs = hsurfs_tot
      ELSE                            ! we are running on multiple PE's
        CALL distribute_field (hsurfs_tot(:,:,1), ie_tot, je_tot,           &
             hsurfs    (:,:,1), ie,     je,     0, izerror)
        CALL distribute_field (hsurfs_tot(:,:,2), ie_tot, je_tot,           &
             hsurfs    (:,:,2), ie,     je,     0, izerror)
      ENDIF

      ! dellaocte memory
      DEALLOCATE( hsurfs_tot, hsurf_tot )

    ENDIF

    !------------------------------------------------------------------------------
    ! 4. Computation of the constant reference atmosphere
    !------------------------------------------------------------------------------

    ! new switch for reference_atmosphere_xxx() to actually calculate hhl
    ! and not assume that it already exists
    ! (e.g., it would exist, if hhl had been read from a grib2-file in full precision):
    calc_hhl = .TRUE.

    SELECT CASE ( refatm%irefatm )
    CASE (1)

      ! preliminary choice to be compatible with the current fast_waves solver
      IF ( itype_fast_waves == 2 ) THEN
        lanalyt_calc_t0p0 = .TRUE.   ! necessary only if irefatm=1
      ELSE
        lanalyt_calc_t0p0 = .FALSE.  ! necessary only if irefatm=1
      END IF

      CALL reference_atmosphere                                                   &
           ( hhl, p0, p0hl, rho0, t0, t0hl, dp0, hsurf, hsurfs, ie, je, ke,       &
             refatm, vcoord, svc1, svc2, r_d, g, lanalyt_calc_t0p0,   &
             calc_hhl, yzerrmsg, izerror)

    CASE (2)
      CALL reference_atmosphere_2                                                 &
           ( hhl, p0, p0hl, rho0, t0, t0hl, dp0, hsurf, hsurfs, ie, je, ke,       &
             refatm, vcoord, svc1, svc2, r_d, g, calc_hhl, yzerrmsg, izerror)

    CASE (3)
      CALL reference_atmosphere_BVconst                                           &
           ( hhl, p0, p0hl, rho0, t0, t0hl, dp0, hsurf, hsurfs, ie, je, ke,       &
             refatm, vcoord, svc1, svc2, r_d, g, calc_hhl, yzerrmsg, izerror)

    CASE default
      CALL model_abort (my_world_id, 10007, &
           'ERROR in reference atmosphere spec., src_artifdata', &
           'src_artifdata, wrong value of irefatm encounterd in '//TRIM(yzroutine)//'()')
    END SELECT

    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 10000+izerror, yzerrmsg, 'reference_atmosphere_x')
    ENDIF

    CALL k_index_of_pressure_levels(  refatm%p0sl, vcoord%sigm_coord, ke, llm,      &
             klv950, klv850, klv800, klv700, klv500, klv400, klv300 )

    ! free memory no longer needed
    DEALLOCATE( hsurfs )

    CALL init_grid_metrics (istata)
    IF (istata /= 0) THEN
      IF (my_cart_id == 0) THEN
        PRINT *, ' Initialization of Grid Metrics failed'
        CALL model_abort (my_world_id, 10007, ' Initialization of Grid Metrics failed', &
                          TRIM(yzroutine)//' in src_artifdata')
      ENDIF
    ENDIF


    ! Set model top height:
    h_top = vcoord%vert_coord(1)

    ! compute height of full levels
    zml(:,:,:) = 0.5_wp*( hhl(:,:,1:ke) + hhl(:,:,2:ke+1) )
    DO i = 1, ie 
      IF (i < ie)  THEN
        zml_u(i,:,:) = 0.5_wp * (zml(i,:,:)+zml(i+1,:,:))
      ELSE
        zml_u(i,:,:) = zml(i,:,:)
      END IF
    ENDDO
    DO i = j, je 
      IF (j < je)  THEN
        zml_v(:,j,:) = 0.5_wp * (zml(:,j,:)+zml(:,j+1,:))
      ELSE
        zml_v(:,j,:) = zml(:,j,:)
      END IF
    ENDDO

    !------------------------------------------------------------------------------
    ! 4. Initialization of the wind fields and thermodynamic variables.
    !
    !    Initial values for t, pp, qv and qc are determined from analytically 
    !    given profiles of the temperature and the relative humidity or by
    !    reading in an external sounding file. The terrain following
    !    coordinates are correctly taken into account.
    !
    !------------------------------------------------------------------------------


    ! itype_artifprofiles:
    !   1 = file only
    !   2 = analytic only
    !   3 = combination of tqv from analytic and uv from file
    !   4 = combination of tqv from file and uv from analytic

    ! initialize the limiting max qv (is set to a lower value for the WK82 test case)
    zmaxqv = HUGE(1.0_wp)
    ! initialize the switch to determine whether the T or Theta is used
    ! for defining the temperature profile:
    ltheta_ini = .FALSE.    


    !*********************************************************************************
    ! Set T, Relhum, T_s, ps, relhum_s and/or U, V from ASCII file:
    !*********************************************************************************

    SELECT CASE (itype_artifprofiles)

    CASE (1, 3, 4)

      !.. if integral_averaging = .true., then values are calculated as 
      !   integral averages over the model layers,
      !   treating the external profile data as linear splines;
      !   otherwise, they are simply linearily interpolated 
      !   from the EXTERNAL profile DATA:

      integral_averaging = .FALSE. !CPS .TRUE. 

      !.. If  lps_from_file = .TRUE. :
      !     Interpolate surface pressure from pressure given in 
      !     radiosonde file.
      !   Else: 
      !     Integrate pressure analytically from radiosonde data
      !     and then interpolate surface pressure.
      !     

      !.. If rasofile_t_is_theta = .TRUE. :
      !     Rasofile contains pot. temp. instead of temp.
      !     and t(:,:,:,nnew) will contain pot. temp.
      !     Therefore, set ltheta_ini = .TRUE.
      
      ltheta_ini = rasofile_t_is_theta

      IF ( rasofile_t_is_theta ) THEN
        
        CALL read_raso(rasofile, zml, hhl, theta_ini(:,:,:), zrhf, &
             u(:,:,:,nnew), v(:,:,:,nnew), &
             ps(:,:,nnew), t_s(:,:,nnew), zrhs, &
             integral_averaging, lps_from_file, rasofile_t_is_theta)

      ELSE

        CALL read_raso(rasofile, zml, hhl, t(:,:,:,nnew), zrhf, &
             u(:,:,:,nnew), v(:,:,:,nnew), &
             ps(:,:,nnew), t_s(:,:,nnew), zrhs, &
             integral_averaging, lps_from_file, rasofile_t_is_theta)

      END IF

    CASE (2)

      ! pure analytic profiles from below are taken!

    CASE default

      CALL model_abort (my_world_id, 10008, 'ERROR in itype_artifprofiles, src_artifdata', &
           'src_artifdata, wrong value of itype_artifprofiles encounterd in '//TRIM(yzroutine)//'()')
      RETURN

    END SELECT


    !*********************************************************************************
    ! Set T, Relhum, T_s, ps and relhum_s analytically (if desired):
    !*********************************************************************************

    IF (itype_artifprofiles == 2 .OR. itype_artifprofiles == 3) THEN

      SELECT CASE (itype_anaprof_tqv)

        !*****************************************************************************
        !.. If you would like to add your own intitialization routine for T and Qv,
        !   the following input/output parameters are necessary:
        !*****************************************************************************
        !
        !   INPUT:
        !     zml(:,:,:)     :  Height of main levels [m]
        !     hsurf(:,:)     :  Height of orography [m]
        !   
        !   OUTPUT:
        !     t(:,:,:,nnew) or theta_ini(:,:,:) : temperature or potential temperature
        !             (if potential temperature is output, then set ltheta_ini = .TRUE.!)
        !     zrhf(:,:,:)    :  Rel. humidity [-]
        !     ps(:,:,nnew)   :  Surface pressure [Pa]
        !     t_s(:,:,nnew)  :  Surface temperature [K]
        !     zrhs(:,:)      :  Rel. humidity at the surface [-] 
        !                       
        !  --> IF E(T) > p, zrh and zrhs will be limited to zrhs_max = p / E(T) later
        !      when calculating the pressure by the iterative method below !!!
        !
        !*****************************************************************************

      CASE (1)

        ! Weisman-Klemp:
        CALL tqv_wk82(zml, hsurf, theta_ini(:,:,:), zrhf, &
             ps(:,:,nnew), t_s(:,:,nnew), zrhs)
        zmaxqv = qv_max_wk
        ltheta_ini = .TRUE.

      CASE (2)

        ! Layered polytrope atmosphere.
        ! Number of layers determined by NAMELIST PARAMETER nlayers_poly:

        CALL tqv_ana_polylayers(zml, hsurf, t(:,:,:,nnew), &
             zrhf, ps(:,:,nnew), t_s(:,:,nnew), zrhs)

        !******************************************************************************
        ! NOTE: THE GRADIENTS TGR_POLY RESP. RHGR_POLY ARE POSITIVE FOR *DECREASING* 
        !       TEMPERATURE RESP. RELHUM WITH HEIGHT!
        !******************************************************************************

      CASE (3)

        ! Layered atmosphere with piecewise constant Brunt-Vaissala-Frequency N.
        ! Number of layers determined by namelist parameter nlayers_poly.
        ! Each layer has constant relhum and N is calculated using virtual temperature.

        CALL tqv_ana_nconstlayers(zml, hsurf, theta_ini(:,:,:), &
             zrhf, ps(:,:,nnew), t_s(:,:,nnew), zrhs)
        ltheta_ini = .TRUE.        

      CASE default

        CALL model_abort (my_world_id, 10008, &
             'ERROR in itype_anaprof_tqv, src_artifdata', &
             'src_artifdata, wrong value of itype_anaprof_tqv '// &
             'encounterd in '//TRIM(yzroutine)//'()')

      END SELECT

    END IF

    !***********************************************************************************
    ! Set U and V analytically (if desired):
    !***********************************************************************************

    IF (itype_artifprofiles == 2 .OR. itype_artifprofiles == 4) THEN

      SELECT CASE (itype_anaprof_uv)

      CASE (1)

        ! Weisman-Klemp tanh-Profile from West to East:
        u (:,:,:,nnew) = u_infty * TANH( (zml_u(:,:,:)-hmin_wk) / (href_wk-hmin_wk) )
        v (:,:,:,nnew) = 0.0_wp

      CASE (2)

        ! Wind profile (from West to East) composed of piecewise 
        ! linear layers, which are not necessarily
        ! continuous. E.g., ramping profile with const. U-layer above.

        CALL u_ana_linearlayers(zml_u, u(:,:,:,nnew))
        v (:,:,:,nnew) = 0.0_wp

        !******************************************************************************
        ! NOTE: THE GRADIENT UGR_LINWIND IS POSITIVE FOR *INCREASING* 
        !       WINDSPEED WITH HEIGHT!
        !******************************************************************************

      CASE (3)

        ! Uniform wind from West to East:
        u (:,:,:,nnew) = u_infty
        v (:,:,:,nnew) = 0.0_wp

      CASE (4)

        ! Uniform wind from South to North:
        u (:,:,:,nnew) = 0.0_wp
        v (:,:,:,nnew) = u_infty

      CASE (5)

        ! Semicircular hodograph, typical for right/leftmoving supercells:
        WRITE (*,*) ' itype_anaprof_uv = 5 not yet implemented ! '
        STOP        

      CASE default

        CALL model_abort (my_world_id, 10008, &
             'ERROR in itype_anaprof_uv, src_artifdata', &
             'src_artifdata, wrong value of itype_anaprof_uv encounterd '// &
             'in '//TRIM(yzroutine)//'()')

      END SELECT
    END IF

    !***********************************************************************************
    ! Generate a thin wind boundary layer with |v_H|=0 at the ground and 
    ! exponent law windprofile:
    !***********************************************************************************
    IF (zo_boundary > 0.0_wp) THEN
      CALL gen_thin_boundary_uv (zml, nnew, ierrstat)
      IF (ierrstat > 0) THEN
        WRITE (*,'(a,i3)') 'ERROR: problem in src_artifdata(): '// &
             'Wind boundary layer interpolation error! ierror = ', ierrstat
        CALL model_abort (my_world_id, 10005, &
             'Wind boundary layer interpolation error', &
             'src_artifdata, setting wind boundary layer')
        RETURN
      END IF
    END IF

    !***********************************************************************************
    ! if desired (NL linitw_followeta=.true.), initialize w in a way that
    ! the flow follows the terrain following coordinate surfaces, otherwise set w = 0:
    !***********************************************************************************
    IF (linitw_followeta) THEN
      CALL init_w_followeta( nnew )
    ELSE
      w(:,:,:,nnew) =  0.0_wp
    END IF


    !------------------------------------------------------------------------------
    ! 6. Initialization of the other thermic variables by hydrostatic balancing
    !    without taking into account coriolis forces and friction, so no geostrophic
    !    balancing or whatsoever.
    !
    !    The initial state should be hydrostatic layered. To get a rather exact
    !    balance of the initial state, the computation of the pressure deviation
    !    is done as in the computation of the pressure gradient in the prognostic
    !    equation for w (=0).
    !
    !------------------------------------------------------------------------------



    ! Computation of t, pp, qv and qc
    ! Problem:
    ! For computing the pressure deviation, the virtual temperature, and 
    ! therefore also the humidities qv and qc are needed. But for computing
    ! qv and qc using the relative humidity, the actual pressure p is needed.
    ! Therefore an iterative method is used here. qv and qc are computed first
    ! with the given pressure deviation (=0 at the beginning) and then the 
    ! resulting deviation is computed. 3 iterations of this kind are 
    ! (hopefully) enough.
    ! The surface pressure, that also is needed for this computation, is taken
    ! either from the sounding file data by hydrostatic polytrope extrapolation 
    ! or as the pressure of the reference atmosphere in case of analytical
    ! profiles of t an rh. In that way, the analytic profiles always get the
    ! surface pressure of the reference atmosphere. Something else could be used here.

    IF ( l2tls ) THEN

      !.. The following is exact for p' T' dynamics. However, this depends also
      !   on irunge_kutta. For now, we apply it here for all cases of the RK-dynamics
      !   regardless of irunge_kutta.

      SELECT CASE (itype_fast_waves)
      CASE (1)
        IF (ltheta_ini) THEN
          ! theta_ini is OPTIONAL; IF present, this is taken instead of temperature,
          ! and temperature is calculated from theta_ini      
          CALL calc_p_hydrostat_psts(ni=ie, nj=je, nk=ke, niter=40, &
               zml=zml, sqrtg=sqrtg_r_w(:,:,1:ke), hsurf=hsurf, psurf=ps(:,:,nnew), &
               t=t(:,:,:,nnew), relhum=zrhf, t0=t0, p0=p0 ,rho0=rho0, &
               piter=pp(:,:,:,nnew), qv=qv(:,:,:), qc=qc(:,:,:), &
               zmaxqv=zmaxqv, r_d=r_d, rvd_m_o=rvd_m_o, &
               itype_fastwaves=itype_fast_waves, &
               theta=theta_ini, ierror=izerror, errmsg=yzerrmsg )
        ELSE
          CALL calc_p_hydrostat_psts(ni=ie, nj=je, nk=ke, niter=40, &
               zml=zml, sqrtg=sqrtg_r_w(:,:,1:ke), hsurf=hsurf, psurf=ps(:,:,nnew), &
               t=t(:,:,:,nnew), relhum=zrhf, t0=t0, p0=p0 ,rho0=rho0, &
               piter=pp(:,:,:,nnew), qv=qv(:,:,:), qc=qc(:,:,:), &
               zmaxqv=zmaxqv, r_d=r_d, rvd_m_o=rvd_m_o, &
               itype_fastwaves=itype_fast_waves, ierror=izerror, errmsg=yzerrmsg )
        END IF
      CASE (2)
        IF (ltheta_ini) THEN
          ! theta_ini is OPTIONAL; IF present, this is taken instead of temperature,
          ! and temperature is calculated from theta_ini
          CALL calc_p_hydrostat_psts(ni=ie, nj=je, nk=ke, niter=40, &
               zml=zml, sqrtg=sqrtg_r_w(:,:,1:ke), hsurf=hsurf, psurf=ps(:,:,nnew), &
               t=t(:,:,:,nnew), relhum=zrhf, t0=t0, p0=p0 ,rho0=rho0, &
               piter=pp(:,:,:,nnew), qv=qv(:,:,:), qc=qc(:,:,:), &
               zmaxqv=zmaxqv, r_d=r_d, rvd_m_o=rvd_m_o, &
               itype_fastwaves=itype_fast_waves, &
               p0hl=p0hl(:,:,1:ke), t0hl=t0hl(:,:,1:ke), wgtfac=wgtfac(:,:,1:ke), &
               theta=theta_ini, ierror=izerror, errmsg=yzerrmsg )
        ELSE
          CALL calc_p_hydrostat_psts(ni=ie, nj=je, nk=ke, niter=40, &
               zml=zml, sqrtg=sqrtg_r_w(:,:,1:ke), hsurf=hsurf, psurf=ps(:,:,nnew), &
               t=t(:,:,:,nnew), relhum=zrhf, t0=t0, p0=p0 ,rho0=rho0, &
               piter=pp(:,:,:,nnew), qv=qv(:,:,:), qc=qc(:,:,:), &
               zmaxqv=zmaxqv, r_d=r_d, rvd_m_o=rvd_m_o, &
               itype_fastwaves=itype_fast_waves, &
               p0hl=p0hl(:,:,1:ke), t0hl=t0hl(:,:,1:ke), wgtfac=wgtfac(:,:,1:ke), &
               ierror=izerror, errmsg=yzerrmsg )
        END IF
      CASE default
        CALL model_abort (my_world_id, 10799, &
             'ERROR, src_artifdata: wrong itype_fastwaves for calc_p_hydrostat_psts()', &
             'src_artifdata, '//TRIM(yzroutine)//'(), hydrostat. pressure init.')
      END SELECT

      IF (izerror /= 0) THEN
        CALL model_abort (my_world_id, izerror, TRIM(yzerrmsg), &
             'gen_ini_data(), calc_p_hydrostat_psts(), calculation of pressure')
      END IF

    ELSE

      IF (ltheta_ini) THEN
        ! theta_ini is OPTIONAL; IF present, this is taken instead of temperature, 
        ! and temperature is calculated from theta_ini      
        CALL calc_p_hydrostat_lf(ni=ie, nj=je, nk=ke, niter=40, &
             zml=zml,hsurf=hsurf, psurf=ps(:,:,nnew), &
             t=t(:,:,:,nnew), relhum=zrhf, t0=t0, p0=p0, dp0=dp0 ,rho0=rho0, &
             piter=pp(:,:,:,nnew), qv=qv(:,:,:), qc=qc(:,:,:), &
             zmaxqv=zmaxqv, r_d=r_d, rvd_m_o=rvd_m_o, &
             theta=theta_ini, ierror=izerror, errmsg=yzerrmsg )                      ! 
      ELSE
        CALL calc_p_hydrostat_lf(ni=ie, nj=je, nk=ke, niter=40, &
             zml=zml, hsurf=hsurf, psurf=ps(:,:,nnew), &
             t=t(:,:,:,nnew), relhum=zrhf, t0=t0, p0=p0, dp0=dp0, rho0=rho0, &
             piter=pp(:,:,:,nnew), qv=qv(:,:,:), qc=qc(:,:,:), &
             zmaxqv=zmaxqv, r_d=r_d, rvd_m_o=rvd_m_o, ierror=izerror, errmsg=yzerrmsg)
      END IF

      IF (izerror /= 0) THEN
        CALL model_abort (my_world_id, izerror, TRIM(yzerrmsg), &
             'gen_ini_data(), calc_p_hydrostat_lf(), calculation of pressure')
      END IF

    END IF

    ! .. If ltheta_ini has not yet been calculated, do it now assuming moist air.
    !    Depending on the namelist switches, it might be needed when
    !    defining temperature disturbances in the initial conditions later on.
    !    NOTE: at this point, pp contains total pressure, not perturbation!
    IF (.NOT.ltheta_ini) THEN
      DO k=1,ke
        DO j=1,je
          DO i=1,ie
            ! moist rd:
            zrdm = rd_moist(qv(i,j,k),qc(i,j,k))
            ! COSMO-approximation of cp:
            zcpm = cp_moist_cosmo(qv(i,j,k),qc(i,j,k),0.0_wp)
            theta_ini(i,j,k) = t(i,j,k,nnew) * (pt00/pp(i,j,k,nnew))**(zrdm/zcpm)
          END DO
        END DO
      END DO
    END IF

    ! .. calc_p_hydrostat_xxx gives pressure, so convert back to perturbation pressure:
    pp(:,:,:,nnew) = pp(:,:,:,nnew) - p0

    ! .. re-compute relative humidity zrhf (for safety):
    DO k=1,ke
      DO j=1,je
        DO i=1,ie
          zrhf(i,j,k) = rh_Tpqv( p0(i,j,k)+pp(i,j,k,nnew), t(i,j,k,nnew), &
                                 qv(i,j,k), qc(i,j,k))
        END DO
      END DO
    END DO

    ! .. control output
    IF (my_cart_id == 0) THEN
      WRITE (*,*)
      WRITE (*,*)  'Initial Profiles on full levels at point (1,1):'
      WRITE (*,'(A3,9A12)') 'k', 'zml', 'T', 'Theta_ini', 'p', &
           'qv', 'qc', 'relhum', 'u', 'v'
      DO k = ke,1,-1
        WRITE(*,'(I3,4F12.3,3F12.8,2F12.5)') k, zml(1,1,k), &
             t(1,1,k,nnew), theta_ini(1,1,k), p0(1,1,k)+pp(1,1,k,nnew), &
             qv(1,1,k), qc(1,1,k), zrhf(1,1,k), &
             u(1,1,k,nnew), v(1,1,k,nnew)
      ENDDO
      WRITE (*,*)

    ENDIF

    ! .. Check Qv about values above 0.1 and issue an error message,
    !    because the COSMO-model is not really designed for such high moistures.
    !    E.g., because c_p is assumed to be that of air, disregarding the
    !    differences bewteen the different hydrometeors.
    DO k=1,ke
      DO j=1,je
        DO i=1,ie
          IF (qv(i,j,k) > 0.1_wp) THEN
            WRITE (*,'(a,3i4,a,es12.5)') 'ERROR: '//TRIM(yzroutine)//': '// &
                 'qv > 0.1 encountered at location ', &
               isubpos(my_cart_id, 1)-nboundlines-1+i, &
               isubpos(my_cart_id, 2)-nboundlines-1+j, k, &
               '  qv(i,j,k) = ',qv(i,j,k) 
            WRITE (*,'(a)') 'IF YOU THINK QV > 0.1 IS OK FOR YOU, ' // &
                 'COMMENT OUT THIS ERROR CHECK IN '//TRIM(YZROUTINE)//'!'
            CALL model_abort (my_world_id, 10789, &
                 'ERROR, src_artifdata: qv > 0.1 encountered', &
                 'src_artifdata, qv > 0.1 in '//TRIM(yzroutine)//'()')
          END IF
        END DO
      END DO
    END DO

    !..  Surface pressure: for consistency the surface pressure is computed
    !    again by extrapolating the pressure from the lowest full level

    CALL calps ( ps(:,:,nnew)   , pp(:,:,ke,nnew), t(:,:,ke,nnew),     &
         qv(:,:,ke), qc(:,:,ke), qrs(:,:,ke),        & 
         rho0(:,:,ke)   , p0(:,:,ke)     , dp0(:,:,ke),        &
         ie, je, rvd_m_o, r_d, 1, ie, 1, je)


    !------------------------------------------------------------------------------
    ! 7. Initialization of soil variables
    !------------------------------------------------------------------------------

    !!! NOTE: If running with soil model (lsoil=.true.), there is an
    !!! extra initialization/correction of some of the soil- and snow variables within the
    !!! routines for the soil model, which is called during the first
    !!! model time step and *after* gen_ini_data! E.g., this relates
    !!! to the snow initialization and soil moisture. 
    !!! If you want to do an idealized run with soil model, you will
    !!! have to carefully check how the soil model modifies your
    !!! below chosen parameters for snow and soil moisture.
    !!! E.g., any initialization value of soil moisture (implicitly it is 0.0) 
    !!! will, depending on your chosen soil type, be set to a minimum of 1 % above
    !!! the socalled air dryness point, which is the soil water content
    !!! which cannot be removed from the soil because of capillary forces.
    !!!
    !!!
    !!! NOTE 1: t_s at this point is the air temperature at the ground deduced
    !!!         from the chosen atmospheric temperature profile.
    !!!
    !!! NOTE 2: in case of periodic BCs, the t_s is alread periodic because of
    !!!         the orography and the atmosphere fields beeing periodic. So no
    !!!         boundary exchange whatsoever is needed.
    

    SELECT CASE (itype_soil_tw)
        
    CASE (1)

        !.. At land points, initialize soil temperature and water saturation
        !   as well as snow parameters and interception store
        !   with constant values t_soil_c, wf_soil_c, t_snow_c, w_snow_c, and
        !   w_i_c. At sea points without ice, initialize t_s with t_water_c. At sea ice points,
        !   initialize t_ice and set t_s with t_ice_c.
        !
        ! NOTE: THE INITIAL SURFACE TEMPERATURE OF SEA OR LAKE POINTS IS THE VARIALBE
        !       t_s(:,:,nnew) AND WILL BE HELD CONSTANT DURING THE SIMULATION.
        !       IF YOU WANT A DISTURBANCE OF THIS TEMPERATURE AT SOME LOCATIONS,
        !       YOU CAN APPLY THE DISTURBANCE TYPE DESIGNED FOR t_s IN THE NAMELIST,
        !       WHICH IS ctype_tempdist = 'hotspot-sfc', OR YOU CAN ORIENT YOURSELF
        !       AT IT TO CREATE YOUR OWN DISTURBANCES IN THE SUBROUTINE set_tempdist_bbc_ts() BELOW.

        CALL init_tw_soil_snow_c()

      
    CASE (2)

      !.. Initialize soil temperature and water saturation
      !   as well as snow parameters and interception store
      !   with 2D fields from ASCII-files
      !   tsurffile, tsoilfile, wfsoilfile, tsnowfile, wsnowfile, wifile.
      !   SAME FORMAT AS FOR OROGRAPHY AND OTHER SOIL PARAMETERS,
      !   SEE ABOVE!
      
      !   "tsurffile"  : "base" temperature t_s of the surface, including all water points and sea ice points.
      !          Might later be overwritten by t_so(0) at land points if using the soil model.
      !   NOTE: if values are < 0.0, then the surface values from the atmosphere remain!

      CALL init_tw_soil_snow_ascii_2d()
      
    CASE default
      
      CALL model_abort (my_world_id, 10009, &
           'ERROR in itype_soil_tw, src_artifdata', &
           'src_artifdata, wrong value of itype_soil_tw encounterd'// &
           ' in '//TRIM(yzroutine)//'()')
      
    END SELECT


    !------------------------------------------------------------------------------
    !    Some cross checks and other computations for internal consistency
    !------------------------------------------------------------------------------

    ! surface temperature of sea ice:
    IF (lseaice) THEN
      WHERE ( .NOT. llandmask .AND. depth_lk <= 0.0_wp .AND. h_ice(:,:,nnew) >= 1e-12_wp )
        t_s(:,:,nnew) = MIN(t_s(:,:,nnew), t0_melt-1.7_wp)
      END WHERE
    END IF

    ! surface temperature of lake ice:
    IF (llake) THEN
      WHERE ( .NOT. llandmask .AND. depth_lk > 0.0_wp .AND. h_ice(:,:,nnew) >= 1e-12_wp )
        t_s(:,:,nnew) = MIN(t_s(:,:,nnew), t0_melt)
      END WHERE
    END IF

    IF (lseaice .OR. llake) THEN
      ! surface temperature of sea water:
      IF (lseaice) THEN
        WHERE ( .NOT. llandmask .AND. depth_lk <= 0.0_wp .AND. h_ice(:,:,nnew) < 1e-12_wp )
          t_s(:,:,nnew) = MAX(t_s(:,:,nnew), t0_melt-1.7_wp)
        END WHERE
      END IF
      ! surface temperature of lake water:
      IF (llake) THEN
        WHERE ( .NOT. llandmask .AND. depth_lk > 0.0_wp .AND. h_ice(:,:,nnew) < 1e-12_wp )
          t_s(:,:,nnew) = MAX(t_s(:,:,nnew), t0_melt)
        END WHERE
      END IF
    ELSE
      ! surface temperature of sea water:
      WHERE ( .NOT. llandmask )
        t_s(:,:,nnew) = MAX(t_s(:,:,nnew), t0_melt-1.7_wp)
      END WHERE
    END IF

    ! snow temperature:
    IF ( lmulti_layer .AND. lmulti_snow ) THEN
      DO k = 0, ke_snow
        ! at water points:
        !    This is not "physically" needed, because there is no
        !    snow at sea or lake points in the model, but t_snow_mult is by convention
        !    set to t_s at water points and is used to detect ice-free water points
        !    in the radiation scheme for the determination of the albedo:
        WHERE (.NOT. llandmask)
          t_snow_mult(:,:,k,nnew) = t_s (:,:,nnew)
        END WHERE
        ! at land points:
        WHERE ( llandmask .AND. w_snow(:,:,nnew) > 0.0_wp )
          t_snow_mult(:,:,k,nnew) = MIN(t_snow_mult(:,:,k,nnew), t0_melt)
        END WHERE
      END DO
    ELSE
      ! at water points:
      !    This is not "physically" needed, because there is no
      !    snow at sea or lake points in the model, but t_snow_mult is by convention
      !    set to t_s at water points and is used to detect ice-free water points
      !    in the radiation scheme for the determination of the albedo:
      WHERE (.NOT. llandmask)
        t_snow(:,:,nnew) = t_s (:,:,nnew)
      END WHERE
      ! at land points:
      WHERE ( llandmask .AND. w_snow(:,:,nnew) > 0.0_wp )
        t_snow(:,:,nnew) = MIN(t_snow(:,:,nnew), t0_melt)
      END WHERE
    END IF

    ! surface temperature of glaciers:
    WHERE ( llandmask .AND. NINT(soiltyp,iintegers) == 1_iintegers )
      t_s(:,:,nnew) = MIN(t_s(:,:,nnew), t0_melt)
    END WHERE

    ! soil water content of glaciers should be 0.0:
    DO k = 1, ke_soil+1
      DO j = 1, je
        DO i = 1, ie
          IF ( llandmask(i,j) .AND. NINT(soiltyp(i,j),iintegers) == 1_iintegers ) THEN
            w_so    (i,j,k,nnew) = 0.0_wp
            w_so_ice(i,j,k,nnew) = 1.0_wp
          END IF
        END DO
      END DO
    END DO

    ! .. If the soil model is running and the soil temperature was not read from an ASCII file, 
    !     enforce at glaciers the "soil" temperature to be the ice surface temperature:
    IF (lsoil .AND. itype_soil_tw == 1) THEN
      IF (lmulti_layer) THEN
        DO k = 0, ke_soil+1
          WHERE (llandmask .AND. NINT(soiltyp,iintegers) == 1_iintegers)
            t_so  (:,:,k,nnew) = t_s(:,:,nnew)
          END WHERE
        END DO
      ELSE
        WHERE (llandmask .AND. NINT(soiltyp,iintegers) == 1_iintegers)
          t_m  (:,:,nnew) =  t_s(:,:,nnew)
        END WHERE
      END IF
      IF (lalloc_t_cl) THEN
        WHERE (llandmask .AND. NINT(soiltyp,iintegers) == 1_iintegers)
          t_cl  (:,:) = t_s(:,:,nnew)
        END WHERE
      END IF
    END IF


    !------------------------------------------------------------------------------
    !    Air Humidity at the surface: qv_s is determined assuming 
    !    constant relative humidity in the lowest model layer and
    !    taking the actual values of t_s:
    !------------------------------------------------------------------------------

    DO j = 1,je
      DO i = 1,ie
        qv_s(i,j,nnew) = qv_Tprelhum( ps(i,j,nnew), t_s(i,j,nnew), &
                                      zrhf(i,j,ke), qc(i,j,ke) )
      ENDDO
    ENDDO


    !------------------------------------------------------------------------------
    ! 6a. Update t_g after t_s has been set or modified
    !------------------------------------------------------------------------------

    ! Compute the temperature at the boundary soil-atmosphere
    IF ( lmulti_layer .AND. lmulti_snow ) THEN
      CALL tgcom ( t_g(:,:,nnew), t_snow_mult(:,:,1,nnew), t_s(:,:,nnew), &
           w_snow(:,:,nnew), llandmask(:,:), ie, je, cf_snow,     &
           1, ie, 1, je )
    ELSE
      CALL tgcom ( t_g(:,:,nnew), t_snow(:,:,nnew), t_s(:,:,nnew),        &
           w_snow(:,:,nnew), llandmask(:,:), ie, je, cf_snow,     &
           1, ie, 1, je )
    ENDIF


    !------------------------------------------------------------------------------
    ! 8. Copy the initial data into timelevel nnow for leapfrog integration
    !------------------------------------------------------------------------------  

    IF ( .NOT.l2tls ) THEN

      u (:,:,:,nnow) = u (:,:,:,nnew)
      v (:,:,:,nnow) = v (:,:,:,nnew)
      w (:,:,:,nnow) = w (:,:,:,nnew)
      t (:,:,:,nnow) = t (:,:,:,nnew)
      pp(:,:,:,nnow) = pp(:,:,:,nnew)
      ps(:,:  ,nnow) = ps(:,:  ,nnew)
      t_s   (:,:,nnow) = t_s (:,: ,nnew)
      t_g   (:,:,nnow) = t_g (:,: ,nnew)
      qv_s  (:,:,nnow) = qv_s(:,: ,nnew)
      IF (lmulti_layer .AND. lmulti_layer) THEN
        t_so(:,:,:,nnow) = t_so(:,:,:,nnew)
        w_so(:,:,:,nnow) = w_so(:,:,:,nnew)
        w_so_ice(:,:,:,nnow) = w_so_ice(:,:,:,nnew)
      ELSE
        t_m(:,:,nnow)  = t_m(:,:,nnew)
        w_g1(:,:,nnow) = w_g1(:,:,nnew)
        w_g2(:,:,nnow) = w_g2(:,:,nnew)
        IF (lalloc_w_g3) THEN
          w_g3(:,:,nnow) = w_g3(:,:,nnew)
        ENDIF
      ENDIF
      IF (lmulti_layer .AND. lmulti_snow) THEN
        t_snow_mult(:,:,:,nnow) = t_snow_mult(:,:,:,nnew)
      ELSE
        t_snow(:,:,nnow) = t_snow (:,: ,nnew)
      ENDIF
      w_snow(:,:,nnow) = w_snow(:,:,nnew)
      w_i(:,:,nnow) = w_i(:,:,nnew)
      IF (lseaice .OR. llake) THEN
        h_ice(:,:,nnow) = h_ice(:,:,nnew)
        t_ice(:,:,nnow) = t_ice(:,:,nnew)
      END IF

      ! loop over tracers
      DO iztrcr = 1, trcr_get_ntrcr()

        ! get pointer to tracer (at nnew)
        CALL trcr_get( izerror, iztrcr, ptr_tlev = nnew, ptr = ztrcr_new )
        IF ( izerror /= 0_iintegers ) THEN
          yzerrmsg = trcr_errorstr( izerror )
          CALL model_abort( my_cart_id, izerror, yzerrmsg, TRIM(yzroutine) )
        ENDIF

        ! get pointer to tracer (at nnow)
        CALL trcr_get( izerror, iztrcr, ptr_tlev = nnow, ptr = ztrcr_now )
        IF ( izerror /= 0_iintegers ) THEN
          yzerrmsg = trcr_errorstr( izerror )
          CALL model_abort( my_cart_id, izerror, yzerrmsg, TRIM(yzroutine) )
        ENDIF
        ! copy new to now 
        ztrcr_now(:,:,:)  = ztrcr_new(:,:,:)

      END DO

    ENDIF

    !------------------------------------------------------------------------------
    !  End of the Subroutine
    !------------------------------------------------------------------------------

  CONTAINS

    ! Check the values of vcoord against other namelist parameters
    SUBROUTINE check_vcoord()

      IMPLICIT NONE

      IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
        WRITE (*,*)  'Subr. check_vcoord() ...'
      END IF

      IF (lspubc) THEN
        ! upper damping layer is active, vcflat has to be smaller than
        ! rdheight in order for the damping layer to function properly!
        IF (rdheight < vcoord%vcflat) THEN
          WRITE (*,*) ' SRC_ARTIFDATA: ERROR in vert. coordinate: '// &
               'vcflat > rdheight ! Break!'
          CALL model_abort (my_world_id, 10005, &
               'ERROR in vert. coordinate specification', &
               'src_artifdata.f90')
          RETURN        
        END IF
      END IF

      ! ..Check strong constraints on vcoord:
      !--------------------------------------
      IF (vcoord%ivctype == 1) THEN
        IF (vcoord%sigm_coord(1) <= 0.0_wp .OR. vcoord%sigm_coord(ke+1) /= 1.0_wp) THEN
          WRITE (*,*) ' SRC_ARTIFDATA: ERROR  *** Wrong values for vcoord: ', &
               vcoord%sigm_coord(1:ke+1),' *** '
          WRITE (*,'(a,i3,a)') '        *** must contain ke+1 increasing values ', &
               'between >0.0 and 1.0 for ivctype = ',vcoord%ivctype,' !!! *** '
          CALL model_abort (my_world_id, 10006, &
               'ERROR in vert. coordinate specification', &
               'src_artifdata.f90')
          RETURN
        END IF
      ELSE IF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        IF (vcoord%vert_coord(ke+1) /= 0.0_wp) THEN
          WRITE (*,*) ' ERROR  *** Wrong value for vcoord(ke+1): ',vcoord%vert_coord(ke+1),' *** '
          WRITE (*,'(a,i3,a)') '        *** must be 0.0 !!! *** '
          CALL model_abort (my_world_id, 10007, &
               'ERROR in vert. coordinate specification', &
               'src_artifdata.f90')
          RETURN
        END IF
      END IF

      ! ..Check settings for the reference atmosphere in relation to vcoord:
      !---------------------------------------------------------------------
      IF (vcoord%ivctype == 1) THEN
        SELECT CASE (refatm%irefatm)
        CASE (1)
          ! Atmosphere with p*dT/dp = const. has a finite height if p*dT/dp > 0
          ! (otherwise there is no finite height). The height is limited 
          ! by T = 0, not by p = 0 in this case!
          IF (refatm%dt0lp > 0.0_wp) THEN
            hmax_refatm = -refatm%t0sl/refatm%dt0lp
            IF (LOG(vcoord%sigm_coord(1)) < hmax_refatm) THEN
              WRITE (*,'(a)') ' SRC_ARTIFDATA: ERROR in vert. coordinate: '// &
                   'vcoord(1) < minimum allowed VALUE ', &
                   EXP(-hmax_refatm),' for reference atmosphere! Break!'
              CALL model_abort (my_world_id, 10011, &
                   'ERROR in vert. coordinate specification', &
                   'src_artifdata.f90')
              RETURN        
            END IF
          END IF
        END SELECT
        ! There are no such constraints for the other reference atmospheres for ivctype = 1!
      ELSE IF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        SELECT CASE (refatm%irefatm)
        CASE (1)
          ! Atmosphere with p*dT/dp = const. has a finite height if p*dT/dp > 0
          ! (otherwise there is no finite height). The height is limited 
          ! by T = 0, not by p = 0 in this case!
          IF (refatm%dt0lp > 0.0_wp) THEN
            ! Here T may become 0 at a certain height (pressure is by definition not 0!)
            hmax_refatm = r_d*refatm%t0sl**2/(2.0_wp*refatm%dt0lp*g)
            IF (vcoord%vert_coord(1) >= hmax_refatm) THEN
              WRITE (*,'(a)') ' SRC_ARTIFDATA: ERROR in vert. coordinate: '// &
                   'vcoord(1) > maximum allowed VALUE ', &
                   hmax_refatm,' for reference atmosphere! Break!'
              CALL model_abort (my_world_id, 10012, &
                   'ERROR in vert. coordinate specification', &
                   'src_artifdata.f90')
              RETURN        
            END IF
          END IF
        CASE (2)
          ! nothing to check, because exponential reference atmosphere does not have
          ! finite height! P and T are always > 0 for any finite height.
        CASE (3)
          ! Atmosphere with const. N has finite height if N < N_tconst = g^2 / (T0 cp)
          ! (otherwise there is no finite height):
          bvref_tconst = SQRT(g**2 / (refatm%t0sl*cp_d))
          IF (refatm%bvref < bvref_tconst) THEN
            hmax_refatm = -g/refatm%bvref**2 * LOG(1.0_wp-(refatm%bvref/bvref_tconst)**2)
            IF (vcoord%vert_coord(1) >= hmax_refatm) THEN
              WRITE (*,'(a)') ' SRC_ARTIFDATA: ERROR in vert. coordinate: '// &
                   'vcoord(1) > maximum allowed VALUE ', &
                   hmax_refatm,' for reference atmosphere! Break!'
              CALL model_abort (my_world_id, 10013, &
                   'ERROR in vert. coordinate specification', &
                   'src_artifdata.f90')
              RETURN        
            END IF
          END IF
        END SELECT
      END IF

    END SUBROUTINE check_vcoord

    ! Subroutine to read external soil parameters from ASCII-Files. Is used only
    ! withtin SR gen_ini_data().
    SUBROUTINE read_soilparams_from_files()

      IMPLICIT NONE

      IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
        WRITE (*,*)  'Subr. read_soilparams_from_files() ...'
      END IF

      CALL read_ascii_field_2d(TRIM(z0file), i_shift_realoro, j_shift_realoro, &
           0.0_wp, gz0, ierrstat)
      IF (my_cart_id == 0) THEN
        IF (ierrstat /= 0 ) THEN
          WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d(), ! Break!'
          CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
               'src_artifdata, reading z0 from file '//TRIM(z0file)// &
               ' with read_ascii_field_2d()')
          RETURN
        END IF
      END IF
      ! The external file contains z0, so multiply by g to get g*z0:
      gz0 = g * gz0

      CALL read_ascii_field_2d(TRIM(frlandfile), i_shift_realoro, j_shift_realoro, &
           0.0_wp, fr_land, ierrstat)
      IF (my_cart_id == 0) THEN
        IF (ierrstat /= 0 ) THEN
          WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
          CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
               'src_artifdata, reading frland from file '//TRIM(frlandfile)// &
               ' with read_ascii_field_2d()')
          RETURN
        END IF
      END IF

      CALL read_ascii_field_2d(TRIM(plcovfile), i_shift_realoro, j_shift_realoro, &
           0.0_wp, plcov, ierrstat)
      IF (my_cart_id == 0) THEN
        IF (ierrstat /= 0 ) THEN
          WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
          CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
               'src_artifdata, reading plcov from file '//TRIM(plcovfile)// &
               ' with read_ascii_field_2d()')
          RETURN
        END IF
      END IF

      CALL read_ascii_field_2d(TRIM(laifile), i_shift_realoro, j_shift_realoro, &
           0.0_wp, lai, ierrstat)
      IF (my_cart_id == 0) THEN
        IF (ierrstat /= 0 ) THEN
          WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
          CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
               'src_artifdata, reading lai from file '//TRIM(laifile)// &
               ' with read_ascii_field_2d()')
          RETURN
        END IF
      END IF

      CALL read_ascii_field_2d(TRIM(soiltypefile), i_shift_realoro, j_shift_realoro, &
           1.0_wp, soiltyp, ierrstat)
      IF (my_cart_id == 0) THEN
        IF (ierrstat /= 0 ) THEN
          WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
          CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
               'src_artifdata, reading soiltyp from file '//TRIM(soiltypefile)// &
               ' with read_ascii_field_2d()')
          RETURN
        END IF
      END IF
        
      IF ( lsoil ) THEN

        CALL read_ascii_field_2d(TRIM(rootdpfile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, rootdp, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
            CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
                 'src_artifdata, reading rootdp from file '//TRIM(rootdpfile)// &
                 ' with read_ascii_field_2d()')
            RETURN
          END IF
        END IF

      END IF

      IF (lforest) THEN

        CALL read_ascii_field_2d(TRIM(forefile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, for_e, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
            CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
                 'src_artifdata, reading for_e from file '//TRIM(forefile)// &
                 ' with read_ascii_field_2d()')
            RETURN
          END IF
        END IF

        CALL read_ascii_field_2d(TRIM(fordfile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, for_d, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
            CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
                 'src_artifdata, reading for_d from file '//TRIM(fordfile)// &
                 ' with read_ascii_field_2d()')
            RETURN
          END IF
        END IF
      ENDIF

      IF (lsso) THEN

        CALL read_ascii_field_2d(TRIM(ssostdhfile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, sso_stdh, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
            CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
                 'src_artifdata, reading sso_stdh from file '//TRIM(ssostdhfile)// &
                 ' with read_ascii_field_2d()')
            RETURN
          END IF
        END IF

        CALL read_ascii_field_2d(TRIM(ssogammafile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, sso_gamma, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
            CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
                 'src_artifdata, reading sso_gamma from file '//TRIM(ssogammafile)// &
                 ' with read_ascii_field_2d()')
            RETURN
          END IF
        END IF

        CALL read_ascii_field_2d(TRIM(ssothetafile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, sso_theta, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
            CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
                 'src_artifdata, reading sso_theta from file '//TRIM(ssothetafile)// &
                 ' with read_ascii_field_2d()')
            RETURN
          END IF
        END IF

        CALL read_ascii_field_2d(TRIM(ssosigmafile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, sso_sigma, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
            CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
                 'src_artifdata, reading sso_sigma from file '//TRIM(ssosigmafile)// &
                 ' with read_ascii_field_2d()')
            RETURN
          END IF
        END IF

      ENDIF

      IF (lseaice .OR. llake) THEN

        CALL read_ascii_field_2d(TRIM(hicefile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, h_ice(:,:,nnew), ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in read_ascii_field_2d()! Break!'
            CALL model_abort (my_world_id, 10005, 'ERROR in read_ascii_field_2d()', &
                 'src_artifdata, reading h_ice from file '//TRIM(hicefile)// &
                 ' with read_ascii_field_2d()')
            RETURN
          END IF
        END IF

      ENDIF

    END SUBROUTINE read_soilparams_from_files


    ! Subroutine to initialize soil temperature and water saturation
    ! as well as snow parameters and interception store
    ! with constant values t_soil_c, wf_soil_c, t_snow_c, w_snow_c, and
    !   w_i_c:
    SUBROUTINE init_tw_soil_snow_c()

      IMPLICIT NONE

      IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
        WRITE (*,*)  'Subr. init_tw_soil_snow_c() ...'
      END IF

      IF (t_surf_c >= 0.0_wp) THEN

        !.. Initialize surface temperature t_s with a constant value differently from the
        !   air temperature at the ground.
        !
        ! NOTE: Might be overwritten later:
        !   - at land points, if lsoil=.true., itype_soil_tw=1 and t_soil_c >= 0
        !   - everywhere if lsoil=.true. and itype_soil_tw=2 (read soil parameters from ASCII files)
        
        t_s(:,:,nnew) = t_surf_c
        
      END IF

      IF (lsoil) THEN
        IF (lmulti_layer) THEN

          ! T in the soil:
          IF (t_soil_c < 0.0_wp) THEN
            DO k = 0, ke_soil+1
              WHERE (llandmask)
                t_so(:,:,k,nnew) = t_s(:,:,nnew)
              END WHERE
            END DO
          ELSE
            DO k = 0, ke_soil+1
              WHERE (llandmask)
                t_so(:,:,k,nnew) = t_soil_c
              END WHERE
            END DO
          ENDIF

          ! Water content of soil layers:
          DO k=1,ke_soil+1
            DO j=1,je
              DO i=1,ie
                IF (llandmask(i,j)) THEN
                !  w_so(i,j,k,nnew) = wf_soil_c * cporv(NINT(soiltyp(i,j),iintegers)) * czmls(k)
                ! CPS BugFix
                 w_so(i,j,k,nnew) = wf_soil_c * cporv(NINT(soiltyp(i,j),iintegers)) * (czhls(k)-czhls(k-1))
                END IF
              END DO
            END DO
            WHERE (llandmask)
              w_so_ice(:,:,k,nnew) = 0.0_wp
            END WHERE
          ENDDO

          ! Weighting function indicating "freshness" of snow in 
          ! upper few cm of snow cover:
          freshsnow(:,:) = 1.0_wp

        ELSE ! .not. lmulti_layer

          ! T in the middle soil layer:
          IF (t_soil_c < 0.0_wp) THEN
            WHERE (llandmask)
              t_m(:,:,nnew) = t_s(:,:,nnew)
            END WHERE
          ELSE
            WHERE (llandmask)
              t_m(:,:,nnew) = t_soil_c
            END WHERE
          END IF

          ! Water saturation of soil layers:
          DO j=1,je
            DO i=1,ie
              IF (llandmask(i,j)) THEN
                w_g1(i,j,nnew) = wf_soil_c * cporv(NINT(soiltyp(i,j),iintegers)) * cdzw13
                w_g2(i,j,nnew) = wf_soil_c * cporv(NINT(soiltyp(i,j),iintegers)) * cdzw23
                w_cl(i,j)      = wf_soil_c * cporv(NINT(soiltyp(i,j),iintegers)) * cdzw33
              END IF
            END DO
          END DO
          IF (lalloc_w_g3) THEN
            DO j=1,je
              DO i=1,ie
                IF (llandmask(i,j)) THEN
                  w_g3(i,j,nnew) = wf_soil_c * cporv(NINT(soiltyp(i,j),iintegers)) * cdzw33
                END IF
              END DO
            END DO
          END IF

        END IF   ! lmulti_layer

        ! Deep soil temperature:
        IF (lalloc_t_cl) THEN
          IF (t_soil_c < 0.0_wp) THEN
            WHERE (llandmask)
              t_cl(:,:) = t_s (:,:,nnew)
            END WHERE
          ELSE
            WHERE (llandmask)
              t_cl(:,:) = t_soil_c
            END WHERE
          END IF
        END IF

        ! Interception storage water on Plants [m H2O]:
        WHERE (llandmask)
          w_i(:,:,nnew) = w_i_c
        END WHERE

      END IF   ! lsoil

      ! Snow temperature:
      IF( lmulti_layer .AND. lmulti_snow ) THEN

        IF (t_snow_c < 0.0_wp) THEN
          DO k = 0, ke_snow
            WHERE (llandmask)
              t_snow_mult(:,:,k,nnew) = t_s (:,:,nnew)
            END WHERE
          END DO
        ELSE
          DO k = 0, ke_snow
            WHERE (llandmask)
              t_snow_mult(:,:,k,nnew) = t_snow_c
            END WHERE
          END DO
        END IF

      ELSE

        IF (t_snow_c < 0.0_wp) THEN
          WHERE (llandmask)
            t_snow(:,:,nnew) = t_s (:,:,nnew)
          END WHERE
        ELSE
          WHERE (llandmask)
            t_snow(:,:,nnew) = t_snow_c
          END WHERE
        ENDIF

      ENDIF

      ! Snow water equivalent [m H2O]:
      WHERE (llandmask)
        w_snow(:,:,nnew) = w_snow_c
      END WHERE

      ! Temperature of water points (might be overwritten at icy points right below):
      IF (t_water_c >= 0.0_wp) THEN
        WHERE (.NOT. llandmask)
          t_s(:,:,nnew) = t_water_c
        END WHERE
        IF (llake) THEN  ! unnecessary?
          WHERE (depth_lk > 0.0_wp)
            t_s(:,:,nnew) = t_water_c
          END WHERE
        END IF
      END IF

      ! Surface temperature of ice at sea- or lake points with ice or at land points with glaciers:
      IF (t_ice_c >= 0.0_wp) THEN
        IF (lseaice .OR. llake) THEN
          ! .. Dynamic ice bodies, like sea ice or lake ice:
          WHERE ( (.NOT. llandmask .OR. depth_lk > 0.0_wp) .AND. h_ice(:,:,nnew) > 0.0_wp)
            t_s  (:,:,nnew) = t_ice_c
            t_ice(:,:,nnew) = t_ice_c
          END WHERE
        END IF
        ! .. Static ice bodies, like glaciers:
        WHERE (llandmask .AND. NINT(soiltyp,iintegers) == 1_iintegers)
          t_s  (:,:,nnew) = t_ice_c
        END WHERE
      END IF

    END SUBROUTINE init_tw_soil_snow_c

    ! Subroutine to initialize soil temperature and water saturation
    ! as well as snow parameters and interception store
    ! with 2D fields from ASCII-files:

    SUBROUTINE init_tw_soil_snow_ascii_2d()

      IMPLICIT NONE

      REAL(KIND=wp)     :: zt_soil(ie,je), zt_snow(ie,je), zwf_soil(ie,je), zt_surf(ie,je), &
                           zw_snow(ie,je)

      IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
        WRITE (*,*)  'Subr. init_tw_soil_snow_ascii_2d() ...'
      END IF

      ! T at the surface: Baseline initialization of t_s from t_surf. Where the values
      !   are < 0.0, the surface values from the atmosphere remain.
      !
      ! NOTE: Should contain the respective surface values of water points and of static ice bodies
      !       like glaciers (soiltyp == 10) and of sea ice points with h_ice > 0.0!
      !
      CALL read_ascii_field_2d(TRIM(tsurffile), i_shift_realoro, j_shift_realoro, &
           0.0_wp, zt_surf, ierrstat)
      IF (my_cart_id == 0) THEN
        IF (ierrstat /= 0 ) THEN
          WRITE (*,*) ' SRC_ARTIFDATA: ERROR in init_tw_soil_snow_ascii_2d()! Break!'
          CALL model_abort (my_world_id, 10035, 'ERROR in init_tw_soil_snow_ascii_2d()', &
               'src_artifdata, reading t_surf from file '//TRIM(tsurffile)// &
               ' with init_tw_soil_snow_ascii_2d()')
          RETURN
        END IF
      END IF

      ! .. set t_s to zt_surf (values from the ASCII-file "tsurffile") at those points
      !    where zt_surf is >= 0.0. Otherwise, the surface values from the
      !    atmosphere remain:
      WHERE (zt_surf >= 0.0_wp)
        t_s(:,:,nnew) = zt_surf
      END WHERE

      ! .. Default initialization:
      zt_soil  = t_s(:,:,nnew)
      zt_snow  = t_s(:,:,nnew)
      zwf_soil = 0.0_wp
      zw_snow  = 0.0_wp

      IF (lsoil) THEN

        ! T in the soil:
        CALL read_ascii_field_2d(TRIM(tsoilfile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, zt_soil, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in init_tw_soil_snow_ascii_2d()! Break!'
            CALL model_abort (my_world_id, 10035, 'ERROR in init_tw_soil_snow_ascii_2d()', &
                 'src_artifdata, reading t_soil from file '//TRIM(tsoilfile)// &
                 ' with init_tw_soil_snow_ascii_2d()')
            RETURN
          END IF
        END IF

        WHERE (zt_soil < 0.0_wp)
          zt_soil(:,:) = t_s(:,:,nnew)
        END WHERE
        
        ! Snow temperature:
        CALL read_ascii_field_2d(TRIM(tsnowfile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, zt_snow, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in init_tw_soil_snow_ascii_2d()! Break!'
            CALL model_abort (my_world_id, 10036, 'ERROR in init_tw_soil_snow_ascii_2d()', &
                 'src_artifdata, reading t_snow from file '//TRIM(tsnowfile)// &
                 ' with init_tw_soil_snow_ascii_2d()')
            RETURN
          END IF
        END IF

        WHERE (zt_snow < 0.0_wp)
          zt_snow(:,:) = t_s(:,:,nnew)
        END WHERE
        
        ! Water saturation of soil layers:
        CALL read_ascii_field_2d(TRIM(wfsoilfile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, zwf_soil, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in init_tw_soil_snow_ascii_2d()! Break!'
            CALL model_abort (my_world_id, 10037, 'ERROR in init_tw_soil_snow_ascii_2d()', &
                 'src_artifdata, reading wf_soil from file '//TRIM(wfsoilfile)// &
                 ' with init_tw_soil_snow_ascii_2d()')
            RETURN
          END IF
        END IF

        IF (lmulti_layer) THEN

          ! T in the soil:
          DO k = 0, ke_soil+1
            t_so(:,:,k,nnew) = zt_soil(:,:)
          END DO

          ! Water content of soil layers:
          DO k=1,ke_soil+1
            DO j=1,je
              DO i=1,ie
                w_so(i,j,k,nnew) = zwf_soil(i,j) * cporv(NINT(soiltyp(i,j),iintegers)) * czmls(k)
              END DO
            END DO
            w_so_ice(:,:,k,nnew) = 0.0_wp
          ENDDO

          ! Weighting function indicating "freshness" of snow in 
          ! upper few cm of snow cover:
          freshsnow(:,:) = 1.0_wp

        ELSE ! .not. lmulti_layer

          ! T in the middle soil layer:
          t_m(:,:,nnew) = zt_soil(:,:)

          ! Water content of soil layers:
          DO j=1,je
            DO i=1,ie
              w_g1(i,j,nnew) = zwf_soil(i,j) * cporv(NINT(soiltyp(i,j),iintegers)) * cdzw13
              w_g2(i,j,nnew) = zwf_soil(i,j) * cporv(NINT(soiltyp(i,j),iintegers)) * cdzw23
              w_cl(i,j)      = zwf_soil(i,j) * cporv(NINT(soiltyp(i,j),iintegers)) * cdzw33
            END DO
          END DO
          IF (lalloc_w_g3) THEN
            DO j=1,je
              DO i=1,ie
                w_g3(i,j,nnew) = zwf_soil(i,j) * cporv(NINT(soiltyp(i,j),iintegers)) * cdzw33
              END DO
            END DO
          END IF

        END IF   ! lmulti_layer

        ! Deep soil temperature:
        IF (lalloc_t_cl) THEN
          t_cl(:,:) = zt_soil(:,:)
        END IF

        ! Snow water equivalent [m H2O]:
        CALL read_ascii_field_2d(TRIM(wsnowfile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, zw_snow, ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in init_tw_soil_snow_ascii_2d()! Break!'
            CALL model_abort (my_world_id, 10038, 'ERROR in init_tw_soil_snow_ascii_2d()', &
                 'src_artifdata, reading w_snow from file '//TRIM(wsnowfile)// &
                 ' with init_tw_soil_snow_ascii_2d()')
            RETURN
          END IF
        END IF
        
        ! Interception storage water on Plants [m H2O]:
        CALL read_ascii_field_2d(TRIM(wifile), i_shift_realoro, j_shift_realoro, &
             0.0_wp, w_i(:,:,nnew), ierrstat)
        IF (my_cart_id == 0) THEN
          IF (ierrstat /= 0 ) THEN
            WRITE (*,*) ' SRC_ARTIFDATA: ERROR in init_tw_soil_snow_ascii_2d()! Break!'
            CALL model_abort (my_world_id, 10038, 'ERROR in init_tw_soil_snow_ascii_2d()', &
                 'src_artifdata, reading w_i from file '//TRIM(wifile)// &
                 ' with init_tw_soil_snow_ascii_2d()')
            RETURN
          END IF
        END IF

      END IF   ! soil

      ! Snow temperature:
      IF( lmulti_layer .AND. lmulti_snow ) THEN

        DO k = 0, ke_snow
          WHERE (llandmask)
            t_snow_mult(:,:,k,nnew) = zt_snow (:,:)
          END WHERE
        END DO
        
      ELSE

        WHERE (llandmask)
          t_snow(:,:,nnew) = zt_snow (:,:)
        END WHERE

      ENDIF

      ! Snow water equivalent [m H2O]:
      w_snow(:,:,nnew) = zw_snow(:,:)

      ! Surface temperature of ice at sea- or lake points with ice:
      IF (lseaice .OR. llake) THEN
        WHERE ( (.NOT. llandmask .OR. depth_lk > 0.0_wp) .AND. h_ice(:,:,nnew) > 0.0_wp)
          t_ice(:,:,nnew) = t_s(:,:,nnew)
        END WHERE
      ENDIF

    END SUBROUTINE init_tw_soil_snow_ascii_2d


  END SUBROUTINE gen_ini_data

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  !
  !+ Module procedure in "src_artifdata" for generating artificial boundary data
  !
  !------------------------------------------------------------------------------

  SUBROUTINE gen_bound_data ( naction )

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   Artificial boundary data are generated, if lartif_data = .TRUE.
    !   Only boundary data that are constant in time are provided, so that
    !   it suffices to copy the initial data to the boundary fields
    !
    ! Method:
    !
    !------------------------------------------------------------------------------
    !
    ! Declarations:
    INTEGER (KIND=iintegers)     :: naction  ! first or following call
    INTEGER (KIND=iintegers)     :: i, j, k, izerror, iztrcr
    REAL(KIND=wp)                :: zsqv
    CHARACTER(len=255)           :: yzerrmsg, yzroutine

    !------------------------------------------------------------------------------
    !- End of header -
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE gen_bound_data
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'gen_bound_data'

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. '//TRIM(yzroutine)//'()', naction, ' ...'
    END IF

    ! update pointers to microphysics tracers (at nnew)
    CALL get_tracers()

    IF (naction == 1) THEN

      !==============================================
      ! Action items for the initial model time step:
      !==============================================

      ! copy the initial data fields to the corresponding boundary data fields
      ! for both boundary levels

      ! fields of the atmosphere
      u_bd (:,:,:,nbd1) = u (:,:,:,nnew)
      u_bd (:,:,:,nbd2) = u (:,:,:,nnew)
      v_bd (:,:,:,nbd1) = v (:,:,:,nnew)
      v_bd (:,:,:,nbd2) = v (:,:,:,nnew)
      IF ( .NOT. lw_freeslip ) THEN
        w_bd (:,:,:,nbd1) = w (:,:,:,nnew)
        w_bd (:,:,:,nbd2) = w (:,:,:,nnew)
      END IF
      t_bd (:,:,:,nbd1) = t (:,:,:,nnew)
      t_bd (:,:,:,nbd2) = t (:,:,:,nnew)
      pp_bd(:,:,:,nbd1) = pp(:,:,:,nnew)
      pp_bd(:,:,:,nbd2) = pp(:,:,:,nnew)

      ! tracers
      DO  iztrcr = 1, trcr_get_ntrcr()

        CALL trcr_get( izerror, iztrcr, ptr_tlev = nnew,     &
                       ptr = ztrcr_new, ptr_bd = ztrcr_bd ) 
        IF ( izerror /= 0_iintegers ) THEN
          yzerrmsg = trcr_errorstr( izerror )
          CALL model_abort( my_cart_id, izerror, yzerrmsg, TRIM(yzroutine) )
        ENDIF

        ztrcr_bd(:,:,:,nbd1) = ztrcr_new
        ztrcr_bd(:,:,:,nbd2) = ztrcr_new

      END DO

      ! soil fields
      IF (lalloc_t_s_bd) THEN
        t_s_bd   (:,:,nbd1) = t_s (:,:,nnew)
        t_s_bd   (:,:,nbd2) = t_s (:,:,nnew)
      ELSE
        ALLOCATE(t_s_bd(ie,je,2))
        t_s_bd = 0.0_wp
        t_s_bd   (:,:,nbd1) = t_s (:,:,nnew)
        t_s_bd   (:,:,nbd2) = t_s (:,:,nnew)                
      END IF
      IF (.NOT. lmulti_layer) THEN
        t_m_bd   (:,:,nbd1) = t_m (:,:,nnew)
        t_m_bd   (:,:,nbd2) = t_m (:,:,nnew)
        w_g1_bd   (:,:,nbd1) = w_g1(:,:,nnew)
        w_g1_bd   (:,:,nbd2) = w_g1(:,:,nnew)
        w_g2_bd   (:,:,nbd1) = w_g2(:,:,nnew)
        w_g2_bd   (:,:,nbd2) = w_g2(:,:,nnew)
        IF (lalloc_w_g3_bd) THEN
          w_g3_bd   (:,:,nbd1) = w_g3(:,:,nnew)
          w_g3_bd   (:,:,nbd2) = w_g3(:,:,nnew)
        ENDIF
      END IF
      IF (lbdclim) THEN
        plcov_bd(:,:,nbd1) = plcov(:,:)
        plcov_bd(:,:,nbd2) = plcov(:,:)
        lai_bd(:,:,nbd1) = lai(:,:)
        lai_bd(:,:,nbd2) = lai(:,:)
        rootdp_bd(:,:,nbd1) = rootdp(:,:)
        rootdp_bd(:,:,nbd2) = rootdp(:,:)
        vio3_bd(:,:,nbd1) = vio3(:,:)
        vio3_bd(:,:,nbd2) = vio3(:,:)
        hmo3_bd(:,:,nbd1) = hmo3(:,:)
        hmo3_bd(:,:,nbd2) = hmo3(:,:)
        t_cl_bd(:,:,nbd1) = t_cl(:,:)
        t_cl_bd(:,:,nbd2) = t_cl(:,:)
        w_cl_bd(:,:,nbd1) = w_cl(:,:)
        w_cl_bd(:,:,nbd2) = w_cl(:,:)
      ENDIF
      IF (lmulti_layer .AND. lmulti_snow) THEN
        t_snow_bd(:,:,nbd1) = t_snow_mult(:,:,1,nnew)
        t_snow_bd(:,:,nbd2) = t_snow_mult(:,:,1,nnew)
      ELSE
        t_snow_bd(:,:,nbd1) = t_snow(:,:,nnew)
        t_snow_bd(:,:,nbd2) = t_snow(:,:,nnew)
      ENDIF
      w_snow_bd(:,:,nbd1) = w_snow(:,:,nnew)
      w_snow_bd(:,:,nbd2) = w_snow(:,:,nnew)
      qv_s_bd  (:,:,nbd1) = qv_s(:,:,nnew)
      qv_s_bd  (:,:,nbd2) = qv_s(:,:,nnew)

      ! the boundary fields of the soil-water contents are initialized with 0.0
      ! together with the allocation 

    ELSEIF (naction == 2) THEN

      !==========================================================
      ! Action items for time dependent boundary data. 
      ! Is called in every time step!
      !==========================================================

      ! Nothing implemented so far!
      ! This means that we have time constant idealized boundary data!
      
    ENDIF


!!! UB: If modeltime < hcond_on, then lcond and lgsp are set .false.,
!!! so that up to this time no condensation or cloud microphysics is active. 
!!! This is intended to help spinning up the flow (e.g., over a mountain)
!!! before convection is triggered.

    IF (hcond_on > 0.0_wp) THEN

      IF (naction == 1) THEN
        lcond_buffer = lcond
        lgsp_buffer = lgsp
      END IF

      IF ( ntstep*dt/3600 < hcond_on ) THEN        
        lcond = .FALSE.
        lgsp = .FALSE.
      ELSE
        lcond = lcond_buffer
        lgsp = lgsp_buffer
        IF ((ntstep-1)*dt/3600 < hcond_on .AND. (lcond .OR. lgsp)) THEN
          ! Limit the relative humidity to 100 % in the very time step of
          ! condensation/microphysics activation, to inhibit strongly saturated
          ! regions from "exploding" cloud developments:

          CALL trcr_get( izerror, idt_qv, ptr_tlev = nnow, ptr = ztrcr_now )
          IF ( izerror /= 0_iintegers ) THEN
            yzerrmsg = trcr_errorstr( izerror )
            CALL model_abort( my_cart_id, izerror, yzerrmsg, 'gen_bound_data' )
          ENDIF

          DO k=1,ke
            DO j=1,je
              DO i=1,ie

                zsqv = qvsat_w(p0(i,j,k)+pp(i,j,k,nnew), t(i,j,k,nnew))
                ! qvsat_w returns -999.99, if E(T) > p, in this case set it to 1.0:
                IF (zsqv < -900.0_wp) zsqv = 1.0_wp
                qv(i,j,k) = MAX(MIN(qv(i,j,k), zsqv), 0.0_wp)

                zsqv = qvsat_w(p0(i,j,k)+pp(i,j,k,nnow), t(i,j,k,nnow))
                ! qvsat_w returns -999.99, if E(T) > p, in this case set it to 1.0:
                IF (zsqv < -900.0_wp) zsqv = 1.0_wp
                ztrcr_now(i,j,k) = MAX(MIN(ztrcr_now(i,j,k), zsqv), 0.0_wp)

              END DO
            END DO
          END DO
        END IF
      END IF

    END IF

    !------------------------------------------------------------------------------
    !  End of the Subroutine
    !------------------------------------------------------------------------------
  END SUBROUTINE gen_bound_data

  !------------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  !
  !+ Module procedure in "src_artifdata" for generating artificial tracers
  !
  !----------------------------------------------------------------------------

  SUBROUTINE gen_trcr_data ( yaction, ierror, yerrmsg )

  !----------------------------------------------------------------------------
  !
  ! Description:
  !   This routine can be used to setup artificial tracer substances which will
  !   be transported, diffused, mixed, etc. along similar to a passive scalar
  !   in the atmosphere
  !
  ! Method:
  !
  !----------------------------------------------------------------------------

  IMPLICIT NONE

  !============================================================================
  !
  ! Parameter list:
  CHARACTER (LEN=*),        INTENT(IN)            ::                      &
    yaction      ! action to be performed

  INTEGER (KIND=iintegers), INTENT(OUT)           ::                      &
    ierror       ! error status

  CHARACTER (LEN=*),        INTENT(OUT)           ::                      &
    yerrmsg      ! error message

  ! Local variables: 

  INTEGER (KIND=iintegers) :: &
    i,j,k,          &  ! loop inidices
    izerr              ! error status

  REAL (KIND=wp),     POINTER :: &
    t_dat(:,:,:) => NULL(),   &  ! tracer data pointer
    t_bd(:,:,:,:)=> NULL()       ! tracer BC pointer

  !----------------------------------------------------------------------------
  !- End of header
  !----------------------------------------------------------------------------
 
  !----------------------------------------------------------------------------
  !- Begin SUBROUTINE gen_trcr_data
  !----------------------------------------------------------------------------

  ierror = 0
  yerrmsg(:) = ' '


  !----------------------------------------------------------------------------
  ! Section 1: Definition of tracers
  !----------------------------------------------------------------------------

  IF (yaction == 'define') THEN
 
    ! init error value
    izerr = 0_iintegers

   ! ! define a new passive tracer
   ! CALL trcr_new( &
   !   ierr           = izerr,                      &
   !   yshort_name    = 'TRCR',                     &
   !   igribparam     = 33,                         &
   !   igribtable     = 2,                          &
   !   yparent        = 'src_artifdata',            &
   !   yunits         = 'kg kg-1',                  &
   !   ystandard_name = 'QIpseudo',                 &
   !   ylong_name     = 'QIpseudo',                 &
   !   itype_adv      = T_ADV_ON,                   &
   !   itype_diff     = T_DIFF_OFF,                 &
   !   itype_turbmix  = T_TURB_OFF,                 &
   !   itype_ini      = T_INI_USER,                 &
   !   itype_lbc      = T_LBC_USER,                 &
   !   itype_relax    = T_RELAX_OFF,                &
   !   itype_damp     = T_DAMP_OFF,                 &
   !   itype_clip     = T_CLP_OFF )
   !
   ! ! check for errors
   ! IF (izerr /= 0_iintegers) THEN
   !   ierror = izerr
   !   yerrmsg = trcr_errorstr( izerr )
   !   RETURN
   ! ENDIF

   ! ! check for errors
   ! IF (izerr /= 0_iintegers) THEN
   !   ierror = izerr
   !   yerrmsg = trcr_errorstr( izerr )
   !   RETURN
   ! ENDIF

  !----------------------------------------------------------------------------
  ! Section 2: Initialization of tracers
  !----------------------------------------------------------------------------

  ELSEIF (yaction == 'init') THEN

    ! init error value
    izerr = 0_iintegers

    ! get pointer to tracer data
   !  CALL trcr_get( izerr, 'TRCR', ptr_tlev=nnew, ptr=t_dat )
   !  IF (izerr /= 0_iintegers) THEN
   !    ierror = izerr
   !    yerrmsg = trcr_errorstr( izerr )
   !    RETURN
   !  ENDIF
   !  CALL trcr_get( izerr, 'TRCR', ptr_bd=t_bd )
   !  IF (izerr /= 0_iintegers) THEN
   !    ierror = izerr
   !    yerrmsg = trcr_errorstr( izerr )
   !    RETURN
   !  ENDIF

   !  ! initialize tracer data
   !  DO k = 1, ke
   !    DO j = 1, je
   !      DO i = 1, ie
   !        t_dat(i,j,k) = 0.0_wp
   !      ENDDO
   !    ENDDO
   !  ENDDO

   !  ! initialize boundary data
   !  DO k = 1, ke
   !    DO j = 1, je
   !      DO i = 1, ie
   !        t_bd(i,j,k,:) = 0.0_wp
   !      ENDDO
   !    ENDDO
   !  ENDDO

  !----------------------------------------------------------------------------
  ! Section 3: All other actions are wrong
  !----------------------------------------------------------------------------

  ELSE

    ierror = 1
    yerrmsg = 'ERROR *** No valid action for gen_trcr_data'
  
  ENDIF

  !----------------------------------------------------------------------------
  !  End of the Subroutine
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------

  END SUBROUTINE gen_trcr_data


  SUBROUTINE add_orography_ana(ierrstat)

    IMPLICIT NONE

    INTEGER(KIND=iintegers), INTENT(out) :: ierrstat

    INTEGER(KIND=iintegers) :: &
         i, j, ii, i_td, j_td, ierrstatloc, kzdims(24)
    REAL(KIND=wp)           :: &
         zdx, zdy, zhmax, zhillcut, tmprotangle, zdax, zday, zdr2, &
         zdi, tmpasym_x, tmpasym_y, zdx_rot, zdy_rot

    REAL(KIND=wp)           :: &
         hsurf_inc(ie,je)            ! Increment for hsurf for single mountains/valleys

    INTEGER(KIND=iintegers) :: izerror, hillsign
    CHARACTER(len=250)      :: yzerrmsg

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. add_orography_ana() ...'
    END IF

    ierrstat    = 0_iintegers
    kzdims(:)   = 0_iintegers
    izerror     = 0_iintegers
    yzerrmsg(:) = ' '
    

    DO ii=1, nhill_max

      IF (lhill(ii)) THEN
        ! Generate a Mountain:

        ! Initialize the hsurf-increment, which is computed below:
        hsurf_inc(:,:) = 0.0_wp

        ! precompute constant parameters:

        ! grid length:
        zdy   = r_earth * dlat * degrad
        ! Absolute value of the height of the hill im m:
        zhmax = ABS(hillheight(ii))
        ! Sign of the hill (is needed below in combine_hsurf)
        hillsign = NINT(SIGN(1.0_wp,hillheight(ii)))
        ! Heights smaller than zhillcut x hillheight are flattened:
        zhillcut = zhmax * zhillcutfact(ii)
        ! Rotation angle of the mountain (clockwise) against the Northern 
        ! direction (angle between the positive
        ! Y-direction and the Y-axis of the mountain)
        tmprotangle = hill_rotangle(ii) * degrad

        ! Set asymmetry parameter, if an asymetric hill type is chosen:
        SELECT CASE (TRIM(ADJUSTL(hill_type(ii))))
        CASE ('bellshaped-asym', 'gauss-asym')
          ! This is an asymetric hill, so use the respective namelist parameters:
          tmpasym_x = hillasym_x(ii)
          tmpasym_y = hillasym_y(ii)
        CASE default
          ! 1.0 leads to symetric hills:
          tmpasym_x = 1.0_wp
          tmpasym_y = 1.0_wp
        END SELECT
        
        SELECT CASE (TRIM(ADJUSTL(hill_type(ii))))

        ! CASE ('yourcase')

          !=====================================================================
          ! 
          ! YOU MIGHT ADD YOUR OWN CASE HERE!
          !
          ! ORIENT YOURSELF AT THE CASES BELOW!
          ! 
          ! NOTE: Here we have to compute the absolute value of an orography
          !       increment for your mountain and store it in the local
          !       field "hsurf_inc". Later, the subroutine "combine_hsurf()"
          !       combines this increment with the already existing orography
          !       in the global field "hsurf", depending on the setting
          !       of the namelist switch "hill_combineaction(ii)".
          !
          !       hill_combineaction(ii) = 1 : add to the previous orography
          !                                    "hsurf" (subtract in the case 
          !                                    of a valley)
          !       hill_combineaction(ii) = 2 : take the maximum of "hsurf_inc"
          !                                    and "hsurf" (minimum in case 
          !                                    of a valley) 
          !       
          !=====================================================================


        CASE ('gauss-2d-simple')

          !=====================================================================
          ! 
          ! NOTE: This case is a SIMPLE EXAMPLE on how to generate your
          ! own orography. 
          !
          ! Gaussian shaped mountain, 2D- mountain ridge of infinite length in y-direction 
          ! with halfwidth hill_width_x.
          !
          ! Relevant namelist parameters:
          ! 
          !   - hill_i(ii)           i-Position of mountain center
          !   - ABS(hillheight(ii))  Mountian height in m (absolute value)
          !   - hill_width_x(ii)     Mountain half width diameter in m
          !
          ! For simplicity, it does *NOT* take into account the namelist
          ! parameters 
          !   hill_j, 'hill_width_y', 'hill_2d', 'hillsideradius_y', 'hill_rotangle', 
          !   'zhillcutfact', any asymmetry-parameters, and maybe others.
          !
          ! To find out, how these are to be implemented, refer to the other cases
          ! below!
          !
          !=====================================================================


          ! global coords of lower left corner of this PE domain:
          ! i-index in total domain
          i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
          ! j-index in total domain
          j_td  = isubpos(my_cart_id, 2) - nboundlines - 1

          DO j = 1, je

            DO i = 1, ie

              ! i-position relative to moutain center in GPs
              zdi  = REAL(i_td-hill_i(ii)+i,wp)

              ! grid length in X-dir. incl. metrical grid stretching factor:
              zdx   = r_earth * dlon * degrad * crlat(j,1)

              !.. Orography increment "hsurf_inc" by this mountain
              !
              ! !!! MUST BE POSITIVE ALSO FOR NEGATIVE HILLHEIGHT(II) !!!
              !
              !  - "zhmax" is computed above and is ABS(hillheight(ii)).
              !    If it is added (hill) or subtracted (valley), or
              !    if MAX() or MIN() operations are applied when combining
              !    it with the global field "hsurf", is determined
              !    later in combine_hsurf() by the sign of hillheight(ii)
              !    and the setting of the namelist switch "hill_combineaction(ii)":

              hsurf_inc(i,j) = zhmax * EXP(-LOG(2.0_wp) * (zdi*zdx/hill_width_x(ii))**2 )

            ENDDO
          ENDDO

        CASE ('bellshaped', 'bellshaped-asym')

          ! Bell shaped mountain:

          IF (lhill_2d(ii)) THEN

            ! 2D- mountain ridge of length hill_width_y with 
            ! halfwidth hill_width_x, oriented in y-Direction with rounded edges:

            ! global coords of lower left corner of this PE domain:
            i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
            j_td  = isubpos(my_cart_id, 2) - nboundlines - 1

            DO j = 1, je
              DO i = 1, ie

                ! compute rectangular distances to the (rotated) main axes
                ! of the hill (zdx_rot, zdy_rot):
                CALL hill_rot_coords(i_td+i, j_td+j, hill_i(ii), hill_j(ii), &
                     tmprotangle, 0.0_wp, zdx_rot, zdy_rot)

                IF (zdx_rot >= 0.0_wp) THEN
                  ! we are on the right flank of the mountain - 
                  ! here the asymmetry-parameter takes effect:
                  zdax = hill_width_x(ii) / tmpasym_x
                ELSE
                  ! we are on the left flank of the mountain:
                  zdax = hill_width_x(ii)
                END IF
                IF (ABS(zdy_rot) <= hill_width_y(ii)) THEN
                  hsurf_inc(i,j) = (zhmax+zhillcut) / (1.0_wp + (zdx_rot/zdax)**2) - zhillcut
                  hsurf_inc(i,j) = MAX(hsurf_inc(i,j), 0.0_wp)
                ELSE
                  IF (hillsideradius_y(ii) > 0.0_wp) THEN
                    zdr2 = ((ABS(zdy_rot)-hill_width_y(ii))/hillsideradius_y(ii))**2 + (zdx_rot/zdax)**2
                    hsurf_inc(i,j) = (zhmax+zhillcut) / (1.0_wp + zdr2) - zhillcut
                    hsurf_inc(i,j) = MAX(hsurf_inc(i,j), 0.0_wp)
                  END IF
                END IF
              ENDDO
            ENDDO

          ELSE

            ! 3D- mountain:

            ! global coords of lower left corner of this PE domain:
            i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
            j_td  = isubpos(my_cart_id, 2) - nboundlines - 1

            DO j = 1, je
              DO i = 1, ie

                ! compute rectangular distances to the (rotated) main axes
                ! of the hill (zdx_rot, zdy_rot):
                CALL hill_rot_coords(i_td+i, j_td+j, hill_i(ii), hill_j(ii), &
                     tmprotangle, 0.0_wp, zdx_rot, zdy_rot)

                IF (zdx_rot >= 0.0_wp) THEN
                  ! we are on the right flank of the mountain - 
                  ! here the asymmetry-parameter takes effect:
                  zdax = hill_width_x(ii) / tmpasym_x
                ELSE
                  ! we are on the left flank of the mountain:
                  zdax = hill_width_x(ii)
                END IF
                IF (zdy_rot >= 0.0_wp) THEN
                  ! we are on the right flank of the mountain - 
                  ! here the asymmetry-parameter takes effect:
                  zday = hill_width_y(ii) / tmpasym_y
                ELSE
                  ! we are on the left flank of the mountain:
                  zday = hill_width_y(ii)
                END IF
                hsurf_inc(i,j) = (zhmax+zhillcut) / &
                     (1.0_wp + (zdx_rot/zdax)**2 + (zdy_rot/zday)**2) - zhillcut
                hsurf_inc(i,j) = MAX(hsurf_inc(i,j), 0.0_wp)
              ENDDO
            ENDDO

          END IF

        CASE ('gauss', 'gauss-asym')

          ! Gaussian shaped mountain

          IF (lhill_2d(ii)) THEN

            ! 2D- mountain ridge of length hill_width_y with 
            ! halfwidth hill_width_x, oriented in y-Direction:

            ! global coords of lower left corner of this PE domain:
            i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
            j_td  = isubpos(my_cart_id, 2) - nboundlines - 1

            DO j = 1, je
              DO i = 1, ie

                ! compute rectangular distances to the (rotated) main axes
                ! of the hill (zdx_rot, zdy_rot):
                CALL hill_rot_coords(i_td+i, j_td+j, hill_i(ii), hill_j(ii), &
                     tmprotangle, 0.0_wp, zdx_rot, zdy_rot)

                IF (zdx_rot >= 0.0_wp) THEN
                  ! we are on the right flank of the mountain - 
                  ! here the asymmetry-parameter takes effect:
                  zdax = hill_width_x(ii) / tmpasym_x
                ELSE
                  ! we are on the left flank of the mountain:
                  zdax = hill_width_x(ii)
                END IF
                IF (ABS(zdy_rot) <= hill_width_y(ii)) THEN
                  hsurf_inc(i,j) = (zhmax+zhillcut) * &
                       EXP(-LOG(2.0_wp) * (zdx_rot/zdax)**2 ) - zhillcut
                  hsurf_inc(i,j) = MAX(hsurf_inc(i,j), 0.0_wp)
                ELSE
                  IF (hillsideradius_y(ii) > 0.0_wp) THEN
                    zdr2 = ((ABS(zdy_rot)-hill_width_y(ii))/hillsideradius_y(ii))**2 + (zdx_rot/zdax)**2
                    hsurf_inc(i,j) = (zhmax+zhillcut) * EXP(-LOG(2.0_wp) * zdr2) - zhillcut
                    hsurf_inc(i,j) = MAX(hsurf_inc(i,j), 0.0_wp)
                  END IF
                END IF
              ENDDO
            ENDDO

          ELSE

            ! A 3D- mountain:

            ! global coords of lower left corner of this PE domain:
            i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
            j_td  = isubpos(my_cart_id, 2) - nboundlines - 1

            DO j = 1, je
              DO i = 1, ie

                ! compute rectangular distances to the (rotated) main axes
                ! of the hill (zdx_rot, zdy_rot):
                CALL hill_rot_coords(i_td+i, j_td+j, hill_i(ii), hill_j(ii), &
                     tmprotangle, 0.0_wp, zdx_rot, zdy_rot)

                IF (zdx_rot >= 0.0_wp) THEN
                  ! we are on the right flank of the mountain - 
                  ! here the asymmetry-parameter takes effect:
                  zdax = hill_width_x(ii) / tmpasym_x
                ELSE
                  ! we are on the left flank of the mountain:
                  zdax = hill_width_x(ii)
                END IF
                IF (zdy_rot >= 0.0_wp) THEN
                  ! we are on the right flank of the mountain - 
                  ! here the asymmetry-parameter takes effect:
                  zday = hill_width_y(ii) / tmpasym_y
                ELSE
                  ! we are on the left flank of the mountain:
                  zday = hill_width_y(ii)
                END IF
                hsurf_inc(i,j) = (zhmax+zhillcut) * &
                     EXP(-LOG(2.0_wp) * ((zdx_rot/zdax)**2 + (zdy_rot/zday)**2) ) - zhillcut
                hsurf_inc(i,j) = MAX(hsurf_inc(i,j), 0.0_wp)
              ENDDO
            ENDDO

          END IF

        CASE default

          WRITE (*,*) 'ERROR: problem in specify_orography_ana(): hill_type = '// &
               TRIM(ADJUSTL(hill_type(ii)))//' not available! Abort!'
          ierrstat = 1
          RETURN

        END SELECT


        ! Finally, combine the actual mountain (hsurf_inc) with the
        ! previous existing orography field (hsurf), using the method
        ! specified by hill_combineaction:
        !   1 = add/subtract hsurf_inc to hsurf (depending on hillsign)
        !   2 = take max/min of hsurf and hsurf_inc (depending on hillsign)
        CALL combine_hsurf(hsurf, hsurf_inc, hillsign, hill_combineaction(ii), ierrstatloc)
        IF (ierrstatloc /= 0) THEN
          ierrstat = 2
          WRITE (*,*) 'ERROR in combine_hsurf(), specify_orography_ana(), src_artifdata.'
          RETURN
        END IF

      END IF

    END DO

    ! enforce periodic boundary conditions for hsurf if required:
    IF (lperi_x .OR. lperi_y .OR. l2dim) THEN
      kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
           ( 0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
           kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh,      &
           lperi_x, lperi_y, l2dim, &
           10000, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
           hsurf )
    END IF

    !----------------------------------------------------------------------------
    ! Internal function(s) / subroutine(s):
    !----------------------------------------------------------------------------

  CONTAINS

    SUBROUTINE combine_hsurf(hfeld, hincfeld, hillsign, action, ierrstat)

      IMPLICIT NONE
      REAL(KIND=wp),     INTENT(inout) :: hfeld(ie,je)
      REAL(KIND=wp),     INTENT(in) :: hincfeld(ie,je)
      INTEGER(KIND=iintegers), INTENT(in) :: action, hillsign
      INTEGER(KIND=iintegers), INTENT(out) :: ierrstat

      ierrstat = 0

      SELECT CASE (action)
      CASE (1)
        hfeld = hfeld + hillsign * hincfeld
      CASE (2)
        IF (hillsign >= 0) THEN
          hfeld = MAX(hfeld, hillsign * hincfeld)
        ELSE
          hfeld = MIN(hfeld, hillsign * hincfeld)
        END IF
      CASE default
        WRITE (*,'(a,i2,a)') 'ERROR in combine_hsurf(), src_artifdata: '// &
             'hill_combineaction = ', action,' not available!'
        ierrstat = 1
        RETURN
      END SELECT

      RETURN
    END SUBROUTINE combine_hsurf

    !----------------------------------------------------------------------------
    ! End of the subroutine
    !----------------------------------------------------------------------------
  END SUBROUTINE add_orography_ana

  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !
  ! Read constant parameter fields (socalled "external data") from simple 
  ! Ascii-files. Is used in SR gen_ini_data() to read orography, z0, landuse data,
  ! soil data, etc., for semi-idealized runs.
  !
  ! FORMAT OF THE FILE:
  !         - 1 header line (arbitrary content)
  !         - 1 line with field dimensions i,j
  !         - 1 long column with the data (first index varies first).
  !
  ! FOR EXAMPLE:
  !
  !>> BEGIN ASCII-FILE:
  ! # roughness length z0 [m]
  ! 461 421
  ! 0.001
  ! 0.002
  ! 0.001
  ! ...
  !<< END ASCII-FILE
  !
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  SUBROUTINE read_ascii_field_2d(dateiname, shift_i, shift_j, &
       feld_default, feld_loc, fehler)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: dateiname
    INTEGER(KIND=iintegers), INTENT(inout) :: fehler
    INTEGER(KIND=iintegers), INTENT(in) :: shift_i, shift_j
    REAL(KIND=wp),     INTENT(in) :: feld_default
    REAL(KIND=wp),     INTENT(out) :: feld_loc(ie,je)

    INTEGER(KIND=iintegers) :: ie_in, je_in, i_tmp, j_tmp, i_offset, j_offset, &
         i, j, freieunit, izerror, mpierror, kzdims(24)
    REAL(KIND=wp)     :: f
    REAL(KIND=wp),     ALLOCATABLE :: feld_tot(:,:)
    CHARACTER (LEN=250) :: yzerrmsg, zeile

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. read_ascii_field_2d() ' // TRIM(ADJUSTL(dateiname)) // '...'
    END IF

    fehler = 0

    IF (my_cart_id == 0 .OR. num_compute == 1) THEN

      CALL get_free_unit (freieunit)

      IF (freieunit == -1) THEN
        WRITE (*,*) 'ERROR: problem in read_ascii_field_2d(): '// &
             'no free file unit available! No file written!'
        fehler=fehler+1
        RETURN
      ENDIF

      OPEN(unit=freieunit, file=TRIM(dateiname), status='old', action='read', iostat=izerror)
      WRITE (*,*) 'READ_ASCII_FIELD_2D: Datei '//TRIM(dateiname)//' is read ...'

      IF (izerror /= 0) THEN
        WRITE (*,*) 'ERROR: problem in read_ascii_field_2d(): error opening '// &
             TRIM(dateiname)//' ! Abort!'
        CLOSE(freieunit)
        CALL release_unit(freieunit)
        fehler=fehler+1
        RETURN
      ENDIF

      ! Format: arbitrary number of header lines, 1 line with field dimension and then one long 
      !         column with the data (first index varies first).

      ! jump over the header lines
      DO
        zeile(:) = ' '
        READ(freieunit, '(a)') zeile
        zeile = ADJUSTL(zeile)
        IF (.NOT.(zeile(1:1) == '#' .OR. zeile(1:1) == '!')) EXIT
      END DO
      READ (zeile, *) ie_in, je_in          ! get field dimensions

    END IF

    IF (num_compute > 1) THEN
      CALL distribute_values(ie_in,    1, 0, imp_integers, icomm_cart, mpierror)
      CALL distribute_values(je_in,    1, 0, imp_integers, icomm_cart, mpierror)
    END IF

    ALLOCATE(feld_tot(ie_tot,je_tot))
    feld_tot = feld_default

    IF (my_cart_id == 0 .OR. num_compute == 1) THEN

      ! Read data point wise and center them around
      ! [ ie_tot/2 + shift_i ; je_tot/2 + shift_j ]

      !?????? The thing with the shift_i and shift_j is not very clear !!!!!!!!!!!!!

      i_offset = (ie_in - ie_tot) / 2 + shift_i
      j_offset = (je_in - je_tot) / 2 + shift_j

      IF (i_offset < 0 .OR. j_offset < 0 .OR. &
           i_offset+ie_tot > ie_in .OR. j_offset+je_tot > je_in) THEN
        WRITE (*,*) 'ERROR: problem in read_ascii_field_2d(): external data set ' // &
             TRIM(dateiname) // ' does not cover the whole model domain! Abort!'
        CALL release_unit(freieunit)
        DEALLOCATE (feld_tot)
        fehler=fehler+1
        RETURN          
      ELSE
        DO j=1,je_in
          DO i=1,ie_in
            READ (freieunit, *) f
            ! Compute corresponding indices in the model grid and if inside the domain,
            ! write the value to the model field:
            i_tmp = i - i_offset
            j_tmp = j - j_offset
            IF (i_tmp > 0 .AND. i_tmp <= ie_tot .AND. j_tmp > 0 .AND. j_tmp <= je_tot) THEN
              feld_tot(i_tmp,j_tmp) = f
            END IF
          END DO
        END DO
      END IF

      CLOSE(freieunit)
      CALL release_unit(freieunit)

    ENDIF

    IF (num_compute > 1) THEN
      CALL distribute_field(feld_tot,ie_tot,je_tot,feld_loc,ie,je, 0, izerror)
      fehler = fehler + izerror
      CALL global_values(fehler,1,'MAX',imp_integers,icomm_cart,0,yzerrmsg, izerror)
      fehler = fehler + izerror
    ELSE
      feld_loc = feld_tot
    END IF

    DEALLOCATE (feld_tot)

    ! enforce periodic boundary conditions or 2dim-exchange if required:
    IF (lperi_x .OR. lperi_y .OR. l2dim) THEN
      izerror     = 0_iintegers
      yzerrmsg(:) = ' '
      kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries                                                  &
           ( 0,  sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
           kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh,      &
           lperi_x, lperi_y, l2dim, &
           10000, .FALSE., ncomm_type, izerror, yzerrmsg,                  &
           feld_loc )
    END IF


    RETURN

  END SUBROUTINE read_ascii_field_2d

  !**********************************************************************************
  !
  ! Read radiosonde data from a text file, interpolate it to the model grid
  ! and distribute it to all processors. Only processor 0 reads the file.
  !
  !
  ! The text file is expected to have a certain format: 
  !   4 header lines, arbitrary number of data lines as follows:
  ! 
  ! #
  ! #
  ! #
  ! P [hPa],    Z [m],     T [K],  Dewp [K], Relhum [%], r [g/kg],  WS [m/s],  WD [deg]
  !  1000.0000      0   300.0000   292.4112    63.1962   14.19878    0.00000   270.0000
  !   988.7517    100   299.1311   292.2245    65.7473   14.19415    0.16660   270.0000
  !   977.5984    200   298.2979   292.0383    68.2801   14.19016    0.33284   270.0000
  !   966.5405    300   297.4824   291.8525    70.8582   14.18678    0.49834   270.0000
  !   955.5781    400   296.6796   291.6672    73.5007   14.18399    0.66274   270.0000
  !   944.7110    500   295.8868   291.4822    76.2191   14.18176    0.82570   270.0000
  !   933.9390    600   295.1022   291.2976    79.0218   14.18008    0.98688   270.0000
  !   923.2617    700   294.3247   291.1134    81.9164   14.17893    1.14594   270.0000
  !   912.6788    800   293.5534   290.9296    84.9094   14.17830    1.30260   270.0000
  !   902.1900    900   292.7876   290.7461    88.0072   14.17817    1.45656   270.0000
  !   891.7948   1000   292.0267   290.5629    91.2161   14.17853    1.60756   270.0000
  !   881.4928   1100   291.2702   290.3800    94.5423   14.17937    1.75536   270.0000
  !   871.2825   1200   290.5169   289.8368    95.7824   13.85286    1.89974   270.0000
  !   861.1619   1300   289.7662   289.0174    95.3386   13.29073    2.04052   270.0000
  !   851.1299   1400   289.0184   288.2001    94.8862   12.74963    2.17751   270.0000
  !   841.1862   1500   288.2734   287.3847    94.4256   12.22870    2.31059   270.0000
  !        ...    ...        ...        ...        ...        ...        ...        ...
  ! 
  ! *NOTE 1*: ONLY THE COLUMNS WITH Z, T, RELHUM, WS AND WD ARE TAKEN, THE REST IS
  ! IGNORED! THE REST MAY SERVE FOR PLOTTING PURPOSES WITH USER PLOT SOFTWARE!
  ! HOWEVER, THERE IS ONE OCCASION WHERE THE COLUMN WITH P IS USED TO INTERPOLATE
  ! THE PRESSURE AT THE TOPOGRAPHY HEIGHT TO SERVE AS STARTING POINT FOR THE
  ! MODELS OWN HYDROSTATIC PRESSURE INITIALIZATION, SEE THE COMMENT ON pflag BELOW!
  ! 
  ! *NOTE 2*: THE COLUMN WITH T MAY ALSO CONTAIN THE POT. TEMPERATURE THETA INSTEAD.
  !           THEN YOU HAVE TO SET THE FLAG thetaflag = .true. BELOW!
  !
  ! 
  !**********************************************************************************
  !
  ! Calling parameters of the routine:
  ! ----------------------------------
  !
  ! rasofile (CHAR):         Name of the text file (input)
  ! zml            :         Model full level heights at mass points (input)
  ! hhl            :         Model half level heights at w points (input)
  ! temperatur     :         Temperature or theta field (output)
  ! relfeuchte     :         Rel. humid. field (output)
  ! ukomp, vkomp   :         U / V- component (output)
  ! psurf          :         Surface pressure (output)
  ! tsurf          :         Surface temperature (output)
  ! rhsurf         :         Rel. humid. at the surface (output)
  !
  ! integral_averaging :     Flag to determine the kind of vertical averaging
  !                          .true. : integral averaging based on linear spline
  !                          .false.: simple linear interpolation
  !
  ! pflag :        :         Flag to determine how to interpret pressure
  !                          values in the radiosonde file:
  !                          .TRUE. : p in the rasofile is accurate and can be trusted.
  !                                   it will be used to determine the surface 
  !                                   pressure and/or to calculate t from theta.
  !                          .false.: the routine calc_p_hydrostat_ana() is used 
  !                                   to integrate the radiosonde pressure from the 
  !                                   lowest level upwards. This will then be
  !                                   used to determine the surface 
  !                                   pressure and/or to calculate t from theta.
  !
  !       HOWEVER: THE LOWEST PRESSURE VALUE IN THE FILE IS ALWAYS TAKEN FOR REAL
  !                EVEN IF PFLAG = .FALSE. !!!
  ! 
  ! thetaflag      :         If .true. indicates that the radiosonde file
  !                          contains potential temp. theta instead of t
  !                          in the temperature column. With this it is
  !                          possible to have theta(z) in the file instead of
  !                          t(z) and still have the correct surface pressure
  !                          and temperature and rel. humid. The output "temperatur" will
  !                          then be theta instead of t!
  !
  !**********************************************************************************

  SUBROUTINE read_raso(rasofile, zml, hhl, temperatur, relfeuchte, &
       ukomp, vkomp, psurf, tsurf, rhsurf, &
       integral_averaging, pflag, thetaflag)

    IMPLICIT NONE

    ! Input/Output parameters:
    CHARACTER(*), INTENT(in) :: rasofile
    ! model full level heights at mass points:
    REAL(KIND=wp),     INTENT(in) :: zml(ie,je,ke)
    ! model half level heights at w points:
    REAL(KIND=wp),     INTENT(in) :: hhl(ie,je,ke+1)

    REAL(KIND=wp),     INTENT(out), DIMENSION(ie,je,ke) :: &
         temperatur, relfeuchte, ukomp, vkomp
    REAL(KIND=wp),     INTENT(out), DIMENSION(ie,je) :: tsurf, rhsurf, psurf
    LOGICAL, INTENT(in) :: integral_averaging, pflag, thetaflag

    ! Local variables:
    REAL(KIND=wp),     ALLOCATABLE, DIMENSION(:) :: &
         pin, zin, tempin, taupin, relfin, spezin, wgin, wrin, uin, vin, qcin, &
         thetain
    INTEGER(KIND=iintegers) :: num_headerlines, freieunit, anz, ios, i, j, k, &
         fehler, ierr, mpierror, index
    REAL(KIND=wp)     :: tvu, tvo, gamma,  &
         hhl_u(ke1), hhl_v(ke1), zml_u(ke), zml_v(ke), tvtmp, hhlgew, &
         speztmp, qctmp
    CHARACTER(len=20) :: zeile
    CHARACTER(len=250) :: errmsg

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. read_raso() ...'
    END IF


    IF (my_cart_id == 0 .OR. num_compute == 1) THEN

      ! Count the number of header- and data lines:
      !============================================

      CALL get_free_unit (freieunit)

      IF (freieunit == -1) THEN
        WRITE (*,*) 'ERROR: problem in read_raso(): no free file unit available! Abort!'
        CALL model_abort (my_world_id, 10001, 'no free file unit available', &
             'read_raso, getting free file unit')
        RETURN
      ENDIF

      OPEN(unit=freieunit, file=TRIM(rasofile), status='old', iostat=ios)
      IF (ios /= 0) THEN
        WRITE (*,*) 'ERROR: problem in read_raso(): no such raso-file! Abort!'
        CALL model_abort (my_world_id, 10002, 'no raso-file '//TRIM(rasofile), &
             'read_raso, opening raso-file')
        RETURN
      END IF

      ! Determine number number of header lines (lines starting with # or !):
      num_headerlines = 0
      DO
        zeile(:) = ' '
        READ(freieunit, '(a)', iostat=ios) zeile
        zeile = ADJUSTL(zeile)
        IF (ios < 0 .OR. .NOT.(zeile(1:1) == '#' .OR. zeile(1:1) == '!')) EXIT
        num_headerlines = num_headerlines + 1_iintegers
      END DO
      ! Add the extra header line, which contains the column description and has
      ! no '#' or '!' at the beginning:
      num_headerlines = num_headerlines + 1_iintegers

      ! Determine number number of data lines:
      !  (NOTE: the above extra header line has already been read!)
      anz = 0_iintegers
      DO
        READ(freieunit, *, iostat=ios)
        IF (ios < 0) EXIT
        anz = anz + 1_iintegers
      END DO

      ! anz is now the number of data lines!
      !=====================================

      REWIND(freieunit)

      ALLOCATE(pin(anz), zin(anz), tempin(anz), taupin(anz), &
           relfin(anz), spezin(anz), wgin(anz), wrin(anz), uin(anz), vin(anz))

      ! Read header lines and forget them:
      DO i=1, num_headerlines
        READ(freieunit,*)
      END DO

      ! Read the data:
      DO i=1, anz
        READ(freieunit,*)                            &
             pin(i), zin(i), tempin(i), taupin(i),   &
             relfin(i), spezin(i), wgin(i), wrin(i)
      END DO

      CLOSE(freieunit)
      CALL release_unit(freieunit)

    END IF

    IF (num_compute > 1) THEN
      ! Distribute the data to all PEs:
      CALL distribute_values(anz,    1, 0, imp_integers, icomm_cart, mpierror)

      IF (my_cart_id /= 0) THEN
        ALLOCATE(pin(anz), zin(anz), tempin(anz), taupin(anz), &
             relfin(anz), spezin(anz), wgin(anz), wrin(anz), uin(anz), vin(anz))
      END IF

      CALL distribute_values(pin,    anz, 0, imp_reals, icomm_cart, mpierror)
      CALL distribute_values(zin,    anz, 0, imp_reals, icomm_cart, mpierror)
      CALL distribute_values(tempin, anz, 0, imp_reals, icomm_cart, mpierror)
      CALL distribute_values(taupin, anz, 0, imp_reals, icomm_cart, mpierror)
      CALL distribute_values(relfin, anz, 0, imp_reals, icomm_cart, mpierror)
      CALL distribute_values(spezin, anz, 0, imp_reals, icomm_cart, mpierror)
      CALL distribute_values(wgin,   anz, 0, imp_reals, icomm_cart, mpierror)
      CALL distribute_values(wrin,   anz, 0, imp_reals, icomm_cart, mpierror)

    END IF

    ! Derive some other quantities from the input profiles:
    uin = wgin * SIN((wrin-180.0_wp)*degrad)
    vin = wgin * COS((wrin-180.0_wp)*degrad)
    ! Convert relative humidity from % to "normal" number:
    relfin = relfin * 1.0E-2_wp
    ! Compute specific humidity from the input mixing ratio, convert
    ! from [g/kg] to [kg/kg] before:
    spezin = spezin * 1.0E-3_wp
    spezin = spezin / (1.0_wp + spezin)

    ! Allocate work vector for qc and set it to 0.0. It may be set to something different
    ! if we do not trust the pressure in the radiosonde file (pflag=.false.)
    ! and if there is somewhere supersaturation in the profile:
    ALLOCATE(qcin(anz))
    qcin = 0.0_wp

    IF (thetaflag) THEN
      ALLOCATE(thetain(anz))
      thetain = tempin
    END IF

    IF (.NOT.pflag) THEN
      !.. In this case, the pressure from the rasofile cannot be trusted,
      !   therefore it is recomputed hydrostatically from the height levels
      !   of the radiosonde data, starting with the lowermost pressure
      !   value in the file. Here, condensation at supersaturated levels
      !   is taken into account, which would be not the case
      !   if pflag=.true.
      IF  (thetaflag) THEN
        CALL calc_p_hydrostat_ana(ni=1, nj=1, nk=anz, niter=30, &
             zml=zin, psurf=(/pin(1)*1.0E2_wp/), &
             t=tempin, relhum=relfin, &
             piter=pin, qv=spezin, qc=qcin, &
             zmaxqv=HUGE(1.0_wp), r_d=r_d, rvd_m_o=rvd_m_o, theta=thetain, &
             ierror=ierr, errmsg=errmsg)
      ELSE
        CALL calc_p_hydrostat_ana(ni=1, nj=1, nk=anz, niter=20, &
             zml=zin, psurf=(/pin(1)*1.0E2_wp/), &
             t=tempin, relhum=relfin, &
             piter=pin, qv=spezin, qc=qcin, &
             zmaxqv=HUGE(1.0_wp), r_d=r_d, rvd_m_o=rvd_m_o, &
             ierror=ierr, errmsg=errmsg)
      END IF
      IF (ierr /= 0) THEN
        CALL model_abort (my_world_id, ierr, TRIM(errmsg), &
             'read_raso, calc_p_hydrostat_ana(), calculation of pressure')
      END IF


    ELSE
      ! Convert pin from hPa -> Pa:
      pin = pin * 1.0E2
      IF (thetaflag) THEN
        ! Compute t from theta (needed later for computing
        ! the surface values):
        tempin = thetain * (pin/pt00)**rdocp
      END IF
    END IF

    ! tempin is the temperature in any case, even if theta has
    ! been read from the file. However, the actual content of the
    ! output variable "temperatur"  will
    ! be determined by thetaflag!


    !==================================================================
    !
    ! Interpolation to model levels:
    !
    !==================================================================

    ! catch errors:
    !==============
    ! 1) height range of rasodata is not sufficient:
    IF (MAXVAL(zin) < MAXVAL(zml)) THEN
      WRITE (*,*) 'ERROR: problem in read_raso(): uppermost data point in raso-file'//&
           TRIM(rasofile)//' is too low! Abort!'
      CALL model_abort (my_world_id, 10002, 'raso-file '//TRIM(rasofile), &
           'read_raso, uppermost data point in raso-file too low')
      RETURN
    END IF
    IF (MINVAL(zin) > MINVAL(hhl(:,:,ke+1))) THEN
      WRITE (*,*) 'ERROR: problem in read_raso(): lowermost data point in raso-file'//&
           TRIM(rasofile)//' is too high! Abort!'
      CALL model_abort (my_world_id, 10002, 'raso-file '//TRIM(rasofile), &
           'read_raso, lowermost data point in raso-file too high')
      RETURN
    END IF

    fehler = 0

    !.. Two kinds of possible interpolation methods:
    !===============================================

    IF (integral_averaging) THEN

      ! ... Integral averaging based on the assumption that the 
      !     profile data are the nodes of a linear spline:

      DO j = 1, je
        DO i = 1, ie
          hhl_u = 0.5_wp * (hhl(i,j,:) + hhl(MIN(i+1,ie),j,:))
          hhl_v = 0.5_wp * (hhl(i,j,:) + hhl(i,MIN(j+1,je),:))
          IF (thetaflag) THEN
            CALL mittel_integral_vec(zin, thetain, anz, hhl(i,j,ke+1:1:-1), &
                 temperatur(i,j,ke:1:-1), ke, ierr)
          ELSE            
            CALL mittel_integral_vec(zin, tempin, anz, hhl(i,j,ke+1:1:-1), &
                 temperatur(i,j,ke:1:-1), ke, ierr)
          END IF
          IF (ierr /= 0) fehler = 1
          CALL mittel_integral_vec(zin, relfin, anz, hhl(i,j,ke+1:1:-1), &
               relfeuchte(i,j,ke:1:-1), ke, ierr)
          IF (ierr /= 0) fehler = 1
          CALL mittel_integral_vec(zin, uin, anz, hhl_u(ke+1:1:-1), &
               ukomp(i,j,ke:1:-1), ke, ierr)
          IF (ierr /= 0) fehler = 1
          CALL mittel_integral_vec(zin, vin, anz, hhl_v(ke+1:1:-1), &
               vkomp(i,j,ke:1:-1), ke, ierr)
          IF (ierr /= 0) fehler = 1
        ENDDO
      ENDDO

      IF (fehler /= 0) THEN
        WRITE (*,*) 'ERROR: problem in read_raso(): error in integral averaging!'
        CALL model_abort (my_world_id, 10002, 'integral averaging error', &
             'read_raso, integral averaging')
      END IF

    ELSE

      ! ... ordinary linear interpolation:

      DO j = 1, je
        DO i = 1, ie
          zml_u = 0.5_wp * (zml(i,j,:) + zml(MIN(i+1,ie),j,:))
          zml_v = 0.5_wp * (zml(i,j,:) + zml(i,MIN(j+1,je),:))
          IF (thetaflag) THEN
            CALL linear_interpol_vec(zin, thetain, anz, zml(i,j,:), &
                 temperatur(i,j,:), ke, ierr)
          ELSE
            CALL linear_interpol_vec(zin, tempin, anz, zml(i,j,:), &
                 temperatur(i,j,:), ke, ierr)
          END IF
          IF (ierr /= 0) fehler = 1
          CALL linear_interpol_vec(zin, relfin, anz, zml(i,j,:), &
               relfeuchte(i,j,:), ke, ierr)
          IF (ierr /= 0) fehler = 1
          CALL linear_interpol_vec(zin, uin, anz, zml_u, &
               ukomp(i,j,:), ke, ierr)
          IF (ierr /= 0) fehler = 1
          CALL linear_interpol_vec(zin, vin, anz, zml_v, &
               vkomp(i,j,:), ke, ierr)
          IF (ierr /= 0) fehler = 1
        END DO
      END DO

      IF (fehler /= 0) THEN
        WRITE (*,*) 'ERROR: problem in read_raso(): linear interpolation error!'
        CALL model_abort (my_world_id, 10003, 'linear interpolation error', &
             'read_raso, linear interpolation')
      END IF

    END IF

    !.. Compute surface values of pressure, temperature and relative humidity:
    !=========================================================================

    DO i = 1, ie
      DO j = 1, je
        
        ! Search for the height interval of the surface height within
        ! the profile data:
        index = -999
        DO k = 1, anz-1
          IF (zin(k) <= hhl(i,j,ke+1) .AND. zin(k+1) > hhl(i,j,ke+1)) THEN
            index = k
            EXIT
          END IF
        END DO
        IF (index == -999) THEN
          WRITE (*,*) 'ERROR: problem in read_raso(): calculation of surface pressure!'
          CALL model_abort (my_world_id, 10002, &
               'model topography outside raso height range, ' // &
               'modify radiosonde data and start again!', &
               'read_raso, surface pressure calculation')
        END IF
        
        hhlgew = (hhl(i,j,ke+1)-zin(index)) / (zin(index+1)-zin(index))
        tvu = tempin(index) * (1.0_wp + rvd_m_o*spezin(index) - qcin(index))
        tvo = tempin(index+1) * (1.0_wp + rvd_m_o*spezin(index+1) - qcin(index+1))
        IF (tvu /= tvo) THEN 
          ! Using polytrope atmosphere w.r.t. Tv (piecewise linear Tv profile):
          gamma = (tvu-tvo) / (zin(index+1)-zin(index))
          psurf(i,j) = pin(index) * (1.0_wp - &
               gamma/tvu*(hhl(i,j,ke+1)-zin(index)))**(g/(r_d*gamma))
          ! Linear interpolation of the specific humidity and 
          ! virtual temperature onto the surface height,
          ! from this compute surface temperature "tsurf" and
          ! surface relative humidity "relhumsurf":
          speztmp = spezin(index) + (spezin(index+1)-spezin(index)) * hhlgew
          qctmp = qcin(index) + (qcin(index+1)-qcin(index)) * hhlgew
          tvtmp = tvu + (tvo - tvu) * hhlgew
          tsurf(i,j) = tvtmp / (1.0_wp + rvd_m_o*speztmp - qctmp)
        ELSE
          ! Virt. temperature is constant, the barometric height formula applies:
          psurf(i,j) = pin(index) * EXP(-g*(hhl(i,j,ke+1)-zin(index))/(r_d*tvu))
          speztmp = spezin(index) + (spezin(index+1)-spezin(index)) * hhlgew
          qctmp = qcin(index) + (qcin(index+1)-qcin(index)) * hhlgew
          tsurf(i,j) = tvu / (1.0_wp + rvd_m_o*speztmp - qctmp)
        END IF
        rhsurf(i,j) = rh_Tpqv( psurf(i,j), tsurf(i,j), speztmp, qctmp)
      END DO
    END DO

    ! Transformation of the u- and v-components from the geographical
    ! reference system to the rotated spherical grid of COSMO:
    IF (lmetr) THEN
      DO k = 1, ke
        CALL uv2uvrot_vec(ukomp(:,:,k), vkomp(:,:,k), rlat*raddeg, rlon*raddeg, &
             pollat, pollon, ie, je)
      END DO
    END IF

    DEALLOCATE(pin, zin, tempin, taupin, &
         relfin, spezin, wrin, wgin, uin, vin, qcin)
    IF (ALLOCATED(thetain)) DEALLOCATE(thetain)

  END SUBROUTINE read_raso


  !**********************************************************************************
  !
  ! Generate a thin boundary layer with |v_H|=0 at the ground and exponent law windprofile,
  ! up to a height specified by namelist parameter "zo_boundary".
  ! The exponent of the power law is taken from the namelist parameter 
  ! "exponent_windprof_boundary".
  !
  !**********************************************************************************

  SUBROUTINE gen_thin_boundary_uv (zml, ntlev, ierror)

    IMPLICIT NONE

    ! model full layer heights at mass points:
    REAL(KIND=wp),           INTENT(in)  :: zml(ie,je,ke)  
    INTEGER(KIND=iintegers), INTENT(in)  :: ntlev
    INTEGER(KIND=iintegers), INTENT(out) :: ierror

    INTEGER(KIND=iintegers) :: i, j, k
    REAL(KIND=wp)           :: hsurf_u, hsurf_v, zml_u(ke), zml_v(ke), zzo_boundary(2), &
         zuo_boundary(1), zvo_boundary(1) 

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. gen_thin_boundary_uv() ...'
    END IF


    DO j=1,je
      DO i=1,ie

        ! Compute hsurf and zml at u-points:
        hsurf_u = 0.5_wp * (hsurf(i,j)+hsurf(MIN(i+1,ie),j))
        zml_u = 0.5_wp * (zml(i,j,:)+zml(MIN(i+1,ie),j,:))
        ! Compute hsurf and zml at v-points:
        hsurf_v = 0.5_wp * (hsurf(i,j)+hsurf(i,MIN(j+1,je)))
        zml_v = 0.5_wp * (zml(i,j,:)+zml(i,MIN(j+1,je),:))
        ! Linearily interpolate u and v to the specified thin boundary layer height:
        zzo_boundary(1) = hsurf_u + zo_boundary
        zzo_boundary(2) = hsurf_v + zo_boundary
        CALL linear_interpol(zml_u(ke:1:-1), u(i,j,ke:1:-1,ntlev), &
             ke, zzo_boundary(1), zuo_boundary(1), 1, ierror)
        IF (ierror /= 0) THEN
          ierror = 1
          RETURN
        END IF
        CALL linear_interpol(zml_v(ke:1:-1), v(i,j,ke:1:-1,ntlev), &
             ke, zzo_boundary(2), zvo_boundary(1), 1, ierror)
        IF (ierror /= 0) THEN
          ierror = 2
          RETURN
        END IF

        ! Replace u, v below zo_boundary with power law profile down to v=0 at the surface:
        DO k=ke,1,-1
          IF (zml_u(k) < zzo_boundary(1)) THEN
            u(i,j,k,ntlev) = zuo_boundary(1) * &
                 ( (zml_u(k) - hsurf_u) / zo_boundary )**exponent_windprof_boundary
          END IF
          IF (zml_v(k) < zzo_boundary(2)) THEN
            v(i,j,k,ntlev) = zvo_boundary(1) * &
                 ( (zml_v(k) - hsurf_v) / zo_boundary )**exponent_windprof_boundary
          END IF

        END DO
      END DO
    END DO

  END SUBROUTINE gen_thin_boundary_uv


  !**********************************************************************************
  !
  ! Subroutine to calculate w from u and v in a way that the flow follows the
  ! terrain following coordinate surfaces.
  !
  ! *NOTE*: In case of steep orography and/or strong horizontal flow, this might
  ! be inconsistent with the hydrostatic pressure initialization!
  !
  !**********************************************************************************

  SUBROUTINE init_w_followeta(ntlev)

    IMPLICIT NONE

    INTEGER(KIND=iintegers), INTENT(in ) :: ntlev

    INTEGER(KIND=iintegers) :: i, j, k, izerror, kzdims(24)
    REAL(KIND=wp)           :: zml, zmlm1, u_below, u_above, v_below, v_above, &
         u_w, v_w, zdx, zdy
    CHARACTER (LEN=250) :: yzerrmsg

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. init_w_followeta() ...'
    END IF

    w(:,:,:,ntlev) =  0.0_wp   
    DO k=2,ke
      DO j=1,je
        DO i=1,ie
          ! only from 2 to ke, since in other parts of the model,
          ! at ke1 (lower boundary) resp. 1 (upper boundary)
          ! the vertical velocity is diagnostically taken to be model levels following
          ! via the lower and upper BC on w.

          !.. model full levels:
          !   mass point below w-point (half level):
          zml   = 0.5_wp*( hhl(i,j,k) + hhl(i,j,k+1) )
          !   mass point above w-point:
          zmlm1 = 0.5_wp*( hhl(i,j,k-1) + hhl(i,j,k) )

          ! Interpolate u and v to the mass grid points:
          IF (i > 1) THEN
            u_below = 0.5_wp * (u(i-1,j,k  ,ntlev) + u(i,j,k  ,ntlev))
            u_above = 0.5_wp * (u(i-1,j,k-1,ntlev) + u(i,j,k-1,ntlev))
          ELSE
            u_below = u(i,j,k  ,ntlev)
            u_above = u(i,j,k-1,ntlev)
          END IF
          IF (j > 1) THEN
            v_below = 0.5_wp * (v(i,j-1,k  ,ntlev) + v(i,j,k  ,ntlev))
            v_above = 0.5_wp * (v(i,j-1,k-1,ntlev) + v(i,j,k-1,ntlev))
          ELSE
            v_below = v(i,j,k  ,ntlev)
            v_above = v(i,j,k-1,ntlev)
          END IF

          ! Linearily interpolate u, v to the staggered w-gridpoints (half levels):
          u_w = u_below + (u_above-u_below) / (zmlm1-zml) * (hhl(i,j,k)-zml)
          v_w = v_below + (v_above-v_below) / (zmlm1-zml) * (hhl(i,j,k)-zml)

          ! grid lengths (dlon, dlat with correct metrical terms):
          zdx   = r_earth * dlon * degrad * crlat(j,1)
          zdy   = r_earth * dlat * degrad

          ! w = V * GRAD(hhl) = u*dh/dx + v*dh/dy :x
          ! centered difference (second order) on the inner gridpoints,
          ! first order differencing on the outer gridpoints:

          u_w = u_w * ( hhl(MIN(i+1,ie),j,k) - hhl(MAX(i-1,1),j,k) ) / &
               ( zdx * ( REAL( MIN(i+1,ie)-MAX(i-1,1), wp) ) )
          v_w = v_w * ( hhl(i,MIN(j+1,je),k) - hhl(i,MAX(j-1,1),k) ) / &
               ( zdy * ( REAL( MIN(j+1,je)-MAX(j-1,1), wp) ) )

          w(i,j,k,ntlev) = u_w + v_w

        END DO
      END DO
    END DO

    ! enforce periodic boundary conditions or 2dim-exchange if required:
    IF (lperi_x .OR. lperi_y .OR. l2dim) THEN
      izerror     = 0_iintegers
      yzerrmsg(:) = ' '
      kzdims(1:24)=(/ke1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries &
           (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
           kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
           lperi_x, lperi_y, l2dim, &
           20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
           w(:,:,:,ntlev) )
    END IF

  END SUBROUTINE init_w_followeta

  !=============================================================================
  !
  ! Layered polytrope atmosphere (temperature and relhum) with linear
  ! temperature gradients (TEMPERATURE, NOT VIRTUAL TEMPERATURE).
  !
  ! The number of layers is determined from namelist parameter nlayers_poly
  ! (up to nlayers_poly_max layers). 
  !
  ! nlayers_poly_max is a constant parameter, available globally within 
  ! module src_artifdata,
  ! defining the maximum vector lengths for the below namelist parameters.
  !
  ! The relevant namelist parameters are the following, assuming nlayers_poly = 5:
  !
  ! h_poly(1)           ! lower boundary of lowest polytrope layer [m]
  ! h_poly(2)           ! upper boundary of lowest polytrope layer [m]
  ! h_poly(3)           ! upper boundary of second polytrope layer [m]
  ! h_poly(4)           ! upper boundary of third polytrope layer  [m]
  ! h_poly(5)           ! upper boundary of fourth polytrope layer [m]
  !                     !    at the same time the lower boundary 
  !                     !    of the fifth layer, wich stretches to the model top.
  !
  ! p_base_poly         ! Pressure at height h_poly(1), i.e., at the lower boundary
  !                     ! of the profile specification
  !
  ! t_poly(1)           ! temperature for the height z = h_poly(1) [K]
  ! t_poly(2)           ! temperature for the height z = h_poly(2) [K]
  ! t_poly(3)           ! temperature for the height z = h_poly(3) [K]
  ! t_poly(4)           ! temperature for the height z = h_poly(4) [K]
  ! t_poly(5)           ! temperature for the height z = h_poly(5) [K]
  !
  ! tgr_poly(1)         ! temperature gradient for the height h_poly(1) <= z < h_poly(2) [K/m]
  ! tgr_poly(2)         ! temperature gradient for the height h_poly(2) <= z < h_poly(3) [K/m]
  ! tgr_poly(3)         ! temperature gradient for the height h_poly(3) <= z < h_poly(4) [K/m]
  ! tgr_poly(4)         ! temperature gradient for the height h_poly(4) <= z < h_poly(5) [K/m]
  ! tgr_poly(5)         ! temperature gradient for the height h_poly(5) <= z             [K/m]
  !
  ! rh_poly(1)          ! relhum for the z = h_poly(1) [-]
  ! rh_poly(2)          ! relhum for the z = h_poly(2) [-]
  ! rh_poly(3)          ! relhum for the z = h_poly(3) [-]
  ! rh_poly(4)          ! relhum for the z = h_poly(4) [-]
  ! rh_poly(5)          ! relhum for the z = h_poly(5) [-]
  !
  ! rhgr_poly(1)        ! relhum gradient for the height h_poly(1) <= z < h_poly(2) [1/m]
  ! rhgr_poly(2)        ! relhum gradient for the height h_poly(2) <= z < h_poly(3) [1/m]
  ! rhgr_poly(3)        ! relhum gradient for the height h_poly(3) <= z < h_poly(4) [1/m]
  ! rhgr_poly(4)        ! relhum gradient for the height h_poly(4) <= z < h_poly(5) [1/m]
  ! rhgr_poly(5)        ! relhum gradient for the height h_poly(5) <= z             [1/m]
  !
  !*****************************************************************************
  ! NOTE: THE GRADIENTS TGR_POLY RESP. RHGR_POLY ARE POSITIVE FOR *DECREASING* TEMPERATURE
  !       RESP. RELHUM WITH HEIGHT!
  !*****************************************************************************
  !
  !*****************************************************************************
  ! NOTE 2: IT IS POSSIBLE TO SPECIFY DISCONTINUOUS TEMPERATURE- AND MOISTURE
  !         PROFILES (DISCONTINUOUS AT THE LAYER BORDER HEIGHTS).
  !*****************************************************************************
  !
  !=============================================================================

  SUBROUTINE tqv_ana_polylayers(zml, hsurf, temp, relhum, psurf, tsurf, rhsurf)

    IMPLICIT NONE
    
    REAL(KIND=wp),     INTENT(in) :: zml(ie,je,ke)          ! model full layer heights at mass points
    REAL(KIND=wp),     INTENT(in) :: hsurf(ie,je)           ! model orography height
    REAL(KIND=wp),     INTENT(out) :: temp(ie,je,ke), relhum(ie,je,ke)
    REAL(KIND=wp),     INTENT(out) :: psurf(ie,je), tsurf(ie,je), rhsurf(ie,je)

    ! Local parameters:
    INTEGER(KIND=iintegers) :: i, j, k, index, nlay, ke_interp, ierror
    INTEGER(KIND=iintegers), PARAMETER :: kmax_interp = 1000

    REAL(KIND=wp),     PARAMETER :: dz_interp = 100.0_wp

    REAL(KIND=wp)     :: tvu, tvo, gamma, &
         zinterp(kmax_interp), tinterp(kmax_interp), relhuminterp(kmax_interp), &
         pinterp(kmax_interp), qvinterp(kmax_interp), qcinterp(kmax_interp), ztmp, &
         hsurfmin, hsurfmax

    CHARACTER(len=250) :: errmsg

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. tqv_ana_polylayers() ...'
    END IF

    nlay = nlayers_poly

    !===========================================================================
    !===========================================================================
    ! Layered polytrope atmosphere for temperature: 
    !===========================================================================
    !===========================================================================

    !===========================================================================
    ! check consistencies:
    !===========================================================================

    hsurfmin = MINVAL(hsurf)

    IF (hsurfmin < h_poly(1)) THEN
      WRITE (*,*) 'ERROR: problem in tqv_ana_polylayers(): h_poly(1) > surface height!'
      CALL model_abort (my_world_id, 10101, &
           'h_poly(1) > lowest model layer, decrease h_poly(1)!', &
           'tqv_ana_polylayers(), definition of h_poly(1)')
      RETURN
    END IF
    IF (nlay > 1) THEN
      IF ( ANY( h_poly(2:nlay) - h_poly(1:nlay-1) < 0.0_wp ) ) THEN
        WRITE (*,*) 'ERROR: problem in tqv_ana_polylayers(): '// &
             'h_poly must be monotonically increasing!'
        CALL model_abort (my_world_id, 10102, &
             'h_poly must be monotonically increasing!', &
             'tqv_ana_polylayers(), definition of h_poly')
        RETURN
      END IF
    END IF

    !===========================================================================
    !.. Set T and relhum in the atmosphere and at the ground:
    !===========================================================================

    !.. 1) Set T and relhum within the first nlay-1 polytrope layers:
    DO k=1, nlay-1
      WHERE ( zml(:,:,:) >= h_poly(k) .AND. zml(:,:,:) < h_poly(k+1) )
        temp(:,:,:)   = t_poly(k)  - tgr_poly(k)  * ( zml(:,:,:) - h_poly(k) )
        relhum(:,:,:) = rh_poly(k) - rhgr_poly(k) * ( zml(:,:,:) - h_poly(k) )
      ENDWHERE
    END DO
    !.. 2) Set T and relhum for the top polytrope layer, i.e., above h_poly(nlay):
    WHERE ( zml(:,:,:) >= h_poly(nlay) )
      temp(:,:,:)   = t_poly(nlay)  - tgr_poly(nlay)  * ( zml(:,:,:) - h_poly(nlay) )
      relhum(:,:,:) = rh_poly(nlay) - rhgr_poly(nlay) * ( zml(:,:,:) - h_poly(nlay) )
    ENDWHERE

    !.. 3) Set T and relhum at the surface:
    DO k=1, nlay-1
      WHERE ( hsurf(:,:) >= h_poly(k) .AND. hsurf(:,:) < h_poly(k+1) )
        tsurf(:,:)  = t_poly(k)  - tgr_poly(k)  * ( hsurf(:,:) - h_poly(k) )
        rhsurf(:,:) = rh_poly(k) - rhgr_poly(k) * ( hsurf(:,:) - h_poly(k) )
      ENDWHERE
    END DO
    WHERE ( hsurf(:,:) >= h_poly(nlay) )
      tsurf(:,:)  = t_poly(nlay)  - tgr_poly(nlay)  * ( hsurf(:,:) - h_poly(nlay) )
      rhsurf(:,:) = rh_poly(nlay) - rhgr_poly(nlay) * ( hsurf(:,:) - h_poly(nlay) )
    ENDWHERE

    !===========================================================================
    !.. Interpolate surface pressure:
    !===========================================================================

    !   1) Helper height vector with distance dz_interp, containing
    !      also the layer border heights h_poly and expanding up to hsurfmax
    !      (important because of the possibility to have discontinuous T/QV-profiles):
    hsurfmax = MAXVAL(hsurf)
    ke_interp = 0
    layers: DO i=1, nlay-1
      inner: DO k=0, FLOOR( (h_poly(i+1)-h_poly(i)) / dz_interp )
        ztmp = h_poly(i)+k*dz_interp
        IF (ztmp < h_poly(i+1) .AND. ztmp <= hsurfmax) THEN
          ke_interp = ke_interp + 1
          IF (ke_interp > kmax_interp) EXIT layers
          zinterp(ke_interp) = ztmp
        END IF
      END DO inner
    END DO layers

    DO k=0, FLOOR( (hsurfmax-h_poly(nlay)) / dz_interp )
      ke_interp = ke_interp + 1
      IF (ke_interp <= kmax_interp) THEN
        zinterp(ke_interp) = h_poly(nlay) + k * dz_interp
      END IF
    END DO
    IF (zinterp(ke_interp) <= hsurfmax .OR. ke_interp < 2) THEN
      ke_interp = ke_interp + 1
      IF (ke_interp <= kmax_interp) THEN
        zinterp(ke_interp) = zinterp(ke_interp-1) + dz_interp
      END IF
    END IF
      
    IF (ke_interp > kmax_interp) THEN
      WRITE (*,*) 'ERROR: problem in tqv_ana_polylayers(): ke_interp > kmax_interp!'
      CALL model_abort (my_world_id, 10102, 'ke_interp > kmax_interp!', &
           'tqv_ana_polylayers(), calculation of surface pressure')
      RETURN
    END IF

    !   2) Set T and relhum for the heights of the helper vector:
    DO k=1, nlay-1
      WHERE ( zinterp(1:ke_interp) >= h_poly(k) .AND. &
           zinterp(1:ke_interp) < h_poly(k+1) )
        tinterp(1:ke_interp)  = t_poly(k) - tgr_poly(k) * (zinterp(1:ke_interp)-h_poly(k))
        relhuminterp(1:ke_interp) = rh_poly(k) - rhgr_poly(k) * &
             (zinterp(1:ke_interp)-h_poly(k))
      ENDWHERE
    END DO
    WHERE ( zinterp(1:ke_interp) >= h_poly(nlay) )
      tinterp(1:ke_interp)  = t_poly(nlay) - tgr_poly(nlay) * &
           (zinterp(1:ke_interp)-h_poly(nlay))
      relhuminterp(1:ke_interp) = rh_poly(nlay) - rhgr_poly(nlay) * &
           (zinterp(1:ke_interp)-h_poly(nlay))
    ENDWHERE
    
    !.. 3) Use piecewise polytrope atmosphere pressure as best estimate 
    !      of pressure profile (lowest level: isothermal):
    CALL calc_p_hydrostat_ana(ni=1, nj=1, nk=ke_interp, niter=20, &
         zml=zinterp(1:ke_interp), psurf=(/p_base_poly/), &
         t=tinterp(1:ke_interp), relhum=relhuminterp(1:ke_interp), &
         piter=pinterp(1:ke_interp), qv=qvinterp(1:ke_interp), qc=qcinterp(1:ke_interp), &
         zmaxqv=HUGE(1.0_wp), r_d=r_d, rvd_m_o=rvd_m_o, &
         ierror=ierror, errmsg=errmsg)

    IF (ierror /= 0) THEN
      CALL model_abort (my_world_id, ierror, TRIM(errmsg), &
           'tqv_ana_polylayers(), calc_p_hydrostat_ana(), calculation of surface pressure')
    END IF

    !.. 4) Interpolate surface pressure from best estimated pressure profile. 
    !      The interpolation is based on analytically calculating the pressure 
    !      in a polytrope atmosphere w.r.t. 
    !      virtual temperature:
    !      (uses virtual temperature to include the effects of moisture, neglects
    !      the temperature dependence of cp, as is done in the rest of COSMO)
    DO i = 1, ie
      DO j = 1, je
        
        ! Search the height interval within the helper vector: 
        index = -999
        DO k = 1, ke_interp-1
          IF (zinterp(k) <= hsurf(i,j) .AND. zinterp(k+1) > hsurf(i,j)) THEN
            index = k
            EXIT
          END IF
        END DO
        IF (index == -999) THEN
          WRITE (*,*) 'ERROR: in tqv_ana_polylayers(): calculation of surface pressure!'        
          WRITE (*,'(a,5(1x,f10.1))') 'ERROR:', hsurf(i,j), MAXVAL(hsurf), MINVAL(hsurf), &
               h_poly(1), h_poly(nlay)
          CALL model_abort (my_world_id, 10002, &
               'model topography outside polylayers-height range, ' // &
               'modify namelist parameters and start again!', &
               'tqv_ana_polylayers(), surface pressure calculation')
        END IF
        
        tvu = tinterp(index) * (1.0_wp + rvd_m_o*qvinterp(index) - qcinterp(index))
        tvo = tinterp(index+1) * (1.0_wp + rvd_m_o*qvinterp(index+1) - qcinterp(index+1))
        IF (tvu /= tvo) THEN 
          ! Polytrope atmosphere w.r.t. Tv (linear Tv-profile) within the layer:
          gamma = (tvu-tvo) / (zinterp(index+1)-zinterp(index))
          psurf(i,j) = pinterp(index) * &
               (1.0_wp - gamma/tvu*(hsurf(i,j)-zinterp(index)))**(g/(r_d*gamma))
        ELSE
          ! Const. virt. temperature, barometric height formula:
          psurf(i,j) = pinterp(index) * EXP(-g*(hsurf(i,j)-zinterp(index))/(r_d*tvu))
        END IF
      END DO
    END DO


    RETURN
  END SUBROUTINE tqv_ana_polylayers


  !=============================================================================
  !
  ! Atmosphere with layers exhibiting constant Brunt-Vaisala-Frequency N 
  ! (N takes water vapor relhum into account).
  !
  ! It is assumed that N refers to a moist subsaturated atmosphere, i.e., it contains
  ! the influence of the water vapor on buoyancy):
  !
  ! N^2 = (g / theta_f) * (d_theta_f / dz)
  !   theta_f = T * (1000 hpa/p)^(R_f/cp_f)
  !   R_f = Gas constant of moist air
  !   cp_f = cp of moist air (FOR COSMO: ASSUMPTION cp_f = cp_d)
  !
  ! The number of layers is determined from namelist parameter nlayers_nconst
  ! (up to nlayers_nconst_max layers). 
  !
  ! nlayers_nconst_max is a constant parameter, available globally within module src_artifdata,
  ! setting the maximum vector lengths for the below namelist parameters.
  !
  ! The relevant namelist parameters are the following, assuming nlayers_nconst = 5:
  !
  ! h_nconst(1)           ! lower boundary of lowest polytrope layer [m]
  ! h_nconst(2)           ! upper boundary of lowest polytrope layer [m]
  ! h_nconst(3)           ! upper boundary of second polytrope layer [m]
  ! h_nconst(4)           ! upper boundary of third polytrope layer  [m]
  ! h_nconst(5)           ! upper boundary of fourth polytrope layer [m]
  !                       !    at the same time the lower boundary 
  !                       !    of the fifth layer, wich stretches to the model top.
  !
  ! p_base_nconst         ! Pressure at height h_nconst(1), i.e., at the lower boundary
  !                       ! of the profile specification
  !
  ! theta0_base_nconst    ! pot. temperature for the height z = h_nconst(1) [K]
  !
  ! N_nconst(1)           ! N for the z = h_nconst(1) [-]
  ! N_nconst(2)           ! N for the z = h_nconst(2) [-]
  ! N_nconst(3)           ! N for the z = h_nconst(3) [-]
  ! N_nconst(4)           ! N for the z = h_nconst(4) [-]
  ! N_nconst(5)           ! N for the z = h_nconst(5) [-]
  !
  ! rh_nconst(1)          ! relhum for the z = h_nconst(1) [-]
  ! rh_nconst(2)          ! relhum for the z = h_nconst(2) [-]
  ! rh_nconst(3)          ! relhum for the z = h_nconst(3) [-]
  ! rh_nconst(4)          ! relhum for the z = h_nconst(4) [-]
  ! rh_nconst(5)          ! relhum for the z = h_nconst(5) [-]
  !
  ! rhgr_nconst(1)        ! relhum gradient for the height h_nconst(1) <= z < h_nconst(2) [1/m]
  ! rhgr_nconst(2)        ! relhum gradient for the height h_nconst(2) <= z < h_nconst(3) [1/m]
  ! rhgr_nconst(3)        ! relhum gradient for the height h_nconst(3) <= z < h_nconst(4) [1/m]
  ! rhgr_nconst(4)        ! relhum gradient for the height h_nconst(4) <= z < h_nconst(5) [1/m]
  ! rhgr_nconst(5)        ! relhum gradient for the height h_nconst(5) <= z             [1/m]
  !
  !*****************************************************************************
  ! NOTE: THE GRADIENT RHGR_NCONST IS POSITIVE FOR *DECREASING* RELHUM WITH HEIGHT!
  !*****************************************************************************
  !
  !*****************************************************************************
  ! NOTE 2: IT IS POSSIBLE TO SPECIFY A DISCONTINUOUS MOISTURE
  !         PROFILE (DISCONTINUOUS AT THE LAYER BORDER HEIGHTS).
  !*****************************************************************************
  !
  !=============================================================================

  SUBROUTINE tqv_ana_nconstlayers(zml, hsurf, theta, relhum, psurf, tsurf, rhsurf)

    IMPLICIT NONE
    
    REAL(KIND=wp),     INTENT(in) :: zml(ie,je,ke)          ! model full layer heights at mass points
    REAL(KIND=wp),     INTENT(in) :: hsurf(ie,je)           ! model orography height
    REAL(KIND=wp),     INTENT(out) :: theta(ie,je,ke), relhum(ie,je,ke)
    REAL(KIND=wp),     INTENT(out) :: psurf(ie,je), tsurf(ie,je), rhsurf(ie,je)

    ! Local parameters:
    INTEGER(KIND=iintegers) :: i, j, k, index, nlay, ke_interp, ierror
    INTEGER(KIND=iintegers), PARAMETER :: kmax_interp = 7001

    REAL(KIND=wp),     PARAMETER :: dz_interp = 10.0_wp

    REAL(KIND=wp)     :: tvu, tvo, gamma, theta0_h_nconst(nlayers_nconst_max), &
         zinterp(kmax_interp), tinterp(kmax_interp), relhuminterp(kmax_interp), &
         pinterp(kmax_interp), qvinterp(kmax_interp), qcinterp(kmax_interp), hhlgew, &
         thetasurf(ie,je), thetainterp(kmax_interp), qvtmp, qctmp, hsurfmin, hsurfmax

    CHARACTER(len=250) :: errmsg

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. tqv_ana_nconstlayers() ...'
    END IF

    nlay = nlayers_nconst

    theta0_h_nconst = 0.0_wp

    !===========================================================================
    !===========================================================================
    ! Atmosphere with layers of const. N: 
    !===========================================================================
    !===========================================================================

    !===========================================================================
    ! check consistencies:
    !===========================================================================

    hsurfmin = MINVAL(hsurf)

    IF (hsurfmin < h_nconst(1)) THEN
      WRITE (*,*) 'ERROR: problem in tqv_ana_nconstlayers(): h_nconst(1) > surface layer!'
      CALL model_abort (my_world_id, 10101, &
           'h_nconst(1) > lowest model layer, decrease h_nconst(1)!', &
           'tqv_ana_nconstlayers(), definition of h_nconst(1)')
      RETURN
    END IF
    IF (nlay > 1) THEN
      IF ( ANY( h_nconst(2:nlay) - h_nconst(1:nlay-1) < 0.0_wp ) ) THEN
        WRITE (*,*) 'ERROR: problem in tqv_ana_nconstlayers(): '// &
             'h_nconst must be monotonically increasing!'
        CALL model_abort (my_world_id, 10102, 'h_nconst must be monotonically increasing!', &
             'tqv_ana_nconstlayers(), definition of h_nconst')
        RETURN
      END IF
    END IF

    !===========================================================================
    !.. Set T and relhum in the atmosphere and at the surface:
    !===========================================================================

    !.. 1) Vector of the theta's at the lower boundaries of the layers:
    !      (are needed for steps 2, 3, and 4 below)
    theta0_h_nconst(1) = theta0_base_nconst
    DO k=2,nlay
      theta0_h_nconst(k) = theta0_h_nconst(k-1) * &
           EXP(N_nconst(k-1)**2/g  * (h_nconst(k)-h_nconst(k-1)))
    END DO

    !.. 2) Set T and relhum within the first nlay-1 layers:
    DO k=1, nlay-1
      WHERE ( zml(:,:,:) >= h_nconst(k) .AND. zml(:,:,:) < h_nconst(k+1) )
        theta(:,:,:)  = theta0_h_nconst(k) * &
             EXP( N_nconst(k)**2/g  * (zml(:,:,:)-h_nconst(k)) )
        relhum(:,:,:) = rh_nconst(k) - rhgr_nconst(k) * ( zml(:,:,:) - h_nconst(k) )
      ENDWHERE
    END DO
    !.. 3) Set T and relhum for the top layer, i.e., above h_nconst(nlay):
    WHERE ( zml(:,:,:) >= h_nconst(nlay) )
      theta(:,:,:)  = theta0_h_nconst(nlay) * &
           EXP( N_nconst(nlay)**2/g  * (zml(:,:,:)-h_nconst(nlay)) )
      relhum(:,:,:) = rh_nconst(nlay) - rhgr_nconst(nlay) * &
           ( zml(:,:,:) - h_nconst(nlay) )
    ENDWHERE

    !.. 4) Set T and relhum at the surface:
    DO k=1, nlay-1
      WHERE ( hsurf(:,:) >= h_nconst(k) .AND. hsurf(:,:) < h_nconst(k+1) )
        thetasurf(:,:)  = theta0_h_nconst(k) * &
             EXP( N_nconst(k)**2/g  * (hsurf(:,:) -h_nconst(k)) )
        rhsurf(:,:) = rh_nconst(k) - rhgr_nconst(k) * ( hsurf(:,:) - h_nconst(k) )
      ENDWHERE
    END DO
    WHERE ( hsurf(:,:) >= h_nconst(nlay) )
      thetasurf(:,:)  = theta0_h_nconst(nlay) * &
           EXP( N_nconst(nlay)**2/g  * (hsurf(:,:)-h_nconst(nlay)) )
      rhsurf(:,:) = rh_nconst(nlay) - rhgr_nconst(nlay) * ( hsurf(:,:) - h_nconst(nlay) )
    ENDWHERE

    !===========================================================================
    !.. Interpolate surface pressure:
    !===========================================================================


    !.. 1) Helper height vector with distance dz_interp:
    !      Simpler solution as with tqv_ana_polylayers() above, because
    !      the T-profile is continuous (only QV might be discontinuous,
    !      but this will be not so critical here):
    hsurfmax = MAXVAL(hsurf)
    ke_interp = MAX( CEILING( (hsurfmax-h_nconst(1)) / dz_interp ) + 2, 2_iintegers)
    IF (ke_interp > kmax_interp) THEN
      WRITE (*,*) 'ERROR: problem in tqv_ana_nconstlayers(): ke_interp > kmax_interp!'
      CALL model_abort (my_world_id, 10102, 'ke_interp > kmax_interp!', &
           'tqv_ana_nconstlayers(), calculation of surface pressure')
      RETURN
    END IF
    layers: DO k=1, ke_interp
      zinterp(k) = h_nconst(1) + (k-1)*dz_interp
    END DO layers

    !.. 2) Set T and relhum for the heights of the helper vector:
    DO k=1, nlay-1
      WHERE ( zinterp(1:ke_interp) >= h_nconst(k) .AND. &
           zinterp(1:ke_interp) < h_nconst(k+1) )
        thetainterp(1:ke_interp)  = theta0_h_nconst(k) * &
             EXP( N_nconst(k)**2/g  * (zinterp(1:ke_interp)-h_nconst(k)) )
        relhuminterp(1:ke_interp) = rh_nconst(k) - rhgr_nconst(k) * &
             (zinterp(1:ke_interp)-h_nconst(k))
      ENDWHERE
    END DO
    WHERE ( zinterp(1:ke_interp) >= h_nconst(nlay) )
      thetainterp(1:ke_interp)  = theta0_h_nconst(nlay) * &
           EXP( N_nconst(nlay)**2/g  * (zinterp(1:ke_interp)-h_nconst(nlay)) )
      relhuminterp(1:ke_interp) = rh_nconst(nlay) - rhgr_nconst(nlay) * &
           (zinterp(1:ke_interp)-h_nconst(nlay))
    ENDWHERE
    
    !.. 3) Use piecewise polytrope atmosphere pressure as best estimate of pressure profile (lowest level: isothermal):
    CALL calc_p_hydrostat_ana(ni=1, nj=1, nk=ke_interp, niter=30, &
         zml=zinterp(1:ke_interp), psurf=(/p_base_nconst/), &
         t=tinterp(1:ke_interp), relhum=relhuminterp(1:ke_interp), &
         piter=pinterp(1:ke_interp), qv=qvinterp(1:ke_interp), qc=qcinterp(1:ke_interp), &
         zmaxqv=HUGE(1.0_wp), r_d=r_d, rvd_m_o=rvd_m_o, &
         theta=thetainterp(1:ke_interp), &
         ierror=ierror, errmsg=errmsg)

    IF (ierror /= 0) THEN
      CALL model_abort (my_world_id, ierror, TRIM(errmsg), &
           'tqv_ana_nconstlayers(), calc_p_hydrostat_ana(), calculation of surface pressure')
    END IF

    !.. 4) Interpolate surface pressure from best estimated pressure profileThe interpolation
    !   is based on analytically calculating the pressure in a polytrope atmosphere w.r.t. 
    !   virtual temperature:
    !     (uses virtual temperature to include the effects of moisture, neglects
    !      the temperature dependence of cp, as is done in the rest of COSMO)
    DO i = 1, ie
      DO j = 1, je
        
        ! Get the index of the height interval within the helper vector:         
        index = FLOOR((hsurf(i,j)-h_nconst(1))/dz_interp) + 1

        tvu = tinterp(index) * (1.0_wp + rvd_m_o*qvinterp(index) - qcinterp(index))
        tvo = tinterp(index+1) * (1.0_wp + rvd_m_o*qvinterp(index+1) - qcinterp(index+1))
        IF (tvu /= tvo) THEN 
          ! Polytrope atmosphere w.r.t. Tv (linear Tv-profile) within the layer:
          gamma = (tvu - tvo) / (zinterp(index+1)-zinterp(index))
          psurf(i,j) = pinterp(index) * &
               (1.0_wp - gamma/tvu*(hsurf(i,j)-zinterp(index)))**(g/(r_d*gamma))
        ELSE
          ! Const. virt. temperature, barometric height formula:
          psurf(i,j) = pinterp(index) * EXP(-g*(hsurf(i,j)-zinterp(index))/(r_d*tvu))
        END IF

        !.. Calculate surface temperature from surface pressure and potential temperature:
        hhlgew = (hsurf(i,j)-zinterp(index)) / (zinterp(index+1)-zinterp(index))
        qvtmp = qvinterp(index) + (qvinterp(index+1)-qvinterp(index)) * hhlgew
        qctmp = qcinterp(index) + (qcinterp(index+1)-qcinterp(index)) * hhlgew
        tsurf(i,j) = thetasurf(i,j) * (psurf(i,j) / pt00 )** &
             (rd_moist(qvtmp,qctmp)/cp_moist_cosmo(qvtmp,qctmp,0.0_wp))

      END DO
    END DO

    RETURN
  END SUBROUTINE tqv_ana_nconstlayers

  !=============================================================================
  !
  ! Wind profile composed of piecewise linear layers (not necessarily continuous!).
  !
  ! The number of layers is determined by namelist parameter nlayers_linwind
  ! (up to nlayers_lindwind_max). 
  !
  ! nlayers_linwind_max is a constant parameter, available globally within module src_artifdata,
  ! setting the maximum vector lengths for the below namelist parameters.
  !
  ! The relevant namelist parameters are the following, assuming nlayers_linwind = 5:
  !
  ! h_linwind(1)           ! lower boundary of lowest linear wind layer [m]
  ! h_linwind(2)           ! upper boundary of lowest linear wind layer [m]
  ! h_linwind(3)           ! upper boundary of second linear wind layer [m]
  ! h_linwind(4)           ! upper boundary of third linear wind layer  [m]
  ! h_linwind(5)           ! upper boundary of fourth linear wind layer [m]
  !                        !    at the same time the lower boundary 
  !                        !    of the fifth layer, wich stretches to the model top.
  !
  ! u_linwind(1)           ! windspeed for the height z = h_linwind(1) [m/s]
  ! u_linwind(2)           ! windspeed for the height z = h_linwind(2) [m/s]
  ! u_linwind(3)           ! windspeed for the height z = h_linwind(3) [m/s]
  ! u_linwind(4)           ! windspeed for the height z = h_linwind(4) [m/s]
  ! u_linwind(5)           ! windspeed for the height z = h_linwind(5) [m/s]
  !
  ! ugr_linwind(1)         ! windspeed gradient for the height h_linwind(1) <= z < h_linwind(2) [1/m]
  ! ugr_linwind(2)         ! windspeed gradient for the height h_linwind(2) <= z < h_linwind(3) [1/m]
  ! ugr_linwind(3)         ! windspeed gradient for the height h_linwind(3) <= z < h_linwind(4) [1/m]
  ! ugr_linwind(4)         ! windspeed gradient for the height h_linwind(4) <= z < h_linwind(5) [1/m]
  ! ugr_linwind(5)         ! windspeed gradient for the height h_linwind(5) <= z                [1/m]
  !
  !***********************************************************************************
  ! NOTE: THE GRADIENT UGR_LINWIND IS POSITIVE FOR *INCREASING* WINDSPEED WITH HEIGHT!
  !***********************************************************************************
  !
  !=============================================================================

  SUBROUTINE u_ana_linearlayers(zml_u, uloc)

    IMPLICIT NONE

    ! model full layer heights at u-points:
    REAL(KIND=wp),     INTENT(in) :: zml_u(ie,je,ke)
    ! resulting U-profile:
    REAL(KIND=wp),     INTENT(out) :: uloc(ie,je,ke)

    ! Local parameters:
    INTEGER(KIND=iintegers) :: k, nlay

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. u_ana_linearlayers() ...'
    END IF


    nlay = nlayers_linwind

    ! Wind profile composed of piecewise linear layers (not necessarily continuous!):

    ! check consistencies:
    IF (MINVAL(zml_u) < h_linwind(1)) THEN
      WRITE (*,*) 'ERROR: problem in u_ana_linearlayers(): h_linwind(1) > surface layer!'
      CALL model_abort (my_world_id, 10103, &
           'h_linwind(1) > lowest model layer, decrease h_linwind(1)!', &
           'u_ana_linearlayers(), definition of h_linwind(1)')
      RETURN
    END IF
    IF ( ANY( h_linwind(2:nlay) - h_linwind(1:nlay-1) < 0.0_wp ) ) THEN
       WRITE (*,*) 'ERROR: problem in u_ana_linearlayers(): '// &
            'h_linwind must be monotonically increasing!'
      CALL model_abort (my_world_id, 10104, 'h_linwind must be monotonically increasing!', &
           'u_ana_linearlayers(), definition of h_linwind')
      RETURN
    END IF


    ! Set T and relhum within the first nlay-1 polytrope layers:
    DO k=1, nlay-1
      WHERE     ( zml_u(:,:,:) >= h_linwind(k) .AND. zml_u(:,:,:) < h_linwind(k+1) )
        uloc(:,:,:)  = u_linwind(k)  + ugr_linwind(k)  * ( zml_u(:,:,:) - h_linwind(k) )
      ENDWHERE
    END DO

    ! Set T and relhum for the top polytrope layer, i.e., above h_linwind(nlay):
    WHERE     ( zml_u(:,:,:) >= h_linwind(nlay) )
      uloc(:,:,:)  = u_linwind(nlay)  + ugr_linwind(nlay)  * &
           ( zml_u(:,:,:) - h_linwind(nlay) )
    ENDWHERE


    RETURN
  END SUBROUTINE u_ana_linearlayers

  !=============================================================================
  !
  ! Calculate the environment state (T, p, qv, qc=0) 
  ! for the Weisman-Klemp (1982) test case.
  ! 
  ! Inputs:  zml    = heights of model full levels at mass points
  !
  ! Outputs: temp   = potential temperature
  !          relhum = profile of relat. humid. (not limited
  !                   to qv_max_wk, because this has to be
  !                   done later by the calling routine)
  !          psurf  = surface pressure
  !          tsurf  = surface temperature
  !          rhsurf = surface rel. humid.
  !
  ! Namelist parameters for the initial profile of the pot. temp.:
  !
  !    hmin_wk          = 0.0_wp         ! [m], base height of the profile
  !    p_base_wk        = 1.0E5_wp       ! [Pa], pressure at height hmin_wk
  !    h_tropo_wk       = 12000._wp      ! [m], height of the tropopause
  !    theta_0_wk       = 300._wp        ! [K], pot. temp. at z=0
  !    theta_tropo_wk   = 343._wp        ! [K], pot. temp. at z=h_tropo_wk
  !    expo_theta_wk    = 1.25_wp
  !    expo_relhum_wk   = 1.25_wp
  !
  ! Data for the initial profile of relative humidity:
  !    rh_min_wk   = 0.25_wp        ! [%] ! rel. hum. above the tropopause level
  !    rh_max_wk   = 1.0_wp         ! [%]
  !  
  ! optionally (if set to a value >= -900):
  !    t_tropo_wk   = -999.99           ! [K], temp. at z=h_tropo_wk   
  !                                *(is 213.0 in the original literature)*
  !
  !    *NOTE* : t_tropo_wk is the tropopause temperature. If t_tropo_wk < -900, 
  !             then the tropopause temperature
  !             is determined exactly in this subroutine. This ensures
  !             that the temperature in the tropopause is constant
  !             with height. If the user sets an estimated value, then this is not
  !             necessarily the case. The original literature, e.g., used an
  !             estimated value of 213 K, which did not lead to a constant tropopause
  !             temperature.
  !
  !=============================================================================

  SUBROUTINE tqv_wk82(zml, hsurf, theta, relhum, psurf, tsurf, rhsurf)

    IMPLICIT NONE

    !----------------------------------------------------

    REAL (KIND=wp),     DIMENSION (1:ie, 1:je, 1:ke), INTENT(in)    :: zml
    REAL (KIND=wp),     DIMENSION (1:ie, 1:je),       INTENT(in)    :: hsurf
    REAL (KIND=wp),     DIMENSION (1:ie, 1:je, 1:ke), INTENT(out)   :: theta, relhum
    REAL (KIND=wp),     DIMENSION (1:ie, 1:je), INTENT(out) ::  tsurf, rhsurf, psurf

    !.. Local variables:

    !.. Max. dimension of a helper height vector, which expands
    !   throughout the entire model atmosphere and which 
    !   will be used for the interpolation of the surface pressure:
    INTEGER(KIND=iintegers), PARAMETER :: kmax_interp = 7001
    !   Height increment of the helper height vector:
    REAL(KIND=wp),     PARAMETER       :: dz_interp   = 10.0_wp
    !   actual dimension (to be determined below):
    INTEGER(KIND=iintegers)            :: ke_interp
    !.. Dimension of the helper height vector below h_tropo_wk:
    !     (to be determined below)
    INTEGER(KIND=iintegers) :: ke_interp_u
    REAL(KIND=wp)           :: dz_interp_u

    REAL (KIND=wp)     :: t_tropo, tvu, tvo, gamma, &
         zinterp(kmax_interp), tinterp(kmax_interp), relhuminterp(kmax_interp), &
         pinterp(kmax_interp), qvinterp(kmax_interp), qcinterp(kmax_interp), &
         thetainterp(kmax_interp), thetasurf(ie,je), hhlgew, qvtmp, qctmp, hsurfmax

    INTEGER(KIND=iintegers) :: i, j, k, index, ierror
    CHARACTER(len=250) :: errmsg

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. tqv_wk82() ...'
    END IF


    !.. Check setting of hmin_wk, since hmin_wk has to be smaller or equal than the lowest
    !   model level height:
    IF (ANY(hsurf < hmin_wk)) THEN
      WRITE (*,*) 'ERROR: problem in tqv_wk82(): hsurf must be larger or equal to hmin_wk!'//&
           'I.e., increase href accordingly!'
      CALL model_abort (my_world_id, 10109, 'hsurf must be larger or equal to hmin_wk!', &
           'tqv_wk82(), src_artifdata.f90')
      RETURN
    END IF

    IF (t_tropo_wk < -900.0_wp) THEN

      !==========================================================================
      !.. Determine exact tropopause temperature from hydrostatic approximation
      !==========================================================================

      !.. helper vectors for determining the tropopause temperature 
      !   (expand up to h_tropo_wk):
      ke_interp_u = MAX( CEILING((h_tropo_wk-hmin_wk)/dz_interp) + 1, 2_iintegers)
      dz_interp_u = MAX( (h_tropo_wk - hmin_wk) / (ke_interp_u-1.0_wp) , 10.0_wp)

      IF (ke_interp_u > kmax_interp) THEN
        WRITE (*,*) 'ERROR: problem in tqv_wk82(): ke_interp_u > kmax_interp!'
        CALL model_abort (my_world_id, 10104, 'ke_interp_u > kmax_interp!', &
             'tqv_wk82(), calculation of tropopause pressure')
      END IF

      DO k = 1, ke_interp_u
        zinterp(k) = hmin_wk +  dz_interp_u * (k-1.0_wp)
        thetainterp(k) = theta_0_wk + &
             (theta_tropo_wk - theta_0_wk)*((zinterp(k)-hmin_wk) / &
             (h_tropo_wk-hmin_wk))**(expo_theta_wk)
        relhuminterp(k)  = rh_max_wk - (rh_max_wk - rh_min_wk)*((zinterp(k)-hmin_wk) / &
             (h_tropo_wk-hmin_wk))**(expo_relhum_wk)
      END DO
      
      !.. Use piecewise polytrope atmosphere pressure as best estimate of pressure profile:
      CALL calc_p_hydrostat_ana(ni=1, nj=1, nk=ke_interp_u, niter=20, &
           zml=zinterp(1:ke_interp_u), psurf=(/p_base_wk/), &
           t=tinterp(1:ke_interp_u), relhum=relhuminterp(1:ke_interp_u), &
           piter=pinterp(1:ke_interp_u), qv=qvinterp(1:ke_interp_u), &
           qc=qcinterp(1:ke_interp_u), &
           zmaxqv=qv_max_wk, r_d=r_d, rvd_m_o=rvd_m_o, &
           theta=thetainterp(1:ke_interp_u), &
           ierror=ierror, errmsg=errmsg)
      
      IF (ierror /= 0) THEN
        CALL model_abort (my_world_id, ierror, TRIM(errmsg), &
             'tqv_wk82(), calc_p_hydrostat_ana(), calculation of tropopause pressure')
      END IF

      !.. Exact temperature of the tropopause (throughout the whole layer, 
      !   not just at z=h_tropo_wk:
      t_tropo  =  tinterp(ke_interp_u)

    ELSE

      !.. Take the user defined value as an estimate of temperature of the tropopause:
      t_tropo = t_tropo_wk

    END IF

    !===========================================================================
    !.. Theta and relative humidity in the atmosphere:
    !===========================================================================
    
    !.. 1) The values of Theta and relhum at the model mass points:
    DO k = 1,ke
      DO j = 1,je
        DO i = 1,ie

          IF ( zml(i,j,k) <= h_tropo_wk) THEN
            theta(i,j,k) = theta_0_wk + &
                 (theta_tropo_wk - theta_0_wk)*((zml(i,j,k)-hmin_wk)/ (h_tropo_wk-hmin_wk))**(expo_theta_wk)
            relhum(i,j,k) = rh_max_wk - (rh_max_wk - rh_min_wk)*((zml(i,j,k)-hmin_wk)/(h_tropo_wk-hmin_wk))**(expo_relhum_wk)
          ELSE
            theta(i,j,k) = theta_tropo_wk*EXP(g/(cp_d*t_tropo)*(zml(i,j,k)-h_tropo_wk))
            relhum(i,j,k)  = rh_min_wk
          ENDIF

        ENDDO
      ENDDO
    ENDDO

    !.. 2) The values of theta and relhum at the ground:
    DO j = 1,je
      DO i = 1,ie
        
        IF ( hsurf(i,j) <= h_tropo_wk) THEN
          thetasurf(i,j) = theta_0_wk + &
               (theta_tropo_wk - theta_0_wk)*((hsurf(i,j)-hmin_wk)/ (h_tropo_wk-hmin_wk))**(expo_theta_wk)
          rhsurf(i,j)  = rh_max_wk - (rh_max_wk - rh_min_wk)*((hsurf(i,j)-hmin_wk)/(h_tropo_wk-hmin_wk))**(expo_relhum_wk)
        ELSE
          thetasurf(i,j) = theta_tropo_wk*EXP(g/(cp_d*t_tropo)*(hsurf(i,j) - h_tropo_wk))
          rhsurf(i,j)    = rh_min_wk
        ENDIF
        
      ENDDO
    ENDDO

    !===========================================================================
    !.. Interpolate surface pressure:
    !===========================================================================

    !.. 1) Helper vectors for interpolating the surface pressure:
    !      (contains h_tropo_wk and expands from hmin_wk up to h_top)
!!! OLD VERSION WITH A CONSTANT PARAMETER ke_interp_u = 701 and ke_interp = 1001
!!$    DO k = 1, ke_interp_u
!!$      zinterp(k) = hmin_wk + (h_tropo_wk-hmin_wk) / (ke_interp_u-1.0_wp) * (k-1.0_wp)
!!$    END DO
!!$    DO k = ke_interp_u+1, ke_interp
!!$      zinterp(k) = h_tropo_wk + (h_top - h_tropo_wk) / &
!!$                     (ke_interp - ke_interp_u - 1.0_wp) * &
!!$                     (k         - ke_interp_u - 1.0_wp)      
!!$    END DO


    !.. 1) Helper height vector with distance dz_interp:
    !      Simpler solution as with tqv_ana_polylayers() above, because
    !      both the T- and QV-profiles are continuous:
    hsurfmax = MAXVAL(hsurf)
    ke_interp = MAX( CEILING( (hsurfmax-h_nconst(1)) / dz_interp ) + 2, 2_iintegers)
    IF (ke_interp > kmax_interp) THEN
      WRITE (*,*) 'ERROR: problem in tqv_wk82(): ke_interp > kmax_interp!'
      CALL model_abort (my_world_id, 10105, 'ke_interp > kmax_interp!', &
           'tqv_wk82(), calculation of surface pressure')
      RETURN
    END IF
    layers: DO k=1, ke_interp
      zinterp(k) = h_nconst(1) + (k-1)*dz_interp
    END DO layers


    !.. 2) Set T and relhum for the heights of the helper vector:
    DO k = 1, ke_interp
      IF ( zinterp(k) <= h_tropo_wk) THEN
        thetainterp(k) = theta_0_wk + &
             (theta_tropo_wk - theta_0_wk)*((zinterp(k)-hmin_wk)/ (h_tropo_wk-hmin_wk))**(expo_theta_wk)
        relhuminterp(k)  = rh_max_wk - (rh_max_wk - rh_min_wk)*((zinterp(k)-hmin_wk)/(h_tropo_wk-hmin_wk))**(expo_relhum_wk)
      ELSE
        thetainterp(k) = theta_tropo_wk*EXP(g/(cp_d*t_tropo)*(zinterp(k)-h_tropo_wk))
        relhuminterp(k)  = rh_min_wk
      ENDIF
    END DO

    !.. 3) Use polytrope atmosphere pressure as best estimate of pressure profile (lowest level: isothermal):
    CALL calc_p_hydrostat_ana(ni=1, nj=1, nk=ke_interp, niter=20, &
         zml=zinterp, psurf=(/p_base_wk/), t=tinterp, relhum=relhuminterp, &
         piter=pinterp, qv=qvinterp, qc=qcinterp, &
         zmaxqv=qv_max_wk, r_d=r_d, rvd_m_o=rvd_m_o, &
         ierror=ierror, errmsg=errmsg, &
         theta=thetainterp)

    IF (ierror /= 0) THEN
      CALL model_abort (my_world_id, ierror, TRIM(errmsg), &
           'tqv_wk82(), calc_p_hydrostat_ana(), calculation of surface pressure')
    END IF

    !.. 4) Interpolate surface pressure from best estimated pressure profile. The interpolation
    !   is based on analytically calculating the pressure in a polytrope atmosphere w.r.t. 
    !   virtual temperature:
    !     (uses virtual temperature to include the effects of moisture, neglects
    !      the temperature dependence of cp, as is done in the rest of COSMO)
    DO i = 1, ie
      DO j = 1, je
        
        ! Search the height interval within the helper vector: 
!!! OLD METHOD, CONSISTENT WITH THE OLD METHOD ABOVE
!!$        index = -999
!!$        DO k = 1, ke_interp-1
!!$          IF (zinterp(k) <= hsurf(i,j) .AND. zinterp(k+1) > hsurf(i,j)) THEN
!!$            index = k
!!$            EXIT
!!$          END IF
!!$        END DO
!!$        IF (index == -999) THEN
!!$          WRITE (*,*) 'ERROR: in tqv_wk82(): calculation of surface pressure!'
!!$          CALL model_abort (my_world_id, 10002, 'model topography outside wk82-height range, ' // &
!!$               'modify namelist parameters and start again!', &
!!$               'tqv_wk82(), surface pressure calculation')
!!$        END IF
        
        ! Get the index of the height interval within the helper vector:         
        index = FLOOR((hsurf(i,j)-h_nconst(1))/dz_interp) + 1

        tvu = tinterp(index) * (1.0_wp + rvd_m_o*qvinterp(index) - qcinterp(index))
        tvo = tinterp(index+1) * (1.0_wp + rvd_m_o*qvinterp(index+1) - qcinterp(index+1))
        IF (tvu /= tvo) THEN 
          ! Polytrope atmosphere w.r.t. Tv (linear Tv-profile) within the layer:
          gamma = (tvu - tvo) / (zinterp(index+1)-zinterp(index))
          psurf(i,j) = pinterp(index) * &
               (1.0_wp - gamma/tvu*(hsurf(i,j)-zinterp(index)))**(g/(r_d*gamma))
        ELSE
          ! Const. virt. temperature, barometric height formula:
          psurf(i,j) = pinterp(index) * EXP(-g*(hsurf(i,j)-zinterp(index))/(r_d*tvu))
        END IF

        !.. Calculate surface temperature from surface pressure and potential temperature:
        hhlgew = (hsurf(i,j)-zinterp(index)) / (zinterp(index+1)-zinterp(index))
        qvtmp = qvinterp(index) + (qvinterp(index+1)-qvinterp(index)) * hhlgew
        qctmp = qcinterp(index) + (qcinterp(index+1)-qcinterp(index)) * hhlgew
        tsurf(i,j) = thetasurf(i,j) * (psurf(i,j) / pt00 )** &
             (rd_moist(qvtmp,qctmp)/cp_moist_cosmo(qvtmp,qctmp,0.0_wp))

      END DO
    END DO

    RETURN
  END SUBROUTINE tqv_wk82
  
!===============================================================================
!
! Routine for setting temperature disturbances in the atmosphere at the
! beginning of the model run. It is called in lmorg.f90 before the
! first time step.
!
!===============================================================================

  SUBROUTINE set_tempdist(nx)

    IMPLICIT NONE

    ! time level (nnow or nnew)
    INTEGER(KIND=iintegers), INTENT(in) :: nx

    INTEGER(KIND=iintegers)      :: i, ii, jj, &
         anzbubs, zbub_index(ntempdist_max), &
         izerror, kzdims(24), nk2d
    REAL(KIND=wp),     ALLOCATABLE :: f_xyz_3d(:,:,:), rhtmp_3d(:,:,:)

    REAL(KIND=wp),     ALLOCATABLE ::  ttmp(:,:,:), qvtmp(:,:,:)

    CHARACTER(len=2)  ::  cii
    CHARACTER (LEN=80)         ::  yzerrmsg

    LOGICAL :: touch_moisture

    REAL(KIND=wp),     POINTER :: &
      zqv(:,:,:) => NULL()     ! QV at nx

    !*******************************************************************
    !                Disturbance on T / QV
    !*******************************************************************

    ! Return if not at first timestep (restart runs):
    IF (ntstep > 0) RETURN

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. set_tempdist() ...'
    END IF

    yzerrmsg(:) = ' '

    kzdims(:) = 0_iintegers

    zbub_index(:) = 0

    ! Retrieve the microphysics tracers
    CALL trcr_get( izerror, idt_qv, ptr_tlev = nx, ptr = zqv )
    IF ( izerror /= 0_iintegers ) THEN
      yzerrmsg = trcr_errorstr( izerror )
      CALL model_abort( my_cart_id, izerror, yzerrmsg, 'set_tempdist' )
    ENDIF

    ! superposition of one or several temperature disturbances:
    anzbubs = 0
    DO ii=1, ntempdist_max
      IF (ltempdist(ii)) THEN
        IF (ANY(reg_tempdist == TRIM(ctype_tempdist(ii)))) THEN
          anzbubs = anzbubs + 1
          zbub_index(anzbubs) = ii
        END IF
      END IF
    END DO


    IF (anzbubs > 0) THEN

      ! ..allocate local working arrays:
      ALLOCATE(f_xyz_3d(ie,je,ke))
      f_xyz_3d = 0.0_wp
      ALLOCATE(rhtmp_3d(ie,je,ke))
      rhtmp_3d = 0.0_wp
      IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
        ALLOCATE(ttmp(ie,je,ke), qvtmp(ie,je,ke))
        ttmp=0.0_wp
        qvtmp=0.0_wp
      END IF

      ! ..loop over all bubbles:
      bubbleloop: DO jj=1, anzbubs

        IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
          ttmp = t(:,:,:,nx)
          qvtmp = zqv(:,:,:)
        END IF

        ii = zbub_index(jj)

        WRITE(cii,'(i2.2)') ii

        IF (my_cart_id == 0) THEN
          WRITE (*,*) ('=',i=1,70)
          WRITE (*,*) '        Release Bubble ', ii, ' , Typ ', &
               TRIM(ctype_tempdist(ii)),' , time-step ', ntstep, ' ...'
          WRITE (*,*) ('=',i=1,70)
        END IF

        IF (ladd_bubblenoise_t(ii)) THEN

          ! Certain disturbance types may be modulated by some white noise
          ! (if ladd_bubblenoise_t(ii) = .true.). 
          ! For this noise, we need a storage field, and for technical
          ! reasons (ease of implementation when integrating all
          ! possible combinations of namelist switches), this is a global
          ! persistent array which needs to be allocated once:
          IF (.NOT. ALLOCATED(bub_tnoisebuffer)) THEN
            ! DO NOT TOUCH THE ALLOCATION INDICES, THEY HAVE TO BE THE SAME 
            ! FOR ALL DISTURBANCE ROUTINES!
            ALLOCATE(bub_tnoisebuffer(ie,je,MAX(ke,ke_soil),1:numbubs_noise))
            !            WRITE (*,*) 'bub_tnoisebuffer allokiert!'
            bub_tnoisebuffer = 0.0_wp
          END IF

          ! To generate a reproducible noise on different numbers of
          ! processors, only PE 0 generates the noise on a global field
          ! and distributes it to the other PEs. This is done in SR gen_bubnoise.
          ! The noise is stored in the field 
          ! bub_tnoisebuffer(:,:,kstart:kend,bub_index_noise(ii)).

          SELECT CASE (TRIM(ctype_tempdist(ii)))

          CASE ('cos', 'mcnider')

            ! We need a 3D noise field:

            IF (bub_index_noise(ii) > 0) THEN
              CALL gen_bubnoise(fnoise=bub_tnoisebuffer(:,:,1:ke,bub_index_noise(ii)), &
                   kdim=ke, seed=505_iintegers + ii)
            ELSE
              WRITE (*,*) 'ERROR: in set_tempdist_t(), '// &
                   'src_artifdata.f90: bub_index_noise(ii) == 0'
              CALL model_abort (my_world_id, 10002+izerror, 'bub_index_noise('//cii// &
                   ') == 0 occured, but is wrong here!', &
                   'set correct bub_index_noise('//cii//') in input_artifctl()!')
            END IF


          CASE ('hotspot')

            ! We need a 2D noise field:

            IF (bub_index_noise(ii) > 0) THEN
              CALL gen_bubnoise(fnoise=bub_tnoisebuffer(:,:,1,bub_index_noise(ii)), &
                   kdim=1, seed=505_iintegers + ii)
            ELSE
              WRITE (*,*) 'ERROR: in set_tempdist_t(), '// &
                   'src_artifdata.f90: bub_index_noise(ii) == 0'
              CALL model_abort (my_world_id, 10002+izerror, 'bub_index_noise('//cii// &
                   ') == 0 occured, but is wrong here!', &
                   'set correct bub_index_noise('//cii//') in input_artifctl()!')
            END IF

          END SELECT

        END IF   ! IF ( ladd_bubblenoise_t(ii) ) ...

        ! Store RH before bubble if (probably) needed later:
        IF (lbub_rhconst(ii)) THEN
          ! storage of relat. humid. before heating:
          rhtmp_3d = rh_Tpqv_3d(p0(:,:,:)+pp(:,:,:,nx), t(:,:,:,nx), &
                                zqv(:,:,:), zqv(:,:,:)*0.0_wp, ie, je, ke)
        END IF

        ! Per default disable QV-change; this might be switched 
        ! on for certain bubbles below,
        ! which then allows to keep the relhum constant during heating,
        ! if lbub_rhconst(ii)=.true.
        touch_moisture = .FALSE.

        ! Indicate if the bubble is 3D, then nk2d = -1, or if the bubble
        ! is 2D (z=const.), then nk2d = height level number of bubble.
        ! It defaults to -1 and may be changed below appropriately for certain
        ! bubble type(s).
        nk2d = -1

        SELECT CASE (TRIM(ctype_tempdist(ii)))

!!!=============================================================================

        CASE ('mcnider')

          touch_moisture = .TRUE.
          nk2d = -1

          CALL f_xyz_mcnider(ii, f_xyz_3D)

          IF (ladd_bubblenoise_t(ii)) THEN
            t(:,:,:,nx) = t(:,:,:,nx) + f_xyz_3D(:,:,:) * &
                 (1.0_wp + bub_dT_bubblenoise(ii) * &
                 bub_tnoisebuffer(:,:,1:ke,bub_index_noise(ii)))
          ELSE
            t(:,:,:,nx) = t(:,:,:,nx) + f_xyz_3D(:,:,:)
          END IF

!!!=============================================================================

        CASE ('hotspot')

          ! gaussian temperature disturbance in the lowest model layer,
          ! general 2D ellipse with the possibility to rotate horizontally.

          touch_moisture = .TRUE.
          nk2d = ke

          CALL f_xy_hotspot(ii, f_xyz_3d(:,:,1))

          IF (ladd_bubblenoise_t(ii)) THEN
            t(:,:,ke,nx) = t(:,:,ke,nx) + bub_dT(ii) * f_xyz_3D(:,:,1) * &
                 (1.0_wp + bub_dT_bubblenoise(ii) * &
                 bub_tnoisebuffer(:,:,1,bub_index_noise(ii)))
          ELSE
            t(:,:,ke,nx) = t(:,:,ke,nx) + bub_dT(ii) * f_xyz_3D(:,:,1)
          END IF


!!!=============================================================================

        CASE ('SK94')

          ! warm, dry bubble for the test case of Skamarock, Klemp (1994) MWR
          touch_moisture = .FALSE.
          nk2d = -1

          IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
            ttmp = theta_ini(:,:,:) * ((p0(:,:,:)+pp(:,:,:,nnew))/pt00)**rdocp
          END IF

          CALL f_xyz_SK94(ii, f_xyz_3d)

          ! calc. the absolute temperature:
          t(:,:,:,nnew) = ( theta_ini(:,:,:) + bub_dT(ii) * f_xyz_3d(:,:,:) ) * &
               ((p0(:,:,:)+pp(:,:,:,nnew))/pt00)**rdocp

!!!=============================================================================

        CASE ('cos')

          touch_moisture = .TRUE.
          nk2d = -1

          CALL f_xyz_cos(ii, f_xyz_3d)

          IF (ladd_bubblenoise_t(ii)) THEN
            t(:,:,:,nx) = t(:,:,:,nx) + bub_dT(ii) * f_xyz_3d(:,:,:) * &
                 (1.0_wp + bub_dT_bubblenoise(ii) * &
                 bub_tnoisebuffer(:,:,1:ke,bub_index_noise(ii)))
          ELSE
            t(:,:,:,nx) = t(:,:,:,nx) + bub_dT(ii) * f_xyz_3d(:,:,:)
          END IF

!!!=============================================================================

        CASE ('squall3D')

          ! Dry disturbance:
          touch_moisture = .FALSE.
          nk2d = -1

          CALL  f_xyz_squall3d(ii, f_xyz_3d)

          t(:,:,:,nx) = t(:,:,:,nx) + bub_dT(ii)* f_xyz_3d(:,:,:)

          WRITE (*,*) '=== SQUALL3D-Testcase: dTmax = ', &
               MAXVAL( bub_dT(ii)* f_xyz_3d(:,:,:) )

!!!=============================================================================

        END SELECT

        IF (touch_moisture) THEN
          IF (lbub_rhconst(ii)) THEN
            ! conservation of relat. humid. during heating:
            IF (nk2d > -1) THEN
              zqv(:,:,nk2d) = RESHAPE( &
                   qv_Tprelhum_3d(p0(:,:,nk2d)+pp(:,:,nk2d,nx), t(:,:,nk2d,nx), &
                   rhtmp_3d(:,:,nk2d), rhtmp_3d(:,:,nk2d)*0.0_wp, ie, je, 1_iintegers), &
                   (/ie, je/) )
            ELSE
              zqv(:,:,:) = qv_Tprelhum_3d(p0(:,:,:)+pp(:,:,:,nx), t(:,:,:,nx), &
                                            rhtmp_3d, rhtmp_3d*0.0_wp, ie, je, ke)
            END IF
          END IF
        END IF

!!!=============================================================================

        IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
          CALL output3d_ascii(t(:,:,:,nx)-ttmp(:,:,:), ke, &
               't-'//TRIM(ctype_tempdist(ii))//'-'//cii, &
               'T-disturbance atmosphere', 'K', 'eta')
          CALL output3d_ascii(zqv(:,:,:)-qvtmp(:,:,:), ke, &
               'qv-'//TRIM(ctype_tempdist(ii))//'-'//cii, &
               'Qv-disturbance atmosphere', '-', 'eta')
        END IF

      END DO bubbleloop

      ! ..clean up local memory:
      DEALLOCATE(f_xyz_3d, rhtmp_3d)

      ! exchange temperature and moisture (for safety and for periodic BCs):
      kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries &
           (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
           kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
           lperi_x, lperi_y, l2dim, &
           20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
           t(:,:,:,nx), zqv(:,:,:) )

      ! Clean up memory after all bubbles have been released:
      IF (.NOT.ANY(ltempdist)) THEN
        IF (ALLOCATED(bub_tnoisebuffer)) THEN
          DEALLOCATE(bub_tnoisebuffer)
          !      WRITE (*,*) 'bub_tnoisebuffer deallokiert!'
        END IF
      END IF

      IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
        DEALLOCATE(ttmp,qvtmp)
      END IF

    END IF   ! IF (anzbubs > 0)

  END SUBROUTINE set_tempdist

  !======================================================================================
  !
  ! Temperature disturbance in the soil (initial condition), makes sense if lsoil=.true.
  ! (is called in lmorg.f90)
  !
  !======================================================================================

  SUBROUTINE set_tempdist_tso(nx)

    IMPLICIT NONE

    ! time level (nnow or nnew)
    INTEGER(KIND=iintegers), INTENT(in) :: nx

    INTEGER(KIND=iintegers)      :: i, ii, jj, &
         anzbubs, zbub_index(ntempdist_max), &
         izerror, kzdims(24)
    REAL(KIND=wp),     ALLOCATABLE :: f_xyz_3d(:,:,:)

    REAL(KIND=wp),     ALLOCATABLE ::  ttmp(:,:,:)
    CHARACTER(len=2)  ::  cii

    CHARACTER (LEN=80)         ::  yzerrmsg

    !*******************************************************************
    !                Disturbance on T in the soil
    !*******************************************************************

    ! Return if not at first timestep (restart runs):
    IF (ntstep > 0) RETURN

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. set_tempdist_tso() ...'
    END IF

    yzerrmsg(:) = ' '

    IF (lsoil) THEN

      kzdims(:) = 0_iintegers

      zbub_index(:) = 0_iintegers

      ! superposition of one or several temperature disturbances in the soil:

      anzbubs = 0
      DO ii=1, ntempdist_max
        IF (ltempdist(ii)) THEN
          IF (ANY(reg_tempdist_tso == TRIM(ctype_tempdist(ii)))) THEN
            anzbubs = anzbubs + 1
            zbub_index(anzbubs) = ii
          END IF
        END IF
      END DO


      IF (anzbubs > 0) THEN

        ! ..allocate local working array:
        ALLOCATE(f_xyz_3d(ie,je,0:ke_soil+1))
        f_xyz_3d = 0.0_wp
        IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
          ALLOCATE(ttmp(ie,je,0:ke_soil+1))
          ttmp=0.0_wp
        END IF

        ! ..loop over all bubbles:
        bubbleloop: DO jj = 1, anzbubs

          IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
            IF (lmulti_layer) THEN
              ttmp = t_so(:,:,:,nx)
            ELSE
              ttmp(:,:,1) = t_s(:,:,nx)
              ttmp(:,:,2) = t_m(:,:,nx)
              ttmp(:,:,3) = t_cl(:,:)
            END IF
          END IF

          ii = zbub_index(jj)

          WRITE(cii,'(i2.2)') ii

          IF (my_cart_id == 0) THEN
            WRITE (*,*) ('=',i=1,70)
            WRITE (*,*) '        Set soil T-disturbance ', ii, ' , Typ ', &
                 TRIM(ctype_tempdist(ii)),' , time-step ', ntstep, ' ...'
            WRITE (*,*) ('=',i=1,70)
          END IF

          IF (ladd_bubblenoise_t(ii)) THEN

            ! Certain disturbance types may be modulated by some white noise
            ! (if ladd_bubblenoise_t(ii) = .true.). 
            ! For this noise, we need a storage field, and for technical
            ! reasons (ease of implementation when integrating all
            ! possible combinations of namelist switches), this is a global
            ! persistent array which needs to be allocated once:
            IF (.NOT. ALLOCATED(bub_tnoisebuffer)) THEN
              ! DO NOT TOUCH THE ALLOCATION INDICES, THEY HAVE TO BE THE SAME 
              ! FOR ALL DISTURBANCE ROUTINES!
              ALLOCATE(bub_tnoisebuffer(ie,je,MAX(ke,ke_soil),1:numbubs_noise))
              !            WRITE (*,*) 'bub_tnoisebuffer allokiert!'
              bub_tnoisebuffer = 0.0_wp
            END IF

            ! To generate a reproducible noise on different numbers of
            ! processors, only PE 0 generates the noise on a global field
            ! and distributes it to the other PEs. This is done in SR gen_bubnoise.
            ! The noise is stored in the global field 
            ! bub_tnoisebuffer(:,:,kstart:kend,bub_index_noise(ii)).

            SELECT CASE (TRIM(ctype_tempdist(ii)))

            CASE ('cos-soil')

              ! We need a 3D noise field:

              IF (bub_index_noise(ii) > 0) THEN
                CALL gen_bubnoise(fnoise=bub_tnoisebuffer(:,:,1:ke_soil,bub_index_noise(ii)), &
                     kdim=ke_soil, seed=505_iintegers + ii)
              ELSE
                WRITE (*,*) 'ERROR: in set_tempdist_tso(), '// &
                     'src_artifdata.f90: bub_index_noise(ii) == 0'
                izerror = 1
                CALL model_abort (my_world_id, 10002+izerror, 'bub_index_noise('//cii// &
                     ') == 0 occured, but is wrong here!', &
                     'set correct bub_index_noise('//cii//') in input_artifctl()!')
              END IF

            CASE ('hotspot-soil')

              ! We need a 2D noise field:

              IF (bub_index_noise(ii) > 0) THEN
                CALL gen_bubnoise(fnoise=bub_tnoisebuffer(:,:,1,bub_index_noise(ii)), &
                     kdim=1, seed=505_iintegers + ii)
              ELSE
                WRITE (*,*) 'ERROR: in set_tempdist_tso(), '// &
                     'src_artifdata.f90: bub_index_noise(ii) == 0'
                CALL model_abort (my_world_id, 10002+izerror, 'bub_index_noise('//cii &
                     //') == 0 occured, but is wrong here!', &
                     'set correct bub_index_noise('//cii//') in input_artifctl()!')
              END IF

            END SELECT

          END IF   ! IF ( ladd_bubblenoise_t(ii) ) ...


          SELECT CASE (ctype_tempdist(ii))

!!!=============================================================================

          CASE ('hotspot-soil')

            ! gaussian temperature disturbance in the uppermost soil level,
            ! general 2D ellipse with the possibility to rotate horizontally.

            CALL f_xy_hotspot(ii, f_xyz_3d(:,:,1))

            IF (lmulti_layer) THEN
              IF (ladd_bubblenoise_t(ii)) THEN
                t_so(:,:,0,nx) = t_so(:,:,0,nx) + bub_dT(ii) * f_xyz_3D(:,:,1) * &
                     (1.0_wp + bub_dT_bubblenoise(ii) * &
                     bub_tnoisebuffer(:,:,1,bub_index_noise(ii)))
              ELSE
                t_so(:,:,0,nx) = t_so(:,:,0,nx) + bub_dT(ii) * f_xyz_3D(:,:,1)
              END IF
              t_so(:,:,1,nx) = t_so(:,:,0,nx)
            ELSE
              IF (ladd_bubblenoise_t(ii)) THEN
                t_s(:,:,nx) = t_s(:,:,nx) + bub_dT(ii) * f_xyz_3D(:,:,1) * &
                     (1.0_wp + bub_dT_bubblenoise(ii) * &
                     bub_tnoisebuffer(:,:,1,bub_index_noise(ii)))
              ELSE
                t_s(:,:,nx) = t_s(:,:,nx) + bub_dT(ii) * f_xyz_3D(:,:,1)
              END IF
            END IF

!!!=============================================================================

          CASE ('cos-soil')

            CALL f_xyz_cos2d_soil(ii, f_xyz_3d(:,:,1:ke_soil))

            IF (lmulti_layer) THEN
              IF (ladd_bubblenoise_t(ii)) THEN
                t_so(:,:,1:ke_soil,nx) = t_so(:,:,1:ke_soil,nx) + &
                     bub_dT(ii) * f_xyz_3d(:,:,1:ke_soil) * &
                     (1.0_wp + bub_dT_bubblenoise(ii) * &
                     bub_tnoisebuffer(:,:,1:ke_soil,bub_index_noise(ii)))
              ELSE
                t_so(:,:,1:ke_soil,nx) = t_so(:,:,1:ke_soil,nx) + &
                     bub_dT(ii) * f_xyz_3d(:,:,1:ke_soil)
              END IF
              t_so(:,:,0,nx)         = t_so(:,:,1,nx)
              t_so(:,:,ke_soil+1,nx) = t_so(:,:,ke_soil,nx)
              IF (lalloc_t_cl) THEN
!!! UB: necessary????
                IF (ladd_bubblenoise_t(ii)) THEN
                  t_cl(:,:)   = t_cl(:,:)   + bub_dT(ii) * f_xyz_3D(:,:,ke_soil) * &
                       (1.0_wp + bub_dT_bubblenoise(ii) * &
                       bub_tnoisebuffer(:,:,ke_soil,bub_index_noise(ii)))
                ELSE
                  t_cl(:,:)   = t_cl(:,:)   + bub_dT(ii) * f_xyz_3D(:,:,ke_soil)
                END IF
              END IF
            ELSE
              IF (ladd_bubblenoise_t(ii)) THEN
                t_s(:,:,nx) = t_s(:,:,nx) + bub_dT(ii) * f_xyz_3D(:,:,1) * &
                     (1.0_wp + bub_dT_bubblenoise(ii) * &
                     bub_tnoisebuffer(:,:,1,bub_index_noise(ii)))
                t_m(:,:,nx) = t_m(:,:,nx) + bub_dT(ii) * f_xyz_3D(:,:,2) * &
                     (1.0_wp + bub_dT_bubblenoise(ii) * &
                     bub_tnoisebuffer(:,:,2,bub_index_noise(ii)))
                t_cl(:,:)   = t_cl(:,:)   + bub_dT(ii) * f_xyz_3D(:,:,3) * &
                     (1.0_wp + bub_dT_bubblenoise(ii) * &
                     bub_tnoisebuffer(:,:,3,bub_index_noise(ii)))
              ELSE
                t_s(:,:,nx) = t_s(:,:,nx) + bub_dT(ii) * f_xyz_3D(:,:,1)
                t_m(:,:,nx) = t_m(:,:,nx) + bub_dT(ii) * f_xyz_3D(:,:,2)
                t_cl(:,:)   = t_cl(:,:)   + bub_dT(ii) * f_xyz_3D(:,:,3)
              END IF
            END IF

!!!=============================================================================

          END SELECT

          IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
            IF (lmulti_layer) THEN
              f_xyz_3d(:,:,:) = t_so(:,:,:,nx) - ttmp
              CALL output3d_ascii(f_xyz_3d(:,:,1:ke_soil+1), ke_soil+1, &
                   'tso-'//TRIM(ctype_tempdist(ii))//'-'//cii, &
                   'T-disturbance soil', 'K', 'eta')
            ELSE
              f_xyz_3d(:,:,1) = t_s(:,:,nx) - ttmp(:,:,1)
              f_xyz_3d(:,:,2) = t_m(:,:,nx) - ttmp(:,:,2)
              f_xyz_3d(:,:,3) = t_cl(:,:)   - ttmp(:,:,3)
              CALL output3d_ascii(f_xyz_3d(:,:,1:3), 3, &
                   'tso-'//TRIM(ctype_tempdist(ii))//'-'//cii, &
                   'T-disturbance soil', 'K', 'eta')
            END IF
          END IF

        END DO bubbleloop

        ! ..clean up local memory:
        DEALLOCATE(f_xyz_3d)
        IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
          DEALLOCATE(ttmp)
        END IF

        ! exchange temperature (for safety and for periodic BCs):
        IF (lmulti_layer) THEN
          IF (.NOT. lalloc_t_cl) THEN
            kzdims(1:24)=(/ke_soil+2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
            CALL exchg_boundaries &
                 (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
                 kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
                 lperi_x, lperi_y, l2dim, &
                 20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
                 t_so(:,:,:,nx) )
          ELSE
            kzdims(1:24)=(/ke_soil+2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
            CALL exchg_boundaries &
                 (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
                 kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
                 lperi_x, lperi_y, l2dim, &
                 20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
                 t_so(:,:,:,nx), t_cl(:,:) )
          END IF
        ELSE
          kzdims(1:24)=(/1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries &
               (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
               kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
               lperi_x, lperi_y, l2dim, &
               20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
               t_s(:,:,nx), t_m(:,:,nx), t_cl(:,:) )
        END IF

        ! Clean up global memory after all bubbles have been released:
        IF (.NOT.ANY(ltempdist)) THEN
          IF (ALLOCATED(bub_tnoisebuffer)) THEN
            DEALLOCATE(bub_tnoisebuffer)
            !      WRITE (*,*) 'bub_tnoisebuffer deallokiert!'
          END IF
        END IF

      END IF   ! IF (anzbubs > 0)

    END IF  ! IF (lsoil)

  END SUBROUTINE set_tempdist_tso


  !======================================================================================
  !
  ! ... Disturbance in the bottom boundary condition for T_s.
  !     Is called at the beginning of every timestep in lmorg.f90 in 
  !     SUBROUTINE initialize_loop().
  !     The disturbance is active over a certain number of timesteps (namelist parameter
  !     "bub_timespan") after time "htempdist", so it is an on/off-behaviour.
  !     However, in principle also more complicated time dependencies are possible,
  !     using "ntstep*dt" as a time measure.
  !
  !     NOTE: makes only sense if soil model is turned off (lsoil=.FALSE.) or for water points
  !           (fr_land < 0.5).
  !
  !======================================================================================

  SUBROUTINE set_tempdist_bbc_ts()

    IMPLICIT NONE

    INTEGER(KIND=iintegers)      :: i, ii, jj, &
         anzbubs, zbub_index(ntempdist_max), &
         izerror, kzdims(24)
    REAL(KIND=wp)     :: f_xyz_2d(ie,je)

    CHARACTER(len=2)  ::  cii

    CHARACTER (LEN=80)         ::  yzerrmsg

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. set_tempdist_bbc_ts() ...'
    END IF

    yzerrmsg(:) = ' '

    IF (.NOT.lsoil) THEN

      kzdims(:) = 0_iintegers

      zbub_index(:) = 0_iintegers

      anzbubs = 0
      DO ii=1, ntempdist_max
        IF (ltempdist(ii) .AND. bub_timespan(ii) > 0) THEN
          IF (ANY(reg_tempdist_bbc_ts == TRIM(ctype_tempdist(ii)))) THEN
            anzbubs = anzbubs + 1
            zbub_index(anzbubs) = ii
          END IF
        END IF
      END DO

      IF (anzbubs > 0) THEN

        ! As a basis for the following time series of t_s, initialize t_s with its
        ! (time constant) boundary values, which itself contain the intitial state.
        ! Later in this routine, disturbances are added in every time step within
        ! the disturbance timespan. At the moment, only time constant disturbances
        ! are implemented, but this may changed by the user. Just implement
        ! your own (time dependent) formula and use, e.g., ntstep*dt as a time
        ! measure.
        ! This initialization also ensures, that if the disturbance 
        ! time has elapsed, t_s attains back its initial state.
        t_s(:,:,nnow) = t_s_bd(:,:,1)
        t_s(:,:,nnew) = t_s_bd(:,:,1)

        t_snow(:,:,nnow) = t_s_bd(:,:,1)
        t_snow(:,:,nnew) = t_s_bd(:,:,1)

        f_xyz_2d = 0.0_wp

        bubbleloop: DO jj = 1, anzbubs

          ii = zbub_index(jj)

          WRITE(cii,'(i2.2)') ii

          IF ( (ntstep_bubble(ii) == ntstep .OR. &
               (bub_zeitzaehler(ii) > 0 .AND. bub_zeitzaehler(ii) < bub_timespan(ii)))) THEN

            bub_zeitzaehler(ii) = bub_zeitzaehler(ii) + 1

            IF (my_cart_id == 0) THEN
              WRITE (*,*) ('=',i=1,70)
              WRITE (*,*) '        Set T_S disturbance ', ii, ' , Typ ', &
                   TRIM(ctype_tempdist(ii)),' , disturbance-time-step ', bub_zeitzaehler(ii), ' ...'
              WRITE (*,*) ('=',i=1,70)
            END IF

            IF (bub_zeitzaehler(ii) == 1) THEN

              IF (ladd_bubblenoise_t(ii)) THEN

                ! Certain disturbance types may be modulated by some white noise
                ! (if ladd_bubblenoise_t(ii) = .true.). 
                ! For this noise, we need a storage field, and for technical
                ! reasons (ease of implementation when integrating all
                ! possible combinations of namelist switches), this is a global
                ! persistent array which needs to be allocated once:
                IF (.NOT. ALLOCATED(bub_tnoisebuffer)) THEN
                  ! DO NOT TOUCH THE ALLOCATION INDICES, THEY HAVE TO BE THE SAME 
                  ! FOR ALL DISTURBANCE ROUTINES!
                  ALLOCATE(bub_tnoisebuffer(ie,je,MAX(ke,ke_soil),1:numbubs_noise))
                  !            WRITE (*,*) 'bub_tnoisebuffer allokiert!'
                  bub_tnoisebuffer = 0.0_wp
                END IF

                ! To generate a reproducible noise on different numbers of
                ! processors, only PE 0 generates the noise on a global field
                ! and distributes it to the other PEs. This is done in SR gen_bubnoise.
                ! The noise is stored in the global field 
                ! bub_tnoisebuffer(:,:,kstart:kend,bub_index_noise(ii)).


                SELECT CASE (TRIM(ctype_tempdist(ii)))

                CASE ('hotspot-sfc')  ! (formerly tsurface-gauss)
                  
                  ! We need a 2D noise field:

                  IF (bub_index_noise(ii) > 0) THEN
                    CALL gen_bubnoise(fnoise=bub_tnoisebuffer(:,:,1,bub_index_noise(ii)), &
                         kdim=1, seed=505_iintegers + ii)
                  ELSE
                    WRITE (*,*) 'ERROR: in set_tempdist_bbc_ts(), '// &
                         'src_artifdata.f90: bub_index_noise(ii) == 0'
                    izerror = 1
                    CALL model_abort (my_world_id, 10002+izerror, 'bub_index_noise('//cii// &
                         ') == 0 occured, but is wrong here!', &
                         'set correct bub_index_noise('//cii//') in input_artifctl()!')
                  END IF

                END SELECT

              END IF   ! IF ( ladd_bubblenoise_t(ii) ) ...

            END IF   ! IF (bub_zeitzaehler(ii) == 1) ...

            SELECT CASE (TRIM(ctype_tempdist(ii)))

!!!=============================================================================

            CASE ('hotspot-sfc')  ! (formerly tsurface-gauss)

              ! gaussian temperature disturbance at the soil-atmosphere interface,
              ! general 2D ellipse with the possibility to rotate horizontally.

              CALL f_xy_hotspot(ii, f_xyz_2d(:,:))

              IF (ladd_bubblenoise_t(ii)) THEN
                t_s(:,:,nnow) = t_s(:,:,nnow) + bub_dT(ii) * f_xyz_2D(:,:) * &
                     (1.0_wp + bub_dT_bubblenoise(ii) * &
                     bub_tnoisebuffer(:,:,1,bub_index_noise(ii)))
              ELSE
                t_s(:,:,nnow) = t_s(:,:,nnow) + bub_dT(ii) * f_xyz_2D(:,:)
              END IF
              t_s(:,:,nnew) = t_s(:,:,nnow)

!!!=============================================================================

            END SELECT

            IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
              CALL output3d_ascii(t_s(:,:,nnow), 1, &
                   'ts-'//TRIM(ctype_tempdist(ii))//'-'//cii, &
                   'T_S-disturbance earth surface', 'K', 'eta')
            END IF

          ELSE IF (bub_zeitzaehler(ii) == bub_timespan(ii)) THEN
            
            ! ... If all disturbances are done, increment time counter by 1 
            !     so that it gets larger than bub_timespan:
            bub_zeitzaehler(ii) = bub_zeitzaehler(ii) + 1
            ! ... and set the disturbance master switch to .FALSE. With this,
            !     it is possible to check later if all disturbances are
            !     finished.
            ltempdist(ii) = .FALSE.

          END IF  ! IF ( (ntstep_bubble(ii) == ntstep .OR. ...

        END DO bubbleloop

        !.. simplified diagnosis of interface temperature, since
        !   the soil model (and hence the snow model) is turned off:
        t_g = t_s
        t_snow = t_s


        ! exchange temperature (for safety and for periodic BCs):
        kzdims(1:24)=(/1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries &
             (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
             kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
             lperi_x, lperi_y, l2dim, &
             20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
             t_s(:,:,nnow), t_snow(:,:,nnow), t_g(:,:,nnow),  &
             t_s(:,:,nnew), t_snow(:,:,nnew), t_g(:,:,nnew) )

        ! Clean up global memory after all bubbles have been released:
        IF (.NOT.ANY(ltempdist)) THEN
          IF (ALLOCATED(bub_tnoisebuffer)) THEN
            DEALLOCATE(bub_tnoisebuffer)
            !      WRITE (*,*) 'bub_tbuffer deallokiert!'
          END IF
        END IF


      END IF  ! IF (anzbubs > 0)

    END IF  ! IF (.NOT.lsoil)

  END SUBROUTINE set_tempdist_bbc_ts


  !======================================================================================
  !
  ! ... Heating rate disturbance in the atmosphere within any time interval during
  !     the simulation starting from "htempdist" for "bub_timespan" timesteps.
  !     Is called at the beginning of every timestep in lmorg.f90 in 
  !     SUBROUTINE initialize_loop().
  !     The disturbance is imposed on the temperature tendency variable
  !     ttens. If desired (namelist switch lbub_rhconst=.true.) then also qvtens 
  !     is disturbed in a way that the relative humidity tendency is 0. 
  !
  !     NOTE: REQUIRES THAT RHO IS COMPUTED RIGHT BEFORE CALLING THIS ROUTINE!
  !
  !======================================================================================

  SUBROUTINE artif_heatrate_dist(nx)

    IMPLICIT NONE

    ! time level (nnow or nnew) from which the actual values of QV and T are taken
    ! to determine QVTENS from TTENS in case of conservation of rel. hum. during heating:
    INTEGER(KIND=iintegers), INTENT(in) :: nx

    INTEGER(KIND=iintegers)      :: i, j, k, ii, jj, &
         anzbubs, zbub_index(ntempdist_max), &
         izerror, kzdims(24), nk2d
    REAL(KIND=wp),     ALLOCATABLE :: f_xyz_3d(:,:,:), rhtmp_3d(:,:,:)

    REAL(KIND=wp),     ALLOCATABLE ::  ttmp(:,:,:), qvtmp(:,:,:)
    CHARACTER(len=2)  ::  cii

    CHARACTER (LEN=80)         ::  yzerrmsg

    LOGICAL :: touch_moisture, alloc_ttmp

    REAL(KIND=wp),     POINTER :: &
      zqv_tens(:,:,:) => NULL(), & ! QV tendency
      zqv(:,:,:)      => NULL()    ! QV at nx

    !*******************************************************************
    !                Disturbance on ttens / qvtens
    !*******************************************************************

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. artif_heatrate_dist() ...'
    END IF

    yzerrmsg(:) = ' '
    izerror     = 0_iintegers

    kzdims(:) = 0_iintegers

    zbub_index(:) = 0_iintegers

    ! Retrieve the microphysics tracers
    CALL trcr_get( izerror, idt_qv, ptr_tlev = nx, ptr = zqv )
    IF ( izerror /= 0_iintegers ) THEN
      yzerrmsg = trcr_errorstr( izerror )
      CALL model_abort( my_cart_id, izerror, yzerrmsg, 'set_tempdist' )
    ENDIF
    CALL trcr_get( izerror, idt_qv, ptr_tens = zqv_tens )
    IF ( izerror /= 0_iintegers ) THEN
      yzerrmsg = trcr_errorstr( izerror )
      CALL model_abort( my_cart_id, izerror, yzerrmsg, 'set_tempdist' )
    ENDIF

    ! superposition of one or several temperature disturbances:

    alloc_ttmp = (ldebug_artif .AND. idbg_artif_level > 3)

    anzbubs = 0
    DO ii=1, ntempdist_max
      IF (ltempdist(ii) .AND. bub_timespan(ii) > 0) THEN
        IF (ANY(reg_heatrate_dist == TRIM(ctype_tempdist(ii)))) THEN
          anzbubs = anzbubs + 1
          zbub_index(anzbubs) = ii
          ! Check if the local storage fields for ttmp and qvtmp are needed
          ! besides the debug mode. This is the case if relhum
          ! is to be conserved via specification of a proper qvtens
          ! consistent with ttens:
          IF (lbub_rhconst(ii)) THEN
            alloc_ttmp = .TRUE.
          END IF
        END IF
      END IF
    END DO

    IF (anzbubs > 0) THEN

      ! ..allocate local working arrays:
      ALLOCATE(f_xyz_3d(ie,je,ke))
      f_xyz_3d = 0.0_wp
      ALLOCATE(rhtmp_3d(ie,je,ke))
      rhtmp_3d = 0.0_wp
      IF (alloc_ttmp) THEN
        ALLOCATE(ttmp(ie,je,ke), qvtmp(ie,je,ke))
        ttmp=0.0_wp
        qvtmp=0.0_wp
      END IF

      ! ..loop over all disturbances:
      bubbleloop: DO jj=1, anzbubs

        ii = zbub_index(jj)

        WRITE(cii,'(i2.2)') ii

        IF (lbub_rhconst(ii) .OR. (ldebug_artif .AND. idbg_artif_level > 3)) THEN
          ttmp = ttens
          qvtmp = zqv_tens
        END IF

        IF ( (ntstep_bubble(ii) == ntstep .OR. &
             (bub_zeitzaehler(ii) > 0 .AND. bub_zeitzaehler(ii) < bub_timespan(ii)))) THEN

          bub_zeitzaehler(ii) = bub_zeitzaehler(ii) + 1

          IF (my_cart_id == 0) THEN
            WRITE (*,*) ('=',i=1,70)
            WRITE (*,*) '        Heating rate disturbance ', ii, ' , Typ ', &
                 TRIM(ctype_tempdist(ii)),' , disturbance-time-step ', bub_zeitzaehler(ii), ' ...'
            WRITE (*,*) ('=',i=1,70)
          END IF

          IF (bub_zeitzaehler(ii) == 1) THEN

            IF (ladd_bubblenoise_t(ii)) THEN

              ! Certain disturbance types may be modulated by some white noise
              ! (if ladd_bubblenoise_t(ii) = .true.). 
              ! For this noise, we need a storage field, and for technical
              ! reasons (ease of implementation when integrating all
              ! possible combinations of namelist switches), this is a global
              ! persistent array which needs to be allocated once:
              IF (.NOT. ALLOCATED(bub_tnoisebuffer)) THEN
                ! DO NOT TOUCH THE ALLOCATION INDICES, THEY HAVE TO BE THE SAME 
                ! FOR ALL DISTURBANCE ROUTINES!
                ALLOCATE(bub_tnoisebuffer(ie,je,MAX(ke,ke_soil),1:numbubs_noise))
                !            WRITE (*,*) 'bub_tnoisebuffer allokiert!'
                bub_tnoisebuffer = 0.0_wp
              END IF

              ! To generate a reproducible noise on different numbers of
              ! processors, only PE 0 generates the noise on a global field
              ! and distributes it to the other PEs. This is done in SR gen_bubnoise.
              ! The noise is stored in the global field 
              ! bub_tnoisebuffer(:,:,kstart:kend,bub_index_noise(ii)).


              SELECT CASE (TRIM(ctype_tempdist(ii)))

              CASE ('cos-hrd', 'mcnider-hrd')

                ! We need a 3D noise field:

                IF (bub_index_noise(ii) > 0) THEN
                  CALL gen_bubnoise(fnoise=bub_tnoisebuffer(:,:,1:ke,bub_index_noise(ii)), &
                       kdim=ke, seed=505_iintegers + ii)
                ELSE
                  WRITE (*,*) 'ERROR: in artif_heatrate_dist(), '// &
                       'src_artifdata.f90: bub_index_noise(ii) == 0'
                  CALL model_abort (my_world_id, 10002+izerror, 'bub_index_noise('//cii// &
                       ') == 0 occured, but is wrong here!', &
                       'set correct bub_index_noise('//cii//') in input_artifctl()!')
                END IF

              CASE ('hotspot-hrd')

                ! We need a 2D noise field:

                IF (bub_index_noise(ii) > 0) THEN
                  CALL gen_bubnoise(fnoise=bub_tnoisebuffer(:,:,1,bub_index_noise(ii)), &
                       kdim=1, seed=505_iintegers + ii)
                ELSE
                  WRITE (*,*) 'ERROR: in artif_heatrate_dist(), '// &
                       'src_artifdata.f90: bub_index_noise(ii) == 0'
                  CALL model_abort (my_world_id, 10002+izerror, 'bub_index_noise('//cii// &
                       ') == 0 occured, but is wrong here!', &
                       'set correct bub_index_noise('//cii//') in input_artifctl()!')
                END IF

              END SELECT

            END IF   ! IF ( ladd_bubblenoise_t(ii) ) ...

          END IF   ! IF (bub_zeitzaehler(ii) == 1) ...


          ! Store RH before bubble if (probably) needed later:
          IF (lbub_rhconst(ii)) THEN
            ! storage of relat. humid. before heating:
            rhtmp_3d = rh_Tpqv_3d(p0(:,:,:)+pp(:,:,:,nx), t(:,:,:,nx), &
                                  zqv(:,:,:), zqv(:,:,:)*0.0_wp, ie, je, ke)
          END IF

          ! Per default disable QV-change; this might be switched on 
          ! for certain bubbles below,
          ! which then allows to keep the relhum constant during heating,
          ! if lbub_rhconst(ii)=.true.
          touch_moisture = .FALSE.

          ! Indicate if the bubble is 3D, then nk2d = -1, or if the bubble
          ! is 2D (z=const.), then nk2d = height level number of bubble.
          ! It defaults to -1 and may be changed below appropriately for certain
          ! bubble type(s).
          nk2d = -1

          SELECT CASE (TRIM(ctype_tempdist(ii)))

!!!=============================================================================

          CASE ('mcnider-hrd')

            touch_moisture = .TRUE.
            nk2d = -1

            CALL f_xyz_mcnider(ii, f_xyz_3D)

            IF (ladd_bubblenoise_t(ii)) THEN
              ttens(:,:,:) = ttens(:,:,:) + bub_heatingrate(ii) * f_xyz_3d(:,:,:) / &
                   MAX(MAXVAL(f_xyz_3d),1.0E-12_wp) * &
                   (1.0_wp + bub_dT_bubblenoise(ii) * &
                   bub_tnoisebuffer(:,:,1:ke,bub_index_noise(ii)))
            ELSE
              ttens(:,:,:) = ttens(:,:,:) + bub_heatingrate(ii) * f_xyz_3d(:,:,:) / &
                   MAX(MAXVAL(f_xyz_3d),1.0E-12_wp)
            END IF

!!!=============================================================================

          CASE ('AS2005_hucmtexas-hrd')

            touch_moisture = .FALSE.
            nk2d = -1

            CALL f_xyz_AS2005(ii, f_xyz_3d)

            ttens(:,:,:) = ttens(:,:,:) + bub_heatingrate(ii) * f_xyz_3d(:,:,:)


!!!=============================================================================

          CASE ('hotspot-hrd')
            
            ! gaussian heating rate disturbance in the lowest model layer,
            ! general 2D ellipse with the possibility to rotate horizontally.

            touch_moisture = .TRUE.
            nk2d = ke

            CALL f_xy_hotspot(ii, f_xyz_3d(:,:,1))


            IF (ladd_bubblenoise_t(ii)) THEN
              ttens(:,:,ke) = ttens(:,:,ke) + bub_heatingrate(ii) * f_xyz_3d(:,:,1) * &
                   (1.0_wp + bub_dT_bubblenoise(ii) * &
                   bub_tnoisebuffer(:,:,1,bub_index_noise(ii)))
            ELSE
              ttens(:,:,ke) = ttens(:,:,ke) + bub_heatingrate(ii) * f_xyz_3d(:,:,1)
            END IF


!!!=============================================================================

          CASE ('cos-hrd')

            touch_moisture = .TRUE.
            nk2d = -1

            CALL f_xyz_cos(ii, f_xyz_3d)

            IF (ladd_bubblenoise_t(ii)) THEN
              ttens(:,:,:) = ttens(:,:,:) + bub_heatingrate(ii) * f_xyz_3d(:,:,:) * &
                   (1.0_wp + bub_dT_bubblenoise(ii) * bub_tnoisebuffer(:,:,1:ke,bub_index_noise(ii)))
            ELSE
              ttens(:,:,:) = ttens(:,:,:) + bub_heatingrate(ii) * f_xyz_3d(:,:,:)
            END IF

!!!=============================================================================

          END SELECT

          IF (touch_moisture) THEN
            IF (lbub_rhconst(ii)) THEN

              ! conservation of relat. humid. during heating at constant pressure:
              ! (compute qv-tendency at constant p and relhum for a given T and T-tendency)
              IF (nk2d > -1) THEN
                DO j=1,je
                  DO i=1,ie
                    zqv_tens(i,j,nk2d) = zqv_tens(i,j,nk2d) + &
                         dqvdt_prh ( &
                             t(i,j,nk2d,nx), ttens(i,j,nk2d)-ttmp(i,j,nk2d), &
                             p0(i,j,nk2d)+pp(i,j,nk2d,nx), rhtmp_3d(i,j,nk2d), &
                             0.0_wp)
                  END DO
                END DO
              ELSE
                DO k=1,ke
                  DO j=1,je
                    DO i=1,ie
                      zqv_tens(i,j,k) = zqv_tens(i,j,k) + &
                           dqvdt_prh ( &
                               t(i,j,k,nx), ttens(i,j,k)-ttmp(i,j,k), &
                               p0(i,j,k)+pp(i,j,k,nx), rhtmp_3d(i,j,k), &
                               0.0_wp)
                    END DO
                  END DO
                END DO
              END IF

            END IF
          END IF

!!!=============================================================================

          IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
            CALL output3d_ascii(ttens(:,:,:)-ttmp, ke, &
                 'ttens-'//TRIM(ctype_tempdist(ii))//'-'//cii, &
                 'TTENS-disturbance atmosphere', 'K/s', 'eta')
            CALL output3d_ascii(zqv_tens(:,:,:)-qvtmp, ke, &
                 'qvtens-'//TRIM(ctype_tempdist(ii))//'-'//cii, &
                 'QVTENS-disturbance atmosphere', '1/s', 'eta')
          END IF

        ELSE IF (bub_zeitzaehler(ii) == bub_timespan(ii)) THEN

          ! ... If all disturbances are done, increment time counter by 1 
          !     so that it gets larger than bub_timespan:
          bub_zeitzaehler(ii) = bub_zeitzaehler(ii) + 1
          ! ... and set the disturbance master switch to .FALSE. With this,
          !     it is possible to check later if all disturbances are
          !     finished.
          ltempdist(ii) = .FALSE.

        END IF

      END DO bubbleloop

      ! ..clean up local memory:
      DEALLOCATE(f_xyz_3d)
      IF (alloc_ttmp) THEN
        DEALLOCATE(ttmp,qvtmp)
      END IF

      ! exchange temperature (for safety and for periodic BCs):
      kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries &
           (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
           kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
           lperi_x, lperi_y, l2dim, &
           20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
           ttens(:,:,:), zqv_tens(:,:,:) )

      ! Clean up memory after all bubbles have been released:
      IF (.NOT.ANY(ltempdist)) THEN
        IF (ALLOCATED(bub_tnoisebuffer)) THEN
          DEALLOCATE(bub_tnoisebuffer)
          !      WRITE (*,*) 'bub_tnoisebuffer deallocated!'
        END IF
      END IF

    END IF

  END SUBROUTINE artif_heatrate_dist


  !======================================================================================
  !
  ! ... Heating rate disturbance in the soil within any time interval during
  !     the simulation starting from "htempdist" for "bub_timespan" timesteps.
  !     Is called at the beginning of every timestep in lmorg.f90 in 
  !     SUBROUTINE initialize_loop().
  !     For t_so, there is no tendency variable in the model, so we just update
  !     t_so in the sense of Euler forward time integration instead of setting a tendency.
  !
  !     NOTE: makes only sense if run with prognostic soil model (lsoil=.true.).
  !           It recognizes the multi-layer soil model
  !
  !======================================================================================

  SUBROUTINE artif_heatrate_dist_tso(nx)

    IMPLICIT NONE

    ! time level (nnow or nnew)
    INTEGER(KIND=iintegers), INTENT(in) :: nx

    INTEGER(KIND=iintegers)      :: i, ii, jj, &
         anzbubs, zbub_index(ntempdist_max), &
         izerror, kzdims(24)
    REAL(KIND=wp),     ALLOCATABLE :: f_xyz_3d(:,:,:)

    REAL(KIND=wp),     ALLOCATABLE ::  ttmp(:,:,:)
    CHARACTER(len=2)  ::  cii

    CHARACTER (LEN=80)         ::  yzerrmsg


    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. artif_heatrate_dist_tso() ...'
    END IF


    yzerrmsg(:) = ' '      

    !*******************************************************************
    !   makes only sense if we run with prognostic soil model:
    !*******************************************************************

    IF (lsoil) THEN

      !*******************************************************************
      !            Disturbance on the heating rate in the soil
      !            during the course of the simulation.
      !            Because of lack of a soil temperature tendency
      !            field, the "tendency" is formulated as a Euler forward
      !            update of the soil temperature itself
      !*******************************************************************

      kzdims(:) = 0_iintegers
      
      zbub_index(:) = 0_iintegers

      anzbubs = 0
      DO ii=1, ntempdist_max
        IF (ltempdist(ii) .AND. bub_timespan(ii) > 0) THEN
          IF (ANY(reg_heatrate_dist_tsO == TRIM(ctype_tempdist(ii)))) THEN
            anzbubs = anzbubs + 1
            zbub_index(anzbubs) = ii
          END IF
        END IF
      END DO
      
      IF (anzbubs > 0) THEN
        
        ! ..allocate local working arrays:
        ALLOCATE(f_xyz_3d(ie,je,0:ke_soil+1))
        f_xyz_3d = 0.0_wp
        IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
          ALLOCATE(ttmp(ie,je,0:ke_soil+1))
          ttmp=0.0_wp
        END IF
            
        ! ..loop over all disturbances:
        bubbleloop: DO jj=1, anzbubs
          
          IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
            IF (lmulti_layer) THEN
              ttmp = t_so(:,:,:,nx)
            ELSE
              ttmp(:,:,1) = t_s(:,:,nx)
              ttmp(:,:,2) = t_m(:,:,nx)
              ttmp(:,:,3) = t_cl(:,:)
            END IF
          END IF

          ii = zbub_index(jj)
          
          WRITE(cii,'(i2.2)') ii

          IF ( (ntstep_bubble(ii) == ntstep .OR. &
               (bub_zeitzaehler(ii) > 0 .AND. bub_zeitzaehler(ii) < bub_timespan(ii)))) THEN

            bub_zeitzaehler(ii) = bub_zeitzaehler(ii) + 1

            IF (my_cart_id == 0) THEN
              WRITE (*,*) ('=',i=1,70)
              WRITE (*,*) '        Soil heating rate disturbance ', ii, ' , Typ ', &
                   TRIM(ctype_tempdist(ii)),' , disturbance-time-step ', bub_zeitzaehler(ii), ' ...'
              WRITE (*,*) ('=',i=1,70)
            END IF

            IF (bub_zeitzaehler(ii) == 1) THEN
              
              IF (ladd_bubblenoise_t(ii)) THEN
                
                ! Certain disturbance types may be modulated by some white noise
                ! (if ladd_bubblenoise_t(ii) = .true.). 
                ! For this noise, we need a storage field, and for technical
                ! reasons (ease of implementation when integrating all
                ! possible combinations of namelist switches), this is a global
                ! persistent array which needs to be allocated once:
                IF (.NOT. ALLOCATED(bub_tnoisebuffer)) THEN
                  ! DO NOT TOUCH THE ALLOCATION INDICES, THEY HAVE TO BE THE SAME 
                  ! FOR ALL DISTURBANCE ROUTINES!
                  ALLOCATE(bub_tnoisebuffer(ie,je,MAX(ke,ke_soil),1:numbubs_noise))
                  !            WRITE (*,*) 'bub_tnoisebuffer allokiert!'
                  bub_tnoisebuffer = 0.0_wp
                END IF

                ! To generate a reproducible noise on different numbers of
                ! processors, only PE 0 generates the noise on a global field
                ! and distributes it to the other PEs. This is done in SR gen_bubnoise.
                ! The noise is stored in the global field bub_tnoisebuffer(:,:,kstart:kend,bub_index_noise(ii)).
                
                SELECT CASE (TRIM(ctype_tempdist(ii)))
                  
                CASE ('cos-soil-hrd')
                  
                  ! We need a 3D noise field:

                  IF (bub_index_noise(ii) > 0) THEN
                    CALL gen_bubnoise(fnoise=bub_tnoisebuffer(:,:,1:ke_soil,bub_index_noise(ii)), &
                         kdim=ke_soil, seed=505_iintegers + ii)
                  ELSE
                    WRITE (*,*) 'ERROR: in artif_heatrate_dist_tso(), '// &
                         'src_artifdata.f90: bub_index_noise(ii) == 0'
                    izerror = 1
                    CALL model_abort (my_world_id, 10002+izerror, 'bub_index_noise('//cii// &
                         ') == 0 occured, but is wrong here!', &
                         'set correct bub_index_noise('//cii//') in input_artifctl()!')
                  END IF
                  
                CASE ('hotspot-soil-hrd')
                  
                  ! We need a 2D noise field:

                  IF (bub_index_noise(ii) > 0) THEN
                    CALL gen_bubnoise(fnoise=bub_tnoisebuffer(:,:,1,bub_index_noise(ii)), &
                         kdim=1, seed=505_iintegers + ii)
                  ELSE
                    WRITE (*,*) 'ERROR: in artif_heatrate_dist_tso(), '// &
                         'src_artifdata.f90: bub_index_noise(ii) == 0'
                    CALL model_abort (my_world_id, 10002+izerror, 'bub_index_noise('//cii// &
                         ') == 0 occured, but is wrong here!', &
                         'set correct bub_index_noise('//cii//') in input_artifctl()!')
                  END IF
                  
                END SELECT
                
              END IF   ! IF ( ladd_bubblenoise_t(ii) ) ...
              
            END IF   ! IF (bub_zeitzaehler(ii) == 1) ...
            
            SELECT CASE (TRIM(ctype_tempdist(ii)))

!!!=============================================================================

            CASE ('hotspot-soil-hrd')

              ! gaussian temperature disturbance in the lowest model layer,
              ! general 2D ellipse with the possibility to rotate horizontally.

              CALL f_xy_hotspot(ii, f_xyz_3d(:,:,1))
              
              IF (lmulti_layer) THEN
                IF (ladd_bubblenoise_t(ii)) THEN
                  t_so(:,:,0,nx) = t_so(:,:,0,nx) + bub_heatingrate(ii) * &
                       f_xyz_3D(:,:,1) * dt * &
                       (1.0_wp + bub_dT_bubblenoise(ii) * &
                       bub_tnoisebuffer(:,:,1,bub_index_noise(ii)))
                ELSE
                  t_so(:,:,0,nx) = t_so(:,:,0,nx) + bub_heatingrate(ii) * &
                       f_xyz_3D(:,:,1) * dt
                END IF
                t_so(:,:,1,nx)         = t_so(:,:,0,nx)
              ELSE
                IF (ladd_bubblenoise_t(ii)) THEN
                  t_s(:,:,nx) = t_s(:,:,nx) + bub_heatingrate(ii) * f_xyz_3D(:,:,1) * dt * &
                       (1.0_wp + bub_dT_bubblenoise(ii) * &
                       bub_tnoisebuffer(:,:,1,bub_index_noise(ii)))
                ELSE
                  t_s(:,:,nx) = t_s(:,:,nx) + bub_heatingrate(ii) * f_xyz_3D(:,:,1) * dt
                END IF
              END IF

!!!=============================================================================

            CASE ('cos-soil-hrd')
              
              CALL f_xyz_cos2d_soil(ii, f_xyz_3d(:,:,1:ke_soil))
              
              IF (lmulti_layer) THEN
                IF (ladd_bubblenoise_t(ii)) THEN
                  t_so(:,:,1:ke_soil,nx) = t_so(:,:,1:ke_soil,nx) + bub_heatingrate(ii) * &
                       f_xyz_3d(:,:,1:ke_soil) * dt * &
                       (1.0_wp + bub_dT_bubblenoise(ii) * &
                       bub_tnoisebuffer(:,:,1:ke_soil,bub_index_noise(ii)))
                ELSE
                  t_so(:,:,1:ke_soil,nx) = t_so(:,:,1:ke_soil,nx) + bub_heatingrate(ii) * &
                       f_xyz_3d(:,:,1:ke_soil) * dt
                END IF
                t_so(:,:,ke_soil+1,nx) = t_so(:,:,ke_soil,nx)
                t_so(:,:,0,nx)         = t_so(:,:,1,nx)
                IF (lalloc_t_cl) THEN
                  IF (ladd_bubblenoise_t(ii)) THEN
                    t_cl(:,:) = t_cl(:,:)  + bub_heatingrate(ii) * &
                         f_xyz_3d(:,:,ke_soil) * dt * &
                         (1.0_wp + bub_dT_bubblenoise(ii) * &
                         bub_tnoisebuffer(:,:,ke_soil,bub_index_noise(ii)))
                  ELSE
                    t_cl(:,:) = t_cl(:,:)  + bub_heatingrate(ii) * &
                         f_xyz_3d(:,:,ke_soil) * dt
                  END IF
                END IF
              ELSE
                IF (ladd_bubblenoise_t(ii)) THEN
                  t_s(:,:,nx) = t_s(:,:,nx) + bub_heatingrate(ii) * f_xyz_3D(:,:,1) * dt * &
                       (1.0_wp + bub_dT_bubblenoise(ii) * &
                       bub_tnoisebuffer(:,:,1,bub_index_noise(ii)))
                  t_m(:,:,nx) = t_m(:,:,nx) + bub_heatingrate(ii) * f_xyz_3D(:,:,2) * dt * &
                       (1.0_wp + bub_dT_bubblenoise(ii) * &
                       bub_tnoisebuffer(:,:,2,bub_index_noise(ii)))
                  t_cl(:,:)   = t_cl(:,:)   + bub_heatingrate(ii) * f_xyz_3D(:,:,3) * dt * &
                       (1.0_wp + bub_dT_bubblenoise(ii) * &
                       bub_tnoisebuffer(:,:,3,bub_index_noise(ii)))
                ELSE
                  t_s(:,:,nx) = t_s(:,:,nx) + bub_heatingrate(ii) * f_xyz_3D(:,:,1) * dt
                  t_m(:,:,nx) = t_m(:,:,nx) + bub_heatingrate(ii) * f_xyz_3D(:,:,2) * dt
                  t_cl(:,:)   = t_cl(:,:)   + bub_heatingrate(ii) * f_xyz_3D(:,:,3) * dt
                END IF
              END IF

!!!=============================================================================

            END SELECT
            
            IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
              IF (lmulti_layer) THEN
                f_xyz_3d(:,:,:) = t_so(:,:,:,nx) - ttmp
                CALL output3d_ascii(f_xyz_3d(:,:,1:ke_soil+1), ke_soil+1, &
                     'tso-'//TRIM(ctype_tempdist(ii))//'-'//cii, &
                     'T-disturbance soil (based on heating rate)', 'K', 'eta')
              ELSE
                f_xyz_3d(:,:,1) = t_s(:,:,nx) - ttmp(:,:,1)
                f_xyz_3d(:,:,2) = t_m(:,:,nx) - ttmp(:,:,2)
                f_xyz_3d(:,:,3) = t_cl(:,:)   - ttmp(:,:,3)
                CALL output3d_ascii(f_xyz_3d(:,:,1:3), 3_iintegers, &
                     'tso-'//TRIM(ctype_tempdist(ii))//'-'//cii, &
                     'T-disturbance soil  (based on heating rate)', 'K', 'eta')
              END IF
            END IF

          ELSE IF (bub_zeitzaehler(ii) == bub_timespan(ii)) THEN
            
            ! ... If all disturbances are done, increment time counter by 1 
            !     so that it gets larger than bub_timespan:
            bub_zeitzaehler(ii) = bub_zeitzaehler(ii) + 1
            ! ... and set the disturbance master switch to .false. With this,
            !     it is possible to check later if all disturbances are
            !     finished.
            ltempdist(ii) = .FALSE.
          
          END IF  ! IF ( (ntstep_bubble(ii) == ntstep .OR. ...

        END DO bubbleloop

        ! ..clean up local memory:
        DEALLOCATE(f_xyz_3d)
        IF (ldebug_artif .AND. idbg_artif_level > 3) THEN
          DEALLOCATE(ttmp)
        END IF

        ! exchange temperature (for safety and for periodic BCs):
        IF (lmulti_layer) THEN
           IF (.NOT. lalloc_t_cl) THEN
            kzdims(1:24)=(/ke_soil+2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
            CALL exchg_boundaries  &
                 (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
                 kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
                 lperi_x, lperi_y, l2dim, &
                 20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
                 t_so(:,:,:,nx) )
          ELSE
            kzdims(1:24)=(/ke_soil+2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
            CALL exchg_boundaries  &
                 (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
                 kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
                 lperi_x, lperi_y, l2dim, &
                 20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
                 t_so(:,:,:,nx), t_cl(:,:) )
          END IF
        ELSE
          kzdims(1:24)=(/1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
          CALL exchg_boundaries &
               (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
               kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
               lperi_x, lperi_y, l2dim, &
               20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
               t_s(:,:,nx), t_m(:,:,nx), t_cl(:,:) )
        END IF

        ! Clean up memory after all bubbles have been released:
        IF (.NOT.ANY(ltempdist)) THEN
          IF (ALLOCATED(bub_tnoisebuffer)) THEN
            DEALLOCATE(bub_tnoisebuffer)
            !      WRITE (*,*) 'bub_tnoisebuffer deallokiert!'
          END IF
        END IF

      END IF  ! IF (anzbubs > 0) 

    END IF  ! IF (lsoil)
    
  END SUBROUTINE artif_heatrate_dist_tso


!!!=============================================================================
!!!=============================================================================

  !=======================================================================
  !
  ! SUBROUTINE set_idealized_surffluxes()
  !
  ! Specification of pre-defined surface sensible and latent heat fluxes.
  ! Some spatial random noise can be added to these fluxes if desired. This
  ! noise is kept constant in time to mimick surface inhomogeneities.
  !
  ! METHOD: specify a suitable combination of the transfer coefficient tch
  !         and the surface temperature t_g and moisture qv_s in order to 
  !         obtain the desired fluxes from the COSMO flux parameterization,
  !         for the actual lowest atmospheric temperature and windspeed.
  !         First, t_g is adjusted using a typical linear dependence
  !         of sensible heat flux from the bulk temperature difference
  !         between the surface and the lowest model level.
  !         Then, the transfer coefficient tch is adjusted to yield
  !         the desired flux. Note that tch depends on qv_s to account
  !         for the virtual temperature, which is used in the flux
  !         parameterization.
  !         If a latent heat flux is desired, qv_s is adjusted afterwards
  !         to yield the desired flux from the above value of tch.
  ! 
  ! NOTE:   The result for the fluxes might not be exactly accurate,
  !         because tch for the sensible heat flux depends on qv_s(nnow),
  !         which might afterwards be changed again to yield the desired
  !         latent heat flux. Normally a fixpoint interation would be needed to
  !         resolve this isssue, but we do without it here and accept the
  !         small flux errors.
  !
  ! USAGE:  This subroutine is called in turbulence_interface.f90.
  !         It works for runs without the soil model (lsoil=.false.)
  !         and for itype_tran=1 or 2. It seems not to work with
  !         itype_turb=3 and itype_tran=3.
  !
  ! RELEVANT NAMELIST PARAMETER:
  !
  !         lsensiflux_fix:        Master switch, turn on sensible heat flux H0
  !
  !           llatentflux_fix:     If lsensiflux_fix=.true., turn on/off latent heat flux
  !
  !           sensiflux_c:         Sensible heat flux in W/m^2
  !
  !           latentflux_LzuS:     Bowen-ratio (latent to sensible), if llatentflux_fix=.true.
  !
  !           H0_rel_noise:        Spatial white noise level on the fluxes, constant in time
  !
  !                        The noise is computed from
  !                  H0_noisy = H0 * ( 1 + H0_rel_noise * whitenoise[-1,1] )
  !
  !           iseed_noise_H0:      INTEGER seed for the random number generator (-999 = system time is used)
  !
  !
  !=======================================================================

!!! timelevel for leapfrog integration for ntstep = 0 ?

  SUBROUTINE set_idealized_surffluxes()

    IMPLICIT NONE

    INTEGER(KIND=iintegers) :: i, j, iu, jv, izerror
    INTEGER(KIND=iintegers), SAVE :: firstcall = 0
    REAL(KIND=wp),     SAVE :: zscal_tdiff
    REAL(KIND=wp),     ALLOCATABLE, SAVE :: h0noise(:,:)
    REAL(KIND=wp)     :: sensifluxtmp(ie,je), zdx, zdy, sswdown
    REAL(KIND=wp)     :: latentfluxtmp(ie,je)
    CHARACTER (LEN=80)         ::  yzerrmsg, yzroutine
    REAL(KIND=wp)     :: zvbke, ztvb, ztch

    REAL (KIND=wp),     POINTER :: qv_new  (:,:,:) => NULL() ! QV at nnew
    REAL (KIND=wp),     POINTER :: qv_now  (:,:,:) => NULL() ! QV at nnow

    yzroutine(:) = ' '
    yzroutine = 'set_constant_surffluxes'

    yzerrmsg(:) = ' '


    IF (lsensiflux_fix) THEN

      IF (llatentflux_fix) THEN
        ! Retrieve the microphysics tracers
        CALL trcr_get( izerror, idt_qv, ptr_tlev = nnew, ptr = qv_new )
        IF ( izerror /= 0_iintegers ) THEN
          yzerrmsg = trcr_errorstr( izerror )
          CALL model_abort( my_cart_id, izerror, yzerrmsg, TRIM(yzroutine) )
        ENDIF
        CALL trcr_get( izerror, idt_qv, ptr_tlev = nnow, ptr = qv_now )
        IF ( izerror /= 0_iintegers ) THEN
          yzerrmsg = trcr_errorstr( izerror )
          CALL model_abort( my_cart_id, izerror, yzerrmsg, TRIM(yzroutine) )
        ENDIF
      END IF

      IF (firstcall .NE. 1) THEN
        IF (my_cart_id == 0) THEN
          WRITE (*,*) ('=',i=1,70)
          WRITE (*,*) '        Setze Waermefluss am Boden auf', sensiflux_c, ' W/m^2'
          WRITE (*,*) ('=',i=1,70)
        END IF
        ! Scaling factor for the surface overheating, linearily varying with the heat flux:
        zscal_tdiff = 10.0_wp / 300.0_wp

        ! If desired, overlay spatial white noise on the fluxes, constant in time,
        ! mimicking some surface inhomogeneities:
        ALLOCATE(h0noise(ie,je))
        IF (H0_rel_noise >= 1.0E-5_wp) THEN

          CALL gen_bubnoise(fnoise=h0noise(:,:), kdim=1, seed=iseed_noise_h0)
          DO j= 1, je
            DO i= 1, ie
              h0noise(i,j) = 1.0_wp + 2.0_wp * (h0noise(i,j) - 0.5_wp) * H0_rel_noise
            END DO
          END DO
        ELSE
          h0noise = 1.0_wp
        END IF

        firstcall = 1
      END IF

      ! Time-constant sensible heat flux:
      sensifluxtmp  = sensiflux_c * h0noise


      ! The surface temperature must be warme/cooler as the adiabatically
      ! scaled temperature from the lowest model level. We use a
      ! linearily scaling temperature difference of 10 K per 300 W/m^2, and
      ! the sign of the difference is dictated by the sign of the heat flux.
      ! In the following loop, the transfer coefficient tch is determined in
      ! such a way that the desired heat flux results from the above temperature
      ! difference.
      DO j=1,je
        jv = MAX(j-1, 1)
        DO i=1,ie
          iu = MAX(i-1, 1)
          t_g(i,j,nnow) = t(i,j,ke,nnow) * ((ps(i,j,nnow))/(p0(i,j,ke)+pp(i,j,ke,nnow)))**rdocp + &
               sensifluxtmp(i,j)*zscal_tdiff
          tch(i,j) = ABS(sensifluxtmp(i,j)) / ( ps(i,j,nnow)/(r_d*t_g(i,j,nnow)*(1.0_wp + rvd_m_o*qv_s(i,j,nnow))) * &
               MAX(0.5_wp*SQRT ( (u(i,j,ke,nnow) + u(iu,j,ke,nnow))**2 + &
               (v(i,j,ke,nnow) + v(i,jv,ke,nnow))**2 ), vel_min) * &
               cp_d * MAX(ABS(sensifluxtmp(i,j))*zscal_tdiff, 1.0E-5_wp) )
        END DO
      END DO

      t_g(:,:,nnew) = t_g(:,:,nnow)
      t_s(:,:,nnow) = t_g(:,:,nnow)
      t_s(:,:,nnew) = t_g(:,:,nnow)
      IF (.NOT. lmulti_layer) THEN
        t_s_bd(:,:,1) = t_g(:,:,nnow)
        t_s_bd(:,:,2) = t_g(:,:,nnow)
      END IF

      ! time-constant latent heat flux using the pre-specified bowen ratio:
      IF (llatentflux_fix) THEN

        latentfluxtmp(:,:) = latentflux_LzuS * sensifluxtmp(:,:)

        ! The specific humidity at ground must be larger/smaller than the
        ! specific humidity in the lowest model level in order
        ! to get a latent heat flux from/to the surface.
        ! In the following loop, it will be set in such a way
        ! that the desired latent flux results from the
        ! above transfer coefficient tch.
        DO j=1,je
          jv = MAX(j-1, 1)
          DO i=1,ie
            iu = MAX(i-1, 1)
            zvbke      = MAX( 0.5_wp*SQRT ( (u(i,j,ke,nnow) + u(iu,j,ke,nnow))**2    &
                 + (v(i,j,ke,nnow) + v(i,jv,ke,nnow))**2 ), vel_min)
            ztvb       = t_s (i,j,nnow)*(1.0_wp + rvd_m_o*qv_s(i,j,nnow))

            ztch       = tch(i,j)*zvbke*ps(i,j,nnow)/(r_d*ztvb)
            IF (ABS(ztch) >= 1.0E-20_wp ) THEN
              qv_s(i,j,nnow) = qv_now(i,j,ke) + latentfluxtmp(i,j) / (lh_v * ztch)
            ELSE
              qv_s(i,j,nnow) = qv_now(i,j,ke)
            END IF
          END DO
        END DO

        qv_s(:,:,nnew) = qv_s(:,:,nnow)
        qv_s_bd(:,:,nnow) = qv_s(:,:,nnow)
        qv_s_bd(:,:,nnew) = qv_s(:,:,nnow)

      END IF

    END IF

  END SUBROUTINE set_idealized_surffluxes


!!!=============================================================================
!!!=============================================================================


  !=============================================================================
  !
  ! Add white noise (horizontally uniformly distributed) to the 
  ! temperature and w field in the boundary layer (lowest 100 hPa).
  ! (only at a single time step at the moment)
  !
  !=============================================================================

  SUBROUTINE add_noise_tw(nx)

    IMPLICIT NONE

    ! time level (nnow or nnew)
    INTEGER(KIND=iintegers), INTENT(in) :: nx

    INTEGER(KIND=iintegers)        :: i, j, k, izerror, kzdims(24)
    REAL(KIND=wp),     ALLOCATABLE :: tmpnoise(:,:,:)
    REAL(KIND=wp)                  :: ztnew
    CHARACTER (LEN=250)            :: yzerrmsg


    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. add_noise_tw() ...'
    END IF

    yzerrmsg(:) = ' '

    IF (ladd_noise_t .AND. ntstep_noise == ntstep) THEN

      IF (my_cart_id == 0) THEN
        WRITE (*,*) ('=',i=1,70)
        WRITE (*,*) '        Generate T-Noise with ', dT_noise, ' K   amplitude'
        WRITE (*,*) '        Generate W-Noise with ', dW_noise, ' m/s amplitude'
        WRITE (*,*) ('=',i=1,70)
      END IF

      ! allocate temporary global noise field 2D:
      ALLOCATE(tmpnoise(ie,je,ke))
      tmpnoise = 0.0_wp      

      CALL gen_bubnoise(fnoise=tmpnoise(:,:,:), kdim=ke, seed=iseed_noise_t)

      ! Add the T-noise in the boundary layer (lowest 100 hPa):
      DO k = ke, 1, -1
        DO j = 1, je
          DO i = 1, ie
            ztnew = ABS((ps(i,j,nx)-p0(i,j,k)-pp(i,j,k,nx)))
            IF (ztnew < 1.0E4_wp) THEN
              ztnew = tmpnoise(i,j,k) * dT_noise * COS(0.5_wp*pi*ztnew)
              t(i,j,k,nx) = t(i,j,k,nx) + ztnew
            END IF
          END DO
        END DO
      END DO

      CALL gen_bubnoise(fnoise=tmpnoise(:,:,:), kdim=ke, seed=iseed_noise_w)

      ! Add the W-noise in the boundary layer (lowest 100 hPa):
      DO k = ke, 1, -1
        DO j = 1, je
          DO i = 1, ie
            ztnew = ABS((ps(i,j,nx)-p0(i,j,k)-pp(i,j,k,nx))*1.0E-4_wp)
            IF (ztnew < 1.0_wp) THEN
              ztnew =  tmpnoise(i,j,k) * dW_noise * COS(0.5_wp*pi*ztnew)
              w(i,j,k,nx) = w(i,j,k,nx) +ztnew
            END IF
          END DO
        END DO
      END DO

      ! ..clean up memory:
      DEALLOCATE(tmpnoise)

      ! For safety and for periodic BCs: boundary exchange of w and t:
      kzdims(1:24)=(/ke1,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
      CALL exchg_boundaries &
           (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je, &
           kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
           lperi_x, lperi_y, l2dim, &
           20000, .FALSE., ncomm_type, izerror, yzerrmsg, &
           w(:,:,:,nx), t(:,:,:,nx) )

    END IF


  END SUBROUTINE add_noise_tw

!==============================================================================
!==============================================================================
!
! Subroutines for defining the spatial shape function of temperature
! disturbances
!
!==============================================================================
!==============================================================================

  !----------------------------------------------------------------------------
  ! Shape function for a cos^2 temperature disturbance 
  ! after Weisman & Klemp (1982)
  !----------------------------------------------------------------------------

  ! The cos^2-disturbance may be elliptic and may be horizontally rotated.

  SUBROUTINE f_xyz_cos(ii, f_xyz)

    IMPLICIT NONE

    ! Index of bubble in bubble namelist vectors:
    INTEGER(KIND=iintegers), INTENT(in)  :: ii
    ! spatial shape function of the bubble:
    REAL(KIND=wp),           INTENT(out) :: f_xyz(ie,je,ke)

    INTEGER(KIND=iintegers)      :: i, j, k, i_td, j_td
    REAL(KIND=wp)                :: zdx_rot, zdy_rot, zfak, ztmpwinkel

    f_xyz(:,:,:) = 0.0_wp

    ztmpwinkel = bub_rotangle(ii) * degrad

    j_td  = isubpos(my_cart_id, 2) - nboundlines - 1
!CDIR NOLOOPCHG
    DO j = 1, je
      j_td = j_td + 1
      i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
!CDIR NOLOOPCHG
      DO i = 1, ie
        i_td = i_td + 1
        ! compute rectangular distances to the (rotated) main axes
        ! of the bubble (zdx_rot, zdy_rot):
        CALL hill_rot_coords(i_td, j_td, bub_centi(ii), bub_centj(ii), &
             ztmpwinkel, 0.0_wp, zdx_rot, zdy_rot)
!CDIR NOLOOPCHG
        DO k = 1, ke
          zfak  = SQRT( (zdx_rot/bub_radx(ii))**2 + &
               (zdy_rot/bub_rady(ii))**2 + &
               (( 0.5_wp*(hhl(i,j,k)+hhl(i,j,k+1)) - bub_centz(ii) )/bub_radz(ii))**2 )
          IF (zfak < 1.0_wp) THEN
            f_xyz(i,j,k) = COS(pi*zfak/2.0_wp)**2
          END IF
        ENDDO
      ENDDO
    ENDDO
    
  END SUBROUTINE f_xyz_cos

  !----------------------------------------------------------------------------
  ! Shape function for McNider & Kopp (1990) temperature disturbance 
  !----------------------------------------------------------------------------

  ! The McNider-disturbances are horizontally circular, 
  ! so rotation does not make sense here.
  ! The McNider-disturbances are formulated terrain following, 
  ! in cotrast to the cos^2-disturbances.

  SUBROUTINE f_xyz_mcnider(ii, f_xyz)

    IMPLICIT NONE

    ! Index of bubble in bubble namelist vectors:
    INTEGER(KIND=iintegers), INTENT(in)  :: ii
    ! spatial shape function of the bubble:
    REAL(KIND=wp),           INTENT(out) :: f_xyz(ie,je,ke)

    INTEGER(KIND=iintegers)      :: i, j, k, i_td, j_td
    REAL(KIND=wp)                :: zA, zlambda_m, ztheta_l, zrho_l, zsigma_theta, &
         zdx, zdy, zfak, zheight

    zA = 3.0_wp;
    zlambda_m = 1.5_wp * bub_zi_mcnider(ii)
    ztheta_l = 303.16_wp
    zrho_l = 1.22_wp


    f_xyz(:,:,:) = 0.0_wp

    j_td  = isubpos(my_cart_id, 2) - nboundlines - 1
!CDIR NOLOOPCHG
    DO j = 1, je
      j_td = j_td + 1
      i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
!CDIR NOLOOPCHG
      DO i = 1, ie
        i_td = i_td + 1
        ! compute rectangular distances to the main axes
        ! of the bubble (zdx, zdy):
        CALL hill_rot_coords(i_td, j_td, bub_centi(ii), bub_centj(ii), &
             0.0_wp, 0.0_wp, zdx, zdy)
        zfak  = EXP(-1.0_wp*(zdx**2 + zdy**2) / (0.5_wp*zlambda_m)**2)
        IF (zfak >= 1.0E-2_wp) THEN
!CDIR NOLOOPCHG
          DO k = 1, ke
            zheight = 0.5_wp*(hhl(i,j,k)+hhl(i,j,k+1)) - hhl(i,j,ke+1)
            IF (zheight <= bub_zmax_mcnider(ii)) THEN
              zsigma_theta = 1.34_wp * MAX(zheight,0.1_wp)**(-1.0_wp/3.0_wp) * &
                   (bub_h0_mcnider(ii)/(zrho_l*cp_d))**(2.0_wp/3.0_wp) * (ztheta_l/g)**(1.0_wp/3.0_wp)
              f_xyz(i,j,k) = zA * zsigma_theta * zfak
            END IF
          ENDDO
        END IF
      ENDDO
    ENDDO

  END SUBROUTINE f_xyz_mcnider

  
  !----------------------------------------------------------------------------
  ! Shape function for the "Squall3D"-testcase temperature disturbance
  ! of G. Bryan (NCAR)
  !----------------------------------------------------------------------------

!!! Settings from Axels code:
!!$  zbub_dT    = 1.5_wp      ![K]
!!$  zbub_radx  = 10000._wp   ![m]
!!$  zbub_radz  = 1400._wp    ![m]
!!$  zbub_centi = 300             ![gridpoints]
!!$  zbub_centj = 40              ![gridpoints]
!!$  zbub_centz = 1400._wp    ![m]

  SUBROUTINE f_xyz_squall3d(ii, f_xyz)

    IMPLICIT NONE

    ! Index of bubble in bubble namelist vectors:
    INTEGER(KIND=iintegers), INTENT(in)  :: ii
    ! spatial shape function of the bubble:
    REAL(KIND=wp),           INTENT(out) :: f_xyz(ie,je,ke)

    INTEGER(KIND=iintegers)      :: i, j, k, i_td, j_td, jj, mpierror
    REAL(KIND=wp)                :: zdx, zdy, zfak, zdi, zdj

    ! variables for squallline initialization:
    REAL(KIND=wp)       :: &
         v_mag               ,&! max magnitude (K) of random variation
         bhrad               ,&! x,y radii of random variation (m)
         bvrad               ,&! z radius of random variation (m)
         ric                 ,&! x location of random variation (grid pts)
         rjc_all(1000)       ,&! postions of perturbations (grid pts)
         rndm_nmbr(1000)     ,&! just some random numbers
         beta, rjc
    INTEGER(KIND=iintegers) :: &
         rjc_nmbr              ! number of centers of slight perturbations  

    REAL(KIND=wp),     ALLOCATABLE :: noisedummy(:)


    f_xyz(:,:,:) = 0.0_wp

    zdy   = r_earth* dlat * degrad
    ! Establish radii for the small, random variations on the perturbation used
    ! to trigger convection.  This current implementation is for one or more
    ! perturbations at a fixed x.

    v_mag             =      0.1_wp  ! max. rel. magnitude of random variation.
    bhrad             =      bub_radx(ii) * 0.25_wp  ! x,y radii of random variation (m)
    bvrad             =      bub_radz(ii) + bub_centz(ii)  ! z radius of random variation (m)
    ric               =      bub_centi(ii) ! x location of random variation (grid pts)

    rjc_nmbr = 0_iintegers
    rjc_all = 0.0_wp
    DO
      rjc_nmbr = rjc_nmbr + 1_iintegers
      IF (rjc_nmbr > 500) THEN
        rjc_nmbr = 500
        EXIT
      END IF
      rjc_all(rjc_nmbr) = bhrad * (2.0_wp*rjc_nmbr) / zdy + nboundlines + 1
      IF (rjc_all(rjc_nmbr) > je_tot - bhrad * 2.0_wp / zdy - nboundlines - 1) THEN
        rjc_nmbr = rjc_nmbr - 1
        EXIT
      END IF
    END DO

    rjc = rjc_all(1)  ! initial setting only


    ! Generate the seed that will be used to produce the
    ! random variations in the initial cold pool or warm bubble.  Then fill the
    ! array with a random number for each grid point in y.
    ! NOTE: this random number series has to be the same on each processor,
    ! so get the numbers on PE 0 and send them to all other nodes.
    IF (my_cart_id == 0) THEN

      CALL seed_random_number( 404_iintegers )
      ! on some compilers, the random series shows problems within the first few
      ! hundred random numbers (the numbers are not really random, but
      ! can be monotonic and reproducibly the same on all processors).
      ! Only after some numbers, the series gets more random type.
      ! Therefore, fetch 5000 dummy random numbers, before doing the
      ! really needed random numbers:
      ALLOCATE(noisedummy(5000))
      CALL RANDOM_NUMBER(noisedummy)
      DEALLOCATE(noisedummy)
      
      CALL RANDOM_NUMBER(rndm_nmbr)
    END IF

    IF (num_compute > 1) THEN
      ! Distribute random numbers to all PEs:
      CALL distribute_values(rndm_nmbr, SIZE(rndm_nmbr,1), 0, imp_reals, icomm_cart, mpierror)
    END IF

    ! Calculate distances from each model grid point to the pre-determined center
    ! points or lines of the thermal perturbation (i.e., the trigger, with or
    ! without random variations).  

    j_td  = isubpos(my_cart_id, 2) - nboundlines - 1
!CDIR NOLOOPCHG
    DO j = 1, je
      j_td = j_td + 1
      zdj  = REAL(j_td,wp) - REAL(je_tot,wp) / 2.0_wp
      i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
      zdx   = r_earth* dlat * degrad * crlat(j,1)
!CDIR NOLOOPCHG
      DO i = 1, ie
        i_td = i_td + 1
        zdi  = REAL(i_td,wp)-REAL(bub_centi(ii),wp)
!CDIR NOLOOPCHG
        DO k = 1, ke
          zfak  = ( zdi*zdx/bub_radx(ii) )**2 + &
               ( ( 0.5_wp*(hhl(i,j,k)+hhl(i,j,k+1)) - bub_centz(ii))/bub_radz(ii) )**2
          zfak  = SQRT(zfak)
          
          IF ( zfak < 1.0_wp ) THEN
            
            ! Calculate for each model grid point whether the point falls within the
            ! radius of a random variation on the temperature perturbation.  If it does,
            ! record the index of the center of the particular random variation.
            
            DO jj = 1, rjc_nmbr
              IF ( ABS( (rjc_all(jj)-REAL(j_td,wp))*zdy ) <= bhrad ) THEN
                rjc = rjc_all(jj)
                EXIT
              END IF
            END DO
            beta = SQRT( ( REAL(i_td-INT(ric),wp) * zdx  / bhrad )**2 &
                 + ( REAL(j_td-INT(rjc),wp) * zdy  / bhrad )**2 &
                 + ( 0.5_wp*(hhl(i,j,k)+hhl(i,j,k+1))  / bvrad )**2 )
            
            IF ( beta < 1.0_wp ) THEN
              f_xyz(i,j,k) = COS(pi*zfak/2.0_wp)**2 + &
                   v_mag*rndm_nmbr(INT(rjc))*(COS(0.5_wp*pi*beta)**2)* &
                   ((p0(i,j,k)+pp(i,j,k,nnew))/refatm%p0sl)**rdocp
            ELSE
              f_xyz(i,j,k) = COS(pi*zfak/2.0_wp)**2
            END IF
          END IF
        ENDDO
      ENDDO
    ENDDO
    
    IF (my_cart_id == 0) THEN
      WRITE (*,*) '=== SQALL3D-Testcase:'
      WRITE (*,*) '===  zdx = ',r_earth * dlat * degrad * COS( (startlat_tot+(je_tot/2)*dlat) * degrad)
      WRITE (*,*) '===  zdy = ',zdy
      WRITE (*,*) '===  bub_centz = ',bub_centz(ii)
      WRITE (*,*) '===  bub_radz  = ',bub_radz(ii)
    ENDIF
    
  END SUBROUTINE f_xyz_squall3d



  SUBROUTINE f_xyz_SK94(ii, f_xyz)

    IMPLICIT NONE

    ! Index of bubble in bubble namelist vectors:
    INTEGER(KIND=iintegers), INTENT(in)  :: ii
    ! spatial shape function of the bubble:
    REAL(KIND=wp),           INTENT(out) :: f_xyz(ie,je,ke)

    INTEGER(KIND=iintegers)      :: i, j, k, i_td
    REAL(KIND=wp)                :: x, zml, zdx, zdi

    REAL(KIND=wp)       :: model_height


    f_xyz(:,:,:) = 0.0_wp

    model_height = hhl(1,1,1) - hhl(1,1,ke+1)

!CDIR NOLOOPCHG
    DO j = 1, je
      zdx   = r_earth* dlon * degrad  * crlat(j,1)
      i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
!CDIR NOLOOPCHG
      DO i = 1, ie
        i_td = i_td + 1
        zdi  = REAL(i_td-bub_centi(ii),wp)
        x = zdx * zdi    ! = x-xc
!CDIR NOLOOPCHG
        DO k = 1, ke      
          zml = 0.5_wp * ( hhl(i,j,k) + hhl(i,j,k+1) )
          
          f_xyz(i,j,k) = SIN( zml * pi / model_height )  &
               &         / ( 1.0_wp + ( x / bub_radx(ii) )**2 )
          
        ENDDO
        
      ENDDO
    ENDDO
    
  END SUBROUTINE f_xyz_SK94

  SUBROUTINE f_xy_hotspot(ii, f_xy)

    IMPLICIT NONE

    ! Index of bubble in bubble namelist vectors:
    INTEGER(KIND=iintegers), INTENT(in)  :: ii
    ! spatial shape function of the bubble:
    REAL(KIND=wp),           INTENT(out) :: f_xy(ie,je)

    INTEGER(KIND=iintegers)      :: i, j, i_td, j_td
    REAL(KIND=wp)                :: zdx_rot, zdy_rot, zfak, ztmpwinkel

    f_xy(:,:) = 0.0_wp

    ztmpwinkel = bub_rotangle(ii) * degrad

    j_td  = isubpos(my_cart_id, 2) - nboundlines - 1
    DO j = 1, je
      j_td = j_td + 1
      i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
      DO i = 1, ie
        i_td = i_td + 1
        ! compute rectangular distances to the (rotated) main axes
        ! of the bubble (zdx_rot, zdy_rot):
        CALL hill_rot_coords(i_td, j_td, bub_centi(ii), bub_centj(ii), &
             ztmpwinkel, 0.0_wp, zdx_rot, zdy_rot)
        zfak  = EXP(-LOG(2.0_wp)*((zdx_rot/bub_radx(ii))**2 + (zdy_rot/bub_rady(ii))**2))
        IF (zfak >= 1.0E-2_wp) THEN
          f_xy(i,j) = zfak
        END IF
      ENDDO
    ENDDO

  END SUBROUTINE f_xy_hotspot
  
  ! Blase, die von Axel im 2005-er Hucmvergleichspaper verwendet wurde:
  ! Blasen haben bisher keine Drehmoeglichkeit!
  ! Heizung muss ueber 10 minuten mit einer Heizrate von 0.03 K/s gemacht werden!
  ! Die Hoehe der max. Heizung ist hier als auf einer Modellflaeche
  ! liegend definiert (bub_centk), so dass die max. Heizrate auch wirklich
  ! auftritt.

  ! T0_PLUME=0.03D0 ! K/s --> namelist-parameter "bub_heatingrate"
  ! X_PLUME=15.0D5  ! cm --> namelist-parameter "bub_centi"
  ! Y_PLUME=0.0D0   ! (im Original nicht vorhanden, da 2D-Simulation l2dim=.true.); namelist-parameter "bub_centj"
  ! Z_PLUME=0.5D5   ! cm --> Unterrand der Blase, ueber namelist-parameter bub_centk als k-index setzbar.
  ! DZ_PLUME=0.5D5  ! cm --> weiter unten als 0.002 / m fest verdrahtet

  SUBROUTINE f_xyz_AS2005(ii, f_xyz)

    IMPLICIT NONE

    ! Index of bubble in bubble namelist vectors:
    INTEGER(KIND=iintegers), INTENT(in)  :: ii
    ! spatial shape function of the bubble:
    REAL(KIND=wp),           INTENT(out) :: f_xyz(ie,je,ke)

    INTEGER(KIND=iintegers)      :: i, j, k, i_td, j_td
    REAL(KIND=wp)                :: zdx, zdy, zfak, zdi, zdj, ztmpwinkel, &
         zbubmax, zfak_horiz, zheight, &
         delx, zbub_zoomfakt


    f_xyz(:,:,:) = 0.0_wp

    ztmpwinkel = bub_rotangle(ii) * degrad

    ! Faktor fuer die Verbreiterung der PLUME, die Standardmaessig 5 km Durchmesser hat:
    ! (2.0 bedeutet einen Durchmesser von 10 km durch lineare Streckung)
    zbub_zoomfakt = 1.0E0_wp

    zdy   = r_earth* dlat * degrad
    j_td  = isubpos(my_cart_id, 2) - nboundlines - 1
!CDIR NOLOOPCHG
    DO j = 1, je
      zdx  = r_earth* dlon * degrad * crlat(j,1)
      j_td = j_td + 1
      zdj  = REAL(j_td-bub_centj(ii),wp)
      i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
!CDIR NOLOOPCHG
      DO i = 1, ie
        i_td = i_td + 1
        zdi  = REAL(i_td-bub_centi(ii),wp)
        ! Schicht mit der maximalen Heizung
        k = ke-bub_centk(ii)+1
        zbubmax = 0.5_wp*(hhl(i,j,k)+hhl(i,j,k+1))
        delx = SQRT((zdi*zdx)**2 +(zdj*zdy)**2) / zbub_zoomfakt
        IF(DELX.LE.0.0E0_wp)                      zfak_horiz=0.75E0_wp
        IF(DELX.GT.0.0E0_wp.AND.DELX.LE.0.5E3_wp) zfak_horiz=0.775E0_wp
        IF(DELX.GT.0.5E3_wp.AND.DELX.LE.1.0E3_wp) zfak_horiz=0.8E0_wp
        IF(DELX.GT.1.0E3_wp.AND.DELX.LE.1.5E3_wp) zfak_horiz=0.85E0_wp
        IF(DELX.GT.1.5E3_wp.AND.DELX.LE.2.0E3_wp) zfak_horiz=0.9E0_wp
        IF(DELX.GT.2.0E3_wp.AND.DELX.LE.2.5E3_wp) zfak_horiz=1.0E0_wp
        IF(DELX.GT.2.5E3_wp)                      zfak_horiz=0.0E0_wp
        
!CDIR NOLOOPCHG
        DO k = 1, ke
!!! disturbance is terrain following
          zheight = 0.5_wp*(hhl(i,j,k)+hhl(i,j,k+1))
          ! IF (zheight >= zbubmax) THEN
          IF (zheight >= zbubmax) THEN
            zfak = EXP(-0.002_wp*(zheight-zbubmax))
            zfak = MIN(zfak,1.0_wp)
            IF (zfak_horiz > 0.0_wp) THEN
              zfak = zfak * zfak_horiz
            ELSE
              zfak = 0.0_wp
            END IF
          ELSE
            zfak = 0.0_wp
          ENDIF
          f_xyz(i,j,k) = zfak
        ENDDO
      ENDDO
    ENDDO
    
  END SUBROUTINE f_xyz_AS2005

  !----------------------------------------------------------------------------
  ! Shape function for a 2D cos^2 temperature disturbance, intended
  ! for a soil temperature disturbance
  !----------------------------------------------------------------------------

  ! The cos^2-disturbance may be elliptic and may be horizontally rotated.

  SUBROUTINE f_xyz_cos2d_soil(ii, f_xyz)

    IMPLICIT NONE

    ! Index of bubble in bubble namelist vectors:
    INTEGER(KIND=iintegers), INTENT(in)  :: ii
    ! spatial shape function of the bubble:
    REAL(KIND=wp),           INTENT(out) :: f_xyz(ie,je,ke_soil)

    INTEGER(KIND=iintegers)      :: i, j, k, i_td, j_td
    REAL(KIND=wp)                :: zdx_rot, zdy_rot, zfak, ztmpwinkel

    f_xyz(:,:,:) = 0.0_wp

    ztmpwinkel = bub_rotangle(ii) * degrad

    j_td  = isubpos(my_cart_id, 2) - nboundlines - 1
!CDIR NOLOOPCHG
    DO j = 1, je
      j_td = j_td + 1
      i_td  = isubpos(my_cart_id, 1) - nboundlines - 1
!CDIR NOLOOPCHG
      DO i = 1, ie
        i_td = i_td + 1
        ! compute rectangular distances to the (rotated) main axes
        ! of the bubble (zdx_rot, zdy_rot):
        CALL hill_rot_coords(i_td, j_td, bub_centi(ii), bub_centj(ii), &
             ztmpwinkel, 0.0_wp, zdx_rot, zdy_rot)
!CDIR NOLOOPCHG
        DO k = 1, ke_soil
          zfak  = SQRT( (zdx_rot/bub_radx(ii))**2 + &
               (zdy_rot/bub_rady(ii))**2 )
          IF (zfak < 1.0_wp) THEN
            f_xyz(i,j,k) = COS(pi*zfak/2.0_wp)**2
          END IF
        ENDDO
      ENDDO
    ENDDO
    
  END SUBROUTINE f_xyz_cos2d_soil

!==============================================================================
!==============================================================================

  !==============================================================================
  !==============================================================================
  !
  ! Initialisation of the random number generator on a parallel machine
  ! with the following properties:
  !
  ! - if initialized with no "iseed_in" (optional integer parameter), then
  !   the resulting random number series (RNS) will be different for each
  !   model run, because seed is determined from the microseconds part of the actual date_and_time().
  !   Additionally, the PE-number is blended into the seed to ensure
  !   that the RNS is also different on each processor, even if called at the exact
  !   same time.
  !
  ! - if "iseed_in" is provided, this is used to generate the same seed on 
  !   a processor with a given PE-number each time the program runs, but
  !   again different seeds on PEs with different PE-numbers.
  !   This enables parallel generation of different random number series
  !   on each processor, which are however the same at each successive program run on
  !   each corresponding PE.
  !
  ! - NOTE: If you would like to, e.g., impose a random but reproducible noise on a 
  !   model field (i.e., INDEPENDENT of the number of PEs) which stays the same for each
  !   successive model run, then it is proposed that you calculate the noisy field
  !   globally on one processor and distribute it afterwards to the single nodes using
  !   the SR distribute_field() from parallel_utilities.f90.
  !   THIS IS DONE IN SR gen_bubnoise() BELOW !
  !
  !==============================================================================
  !==============================================================================

  !.. PGI-friendly version:
  SUBROUTINE seed_random_number(iseed_in) 
    
    IMPLICIT NONE 
    
    !.. local vars
    INTEGER(KIND=iintegers), OPTIONAL, INTENT(in) :: iseed_in

!!! UB>> Older settings with detailed type specifications seem not
!!!      to be necessary any more (aside from pgi-compiler, which
!!!      strangely seems to need integer*8 data types - 
!!!      see #ifdef __PGI_FORTRAN__ below)
!!$
!!$    INTEGER*4 :: i
!!$    INTEGER*4 :: k
!!$    INTEGER*4 :: i1
!!$#ifdef __SX__
!!$    INTEGER (KIND=iintegers) :: zeit(8), i2, iseed
!!$#else
!!$    INTEGER*8 :: zeit(8), i2, iseed
!!$#endif
!!$    INTEGER*4, ALLOCATABLE :: seed(:)

    INTEGER :: i
    INTEGER :: k
    INTEGER :: i1
#ifdef __PGI_FORTRAN__
    INTEGER*8 :: zeit(8), i2, iseed
#else
    INTEGER :: zeit(8), i2, iseed
#endif
    INTEGER, ALLOCATABLE :: seed(:)
    
    ! for pgi-compiler, the system_clock starts at 0 when system_clock is first called during a program.
    ! So by default, it measures a time *difference*, anticipating that the user
    ! only wants to time his program. So, at the first call, the returned time is always 0.
    ! Only for the subsequent calls, the time increases. What a nonsense!

    ! Unfortunately, this is unusable for the purpose of initializing the random number generator with
    ! a different seed for every program run, since this will all times lead
    ! to the same result, independent of a certain
    ! random or varying component.
    ! This is different from other compilers, where system_clock() delivers the elapsed time
    ! since 1.1.1970 in milliseconds, modulo HUGE(int).
    !
    ! So, we do it differently:
    ! First, we use HUGE() to determine the maxint value i2:
    i2 = HUGE(i1)
    IF (.NOT.PRESENT(iseed_in)) THEN
      ! and then we use DATE_AND_TIME() to get the milliseconds part of the actual time,
      ! which later will serve as the varying component from program run to program run:
      CALL DATE_AND_TIME(values=zeit)
      iseed=zeit(8)
    ELSE IF (iseed_in == -999_iintegers) THEN
      ! and then we use DATE_AND_TIME() to get the milliseconds part of the actual time,
      ! which later will serve as the varying component from program run to program run:
      CALL DATE_AND_TIME(values=zeit)
      iseed=zeit(8)
    ELSE
      ! or, if it is desired, we use just iseed_in, which leeds to the same random numbers 
      ! everytime:
      iseed = iseed_in
    END IF

    ! get length of seed vector:
    CALL RANDOM_SEED(SIZE=k)
    ALLOCATE( seed(k) )
    seed = 0

    ! However, in any case we want to have a different random number series on each task,
    ! so this is achieved by merging in my_cart_id into the seed.
    ! The seed itself is constructed in a way that it is (multiply) folded into
    ! the number range of integer data type, in order to break somehow the monotonicity
    ! in the seed vector. Monotonicity in the seed leads to a number series, from
    ! which the first 100 elements or so are not random but very close to 0, and only
    ! afterwards convert to more random behaviour.
    DO i=1,k
      seed(i) = INT(MOD(i2/11.0_wp*13.0_wp*((MOD(i,5)+i)*i) + i2*(iseed/1000.0_wp) + &
           my_cart_id*i2*0.95_wp/num_compute, i2*1.0_wp ), KIND=KIND(i2))
    END DO
    CALL RANDOM_SEED(PUT=seed)

    IF (k >= 4) THEN
      WRITE(*,'(a,i3,a,4(1x,i14))') '    SEED_RANDOM_NUMBER (first 4 of ',k,'):', seed(1:4)
    ELSE
      WRITE(*,*) '    SEED_RANDOM_NUMBER : ', seed
    END IF

    DEALLOCATE( seed )    
    RETURN
  END SUBROUTINE seed_random_number

  !.. old version, not working with PGI-Compiler:
  SUBROUTINE seed_random_number_old(iseed_in) 
    
    IMPLICIT NONE 
    
    !.. local vars
    INTEGER(KIND=iintegers), OPTIONAL, INTENT(in) :: iseed_in
    INTEGER(KIND=iintegers) :: i,k,i1,i2,iseed
    INTEGER(KIND=iintegers), ALLOCATABLE :: seed(:)
    REAL(KIND=wp)           :: r 
    
    IF (.NOT.PRESENT(iseed_in)) THEN
      CALL SYSTEM_CLOCK(iseed,i1,i2)
    ELSE
      iseed = iseed_in
    END IF

    CALL RANDOM_SEED(SIZE=k)
    ALLOCATE( seed(k) )
    DO i=1,k
      seed(i) = i*(iseed+1) + (my_cart_id+2)**2
    ENDDO
    CALL RANDOM_SEED(PUT=seed)
    
    IF (.NOT.PRESENT(iseed_in)) THEN
      CALL SYSTEM_CLOCK(iseed,i1,i2)
    ELSE
      iseed = iseed_in + 2
    END IF
    
    DO i=1,k
      CALL RANDOM_NUMBER(r)
      seed(i) = INT(r*(iseed+1)) + (my_cart_id+5)**2
    ENDDO
    CALL RANDOM_SEED(PUT=seed)
    WRITE(*,*) '    SEED_RANDOM_NUMBER: ', seed
    
    DEALLOCATE( seed )    
    RETURN      
  END SUBROUTINE seed_random_number_old


  !==============================================================================
  !==============================================================================
  !
  ! Fetch a REPRODUCIBLE 3D noise field (fnoise), given an (optional) seed or not,
  ! i.e., the noise is independent of the number of PEs. It is fetched
  ! only on one node an then divided to all other nodes using the SR 
  ! distribute_field() from parallel_utilities.f90.
  ! 
  ! This routine has to be called from all compute PEs in order to
  ! properly distribute the noise field.
  !
  ! The output field fnoise has shape (1:ie,1:je,kstart:kend).
  !
  !==============================================================================
  !==============================================================================

  SUBROUTINE gen_bubnoise(fnoise, kdim, seed)

    IMPLICIT NONE

    INTEGER(KIND=iintegers), INTENT(in)           :: kdim
    INTEGER(KIND=iintegers), INTENT(in), OPTIONAL :: seed
    REAL(KIND=wp),     INTENT(out)                :: fnoise(ie,je,1:kdim)

    REAL(KIND=wp),     ALLOCATABLE :: noisedummy(:), bub_glob_tnoisebuffer(:,:)

    INTEGER(KIND=iintegers) :: i, j, k, kzdims(24), izerror
    REAL(KIND=wp)           :: t_zufall(ie_tot)
    CHARACTER (LEN=250)     ::  yzerrmsg

    yzerrmsg(:) = ' '

    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 0) THEN
      WRITE (*,*)  'Subr. gen_bubnoise() ...'
    END IF

    IF (my_cart_id == 0) THEN
                  
      ! initialise the random number generator:
      IF (PRESENT(seed)) THEN
        CALL seed_random_number( seed )
      ELSE
        CALL seed_random_number( )
      END IF

      ! on some compilers, the random series shows problems within the first few
      ! hundred random numbers (the numbers are not really random, but
      ! can be monotonic and reproducibly the same on all processors).
      ! Only after some numbers, the series gets more random type.
      ! Therefore, fetch 5000 dummy random numbers, before doing the
      ! really needed random numbers:
                  
      ALLOCATE(noisedummy(5000))
      CALL RANDOM_NUMBER(noisedummy)
      DEALLOCATE(noisedummy)
                  
      ! allocate temporary global noise field 2D:
      ALLOCATE(bub_glob_tnoisebuffer(ie_tot, je_tot))
      bub_glob_tnoisebuffer = 0.0_wp
                  
    ELSE

      ! dummy allocation, necessary for distribute_field() below at least for gfortran and OpenMPI:
      ALLOCATE(bub_glob_tnoisebuffer(1, 1))
      bub_glob_tnoisebuffer = 0.0_wp

    END IF
                
    ! Fetch white noise in the Interval [-1;1] to the global noise field.
    DO k = 1, kdim
                  
      IF (my_cart_id == 0) THEN
        DO j= 1, je_tot
          CALL RANDOM_NUMBER(t_zufall(:))
          DO i= 1, ie_tot
            bub_glob_tnoisebuffer(i,j) = 2.0_wp * (t_zufall(i) - 0.5_wp)
          END DO
        END DO
      END IF
      
      ! distribute global noise field to all subdomains:
      IF ( num_compute > 1 ) THEN
        CALL distribute_field(bub_glob_tnoisebuffer(:,:),ie_tot,je_tot, &
             fnoise(:,:,k),ie,je,0,izerror)
        IF (izerror /= 0) THEN
          WRITE (*,*) 'ERROR: in gen_bubnoise(), src_artifdata.f90: distribution of global noise field'
          CALL model_abort (my_world_id, 10002+izerror, 'distribution of noise field failed!', &
               'gen_bubnoise(), src_artifdata.f90, distribution of global noise field')
        END IF
      ELSE
        fnoise(:,:,k) = bub_glob_tnoisebuffer(:,:)
      END IF
      
    END DO
    
    ! clean up memory:
    IF (my_cart_id == 0) DEALLOCATE(bub_glob_tnoisebuffer)

    ! exchange noise field for bubble number ii:
    kzdims(1:24)=(/kdim,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL exchg_boundaries                                                &
         (0, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,           &
         kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
         lperi_x, lperi_y, l2dim, &
         20000, .FALSE., ncomm_type, izerror, yzerrmsg,              &
         fnoise(:,:,1:kdim) )

    IF (izerror /= 0) THEN
      WRITE (*,*) 'ERROR: in gen_bubnoise(), src_artifdata.f90: exchg_boundaries'
      WRITE (*,*) 'ERROR: ', TRIM(yzerrmsg)
      CALL model_abort (my_world_id, 10003+izerror, 'exchg_boundaries failed!', &
           'gen_bubnoise(), src_artifdata.f90, exchg_boundaries')
    END IF

  END SUBROUTINE gen_bubnoise

!=====================================================================================
!=====================================================================================
!
! Iteratively calculate pressure piter on the levels zml, given surface
! height hsurf, surface pressure psurf, also given the temperature t 
! (or potential temp. theta, if present), relative humidity relhum,  
! base state pressure p0 and base state density rho0, each on the levels zml.
!
! qv resp. qc is diagnosed from t / theta and the relative humidity in a way 
! as to condense all excess water at supersaturated voxels to qc.
!
! This routine requires monotonically decreasing height levels and distinguishes
! between orography height hsurf and model levels zml above it. It is especially for
! the pressure integration at the end of gen_ini_data().
!
! t and relhum at the surface (hsurf) are assumed to be the same values
! as in the layer above (zml(:,:,nk)).
!
!!! The SR calc_p_hydrostat_psts() leads
!!! to a solution which is exactly compatible with the
!!! p'T'-dynamics of the model.
!
! Changes:
!
!  - 6.8.2012: Added support for new fast waves solver (itype_fastwaves=2) (UB)
!
!==============================================================================
!==============================================================================

  SUBROUTINE calc_p_hydrostat_psts(ni,nj,nk,niter, &
       zml, sqrtg, hsurf, psurf, t, relhum, t0, p0, rho0, &
       piter, qv, qc, zmaxqv, r_d, rvd_m_o, itype_fastwaves, &
       ierror, errmsg, &
       t0hl, p0hl, wgtfac, theta)

    IMPLICIT NONE

    !.. Input/output variables:
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk, niter, itype_fastwaves
    REAL(KIND=wp),     INTENT(in) :: zml(ni,nj,nk), sqrtg(ni,nj,nk), hsurf(ni,nj), psurf(ni,nj), &
         relhum(ni,nj,nk)
    REAL(KIND=wp),     INTENT(out) :: qv(ni,nj,nk), qc(ni,nj,nk), piter(ni,nj,nk)
    !.. p, rho of the reference atmosphere at full levels:
    REAL(KIND=wp),     INTENT(in) :: t0(ni,nj,nk), p0(ni,nj,nk), rho0(ni,nj,nk)
    !.. at half levels needed only for itype_fastwaves=2:
    REAL(KIND=wp),     INTENT(in) :: zmaxqv, r_d, rvd_m_o
    REAL(KIND=wp),     INTENT(inout) :: t(ni,nj,nk)
    INTEGER(kind=iintegers), INTENT(out) :: ierror
    CHARACTER(len=250), INTENT(out) :: errmsg
    REAL(KIND=wp),     INTENT(in), OPTIONAL :: t0hl(ni,nj,nk), p0hl(ni,nj,nk), wgtfac(ni,nj,nk)
    REAL(KIND=wp),     INTENT(in), OPTIONAL :: theta(ni,nj,nk)

    !.. Local variables:
    REAL(KIND=wp)     :: piterold(ni,nj,nk), relhum_lim(ni,nj,nk)
    REAL(KIND=wp)     :: eps_pp, zpa, ztvw, ztvw1, zesat, zsqv, zrhogdh, zrdm, zcpm
    INTEGER(KIND=iintegers) ::  i, j, k, iter
    LOGICAL :: ltheta_ini

    !.. For itype_fastwaves == 2:
    REAL(KIND=wp)     :: one_m_wgtfac, Tp_k_avg, pp_k_avg, T_hl, p_hl, rho_inv_k_avg, &
                         p0_p_k_avg, T_T0_k_avg, p_inv_k_avg, T0_inv_k_avg, q_x_k_avg, buoy


    ierror = 0
    errmsg(:) = ' '

    eps_pp = 1.0E-4_wp  ! required absolute accuracy of the pressure iteration, [Pa]

    !.. Check input parameters:

    IF (itype_fastwaves < 1 .OR. itype_fastwaves > 2) THEN
      ierror = 10010
      errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_psts(): Wrong itype_fastwaves!'
      RETURN
    END IF

    IF (itype_fastwaves == 2) THEN
      IF (.NOT.PRESENT(t0hl)) THEN
        ierror = 10011
        errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_psts(): '// &
                 'itype_fastwaves == 2, but t0hl not present in input parameters!'
        RETURN
      END IF
      IF (.NOT.PRESENT(p0hl)) THEN
        ierror = 10011
        errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_psts(): '// &
                 'itype_fastwaves == 2, but p0hl not present in input parameters!'
        RETURN
      END IF
      IF (.NOT.PRESENT(wgtfac)) THEN
        ierror = 10011
        errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_psts(): '// &
                 'itype_fastwaves == 2, but wgtfac not present in input parameters!'   
        RETURN
      END IF
    END IF

    IF (PRESENT(theta)) THEN
      ltheta_ini = .TRUE.
      ! Initialize qv and qc with their "dry" values so that zrdm and zcpm
      ! have a defined value in the first iteration:
      qv(:,:,:) = 0.0_wp
      qc(:,:,:) = 0.0_wp
    ELSE
      ltheta_ini = .FALSE.
      ! Initialize qv and qc with their "dry" values so that qv_Tprelhum()
      ! gets defined values in the first "moist" iteration:
      qv(:,:,:) = 0.0_wp
      qc(:,:,:) = 0.0_wp
    END IF


    ! Instead of iteration until convergence into some max. absolute residuum
    ! we choose to do a fixed number of iterations, since this produces
    ! reproducible results on vector machines.
    ! However, after the iteration it is checked whether the iteration
    ! has converged to the specified absolute error bounds resp.
    ! iteration increment on all gridpoints.

    ! initial pp at the surface, will be iteratively corrected below:
    piter(:,:,:)  = p0(:,:,:)
    piter(:,:,nk) = psurf(:,:)

    DO iter = 1, niter

      ! New variable for limited relhum to max. allowed value p / E(T):
      ! Re-initialize it to relhum for each iteration
      ! to avoid drifts of relhum_lim during the iteration:
      relhum_lim(:,:,:) = relhum(:,:,:)

      ! Store the initial pp at the beginning of the iteration step:
      piterold(:,:,:) = piter(:,:,:)

      ! pressure deviation on the lowest full level
      ! and other initializations
      DO j = 1,nj
        DO i = 1,ni
          ! decomposition of the relative humidity into qv, qc using the 
          ! momentary pressure zpa = p0+pp, depending on the settings of lcond:
          zpa          = piter(i,j,nk)
          IF (ltheta_ini) THEN
            ! "moist" r_d:
            zrdm      = rd_moist(qv(i,j,nk),qc(i,j,nk))
            ! COSMO-approximation of cp:
            zcpm      = cp_moist_cosmo(qv(i,j,nk),qc(i,j,nk),0.0_wp)
            t(i,j,nk) = theta(i,j,nk)*(zpa/pt00)**(zrdm/zcpm)
          END IF
!!$          IF (t(i,j,nk) > 233.16) THEN
            zesat       = esat_w(t(i,j,nk))
!!$          ELSE
!!$            zesat       = esat_i(t(i,j,nk))
!!$          ENDIF
          ! Limit relhum to its maximum possible value p / E(T):
          relhum_lim(i,j,nk) = MIN( zpa/zesat , relhum(i,j,nk) )
          ! Compute actual Qv:
          qv(i,j,nk) = qv_Tprelhum( zpa, t(i,j,nk), relhum_lim(i,j,nk), qc(i,j,nk) )
          IF (lcond .AND. relhum_lim(i,j,1) > 1.0_wp) THEN
            ! condensation is allowed and physically can happen
            ! and relhum > 1.0, so convert qv -> qc to limit relhum to 1.0:
            zsqv = qvsat_w( zpa, t(i,j,nk) )
            qv(i,j,nk) = MIN ( zsqv       ,  MIN(qv(i,j,nk), zmaxqv) )
            qc(i,j,nk) = MAX ( 0.0_wp ,  MIN(qv(i,j,nk), zmaxqv) - zsqv )
          ELSE
            ! condensation is not allowed or cannot happen physically at 
            ! that pressure and temperature, so just impose the limit zmaxqv:
            qv(i,j,nk) = MIN(qv(i,j,nk), zmaxqv)
            qc(i,j,nk) = 0.0_wp
          END IF

          ! pressure on the full level nk by isothermal extrapolation from the ground:
          ztvw         = t(i,j,nk) &
               * (1.0_wp + rvd_m_o*qv(i,j,nk) - qc(i,j,nk))
          piter(i,j,nk) = psurf(i,j) * EXP((hsurf(i,j)-zml(i,j,nk))*g/(ztvw*r_d))            

        ENDDO
      ENDDO

      ! pressure deviation on the full level k-1
      DO k = nk,2, -1

        DO j = 1,nj
          DO i = 1,ni

            ! decomposition of the relative humidity into qv, qc using the 
            ! momentary pressure zpa = p0+pp
            zpa           = piter(i,j,k-1)
            IF (ltheta_ini) THEN
              ! "moist" r_d:
              zrdm       = rd_moist(qv(i,j,k-1),qc(i,j,k-1))
              ! COSMO-approximation of cp:
              zcpm       = cp_moist_cosmo(qv(i,j,k-1),qc(i,j,k-1),0.0_wp)
              t(i,j,k-1) = theta(i,j,k-1)*(zpa/pt00)**(zrdm/zcpm)
            END IF
!!$            IF (t(i,j,k-1) > 233.16) THEN
              zesat       = esat_w(t(i,j,k-1))
!!$            ELSE
!!$              zesat       = esat_i(t(i,j,k-1))
!!$            ENDIF
            ! Limit relhum to its maximum possible value p / E(T):
            relhum_lim(i,j,k-1) = MIN( zpa/zesat , relhum(i,j,k-1) )
            ! Compute actual Qv:
            qv(i,j,k-1) = qv_Tprelhum( zpa, t(i,j,k-1), relhum_lim(i,j,k-1), qc(i,j,k-1) )
            IF (lcond .AND. relhum_lim(i,j,k-1) > 1.0_wp) THEN
              ! condensation is allowed and physically can happen
              ! and relhum > 1.0, so convert qv -> qc to limit relhum to 1.0:
              zsqv = qvsat_w( zpa, t(i,j,k-1) )
              qv(i,j,k-1) = MIN ( zsqv       ,  MIN(qv(i,j,k-1), zmaxqv) )
              qc(i,j,k-1) = MAX ( 0.0_wp ,  MIN(qv(i,j,k-1), zmaxqv) - zsqv )
            ELSE
              ! condensation is not allowed, just set the limit to zmaxqv:
              qv(i,j,k-1) = MIN(qv(i,j,k-1), zmaxqv)
              qc(i,j,k-1) = 0.0_wp
            END IF

            IF (itype_fastwaves == 1) THEN

              ! virtual temperature at levels k and k-1, so that we can use "dry" R_d instead of rd_moist:
              ztvw1 = t(i,j,k-1)*(1.0_wp + rvd_m_o*qv(i,j,k-1) - qc(i,j,k-1))
              ztvw  = t(i,j,k  ) * (1.0_wp + rvd_m_o*qv(i,j,k) - qc(i,j,k))

              zrhogdh = 0.25_wp / sqrtg(i,j,k) * g * (rho0(i,j,k-1)+rho0(i,j,k))

              ! This expression for the pressure has been derived by solving the
              ! discretized r.h.s. of the vertical equation of motion for pp.
              ! Provides exact hydrostatic balance for p'T'-dynamics (itheta_adv=0)
              zpa   = piter(i,j,k) - p0(i,j,k)
              piter(i,j,k-1) = p0(i,j,k-1) + ( zpa + zrhogdh*(2._wp-t0(i,j,k-1)/ztvw1 &
                   - t0(i,j,k)/ztvw*(1._wp + zpa/p0(i,j,k)))) &
                   / (1._wp + zrhogdh*t0(i,j,k-1)/(ztvw1*p0(i,j,k-1)))

            ELSE IF (itype_fastwaves == 2) THEN

              one_m_wgtfac = 1.0_wp - wgtfac(i,j,k)

              Tp_k_avg = wgtfac(i,j,k) * ( t(i,j,k)   - t0(i,j,k)   )  +  &
                   &      one_m_wgtfac * ( t(i,j,k-1) - t0(i,j,k-1) )

              pp_k_avg = wgtfac(i,j,k) * ( piter(i,j,k)   - p0(i,j,k)   )  +  &
                   &      one_m_wgtfac * ( piter(i,j,k-1) - p0(i,j,k-1) )

              T_hl = t0hl(i,j,k) + Tp_k_avg
              p_hl = p0hl(i,j,k) + pp_k_avg

              zrdm       = rd_moist(qv(i,j,k-1),qc(i,j,k-1))
              rho_inv_k_avg   = ( zrdm * T_hl ) / p_hl

              p0_p_k_avg   = 1.0_wp / ( 1.0_wp + pp_k_avg / p0hl(i,j,k) )
              T_T0_k_avg   = 1.0_wp + Tp_k_avg / t0hl(i,j,k)
              p_inv_k_avg  = 1.0_wp / ( p0hl(i,j,k) + pp_k_avg )
              T0_inv_k_avg = 1.0_wp / t0hl(i,j,k)

              q_x_k_avg = wgtfac(i,j,k) * ( rvd_m_o*qv(i,j,k)   - qc(i,j,k)   )  +  &
                           one_m_wgtfac * ( rvd_m_o*qv(i,j,k-1) - qc(i,j,k-1) )

              buoy =  g * ( p0_p_k_avg * T0_inv_k_avg * Tp_k_avg         &
                            - p_inv_k_avg * pp_k_avg                     &
                            + p0_p_k_avg * T_T0_k_avg * q_x_k_avg  )

              zpa = p0(i,j,k) - p0(i,j,k-1)

              piter(i,j,k-1) = piter(i,j,k) - zpa + buoy / ( rho_inv_k_avg * sqrtg(i,j,k) )

            END IF

          ENDDO
        ENDDO

      ENDDO

      ! Catch a badly diverged iteration, which can happen in case
      !   of initialization with very stable N-const profiles (itype_anaprof_tqv = 3):
      IF ( ANY(piter < 0.0_wp .OR. t < 0.0_wp .OR. t > 1000.0_wp) ) THEN
        !.. For meaningful error message, print first erroneous point:
        DO j = 1,nj
          DO i = 1,ni
            IF ( ANY( piter(i,j,:) < 0.0_wp .OR. t(i,j,:) < 0.0_wp .OR. t(i,j,:) > 1000.0_wp ) ) THEN

              WRITE (*,'(a)') 'ERROR: Subroutine calc_p_hydrostat_psts() failed during pressure iteration'
              WRITE (*,'(a,i4,1x,i4)') '    at location (i,j) ', &
                   isubpos(my_cart_id, 1)-nboundlines-1+i, &
                   isubpos(my_cart_id, 2)-nboundlines-1+j
              IF (ltheta_ini) THEN
                WRITE (*,'(a)') '   k   zml(k) theta(k)    t(k) relhum(k)         p(k)        '// &
                     'qv(k)        qc(k)'
                DO k=1,nk
                  WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5))') &
                       k,zml(i,j,k),theta(i,j,k),t(i,j,k),relhum(i,j,k),piter(i,j,k),qv(i,j,k),qc(i,j,k)
                END DO
                ierror = 10081
                errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_psts(): Problem in iteration of '// &
                     'hydrostatic initialisation of pressure field! '// &
                     'Check namelist parameters for stability (N) and/or try to reduce the model top height!'
              ELSE
                WRITE (*,'(a)') '   k   zml(k)    t(k) relhum(k)           p(k)        '// &
                     'qv(k)        qc(k)'
                DO k=1,nk
                  WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5))') &
                       k,zml(i,j,k),t(i,j,k),relhum(i,j,k),piter(i,j,k),qv(i,j,k),qc(i,j,k)
                END DO
                ierror = 10082
                errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_psts(): Problem in iteration of '// &
                     'hydrostatic initialisation of pressure field! '// &
                     'Initial T profile seems to be problematic!'
              END IF

              RETURN

            END IF
          END DO
        END DO
      END IF

    ENDDO   ! of the iteration

    !.. If necessary, re-calculate T from Theta and iterated pressure:
    IF (ltheta_ini) THEN
      DO k = 1,nk
        DO j = 1,nj
          DO i = 1,ni
            ! "moist" r_d:
            zrdm     = rd_moist(qv(i,j,k),qc(i,j,k))
            ! COSMO-approximation of cp:
            zcpm     = cp_moist_cosmo(qv(i,j,k),qc(i,j,k),0.0_wp)
            t(i,j,k) = theta(i,j,k)*(piter(i,j,k)/pt00)**(zrdm/zcpm)
          END DO
        END DO
      END DO
    END IF

    !.. Debug output of the iterated pressure profile:
    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 4) THEN
      IF (ltheta_ini) THEN
        WRITE (*,'(a)') 'Subroutine calc_p_hydrostat_psts() after '// &
             'pressure iteration at point (1,1):'
        WRITE (*,'(a)') 'k     zml(k) theta(k) relhum(k)     t(k)  piter(k)      '// &
             'qv(k)        qc(k) relhum_out(k)'
        DO k=1,nk
          WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f8.1,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
               k,zml(1,1,k),theta(1,1,k),relhum(1,1,k),t(1,1,k), &
               piter(1,1,k),qv(1,1,k),qc(1,1,k), &
               rh_Tpqv(piter(1,1,k),t(1,1,k),qv(1,1,k),qc(1,1,k))
        END DO

      ELSE
        WRITE (*,'(a)') 'Subroutine calc_p_hydrostat_psts() after '// &
             'pressure iteration at point (1,1):'
        WRITE (*,'(a)') 'k     zml(k)     t(k) relhum(k)     p(k)      '// &
             'qv(k)        qc(k) relhum_out(k)'
        DO k=1,nk
          WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
               k,zml(1,1,k),t(1,1,k),relhum(1,1,k),piter(1,1,k),qv(1,1,k),qc(1,1,k), &
               rh_Tpqv(piter(1,1,k),t(1,1,k),qv(1,1,k),qc(1,1,k))
        END DO
      END IF
    END IF

    ! Convergence check: 
    DO j = 1,nj
      DO i = 1,ni
        IF ( ANY( ABS(piter(i,j,:) - piterold(i,j,:)) > eps_pp ) ) THEN
          WRITE (*,'(a)') 'ERROR: in calc_p_hydrostat_psts(), src_artifdata.f90: '//&
               'no convergence in iteration of hydrostatic initialisation!'
          WRITE (*,'(a,i4,1x,i4,a,f12.4,a,f14.7)') 'at location (i,j) ', &
               isubpos(my_cart_id, 1)-nboundlines-1+i, &
               isubpos(my_cart_id, 2)-nboundlines-1+j, &
               '    psurf(i,j) = ', psurf(i,j), '   dp_last_max = ', &
               MAXVAL(ABS(piterold(i,j,:)-piter(i,j,:)))
          IF (ltheta_ini) THEN
            WRITE (*,'(a)') 'k     zml(k) theta(k) relhum(k)     t(k) piter(k)      '// &
                 'qv(k)        qc(k) relhum_out(k)'
            DO k=1,nk
              WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f8.1,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
                   k,zml(i,j,k),theta(i,j,k),relhum(i,j,k),t(i,j,k),&
                   piter(i,j,k),qv(i,j,k),qc(i,j,k), &
                   rh_Tpqv(piter(i,j,k),t(i,j,k),qv(i,j,k),qc(i,j,k))
            END DO
          ELSE
            WRITE (*,'(a)') 'k     zml(k)     t(k) relhum(k) piter(k)      '// &
                 'qv(k)        qc(k) relhum_out(k)'
            DO k=1,nk
              WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
                   k,zml(i,j,k),t(i,j,k),relhum(i,j,k),piter(i,j,k),qv(i,j,k),qc(i,j,k), &
                   rh_Tpqv(piter(i,j,k),t(i,j,k),qv(i,j,k),qc(i,j,k))
            END DO
          END IF
          ierror = 10008
          errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_psts(): No convergence in iteration of '// &
                   'hydrostatic initialisation of pressure field!'
          RETURN
        END IF
      END DO
    END DO

    !.. Issue a warning message if the relat. humid. has been altered above:
    IF (ANY(relhum_lim(:,:,:) /= relhum(:,:,:))) THEN
      WRITE (*,'(a)') REPEAT('*',70)
      WRITE (*,'(a)') 'WARNING: calc_p_hydrostat_psts: Rel. humid. '// &
           'limited to max. allowed VALUE p / E(T)'
      WRITE (*,'(a)') REPEAT('*',70)
    END IF

  END SUBROUTINE calc_p_hydrostat_psts

!==============================================================================
!==============================================================================
!
! Subroutine similar to calc_p_hydrostat_psts() above, 
! but calc_p_hydrostat_lf() leads
! to a solution which is exactly compatible with the
! Leapfrog-dynamics of the model instead of RK dynamics.
!
!==============================================================================
!==============================================================================

  SUBROUTINE calc_p_hydrostat_lf(ni,nj,nk,niter, &
       zml, hsurf, psurf, t, relhum, t0, p0, dp0, rho0, &
       piter, qv, qc, zmaxqv, r_d, rvd_m_o, ierror, errmsg, &
       theta)

    IMPLICIT NONE

    !.. Input/output variables:
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk, niter
    REAL(KIND=wp),     INTENT(in) :: zml(ni,nj,nk), hsurf(ni,nj), psurf(ni,nj), &
         relhum(ni,nj,nk), t0(ni,nj,nk), p0(ni,nj,nk), dp0(ni,nj,nk), rho0(ni,nj,nk)
    REAL(KIND=wp),     INTENT(out) :: qv(ni,nj,nk), qc(ni,nj,nk), piter(ni,nj,nk)
    REAL(KIND=wp),     INTENT(in) :: zmaxqv, r_d, rvd_m_o
    REAL(KIND=wp),     INTENT(inout) :: t(ni,nj,nk)
    INTEGER(kind=iintegers), INTENT(out) :: ierror
    CHARACTER(len=250), INTENT(out) :: errmsg
    REAL(KIND=wp),     INTENT(in), OPTIONAL :: theta(ni,nj,nk)

    !.. Local variables:
    REAL(KIND=wp)     :: piterold(ni,nj,nk), relhum_lim(ni,nj,nk), &
         ztvdt(ni,nj,2), zt0dp0t(ni,nj,2)
    REAL(KIND=wp)     :: eps_pp, zpa, ztvw, zesat, zsqv, zleft, zrhs, zrdm, zcpm
    INTEGER(KIND=iintegers) ::  i, j, k, iter, klow, kup
    LOGICAL :: ltheta_ini


    ierror = 0
    errmsg(:) = ' '

    IF (PRESENT(theta)) THEN
      ltheta_ini = .TRUE.
      ! Initialize qv and qc with their "dry" values so that zrdm and zcpm
      ! have a defined value in the first iteration:
      qv(:,:,:) = 0.0_wp
      qc(:,:,:) = 0.0_wp
    ELSE
      ltheta_ini = .FALSE.
      ! Initialize qv and qc with their "dry" values so that qv_Tprelhum()
      ! gets defined values in the first "moist" iteration:
      qv(:,:,:) = 0.0_wp
      qc(:,:,:) = 0.0_wp
    END IF

    ! Instead of iteration until convergence into some max. absolute residuum
    ! we choose to do a fixed number of iterations, since this produces
    ! reproducible results on vector machines.
    ! However, after the iteration it is checked whether the iteration
    ! has converged to the specified absolute error bounds resp.
    ! iteration increment on all gridpoints.

    eps_pp = 1.0E-4_wp  ! required absolute accuracy of the pressure iteration, [Pa]

    ! initial pp at the surface, will be iteratively corrected below:
    piter(:,:,:)  = p0(:,:,:)
    piter(:,:,nk) = psurf(:,:)
    

    DO iter = 1, niter

      ! organizational indices
      kup  = 2
      klow = 1

      ! Store the initial pp at the beginning of the iteration step:
      piterold(:,:,:) = piter(:,:,:)

      ! New variable for limited relhum to max. allowed value p / E(T):
      ! Re-initialize it to relhum for each iteration
      ! to avoid drifts of relhum_lim during the iteration:
      relhum_lim(:,:,:) = relhum(:,:,:)

      ! pressure deviation on the lowest full level
      ! and other initializations
      DO j = 1,nj
        DO i = 1,ni
          ! decomposition of the relative humidity into qv, qc using the 
          ! momentary pressure zpa = p0+pp, depending on the settings of lcond:
          zpa          = piter(i,j,nk)
          IF (ltheta_ini) THEN
            ! "moist" r_d:
            zrdm      = rd_moist(qv(i,j,nk),qc(i,j,nk))
            ! COSMO-approximation of cp
            zcpm      = cp_moist_cosmo(qv(i,j,nk),qc(i,j,nk),0.0_wp)
            t(i,j,nk) = theta(i,j,nk)*(zpa/pt00)**(zrdm/zcpm)
          END IF
!!$          IF (t(i,j,nk) > 233.16) THEN
            zesat       = esat_w(t(i,j,nk))
!!$          ELSE
!!$            zesat       = esat_i(t(i,j,nk))
!!$          ENDIF
          ! Limit relhum to its maximum possible value p / E(T):
          relhum_lim(i,j,nk) = MIN( zpa/zesat , relhum(i,j,nk) )
          ! Compute actual Qv:
          qv(i,j,nk) = qv_Tprelhum( zpa, t(i,j,nk), relhum_lim(i,j,nk), qc(i,j,nk) )
          IF (lcond .AND. relhum_lim(i,j,1) > 1.0_wp) THEN
            ! condensation is allowed and physically can happen
            ! and relhum > 1.0, so convert qv -> qc to limit relhum to 1.0:
            zsqv = qvsat_w( zpa, t(i,j,nk) )
            qv(i,j,nk) = MIN ( zsqv       ,  MIN(qv(i,j,nk), zmaxqv) )
            qc(i,j,nk) = MAX ( 0.0_wp ,  MIN(qv(i,j,nk), zmaxqv) - zsqv )
          ELSE
            ! condensation is not allowed or cannot happen physically at 
            ! that pressure and temperature, so just impose the limit zmaxqv:
            qv(i,j,nk) = MIN(qv(i,j,nk), zmaxqv)
            qc(i,j,nk) = 0.0_wp
          END IF

          ! pressure on the full level nk by isothermal extrapolation from the ground:
          ztvw         = t(i,j,nk) &
               * (1.0_wp + rvd_m_o*qv(i,j,nk) - qc(i,j,nk))
          piter(i,j,nk) = psurf(i,j) * EXP((hsurf(i,j)-zml(i,j,nk))*g/(ztvw*r_d))

          ! contribution of the virtual temperature to the buoyancy term
          ztvdt(i,j,klow)   = (ztvw - t0(i,j,nk)) / t(i,j,nk)
          
          ! coefficient of the pressure contribution to the buoyancy term
          zt0dp0t(i,j,klow) = 0.5_wp / ( r_d * rho0(i,j,nk) * t(i,j,nk) )

        ENDDO
      ENDDO

      ! pressure deviation on the full level k-1
      DO k = nk,2, -1

        DO j = 1,nj
          DO i = 1,ni

            ! decomposition of the relative humidity into qv, qc using the 
            ! momentary pressure zpa = p0+pp
            zpa           = piter(i,j,k-1)
            IF (ltheta_ini) THEN
              ! "moist" r_d:
              zrdm       = rd_moist(qv(i,j,k-1),qc(i,j,k-1))
              ! COSMO-approximation of cp:
              zcpm       = cp_moist_cosmo(qv(i,j,k-1),qc(i,j,k-1),0.0_wp)
              t(i,j,k-1) = theta(i,j,k-1)*(zpa/pt00)**(zrdm/zcpm)
            END IF
!!$            IF (t(i,j,k-1) > 233.16) THEN
              zesat       = esat_w(t(i,j,k-1))
!!$            ELSE
!!$              zesat       = esat_i(t(i,j,k-1))
!!$            ENDIF
            ! Limit relhum to its maximum possible value p / E(T):
            relhum_lim(i,j,k-1) = MIN( zpa/zesat , relhum(i,j,k-1) )
            ! Compute actual Qv:
            qv(i,j,k-1) = qv_Tprelhum( zpa, t(i,j,k-1), relhum_lim(i,j,k-1), qc(i,j,k-1) )
            IF (lcond .AND. relhum_lim(i,j,k-1) > 1.0_wp) THEN
              ! condensation is allowed and physically can happen
              ! and relhum > 1.0, so convert qv -> qc to limit relhum to 1.0:
              zsqv = qvsat_w( zpa, t(i,j,k-1) )
              qv(i,j,k-1) = MIN ( zsqv       ,  MIN(qv(i,j,k-1), zmaxqv) )
              qc(i,j,k-1) = MAX ( 0.0_wp ,  MIN(qv(i,j,k-1), zmaxqv) - zsqv )
            ELSE
              ! condensation is not allowed, just set the limit to zmaxqv:
              qv(i,j,k-1) = MIN(qv(i,j,k-1), zmaxqv)
              qc(i,j,k-1) = 0.0_wp
            END IF

            ! virtual temperature at levels k and k-1
            ztvw = t(i,j,k-1)*(1.0_wp + rvd_m_o*qv(i,j,k-1) - qc(i,j,k-1))
            ztvdt(i,j,kup) = (ztvw - t0(i,j,k-1)) / t(i,j,k-1)

            ! coefficient of the pressure contribution to the buoyancy term
            zt0dp0t(i,j,kup) = 0.5_wp / ( r_d * rho0(i,j,k-1) * t(i,j,k-1) )

            ! pressure deviation on the full level k-1
            zpa   = piter(i,j,k) - p0(i,j,k)
            zleft         = 1.0_wp + zt0dp0t(i,j,kup)*dp0(i,j,k)
            zrhs  = zpa*( 1.0_wp - zt0dp0t(i,j,klow)*dp0(i,j,k-1) ) &
                 + 0.5_wp*( ztvdt(i,j,klow)*dp0(i,j,k-1)  &
                        +ztvdt(i,j,kup )*dp0(i,j,k  ) )
            piter(i,j,k-1) = p0(i,j,k-1) + zrhs / zleft

          ENDDO
        ENDDO

        ! changing the organizational indices
        klow = 3 - klow
        kup  = 3 - kup

      ENDDO

      ! Catch a badly diverged iteration, which can easily happen in case
      !   of initialization with N-const profiles (itype_anaprof_tqv = 3):
      IF ( ANY(piter < 0.0_wp .OR. t < 0.0_wp .OR. t > 1000.0_wp) ) THEN
        !.. For meaningful error message, print first erroneous point:
        DO j = 1,nj
          DO i = 1,ni
            IF ( ANY( piter(i,j,:) < 0.0_wp .OR. t(i,j,:) < 0.0_wp .OR. t(i,j,:) > 1000.0_wp ) ) THEN

              WRITE (*,'(a)') 'ERROR: Subroutine calc_p_hydrostat_lf() failed during pressure iteration'
              WRITE (*,'(a,i4,1x,i4)') '    at location (i,j) ', &
                   isubpos(my_cart_id, 1)-nboundlines-1+i, &
                   isubpos(my_cart_id, 2)-nboundlines-1+j
              IF (ltheta_ini) THEN
                WRITE (*,'(a)') '   k   zml(k) theta(k)    t(k) relhum(k)         p(k)        '// &
                     'qv(k)        qc(k)'
                DO k=1,nk
                  WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5))') &
                       k,zml(i,j,k),theta(i,j,k),t(i,j,k),relhum(i,j,k),piter(i,j,k),qv(i,j,k),qc(i,j,k)
                END DO
                ierror = 10083
                errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_lf(): Problem in iteration of '// &
                     'hydrostatic initialisation of pressure field! '// &
                     'Check namelist parameters for stability (N) and/or try to reduce the model top height!'
              ELSE
                WRITE (*,'(a)') '   k   zml(k)    t(k) relhum(k)           p(k)        '// &
                     'qv(k)        qc(k)'
                DO k=1,nk
                  WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5))') &
                       k,zml(i,j,k),t(i,j,k),relhum(i,j,k),piter(i,j,k),qv(i,j,k),qc(i,j,k)
                END DO
                ierror = 10084
                errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_lf(): Problem in iteration of '// &
                     'hydrostatic initialisation of pressure field! '// &
                     'Initial T profile seems to be problematic!'
              END IF

              RETURN

            END IF
          END DO
        END DO
      END IF

    ENDDO   ! of the iteration

    !.. If necessary, re-calculate T from Theta and iterated pressure:
    IF (ltheta_ini) THEN
      DO k = 1,nk
        DO j = 1,nj
          DO i = 1,ni
            ! "moist" r_d:
            zrdm     = rd_moist(qv(i,j,k),qc(i,j,k))
            ! COSMO-approximation of cp:
            zcpm     = cp_moist_cosmo(qv(i,j,k),qc(i,j,k),0.0_wp)
            t(i,j,k) = theta(i,j,k)*(piter(i,j,k)/pt00)**(zrdm/zcpm)
          END DO
        END DO
      END DO
    END IF

    !.. Debug output of the iterated pressure profile:
    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 4) THEN
      IF (ltheta_ini) THEN
        WRITE (*,'(a)') 'Subroutine calc_p_hydrostat_lf() after pressure '// &
             'iteration at point (1,1):'
        WRITE (*,'(a)') 'k     zml(k) theta(k) relhum(k)     t(k)  piter(k)      '// &
             'qv(k)        qc(k) relhum_out(k)'
        DO k=1,nk
          WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f8.1,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
               k,zml(1,1,k),theta(1,1,k),relhum(1,1,k),t(1,1,k), &
               piter(1,1,k),qv(1,1,k),qc(1,1,k), &
               rh_Tpqv(piter(1,1,k),t(1,1,k),qv(1,1,k),qc(1,1,k))
        END DO

      ELSE
        WRITE (*,'(a)') 'Subroutine calc_p_hydrostat_lf() after pressure '// &
             'iteration at point (1,1):'
        WRITE (*,'(a)') 'k     zml(k)     t(k) relhum(k)     p(k)      '// &
             'qv(k)        qc(k) relhum_out(k)'
        DO k=1,nk
          WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
               k,zml(1,1,k),t(1,1,k),relhum(1,1,k),piter(1,1,k),qv(1,1,k),qc(1,1,k), &
               rh_Tpqv(piter(1,1,k),t(1,1,k),qv(1,1,k),qc(1,1,k))
        END DO
      END IF
    END IF

    ! Convergence check:
    DO j = 1,nj
      DO i = 1,ni
        IF ( ANY( ABS(piter(i,j,:) - piterold(i,j,:)) > eps_pp ) ) THEN
          WRITE (*,'(a)') 'ERROR: in calc_p_hydrostat_lf(), src_artifdata.f90: '//&
               'no convergence in iteration of hydrostatic initialisation!'
          WRITE (*,'(a,i4,1x,i4,a,f12.4,a,f14.7)') 'at location (i,j) ', &
               isubpos(my_cart_id, 1)-nboundlines-1+i, &
               isubpos(my_cart_id, 2)-nboundlines-1+j, &
               '    psurf(i,j) = ', psurf(i,j), '   dp_last_max = ', &
               MAXVAL(ABS(piterold(i,j,:)-piter(i,j,:)))
          IF (ltheta_ini) THEN
            WRITE (*,'(a)') 'k     zml(k) theta(k) relhum(k)     t(k) piter(k)      '// &
                 'qv(k)        qc(k) relhum_out(k)'
            DO k=1,nk
              WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f8.1,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
                   k,zml(i,j,k),theta(i,j,k),relhum(i,j,k),t(i,j,k),&
                   piter(i,j,k),qv(i,j,k),qc(i,j,k), &
                   rh_Tpqv(piter(i,j,k),t(i,j,k),qv(i,j,k),qc(i,j,k))
            END DO
          ELSE
            WRITE (*,'(a)') 'k     zml(k)     t(k) relhum(k) piter(k)      '// &
                 'qv(k)        qc(k) relhum_out(k)'
            DO k=1,nk
              WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
                   k,zml(i,j,k),t(i,j,k),relhum(i,j,k),piter(i,j,k),qv(i,j,k),qc(i,j,k), &
                   rh_Tpqv(piter(i,j,k),t(i,j,k),qv(i,j,k),qc(i,j,k))
            END DO
          END IF
          ierror = 10009
          errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_lf(): No convergence in iteration of '// &
                   'hydrostatic initialisation of pressure field!'
          RETURN
        END IF
      END DO
    END DO
    
    !.. Issue a warning message if the relat. humid. has been altered above:
    IF (ANY(relhum_lim(:,:,:) /= relhum(:,:,:))) THEN
      WRITE (*,'(a)') REPEAT('*',70)
      WRITE (*,'(a)') 'WARNING: calc_p_hydrostat_lf: Rel. humid. limited to '// &
           'max. allowed value p / E(T)'
      WRITE (*,'(a)') REPEAT('*',70)
    END IF

  END SUBROUTINE calc_p_hydrostat_lf


!=====================================================================================
!=====================================================================================
!
! Subroutine for computing the analytic hydrostatic pressure profile
! for a given temperature t (or pot. temperature theta) profile.
! The discrete points in the t- (or theta-) profile are
! interpreted as a linear spline (polygone) w.r.t. t, so that t in the layer
! between two points is treated as a polytrope layer (linear
! dependence on height).
!
! The resulting solution for the pressure profile is
! perhaps the most accurate analytic solution, but is not
! necessarily equal to the numeric representation of 
! the hydrostatic balance used in the model dynamics!
! IN CONTRAST, the SR calc_p_hydrostat_psts() leads
! to a solution which is exactly compatible with the
! p'T'-dynamics of the model, and calc_p_hydrostat_lf()
! is compatible with leapfrog dynamics.
!
! - psurf is the pressure at the lowest height level zml(:,:,1)
!   and is the base for upwards pressure integration.
!
! - If theta (optional argument) is provided on input,
!   then this is used as the basis for the temperature
!   profile rather than directly using t, and t is
!   calculated and given back to the calling routine.  
!
! - The pressure integration takes into account moisture and
!   perhaps cloud droplets if relhum > 1 at some voxels.
! 
! - At the same time, relhum (and consequently Qv)  is 
!   limited to its upper bound of p / E(T).
!
! - Requires monotonically increasing height levels !!!
!
!
!=====================================================================================
!=====================================================================================

  SUBROUTINE calc_p_hydrostat_ana(ni, nj, nk, niter, &
       zml, psurf, t, relhum, &
       piter, qv, qc, zmaxqv, r_d, rvd_m_o, ierror, errmsg, &
       theta)

    IMPLICIT NONE

    !.. Input/output variables:
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk, niter
    REAL(KIND=wp),     INTENT(in) :: zml(ni,nj,nk), psurf(ni,nj), &
         relhum(ni,nj,nk)
    REAL(KIND=wp),     INTENT(out) :: qv(ni,nj,nk), qc(ni,nj,nk), piter(ni,nj,nk)
    REAL(KIND=wp),     INTENT(in) :: zmaxqv, r_d, rvd_m_o
    REAL(KIND=wp),     INTENT(inout) :: t(ni,nj,nk)
    INTEGER(kind=iintegers), INTENT(out) :: ierror
    CHARACTER(len=250), INTENT(out) :: errmsg
    REAL(KIND=wp),     INTENT(in), OPTIONAL :: theta(ni,nj,nk)

    !.. Local variables:
    REAL(KIND=wp)     :: piterold(ni,nj,nk), relhum_lim(ni,nj,nk)
    REAL(KIND=wp)     :: eps_pp, zpa, ztv1, ztv2, zesat, zsqv, zgamma, zrdm, zcpm
    INTEGER(KIND=iintegers) ::  i, j, k, iter
    LOGICAL :: ltheta_ini

    ierror = 0
    errmsg(:) = ' '

    eps_pp = 1.0E-4_wp     ! Pa, required abs. accuracy of the below pressure iteration

    IF (PRESENT(theta)) THEN
      ltheta_ini = .TRUE.
      ! Initialize qv and qc with their "dry" values so that zrdm and zcpm
      ! have a defined value in the first "moist" iteration:
      qv(:,:,:) = 0.0_wp
      qc(:,:,:) = 0.0_wp
    ELSE
      ltheta_ini = .FALSE.
      ! Initialize qv and qc with their "dry" values so that qv_Tprelhum()
      ! gets defined values in the first "moist" iteration:
      qv(:,:,:) = 0.0_wp
      qc(:,:,:) = 0.0_wp
    END IF

    !==========================================================================
    !.. First: Integrate atmospheric pressure assuming dry air. 
    !   This will be the starting point for
    !   a fixpoint iteration to include also moisture and condensation 
    !   in supersaturated voxels.
    !==========================================================================

    !   In case of ltheta_ini = .true., the temperature profile is not known a priori,
    !   because we would need the pressure to compute it from the theta-profile.
    !   As a starting point, the temperature will simply be that of the ICAO 
    !   polytrope atmosphere with a constant temperature lapse rate of 0.0065 K/m 
    !   up to 11 km height and a base temperature of 293.16 K. Then, the
    !   pressure estimate for the dry base will be that of the ICAO standard
    !   atmosphere. The below iteration will compute the correct temperature-
    !   and pressure profile afterwards.
    !    
    IF (ltheta_ini) THEN
      t(:,:,:)  = 293.16_wp - 0.0065_wp * MIN(zml(:,:,:), 11000.0_wp)
    END IF

    ! Pressure at the lowest level:
    piter(:,:,1) = psurf

    ! The other levels are polytrope layers assuming const. temp. gradient:
    DO k = 2, nk
      DO j=1,nj
        DO i=1,ni

          ztv1 = t(i,j,k-1)
          ztv2 = t(i,j,k)
          IF (ABS(ztv2-ztv1) > 1.0E-6_wp) THEN
            zgamma = (ztv1-ztv2) / (zml(i,j,k)-zml(i,j,k-1))
            piter(i,j,k) = piter(i,j,k-1) * &
                 (1.0_wp-(zml(i,j,k)-zml(i,j,k-1))*zgamma/ztv1)**(g/(r_d*zgamma))
          ELSE
            piter(i,j,k) = piter(i,j,k-1) * EXP(-g*(zml(i,j,k)-zml(i,j,k-1))/(r_d*ztv1))
          END IF

        END DO
      END DO
    END DO

    !.. Debug output of the iterated dry pressure profile:
    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 4) THEN
      WRITE (*,*) 'Subroutine calc_p_hydrostat_ana() after dry '// &
           'pressure iteration at point (1,1):'
      WRITE (*,*) 'k     zml(k)     p(k) for the dry atmosphere'
      DO k=1,nk
        WRITE (*,'(i4,1x,f8.1,1x,f12.2)') k,zml(1,1,k),piter(1,1,k)
      END DO
    END IF

    !==========================================================================
    ! Now humidity and clouds come into play:
    !
    ! The total density depends on qv, qc and therefore also on piter (through qv),
    ! but piter itself depends on qv via the hydrostatic approximation --> 
    ! iterative solution of this implicit equation for piter necessary!
    ! For this, each layer is again assumed to be a polytrope layer,
    ! this time with a constant virtual temperature gradient:
    !
    ! For the sake of reproducible results, we do a fixed number of iterations
    ! instead of iterating until convergence to a certain accuracy.
    !==========================================================================

    DO iter = 1, niter

      ! New variable for limited relhum to max. allowed value p / E(T):
      ! Re-initialize it to relhum for each iteration
      ! to avoid drifts of relhum_lim during the iteration:
      relhum_lim(:,:,:) = relhum(:,:,:)

      piterold(:,:,:) = piter(:,:,:)

      DO j=1,nj
        DO i=1,ni

          ! Compute t, qv, qc for the lowest level:

          ! decomposition of the relative humidity into qv, qc using the 
          ! momentary pressure zpa
          zpa          = piter(i,j,1)
          ! saturation specific vapor content (Spezifische Feuchte bei Saettigung):
          IF (ltheta_ini) THEN
            ! "moist" r_d:
            zrdm     = rd_moist(qv(i,j,1),qc(i,j,1))
            ! "moist" cp:
!!$            zcpm = cp_moist(qv(i,j,1),qc(i,j,1),0.0_wp)
            ! COSMO-approximation of "moist" cp:
            zcpm     = cp_moist_cosmo(qv(i,j,1),qc(i,j,1),0.0_wp)
            t(i,j,1) = theta(i,j,1)*(zpa/pt00)**(zrdm/zcpm)
          END IF
!!$          IF (t(i,j,1) > 233.16) THEN
            zesat       = esat_w(t(i,j,1))
!!$          ELSE
!!$            zesat       = esat_i(t(i,j,1))
!!$          ENDIF
          ! Limit relhum to its maximum possible value p / E(T):
          relhum_lim(i,j,1) = MIN( zpa/zesat , relhum(i,j,1) )
          ! Compute actual Qv:
          qv(i,j,1) = qv_Tprelhum( zpa, t(i,j,1), relhum_lim(i,j,1), qc(i,j,1) )
          IF (lcond .AND. relhum_lim(i,j,1) > 1.0_wp) THEN
            ! condensation is allowed and physically can happen
            ! and relhum > 1.0, so convert qv -> qc to limit relhum to 1.0:
            zsqv = qvsat_w( zpa, t(i,j,1) )
            qv(i,j,1) = MIN ( zsqv       ,  MIN(qv(i,j,1), zmaxqv) )
            qc(i,j,1) = MAX ( 0.0_wp ,  MIN(qv(i,j,1), zmaxqv) - zsqv )
          ELSE
            ! condensation is not allowed or cannot happen physically at 
            ! that pressure and temperature, so just impose the limit zmaxqv:
            qv(i,j,1) = MIN(qv(i,j,1), zmaxqv)
            qc(i,j,1) = 0.0_wp
          END IF

        END DO
      END DO

      ! Again, the other levels are polytrope layers assuming const. virt. temp. gradient:
      DO k = 2, nk

        DO j=1,nj
          DO i=1,ni

            !.. virtual temperature ztv1 at the bottom of the layer:
            ztv1 = t(i,j,k-1) * (1.0_wp + rvd_m_o*qv(i,j,k-1) - qc(i,j,k-1))

            !.. virtual temperature ztv2 at the top of the layer:
            zpa          = piter(i,j,k)
            IF (ltheta_ini) THEN
              ! "moist" r_d:
              zrdm     = rd_moist(qv(i,j,k),qc(i,j,k))
              ! "moist" cp:
!!$              zcpm = cp_moist(qv(i,j,k),qc(i,j,k),0.0_wp)
              ! COSMO-approximation o "moist" cp:
              zcpm     = cp_moist_cosmo(qv(i,j,k),qc(i,j,k),0.0_wp)
              t(i,j,k) = theta(i,j,k)*(zpa/pt00)**(zrdm/zcpm)
            END IF
!!$            IF (t(i,j,k) > 233.16) THEN
              zesat       = esat_w(t(i,j,k))
!!$            ELSE
!!$              zesat       = esat_i(t(i,j,k))
!!$            ENDIF
            ! Limit relhum to its maximum possible value p / E(T):
            relhum_lim(i,j,k) = MIN( zpa/zesat , relhum(i,j,k) )
            ! Compute actual Qv:
            qv(i,j,k) = qv_Tprelhum( zpa, t(i,j,k), relhum_lim(i,j,k), qc(i,j,k) )

            IF (lcond .AND. relhum_lim(i,j,k) > 1.0_wp) THEN
              ! condensation is allowed and physically can happen
              ! and relhum > 1.0, so convert qv -> qc to limit relhum to 1.0:
              zsqv = qvsat_w( zpa, t(i,j,k) )
              qv(i,j,k) = MIN ( zsqv       ,  MIN(qv(i,j,k), zmaxqv) )
              qc(i,j,k) = MAX ( 0.0_wp ,  MIN(qv(i,j,k), zmaxqv) - zsqv )
            ELSE
              ! condensation is not allowed, just set the limit to zmaxqv:
              qv(i,j,k) = MIN(qv(i,j,k), zmaxqv)
              qc(i,j,k) = 0.0_wp
            END IF

!!$            IF (i==1 .AND. j==1 .AND. k==nk .AND. ldebug_artif .AND. &
!!$                 idbg_artif_level > 4) THEN
!!$              WRITE (*,'(a,i4,5(1x,es11.4))') 'calc_p_hydrostat_ana(): '// &
!!$                   'ITERATION CHECK VALUES AT TOP: ', &
!!$                   iter, zpa, zesat, zpa/zesat, relhum_lim(i,j,k), qv(i,j,k)
!!$            END IF

            ztv2   = t(i,j,k) * (1.0_wp + rvd_m_o*qv(i,j,k) - qc(i,j,k))

            !.. pressure at the top of the polytrope layer:
            IF (ABS(ztv2-ztv1) > 1.0E-6_wp) THEN
              zgamma = (ztv1-ztv2) / (zml(i,j,k)-zml(i,j,k-1))
              piter(i,j,k) = piter(i,j,k-1) * &
                   (1.0_wp-(zml(i,j,k)-zml(i,j,k-1))*zgamma/ztv1)**(g/(r_d*zgamma))
            ELSE
              piter(i,j,k) = piter(i,j,k-1) * EXP(-g*(zml(i,j,k)-zml(i,j,k-1))/(r_d*ztv1))
            END IF

          END DO
        END DO

      END DO

      ! Catch a badly diverged iteration, which can easily happen in case
      !   of initialization with N-const profiles (itype_anaprof_tqv = 3):
      IF ( ANY(piter < 0.0_wp .OR. t < 0.0_wp .OR. t > 1000.0_wp) ) THEN
        !.. For meaningful error message, print first erroneous point:
        DO j = 1,nj
          DO i = 1,ni
            IF ( ANY( piter(i,j,:) < 0.0_wp .OR. t(i,j,:) < 0.0_wp .OR. t(i,j,:) > 1000.0_wp ) ) THEN

              WRITE (*,'(a)') 'ERROR: Subroutine calc_p_hydrostat_ana() failed during pressure iteration'
              WRITE (*,'(a,i4,1x,i4)') '    at location (i,j) ', &
                   isubpos(my_cart_id, 1)-nboundlines-1+i, &
                   isubpos(my_cart_id, 2)-nboundlines-1+j
              IF (ltheta_ini) THEN
                WRITE (*,'(a)') '   k   zml(k) theta(k)    t(k) relhum(k)         p(k)        '// &
                     'qv(k)        qc(k)'
                DO k=1,nk
                  WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5))') &
                       k,zml(i,j,k),theta(i,j,k),t(i,j,k),relhum(i,j,k),piter(i,j,k),qv(i,j,k),qc(i,j,k)
                END DO
                ierror = 10085
                errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_ana(): Problem in iteration of '// &
                     'hydrostatic initialisation of pressure field! '// &
                     'Check namelist parameters for stability (N) and/or try to reduce the model top height!'
              ELSE
                WRITE (*,'(a)') '   k   zml(k)    t(k) relhum(k)           p(k)        '// &
                     'qv(k)        qc(k)'
                DO k=1,nk
                  WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5))') &
                       k,zml(i,j,k),t(i,j,k),relhum(i,j,k),piter(i,j,k),qv(i,j,k),qc(i,j,k)
                END DO
                ierror = 10086
                errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_ana(): Problem in iteration of '// &
                     'hydrostatic initialisation of pressure field! '// &
                     'Initial T profile seems to be problematic!'
              END IF

              RETURN

            END IF
          END DO
        END DO
      END IF

    END DO

    !.. If necessary, re-calculate T from Theta and iterated pressure:
    IF (ltheta_ini) THEN
      DO k = 1,nk
        DO j = 1,nj
          DO i = 1,ni
            ! "moist" r_d:
            zrdm     = rd_moist(qv(i,j,k),qc(i,j,k))
            ! COSMO-approximation of cp:
            zcpm     = cp_moist_cosmo(qv(i,j,k),qc(i,j,k),0.0_wp)
            t(i,j,k) = theta(i,j,k)*(piter(i,j,k)/pt00)**(zrdm/zcpm)
          END DO
        END DO
      END DO
    END IF

    !.. Debug output of the iterated pressure profile:
    IF (my_cart_id == 0 .AND. ldebug_artif .AND. idbg_artif_level > 4) THEN
      IF (ltheta_ini) THEN
        WRITE (*,'(a)') 'Subroutine calc_p_hydrostat_ana() after pressure '// &
             'iteration at point (1,1):'
        WRITE (*,'(a)') 'k     zml(k) theta(k) relhum(k)     t(k)  piter(k)      '// &
             'qv(k)        qc(k) relhum_out(k)'
        DO k=1,nk
          WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f8.1,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
               k,zml(1,1,k),theta(1,1,k),relhum(1,1,k),t(1,1,k), &
               piter(1,1,k),qv(1,1,k),qc(1,1,k), &
               rh_Tpqv(piter(1,1,k),t(1,1,k),qv(1,1,k),qc(1,1,k))
        END DO

      ELSE
        WRITE (*,'(a)') 'Subroutine calc_p_hydrostat_ana() after pressure '// &
             'iteration at point (1,1):'
        WRITE (*,'(a)') 'k     zml(k)     t(k) relhum(k)     p(k)      '// &
             'qv(k)        qc(k) relhum_out(k)'
        DO k=1,nk
          WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
               k,zml(1,1,k),t(1,1,k),relhum(1,1,k),piter(1,1,k),qv(1,1,k),qc(1,1,k), &
               rh_Tpqv(piter(1,1,k),t(1,1,k),qv(1,1,k),qc(1,1,k))
        END DO
      END IF
    END IF

    ! Convergence check:
    DO j = 1,nj
      DO i = 1,ni
        IF ( ANY( ABS(piter(i,j,:) - piterold(i,j,:)) > eps_pp ) ) THEN
          WRITE (*,'(a)') 'ERROR: in calc_p_hydrostat_ana(), src_artifdata.f90: '//&
               'no convergence in iteration of hydrostatic initialisation!'
          WRITE (*,'(a,i4,1x,i4,a,f12.4,a,f14.7)') 'at location (i,j) ', &
               isubpos(my_cart_id, 1)-nboundlines-1+i, &
               isubpos(my_cart_id, 2)-nboundlines-1+j, &
               '    psurf(i,j) = ', psurf(i,j), '   dp_last_max = ', &
               MAXVAL(ABS(piterold(i,j,:)-piter(i,j,:)))
          IF (ltheta_ini) THEN
            WRITE (*,'(a)') 'k     zml(k) theta(k) relhum(k)     t(k) piter(k)      '// &
                 'qv(k)        qc(k) relhum_out(k)'
            DO k=1,nk
              WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f8.1,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
                   k,zml(i,j,k),theta(i,j,k),relhum(i,j,k),t(i,j,k), &
                   piter(i,j,k),qv(i,j,k),qc(i,j,k), &
                   rh_Tpqv(piter(i,j,k),t(i,j,k),qv(i,j,k),qc(i,j,k))
            END DO
          ELSE
            WRITE (*,'(a)') 'k     zml(k)     t(k) relhum(k) piter(k)      '// &
                 'qv(k)        qc(k) relhum_out(k)'
            DO k=1,nk
              WRITE (*,'(i4,1x,f8.1,1x,f8.1,1x,f8.5,1x,f12.2,2(1x,es12.5),1x,f8.5)') &
                   k,zml(i,j,k),t(i,j,k),relhum(i,j,k),piter(i,j,k),qv(i,j,k),qc(i,j,k), &
                   rh_Tpqv(piter(i,j,k),t(i,j,k),qv(i,j,k),qc(i,j,k))
            END DO
          END IF

          ierror = 10010
          errmsg = 'ERROR src_artifdata.f90, calc_p_hydrostat_ana(): No convergence in iteration of '// &
                   'hydrostatic initialisation of pressure field!'

          RETURN
        END IF
      END DO
    END DO

    !.. Issue a warning message if the relat. humid. has been altered above:
    IF (ANY(relhum_lim(:,:,:) /= relhum(:,:,:))) THEN
      WRITE (*,'(a)') REPEAT('*',70)
      WRITE (*,'(a)') 'WARNING: calc_p_hydrostat_ana: Rel. humid. limited '// &
           'to max. allowed value p / E(T)'
      WRITE (*,'(a)') REPEAT('*',70)
    END IF

  END SUBROUTINE calc_p_hydrostat_ana

!=====================================================================================
!=====================================================================================

!=====================================================================================
!
! Functions for computing moisture quantities:
!
!=====================================================================================

  ! saturation vapor pressure w.r.t. plain water surface (scalar version):
  REAL(KIND=wp)     FUNCTION esat_w(temp)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: temp
    esat_w = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
  END FUNCTION esat_w

  ! saturation vapor pressure w.r.t. plain ice surface (scalar version):
  REAL(KIND=wp)     FUNCTION esat_i(temp)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: temp
    esat_i = b1 * EXP( b2i * (temp-b3)/(temp-b4i) )
  END FUNCTION esat_i

  ! Derivative of saturation vapor pressure w.r.t. plain water surface (scalar version):
  REAL(KIND=wp)     FUNCTION desat_w_dT(temp)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: temp
    desat_w_dT = b234w/(temp-b4w)**2 * b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
  END FUNCTION desat_w_dT
  
  ! saturation specific humidity w.r.t. plain water surface (scalar version):
  ! RETURNS -999.99, IF THE SATURATION VAPOR PRESSURE IS > TOTAL PRESSURE!
  !   (IN THIS CASE, QVSAT IS NO LONGER A MEANINGFUL QUANTITY,
  !    BECAUSE SATURATION CANNOT OCCUR AT THAT PRESSURE AND TEMPERATURE)
  REAL(KIND=wp)     FUNCTION qvsat_w(p, temp)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: p, temp
    REAL(KIND=wp)     :: zesat_w
    zesat_w  = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
    IF (zesat_w <= p) THEN
      qvsat_w = rdv * zesat_w / ( p - o_m_rdv * zesat_w )
    ELSE
      qvsat_w = -999.99_wp
    END IF
  END FUNCTION qvsat_w

  ! saturation specific humidity w.r.t. plain ice surface (scalar version):
  ! RETURNS -999.99, IF THE SATURATION VAPOR PRESSURE IS > TOTAL PRESSURE!
  !   (IN THIS CASE, QVSAT IS NO LONGER A MEANINGFUL QUANTITY,
  !    BECAUSE SATURATION CANNOT OCCUR AT THAT PRESSURE AND TEMPERATURE)
  REAL(KIND=wp)     FUNCTION qvsat_i(p, temp)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: p, temp
    REAL(KIND=wp)     :: zesat_i
    zesat_i  = b1 * EXP( b2i * (temp-b3)/(temp-b4i) )
    IF (zesat_i <= p) THEN
      qvsat_i = rdv * zesat_i / ( p - o_m_rdv * zesat_i )
    ELSE
      qvsat_i = -999.99_wp
    END IF
  END FUNCTION qvsat_i

  ! saturation vapor pressure w.r.t. plain water surface (3D version):
  FUNCTION esat_w_3d(temp, ni, nj, nk)
    IMPLICIT NONE
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk
    REAL(KIND=wp),     INTENT(in) :: temp(ni,nj,nk)
    REAL(KIND=wp)                 :: esat_w_3d(ni,nj,nk)
    esat_w_3d = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
  END FUNCTION esat_w_3d

  ! saturation vapor pressure w.r.t. plain ice surface (3D version):
  FUNCTION esat_i_3d(temp, ni, nj, nk)
    IMPLICIT NONE
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk
    REAL(KIND=wp),     INTENT(in) :: temp(ni,nj,nk)
    REAL(KIND=wp)                 :: esat_i_3d(ni,nj,nk)
    esat_i_3d = b1 * EXP( b2i * (temp-b3)/(temp-b4i) )
  END FUNCTION esat_i_3d

  ! Derivative of saturation vapor pressure w.r.t. plain water surface (scalar version):
  FUNCTION desat_w_dT_3d(temp, ni, nj, nk)
    IMPLICIT NONE
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk
    REAL(KIND=wp),     INTENT(in) :: temp(ni,nj,nk)
    REAL(KIND=wp)                 :: desat_w_dT_3d(ni,nj,nk)
    desat_w_dT_3d = b234w/(temp-b4w)**2 * b1*EXP( b2w*(temp-b3)/(temp-b4w) )
  END FUNCTION desat_w_dT_3d
  
  ! saturation specific humidity w.r.t. plain water surface (3D version):
  ! RETURNS -999.99, IF THE SATURATION VAPOR PRESSURE IS > TOTAL PRESSURE!
  !   (IN THIS CASE, QVSAT IS NO LONGER A MEANINGFUL QUANTITY,
  !    BECAUSE SATURATION CANNOT OCCUR AT THAT PRESSURE AND TEMPERATURE)
  FUNCTION qvsat_w_3d(p, temp, ni, nj, nk)
    IMPLICIT NONE
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk
    REAL(KIND=wp),     INTENT(in) :: p(ni,nj,nk), temp(ni,nj,nk)
    REAL(KIND=wp)                 :: qvsat_w_3d(ni,nj,nk)
    REAL(KIND=wp)     :: zesat_w(ni,nj,nk)
    zesat_w  = esat_w_3d(temp,ni,nj,nk)
    WHERE (zesat_w <= p) 
      qvsat_w_3d = rdv * zesat_w / ( p - o_m_rdv * zesat_w )
    ELSEWHERE
      qvsat_w_3d = -999.99_wp
    END WHERE
  END FUNCTION qvsat_w_3d

  ! saturation specific humidity w.r.t. plain ice surface (3D version):
  ! RETURNS -999.99, IF THE SATURATION VAPOR PRESSURE IS > TOTAL PRESSURE!
  !   (IN THIS CASE, QVSAT IS NO LONGER A MEANINGFUL QUANTITY,
  !    BECAUSE SATURATION CANNOT OCCUR AT THAT PRESSURE AND TEMPERATURE)
  FUNCTION qvsat_i_3d(p, temp, ni, nj, nk)
    IMPLICIT NONE
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk
    REAL(KIND=wp),     INTENT(in) :: p(ni,nj,nk), temp(ni,nj,nk)
    REAL(KIND=wp)                 :: qvsat_i_3d(ni,nj,nk)
    REAL(KIND=wp)     :: zesat_i(ni,nj,nk)
    zesat_i  = esat_i_3d(temp,ni,nj,nk)
    WHERE (zesat_i <= p) 
      qvsat_i_3d = rdv * zesat_i / ( p - o_m_rdv * zesat_i )
    ELSEWHERE
      qvsat_i_3d = -999.99_wp
    END WHERE
  END FUNCTION qvsat_i_3d

  ! Specific humidity as function of T, p, and relhum 
  !   (and qcrs = sum of all hydrometeor contents):
  ! NOTE: on input, relhum has to be smaller than p / E(T)!
  REAL(KIND=wp)     FUNCTION qv_Tprelhum(p, temp, relhum, qcrs)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: p, temp, relhum, qcrs
    REAL(KIND=wp)     :: zesat_w, coeff
    zesat_w  = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
    coeff = relhum * zesat_w * rdv / p
    qv_Tprelhum = coeff * (1.0_wp + qcrs) / (1.0_wp - rvd_m_o*coeff)
  END FUNCTION qv_Tprelhum

  ! Specific humidity as function of T, p, and relhum (3D-version)
  !   (and qcrs = sum of all hydrometeor contents):
  ! NOTE: on input, relhum has to be smaller than p / E(T)!
  FUNCTION qv_Tprelhum_3d(p, temp, relhum, qcrs, ni, nj, nk)
    IMPLICIT NONE
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk
    REAL(KIND=wp),     INTENT(in), DIMENSION(ni,nj,nk) :: p, temp, relhum, qcrs
    REAL(KIND=wp),                 DIMENSION(ni,nj,nk) :: qv_Tprelhum_3d
    REAL(KIND=wp),                 DIMENSION(ni,nj,nk) :: zesat_w, coeff
    zesat_w  = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
    coeff   = relhum * zesat_w * rdv / p
    qv_Tprelhum_3d = coeff * (1.0_wp + qcrs) / (1.0_wp - rvd_m_o*coeff)
  END FUNCTION qv_Tprelhum_3d


  ! Rel. humidity as function of T, p, and qv
  !   (and qcrs = sum of all hydrometeor contents):
  REAL(KIND=wp)     FUNCTION rh_Tpqv(p, temp, qv, qcrs)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: p, temp, qv, qcrs
    REAL(KIND=wp)     :: zesat_w
    zesat_w  = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
    rh_Tpqv = p * qv / (rdv * zesat_w * (1.0_wp+rvd_m_o*qv-qcrs) )
  END FUNCTION rh_Tpqv

  ! Rel. humidity as function of T, p, and qv (3D-version)
  !   (and qcrs = sum of all hydrometeor contents):
  FUNCTION rh_Tpqv_3d(p, temp, qv, qcrs, ni, nj, nk)
    IMPLICIT NONE
    INTEGER(KIND=iintegers), INTENT(in) :: ni, nj, nk
    REAL(KIND=wp),     INTENT(in), DIMENSION(ni,nj,nk) :: p, temp, qv, qcrs
    REAL(KIND=wp),                 DIMENSION(ni,nj,nk) :: rh_Tpqv_3d
    REAL(KIND=wp),                 DIMENSION(ni,nj,nk) :: zesat_w
    zesat_w     = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
    rh_Tpqv_3d = p * qv / (rdv * zesat_w * (1.0_wp+rvd_m_o*qv-qcrs) )
  END FUNCTION rh_Tpqv_3d


  ! Gas constant of moist air containing hydrometeors (qcrs is the sum of
  ! the specific hydrometeor contents):
  REAL(KIND=wp)     FUNCTION rd_moist(qv,qcrs)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: qv, qcrs

    rd_moist = r_d * (1.0_wp + rvd_m_o*qv - qcrs)

  END FUNCTION rd_moist

  ! Cp of moist air containing hydrometeors (ql is the sum of
  ! the liquid specific hydrometeor contents, qi the sum of the specific ice
  ! hydrometeor contents):
  REAL(KIND=wp)     FUNCTION cp_moist(qv,ql,qi)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: qv, ql, qi
    REAL(KIND=wp)                 :: qd

    ! Dry air content:
    qd = 1.0_wp - qv - ql -qi

    ! cp:
!    cp_moist = qd*cp_d + qv*cp_v + ql*cp_l + qi*cp_i
    cp_moist = qd*cp_d + qv*(cp_d*(1.0_wp+rcpv)) + ql*(cp_d*(1.0_wp+rcpl)) + qi*2060.0_wp

  END FUNCTION cp_moist

  ! COSMO-APPROXIMATION: CP IS APPROXIMATED TO BE THAT OF DRY AIR
  REAL(KIND=wp)     FUNCTION cp_moist_cosmo(qv,ql,qi)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: qv, ql, qi
    REAL(KIND=wp)                 :: qd

    ! Dry air content:
    qd = 1.0_wp - qv - ql -qi

    ! cp:
    cp_moist_cosmo = qd*cp_d + qv*cp_d + ql*cp_d + qi*cp_d

  END FUNCTION cp_moist_cosmo

  ! Qvtens = dqv/dt at constant p, relhum and qcrs (total hydrometeor content) for a given
  !   temperature temp and temperature tendency dTdt = dtemp/dt:
  REAL(KIND=wp)     FUNCTION dqvdt_prh (temp, dTdt, p, relhum, qcrs)
    IMPLICIT NONE
    REAL(KIND=wp),     INTENT(in) :: p, temp, relhum, qcrs, dTdt
    REAL(KIND=wp)     :: zesat_w, zdesat_w, coeff, qdv, fakt1, fakt2

    zesat_w = b1 * EXP( b2w * (temp-b3)/(temp-b4w) )
    zdesat_w = b234w/(temp-b4w)**2 * zesat_w
    fakt1 = relhum * rdv / p
    qdv = 1.0_wp - qcrs
    fakt2 = 1.0_wp / (1.0_wp - rvd_m_o*zesat_w*fakt1)
    coeff = fakt1*qdv*fakt2 + zesat_w*fakt1*fakt1*qdv*fakt2*fakt2*rvd_m_o
    
    dqvdt_prh = coeff * zdesat_w * dTdt

  END FUNCTION dqvdt_prh

!=================================================================================
!=================================================================================

!=================================================================================
!
! Functions for coordinate operations in spherical coordinates:
!
!=================================================================================

  ! Distances relative to hill/bubble main axes (these can be rotated by <rotangle> 
  ! relative to the rotated North direction) for coordinates
  ! given at the mass points. The coordinates must be global (i,j) coordinates, 
  ! not these for the local PE!
  ! 
  ! 1) If lmetr=.false., these distances are just the normal Euklidian distances
  !    from the respective point to the main axes.
  ! 
  ! 2) Else, the main axes are along great circles and the distances are 
  !    also measured along great circles.
  ! 
  SUBROUTINE hill_rot_coords(i_td, j_td, ic_td, jc_td, rotangle, height, rx, ry)

    IMPLICIT NONE

    !.. Input/Output parameters:
    !---------------------------
    !   global index of model grid point:
    INTEGER(KIND=iintegers), INTENT(in) :: i_td, j_td
    !   global index of center of bubble/hill:
    REAL(KIND=wp),     INTENT(in) :: ic_td, jc_td
    !   rotation angle of main hill/bubble y-axis clockwise relative to north in rad:
    REAL(KIND=wp),     INTENT(in)       :: rotangle
    !   height level where the arc lengths are referenced to in m:
    REAL(KIND=wp),     INTENT(in)       :: height
    !   arc lengths along great circles perpendicular to the 
    !   hill/bubble main axes (also great circles):
    REAL(KIND=wp),     INTENT(out)      :: rx, ry


    !.. Local variables:
    !---------------------------
    !   angles and arc lengths on the unit sphere in rad:
    REAL(KIND=wp)     ::  d, delta, &
         zlon, zlat, zlon_c, zlat_c, tmp, zdi, zdj, zdx, zdy
#ifdef __SX__
    REAL(KIND=wp)     ::  cangle, cos_cangle, cos_arg
#endif

    IF (lmetr) THEN

      zlon   = ( startlon_tot + (i_td -1.0_wp)*dlon ) * degrad
      zlat   = ( startlat_tot + (j_td -1.0_wp)*dlat ) * degrad
      zlon_c = ( startlon_tot + (ic_td-1.0_wp)*dlon ) * degrad
      zlat_c = ( startlat_tot + (jc_td-1.0_wp)*dlat ) * degrad

#ifdef __SX__

!!$--------------- inlined code from geo_dist -----------------------------------
      d =  ACOS( SIN(zlat)*SIN(zlat_c)+COS(zlat)*COS(zlat_c)*COS(zlon-zlon_c) )
!!$--------------- end inlined code ---------------------------------------------

!!$--------------- inlined code from geo_heading --------------------------------
      cos_cangle =  SIN(zlat)*SIN(zlat_c)+COS(zlat)*COS(zlat_c)*COS(zlon-zlon_c) 
      cangle = ACOS(cos_cangle)
      IF (ABS(cangle) < 1.0E-20_wp) cangle = 1.0E-20_wp

      cos_arg = (SIN(zlat)-SIN(zlat_c)*cos_cangle) / (COS(zlat_c)*SIN(cangle))
      cos_arg = MAX(MIN(cos_arg, 1.0_wp), -1.0_wp)

      delta = ACOS(cos_arg)
      IF (zlon_c > zlon) delta = 2.0_wp*pi - delta
!!$--------------- end inlined code ---------------------------------------------

#else

      d     =  geo_dist   ( zlon, zlat, zlon_c, zlat_c, 0.0_wp ) / r_earth
      delta =  geo_heading( zlon, zlat, zlon_c, zlat_c )

#endif

      !.. take rotation into account:
      delta = delta - rotangle

      rx = ASIN(SIN(d)*SIN(delta))
      tmp = 1.0_wp - SIN(delta)*SIN(rx)*SIN(d)
      IF (ABS(tmp) < 1.0E-20_wp) THEN
        ry = 0.0_wp
      ELSE
        ry = ACOS(COS(rx)*COS(d) / tmp)
        IF (zlat_c > zlat) ry = -ry
      END IF
      rx = (r_earth + height) * rx
      ry = (r_earth + height) * ry

    ELSE

      ! grid lengths (without metrical term because lmetr = .false.):
      zdx   = r_earth * dlon * degrad
      zdy   = r_earth * dlat * degrad
      ! index distance to hill/bubble main axes (unrotated hill/bubble):
      zdi = i_td - ic_td
      zdj = j_td - jc_td
      ! physical distance to hill/bubble main axes (rotated hill/bubble):
      rx = zdi*COS(rotangle)*zdx - zdj*SIN(rotangle)*zdy
      ry = zdi*SIN(rotangle)*zdx + zdj*COS(rotangle)*zdy

    END IF

    RETURN
  END SUBROUTINE hill_rot_coords

  FUNCTION geo_dist(zlon, zlat, zlon_c, zlat_c, height) RESULT (dist)

    IMPLICIT NONE

    !   coordinates of the target grid point for the great circle in rad:
    REAL(KIND=wp),     INTENT(in)       :: zlon, zlat
    !   coordinates of the start grid point for the great circle in rad:
    REAL(KIND=wp),     INTENT(in)       :: zlon_c, zlat_c
    !   height level where the arc lengths are referenced to in m:
    REAL(KIND=wp),     INTENT(in)       :: height

    !   great circle distance between the two points in m:
    REAL(KIND=wp)                       :: dist

    dist = (r_earth+height) * ACOS( &
         SIN(zlat)*SIN(zlat_c)+COS(zlat)*COS(zlat_c)*COS(zlon-zlon_c) )

  END FUNCTION geo_dist

  FUNCTION geo_heading(zlon, zlat, zlon_c, zlat_c) RESULT (truecourse)

    IMPLICIT NONE

    !   coordinates of the target grid point for the great circle in rad:
    REAL(KIND=wp),     INTENT(in)       :: zlon, zlat
    !   coordinates of the start grid point for the great circle in rad:
    REAL(KIND=wp),     INTENT(in)       :: zlon_c, zlat_c
    
    !   geogr. direction from the start point to the target point
    !   at the start point in rad:
    REAL(KIND=wp)                       :: truecourse

    REAL(KIND=wp)     :: cangle, cos_cangle, cos_arg

    cos_cangle =  SIN(zlat)*SIN(zlat_c)+COS(zlat)*COS(zlat_c)*COS(zlon-zlon_c) 
    cangle = ACOS(cos_cangle)
    IF (ABS(cangle) < 1.0E-20_wp) cangle = 1.0E-20_wp

    cos_arg = (SIN(zlat)-SIN(zlat_c)*cos_cangle) / (COS(zlat_c)*SIN(cangle))
    cos_arg = MAX(MIN(cos_arg, 1.0_wp), -1.0_wp)

    truecourse = ACOS(cos_arg)
    IF (zlon_c > zlon) truecourse = 2.0_wp*pi - truecourse

  END FUNCTION geo_heading


!=================================================================================
!=================================================================================

!=================================================================================
!
! Functions for computing p, T, theta, rho as function of z of atmosphere
! layers with const. Brunt-Vaisala-frequency N (dry air only).
!
!=================================================================================

  !=================================================================================
  !
  !   Pot. temp. theta at height z of a gas with const. Brunt-Vaisala-Frequency n_bv
  !   as function of base height z0 and base pot. temp. theta0 at height z0
  !
  !   UNTESTED YET!
  !
  !=================================================================================

  FUNCTION  theta_z_N_const( z, z0, n_bv, theta0 ) &
       RESULT (theta)

    USE data_constants  , ONLY :   &
         g               ! acceleration due to gravity

    IMPLICIT NONE
    REAL (KIND=wp),     INTENT(in) :: z, z0, n_bv, theta0
    REAL (KIND=wp)     :: theta

    theta = theta0 * EXP(n_bv**2 / g * (z-z0))

  END FUNCTION theta_z_N_const

  !=================================================================================
  !
  !   Temp. T at height z of air with const. Brunt-Vaisala-Frequency n_bv
  !   as function of base height z0 and base temp. T0 at height z0
  !
  !   NOTE: BECAUSE WE APPROXIMATE CP OF MOIST AIR BY CP_D IN THE 
  !       COSMO-MODEL, IT IS ALSO APPLICABLE FOR MOIST AIR (AIR + 
  !       WATER VAPOR) HERE.
  !
  !   UNTESTED YET!
  !
  !=================================================================================

  FUNCTION  T_z_N_const_dry( z, z0, n_bv, T0 ) &
       RESULT (Tz)

    USE data_constants  , ONLY :   &
         g,            & ! acceleration due to gravity
         cp_d            ! specific heat of dry air at constant pressure

    IMPLICIT NONE
    REAL (KIND=wp),     INTENT(in) :: z, z0, n_bv, T0
    REAL (KIND=wp)     :: Tz
    REAL (KIND=wp)     :: coeff

    coeff = g**2 / (n_bv**2 * cp_d)
    Tz = coeff + (T0 - coeff) * EXP(n_bv**2 / g * (z-z0))

  END FUNCTION T_z_N_const_dry

  !=================================================================================
  !
  !   Pressure p at height z of dry air with const. Brunt-Vaisala-Frequency n_bv
  !   as function of base height z0, base pressure p0 at height z0 and 
  !   base temp. T0 at height z0
  !
  !   UNTESTED YET!
  !
  !=================================================================================

  FUNCTION  p_z_N_const_dry_t( z, z0, n_bv, p0, T0 ) &
       RESULT (pz)

    USE data_constants  , ONLY :   &
         g,            & ! acceleration due to gravity
         r_d,          & ! gas constant for dry air
         cp_d            ! specific heat of dry air at constant pressure

    IMPLICIT NONE
    REAL (KIND=wp),     INTENT(in) :: z, z0, n_bv, p0, T0
    REAL (KIND=wp)     :: pz
    REAL (KIND=wp)     :: coeff

    coeff = g**2 / (T0 * n_bv**2 * cp_d)
    pz = p0 * EXP( (cp_d/r_d) * LOG(1.0_wp + coeff*(EXP(-n_bv**2/g*(z-z0)) - 1.0_wp) ) )

  END FUNCTION p_z_N_const_dry_t

  !=================================================================================
  !
  !   Pressure p at height z of dry air with const. Brunt-Vaisala-Frequency n_bv
  !   as function of base height z0, base pressure p0 at height z0 and 
  !   base pot. temp. theta0 at height z0
  !
  !   p00:   reference pressure for pot. temperature, usually 1.0E5 Pa
  !
  !   UNTESTED YET!
  !
  !=================================================================================

  FUNCTION  p_z_N_const_dry_theta( z, z0, n_bv, p0, p00, theta0 ) &
       RESULT (pz)

    USE data_constants  , ONLY :   &
         g,            & ! acceleration due to gravity
         r_d,          & ! gas constant for dry air
         cp_d            ! specific heat of dry air at constant pressure

    IMPLICIT NONE
    REAL (KIND=wp),     INTENT(in) :: z, z0, n_bv, p0, p00, theta0
    REAL (KIND=wp)     :: pz
    REAL (KIND=wp)     :: coeff, T0

    T0 = theta0 * (p0/p00)**(r_d/cp_d)
    coeff = g**2 / (T0 * n_bv**2 * cp_d)
    pz = p0 * EXP( (cp_d/r_d) * LOG(1.0_wp + coeff*(EXP(-n_bv**2/g*(z-z0)) - 1.0_wp) ) )

  END FUNCTION p_z_N_const_dry_theta

  !=================================================================================
  !
  !   Density rho at height z of dry air with const. Brunt-Vaisala-Frequency n_bv
  !   as function of base height z0, base pressure p0 at height z0 and 
  !   base temp. T0 at height z0
  !
  !   UNTESTED YET!
  !
  !=================================================================================

  FUNCTION  rho_z_N_const_dry_t( z, z0, n_bv, p0, T0 ) &
       RESULT (rhoz)

    USE data_constants  , ONLY :   &
         g,            & ! acceleration due to gravity
         r_d,          & ! gas constant for dry air
         cp_d            ! specific heat of dry air at constant pressure

    IMPLICIT NONE
    REAL (KIND=wp),     INTENT(in) :: z, z0, n_bv, p0, T0
    REAL (KIND=wp)     :: rhoz
    REAL (KIND=wp)     :: coeff, pz, Tz

    !.. No futher subroutine- or function calls for Tz and pz for better
    !   vectorization properties!

    coeff = g**2 / (T0 * n_bv**2 * cp_d)
    Tz = T0 * (coeff + (1.0_wp - coeff) * EXP(n_bv**2 / g * (z-z0)) )
    pz = p0 * EXP( (cp_d/r_d) * LOG(1.0_wp + coeff*(EXP(-n_bv**2/g*(z-z0)) - 1.0_wp) ) )
    rhoz = pz / (r_d * tz)

  END FUNCTION rho_z_N_const_dry_t

  !=================================================================================
  !
  !   Density rho at height z of dry air with const. Brunt-Vaisala-Frequency n_bv
  !   as function of base height z0, base pressure p0 at height z0 and 
  !   base pot. temp. theta0 at height z0
  !
  !   p00:   reference pressure for pot. temperature, usually 1.0E5 Pa
  !
  !   UNTESTED YET!
  !
  !=================================================================================

  FUNCTION  rho_z_N_const_dry_theta( z, z0, n_bv, p0, p00, theta0 ) &
       RESULT (rhoz)

    USE data_constants  , ONLY :   &
         g,            & ! acceleration due to gravity
         r_d,          & ! gas constant for dry air
         cp_d            ! specific heat of dry air at constant pressure

    IMPLICIT NONE
    REAL (KIND=wp),     INTENT(in) :: z, z0, n_bv, p0, p00, theta0
    REAL (KIND=wp)     :: rhoz
    REAL (KIND=wp)     :: coeff, pz, Tz, T0

    !.. No futher subroutine- or function calls for Tz and pz for better
    !   vectorization properties!

    T0 = theta0 * (p0/p00)**(r_d/cp_d)
    coeff = g**2 / (T0 * n_bv**2 * cp_d)
    Tz = T0 * (coeff + (1.0_wp - coeff) * EXP(n_bv**2 / g * (z-z0)) )
    pz = p0 * EXP( (cp_d/r_d) * LOG(1.0_wp + coeff*(EXP(-n_bv**2/g*(z-z0)) - 1.0_wp) ) )
    rhoz = pz / (r_d * tz)

  END FUNCTION rho_z_N_const_dry_theta

  !=====================================================================================
  !=====================================================================================

  !=====================================================================================
  !
  ! Functions of Michael Baldauf for atmosphere with const. N:
  !
  !=====================================================================================
  
  REAL (KIND=wp)     FUNCTION  dry_atmosph_N_const_temp( z, n_bv, T0 )

    ! absolute temperature T in height z for a dry atmosphere with constant Brunt-Vaisala frequency n_bv
    ! with temperature T0 in z=0 

    USE data_constants  , ONLY :   &
      g,            & ! acceleration due to gravity
      r_d,          & ! gas constant for dry air
      cp_d            ! specific heat of dry air at constant pressure

    IMPLICIT NONE

    REAL (KIND=wp),     INTENT(in) :: z, n_bv, T0
    REAL (KIND=wp)     :: rcoeff

    rcoeff = T0 * n_bv**2 * cp_d / g**2 - 1.0_wp

    dry_atmosph_N_const_temp = g**2 / n_bv**2 / cp_d * ( 1.0_wp + rcoeff * EXP( n_bv**2 / g * z ) )

  END FUNCTION dry_atmosph_N_const_temp


  REAL (KIND=wp)     FUNCTION  dry_atmosph_N_const_p(z, n_bv, T0, p0 )

    ! pressure p in height z for a dry atmosphere with constant Brunt-Vaisala frequency n_bv
    ! with temperature T0 und pressure p0 in z=0 

    USE data_constants  , ONLY :   &
      g,            & ! acceleration due to gravity
      r_d,          & ! gas constant for dry air
      cp_d            ! specific heat of dry air at constant pressure

    IMPLICIT NONE

    REAL (KIND=wp),     INTENT(in) :: z, n_bv, T0, p0
    REAL (KIND=wp)     :: rcoeff

    rcoeff = T0 * n_bv**2 * cp_d / g**2 - 1.0_wp

    dry_atmosph_N_const_p =  &
      &          p0 * EXP( cp_d/r_d * ( LOG( 1.0_wp + rcoeff * EXP( n_bv**2 / g * z ) )   &
      &                               - LOG( 1.0_wp + rcoeff )    - n_bv**2 / g * z ) ) 

  END FUNCTION dry_atmosph_N_const_p


  REAL (KIND=wp)     FUNCTION  dry_atmosph_N_const_temp_rho(z, n_bv, T0, p0)

    ! density rho in height z for a dry atmosphere with constant Brunt-Vaisala frequency n_bv
    ! with temperature T0 und pressure p0 in z=0 

    USE data_constants  , ONLY :   &
      r_d            ! gas constant for dry air

    IMPLICIT NONE

    REAL (KIND=wp),     INTENT(in) :: z, n_bv, T0, p0

    dry_atmosph_N_const_temp_rho = dry_atmosph_N_const_p   (z, n_bv, T0, p0) / r_d / &
      &                            dry_atmosph_N_const_temp(z, n_bv, T0)

  END FUNCTION dry_atmosph_N_const_temp_rho

  ! END functions of Michael Baldauf
  !
  !=====================================================================================
  !=====================================================================================

  !=========================================================================
  !
  ! Output of a 3D model field with vertical dimension nk 
  ! (horizontal DIMENSION has to be ie, je) to a 3D-ASCII-file: 
  !
  ! Method: Has to be called on every processor. The global field
  !         is collected from all processors to the root processor
  !         written to a file. On systems other than the NEC, a simple
  !         text file is written, consisting of 1 header line composed
  !         of the "fieldcomment" and "fieldunit" (string provided by the user),
  !         a line containing the field dimensions (3 integer values),
  !         and a long column of data (first index varies first).
  !
  !         Filename convention: <fieldname>_YYYYMMDDHH_DDMMHHSS_<ylevtyp>.dat
  !
  !         where: YYYYMMDDHH  =  model start time
  !                DDMMHHSS    =  forecast time
  !                <fieldname> =  string provided by the user on SR call 
  !                               (e.g., "w" for vertical speed)
  !                <ylevtyp>   =  string provided by the user on SR call 
  !                               (e.g., "eta" for eta-levels)
  !
  !         To compensate for the very inefficient formatted output on a NEC,
  !         on this machine a Fortran binary file (suffix .bin instead of .dat) is written, which
  !         can be converted to the abovementioned ASCII format by using
  !         the small program "bin2ascii_convrates3d.f90" by Ulrich Blahak:
  !
  !         $> bin2ascii_convrates3d file.bin
  !
  !         This will generate the .dat-file from the .bin-file.
  !
  ! Author:  Ulrich Blahak, ulrich.blahak@dwd.de
  !
  ! Example of SR call:
  !
  !         call output3d_ascii(w(:,:,:,nnow), 51, "w", "vertical velocity", "m/s", "eta")
  !
  !=========================================================================

  SUBROUTINE output3d_ascii(zfield, nk, fieldname, fieldcomment, fieldunit, ylevtyp)

    IMPLICIT NONE

    !.. Input/Output parameters:

    INTEGER(KIND=iintegers), INTENT(in) :: nk
    REAL(KIND=wp),     INTENT(in) :: zfield(ie,je,nk)
    CHARACTER(len=*), INTENT(in) :: fieldname, fieldcomment, fieldunit, ylevtyp

    !.. Local variables:

    REAL(KIND=wp),     ALLOCATABLE :: zfield_glob(:,:)
    CHARACTER(len=250) :: ausdateiname
    INTEGER(KIND=iintegers) :: i , j, k, ios, funit

    INTEGER (KIND=iintegers)              ::    &
         mfor_s,    & ! Forecast range (seconds)
         mfor_m,    & ! Forecast range (minutes)
         mfor_h,    & ! Forecast range (hours)
         mfor_d       ! Forecast range (days)
    REAL (KIND=wp)     :: zforrange
    CHARACTER (LEN= 8)                    ::    &
         yforrange      ! Forecast range as character

#ifdef __SX__
    REAL(KIND=irealgrib), ALLOCATABLE :: zfield_glob_32bit(:,:)
    CHARACTER (len=300) :: outstring
    REAL(KIND=irealgrib) :: d32bit = 0.0_irealgrib
    CHARACTER(LEN=1) :: dimen
    CHARACTER(LEN=5) :: cnrbuf
#endif

    IF (my_cart_id == 0) THEN

      ! 1) Construct filename for the ASCII output file:
      ausdateiname = REPEAT(' ',LEN(ausdateiname))

      zforrange = REAL(ntstep, wp)*dt
      mfor_d    =  INT ( zforrange/86400.0_wp, iintegers)
      mfor_h    =  INT ((zforrange                                 &
           - REAL(mfor_d, wp)*86400.0_wp)/3600.0_wp, iintegers)
      mfor_m    =  INT ((zforrange                                 &
           - REAL(mfor_d, wp)*86400.0_wp                           &
           - REAL(mfor_h, wp)* 3600.0_wp)/  60.0_wp, iintegers)
      mfor_s    = NINT ( zforrange                                 &
           - REAL(mfor_d, wp)*86400.0_wp                           &
           - REAL(mfor_h, wp)* 3600.0_wp                           &
           - REAL(mfor_m, wp)*   60.0_wp, iintegers)
      WRITE (yforrange,'(4(I2.2))') mfor_d, mfor_h, mfor_m, mfor_s
      
      ausdateiname = TRIM(root%ydir)//'/'//TRIM(ADJUSTL(fieldname))//'_'//ydate_ini//'_'//yforrange//'_'//ylevtyp

      WRITE (*,*) ausdateiname

      CALL get_free_unit (funit)
      IF (funit == -1) THEN
        CALL model_abort (my_world_id, 10071, &
             'ERROR: problem in output3d_ascii(): no free file unit available! Abort!', &
             'output3d_ascii, opening output file')
      END IF

#ifdef __SX__

      OPEN(unit=funit, file=TRIM(ausdateiname)//'.bin', status='replace', &
           form='unformatted', iostat=ios)
      IF (ios /= 0) THEN
        CALL model_abort (my_world_id, 10072, &
             'ERROR: problem in output3d_ascii(): error opening '//TRIM(ausdateiname)//'.bin', &
             'output3d_ascii, opening output file')
      ENDIF
      
      outstring = REPEAT(' ',LEN(outstring))
      WRITE (outstring, '(a)') &
           '# '//TRIM(ADJUSTL(fieldcomment))//' ['//TRIM(ADJUSTL(fieldunit))//']'
      cnrbuf(:) = ' '
      WRITE (cnrbuf,'(i5.5)') LEN_TRIM(outstring)
      WRITE (funit) cnrbuf
      WRITE (funit) TRIM(outstring)
      
      dimen = ACHAR(3)
      WRITE (funit) dimen
      WRITE (funit) ie_tot, je_tot, nk

#else

      OPEN(unit=funit, file=TRIM(ausdateiname)//'.dat', status='replace', iostat=ios)
      IF (ios /= 0) THEN
        CALL model_abort (my_world_id, 10072, &
             'ERROR: problem in output3d_ascii(): error opening '//TRIM(ausdateiname)//'.dat', &
             'output3d_ascii, opening output file')
      ENDIF

      WRITE (funit, '(a)') &
           '# '//TRIM(ADJUSTL(fieldcomment))//' ['//TRIM(ADJUSTL(fieldunit))//']'
      WRITE (funit, '(i4,1x,i4,1x,i4)') ie_tot, je_tot, nk
#endif

    END IF

    IF (num_compute > 1) THEN
      ! Collect field on the root node and do the output:
      ALLOCATE(zfield_glob(ie_tot,je_tot))
      zfield_glob = 0.0_wp
#ifdef __SX__
      ALLOCATE(zfield_glob_32bit(ie_tot,je_tot))
      zfield_glob_32bit = 0.0_irealgrib
#endif
      DO k=1,nk
        CALL gather_field (zfield(:,:,k), ie, je, zfield_glob, ie_tot, je_tot, 0, ios)
        IF (my_cart_id == 0) THEN
#ifdef __SX__
          zfield_glob_32bit = REAL(zfield_glob(:,:), KIND=KIND(d32bit))
          WRITE (funit) zfield_glob_32bit
#else
          DO j=1,je_tot
            DO i=1,ie_tot
              IF (ABS(zfield_glob(i,j)) >= 1.0E-30_wp) THEN
                WRITE (funit, '(es12.5)') zfield_glob(i,j)
              ELSE
                WRITE (funit, '(i1)') 0_iintegers
              END IF
            END DO
          END DO
#endif
        END IF
      END DO
      DEALLOCATE(zfield_glob)
#ifdef __SX__
      DEALLOCATE(zfield_glob_32bit)
#endif
    ELSE
      IF (my_cart_id == 0) THEN
        DO k=1,nk
#ifdef __SX__
          WRITE (funit) REAL(zfield(:,:,k), KIND=KIND(d32bit))
#else
          DO j=1,je
            DO i=1,ie
              IF (ABS(zfield(i,j,k)) >= 1.0E-30_wp) THEN
                WRITE (funit, '(es12.5)') zfield(i,j,k)
              ELSE
                WRITE (funit, '(i1)') 0_iintegers
              END IF
            END DO
          END DO
#endif
        END DO
      END IF
    END IF

    IF (my_cart_id == 0) THEN
      CLOSE (funit)
      CALL release_unit(funit)
    END IF

  END SUBROUTINE output3d_ascii

!==============================================================================



END MODULE src_artifdata
