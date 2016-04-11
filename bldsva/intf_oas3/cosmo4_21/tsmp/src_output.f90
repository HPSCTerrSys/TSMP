!+ Source Module for writing Grib files
!------------------------------------------------------------------------------

MODULE src_output

!------------------------------------------------------------------------------
!
! Description:
!   This module contains subroutines necessary for writing the result data
!   of the LM. It uses also routines from the module "io_utilities".
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.2        1998/03/30 Ulrich Schaettler
!  Regarding whether digital filtering has been performed.
! 1.7        1998/07/16 Guenther Doms
!  Use of additional global fields to perform output of non-global fields.
! 1.8        1998/08/03 Ulrich Schaettler
!  Use grib parameters from module data_io.f90.
! 1.9        1998/09/16 Guenther Doms
!  Use of parameters 'nincmxt' and 'nincmxu' (replacing 'nincmxn') from
!  data module 'data_runcontrol.f90'.
! 1.10       1998/09/29 Ulrich Schaettler
!  Eliminate dependency from routine remark.
! 1.14       1998/10/26 Ulrich Schaettler
!  Changed igds to igds_out.
! 1.17       1998/11/17 Ulrich Schaettler
!  Changes for reading and writing ready files and constant variables.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.24       1999/03/01 Guenther Doms
!  Inclusion of the new prognostic 3-D array 'qi'.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and prepared use of MPE_IO;
!  Subroutine tautsp has been put to module utilities.f90
! 1.30       1999/06/24 Matthias Raschendofer
!  Use 6 additional fields: t, t_g, qv_s, tfm, tfh, tke, gz0.
!  Use 1 additional parameter form module data_runcontrol: ntke
! 1.32       1999/08/24 Guenther Doms
!  Use of postprocessing utilities from new routine 'pp_utilities'
! 1.33       1999/10/14 Guenther Doms
!  Use of postprocessing utility 'caliq' from 'pp_utilities'.
! 1.34       1999/12/10 Ulrich Schaettler
!  Use new Namelist variables and include all module procedures in this file
! 1.38       2000/04/06 Christoph Schraff
!  Correction for mean values when writing analyses (in 'output_grib_data').
! 1.39       2000/05/03 Ulrich Schaettler 
!  Changed names for variables concerned to latitude or longitude.
!  Introduced possibility for database IO.
!  Prepared output for different nests (for interactive nesting).
! 1.40       2000/05/23 Ulrich Schaettler
!  No interpolation to masspoints for u,v for p- and z-interpolation.
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists;
!  Adapted input to new organization of I/O;
!  Introduced 2D version of p- and z-interpolation
! 2.10       2001/07/24 Ulrich Schaettler
!  Corrected declaration of my_iee in subroutine output_grib_data
! 2.11       2001/09/28 Ulrich Schaettler
!  Eliminated extensive use of LEN_TRIM (which was not performant on e.g. NEC)
!  Corrected a bug when gribbing the pole of the rotation grid (igds_out(21))
! 2.14       2002/02/15 Ulrich Schaettler
!  Modifications for computation of the height of the snow-fall limit
!  Correction in the grib-coding of the upper right corner (igds_out(11))
! 2.17       2002/05/08 Ulrich Schaettler
!  Modifications to perform I/O-communications in irealgrib-format
! 2.19       2002/10/24 Ulrich Schaettler
!  Corrected bugs in writing variables from multi-layer soil model and
!  in calculating the gds-values in case of luvmasspoints=.TRUE.
! 3.5        2003/09/02 Ulrich Schaettler
!  Include output for zenith delay (routine calztd from pp_utilities)
!  Include output for ICW, ICI (integrated cloud water; integrated cloud ice)
!  Corrected output for calrelhum, calomega
! 3.6        2003/12/11 Ulrich Schaettler
!  Adaptations for multi-layer soil model
! 3.7        2004/02/18 Ulrich Schaettler
!  Adaptations for output of synthetic satellite images as GRIB fields;
!  Possibility for specifying additional GRIB tables
!  Renamed phi to rlat
! 3.8        2004/03/23 Ulrich Schaettler
!  Bug-Fix in the treatment of the linked-list for Namelist group /GRIBOUT/
! 3.13       2004/12/03 Ulrich Schaettler
!  Put KIND-parameters for Grib-library to data_parameters;
!  Changed W_ICE to W_SO_ICE (Reinhold Schrodin)
!  Adaptations for output of new variables (graupel scheme, 3D turbulence)
!  Possibility to write 3dimensional variables on p- and z-levels
!                            (Thorsten Reinhardt, Jochen Foerstner)
!  Bug correction for writing gds for u and v in case of p- or z-levels
!                            (Emanuele Zala)
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL; Implemented NL variable nunit_of_time
! 3.16       2005/07/22 Ulrich Schaettler
!  Adapted length of pathname for output-directory
!  Calculate new fields for output of radar images
! 3.18       2006/03/03 Ulrich Schaettler
!  Introduction of writing NetCDF and Restart files
!  Changed treatment of ASCII files for introducing restart possibility
!  Added output variables TQR, TQS, TQG, RELHUM_2M
!  Changes to introduce new type of coding the vertical coordinate parameters
!    (introduced by MeteoSwiss for SLEVE vertical coordinate)
!  Introduction of ldwd_grib_use: if .TRUE., special non-standard Grib Code
!    settings used at DWD are done (setting of ipds(4), tri, reference time
!    for some analysis products)
!  Introduction of possibility to write a subdomain field
!  Determination of grib record length now by return value idims_out(19) from
!    routine grbex1 (instead of idims_out(18)*iwlength)
!  New routines for NetCDF: write_nc_gdefs, write_nc_vdefs
!  New routines smooth_pmsl, smooth_geopot for extreme smoothing pmsl, geopot
!    in mountainous terrain (introduced by Bundeswehr
! 3.19       2006/04/25 Ulrich Schaettler
!  Corrections in the NetCDF output
! 3.21       2006/12/04 Ulrich Schaettler
!  changes in write_nc_gdefs to meet netCDF CF standards (B. Rockel)
!  Correction for specifying soil types for NetCDF output
!  land/sea masks in netCDF output included
!  polgam introduced
!  Save the vertical coordinate parameters in restart files (U. Schaettler)
!  Use nnow for output
!  Adaptations for Ensemble Prediction output (C. Gebhardt)
!  Additional output of several convective indices (D. Leuenberger)
!  Time integrated analysis increment fields introduced (C. Schraff)
! V3_23        2007/03/30 Jochen Foerstner, Michael Baldauf, Ulrich Schaettler
!  Corrected computation of q_sedim (must be done with qi(:,:,:,nnew)
!  Added call to SR calc_ceiling
!  Introduced idbg_level for verbosity of output
! V3_24        2007/04/26 Ulrich Schaettler
!  Eliminated nincmxu, nincmxt and introduced control as for other increments
! V3_25        2007/05/21 Ulrich Schaettler
!  Corrections for writing synthetic satellite images for MSG2
!  Modifications for writing lat/lon values to NetCDF files
! V3_26        2007/05/29 Ulrich Schaettler
!  More debug level output
! V3_27        2007/06/29 Ulrich Schaettler
!  Additional correction for flexible dt in makepds for writing analysis
! V4_1         2007/12/04 Ulrich Schaettler
!  Corrected settings of igds_out for ivctype=3
!  Bug fix for re-initializing rain rates after restart files (Uwe Boehm)
!  Introduced output for SDI (supercell detection indices)
!  Introduced additional clipping of variables, if values are around 0
!  (only for grib-output)
! V4_3         2008/02/25 Ulrich Schaettler
!  Omit grib warnings in case of climate simulations
! V4_4         2008/07/16 Ulrich Schaettler
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  Adapted interface of get_timings
! V4_5         2008/09/10 Guenther Zaengl
!  Adaptations for new reference atmosphere
! V4_8         2009/02/16 Ulrich Schaettler, Guenther Zaengl
!  Corrections for NetCDF Input and Grib Output
!  Use noutlevels as maximum number of output levels
!  Adapted interface to SR cal_conv_ind to changes in pp_utilities
!  Adapted GDS encoding and restarts for new reference atmosphere
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation
!  Bug fix for grib encoding of LM output (affects only sleve coordinate and
!  new reference atmosphere)
!  Add l_ke_in_gds to partly replace ldwd_grib_use
! V4_9         2009/07/16 Ulrich Schaettler, Burkhardt Rockel
!  Corrections in closing YUCHKDAT
!  Vectorization of p_int and z_int routines
!  Adaptations for new reference atmosphere in netCDF output
!  Include ivctype=2 and ivctype=3 in netCDF output
! V4_10        2009/09/11 MCH
!  Computation and output of BRN and HPBL
!  Added compiler directive to use option _on_adb for NEC
! V4_11        2009/11/30 Guenther Zaengl, Lucio Torrisi
!  Adaptations for output using irefatm=3 (const. Brunt-Vaisala frequency)
!  Adaptations for output of additional fluxes
! V4_12        2010/05/11 Michael Baldauf, Ulrich Schaettler
!  Introduced output for potential and (relative) Vorticity
!  Moved subroutine calc_sdi from pp_utilities to here because of formal reasons
!  Renamed t0 to t0_melt because of conflicting names
!  Added more convective indices fields for output; 
!  Introduced unit_of_time=15 for 10 minutes output (Oli Fuhrer)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler, Oliver Fuhrer
!  Eliminated tgrlat, acrlat from interface to SR potential_vorticity_rho
!  Use of nnow instead of itl (for itimelevel) in src_output for some variables
!   VORTIC_U,_V,_W, POT_VORTIC,  T_WATER, and SR caliq, calztd, calomega
! V4_17        2011/02/24 Ulrich Blahak
!  - Adapted interface of exchg_boundaries; 
!  - corrected kzdims(1:20) -> kzdims(1:24);
!  - eliminated my_peri_neigh; 
!  - implemented linear interpolation as alternative to
!    cubic tension splines in SRs z_int() and pint() -- this uses 
!    the new interpolation routine lininterp2D_xinter1D_vec() 
!    from utilities.f90 and may be activated by the new namelist
!    parameter itype_vertint in namelist group gribout (1=cubic spline, 2=linear);
!  - surface values of U, V, W, and T
!    for *linear* z_int()-interpolation now depend on lnosurffluxes_m/h (free-slip b.c or not)
!    FOR NOW: THIS IS NOT DONE FOR P-LEVELS OUTPUT AND Z-LEVELS CUBIC 
!    SPLINES TO PRESERVE "OLD" BEHAVIOUR AND TO AVOID PROBLEMS
!    WITH OVERSHOOTS IN CUBIC SPLINE INTERPOLATION.
!  - surface value of pressure for z_int()-interpolation is now
!    taken to be the surface pressure ps instead of p on the
!    lowest main levels. This changes slightly the results
!    of z-interpolated pressure near the surface.
!  - changed variable name "result"
!    to "results", since "result" is a fortran key word; 
!  - added debug output on processed variables for constant fields output.
! V4_18        2011/05/26 Ulrich Schaettler
!  Introduced conditional compilation for synthetic satellite images
!  Moved NL parameter yform_write to the group(s) /GRIBOUT/ to be able to
!     specify the output format differently for every group. 
!  Adapted NetCDF I/O to deal with 3D external parameter field for sectors of
!   the horizon (for topographical corrections) and its attributes (Anne Roches)
!  Adapted NetCDF I/O to deal with synthetic satellite data (but only MSG)
!  Introduced 4 additional fields for each group of products for syn sat data
!   (Anne Roches et al.)
!  More general if-clauses for SR caliq (Michael Baldauf)
!  Bug fixes in calls to SR calc_Theta_Tppp and potential_vorticity_rho
!   (reported by Jean-Marie Bettems)
!  Exchange of boundaries for AUMFL_S, AVMFL_S, if luvmasspoint=.TRUE.
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for NetCDF and GRIBOUT
!  Check inconsistent RTTOV- and OUTPUT-settings for synthetic satellite images
!     (Robin Faulwetter)
! V4_20        2011/08/31 Ulrich Schaettler
!  Replaced variablename namelist (which is Fortran Keyword) by outblock
!  Bug fix for computing total pressure on highest half level (J. Schmidli)
! V4_21        2011/12/06 Ulrich Blahak
!  Bugfixes p_int for the case itype_vertint=2 (linear vertical interpolation):
!   - monotonically increasing (dummy) p-values are also required below the
!     surface for routine lininterp2D_xinter1D_vec(). 
!   - error in field dimension when calling lininterp2D_xinter1D_vec
!     (-> model crashes) was corrected.
!  Bug in calling sequences of radar_lm_ray: The routine requires
!   hydrometeor densities (kg/m**3), not specific values (kg/kg), so added 
!   necessary multiplications with rho.
!  Initialized variable izerror in SR calc_sdi (Oli Fuhrer)
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
    ireals,    & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib, & ! KIND-type parameter for real variables in the grib library
    iwlength,  & ! length of an integer word in byte
    intgribf,  & ! KIND-type parameter for fortran files in the grib library
    intgribc     ! KIND-type parameter for C files in the grib library

USE data_fields, ONLY :  &
    hhl,          & ! geometrical height of model half levels
    hsurf,        & ! height of surface topography
    llandmask,    & ! landpoint mask
    rlat,         & ! geographical latitude
    rain_gsp,     & ! amount of rain from grid-scale precip. (sum)  (kg/m2)
    snow_gsp,     & ! amount of snow from grid-scale precip. (sum)  (kg/m2)
    grau_gsp,     & ! amount of graupel from grid-scale prec. (sum) (kg/m2)
    rain_con,     & ! amount of rain from convective precip. (sum)  (kg/m2)
    snow_con,     & ! amount of snow from convective precip. (sum)  (kg/m2)
    pp,           & ! deviation from the reference pressure 
    ps,           & ! surface pressure 
    dp0,          & ! pressure thickness of model layers
    rho0,         & ! base-state density
    rho,          & ! total air density
    p0,           & ! base-state pressure of full model levels
    p0hl,         & ! base-state pressure of half model levels
    t0,           & ! base state temperature
    clc_sgs,      & ! subgrid-scale stratiform cloud cover  
    clc_con,      & ! convective cloud cover
    top_con,      & ! level index of convective cloud top 
    bas_con,      & ! level index of convective cloud base
    pptens          ! pressure tendency

USE data_fields, ONLY :  &
    qv,           & ! specific humidity
    qc,           & ! specific cloud water content
    qi,           & ! specific cloud ice content
    qr,           & ! specific rain content                         (kg/kg)
    qs,           & ! specific snow content                         (kg/kg)
    qg,           & ! specific graupel content                      (kg/kg)
    qrs,          & ! specific precip. water content                (kg/kg)
    sqrtg_r_s,    & ! 1 / square root of G at scalar points         ( 1/m )
    crlat,        & ! cosine of transformed latitude
    u,            & ! zonal velocity
    v,            & ! meridional velocity
    w ,           & ! vertical velocity
    t ,           & ! temperature
    t_g,          & ! weighted surface temperature                   (  k  )
    t_2m,         & ! temperature in 2m                             (  K  )
    qv_2m,        & ! specific water vapor content                  (kg/kg)
    w_snow,       & ! water content of snow                         (m H2O)
    p_anai,       & ! deviation from the reference pressure         ( Pa  )
    qv_anai,      & ! specific water vapor content                  (kg/kg)
    qc_anai,      & ! specific cloud water content (via saturation adjustm)
    synme7,       & !
    synmsg,       & !
    fc,           & ! coriolis-parameter                            ( 1/s )
    fccos,        & ! horizontal coriolis-parameter                 ( 1/s )
    acrlat,       & ! 1 / ( crlat * radius of the earth )           ( 1/m )
    tgrlat,       & ! tangens of transformed latitude                 --
    sqrtg_r_s,    & ! 1 / square root of G at scalar points       ( 1/m )
    dzeta_dlam,   & ! d zeta / d lambda (for constant phi,    z)
                    ! at the scalar position                      ( 1   )
    dzeta_dphi      ! d zeta / d phi    (for constant lambda, z)
                    ! at the scalar position                      ( 1   )


USE data_modelconfig,ONLY : &
    czmls,        & ! depth of the main soil layers in m
    czhls,        & ! depth of the half soil layers in m
    msoilgrib,    & ! grib coded depth of main soil levels in centimeters
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    dt,           & ! long time-step
    istartpar,    & ! start index for computations in the parallel program
    jstartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    ie,           & ! number of grid points in zonal direction
    ie_tot,       & ! number of grid points in zonal direction total
    je,           & ! number of grid points in meridional direction
    je_tot,       & ! number of grid points in meridional direction total
    ie_max,       & ! Max. of ie on all processors
    je_max,       & ! Max. of je on all processors
    ke,           & ! number of grid points in vertical direction
    ke_tot,       & ! number of grid points in vertical direction total
    ke1,          & ! KE+1
    ke_soil,      & ! number of layers in the multi-layer soil model
    p0sl,         & ! reference pressure at sea level
    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    polgam,       & ! angle between the north poles of the systems
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    startlat_tot    ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)

USE data_modelconfig,ONLY : &
    hhlr ,        & ! height-coordinate refering to MSL (z=0)      half level
    vcoord,       & ! vertical coordinate of LM                    half level
    ivctype,      & ! type of vertical coordinate
    svc1,         & ! vertical decay rate for large-scale topo part of SLEVE
                    ! coordinate (in meter)
    svc2,         & ! vertical decay rate for small-scale topo part of SLEVE
                    ! coordinate (in meter)
    nfltvc,       & ! number of filter applications in topo-splitting of
                    ! SLEVE coordinate
    vcflat,       & ! vertical coordinate for changing back to z-system
    irefatm,      & ! 1: old reference atm. based on dT/dlnp = const
                    ! 2: new reference atm. with asymptotically isothermal stratosphere
                    ! 3: constant Brunt-Vaisala frequency (idealized runs only)
    t0sl,         & ! reference temperature at sea level
    dt0lp,        & ! d (t0) / d (ln p0) (for irefatm = 1)
    delta_t,      & ! temp. difference between sea level and stratosphere (for irefatm = 2)
    h_scal,       & ! scale height for irefatm = 2
    bvref,        & ! BV-freq. for irefatm=3
    klv850,       & !
    eddlon,       & ! 1 / dlon
    eddlat          ! 1 / dlat


USE data_constants, ONLY : &
    pi,           & ! circle constant
    r_d,          & ! gas constant for dry air
    r_earth,      & ! mean radius of the earth (m)
    g,            & ! gravity acceleration
    rdv,          & ! r_d / r_v
    rvd_m_o,      & ! r_v/r_d - 1
    o_m_rdv,      & ! 1 - r_d/r_v
    cpdr,         & ! 1.0 / cp_d
    t0_melt,      & ! melting temperature
    rho_w,        & ! density of liquid water
    rho_ice,      & ! density of ice          (kg/m^3)
    K_w,          & ! dielectric constant for water
    K_ice,        & ! dielectric constant for ice
    b1,           & !
    b2w,          & !
    b3,           & !
    b4w,          & !
    cp_d,         & ! specific heat capacity of dry air
    lh_v            ! latent heat of condensation

USE data_runcontrol, ONLY : &
    nlastmxu,     & ! last step when vbmax was "nullified"
    nlastmxt,     & ! last step when tmin, tmax were "nullified"
    nnew,         & ! corresponds to ntstep + 1
    nnow,         & ! corresponds to ntstep
    nstart,       & ! first time step of the forecast
    ntstep,       & ! actual time step
    ntke,         & ! actual TKE-time step, corresponds to ntstep
    nvers,        & ! version number of experiment for documentation
    l2tls,        & ! time integration by two timelevel RK-scheme (.TRUE.)
                    ! else with split-explicit scheme (only for l2tls=FALSE!)
    leps,         & ! switch ensemble mode on/off
    iepsmem,      & ! ID of ensemble member (EPS)
    iepstot,      & ! total number ensemble members (EPS)
    iepstyp,      & ! ID of ensemble generation type (EPS)
    itype_turb,   & ! type of turbulent diffusion parametrization
    lprogprec,    & ! forecast with prognostic rain and snow (qr, qs)
    lprog_qi,     & ! if .TRUE., running with cloud ice
    itype_gscp,   & ! type of grid-scale precipitation physics
    lmetr,        & ! lartif_data=.TRUE.:  with metric terms
                    !            =.FALSE.: or without metric terms
    ldfi,         & ! whether digital filtering or not
    lmulti_layer, & ! run multi-layer soil model
    nhori,        & ! number of sectors for the horizont array by the topographic
                    ! correction of the radiation
    lradtopo,     & ! if .TRUE., calculate topographic correction of radiation
    lcori_deep,   & ! if =.TRUE.: take cos(phi) coriolis terms into account
    itype_calendar,&! for specifying the calendar used
    psm0,         & ! initial value for mean surface pressure ps
    dsem0,        & ! initial value for mean dry static energy
    msem0,        & ! initial value for mean moist static energy
    kem0,         & ! initial value for mean kinetic energy
    qcm0,         & ! initial value for mean cloudwater content
    yakdat1,      & ! actual date (ydate_ini+ntstep/dt) in the form
                    ! ddmmyyyyhh (day, month, year, hour)
    luse_rttov      ! if .true. calculate synthetic satellite images

USE data_runcontrol, ONLY : &
    ltime,        & ! detailed timings of the program are given
    idbg_level,   & ! to control the verbosity of debug output
    ldebug_io ,   & ! if .TRUE., debug output for I/O
    lprintdeb_all,& ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output
    luse_rttov,   & ! if rttov-library is used
    lartif_data,  & ! forecast with self-defined artificial data
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim,        & ! 2 dimensional runs
    nbl_exchg       !


USE data_parallel,      ONLY :  &
    lasync_io,      & ! if .TRUE.: the model runs with extra PEs for
                      ! asynchronous IO
    my_cart_id,     & ! rank of this subdomain in the cartesian communicator
    num_compute,    & ! number of compute PEs
    icomm_cart,     & ! communicator that belongs to the cartesian grid
    imp_reals,      & ! determines the correct REAL type used in the model
                      ! for MPI
    imp_grib,       & ! determines the REAL type for the GRIB library
    imp_integers,   & ! determines the correct INTEGER type used in the model
                      ! for MPI
    nboundlines,    & ! number of boundary lines of the domain for which
                      ! no forecast is computed = overlapping boundary
                      ! lines of the subdomains
    my_cart_neigh,  & ! neighbors of this subdomain in the cartesian grid
    iexch_req,      & ! stores the sends requests for the neighbor-exchange
                      ! that can be used by MPI_WAIT to identify the send
    ldatatypes,     & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,  & ! if .TRUE.: use additional barriers for determining the
                      ! load-imbalance
    ncomm_type,     & ! type of communication
    sendbuf,        & ! sending buffer for boundary exchange:
                      ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen,    & ! length of one column of sendbuf
    nproc, realbuf, intbuf

USE data_io,        ONLY : &
  nhour_restart,     & ! start-, stop-, inc of writing restart files (tstep)
  nzmxid,            & ! maximum number of NetCDF variabe IDs
  ydate_ini,         & ! start of the forecast 
  ymode_write,       & ! mode for opening the (write) Grib files
  yform_read,        & ! format of the (read) files
  ntrans_out,        & ! Unit Number for writing ready-Files during output
  nuchkdat,          & ! Unit number for checking the I/O data
  yuchkdat,          & ! checking the I/O data
  ytrans_out,        & ! directory for writing ready files
  nsma_stat,         & ! status for soil moisture analysis
  npds,              & ! Dimension for product definition section (pds)
  ngds,              & ! Dimension for grid description section (gds)
  nbms,              & ! Dimension for bit map section (bms)
  nbds,              & ! Dimension for binary data section
  ndsup,             & ! Dimension for dsup
  ndims,             & ! Dimension for idims (contains all dimensions)
  lfd,               & !
  lbm,               & !
  lds,               & !
  ntrip,             & ! maximum number of timing triples
  noutlevels,        & ! maximum actual existing number of output levels
  nvar,              & ! maximum number of variables in LM variable table
  max_gribtabs,      & ! maximum number of GRIB tables in LM variable table
  iednr,             & ! grib edition number
  undefgrib,         & ! value for "undefined" in the grib routines
  undefncdf,         & ! value for "undefined" in the netcdf routines
  undef,             & ! the same as undefgrib but with other KIND-Parameter
  ncenter,           & ! originating center identification
  lst_gribtabs,      & ! IDs of GRIB tables use
  nprocess_ini_in,   & ! process gener. identification for initial (analysis)
  nprocess_bd_in       ! and for boundary (forecasts) data from input data

USE data_io,        ONLY : &
! Global arrays
  iblock,            & ! array for gribed data
  idims_out,         & ! array for all dimensions
  ibmap,             & ! array for
  ipds,              & ! product definition section
  igds_out,          & ! grid description section
  ibms,              & ! bit map section
  ibds,              & ! binary data section
  dsup,              & ! Parameter for grib routines
  ds_grib,           & ! array for unpacked data
  ds_real,           & ! array for unpacked data

! Global types
  pp_nl,             & ! structure for gribout namelist
  var                  ! array for LM variable table

USE data_io,        ONLY : &
  ldwd_grib_use,     & ! use some DWD specific Grib settings
  l_ke_in_gds,       & ! explicit GDS entry for number of model levels
  lbdclim,           & ! boundary data in climate model     ! PIK  (D.Hauffe)
                       !  (in climate mode also some external parameters have
                       !  to be updated, which are held constant in forecast
                       !  mode; e.g. plant cover, root depth)
  idims_id_out,      & ! array for the IDs of the dimensions of netCDF 
                       ! formatted output
  yncglob_institution,   & ! originating center name
  yncglob_title,         & ! title string for the output
  yncglob_source,        & ! program name and version
  yncglob_project_id,    & ! identification of the project of simulation
  yncglob_experiment_id, & ! identification of the experiment of simulation
  yncglob_contact,       & ! contact e.g. email address
  yncglob_references,    & ! URL, report etc.
  ncglob_realization       ! number of the realization of the experiment


#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
USE data_satellites,          ONLY :  &
    sat_compute, num_sensors
#endif

!------------------------------------------------------------------------------

USE utilities,                ONLY :  &
    smoother, dfilt4, dfilt8, tautsp, tautsp2D,                     &
    phirot2phi, rlarot2rla, get_utc_date, &
    lininterp2D_xinter1D_vec

USE pp_utilities,             ONLY :  &
    calpmsl,    calprsum,   caltopdc,   calhzero,   calsnowlmt,     &
    calcldepth, calclmod,   calomega,   calrelhum,  caliq,  calztd, &
    radar_lm_ray, cal_conv_ind, calc_ceiling,                       &
    calc_bulk_richardson, calc_pbl_brn,                             &
    potential_vorticity_rho

USE numeric_utilities,             ONLY :  &
    curl, calc_Theta_Tppp, mean_over_box, mean_cov_over_box, vert_avg

USE environment,              ONLY :  &
    model_abort, get_free_unit, release_unit, exchg_boundaries

USE parallel_utilities,       ONLY :  &
    gather_values, combine_subarrays, distribute_values, gather_field

USE io_utilities,             ONLY :  &
    open_file, close_file, write_grib, write_netcdf, write_restart, &
    make_fn, check_record

USE time_utilities,     ONLY: get_timings, i_computations_O, i_gather_data, &
                                           i_write_data

USE src_artifdata,      ONLY:  lnosurffluxes_m, lnosurffluxes_h

!------------------------------------------------------------------------------

#ifdef NETCDF
USE netcdf,           ONLY :   &
  nf90_def_dim,            &
  nf90_def_var,            &
  nf90_enddef,             &
  nf90_redef,              &
  nf90_put_att,            &
  nf90_put_var,            &
  nf90_noerr,              &
  nf90_strerror,           &
  NF90_CHAR,               &
  NF90_DOUBLE,             &
  NF90_FLOAT,              &
  NF90_GLOBAL,             &
  NF90_UNLIMITED
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! string variable to hold grid information
  CHARACTER (LEN=200)     grid_mapping

! for smoothing fi and pmsl, the global hsurf-field is needed
  REAL (KIND=ireals), ALLOCATABLE  :: hsurf_tot(:,:)

! for NetCDF output of MSG-variables
INTEGER (KIND=iintegers), PARAMETER  ::          &
  nmsgchan = 8

!==============================================================================
! Module procedures
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in src_output for initializing the output organization
!------------------------------------------------------------------------------

SUBROUTINE init_output (root, n_num_now)

!------------------------------------------------------------------------------
!
! Description:
!  The routine init_output initializes organizational variables of the 
!  output routines dealing with the grib code (such as the dimensions of
!  the different grib code sections). Also the grid description section
!  is initialized (except the location of the lower left grid point, because
!  it depends on the variable (U,V or other)).
!  
! Method:
!
!------------------------------------------------------------------------------

! Pointers with intent(in):
TYPE (pp_nl), POINTER                 :: root

INTEGER (KIND=iintegers), INTENT(IN)  ::    &
  n_num_now        ! current nest number

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)     :: i,i1,i2,i3,k,n, niostat, ierrstat, itable,  &
                                nuedat, nzrecords, kbot, ktop, izerror,     &
                                iorg_data(3,0:num_compute-1), izdebug

LOGICAL                      :: lzcheck, lzopen, lzwriteready

CHARACTER  (LEN=260)         :: yname
CHARACTER  (LEN=25)          :: yroutine
CHARACTER  (LEN=80)          :: yzerrmsg
CHARACTER  (LEN= 3)          :: yzhead
CHARACTER  (LEN= 3)          :: c_nnum

! Local arrays:
REAL (KIND=ireals)      ::     &
  zvarlev   (ie,je,0:noutlevels),      & ! variable for vertical interpolation
  zprocarray_real(ie_max,je_max,num_compute)
                                 ! has a two dimensional array for all 
                                 ! nodes for gathering values from the nodes

REAL (KIND=irealgrib)   ::     &
  zprocarray_grib(ie_max,je_max,num_compute) 
                                 ! has a two dimensional array for all
                                 ! nodes for gathering values from the nodes

INTEGER (KIND=iintegers)::     &
  ivar_id(nzmxid)                ! NetCDF-ID of each variable in the list

! Local Pointers:
TYPE (pp_nl), POINTER   :: now

!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

yroutine = 'init_output'

! Initialize, whether additional debug output shall be done
IF (ldebug_io) THEN
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF
ELSE
  izdebug = 0
ENDIF

! Set lfd, lds and lbm
lds = ie_tot * je_tot
lbm = 1875
lfd = lds * 4 / iwlength + 2000   ! the "4" means 4 bytes per word to allow for 24/32 bit grib output
                                  ! the "2000" just is safety
! Allocate GRIB arrays
ALLOCATE (iblock(lfd), ibmap(lbm), STAT=ierrstat)
ALLOCATE (ds_real(lds), ds_grib(lds), dsup(ndsup), STAT=ierrstat)

! Initializations for the grib library
!  moving arraydimensions into idims
!  declaration dimensions
idims_out( 1) = npds
idims_out( 2) = ngds
idims_out( 3) = nbms
idims_out( 4) = nbds
idims_out( 5) = lbm
idims_out( 6) = ndsup
idims_out( 8) = lfd

!  real dimensions
idims_out(11) = 47
IF (ldwd_grib_use .AND. (.NOT.l_ke_in_gds)) THEN
  ! old style of coding the vertical coordinate parameters
  IF (ivctype == 3) THEN   ! SLEVE coordinates
    idims_out(12) = 25 + 4 + ke1 + 4
  ELSE                     ! NON SLEVE coordinates
    ! should logically be ke1+1, but is kept for backward compatibility
    idims_out(12) = 25 + 4 + ke1
  ENDIF
  IF (irefatm == 2) idims_out(12) = 25 + 4 + ke1 + 6
  IF (irefatm == 3) idims_out(12) = 25 + 4 + ke1 + 5
ELSE 
  ! new style of coding the vertical coordinate parameters
  IF (ivctype == 3) THEN   ! SLEVE coordinates
    idims_out(12) = 25 + 6 + ke1 + 3
  ELSE                     ! NON SLEVE coordinates
    idims_out(12) = 25 + 6 + ke1
  ENDIF
  IF (irefatm == 2) idims_out(12) = 25 + 6 + ke1 + 5
  IF (irefatm == 3) idims_out(12) = 25 + 6 + ke1 + 4
ENDIF
idims_out(13) = 3
idims_out(14) = 5
idims_out(16) = 0

! idims_out(7,15,17) depend on the special namelist group and are set later

! Initializations for different nests
c_nnum   = '   '
IF ( n_num_now > 1 ) THEN
  IF ( n_num_now < 10 ) THEN
    WRITE ( c_nnum, '( A1, I1, A1 )' ) '.', n_num_now, ' '
  ELSE
    WRITE ( c_nnum, '( A1, I2.2 )'   ) '.', n_num_now
  END IF
END IF

! Initialization for ivar_id (should be set also for non netcdf output)
ivar_id(:) = 0_iintegers

!------------------------------------------------------------------------------
! Section 2: Set the grid description section
!------------------------------------------------------------------------------

CALL makegds

!------------------------------------------------------------------------------
! Section 3: Open YUCHKDAT, if it has not been opened already
!------------------------------------------------------------------------------

! check whether a gribout namelists wants to write to YUCHKDAT
! define startpoint of the list
now => root
lzcheck = .FALSE.

gribout_loop: DO
  IF ( now%n_num == n_num_now ) THEN
    IF (now%lcheck .EQV. .TRUE.) THEN
      lzcheck = .TRUE.
    ENDIF
  ENDIF

  IF (ASSOCIATED(now%next)) THEN
    now => now%next
  ELSE
    EXIT gribout_loop
  ENDIF
ENDDO gribout_loop

! open file YUCHKDAT
IF ( (lzcheck .EQV. .TRUE.) .AND. (my_cart_id == 0) ) THEN
  OPEN(nuchkdat, FILE=yuchkdat, FORM=  'FORMATTED', STATUS='UNKNOWN',  &
                 POSITION='APPEND', IOSTAT=niostat)
  IF(niostat /= 0) THEN
    yzerrmsg = ' ERROR    *** Error while opening file YUCHKDAT *** '
    ierrstat = 2005
    CALL model_abort (my_cart_id, ierrstat, yzerrmsg, yroutine)
  ENDIF
ENDIF

!------------------------------------------------------------------------------
! Section 4: Write file with constant fields, if required
!------------------------------------------------------------------------------

! All namelist groups are inquired, whether constant fields shall be written.
! If no further output is done in step 0, a ready-file is written (checked
! with lzwriteready)

now => root
lzwriteready = .FALSE.

gribout_loop_2: DO
  IF ( now%n_num == n_num_now ) THEN
    IF (now%lwrite_const .EQV. .TRUE.) THEN

      ! set the namelist-dependent idims_out values
      idims_out( 7) = now%ie_out_tot * now%je_out_tot
      idims_out(15) = now%ie_out_tot * now%je_out_tot
      idims_out(17) = now%ie_out_tot * now%je_out_tot

      ! Create the filename
      yzhead = 'lf'//now%ydomain
      CALL make_fn (yzhead, ydate_ini, now%ytunit,'c', ntstep, dt,          &
                    now%lhour, now%ydir, yname, izdebug, izerror)

      ! Add nest specification to filename
      IF (now%n_num > 1) THEN
        yname = yname(1:LEN_TRIM(yname)) // c_nnum
      ENDIF

      ! In case of netcdf, add extension '.nc' to filename
      IF (now%yform_write == 'ncdf') THEN
        yname = yname(1:LEN_TRIM(yname)) // '.nc'
      ENDIF

      ! Add optional suffix to filename
      IF ( LEN_TRIM(now%ysuffix) /= 0 ) THEN
        yname = yname(1:LEN_TRIM(yname)) // TRIM(now%ysuffix)
      ENDIF

      ! open file
      IF (now%yform_write == 'bina') THEN
        ! get a free unit-number for Fortran OPEN
        CALL get_free_unit (nuedat)
      ENDIF
      CALL open_file (nuedat, yname, ymode_write, now%yform_write, icomm_cart,  &
                      my_cart_id, num_compute, lasync_io, idbg_level,           &
                      yzerrmsg, izerror)
      IF (izerror /= 0) THEN
         CALL model_abort (my_cart_id, 2031, yzerrmsg, yroutine)
      ENDIF

      ! gridpoints, simple packing, floating point data
      ibds(2)   = 0

      ! nrbit, number of bits
      ibds(5)   = now%nrbit

      ! no bitmap
      ibms(3)   = -2

#ifdef NETCDF
      IF (now%yform_write == 'ncdf') THEN
        ! Write global headers for netcdf file
        CALL write_nc_gdefs (nuedat, now, icomm_cart, num_compute, 'c',  &
                             yzerrmsg, izerror)
        IF (izerror /= 0) THEN
          CALL model_abort (my_cart_id, 8052, yzerrmsg, yroutine)
        ENDIF

        CALL write_nc_vdefs (nuedat, now%nyvar_c, now%ilist_c, ivar_id,  &
                             now%luvmasspoint, icomm_cart, num_compute,  &
                             'c', yzerrmsg, izerror)
        IF (izerror /= 0) THEN
          CALL model_abort (my_cart_id, 8053, yzerrmsg, yroutine)
        ENDIF
      ENDIF
#endif

      ! Write the headline in YUCHKDAT for the file with constant variables
      IF ( (now%lcheck) .AND. (my_cart_id == 0) ) THEN
        WRITE (nuchkdat,'(A,I7)')                                           &
                            'Check the constant data: '
        WRITE (nuchkdat,'(A,A)')                                            &
           '    File:   ',yname
        WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                                 &
           '    ie_tot =',ie_tot,'   je_tot =',je_tot,'   ke_tot =',ke_tot
        WRITE (nuchkdat,'(A)') '    '
        WRITE (nuchkdat,'(A,A)')                                            &
          '     var       ee    lev         min      ',                     &
          'imin   jmin          max      imax   jmax         mean  '
      ENDIF

      ! loop over all constant variables that should be written
      nzrecords = 0
      DO n = 1, now%nyvar_c
        ! location in the variable table
        i1 = now%ilist_c(1,n)
        i2 = now%ilist_c(2,n)
        i3 = now%ilist_c(3,n)

        IF (izdebug >= 5) THEN
          WRITE (*,'(3a,i4,a)') '  src_output: processing ', &
               TRIM(ADJUSTL(var(i1,i2,i3)%name)),' on PE ',my_cart_id,' ...'
        END IF

        SELECT CASE (var(i1,i2,i3)%rank)
        CASE(3)
          kbot = LBOUND(var(i1,i2,i3)%p3,3)
          ktop = UBOUND(var(i1,i2,i3)%p3,3)
          zvarlev(:,:,kbot:ktop) = var(i1,i2,i3)%p3(:,:,kbot:ktop)


          DO k=kbot,ktop
            nzrecords = nzrecords + 1
            CALL output_data (nuedat, nzrecords, i1,i2,i3, k, ktop,           &
                      zvarlev(1:ie,1:je,k), zprocarray_grib, zprocarray_real, &
                      now, .FALSE., 'c', 0.0_ireals, .FALSE., ivar_id(n),     &
                      iorg_data)
          ENDDO

        CASE(2)
          IF (now%yvarc(n)(1:LEN_TRIM(now%yvarc(n))) == 'FIS' ) THEN
            zvarlev(1:ie,1:je,1) = hsurf(1:ie,1:je) * g
          ELSE
            zvarlev(1:ie,1:je,1) = var(i1,i2,i3)%p2(1:ie,1:je)
          ENDIF
          nzrecords   = nzrecords+1
          CALL output_data (nuedat, nzrecords, i1,i2,i3, 1, 1,                &
               zvarlev(:,:,1), zprocarray_grib, zprocarray_real, now,         &
               .FALSE., 'c', 0.0_ireals, .FALSE., ivar_id(n), iorg_data)
        END SELECT
      ENDDO

      ! Flush the output buffers and close the file
      CALL output_data (nuedat, -1, -1,-1,-1, -1, -1,                         &
                  zvarlev(:,:,1), zprocarray_grib, zprocarray_real, now,      &
                  .TRUE., 'c', 0.0_ireals, .FALSE., -1, iorg_data)
      CALL close_file (nuedat, now%yform_write, icomm_cart, my_cart_id,       &
                       num_compute, lasync_io, idbg_level, yzerrmsg, izerror)
      IF (izerror /= 0) THEN
         CALL model_abort (my_cart_id, 2032, yzerrmsg, yroutine)
      ENDIF

      IF (now%yform_write == 'bina') THEN
        ! release the unit-number again
        CALL release_unit (nuedat)
      ENDIF

      ! Indicate that a ready file is needed
      lzwriteready = .TRUE.

      ! Write a blank line to YUCHKDAT
      IF ( (now%lcheck) .AND. (my_cart_id == 0) ) THEN
        WRITE (nuchkdat,'(A)') '   '
        WRITE (nuchkdat,'(A)') '   '
      ENDIF
    ENDIF
  ENDIF

  IF (ASSOCIATED(now%next)) THEN
    now => now%next
  ELSE
    EXIT gribout_loop_2
  ENDIF
ENDDO gribout_loop_2

! Check for further output in step 0
IF ( (my_cart_id == 0) .AND. (lzwriteready) ) THEN
  now => root
  gribout_loop_3: DO
    IF ( now%n_num == n_num_now ) THEN
      IF (now%ngrib(now%nextstep) == 0) THEN
        ! further output will be written: no ready file is necessary here
        lzwriteready = .FALSE.
      ENDIF
    ENDIF

    ! Add nest specification to filename
    IF (now%n_num > 1) THEN
      yname = yname(1:LEN_TRIM(yname)) // c_nnum
    ENDIF

    ! Add optional suffix to filename
    IF ( LEN_TRIM(now%ysuffix) /= 0 ) THEN
      yname = yname(1:LEN_TRIM(yname)) // TRIM(now%ysuffix)
    ENDIF

    IF (ASSOCIATED(now%next)) THEN
      now => now%next
    ELSE
      EXIT gribout_loop_3
    ENDIF
  ENDDO gribout_loop_3

  IF ( (lzwriteready) .AND. (ytrans_out /= '   ') ) THEN
    ! Create the filename LMF_forecasttime
    yzhead = 'LMF'
    CALL make_fn (yzhead, yakdat1, 'f',' ', ntstep, dt, .TRUE.,     &
                  ytrans_out, yname, izdebug, izerror)

    ! Write the file
    OPEN  (ntrans_out, FILE=yname, FORM='FORMATTED')
    WRITE (ntrans_out, '(A)') 'ready'
    CLOSE (ntrans_out)
  ENDIF
ENDIF

! close file nuchkdat
IF ( (now%lcheck) .AND. (my_cart_id == 0) ) THEN
  CLOSE (nuchkdat, STATUS='KEEP')
ENDIF

! Deallocate arrays for IO
DEALLOCATE (iblock, ibmap, ds_real, ds_grib, dsup)

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE init_output

!==============================================================================
!+ Module procedure in src_output for organizing the output
!------------------------------------------------------------------------------

SUBROUTINE organize_output (outblock, yextension, numlist, ylist, ilist)

!------------------------------------------------------------------------------
!
! Description:
!  The routine organize_output is called for every namelist output group and
!  for every of the three output lists (model variables, pressure level 
!  variables and height level variables). In case of pressure level variables 
!  or height level variables the routines p_int and z_int, resp., are called 
!  for vertical interpolation.
!  
!  Parallelization for the output is by layers that should be written
!  to the grib file. Every PE gets a layer and packs it into grib format
!  (in routine output_data).
!
! Method:
!  - Initializations (for the grib library)
!  - Opening the output grib file
!  - Scanning through the list (loop over all variables)
!  - Closing the output grib file
!
! Output files:
!  Output grib files for model variables (without extension), for 
!  pressure level variables (with extension 'p') and for height level
!  variables (with extension 'z').
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
TYPE(pp_nl),              INTENT(IN)     ::    &
  outblock         ! pointer to the namelist group

! Scalar arguments with intent(in):
CHARACTER (LEN=1),        INTENT(IN)     ::    &
  yextension       ! indicates model variables (''), p-('p') or z-levels ('z')

INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  numlist          ! number of elements in ylist

CHARACTER (LEN=10),       INTENT(IN)     ::    &
  ylist(numlist)   ! list of variables for output

INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  ilist(3,numlist) ! number of elements in ylist

!------------------------------------------------------------------------------
!
! Local scalars:
  INTEGER (KIND=iintegers) :: i1,i2,i3, i,j,k, n, ktop, kbot, itl, nentry,  &
                              itimelevel, isens, iorg_data(3,0:num_compute-1)
  INTEGER (KIND=iintegers) :: klev, nuedat, niostat, ierrstat, izlen, nzrecords
  CHARACTER (LEN=250)      :: yname
  CHARACTER (LEN= 25)      :: yroutine
  CHARACTER (LEN= 80)      :: yerrmsg
  CHARACTER (LEN=  3)      :: c_nnum
  CHARACTER (LEN= 10)      :: yzdat1
  CHARACTER (LEN= 22)      :: yzdat2


! Local arrays:
REAL (KIND=irealgrib)   ::        &
  zprocarray_grib(ie_max,je_max,num_compute) 
                                 ! has a two dimensional array for all 
                                 ! nodes for gathering values from the nodes

REAL (KIND=ireals)   ::     &
  zvarlev  (ie,je,0:noutlevels),&! variable for vertical interpolation
  zprocarray_real(ie_max,je_max,num_compute), &
                                 ! has a two dimensional array for all 
                                 ! nodes for gathering values from the nodes
  slev     (0:noutlevels)        ! stores the z- or the p-levels

REAL (KIND=ireals)   ::     & !
  zenith_t (ie,je),         & ! Arrays for computing the zenith total (dry,
  zenith_w (ie,je),         & ! hydrostatic) delay
  zenith_h (ie,je),         & ! 
  zcape_mu (ie,je),         & ! Arrays for most unstable CAPE
  zcin_mu  (ie,je),         & ! ... and CIN
  zcape_ml (ie,je),         & ! Arrays for mixed layer CAPE
  zcin_ml  (ie,je),         & ! ... and CIN
  zcape_3km(ie,je),         & ! ... and CAPE 3KM
  zlcl_ml  (ie,je),         & ! ... and LCL
  zlfc_ml  (ie,je),         & ! ... and LFC    
  zhelp2d  (ie,je),         & !
  zbrn  (ie,je,ke),         & ! Array for Bulk Richardson Number
  zhelp1(ie,je,ke),         & !
  zhelp2(ie,je,ke),         & !
  zhelp3(ie,je,ke),         & !
  zhelp4(ie,je,ke)            !

REAL (KIND=ireals), PARAMETER   ::     & !
  missing_value = -999.9_ireals

REAL (KIND=ireals)              ::     & !
  zacthour

LOGICAL              ::     &
  lzenith,                  & ! indicates whether to compute zenith delay
  lzrestart,                & ! indicates whether restart-files are written
  lconvind_mu,              & ! indicates whether CAPE_MU/CIN_MU has been computed
  lconvind_ml,              & ! indicates whether CAPE_ML/CIN_ML has been computed
  l_brn,                    & ! indicates whether BULK RICHARDSON NUMBER has been computed
  lwarn,                    & ! to indicate whether a SR issues some warnings
  l_calc_sdi_already_comp     ! prevents from a double computing of Subr. 'calc_sdi'

INTEGER (KIND=iintegers) :: kzdims(24)

INTEGER (KIND=iintegers) :: &
  ivar_id(nzmxid), & ! NetCDF-ID of each variable in the output list
  nzjulianday, izvctype_write, izdebug, izerror

CHARACTER (LEN=3)        :: &
  yzhead           ! characterizes the special kind of the data

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!  Section 1a: Initializations
!------------------------------------------------------------------------------

  yroutine = 'organize_output'
  ierrstat = 0
  izerror  = 0
  lzenith  = .FALSE.
  lconvind_mu  = .FALSE.
  lconvind_ml  = .FALSE.
  l_brn = .FALSE.
  lwarn = .FALSE.
  l_calc_sdi_already_comp = .FALSE.

  ! Initialize, whether additional debug output shall be done
  IF (ldebug_io) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

  IF ((yextension == 'o') .OR. (yextension == 'n')) THEN
    lzrestart = .TRUE.
  ELSE
    lzrestart = .FALSE.
  ENDIF

  ! Set lfd, lds and lbm
  lds = ie_tot * je_tot
  lbm = 1875
  lfd = lds * 4 / iwlength + 2000   ! the "4" means 4 bytes per word to allow for 24/32 bit grib output
                                    ! the "2000" just is safety

  ! Initializations for the grib library
  !  moving arraydimensions into idims
  !  declaration dimensions
  idims_out( 1) = npds
  idims_out( 2) = ngds
  idims_out( 3) = nbms
  idims_out( 4) = nbds
  idims_out( 5) = lbm
  idims_out( 6) = ndsup
  idims_out( 7) = outblock%ie_out_tot * outblock%je_out_tot
  idims_out( 8) = lfd
  
  !  real dimensions
  idims_out(11) = 47
  IF (ldwd_grib_use .AND. (.NOT.l_ke_in_gds)) THEN
    ! old style of coding the vertical coordinate parameters
    IF (ivctype == 3) THEN   ! SLEVE coordinates
      idims_out(12) = 25 + 4 + ke1 + 4
    ELSE                     ! NON SLEVE coordinates
      ! should logically be ke1+1, but is kept for backward compatibility
      idims_out(12) = 25 + 4 + ke1
    ENDIF
    IF (irefatm == 2) idims_out(12) = 25 + 4 + ke1 + 6
    IF (irefatm == 3) idims_out(12) = 25 + 4 + ke1 + 5
  ELSE
    ! new style of coding the vertical coordinate parameters
    IF (ivctype == 3) THEN   ! SLEVE coordinates
      idims_out(12) = 25 + 6 + ke1 + 3
    ELSE                     ! NON SLEVE coordinates
      idims_out(12) = 25 + 6 + ke1
    ENDIF
    IF (irefatm == 2) idims_out(12) = 25 + 6 + ke1 + 5
    IF (irefatm == 3) idims_out(12) = 25 + 6 + ke1 + 4
  ENDIF
  idims_out(13) = 3
  idims_out(14) = 5
  idims_out(15) = outblock%ie_out_tot * outblock%je_out_tot
  idims_out(16) = 0
  idims_out(17) = outblock%ie_out_tot * outblock%je_out_tot

  ! Allocate GRIB arrays
  ALLOCATE (iblock(lfd), ibmap(lbm), STAT=ierrstat)
  ALLOCATE (ds_real(lds), ds_grib(lds), dsup(ndsup)   , STAT=ierrstat)

  c_nnum   = '   '
  IF ( outblock%n_num > 1 ) THEN
    IF ( outblock%n_num < 10 ) THEN
      WRITE ( c_nnum, '( A1, I1, A1 )' ) '.', outblock%n_num, ' '
    ELSE
      WRITE ( c_nnum, '( A1, I2.2 )'   ) '.', outblock%n_num
    END IF
  END IF

  slev(:)  = 0.0_ireals

  ! gridpoints, simple packing, floating point data
  ibds(2)   = 0
 
  ! nrbit, number of bits
  ibds(5)   = outblock%nrbit

  ! no bitmap
  ibms(3)   = -2

  ! determine the timelevel for output
  IF (lzrestart) THEN
    IF (yextension == 'o') THEN
      IF (.NOT. l2tls) THEN
        itimelevel = nnow   ! for leapfrog
      ELSE
        itimelevel = nnew   ! for Runge-Kutta
      ENDIF
    ELSEIF (yextension == 'n') THEN
      itimelevel = nnew
    ENDIF
  ELSE
    ! use nnow for output (this was nnew before)
    itimelevel = nnow
  ENDIF

  ! Initialization for ivar_id (should be set also for non netcdf output)
  ivar_id(:) = 0_iintegers

!------------------------------------------------------------------------------
!  Section 1b: Gather hsurf field to all PEs, if necessary
!------------------------------------------------------------------------------

  IF (outblock%l_fi_ps_smooth) THEN
    ALLOCATE (hsurf_tot(ie_tot,je_tot), STAT=ierrstat)

    IF (num_compute > 1) THEN
      CALL gather_field (hsurf, ie,je, hsurf_tot, ie_tot,je_tot, -1, ierrstat)
    ELSE
      hsurf_tot(:,:) = hsurf(:,:)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Open the grib file
!------------------------------------------------------------------------------

  ! creating filename
  IF (outblock%lanalysis) THEN
    yzhead = 'la'//outblock%ydomain
  ELSE
    IF (lzrestart) THEN
      yzhead = 'lr'//outblock%ydomain
    ELSE
      yzhead = 'lf'//outblock%ydomain
    ENDIF
  ENDIF

  IF ((outblock%yform_write == 'bina') .AND. (yextension=='o' .OR. yextension=='n')) THEN
    ! The date of the next time step has to be determined to get the proper
    ! file name also for ytunit='d'
    CALL get_utc_date(ntstep+1, ydate_ini, dt, itype_calendar, yzdat1,       &
                      yzdat2, nzjulianday, zacthour)
    CALL make_fn (yzhead, yzdat1, outblock%ytunit, yextension, ntstep+1, dt, &
                  outblock%lhour, outblock%ydir, yname, izdebug, ierrstat)
  ELSE
    IF (izdebug > 10) THEN
      PRINT *, ' calling make_fn with date/unit: ', yakdat1, ' ', outblock%ytunit
    ENDIF
    CALL make_fn (yzhead, yakdat1, outblock%ytunit, yextension, ntstep, dt,  &
                  outblock%lhour, outblock%ydir, yname, izdebug, ierrstat)
  ENDIF

  ! Add nest specification to filename
  IF (outblock%n_num > 1) THEN
    yname = yname(1:LEN_TRIM(yname)) // c_nnum
  ENDIF

  ! In case of netcdf, add extension '.nc' to filename
  IF (outblock%yform_write == 'ncdf') THEN
    yname = yname(1:LEN_TRIM(yname)) // '.nc'
  ENDIF

  ! Add optional suffix to filename
  IF ( LEN_TRIM(outblock%ysuffix) /= 0 ) THEN
    yname = yname(1:LEN_TRIM(yname)) // TRIM(outblock%ysuffix)
  ENDIF

  IF (outblock%yform_write == 'bina') THEN
    ! get a free unit-number for Fortran OPEN
    CALL get_free_unit (nuedat)
  ENDIF
  CALL open_file (nuedat, yname, ymode_write, TRIM(outblock%yform_write),      &
                  icomm_cart, my_cart_id, num_compute, lasync_io, idbg_level,  &
                  yerrmsg, ierrstat)
  IF (ierrstat /= 0) THEN
     CALL model_abort (my_cart_id, 2033, yerrmsg, yroutine)
  ENDIF
  IF ( (outblock%yform_write == 'bina') .AND. (my_cart_id == 0) .AND. (yextension == 'o')) THEN
    ! write the initial values for the meanvalues
    WRITE (nuedat,IOSTAT=niostat) psm0, dsem0, msem0, kem0, qcm0
    ! write the vertical coordinate parameters
    IF     (irefatm == 1) THEN
      izvctype_write = ivctype
      WRITE (nuedat,IOSTAT=niostat) izvctype_write, p0sl, t0sl, dt0lp, vcflat, vcoord
    ELSEIF (irefatm == 2) THEN
      izvctype_write = ivctype+100
      WRITE (nuedat,IOSTAT=niostat) izvctype_write, p0sl, t0sl, dt0lp, vcflat, vcoord
      WRITE (nuedat,IOSTAT=niostat) delta_t, h_scal
    ELSEIF (irefatm == 3) THEN
      izvctype_write = ivctype+200
      WRITE (nuedat,IOSTAT=niostat) izvctype_write, p0sl, t0sl, dt0lp, vcflat, vcoord
      WRITE (nuedat,IOSTAT=niostat) bvref
    ENDIF
  ENDIF

#ifdef NETCDF
  IF (outblock%yform_write == 'ncdf') THEN
    ! Write global headers for netcdf file
    CALL write_nc_gdefs (nuedat, outblock, icomm_cart, num_compute,       &
                         yextension, yerrmsg, ierrstat)
    IF (ierrstat /= 0) THEN
      CALL model_abort (my_cart_id, 8052, yerrmsg, yroutine)
    ENDIF

    CALL write_nc_vdefs (nuedat, numlist, ilist, ivar_id,                 &
                         outblock%luvmasspoint, icomm_cart, num_compute,  &
                         yextension, yerrmsg, ierrstat)
    IF (ierrstat /= 0) THEN
      CALL model_abort (my_cart_id, 8053, yerrmsg, yroutine)
    ENDIF
  ENDIF
#endif

  ! Write the headline in YUCHKDAT for this file
  IF ( (outblock%lcheck) .AND. (my_cart_id == 0) ) THEN
    OPEN(nuchkdat, FILE=yuchkdat, FORM=  'FORMATTED', STATUS='UNKNOWN',  &
                   POSITION='APPEND', IOSTAT=niostat)
    IF(niostat /= 0) THEN
      yerrmsg = ' ERROR    *** Error while opening file YUCHKDAT *** '
      ierrstat = 2005
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    ENDIF

    WRITE (nuchkdat,'(A,I7)')                                               &
                        'Check the data in output file for step: ', ntstep
    WRITE (nuchkdat,'(A,A)')                                                &
           '    File:   ',yname
    WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                                     &
           '    ie_tot =',ie_tot,'   je_tot =',je_tot,'   ke_tot =',ke_tot
    WRITE (nuchkdat,'(A)') '    '
    WRITE (nuchkdat,'(A,A)')                                                &
          '     var       ee    lev         min      ',                     &
          'imin   jmin          max      imax   jmax         mean  '
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Look for output variables in LM variable table
!------------------------------------------------------------------------------

  nzrecords = 0

  ! loop over all variables that should be written and loop over all variables
  ! in the LM variable table until equal elements are found
  write_loop: DO n = 1, numlist
    ! indices of field in variable table
    i1 = ilist(1,n)
    i2 = ilist(2,n)
    i3 = ilist(3,n)

    izlen = LEN_TRIM(ylist(n))
    IF (ylist(n)(1:izlen) == 'TKE') THEN
      itl = ntke
    ELSE
      itl = itimelevel
    ENDIF

    IF ( izdebug >= 5 ) THEN
      PRINT *, ' src_output: processing ', ylist(n)(1:izlen)
    ENDIF

    SELECT CASE (var(i1,i2,i3)%rank)
      ! pack data depending on the rank

    CASE(4)
      ! vertical interpolation, if necessary
      IF     (yextension == 'p') THEN
        CALL p_int (outblock, i1,i2,i3, n, itl, zvarlev(:,:,1:outblock%kepin))
        kbot = 1
        ktop = outblock%kepin
        slev(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
      ELSEIF (yextension == 'z') THEN
        CALL z_int (outblock, i1,i2,i3, n, itl, zvarlev(:,:,1:outblock%kezin))
        kbot = 1
        ktop = outblock%kezin
        slev(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)
      ELSE
        ! calculate additional non-global 4-d fields if required
        IF (ylist(n)(1:izlen) == 'TKE' .AND. (.NOT. lzrestart) ) THEN
          kbot =    1
          ktop = ke+1
          SELECT CASE( itype_turb )
          CASE( 5:8 )
            zvarlev(:,:,kbot:ktop-1) = var(i1,i2,i3)%p4(:,:,kbot:ktop-1,itl)
            zvarlev(:,:,ktop) = &
              0.5_ireals * (var(i1,i2,i3)%p4(:,:,ktop,itl))**2
          CASE default
            zvarlev(:,:,kbot:ktop) = &
              0.5_ireals * (var(i1,i2,i3)%p4(:,:,kbot:ktop,itl))**2
          END SELECT
        ELSEIF ( (ylist(n)(1:izlen) == 'U') .AND.    &
                    (outblock%luvmasspoint) .AND. (.NOT. lzrestart) ) THEN
          ! determine first and last level
          kbot = LBOUND (var(i1,i2,i3)%p4,3)
          ktop = UBOUND (var(i1,i2,i3)%p4,3)
          DO k = kbot, ktop
            DO j = 1, je
              DO i = 2, ie
                zvarlev(i,j,k) = 0.5 * (var(i1,i2,i3)%p4(i-1,j,k,itl) + &
                                        var(i1,i2,i3)%p4(i  ,j,k,itl))
              ENDDO
              zvarlev(1,j,k) = zvarlev(2,j,k)
            ENDDO
          ENDDO
        ELSEIF ( (ylist(n)(1:izlen) == 'V') .AND.    &
                    (outblock%luvmasspoint) .AND. (.NOT. lzrestart) ) THEN
          ! determine first and last level
          kbot = LBOUND (var(i1,i2,i3)%p4,3)
          ktop = UBOUND (var(i1,i2,i3)%p4,3)
          DO k = kbot, ktop
            DO j = 2, je
              DO i = 1, ie
                zvarlev(i,j,k) = 0.5 * (var(i1,i2,i3)%p4(i,j-1,k,itl) + &
                                        var(i1,i2,i3)%p4(i,j  ,k,itl))
              ENDDO
            ENDDO
            zvarlev(:,1,k) = zvarlev(:,2,k)
          ENDDO
        ELSE
          ! determine first and last level
          IF ( (ylist(n)(1:izlen) == 'T_SO') .AND. (outblock%yform_write == 'ncdf') ) THEN
            kbot = 1    !LBOUND (var(i1,i2,i3)%p4,3)+1
          ELSE
            kbot = LBOUND (var(i1,i2,i3)%p4,3)
          ENDIF
          ktop = UBOUND (var(i1,i2,i3)%p4,3)
          zvarlev(:,:,kbot:ktop) = var(i1,i2,i3)%p4(:,:,kbot:ktop,itl)
        END IF
      ENDIF

      DO k = kbot, ktop
        nzrecords = nzrecords + 1
        CALL output_data (nuedat, nzrecords, i1,i2,i3, k, ktop,            &
             zvarlev(1:ie,1:je,k), zprocarray_grib, zprocarray_real,       &
             outblock, .FALSE., yextension, slev(k), lzrestart, ivar_id(n),&
             iorg_data)
      ENDDO

    CASE(3)
      IF     (yextension == 's') THEN
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
        IF     (ylist(n)(1:izlen) == 'SYNME7') THEN
          ! Look for entry in sat_compute
          nentry = -1
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:8)=='METEOSAT') .AND.    &
                 (sat_compute(isens)%nsat_id        == 7        ) .AND.    &
                 (sat_compute(isens)%ysensor        =='MVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          kbot = 1
          IF (luse_rttov .AND. ALLOCATED(synme7) .AND. (nentry > 0)) THEN
             ktop = UBOUND(synme7,3)
             zvarlev(:,:,kbot:ktop) = synme7(:,:,kbot:ktop)
             slev(:) = REAL(nentry)
          ELSE
             ktop = 0
          ENDIF
        ELSEIF (ylist(n)(1:izlen) == 'SYNMSG') THEN
          ! Look for entry in sat_compute
          nentry = -1
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          kbot = 1
          IF (luse_rttov .AND. ALLOCATED(synmsg) .AND. (nentry > 0)) THEN
             ktop = UBOUND(synmsg,3)
             zvarlev(:,:,kbot:ktop) = synmsg(:,:,kbot:ktop)
             slev(:) = REAL(nentry)
          ELSE
             ktop = 0
          ENDIF
        ELSEIF (ylist(n)(1:izlen) == 'MSG_TB') THEN
          kbot = 1
          ktop = nmsgchan
          zvarlev(:,:,kbot:ktop) = synmsg(:,:,1:29:4)
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          slev(:) = REAL(nentry)
        ELSEIF (ylist(n)(1:izlen) == 'MSG_TBC') THEN
          kbot = 1
          ktop = nmsgchan
          zvarlev(:,:,kbot:ktop) = synmsg(:,:,2:30:4)
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          slev(:) = REAL(nentry)
        ELSEIF (ylist(n)(1:izlen) == 'MSG_RAD') THEN
          kbot = 1
          ktop = nmsgchan
          zvarlev(:,:,kbot:ktop) = synmsg(:,:,3:31:4)
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          slev(:) = REAL(nentry)
        ELSEIF (ylist(n)(1:izlen) == 'MSG_RADC') THEN
          kbot = 1
          ktop = nmsgchan
          zvarlev(:,:,kbot:ktop) =  synmsg(:,:,4:32:4)
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          slev(:) = REAL(nentry)
        ENDIF
#endif
      ELSEIF (yextension == 'p') THEN
        ! vertical interpolation on p-levels, if necessary
        CALL p_int (outblock, i1,i2,i3, n, itl, zvarlev(:,:,1:outblock%kepin))
        kbot = 1
        ktop = outblock%kepin
        slev(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
      ELSEIF (yextension == 'z') THEN
        ! vertical interpolation on z-levels, if necessary
        CALL z_int (outblock, i1,i2,i3, n, itl, zvarlev(:,:,1:outblock%kezin))
        kbot = 1
        ktop = outblock%kezin
        slev(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)
      ELSE
        !Calculate additional non-global 3-d fields if requird
        IF ( ylist(n)(1:izlen) == 'P' ) THEN
          kbot = 1
          ktop = ke
          zvarlev(:,:,kbot:ktop) = p0(:,:,kbot:ktop) + pp(:,:,kbot:ktop,itl)
        ELSEIF (ylist(n)(1:izlen) == 'OMEGA') THEN
          kbot = 1
          ktop = ke
          CALL calomega (zvarlev(:,:,kbot:ktop), pp(:,:,:,nnew),          &
                         pp(:,:,:,itl ), pptens(:,:,:), w(:,:,:,itl),     &
                         rho0(:,:,:), ie, je, ke, dt, g )
        ELSEIF (ylist(n)(1:izlen) == 'CLC') THEN
          kbot = 1
          ktop = ke
          zvarlev(:,:,kbot:ktop) = clc_sgs(:,:,kbot:ktop) &
             +  clc_con(:,:,kbot:ktop)*(1.0_ireals -  clc_sgs(:,:,kbot:ktop))
        ELSEIF (ylist(n)(1:izlen) == 'Q_SEDIM') THEN
          kbot = 1
          ktop = ke
          IF (lprog_qi) THEN
            zvarlev(:,:,kbot:ktop) = qrs(:,:,:) - qi(:,:,:,nnew)
          ELSE
            zvarlev(:,:,kbot:ktop) = qrs(:,:,:)
          ENDIF
        ELSEIF (ylist(n)(1:izlen) == 'RELHUM') THEN
          kbot = 1
          ktop = ke
          CALL calrelhum(zvarlev(:,:,kbot:ktop), t(:,:,:,itl), pp(:,:,:,itl),&
                         p0(:,:,:), qv(:,:,:,itl), ie, je, ke,               &
                         b1, b2w, b3, b4w, rdv, o_m_rdv )
        ELSEIF (ylist(n)(1:izlen) == 'BRN') THEN
          kbot = 1
          ktop = ke
          CALL calc_bulk_richardson(zvarlev(:,:,kbot:ktop),t(:,:,:,itl),   &
                qv(:,:,:,itl), u(:,:,:,itl), v(:,:,:,itl),                 &
                p0(:,:,:)+pp(:,:,:,itl), hsurf(:,:), ps(:,:,itl),          &
                t_2m(:,:), qv_2m(:,:),hhl(:,:,:),                          &
                ie, je, ke, cp_d, r_d,  rvd_m_o, g)
          zbrn(:,:,:) = zvarlev(:,:,kbot:ktop)
          l_brn = .TRUE.
        ELSEIF (ylist(n)(1:izlen) == 'FI_ANAI') THEN
          kbot = 1
          ktop = ke
          zvarlev(:,:,kbot:ktop) = p_anai(:,:,kbot:ktop) / rho(:,:,kbot:ktop)
        ELSEIF (ylist(n)(1:izlen) == 'DBZ') THEN
          kbot = 1
          ktop = ke
          IF (itype_gscp == 3) THEN
            CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,     &
              klv850, my_cart_id, itype_gscp, t(:,:,:,itl), qc(:,:,:,itl)*rho, qr(:,:,:,itl)*rho,   &
              qi(:,:,:,itl)*rho, qs(:,:,:,itl)*rho, z_radar = zvarlev(:,:,kbot:ktop))
          ELSEIF (itype_gscp == 4) THEN
            CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,     &
              klv850, my_cart_id, itype_gscp, t(:,:,:,itl), qc(:,:,:,itl)*rho, qr(:,:,:,itl)*rho,   &
              qi(:,:,:,itl)*rho, qs(:,:,:,itl)*rho, q_grau  = qg(:,:,:,itl)*rho,         &
              z_radar = zvarlev(:,:,kbot:ktop) )
          ENDIF
        ELSEIF ( ylist(n)(1:izlen) == 'VORTIC_U' ) THEN
          kbot = 1
          ktop = ke
          CALL curl (ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,    &
                     sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,               &
                     u(:,:,:,itl ), v(:,:,:,itl ), w(:,:,:,itl ),            &
                     .TRUE., zvarlev(:,:,1:ke), zhelp2(:,:,:), zhelp3(:,:,:))
        ELSEIF ( ylist(n)(1:izlen) == 'VORTIC_V' ) THEN
          kbot = 1
          ktop = ke
          CALL curl (ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,    &
                     sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,               &
                     u(:,:,:,itl ), v(:,:,:,itl ), w(:,:,:,itl ),            &
                     .TRUE., zhelp1(:,:,:), zvarlev(:,:,1:ke), zhelp3(:,:,:))
        ELSEIF ( ylist(n)(1:izlen) == 'VORTIC_W' ) THEN
          kbot = 1
          ktop = ke
          CALL curl (ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,    &
                     sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,               &
                     u(:,:,:,itl ), v(:,:,:,itl ), w(:,:,:,itl ),            &
                     .TRUE., zhelp1(:,:,:), zhelp2(:,:,:), zvarlev(:,:,1:ke))
        ELSEIF ( ylist(n)(1:izlen) == 'POT_VORTIC' ) THEN
          kbot = 1
          ktop = ke
          CALL calc_Theta_Tppp( t(:,:,:,itl ), pp(:,:,:,itl ), p0,           &
                                ie, je, ke, r_d, cp_d, zhelp4)
          CALL curl (ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,    &
                     sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,               &
                     u(:,:,:,itl ), v(:,:,:,itl ), w(:,:,:,itl ),            &
                     .FALSE., zhelp1, zhelp2, zhelp3)

          ! coriolis parameter with cosine only available for deep atmosphere
          IF ( .NOT. lcori_deep ) THEN
            zhelp2d(:,:) = 0.0_ireals
          ELSE
            zhelp2d(:,:) = fccos(:,:)
          ENDIF

          CALL potential_vorticity_rho( ie, je, ke, eddlon, eddlat, r_earth, &
                     fc, zhelp2d, sqrtg_r_s, dzeta_dlam, dzeta_dphi,         &
                     zhelp1, zhelp2, zhelp3, lmetr, zhelp4,                  &
                     u, v, w, zvarlev(:,:,1:ke))
          zvarlev(:,:,1:ke) = zvarlev(:,:,1:ke) / rho(:,:,1:ke)
        ELSE
          kbot = LBOUND(var(i1,i2,i3)%p3,3)
          ktop = UBOUND(var(i1,i2,i3)%p3,3)
          ! decision, if mlf or slf
          IF (ktop <= 3) THEN
            ! These are just 2D arrays with a timelevel
            kbot = 1
            ktop = 1
            zvarlev(:,:,1) = var(i1,i2,i3)%p3(1:ie,1:je,itl)
          ELSE
            zvarlev(:,:,kbot:ktop) = var(i1,i2,i3)%p3(:,:,kbot:ktop)
          ENDIF
        ENDIF
      ENDIF

      ! Distribute the multidimensional fields to the PEs
      DO k=kbot,ktop
        nzrecords = nzrecords + 1
        CALL output_data (nuedat, nzrecords, i1,i2,i3, k, ktop,            &
             zvarlev(1:ie,1:je,k), zprocarray_grib, zprocarray_real,       &
             outblock, .FALSE., yextension, slev(k), lzrestart, ivar_id(n),&
             iorg_data)
      ENDDO

    CASE(2)
      ! Calculate additional non-global output fields if required
      IF ( ylist(n)(1:izlen) == 'PMSL' ) THEN
        CALL calpmsl( zvarlev(:,:,1), ps(:,:,itl), t(:,:,ke,itl),       &
               rho0(:,:,ke), dp0(:,:,ke), hsurf, ie, je, g, r_d )
      ELSEIF ( ylist(n)(1:izlen) == 'TOT_PREC' ) THEN
        IF (itype_gscp == 4) THEN
          CALL calprsum( zvarlev(:,:,1), rain_gsp, snow_gsp + grau_gsp, &
                           rain_con, snow_con, ie, je )
        ELSE
          CALL calprsum( zvarlev(:,:,1), rain_gsp, snow_gsp, rain_con,  &
                         snow_con, ie, je )
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'HTOP_DC' ) THEN
        CALL caltopdc( zvarlev(:,:,1), t(:,:,:,itl), p0, pp(:,:,:,itl), &
                       qv(:,:,:,itl), hhl, hhlr, ie, je, ke,            &
                       b1, b2w, b3, b4w, rdv, o_m_rdv, rvd_m_o, cpdr, g )
      ELSEIF ( ylist(n)(1:izlen) == 'FIS' ) THEN
        zvarlev(1:ie,1:je,1) = hsurf(1:ie,1:je) * g
      ELSEIF ( ylist(n)(1:izlen) == 'HTOP_CON'         .OR.      &
               ylist(n)(1:izlen) == 'HTOP_SC' ) THEN
        DO j= 1, je
          DO i = 1, ie 
            zvarlev(i,j,1) = 0.0_ireals
            klev = NINT( top_con(i,j) )
            IF(klev > 0) THEN
              zvarlev(i,j,1) = 0.5*(hhl(i,j,klev)+hhl(i,j,klev+1))
            ENDIF
          ENDDO 
        ENDDO
      ELSEIF ( ylist(n)(1:izlen) == 'HBAS_CON'         .OR.      &
               ylist(n)(1:izlen) == 'HBAS_SC' ) THEN
        DO j= 1, je
          DO i = 1, ie
            zvarlev(i,j,1) = 0.0_ireals
            klev = NINT( bas_con(i,j) )
            IF(klev > 0) THEN
              zvarlev(i,j,1) = hhl(i,j,klev)
            ENDIF
          ENDDO 
        ENDDO
      ELSEIF ( ylist(n)(1:izlen) == 'HZEROCL' ) THEN
        CALL calhzero( zvarlev(:,:,1), t(:,:,:,itl), hhl, hhlr, ie, je, ke, t0_melt)
      ELSEIF ( ylist(n)(1:izlen) == 'SNOWLMT' ) THEN
        CALL calsnowlmt( zvarlev(:,:,1), t(:,:,:,itl), pp(:,:,:,itl),      &
              p0(:,:,:), qv(:,:,:,itl), hhl, hhlr, ie, je, ke, t0_melt, 1.3_ireals)
      ELSEIF ( ylist(n)(1:izlen) == 'CLCT_MOD' ) THEN
        CALL calclmod(zvarlev(:,:,1), clc_sgs, clc_con, p0hl, pi, ie, je, ke)
      ELSEIF ( ylist(n)(1:izlen) == 'CLDEPTH' ) THEN
        CALL calcldepth( zvarlev(:,:,1), clc_sgs, clc_con, dp0, ie, je, ke )
      ELSEIF ( ylist(n)(1:izlen) == 'TQV' ) THEN
        CALL caliq( zvarlev(:,:,1), rho, hhl, qv(:,:,:,itl ), ie, je, ke )
      ELSEIF ( ylist(n)(1:izlen) == 'TQC' ) THEN
        CALL caliq( zvarlev(:,:,1), rho, hhl, qc(:,:,:,itl ), ie, je, ke )
      ELSEIF ( ylist(n)(1:izlen) == 'TQR' ) THEN
        IF ( ALLOCATED(qr) ) THEN
          CALL caliq( zvarlev(:,:,1), rho, hhl, qr(:,:,:,itl ), ie, je, ke )
        ELSE
          zvarlev(:,:,1) = 0.0_ireals
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'TQS' ) THEN
        IF ( ALLOCATED(qs) ) THEN
          CALL caliq( zvarlev(:,:,1), rho, hhl, qs(:,:,:,itl ), ie, je, ke )
        ELSE
          zvarlev(:,:,1) = 0.0_ireals
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'TQG' ) THEN
        IF ( ALLOCATED(qg) ) THEN
          CALL caliq( zvarlev(:,:,1), rho, hhl, qg(:,:,:,itl ), ie, je, ke )
        ELSE
          zvarlev(:,:,1) = 0.0_ireals
        ENDIF
      ELSEIF ( (ylist(n)(1:izlen) == 'ZTD') .OR.                     &
               (ylist(n)(1:izlen) == 'ZWD') .OR.                     &
               (ylist(n)(1:izlen) == 'ZHD') ) THEN
        IF (.NOT. lzenith) THEN
          CALL calztd( zenith_t(:,:), zenith_w(:,:), zenith_h(:,:),        &
                       rho, hhl, qv(:,:,:,itl ), ps(:,:,itl ),             &
                       t(:,:,ke,itl ), hsurf, rlat, pi, ie, je, ke)
          lzenith = .TRUE.
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'TWATER' ) THEN
        zvarlev(:,:,2:ke+1) = qv (:,:,1:ke,itl ) + qc(:,:,1:ke,itl )
        IF (lprog_qi) THEN
          zvarlev(:,:,2:ke+1) = zvarlev(:,:,2:ke+1) + qi(:,:,1:ke,itl )
        ENDIF
        IF (ALLOCATED(qr)) THEN
          zvarlev(:,:,2:ke+1) = zvarlev(:,:,2:ke+1) + qr(:,:,1:ke,itl )
        ENDIF
        IF (ALLOCATED(qs)) THEN
          zvarlev(:,:,2:ke+1) = zvarlev(:,:,2:ke+1) + qs(:,:,1:ke,itl )
        ENDIF
        IF (ALLOCATED(qg)) THEN
          zvarlev(:,:,2:ke+1) = zvarlev(:,:,2:ke+1) + qg(:,:,1:ke,itl )
        ENDIF
        CALL caliq( zvarlev(:,:,1), rho, hhl, zvarlev(:,:,2:ke+1), ie,je,ke)
      ELSEIF ( ylist(n)(1:izlen) == 'TQI') THEN
        IF (lprog_qi) THEN
          CALL caliq( zvarlev(:,:,1), rho, hhl, qi(:,:,:,itl ), ie,je,ke)
        ELSE
          CYCLE write_loop
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'TQV_ANAI') THEN
        ! for analysis increments, using rho(nnow) is sufficiently accurate
        CALL caliq( zvarlev(:,:,1), rho, hhl, qv_anai(:,:,:), ie,je,ke)
      ELSEIF ( ylist(n)(1:izlen) == 'TQC_ANAI') THEN
        CALL caliq( zvarlev(:,:,1), rho, hhl, qc_anai(:,:,:), ie,je,ke)
      ELSEIF ( ylist(n)(1:izlen) == 'PMSL_ANAI') THEN
        CALL calpmsl( zvarlev(:,:,1), ps(:,:,itl), t(:,:,ke,itl),          &
               rho0(:,:,ke), dp0(:,:,ke), hsurf, ie, je, g, r_d )
        zvarlev(:,:,1) = p_anai(:,:,ke) * zvarlev(:,:,1)                   &
                                        /(p0(:,:,ke) + pp(:,:,ke,itl))
      ELSEIF ( (ylist(n)(1:izlen) == 'AUMFL_S') .AND.     &
               (.NOT. lzrestart) .AND. outblock%luvmasspoint) THEN
        ! Exchange aumfl_s first
        kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                &
          (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
          ie, je, kzdims, jstartpar, jendpar,                                &
          nbl_exchg, nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,    &
          20000+ntstep, .FALSE.   , ncomm_type, izerror, yerrmsg,            &
          var(i1,i2,i3)%p2(:,:))
        DO j = 1, je
          DO i = 2, ie
            zvarlev(i,j,1) = 0.5 * (var(i1,i2,i3)%p2(i-1,j) +      &
                                    var(i1,i2,i3)%p2(i  ,j))
          ENDDO
          zvarlev(1,j,1) = zvarlev(2,j,1)
        ENDDO
      ELSEIF ( (ylist(n)(1:izlen) == 'AVMFL_S') .AND.     &
               (.NOT. lzrestart) .AND. outblock%luvmasspoint) THEN
        ! Exchange avmfl_s first
        kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                &
          (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
          ie, je, kzdims, jstartpar, jendpar,                                &
          nbl_exchg, nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,    &
          20000+ntstep, .FALSE.   , ncomm_type, izerror, yerrmsg,            &
          var(i1,i2,i3)%p2(:,:))
        DO j = 2, je
          DO i = 1, ie
            zvarlev(i,j,1) = 0.5 * (var(i1,i2,i3)%p2(i,j-1) +      &
                                    var(i1,i2,i3)%p2(i,j  ))
          ENDDO
        ENDDO
        zvarlev(:,1,1) = zvarlev(:,2,1)
      ELSEIF (ylist(n)(1:izlen) == 'DBZ_850') THEN
        IF (itype_gscp == 3) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,     &
            klv850, my_cart_id, itype_gscp, t(:,:,:,itl), qc(:,:,:,itl)*rho, qr(:,:,:,itl)*rho,   &
            qi(:,:,:,itl)*rho, qs(:,:,:,itl)*rho, z_radar_850 = zvarlev(:,:,1) )
        ELSEIF (itype_gscp == 4) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,     &
            klv850, my_cart_id, itype_gscp, t(:,:,:,itl), qc(:,:,:,itl)*rho, qr(:,:,:,itl)*rho,   &
            qi(:,:,:,itl)*rho, qs(:,:,:,itl)*rho, q_grau  = qg(:,:,:,itl)*rho,         &
            z_radar_850 = zvarlev(:,:,1) )
        ENDIF
      ELSEIF (ylist(n)(1:izlen) == 'DBZ_CMAX') THEN
        IF (itype_gscp == 3) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,     &
            klv850, my_cart_id, itype_gscp, t(:,:,:,itl), qc(:,:,:,itl)*rho, qr(:,:,:,itl)*rho,   &
            qi(:,:,:,itl)*rho, qs(:,:,:,itl)*rho, z_radar_cmax = zvarlev(:,:,1) )
        ELSEIF (itype_gscp == 4) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,     &
            klv850, my_cart_id, itype_gscp, t(:,:,:,itl), qc(:,:,:,itl)*rho, qr(:,:,:,itl)*rho,   &
            qi(:,:,:,itl)*rho, qs(:,:,:,itl)*rho, q_grau  = qg(:,:,:,itl)*rho,         &
            z_radar_cmax = zvarlev(:,:,1) )
        ENDIF
      ELSEIF (ylist(n)(1:izlen) == 'SWISS00') THEN
        CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:,itl),                    &
             u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),            &
             p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w, b3, &
             b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,               &
             missing_value, izdebug, lwarn, ierrstat, yerrmsg,             &
             swiss00 = zvarlev(:,:,1))
      ELSEIF (ylist(n)(1:izlen) == 'SWISS12') THEN
        CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:,itl),                    &
             u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),            &
             p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w, b3, &
             b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,               &
             missing_value, izdebug, lwarn, ierrstat, yerrmsg,             &
             swiss12 = zvarlev(:,:,1))
      ELSEIF (ylist(n)(1:izlen) == 'SI') THEN
        CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:,itl),                    &
             u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),            &
             p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w, b3, &
             b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,               &
             missing_value, izdebug, lwarn, ierrstat, yerrmsg,             &
             si= zvarlev(:,:,1))
      ELSEIF (ylist(n)(1:izlen) == 'SLI') THEN
        CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:,itl),                    &
             u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),            &
             p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w, b3, &
             b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,               &
             missing_value, izdebug, lwarn, ierrstat, yerrmsg,             &
             sli= zvarlev(:,:,1))
      ELSEIF ( (ylist(n)(1:izlen) == 'CAPE_MU')  .OR. &
               (ylist(n)(1:izlen) == 'CIN_MU') ) THEN
        IF (.NOT. lconvind_mu) THEN
           CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:,itl),                 &
                u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),         &
                p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w,  &
                b3, b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,        &
                missing_value, izdebug, lwarn, ierrstat, yerrmsg,          &
                cape_mu=zcape_mu, cin_mu=zcin_mu)
           lconvind_mu = .TRUE.
        ENDIF
      ELSEIF ( (ylist(n)(1:izlen) == 'CAPE_ML')  .OR. &
               (ylist(n)(1:izlen) == 'CIN_ML')   .OR. &
               (ylist(n)(1:izlen) == 'LCL_ML')   .OR. &
               (ylist(n)(1:izlen) == 'LFC_ML')   .OR. &
               (ylist(n)(1:izlen) == 'CAPE_3KM') ) THEN
        IF (.NOT. lconvind_ml) THEN
           CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:,itl),                 &
                u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),         &
                p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w,  &
                b3, b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,        &
                missing_value, izdebug, lwarn, ierrstat, yerrmsg,          &
                cape_ml=zcape_ml, cin_ml=zcin_ml, lcl_ml=zlcl_ml,          &
                lfc_ml=zlfc_ml, cape_3km=zcape_3km)
           lconvind_ml = .TRUE.
        ENDIF
      ELSEIF (ylist(n)(1:izlen) == 'HPBL') THEN
        IF (.NOT. l_brn) THEN
          kbot = 1
          ktop = ke
          CALL calc_bulk_richardson(zvarlev(:,:,kbot:ktop),t(:,:,:,itl),  &
               qv(:,:,:,itl), u(:,:,:,itl), v(:,:,:,itl),                 &
               p0(:,:,:)+pp(:,:,:,itl), hsurf(:,:), ps(:,:,itl),          &
               t_2m(:,:), qv_2m(:,:),hhl(:,:,:),                          &
               ie, je, ke, cp_d, r_d,  rvd_m_o, g)
          l_brn = .TRUE.
          zbrn(:,:,:) = zvarlev(:,:,kbot:ktop)
        ENDIF
        CALL  calc_pbl_brn(t(:,:,:,itl), qv(:,:,:,itl),                    &
                p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:),hsurf(:,:),            &
                zbrn(:,:,:), ie, je, ke, cp_d, r_d, rvd_m_o, missing_value,&
                zvarlev(:,:,1))
        zvarlev(:,:,0)=zvarlev(:,:,1)
      ELSEIF (ylist(n)(1:izlen) == 'CEILING') THEN
        CALL calc_ceiling( zvarlev(:,:,1), clc_sgs, hhl, ie, je, ke )
      ELSEIF ( (ylist(n)(1:izlen) == 'SDI_1') .OR.          &
               (ylist(n)(1:izlen) == 'SDI_2') )  THEN
        IF (.NOT. l_calc_sdi_already_comp ) THEN
          CALL calc_sdi( zvarlev(:,:,1), zvarlev(:,:,2) )
          l_calc_sdi_already_comp = .TRUE.
        END IF
      ELSE
        IF (ASSOCIATED(var(i1,i2,i3)%p2)) THEN
          zvarlev(1:ie,1:je,1) = var(i1,i2,i3)%p2(1:ie,1:je)
        ELSE
          PRINT *, ' *** ERROR: Trying to output unassociated variable: ', &
                     ylist(n)(1:izlen)
        ENDIF
      ENDIF


      nzrecords   = nzrecords+1
      IF     ( ylist(n)(1:izlen) == 'ZTD' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zenith_t(:,:), zprocarray_grib, zprocarray_real, &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)
      ELSEIF ( ylist(n)(1:izlen) == 'ZWD' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zenith_w(:,:), zprocarray_grib, zprocarray_real, &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)
      ELSEIF ( ylist(n)(1:izlen) == 'ZHD' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zenith_h(:,:), zprocarray_grib, zprocarray_real, &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)
      ELSEIF ( ylist(n)(1:izlen) == 'CAPE_MU' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zcape_mu(:,:), zprocarray_grib, zprocarray_real, &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)

      ELSEIF ( ylist(n)(1:izlen) == 'CIN_MU' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zcin_mu(:,:), zprocarray_grib, zprocarray_real,  &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)

      ELSEIF ( ylist(n)(1:izlen) == 'CAPE_ML' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zcape_ml(:,:), zprocarray_grib, zprocarray_real, &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)

      ELSEIF ( ylist(n)(1:izlen) == 'CIN_ML' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zcin_ml(:,:), zprocarray_grib, zprocarray_real,  &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)
      ELSEIF ( ylist(n)(1:izlen) == 'LCL_ML' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zlcl_ml(:,:), zprocarray_grib, zprocarray_real,  &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)

      ELSEIF ( ylist(n)(1:izlen) == 'LFC_ML' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zlfc_ml(:,:), zprocarray_grib, zprocarray_real,  &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)

      ELSEIF ( ylist(n)(1:izlen) == 'CAPE_3KM' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zcape_3km(:,:), zprocarray_grib, zprocarray_real,  &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)
      ELSEIF ( ylist(n)(1:izlen) == 'SDI_1' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zvarlev(:,:,1), zprocarray_grib, zprocarray_real,&
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)

      ELSEIF ( ylist(n)(1:izlen) == 'SDI_2' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zvarlev(:,:,2), zprocarray_grib, zprocarray_real,&
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)

      ELSE
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1,                &
                          zvarlev(:,:,1), zprocarray_grib, zprocarray_real,&
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), iorg_data)
      ENDIF

    END SELECT
  ENDDO write_loop

!CPS#ifdef COUP_OAS_COS
!CPS   reset_loop: DO n = 1, numlist
!CPS      ! indices of field in variable table
!CPS      i1 = ilist(1,n)
!CPS      i2 = ilist(2,n)
!CPS      i3 = ilist(3,n)

!CPS      SELECT CASE (var(i1,i2,i3)%rank)
        ! pack data depending on the rank

!CPS      CASE(2)
        ! reset summation variables
!CPS        IF ( (var(i1,i2,i3)%ntri == 3)) THEN
!CPS          IF (ASSOCIATED(var(i1,i2,i3)%p2)) THEN
!CPS            var(i1,i2,i3)%p2(:,:) = 0.0_ireals
!CPS          ENDIF
!CPS        ENDIF

!CPS      END SELECT
!CPS    ENDDO reset_loop
!CPS#else
  IF (lbdclim .AND. (.NOT. lzrestart)) THEN
    ! Uwe Boehm 09.08.2007
    ! reset all summation variables in an additional loop and NOT in the write_loop
    ! to avoid reset of total precipitation components before adding them up
    reset_loop: DO n = 1, numlist
      ! indices of field in variable table
      i1 = ilist(1,n)
      i2 = ilist(2,n)
      i3 = ilist(3,n)

      SELECT CASE (var(i1,i2,i3)%rank)
        ! pack data depending on the rank

      CASE(2)
        ! reset summation variables
        IF ( (var(i1,i2,i3)%ntri == 3) .OR. (var(i1,i2,i3)%ntri == 4) ) THEN
          IF (ASSOCIATED(var(i1,i2,i3)%p2)) THEN
            var(i1,i2,i3)%p2(:,:) = 0.0_ireals
          ENDIF
        ENDIF

      END SELECT
    ENDDO reset_loop
  ENDIF
!CPS#endif
!------------------------------------------------------------------------------
! Section 4: Flush output buffers and close grib file
!------------------------------------------------------------------------------

  CALL output_data (nuedat, -1, -1,-1,-1, -1, -1,           zvarlev(:,:,1),&
                    zprocarray_grib, zprocarray_real, outblock, .TRUE.,    &
                    yextension, 0.0_ireals, lzrestart, -1, iorg_data)
  CALL close_file (nuedat, TRIM(outblock%yform_write), icomm_cart, my_cart_id, &
                   num_compute, lasync_io, idbg_level, yerrmsg, ierrstat)

  IF (ierrstat /= 0) THEN
     CALL model_abort (my_cart_id, 2034, yerrmsg, yroutine)
  ENDIF

  IF (outblock%yform_write == 'bina') THEN
    ! release the unit-number again
    CALL release_unit (nuedat)
  ENDIF

  ! Write a blank line to YUCHKDAT
  IF ( (outblock%lcheck) .AND. (my_cart_id == 0) ) THEN
    WRITE (nuchkdat,'(A)') '   '
    WRITE (nuchkdat,'(A)') '   '
  ENDIF

  ! close file nuchkdat
  IF ( (outblock%lcheck) .AND. (my_cart_id == 0) ) THEN
    CLOSE (nuchkdat, STATUS='KEEP')
  ENDIF

  ! Deallocate arrays for IO
  DEALLOCATE (iblock, ibmap, ds_real, ds_grib, dsup)

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE organize_output

!==============================================================================
!+ distributes records to PEs and packs them into grib format
!------------------------------------------------------------------------------

SUBROUTINE output_data (nuedat, irec, i1,i2,i3, k, klevels,                 &
                        array2d_real, procarray_grib, procarray_real,       &
                        outblock, lflush, yextension, slev, lrestart,       &
                        incdf_var_id, my_orgdata)

!------------------------------------------------------------------------------
!
! Description:
!  output_data distributes records to the PEs for packing these into 
!  grib format. First, the records are only gathered from all PEs and
!  every PE stores one record into the variable procarray_xxx. Only if every
!  PE has got a record, the data are packed and written to disk (in the
!  routine "write_xxxx"). If some PEs have got no record because no more
!  records are left, the output buffers (variable procarray_xx) are "flushed".
!
! Method:
!  output_data is called for every record that is processed. The PE that
!  gets a special record, saves the characteristics of this record for the
!  output step later on.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  nuedat,  & ! descriptor of the grib file
  irec,    & ! number of record to be processed
  i1,i2,i3,& ! location of the variable in the LM variable table
  k,       & ! number of the actual level
  klevels    ! number of levels this variable has

TYPE(pp_nl),              INTENT(IN)     ::    &
  outblock        ! pointer to the namelist group

LOGICAL                 , INTENT (IN)    ::    &
  lflush,  & ! for flushing the output buffers
  lrestart   ! whether restart-files are written or not

CHARACTER (LEN= 1)      , INTENT (IN)    ::    &
  yextension ! to check which output list is processed

REAL (KIND=ireals)      , INTENT (IN)    ::    &
  slev       ! level for vertical interpolated fields
             ! has to be present for p- and z-levels

! Array arguments with intent(inout):
REAL (KIND=ireals)      , INTENT (INOUT) ::    &
  array2d_real (ie,je)  , & ! values of the variable to be processed
  procarray_real(ie_max,je_max,num_compute)

REAL (KIND=irealgrib)   , INTENT (INOUT) ::    &
  procarray_grib(ie_max,je_max,num_compute)

INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  incdf_var_id              ! NetCDF-ID of each variable in the output list
                            ! only PE 0 has a reasonable value here

INTEGER (KIND=iintegers), INTENT(INOUT)  ::    &
  my_orgdata(3,0:num_compute-1) ! necessary only for PE 0 to save information

!------------------------------------------------------------------------------
!
! Local parameters:

!
! Local scalars:
INTEGER (KIND=intgribf)         :: iz_ps=1
INTEGER (KIND=intgribf)         :: ierrf, tri, my_k_f, my_iee
INTEGER (KIND=iintegers)        :: npe, iz_lfd, izerror, izlen, i,j,         &
                                   irecord_len, implcode, nlastout,          &
                                   ij_out, ij_ful

REAL    (KIND=ireals)           :: zstartlon_tot, zendlon_tot, zstartlat_tot,&
                                   zendlat_tot, zavgfactor
REAL    (KIND=ireals)           :: zbias, zgribfactor, array2dreal(ie_max, je_max)
REAL    (KIND=irealgrib)        :: array2d_grib(ie_max,je_max),              &
                                   ds_out(ie_tot*je_tot), undefsub

CHARACTER (LEN=25)              :: yroutine
CHARACTER (LEN=80)              :: yerrmsg
CHARACTER (LEN=10)              :: yzname 

INTEGER (KIND=iintegers), SAVE  :: my_i1, my_i2, my_i3, my_k, my_irec
LOGICAL                 , SAVE  :: loutput
REAL    (KIND=ireals)   , SAVE  :: rmy_slev

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

  yroutine = 'output_data'
  yzname   = '          '
  izerror  = 0_iintegers
  iz_lfd   = INT (lfd, iintegers)
  IF (outblock%yform_write /= 'ncdf') THEN
    undef     = REAL(undefgrib, ireals)
    undefsub  = undefgrib
  ELSE
    undef     = REAL(undefncdf, ireals)
    undefsub  = undefncdf
  ENDIF

  IF (irec == 1) THEN
    my_orgdata(:,:) = 0_iintegers
    loutput         = .FALSE.
  ENDIF

  IF (.NOT. lrestart) THEN
    IF (outblock%nextstep == 1) THEN
      nlastout = 0
      ! Adaptation, if summation and meanvalues are done between output steps
      ! and the first output step is not 0
      IF (outblock%ngrib( outblock%nextstep) > 0 ) THEN
          nlastout = outblock%ngrib( outblock%nextstep ) -                             &
                    (outblock%ngrib( outblock%nextstep+1 )-outblock%ngrib( outblock%nextstep ))
        !US should it be:  nlastout = outblock%ngrib( outblock%nextstep) ???????
      ENDIF
    ELSE
      nlastout = outblock%ngrib( outblock%nextstep-1 )
    ENDIF
  ELSE
    nlastout = outblock%nextstep - nhour_restart(3) * NINT (3600.0_ireals / dt)
  ENDIF

!------------------------------------------------------------------------------
!  Section 2: If this is not a call to only flush the output data,
!             gather the field on the next free PE
!------------------------------------------------------------------------------

  IF (.NOT.lflush) THEN

    ! Get the number of the PE to deal with that slice
    npe     = MOD(irec-1,num_compute)

    ! Save necessary values for later processing
    IF (my_cart_id == npe ) THEN
      my_irec   = irec
      my_i1     = i1
      my_i2     = i2
      my_i3     = i3
      my_k      = k
      rmy_slev  = slev
      loutput   = .TRUE.
    ENDIF

    ! Processor 0 has to save some organizational data
    IF (my_cart_id == 0) THEN
      my_orgdata(1,npe) = k
      my_orgdata(2,npe) = klevels
      my_orgdata(3,npe) = incdf_var_id
    ENDIF

    ! scale the fields with time range indicator 3 
    IF ((.NOT. lrestart) .AND. (var(i1,i2,i3)%ntri == 3)) THEN
      IF (ntstep == 0) THEN
        zavgfactor = 1.0_ireals
      ELSE
        IF (lbdclim) THEN
          ! averaging is done between output steps
          zavgfactor = 1.0_ireals / REAL  (ntstep - nlastout, ireals)
        ELSE
          ! averaging is done between beginning of forecast and
          ! actual output step
          IF (l2tls) THEN
            zavgfactor = 1.0_ireals / REAL  (ntstep+1, ireals)
          ELSE
            zavgfactor = 1.0_ireals / REAL  (ntstep, ireals)
          ENDIF
!CPS#ifdef COUP_OAS_COS
!CPS          ! averaging is done between output steps
!CPS          zavgfactor = 1.0_ireals / REAL  (ntstep - nlastout, ireals)
!CPS#endif
        ENDIF
      ENDIF
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          array2d_real(i,j) = array2d_real(i,j) * zavgfactor
        ENDDO
      ENDDO
    ENDIF

    ! Scale field with factor and bias (only for grib)
    ! transform it to single precision (for grib and Netcdf)
    zbias       = var(i1,i2,i3)%bias
    zgribfactor = var(i1,i2,i3)%factor
    IF (outblock%yform_write == 'grb1') THEN
      ! Check for undefined values in case of NetCDF-Input and Grib output
      IF (yform_read == 'ncdf') THEN
        ! this is a very pragmatic solution, which might not be satisfactory
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            IF (array2d_real(i,j) == REAL(undefncdf, ireals)) THEN
              array2d_real(i,j) = 0.0_ireals
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          ! Do an additional clipping, if values are too small
          IF (ABS(array2d_real(i,j)) < 1.0E-15_ireals) THEN
            array2d_grib(i,j) =                                             &
              REAL (((0.0_ireals        + zbias) * zgribfactor), irealgrib)
          ELSE
            array2d_grib(i,j) =                                             &
              REAL (((array2d_real(i,j) + zbias) * zgribfactor), irealgrib)
          ENDIF
        ENDDO
      ENDDO
    ELSEIF (outblock%yform_write == 'ncdf') THEN
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          ! No scaling for NetCDF
          array2d_grib(i,j) = REAL (array2d_real(i,j), irealgrib)
        ENDDO
      ENDDO

      ! take care of special fields
      IF (var(i1,i2,i3)%lsm == 'l') THEN
        WHERE (.NOT. llandmask(1:ie,1:je)) array2d_grib(1:ie,1:je) = undefncdf
      ENDIF
      IF (var(i1,i2,i3)%lsm == 's') THEN
        WHERE (llandmask(1:ie,1:je)) array2d_grib(1:ie,1:je) = undefncdf
      ENDIF
      IF (TRIM(var(i1,i2,i3)%name) == 'SNOWLMT' .OR. &
          TRIM(var(i1,i2,i3)%name) == 'HZEROCL') THEN
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            IF (array2d_grib(i,j) == -999.0_irealgrib) THEN
              array2d_grib(i,j) = undefncdf
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF (TRIM(var(i1,i2,i3)%name) == 'HBAS_CON' .OR. &
          TRIM(var(i1,i2,i3)%name) == 'HTOP_CON') THEN
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            IF (array2d_grib(i,j) == 0.0_irealgrib) THEN
              array2d_grib(i,j) = undefncdf
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF (TRIM(var(i1,i2,i3)%name) == 'T_SNOW') THEN
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            IF (w_snow(i,j,nnow) == 0.0_ireals) THEN
              array2d_grib(i,j) = undefncdf
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF (TRIM(var(i1,i2,i3)%name) == 'Z0') THEN
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            array2d_grib(i,j) = array2d_grib(i,j) * 0.10197
          ENDDO
        ENDDO
      ENDIF

    ELSEIF (outblock%yform_write == 'bina') THEN
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          ! no scaling and converting in case of restart files,
          ! because this is not reproducible!
          array2dreal(i,j) = array2d_real(i,j)
        ENDDO
      ENDDO
    ENDIF

    IF ( ldebug_io .AND. (idbg_level >= 10) .AND. (my_cart_id == npe)) THEN
      PRINT *, ' src_output: gathering  ', TRIM(var(my_i1, my_i2, my_i3)%name), my_k, &
               ' to PE ', npe
    ENDIF

    IF (ltime) CALL get_timings (i_computations_O, ntstep, dt, izerror)

    ! Gather the data
    IF (.NOT. lrestart) THEN
      IF (num_compute > 1) THEN
        CALL gather_values (array2d_grib, procarray_grib, ie_max, je_max, &
               num_compute, imp_grib, npe, icomm_cart, yerrmsg, izerror)
      ELSE
        procarray_grib(:,:,1) = array2d_grib(:,:)
      ENDIF
    ELSE
      IF (num_compute > 1) THEN
        CALL gather_values (array2dreal, procarray_real, ie_max, je_max,     &
               num_compute, imp_reals, npe, icomm_cart, yerrmsg, izerror)
      ELSE
        procarray_real(:,:,1) = array2dreal(:,:)
      ENDIF
    ENDIF

    IF (ltime) CALL get_timings (i_gather_data, ntstep, dt, izerror)

  ENDIF

!------------------------------------------------------------------------------
!  Section 3: If lflush is .TRUE. or all PEs have gotten data, do the output
!------------------------------------------------------------------------------

IF ( lflush .OR. MOD(irec-1,num_compute) == num_compute-1) THEN

   !---------------------------------------------------------------------------
   !  Section 3.1: All PEs that have gotten data must combine the subarrays
   !               and convert data to GRIB
   !---------------------------------------------------------------------------

   IF ( loutput ) THEN

     IF (outblock%yform_write /= 'ncdf') THEN
 
       ! complete gds

       ! number of gridpoints (depend on namelist group)
       igds_out( 5) = outblock%ie_out_tot
       igds_out( 6) = outblock%je_out_tot

       ! calculation of the left bottom and upper right corner
       ! these depend also on the namelist group
       zstartlon_tot = startlon_tot + REAL (outblock%i_out_start-1,ireals)*dlon
       zstartlat_tot = startlat_tot + REAL (outblock%j_out_start-1,ireals)*dlat

       ! Careful: In case of p- or z-levels, u and v are interpolated to 
       !          to the mass grid point. The gds has not to be "staggered"
       !          in this case
       yzname = var(my_i1,my_i2,my_i3)%name    ! to shorten the following lines
       izlen  = LEN_TRIM(yzname)
       IF( ( (yzname(1:izlen) == 'U')            .OR.      &
             (yzname(1:izlen) == 'AUMFL_S') )    .AND.     &
           (.NOT. outblock%luvmasspoint)         .AND.     &
           (yextension /= 'p') .AND. (yextension /= 'z') ) THEN
         zstartlon_tot = zstartlon_tot + 0.5_ireals * dlon
       ENDIF
       zendlon_tot = zstartlon_tot + (outblock%ie_out_tot-1)*dlon

       IF( ( (yzname(1:izlen) == 'V')            .OR.      &
             (yzname(1:izlen) == 'AVMFL_S') )    .AND.     &
           (.NOT. outblock%luvmasspoint)         .AND.     &
           (yextension /= 'p') .AND. (yextension /= 'z') ) THEN
         zstartlat_tot = zstartlat_tot + 0.5_ireals * dlat
       ENDIF
       zendlat_tot = zstartlat_tot + (outblock%je_out_tot-1)*dlat

       ! All longitude values have to be limited to the range (-180.0,+180.0)
       IF (zstartlon_tot > 180.0_ireals) THEN
         zstartlon_tot = zstartlon_tot - 360.0_ireals
       ENDIF
       IF (zendlon_tot > 180.0_ireals) THEN
         zendlon_tot = zendlon_tot - 360.0_ireals
       ENDIF

       igds_out( 7) = NINT(zstartlat_tot * 1000.0_ireals)
       igds_out( 8) = NINT(zstartlon_tot * 1000.0_ireals)
       igds_out(10) = NINT(zendlat_tot   * 1000.0_ireals)
       igds_out(11) = NINT(zendlon_tot   * 1000.0_ireals)

       ! create pds
       call makepds (my_i1,my_i2,my_i3, my_k, tri, outblock%nunit_of_time, &
                     nlastout, outblock%nprocess_ini_out,                  &
                     outblock%nprocess_bd_out, outblock%lanalysis,         &
                     yextension, rmy_slev, lrestart)
     ENDIF

     IF (.NOT. lrestart) THEN
       ! combine the subarrays in the correct order
       CALL combine_subarrays (procarray_grib, ds_grib)

       IF (outblock%l_fi_ps_smooth) THEN
         ! Special smoothing in mountainous terrain and
         ! apply a digital smoother for selected fields
         IF (     (var(my_i1,my_i2,my_i3)%name == 'PMSL     ')             &
             .OR. (var(my_i1,my_i2,my_i3)%name == 'PMSL_ANAI')) THEN
           CALL smooth_pmsl (ds_grib, hsurf_tot, ie_tot, je_tot )
           CALL smoother(ds_grib, ie_tot, je_tot, 4,20 )
         ENDIF
         IF(var(my_i1,my_i2,my_i3)%name == 'FI      ') THEN
           CALL smooth_geopot (ds_grib, hsurf_tot, ie_tot, je_tot )
           CALL smoother(ds_grib, ie_tot, je_tot, 4,20 )
         ENDIF
       ENDIF

       ! Apply a digital smoother for selected fields
       IF ((     (var(my_i1,my_i2,my_i3)%name == 'PMSL     ')              &
            .OR. (var(my_i1,my_i2,my_i3)%name == 'PMSL_ANAI'))             &
                                     .AND. outblock%l_z_filter) THEN
         CALL smoother(ds_grib, ie_tot, je_tot, 4,20 )
       ENDIF
       IF(yextension == 'p' .AND. outblock%l_p_filter) THEN
         CALL smoother(ds_grib, ie_tot, je_tot, 4,20 )
       ENDIF
       IF(yextension == 'z' .AND. outblock%l_z_filter) THEN
         CALL smoother(ds_grib, ie_tot, je_tot, 4,20 )
       ENDIF

       ! Limit certain variables to the allowed range (again)
       SELECT CASE (var(my_i1,my_i2,my_i3)%name)
       CASE ('RELHUM    ')
         ds_grib(:) = MAX (0.0_irealgrib, MIN(100.0_irealgrib,ds_grib(:)))
       CASE ('QV        ','QC        ','QI        ',    &
             'QR        ','QS        ','QG        ')
         ds_grib(:) = MAX (0.0_irealgrib, ds_grib(:))
       END SELECT

       ! Now cut out the proper (sub-)field which was chosen with
       ! slon, slat, elon, elat
       ij_out = 0
       DO j = outblock%j_out_start, outblock%j_out_end
         DO i = outblock%i_out_start, outblock%i_out_end
           ij_out = ij_out+1
           ij_ful = (j-1) * ie_tot + i
           ds_out(ij_out) = ds_grib(ij_ful)
         ENDDO
       ENDDO

       IF (outblock%yform_write == 'grb1') THEN
#ifdef GRIBDWD
         ! degrib the level
         CALL grbex1(iednr, iz_ps, undefgrib, ndims, idims_out, ipds,     &
                  igds_out, ibms, ibds, ibmap, dsup, ds_out, iblock, ierrf)
         IF (ierrf /= 0) THEN
           yerrmsg = 'error in grbex1'
           CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
         ENDIF

         ! length of GRIB record in bytes
         irecord_len = idims_out(19)
#endif
       ELSE
 
         ! length of netcdf record in words
         irecord_len = outblock%ie_out_tot * outblock%je_out_tot
       ENDIF

     ELSE
       ! this is for restart-output
       CALL combine_subarrays (procarray_real, ds_real)

       ! length of restart record in words
       irecord_len = outblock%ie_out_tot * outblock%je_out_tot
     ENDIF

   ELSE

      irecord_len = 0_iintegers

   ENDIF

   !---------------------------------------------------------------------------
   !  Section 3.2: Check data and write output to disk
   !---------------------------------------------------------------------------

   ! check the data, if wanted
   IF (outblock%lcheck .AND. (.NOT. lrestart)) THEN
      my_k_f  = INT (my_k,  intgribf)
      my_iee  = INT (my_i2, intgribf)
      IF (.NOT. loutput) THEN
        my_i1 = 1; my_i2 = 1; my_i3 = 1           ! for safety reasons
      ENDIF
      CALL check_record (ds_grib, 1, ie_tot, 1, je_tot, 1, 1,               &
               outblock%i_out_start, outblock%i_out_end,                    &
               outblock%j_out_start, outblock%j_out_end, 1, 1,              &
               undefsub,  var(my_i1,my_i2,my_i3)%name,                      &
               my_iee, my_k_f, loutput, nuchkdat, num_compute,              &
               icomm_cart, my_cart_id, yerrmsg, izerror)
   ELSEIF (lrestart) THEN
      my_k_f  = INT (my_k,  intgribf)
      my_iee  = INT (my_i2, intgribf)
      IF (.NOT. loutput) THEN
        my_i1 = 1; my_i2 = 1; my_i3 = 1           ! for safety reasons
      ENDIF
      ds_grib(:) = REAL (ds_real(:), irealgrib)
      CALL check_record (ds_grib, 1, ie_tot, 1, je_tot, 1, 1, 1, ie_tot, 1, &
               je_tot, 1, 1, undefsub,  var(my_i1,my_i2,my_i3)%name,        &
               my_iee, my_k_f, loutput, nuchkdat, num_compute,              &
               icomm_cart, my_cart_id, yerrmsg, izerror)
   ENDIF

   IF (ltime) CALL get_timings (i_computations_O, ntstep, dt, izerror)

   SELECT CASE (outblock%yform_write)

   CASE ('grb1')

     CALL write_grib    (nuedat, iblock, irecord_len, iz_lfd, icomm_cart, &
                         num_compute, lflush, outblock%ydbtype,           &
                         lasync_io, yerrmsg, izerror)

   CASE ('ncdf')

     CALL write_netcdf  (nuedat, ds_out, outblock%ie_out_tot,             &
                         outblock%je_out_tot, irecord_len, my_orgdata,    &
                         icomm_cart, my_cart_id, num_compute, imp_grib,   &
                         lasync_io, yerrmsg, izerror)

   CASE ('bina')

     CALL write_restart (nuedat, ds_real, ie_tot, je_tot, irecord_len,    &
                         ipds, npds, igds_out, ngds, icomm_cart,          &
                         my_cart_id, num_compute, imp_reals, lasync_io,   &
                         yerrmsg, izerror)

   END SELECT

   IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 2035, yerrmsg, yroutine)
   ENDIF

   IF (ltime) CALL get_timings (i_write_data, ntstep, dt, izerror)

   IF ( ldebug_io .AND. (idbg_level >= 10)) THEN
     PRINT *, ' src_output: saved ', TRIM(var(my_i1, my_i2, my_i3)%name), my_k, ' to disk'
   ENDIF

   ! Reset organizational variables
   my_orgdata(:,:) = 0_iintegers
   loutput         = .FALSE.

ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE output_data

!==============================================================================
!+ Module procedure in src_output to create the grid definition block
!------------------------------------------------------------------------------

SUBROUTINE makegds

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure creates the grid definition block according to the
!   description of putgd1 (dwd) and the WMO.
!
! Method:
!
!==============================================================================
!
! Local scalars:
INTEGER  (KIND=iintegers)          :: k
REAL     (KIND=ireals)             :: pollon_grib
!

#ifdef GRIBDWD
INTEGER (KIND=intgribf), EXTERNAL  :: IREFTS

!- End of header
!==============================================================================

igds_out(:) = -999999

! length of igds_out in bytes
igds_out(1) = 42 + (4 + ke1) * 4

IF (ldwd_grib_use .AND. (.NOT.l_ke_in_gds)) THEN
  ! old style of coding the vertical coordinate parameters
  IF (ivctype == 3) THEN   ! SLEVE coordinates
    igds_out(1) = 42 + (8 + ke1) * 4
  ELSE                     ! NON SLEVE coordinates
    ! should logically be (5+ke1), but is kept for backward compatibility
    igds_out(1) = 42 + (4 + ke1) * 4
  ENDIF
  IF (irefatm == 2) igds_out(1) = 42 + (10 + ke1) * 4
  IF (irefatm == 3) igds_out(1) = 42 + (9 + ke1) * 4
ELSE
  ! new style of coding the vertical coordinate parameters
  IF (ivctype == 3) THEN   ! SLEVE coordinates
    igds_out(1) = 42 + (9 + ke1) * 4
  ELSE                     ! NON SLEVE coordinates
    igds_out(1) = 42 + (6 + ke1) * 4
  ENDIF
  IF (irefatm == 2) igds_out(1) = 42 + (11 + ke1) * 4
  IF (irefatm == 3) igds_out(1) = 42 + (10 + ke1) * 4
ENDIF

! number of the vertical coordinate parameters
IF (ldwd_grib_use .AND. (.NOT.l_ke_in_gds)) THEN
  ! old style of coding the vertical coordinate parameters
  IF (ivctype == 3) THEN   ! SLEVE coordinates
    igds_out(2) = 8 + ke1
  ELSE                     ! NON SLEVE coordinates
    ! should logically be 5+ke1, but is kept for backward compatibility
    igds_out(2) = 4 + ke1
  ENDIF
  IF (irefatm == 2) igds_out(2) = 10 + ke1
  IF (irefatm == 3) igds_out(2) = 9 + ke1
ELSE
  ! new style of coding the vertical coordinate parameters
  IF (ivctype == 3) THEN   ! SLEVE coordinates
    igds_out(2) = 9 + ke1
  ELSE                     ! NON SLEVE coordinates
    igds_out(2) = 6 + ke1
  ENDIF
  IF (irefatm == 2) igds_out(2) = 11 + ke1
  IF (irefatm == 3) igds_out(2) = 10 + ke1
ENDIF

! location of the list of vertical coordinate parameters in bytes
igds_out(3) = 43

! data representation type
igds_out(4) = 10

! calculation of the left bottom corner and
! left bottom and right upper corner in millidegrees:
! This depends on the variable (U- and V-related variables are shifted in
! the Arakawa C-grid) and is determined during the output step.

igds_out( 9) = 8       ! this was 0 before; but due to Grib 1, Code Table 7
                       ! bit no. 5 must be set to 1, to indicate that u,v
                       ! components are relative to defined (rotated) grid

! increments
igds_out(12) = 0
igds_out(13) = 0


igds_out(14) = 64
igds_out(15:19) = 0

! coordinates of the pole: in the grib code the southern pole has to be
! specified: for the latitude this is the negative value, for the
! longitude it is the value + 180.0 and then limiting the value again
! to -180.0 ... 180.0
pollon_grib = pollon + 180.0_ireals
IF (pollon_grib > 180.0_ireals) THEN
  pollon_grib = pollon_grib - 360.0_ireals
ENDIF

igds_out(20) = NINT(-pollat      * 1000.0_ireals)
igds_out(21) = NINT( pollon_grib * 1000.0_ireals)
!US igds_out(22) = IREFTS(0.0_irealgrib)
igds_out(22) = IREFTS(REAL(polgam, irealgrib))

! vertical coordinate parameters
IF (ldwd_grib_use .AND. (.NOT.l_ke_in_gds)) THEN
  ! old style of coding: first 4 values for the reference atmosphere,
  ! then the vertical coordinate parameters and last eventually
  ! additional parameters for the SLEVE coordinate
  igds_out(26) = IREFTS( REAL (p0sl,   irealgrib) )
  igds_out(27) = IREFTS( REAL (t0sl,   irealgrib) )
  igds_out(28) = IREFTS( REAL (dt0lp,  irealgrib) )
  igds_out(29) = IREFTS( REAL (vcflat, irealgrib) )

  DO k = 1, ke+1
    igds_out(29 + k) = IREFTS( REAL (vcoord(k), irealgrib) )
  ENDDO

  IF (ivctype == 3) THEN
    ! SLEVE coordinate: 3 new parameters are written to GDS
    igds_out(29+ke1+1) = ivctype
    igds_out(29+ke1+2) = IREFTS (REAL(svc1,    irealgrib))
    igds_out(29+ke1+3) = IREFTS (REAL(svc2,    irealgrib))
    igds_out(29+ke1+4) = IREFTS (REAL(nfltvc,  irealgrib))
  ENDIF

  IF (irefatm == 2) THEN ! Write parameters for new reference atmosphere
    igds_out(29+ke1+1) = ivctype+100
    IF (ivctype /= 3) THEN
      igds_out(29+ke1+2) = IREFTS (0.0_irealgrib)
      igds_out(29+ke1+3) = IREFTS (0.0_irealgrib)
      igds_out(29+ke1+4) = IREFTS (0.0_irealgrib)
    ENDIF
    igds_out(29+ke1+5) = IREFTS (REAL(delta_t, irealgrib))
    igds_out(29+ke1+6) = IREFTS (REAL(h_scal,  irealgrib))
  ELSE IF (irefatm == 3) THEN ! Write parameters for constant-BV reference atmosphere
    igds_out(29+ke1+1) = ivctype+200
    IF (ivctype /= 3) THEN
      igds_out(29+ke1+2) = IREFTS (0.0_irealgrib)
      igds_out(29+ke1+3) = IREFTS (0.0_irealgrib)
      igds_out(29+ke1+4) = IREFTS (0.0_irealgrib)
    ENDIF
    igds_out(29+ke1+5) = IREFTS (REAL(bvref, irealgrib))
  ENDIF
ELSE
  ! new style of coding: the first 2 values are the vertical coordinate
  ! type and the number of main levels used. Then 4 values for the 
  ! reference atmosphere, the vertical coordinate parameters and 
  ! last eventually additional parameters for the SLEVE coordinate
  IF (irefatm == 1) THEN
    igds_out(26) = IREFTS( REAL (ivctype,irealgrib) )
  ELSE IF (irefatm == 2) THEN
    igds_out(26) = IREFTS( REAL (ivctype+100,irealgrib) )
  ELSE IF (irefatm == 3) THEN
    igds_out(26) = IREFTS( REAL (ivctype+200,irealgrib) )
  ENDIF
  igds_out(27) = IREFTS( REAL (ke,     irealgrib) )
  igds_out(28) = IREFTS( REAL (p0sl,   irealgrib) )
  igds_out(29) = IREFTS( REAL (t0sl,   irealgrib) )
  igds_out(30) = IREFTS( REAL (dt0lp,  irealgrib) )
  igds_out(31) = IREFTS( REAL (vcflat, irealgrib) )

  DO k = 1, ke+1
    igds_out(31 + k) = IREFTS( REAL (vcoord(k), irealgrib) )
  ENDDO

  IF (ivctype == 3) THEN
    ! SLEVE coordinate: 3 new parameters are written to GDS
    igds_out(31+ke1+1) = IREFTS (REAL(svc1,  irealgrib))
    igds_out(31+ke1+2) = IREFTS (REAL(svc2,  irealgrib))
    igds_out(31+ke1+3) = IREFTS (REAL(nfltvc,irealgrib))
  ENDIF

  IF (irefatm == 2) THEN ! Write parameters for new reference atmosphere
    IF (ivctype /= 3) THEN
      igds_out(31+ke1+1) = IREFTS (0.0_irealgrib)
      igds_out(31+ke1+2) = IREFTS (0.0_irealgrib)
      igds_out(31+ke1+3) = IREFTS (0.0_irealgrib)
    ENDIF
    igds_out(31+ke1+4) = IREFTS (REAL(delta_t, irealgrib))
    igds_out(31+ke1+5) = IREFTS (REAL(h_scal,  irealgrib))
  ELSE IF (irefatm == 3) THEN ! Write parameters for constant-BV reference atmosphere
    IF (ivctype /= 3) THEN
      igds_out(31+ke1+1) = IREFTS (0.0_irealgrib)
      igds_out(31+ke1+2) = IREFTS (0.0_irealgrib)
      igds_out(31+ke1+3) = IREFTS (0.0_irealgrib)
    ENDIF
    igds_out(31+ke1+4) = IREFTS (REAL(bvref, irealgrib))
  ENDIF
ENDIF
#endif

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE makegds

!==============================================================================
!+ Module procedure in src_output to create the product definition block
!------------------------------------------------------------------------------

SUBROUTINE makepds ( n1,n2,n3, nlevel, ntri, nuot, nlastout, nprocess_ini_out,&
                     nprocess_bd_out, lanalysis, yextension, slev, lrestart)

!------------------------------------------------------------------------------
!
! Description:
!   This module procedure creates the product definition block according to the 
!   description of putpd1 (dwd) and the WMO.
!
! Method:
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER  (KIND=iintegers) , INTENT(IN) :: n1, n2, n3, nlevel, nlastout, nuot
INTEGER  (KIND=iintegers) , INTENT(IN) :: nprocess_ini_out, nprocess_bd_out
REAL     (KIND=ireals)    , INTENT(IN) :: slev 
CHARACTER (LEN=1),          INTENT(IN) :: yextension 
LOGICAL                   , INTENT(IN) :: lanalysis, lrestart

INTEGER  (KIND=intgribf)  , INTENT(OUT):: ntri
!
! Local scalars:
INTEGER  (KIND=iintegers)              :: ndgb, mmgb, jjgb, nhgb, jcgb
REAL     (KIND=ireals)                 :: hm, zacthour
INTEGER  (KIND=iintegers)              :: jj,mm,dd,hh,mi, nzentry
INTEGER  (KIND=iintegers)              :: timestep, nzjulianday
CHARACTER (LEN=  8)                    :: ydate
CHARACTER (LEN= 10)                    :: yzdat1, ytime
CHARACTER (LEN= 22)                    :: yzdat2

!
!- End of header
!==============================================================================

  timestep = ntstep

  IF ( yextension == 's') THEN
    nzentry = INT(slev)
  ENDIF

  ! Initialisation
  ipds(:)  = -999999
   
  ! grib table nr. of version
!  SELECT CASE (n3)
!  CASE (1)
!    ipds(2) =   2_intgribf
!  CASE (2)
!    ipds(2) = 201_intgribf
!  CASE (3)
!    ipds(2) = 202_intgribf
!  CASE (4)
!    ipds(2) = 203_intgribf
!  CASE (5)
!    ipds(2) = 205_intgribf
!  END SELECT
   ipds(2)  = lst_gribtabs(n3)
   
  ! centre identification
  ipds(3)  = ncenter
   
  ! generating process identification
  IF (ldwd_grib_use) THEN
    IF (lanalysis) THEN
      IF (nprocess_ini_out == -999999) THEN
        ipds(4) =  nprocess_ini_in
      ELSE
        ipds(4) =  nprocess_ini_out
      ENDIF
    ELSE
      IF (nprocess_bd_out == -999999) THEN
        ipds(4) =  nprocess_bd_in
      ELSE
        ipds(4) =  nprocess_bd_out
      ENDIF
    ENDIF
  ENDIF
   
  ! number of used grid
  ipds(5)  = 255_intgribf
   
  ! flag, that indicates whether GDS and BMS follows
  ipds(6)  = 128_intgribf
   
  ! number of element
  ipds(7)  = INT (n2, intgribf)
   
  IF ( yextension == 's') THEN
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    ! satellite channels on output
    ! leveltyp and value of level
    ipds(8)  = 222_intgribf
    ipds(9)  = 0_intgribf
    ipds(10) = sat_compute(nzentry)%ngrib_chan(nlevel)
#endif
  ELSEIF ( yextension == 'p' ) THEN
    ! pressure level on output
    ! leveltyp and value of level
    ipds(8)  = 100_intgribf
    ipds(9)  = 0_intgribf
    ipds(10) = NINT(slev*0.01_ireals)
  ELSE IF ( yextension == 'z') THEN 
    ! z- levels on output
    ipds(8)  = 103_intgribf
    ipds(9)  = 0_intgribf
    ipds(10) = NINT(slev)
  ELSE
    ! model-levels on output
    ! leveltyp
    ipds(8)  = var(n1,n2,n3)%levtyp
   
    ! upper and lower boundary of the level
      
    SELECT CASE(var(n1,n2,n3)%levtyp)
    CASE (1,2,3,4,8,100,102,103,107,113,115,117,119,125,160,200,201)
      ipds(9)  = 0_intgribf
      ipds(10) = 0_intgribf
    CASE (109)
      ipds(9)  = 0_intgribf
      ipds(10) = nlevel
    CASE (105,111)
      IF ( (var(n1,n2,n3)%name == 'ALHFL_PL  ') .OR. &
           (var(n1,n2,n3)%name == 'W_SO      ') .OR. &
           (var(n1,n2,n3)%name == 'T_SO      ') .OR. &
           (var(n1,n2,n3)%name == 'W_SO_ICE  ') ) THEN
        ipds(9)  = 0.0_ireals        ! levtop
        ipds(10) = msoilgrib(nlevel) ! depth of main soil level in cm
      ELSE
        ipds(9)  = 0_intgribf
        ipds(10) = var(n1,n2,n3)%levbot
      ENDIF
    CASE (112)
      ipds(9)  = var(n1,n2,n3)%levtop
      ipds(10) = var(n1,n2,n3)%levbot
    CASE DEFAULT
      ipds(9)  = INT (nlevel  , intgribf)
      ipds(10) = INT (nlevel+1, intgribf)
    END SELECT 
  ENDIF
   
  ! Reference time of data
  IF (lanalysis) THEN
    IF (ldwd_grib_use) THEN
      ! When nudging is active the output fields are treated as analyses
      ! and the reference time is the actual forecast time
      ! When using a flexible dt, the hour has eventually to be updated
      ! to the nearest full hour

      CALL get_utc_date(ntstep, ydate_ini, dt, itype_calendar, yzdat1,    &
                        yzdat2, nzjulianday, zacthour)
      READ(yzdat1,'(5I2)') jcgb, jjgb, mmgb, ndgb, nhgb
      nhgb = NINT(zacthour, iintegers)
    ELSE
      ! Local modification for aLMo verification.     gdm, May 4, 2004.
      ! Use ydate_ini as reference time for ACCUM/MIN/MAX/MEAN fields
      ! for aLMo lanalysis
      SELECT CASE(var(n1,n2,n3)%ntri)
      CASE(0,1)
        ! When nudging is active the output fields are treated as analyses
        ! and the reference time is the actual forecast time
        READ(yakdat1  ,'(5I2)') jcgb, jjgb, mmgb, ndgb, nhgb
      CASE DEFAULT
        READ(ydate_ini,'(5I2)') jcgb, jjgb, mmgb, ndgb, nhgb
      END SELECT
    ENDIF
  ELSE
    READ(ydate_ini,'(5I2)') jcgb, jjgb, mmgb, ndgb, nhgb
  ENDIF

  IF(jjgb == 0) THEN
    ipds(11) = 100_intgribf
    ipds(22) = INT (jcgb, intgribf)
  ELSE  
    ipds(11) = INT (jjgb, intgribf)
    ipds(22) = INT (jcgb, intgribf) + 1
  ENDIF
  ipds(12) = INT (mmgb, intgribf)
  ipds(13) = INT (ndgb, intgribf)
  ipds(14) = INT (nhgb, intgribf)
  ipds(15) = 0_intgribf    ! hier muessen die Minuten rein
   
  ! Indicator of unit of time range
  ! (hm gives the seconds for this range
  ipds(16) = INT (nuot, intgribf)
  SELECT CASE (nuot)
  CASE ( 0)                  !  1 minute
    hm  =    60.0_ireals
  CASE ( 1)                  !  1 hour
    hm  =  3600.0_ireals
  CASE ( 2)                  !  1 day
    hm  = 86400.0_ireals
  CASE (10)                  !  3 hours
    hm  = 10800.0_ireals
  CASE (11)                  !  6 hours
    hm  = 21600.0_ireals
  CASE (12)                  ! 12 hours
    hm  = 43200.0_ireals
  CASE (13)                  ! 15 minutes
    hm  =   900.0_ireals
  CASE (14)                  ! 30 minutes
    hm  =  1800.0_ireals
  CASE (15)                  ! 10 minutes
    hm  =   600.0_ireals
  CASE DEFAULT
    PRINT *, ' ERROR *** Wrong value for unit-of-time:  ', nuot
    CALL model_abort (my_cart_id, 2015, 'wrong unit-of-time', 'makepds')
  END SELECT

  ! Time Range Indicator
  IF( timestep == 0 ) THEN
    IF (ldfi) THEN
      ntri = 1_intgribf   ! this means initialized analyses. 
    ELSE
      ntri = 0_intgribf   ! this means uninitialized analyses. 
    ENDIF
  ELSE
    IF (ldwd_grib_use) THEN
      IF (lanalysis) THEN ! analysis cannot be written in time step 0
        ntri = 13_intgribf
      ELSE
        ntri = var(n1,n2,n3)%ntri
      ENDIF
    ELSE
      ntri = var(n1,n2,n3)%ntri
    ENDIF
  ENDIF
   
  ! Time p1 and p2:  
  !  p1:       actual time for which a forecast product is valid
  !            last deletion of min, max-values (tri=2)
  !       or:  begin of averaging period        (tri=3)
  !  p2:       actual time for which a forecast product is valid
  !       or:  end of averaging period          (tri=3)
  SELECT CASE(ntri)
  CASE(0,1)
    IF (ldwd_grib_use) THEN
      ipds(17) = NINT (timestep*dt/hm)
      ipds(18) = 0_intgribf
    ELSE
      IF (lanalysis) THEN
        ipds(17) = 0_intgribf
        ipds(18) = 0_intgribf
      ELSE
        ipds(17) = NINT (timestep*dt/hm)
        ipds(18) = 0_intgribf
      ENDIF
    ENDIF
  CASE(2)
    IF(var(n1,n2,n3)%name == 'VMAX_10M'  .OR. &
       var(n1,n2,n3)%name == 'VGUST_DYN' .OR. &
       var(n1,n2,n3)%name == 'VGUST_CON') THEN
      ipds(17) = NINT (nlastmxu * dt/hm, intgribf)
    ELSE   ! tmin_2m, tmax_2m
      ipds(17) = NINT (nlastmxt * dt/hm, intgribf)
    ENDIF
    ipds(18) = NINT (timestep*dt/hm, intgribf)
  CASE(3)
    IF (lbdclim) THEN
      ipds(17) = NINT (nlastout*dt/hm, intgribf)
    ELSE
      ipds(17) = 0_intgribf
    ENDIF
!CPS#ifdef COUP_OAS_COS
!CPS    ipds(17) = NINT (nlastout*dt/hm, intgribf)   !CPS sanity
!CPS#endif
    ipds(18) = NINT (timestep*dt/hm, intgribf)
  CASE(4)
    IF (lbdclim) THEN
      ipds(17) = NINT (nlastout*dt/hm, intgribf)
    ELSE
      ipds(17) = 0_intgribf
    ENDIF
    ipds(18) = NINT (timestep*dt/hm, intgribf)
  CASE(10)
    ipds(18) = NINT (timestep*dt/hm, intgribf)
  CASE(13)
    ! assimilation runs at DWD
    ipds(17) = 0_intgribf
    ipds(18) = 0_intgribf

    IF (var(n1,n2,n3)%ntri >= 2) THEN
      ipds(18) = NINT (timestep*dt/hm, intgribf)
    ENDIF

    IF (var(n1,n2,n3)%name(1:8) == 'VMAX_10M' .OR. &
        var(n1,n2,n3)%name(1:9) == 'VGUST_DYN' .OR. &
        var(n1,n2,n3)%name(1:9) == 'VGUST_CON') THEN
      ipds(18) = ipds(18)-NINT (nlastmxu * dt/hm, intgribf)
    ENDIF
    IF ( (var(n1,n2,n3)%name(1:7) == 'TMIN_2M')   .OR.     &
         (var(n1,n2,n3)%name(1:7) == 'TMAX_2M') ) THEN
      ipds(18) = ipds(18)-NINT (nlastmxt * dt/hm, intgribf)
    ENDIF
  END SELECT

  ! Check, if ipds(17,18) are between 0 and 254
  ! (but for restart-files in climate mode this will be violated, 
  !  therefore omit the warnings; the same might be true for very
  !  high resolution artificial cases with a high output frequency)
  IF (.NOT. lbdclim .AND. .NOT. lartif_data) THEN
    IF (.NOT. lrestart) THEN
      IF ( (ipds(17) < 0_intgribf) .OR. (ipds(17) > 254_intgribf) ) THEN
        PRINT *, ' WARNING: *** ipds(17) = ', ipds(17),      &
                 ' is outside the valid range *** ', n1, n2, n3, nuot, ntri, hm
      ENDIF
      IF ( (ipds(18) < 0_intgribf) .OR. (ipds(18) > 254_intgribf) ) THEN
        PRINT *, ' WARNING: *** ipds(18) = ', ipds(18),      &
                 ' is outside the valid range *** ', n1, n2, n3, nuot, ntri, hm
      ENDIF
    ENDIF
  ENDIF

  ipds(19) = ntri
  ipds(20) = 0_intgribf
  ipds(21) = 0_intgribf
  ipds(24) = 0_intgribf
   
  ! local use area
  IF (leps) THEN
    idims_out(11) = 54_intgribf
    ipds( 1)    = 66_intgribf
    ipds(37)    = 253_intgribf
    ipds(48:49) = 0_intgribf
    ipds(53:54) = 0_intgribf
  ELSE
    ipds(37)    = 254_intgribf
  ENDIF
  ipds(38:40) = 0_intgribf

  IF ( yextension == 's') THEN
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
    ipds(41) = sat_compute(nzentry)%ngrib_aees(nlevel)
#endif
  ELSE
    SELECT CASE(nsma_stat)
      CASE(21,23,27)
        ipds(41)    = INT (nsma_stat, intgribf)
      CASE DEFAULT
        ipds(41)    = 0_intgribf
    END SELECT
  ENDIF

  ! write actual date and time to local use area of pds
  CALL DATE_AND_TIME(ydate,ytime)
  READ(ydate,'(I4,2I2)') jj,mm,dd
  READ(ytime,'(2I2,6X)') hh,mi
  ! According to DWD standard, ipds(42) contains the offset to the year 1900
  ipds(42)    = INT (jj-1900, intgribf)
  ipds(43)    = INT (mm, intgribf)
  ipds(44)    = INT (dd, intgribf)
  ipds(45)    = INT (hh, intgribf)
  ipds(46)    = INT (mi, intgribf)

  ipds(47)    = INT (nvers, intgribf)

  IF (leps) THEN
    ipds(50)   = INT (iepstyp, intgribf)
    ipds(51)   = INT (iepstot, intgribf)
    ipds(52)   = INT (iepsmem, intgribf)
  ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE makepds

!==============================================================================
!+ Module procedure in src_output for the p-interpolation
!------------------------------------------------------------------------------

SUBROUTINE p_int (outblock, i1,i2,i3, nlist, itlev, results)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine interpolates variables given in the namelist from
!   model leves to pressure levels. The result of the interpolation
!   is given back to the calling procedure in the three dimensional variable
!   "results". 
!
! Method:
!   Column wise interpolation with tension splines.
!   For the interpolation of FI and QV the logarithm of the pressure is
!   used.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
TYPE(pp_nl),              INTENT(IN) ::     &
  outblock       ! pointer to the namelist group

INTEGER (KIND=iintegers), INTENT(IN) ::     &
  i1,i2,i3,      & ! location of the variable to be processed in the LM 
                 ! variable table
  nlist,       & ! location of the variable in the output list
  itlev          ! time level of the variables

! Array arguments with intent(out):
REAL (KIND=ireals)      , INTENT(OUT)::     &
  results(ie,je,outblock%kepin)

!------------------------------------------------------------------------------
!
! Local parameters
REAL(KIND=ireals),PARAMETER    :: gamma   = 5.5_ireals    ! tension factor
REAL(KIND=ireals),PARAMETER    :: delpchk = 1.0_ireals

! Local scalars:
INTEGER (KIND=iintegers)       :: i, j, k, kint(ie,je)
INTEGER (KIND=iintegers)       :: ierrstat, ierr, izlen
INTEGER (KIND=iintegers)       :: nldim(ie,je)

REAL(KIND=ireals)              :: zt0s, zalnp, zpexp

CHARACTER (LEN=25)             :: yroutine
CHARACTER (LEN=80)             :: yerrmsg
CHARACTER (LEN=10)             :: yzname

! Local arrays:
REAL(KIND=ireals)              :: fmfl(ie,je,ke)
REAL(KIND=ireals)              :: fexp(ie,ke+4),       &
                                  pexp(ie,ke+4)
REAL(KIND=ireals)              :: zdelp(ie,je)
REAL(KIND=ireals)              :: ztstar(ie,je), zalpha(ie,je), zt0(ie,je)

! Output from tautsp
REAL(KIND=ireals)              :: s_vec    (ie,(ke+4)*6), &
                                  break_vec(ie,(ke+4)*3), &
                                  coef_vec (ie,4,(ke+4)*3)

! Output from spline
REAL(KIND=ireals)              :: fpfls(outblock%kepin),  &
                                  pfls (outblock%kepin)

! Switch to choose between vertical interpolation methods:
! Local value; there is a global namelist parameter itype_vertint
! in each namelist of group GRIBOUT.
INTEGER(KIND=iintegers) :: zitype_vertint

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Set variable that has to be interpolated
!------------------------------------------------------------------------------

yroutine = 'p_int'
yzname   = '          '
yzname(1:LEN_TRIM(outblock%yvarpl(nlist))) =                            &
         outblock%yvarpl(nlist)(1:LEN_TRIM(outblock%yvarpl(nlist)))
izlen = LEN_TRIM(yzname)

! Spline interpolation is the default value, which may be changed for single fields below
! (e.g., if you would like to have splineinterpolation for all fiels but not for QC ...)
zitype_vertint = outblock%itype_vertint

! calculation of u and v on masspoint
! -----------------------------------

  IF     (yzname(1:izlen) == 'U' ) THEN
    fmfl(2:ie,:,1:ke) = 0.5_ireals * (var(i1,i2,i3)%p4(2:ie,:,:,itlev) &
                                    + var(i1,i2,i3)%p4(1:ie-1,:,:,itlev))  
    fmfl(1,:,1:ke)    = fmfl(2,:,1:ke)
  ELSEIF (yzname(1:izlen) == 'V' ) THEN
!CDIR COLLAPSE
    fmfl(:,2:je,1:ke) = 0.5_ireals * (var(i1,i2,i3)%p4(:,2:je,:,itlev) &
                                 + var(i1,i2,i3)%p4(:,1:je-1,:,itlev))
    fmfl(:,1,1:ke)    = fmfl(:,2,1:ke)
  ELSEIF (yzname(1:izlen) == 'TKVM' .OR. yzname(1:izlen) == 'TKVH') THEN
!CDIR COLLAPSE
    fmfl(:,:,2:ke)    = var(i1,i2,i3)%p3(:,:,2:ke)
    fmfl(:,:,1)       = fmfl(:,:,2)
  ELSEIF (yzname(1:izlen) == 'FI' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    =  0.5_ireals * ( hhl(:,:,1:ke) + hhl(:,:,2:ke+1) ) * g
  ELSEIF (yzname(1:izlen) == 'Q_SEDIM') THEN
    IF (lprog_qi) THEN
!CDIR COLLAPSE
      fmfl(:,:,1:ke) = qrs(:,:,:) - qi(:,:,:,itlev)
    ELSE
!CDIR COLLAPSE
      fmfl(:,:,1:ke) = qrs(:,:,:)
    ENDIF
  ELSEIF (yzname(1:izlen) == 'RELHUM') THEN
    CALL calrelhum(fmfl(:,:,:), t(:,:,:,itlev), pp(:,:,:,itlev), p0(:,:,:), &
                   qv(:,:,:,itlev),ie, je, ke, b1, b2w, b3, b4w, rdv, o_m_rdv)
  ELSEIF (yzname(1:izlen) == 'OMEGA') THEN
    CALL calomega (fmfl(:,:,:), pp(:,:,:,nnew), pp(:,:,:,itlev), pptens(:,:,:),&
                   w(:,:,:,itlev), rho0(:,:,:), ie, je, ke, dt, g )
  ELSEIF (yzname(1:izlen) == 'FI_ANAI' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = p_anai(:,:,1:ke) / rho(:,:,1:ke)
  ELSEIF (yzname(1:izlen) == 'DBZ') THEN
    IF (itype_gscp == 3) THEN
      CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,         &
        klv850, my_cart_id, itype_gscp, t(:,:,:,itlev), qc(:,:,:,itlev)*rho, qr(:,:,:,itlev)*rho,   &
        qi(:,:,:,itlev)*rho, qs(:,:,:,itlev)*rho, z_radar = fmfl )
    ELSEIF (itype_gscp == 4) THEN
      CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,         &
        klv850, my_cart_id, itype_gscp, t(:,:,:,itlev), qc(:,:,:,itlev)*rho, qr(:,:,:,itlev)*rho,   &
        qi(:,:,:,itlev)*rho, qs(:,:,:,itlev)*rho, q_grau  = qg(:,:,:,itlev)*rho,       &
        z_radar = fmfl )
    ENDIF
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
    ! Reflectivity interpolation should be done in linear space, not in logarithmic (as alternative, interpolation
    ! in the space of rain rate -- Z^0.6666 --- would also be desireable):
    WHERE (fmfl(:,:,:) >= -90.0)
      fmfl(:,:,:) = 10.0_ireals ** (0.1_ireals * fmfl(:,:,:))
    ELSEWHERE
      fmfl(:,:,:) = 0.0_ireals
    END WHERE
  ELSE
    SELECT CASE( var(i1,i2,i3)%rank )
    CASE(4)
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = var(i1,i2,i3)%p4(:,:,1:ke,itlev)
    CASE(3)
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = var(i1,i2,i3)%p3(:,:,1:ke)
    END SELECT
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Interpolation
!------------------------------------------------------------------------------

! slicewise interpolation
! -----------------------

  DO j=jstartpar,jendpar

    ! Calculation of the pressure at full model levels.
    ! Variables defined on half levels are first averaged to
    ! full levels.
    ! For water vapour (QV) and geopotential height (FI), the 
    ! interpolation on constatnt pressure levels is logarithmic 
    ! with respect to pressure.

    SELECT CASE(var(i1,i2,i3)%levtyp)
    CASE(110)
      IF (yzname(1:izlen) == 'QV'       .OR.                     &
          yzname(1:izlen) == 'FI'       .OR.                     &
          yzname(1:izlen) == 'QV_ANAI'  .OR.                     &
          yzname(1:izlen) == 'FI_ANAI' ) THEN
        DO k = 1, ke
          DO i = 1, ie
            pexp(i,k) = LOG(p0(i,j,k) + pp(i,j,k,itlev))
          ENDDO
        ENDDO
        pfls(1:outblock%kepin) = LOG(outblock%plev(1:outblock%kepin))
      ELSE
        DO k = 1, ke
          DO i = 1, ie
            pexp(i,k)  = p0(i,j,k) + pp(i,j,k,itlev)
          ENDDO
        ENDDO
        pfls(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
      ENDIF
    CASE(109)
      DO k = 2, ke
        DO i = 1, ie
          pexp(i,k) = p0hl(i,j,k)  &
                  + 0.5_ireals*(pp(i,j,k,itlev)+pp(i,j,k-1,itlev))
        ENDDO
      ENDDO
      DO i = 1, ie
        pexp(i,1) = p0hl(i,j,1) + pp(i,j,1,itlev)
             ! best approximation we can get there
             ! before, this was used, but factor 0.5 is wrong
             !                + 0.5_ireals *pp(i,j,1,itlev)
      ENDDO
      pfls(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
    CASE DEFAULT
               yerrmsg = 'wrong leveltyp of input field'
               ierrstat = 2004
               CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    END SELECT


    ! calculate pressure difference of lowest constant pressure level
    ! and surface pressure

    zdelp(:,j) = outblock%plev(outblock%kepin) - ps(:,j,itlev)

    ! vertical interpolation colum wise
    ! ---------------------------------

    ! copy slice variable to colum variable and set lowest modellevel
    DO k = 1, ke
      DO i=istartpar,iendpar
        fexp(i,k)  = fmfl(i,j,k)
      ENDDO
    ENDDO
    DO i=istartpar,iendpar
      fexp(i,ke+1)  = fexp(i,ke)
      pexp(i,ke+1)  = ps(i,j,itlev)

      IF (zdelp(i,j) > 0.0) THEN
        kint(i,j) = ke + 4
      ELSE
        kint(i,j) = ke + 1
      ENDIF
      nldim(i,j) = (ke+4)*3
    ENDDO

    IF (TRIM(yzname) == 'T') THEN
      ! calculation of the soil temperatur (ztstar) from the temperatur 
      ! of the lowest modell level
      DO i=istartpar,iendpar
        ztstar(i,j) = t(i,j,ke,itlev) + 0.0065_ireals                       &
                      * 0.5_ireals*(hhl(i,j,ke)-hhl(i,j,ke+1))
        zalpha(i,j) = 0.0065*r_d/g
        zt0   (i,j) = ztstar(i,j)  + 0.0065 * hsurf(i,j)
      ENDDO

      DO i=istartpar,iendpar
        IF (ztstar(i,j) > 315.0 .AND. zt0(i,j) > 315.0) THEN
          ztstar(i,j) = 0.5 * ( 315.0 + ztstar(i,j))
        ELSEIF (ztstar(i,j) < 255.0) THEN
          ztstar(i,j) = 0.5 * ( 255.0 + ztstar(i,j))
        ENDIF
        IF(hsurf(i,j) > 2500_ireals) THEN
          zt0s   = MIN(zt0(i,j),298.0_ireals)
          zalpha(i,j) = r_d * (zt0s - ztstar(i,j))/(hsurf(i,j)*g)
          IF(zt0s < ztstar(i,j)) zalpha(i,j) = 0_ireals
        ELSEIF((hsurf(i,j) <= 2500_ireals)                        &
          .AND.(hsurf(i,j) >= 2000_ireals)) THEN
          zt0s   = 0.002_ireals * (2500_ireals-hsurf(i,j) * zt0(i,j)  &
                   +(hsurf(i,j)-2000_ireals)*MIN(zt0(i,j),298.0_ireals))
          zalpha(i,j) = r_d * (zt0s - ztstar(i,j))/(hsurf(i,j)*g)
          IF(zt0s < ztstar(i,j)) zalpha(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDIF

    IF (TRIM(yzname) == 'FI') THEN
      ! calculation of the soil temperatur (ztstar) from the temperatur 
      ! of the lowest modell level
      DO i=istartpar,iendpar
        ztstar(i,j) = t(i,j,ke,itlev) + 0.0065_ireals                       &
                      * 0.5_ireals*(hhl(i,j,ke)-hhl(i,j,ke+1))
        zalpha(i,j) = 0.0065*r_d/g
        zt0   (i,j) = ztstar(i,j)  + 0.0065 * hsurf(i,j)
      ENDDO

      DO i=istartpar,iendpar
        IF (ztstar(i,j) <= 290.5 .AND. zt0(i,j) > 290.5) THEN
          zalpha(i,j) = r_d*(290.5 - ztstar(i,j))/(hsurf(i,j)*g)
        ELSEIF (ztstar(i,j) > 290.5 .AND. zt0(i,j) > 290.5) THEN
          ztstar(i,j) = 0.5*(290.5 + ztstar(i,j))
          zalpha(i,j) = 0.0
        ELSEIF (ztstar(i,j) < 255.0) THEN
          ztstar(i,j) = 0.5*(255.0 + ztstar(i,j))
        ENDIF
      ENDDO
    ENDIF


    IF (yzname(1:izlen) == 'T' ) THEN
      DO i=istartpar,iendpar
        IF (zdelp(i,j) > 0.0) THEN
          ! extrapolation of the model data if the required pressure level
          ! is below the lowest modellevel
          DO k=ke+1,ke+4
            pexp(i,k) = ps(i,j,itlev) &
                    + (zdelp(i,j) + delpchk)/4.0_ireals* REAL (k-ke, ireals)
            zalnp   = zalpha(i,j) * LOG(pexp(i,k)/ps(i,j,itlev))
            fexp(i,k) = ztstar(i,j) *                                     &
                      ( 1.0_ireals + zalnp + 0.5*zalnp**2 + 0.1667*zalnp**3)
          ENDDO
        ENDIF
      ENDDO
    ELSEIF (yzname(1:izlen) == 'FI' ) THEN
      DO i=istartpar,iendpar
        IF (zdelp(i,j) > 0.0) THEN
          DO k=ke+1,ke+4
            zpexp   = ps(i,j,itlev) &
                     + (zdelp(i,j) + delpchk)/4.0_ireals* REAL (k-ke, ireals)
            zalnp   = zalpha(i,j) * LOG(zpexp/ps(i,j,itlev))
            fexp(i,k) = hsurf(i,j) * g                                   &
                          -r_d*ztstar(i,j) * LOG(zpexp/ps(i,j,itlev))*   &
                                (1.0_ireals+zalnp                        &
                                +0.5_ireals*zalnp**2                     &
                           +0.1667_ireals*zalnp**3)
            pexp(i,k) = LOG(zpexp)
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO i=istartpar,iendpar
        IF (zdelp(i,j) > 0.0) THEN
          DO k=ke+1,ke+4
            pexp(i,k) = ps(i,j,itlev) &
                    + (zdelp(i,j) + delpchk)/4.0_ireals* REAL (k-ke, ireals)
            fexp(i,k) = fexp(i,ke)
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    SELECT CASE (zitype_vertint)

    CASE (1)

  ! Spline interpolation over (i,k)-slices
  ! --------------------------------------

      CALL tautsp2D (pexp(:,:), fexp(:,:), kint(:,j) ,ie, istartpar,     &
           iendpar, ke+4, gamma, s_vec, break_vec, coef_vec, nldim, ierr)
      IF (ierr == 2) THEN
        yerrmsg = 'wrong input in tautsp'
        ierrstat = 1004
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF
      IF (ierr /= 0) THEN
        yerrmsg = ' ERROR *** Error in tautsp *** '
        ierrstat = 1005
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF
      
      CALL spline2D(break_vec,coef_vec,j,nldim,pfls,outblock%kepin,    &
           results, outblock%yvarpl(nlist))
      
    CASE (2)

      ! .. provide monotonically increasing dummy values for pexp(:,ke+1:ke+4)
      !    below the surface, because lininterp2D_xinter1D_vec()
      !    requires that all values of pexp are prescribed:
      DO i=istartpar,iendpar
        IF (zdelp(i,j) <= 0.0) THEN
          DO k=ke+1,ke+4
            pexp(i,k) = ps(i,j,itlev) &
                    + delpchk/4.0_ireals* REAL (k-ke, ireals)
            fexp(i,k) = fexp(i,ke)
          ENDDO
        ENDIF
      ENDDO

      ! Linear Interpolation with respect to pressure:
      ierr = 0
      CALL lininterp2D_xinter1D_vec(pexp(istartpar:iendpar,:), &
           fexp(istartpar:iendpar,:), iendpar-istartpar+1, ke+4, &
           pfls(:), results(istartpar:iendpar,j,:), outblock%kepin, ierr)
      
      IF (ierr /= 0) THEN
        yerrmsg = ' ERROR *** Error in lininterp2D_xinter1D_vec *** '
        ierrstat = 1006
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

    CASE default
      
      yerrmsg = ' ERROR *** Wrong value for zitype_vertint *** '
      ierrstat = 1007
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      
    END SELECT

    !.. UB: Correct interpolation artifacts:
    !   (MAY BE COMMENTED IN BY SPECIALIZED USERS, BUT NOT OPERATIONNALY AT DWD!)
!!$    IF (TRIM(yzname) == 'RELHUM') THEN
!!$      results(istartpar:iendpar,j,:) = MAX (0.0_ireals, MIN(100.0_ireals,results(istartpar:iendpar,j,:)))
!!$    ELSEIF (TRIM(yzname) == 'QC') THEN
!!$      results(istartpar:iendpar,j,:) = MAX (0.0_ireals, results(istartpar:iendpar,j,:) )
!!$    END IF

  ENDDO ! loop over j

  ! In case of radar reflectivity, transform back to log space after interpolation:
  IF (yzname(1:izlen) == 'DBZ') THEN
    WHERE (results(istartpar:iendpar,jstartpar:jendpar,:) >= 1e-9_ireals)
      results(istartpar:iendpar,jstartpar:jendpar,:) = &
           10.0_ireals * LOG10(results(istartpar:iendpar,jstartpar:jendpar,:))
    ELSEWHERE
      results(istartpar:iendpar,jstartpar:jendpar,:) = -99.99
    END WHERE
  END IF

  !.. UB: Set data below the surface to missing values: 
  !   (MAY BE COMMENTED IN BY SPECIALIZED USERS, BUT NOT OPERATIONNALY AT DWD!)
!!$  SELECT CASE (yzname(1:izlen))
!!$  CASE ('DBZ')
!!$    DO i=istartpar,iendpar
!!$      DO j=jstartpar,jendpar
!!$        DO k=1,outblock%kepin
!!$          IF (outblock%plev(k) > ps(i,j,itlev)) results(i,j,k) = -99.99
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO
!!$  CASE default
!!$    DO i=istartpar,iendpar
!!$      DO j=jstartpar,jendpar
!!$        DO k=1,outblock%kepin
!!$          IF (outblock%plev(k) > ps(i,j,itlev)) results(i,j,k) = 0.0_ireals
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO
!!$  END SELECT

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE p_int

!==============================================================================
!+ Module procedure in src_output for the z-interpolation
!------------------------------------------------------------------------------

SUBROUTINE z_int (outblock, i1,i2,i3, nlist, itlev, results)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine interpolates a number of 3-d variables, given by namelist-
!   input to this routine, from model layers to fixed-height levels (z-levles)
!
! Method:
!   Colum wise interpolation with tension splines.
!   For the interpolation of FI and QV the logarithm of the pressure is
!   used.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
TYPE(pp_nl),              INTENT(IN) ::     &
  outblock       ! pointer to the namelist group

INTEGER (KIND=iintegers), INTENT(IN) ::     &
  i1,i2,i3,      & ! location of the variable to be processed in the LM
                 ! variable table
  nlist,       & ! location of the variable in the output list
  itlev          ! time level of the variables

! Array arguments with intent(out):
REAL (KIND=ireals)      , INTENT(OUT)::     &
  results(ie,je,outblock%kezin)

!------------------------------------------------------------------------------

! Local parameters
REAL(KIND=ireals),PARAMETER    :: gamma   = 5.5_ireals    ! tension factor

INTEGER (KIND=iintegers)       :: i, j, k, kint(ie)
INTEGER (KIND=iintegers)       :: ierrstat, ierr, izlen
INTEGER (KIND=iintegers)       :: nldim(ie,1)

CHARACTER (LEN=25)             :: yroutine
CHARACTER (LEN=80)             :: yerrmsg
CHARACTER (LEN=10)             :: yzname

! Local arrays:
REAL(KIND=ireals)              :: fmfl(ie,je,ke+1),  &
                                  zmfl(ie,   ke+1)
REAL(KIND=ireals)              :: fexp(ie,ke+4),  &
                                  zexp(ie,ke+4)

! Output from tautsp
REAL(KIND=ireals)              :: s_vec    (ie,(ke+4)*6), &
                                  break_vec(ie,(ke+4)*3), &
                                  coef_vec (ie,4,(ke+4)*3)

! Output from spline
REAL(KIND=ireals)              :: fzfls(outblock%kezin),  &
                                  zfls (outblock%kezin)

! Switch to choose between vertical interpolation methods:
! Local value; there is a global namelist parameter itype_vertint
! in each namelist of group GRIBOUT.
INTEGER(KIND=iintegers) :: zitype_vertint
!
! Switch to choose if values of, U, V, W, and T are extrapolated
! to the surface based on no-slip- and skin-condition
! or if a constant extrapolation is done.
! The former is not so good for cubic spline interpolation
! because of oszillations with unwanted results, so we 
! allow it only for linear interpolation.
LOGICAL :: zlsurfextrapol_noslip

!------------------------------------------------------------------------------
! Section 1: general and grib-io preparations
!------------------------------------------------------------------------------

yroutine = 'z_int'
yzname   = '          '
yzname(1:LEN_TRIM(outblock%yvarzl(nlist))) =                            &
         outblock%yvarzl(nlist)(1:LEN_TRIM(outblock%yvarzl(nlist)))
izlen = LEN_TRIM(yzname)

! The default value is taken from the namelist value, 
! which may be "hard" changed for single fields below
! (e.g., if you would not like to have splineinterpolation for fiels like QC ...)
zitype_vertint = outblock%itype_vertint

! Allow  no-slip- and skin-condition-extrapolation to the
! surface only for linear interpolation:
IF (zitype_vertint == 2) THEN
  zlsurfextrapol_noslip = .TRUE.
ELSE
  zlsurfextrapol_noslip = .FALSE.
END IF

!------------------------------------------------------------------------------
! Section 1: Set variable that has to be interpolated
!------------------------------------------------------------------------------

  ! calculation of u and v on masspoint
  ! -----------------------------------

  IF     (yzname(1:izlen) == 'U' ) THEN
    fmfl(2:ie,:,1:ke) = 0.5_ireals * (var(i1,i2,i3)%p4(2:ie,:,1:ke,itlev)    &
                                    + var(i1,i2,i3)%p4(1:ie-1,:,1:ke,itlev))
!CDIR COLLAPSE
    fmfl(1,:,1:ke)    = fmfl(2,:,1:ke)
    IF (lnosurffluxes_m .OR. .NOT.zlsurfextrapol_noslip) THEN
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = fmfl(:,:,ke  )
    ELSE
      ! For interpolation: Surface velocity = 0
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = 0.0
    END IF
  ELSEIF (yzname(1:izlen) == 'V' ) THEN
!CDIR COLLAPSE
    fmfl(:,2:je,1:ke) = 0.5_ireals * (var(i1,i2,i3)%p4(:,2:je,1:ke,itlev)    &
                                    + var(i1,i2,i3)%p4(:,1:je-1,1:ke,itlev))
    fmfl(:,1,1:ke)    = fmfl(:,2,1:ke)
    IF (lnosurffluxes_m .OR. .NOT.zlsurfextrapol_noslip) THEN
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = fmfl(:,:,ke  )
    ELSE
      ! For interpolation: Surface velocity = 0
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = 0.0
    END IF
  ELSEIF (yzname(1:izlen) == 'P' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = p0(:,:,1:ke) + pp(:,:,1:ke,itlev)
! For interpolation: Surface pressure at the ground!
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = ps(:,:,itlev)
  ELSEIF (yzname(1:izlen) == 'W' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = var(i1,i2,i3)%p4(:,:,1:ke,itlev)
    IF (lnosurffluxes_m .OR. .NOT.zlsurfextrapol_noslip) THEN
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = fmfl(:,:,ke)
    ELSE
      ! For interpolation: Surface vertical velocity = 0
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = 0.0
    END IF
  ELSEIF (yzname(1:izlen) == 'T' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = var(i1,i2,i3)%p4(:,:,1:ke,itlev)
    IF (lnosurffluxes_h .OR. .NOT.zlsurfextrapol_noslip) THEN
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = fmfl(:,:,ke)
    ELSE
      ! For interpolation: Surface temperature = interfacial temperature
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = t_g(:,:,itlev)
    END IF
  ELSEIF (yzname(1:izlen) == 'TKVM' .OR. yzname(1:izlen) == 'TKVH') THEN
!CDIR COLLAPSE
    fmfl(:,:,2:ke)    = var(i1,i2,i3)%p3(:,:,2:ke)
!CDIR COLLAPSE
    fmfl(:,:,1)       = fmfl(:,:,2)
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'Q_SEDIM') THEN
    IF (lprog_qi) THEN
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = qrs(:,:,:) - qi(:,:,:,itlev)
    ELSE
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = qrs(:,:,:)
    ENDIF
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'RELHUM') THEN
    CALL calrelhum(fmfl, t(:,:,:,itlev), pp(:,:,:,itlev), p0(:,:,:), &
                   qv(:,:,:,itlev),ie, je, ke, b1, b2w, b3, b4w, rdv, o_m_rdv)
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'OMEGA') THEN
    CALL calomega (fmfl, pp(:,:,:,nnew), pp(:,:,:,itlev), pptens(:,:,:), &
                   w(:,:,:,itlev), rho0(:,:,:), ie, je, ke, dt, g )
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'FI_ANAI' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = p_anai(:,:,1:ke) / rho(:,:,1:ke)
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'DBZ') THEN
    IF (itype_gscp == 3) THEN
      CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,         &
        klv850, my_cart_id, itype_gscp, t(:,:,:,itlev), qc(:,:,:,itlev)*rho, qr(:,:,:,itlev)*rho,   &
        qi(:,:,:,itlev)*rho, qs(:,:,:,itlev)*rho, z_radar = fmfl )
    ELSEIF (itype_gscp == 4) THEN
      CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice,         &
        klv850, my_cart_id, itype_gscp, t(:,:,:,itlev), qc(:,:,:,itlev)*rho, qr(:,:,:,itlev)*rho,   &
        qi(:,:,:,itlev)*rho, qs(:,:,:,itlev)*rho, q_grau  = qg(:,:,:,itlev)*rho,       &
        z_radar = fmfl )
    ENDIF
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
    ! Reflectivity interpolation should be done in linear space, not in logarithmic (as alternative, interpolation
    ! in the space of rain rate -- Z^0.6666 --- would also be desireable):
    WHERE (fmfl(:,:,:) >= -90.0)
      fmfl(:,:,:) = 10.0_ireals ** (0.1_ireals * fmfl(:,:,:))
    ELSEWHERE
      fmfl(:,:,:) = 0.0_ireals
    END WHERE
  ELSE
    SELECT CASE( var(i1,i2,i3)%rank )
    CASE(4)
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = var(i1,i2,i3)%p4(:,:,1:ke,itlev)
    CASE(3)
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = var(i1,i2,i3)%p3(:,:,1:ke)
    END SELECT
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke  )
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Interpolation
!------------------------------------------------------------------------------

  ! slicewise interpolation
  ! -----------------------

  DO j=jstartpar,jendpar

    ! Set the height of the modell levels (zmfl) and the height of
    ! the constant z-levels to interpolate to (zfls).

    SELECT CASE(var(i1,i2,i3)%levtyp)
    CASE(110)
      zmfl(:,1:ke) = 0.5_ireals * ( hhl(:,j,1:ke) + hhl(:,j,2:ke+1) )
      zmfl(:,ke+1) = hsurf(:,j)
      zfls(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)
    CASE(109)
      zmfl(:,1:ke+1) = hhl(:,j,1:ke+1)
      zfls(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)
    CASE DEFAULT
      yerrmsg = 'wrong leveltyp of input field'
      ierrstat = 2004
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    END SELECT

    ! To handle cases where the lowest z-level is below the orography, an
    ! additional data point at z = - 100mf  will be created, assuming a
    ! constant profile.
    ! The interpolated variable value will be set to undefined.

    ! vertical interpolation colum wise
    ! ---------------------------------

    DO i=istartpar,iendpar

      ! copy slice variable to colum variable (in reverse order, because
      ! data points must strictly increase in tautsp) and add an extra
      ! data/function point in any case

      kint (i) = ke + 2
!US   zexp (i,1) = - 100.0_ireals
      zexp (i,1) = MIN (-100.0_ireals, zmfl(i,ke+1) - 20.0_ireals)
      fexp (i,1) = fmfl(i,j,ke+1)
      nldim(i,1) = (ke+4)*3
    ENDDO

    DO k = 1, ke1
      DO i=istartpar,iendpar
        zexp(i,k+1)  = zmfl(i,  ke+2-k)
        fexp(i,k+1)  = fmfl(i,j,ke+2-k)
      ENDDO

    ENDDO ! loop over i


  ! Spline interpolation over total domain or slicewise
  ! ---------------------------------------------------

    SELECT CASE (zitype_vertint)

    CASE (1)
       ! Interpolation durch kubische TauT-Splines (Tension Splines):

      CALL tautsp2D (zexp(:,:), fexp(:,:), kint(:), ie, istartpar,     &
           iendpar, ke+4, gamma, s_vec, break_vec, coef_vec, nldim, ierr)
      IF (ierr == 2) THEN
        yerrmsg = 'wrong input in tautsp'
        ierrstat = 1004
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF
      IF (ierr /= 0) THEN
        yerrmsg = ' ERROR *** Error in tautsp *** '
        ierrstat = 1005
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

      CALL spline2D(break_vec,coef_vec,j,nldim,zfls,outblock%kezin,    &
           results,outblock%yvarzl(nlist))

    CASE (2)

      ! Linear Interpolation with respect to height:
       ierr = 0
       CALL lininterp2D_xinter1D_vec(zexp(istartpar:iendpar,1:ke+2), &
            fexp(istartpar:iendpar,1:ke+2), iendpar-istartpar+1, ke+2, &
            zfls(:), results(istartpar:iendpar,j,:), outblock%kezin, ierr)

       IF (ierr /= 0) THEN
         yerrmsg = ' ERROR *** Error in lininterp2D_xinter1D_vec *** '
         ierrstat = 1006
         CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
       ENDIF

     CASE default
       
       yerrmsg = ' ERROR *** Wrong value for zitype_vertint *** '
       ierrstat = 1007
       CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
       
     END SELECT
    
     !.. UB: Correct interpolation artifacts:
     !   (MAY BE COMMENTED IN BY SPECIALIZED USERS, BUT NOT OPERATIONNALY AT DWD!)
!!$     IF (TRIM(yzname) == 'RELHUM') THEN
!!$       results(:,j,:) = MAX (0.0_ireals, MIN(100.0_ireals,results(:,j,:)))
!!$     ELSEIF (TRIM(yzname) == 'QC') THEN
!!$       results(:,j,:) = MAX (0.0_ireals, results(:,j,:) )
!!$     END IF

  ENDDO ! loop over j

  ! In case of radar reflectivity, transform back to log space after interpolation:
  IF (yzname(1:izlen) == 'DBZ') THEN
    WHERE (results(istartpar:iendpar,jstartpar:jendpar,:) >= 1e-9_ireals)
      results(istartpar:iendpar,jstartpar:jendpar,:) = &
           10.0_ireals * LOG10(results(istartpar:iendpar,jstartpar:jendpar,:))
    ELSEWHERE
      results(istartpar:iendpar,jstartpar:jendpar,:) = -99.99
    END WHERE
  END IF

  !.. UB: Set data below the surface to missing values:
  !   (MAY BE COMMENTED IN BY SPECIALIZED USERS, BUT NOT OPERATIONNALY AT DWD!)
!!$  SELECT CASE (yzname(1:izlen))
!!$  CASE ('DBZ')
!!$    DO i=istartpar,iendpar
!!$      DO j=jstartpar,jendpar
!!$        DO k=1,outblock%kezin
!!$          IF (outblock%zlev(k) < hsurf(i,j)) results(i,j,k) = -99.99
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO
!!$  CASE default
!!$    DO i=istartpar,iendpar
!!$      DO j=jstartpar,jendpar
!!$        DO k=1,outblock%kezin
!!$          IF (outblock%zlev(k) < hsurf(i,j)) results(i,j,k) = 0.0_ireals
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO    
!!$  END SELECT


!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE z_int

!==============================================================================
!+ calculates the function values for given coefficients
!------------------------------------------------------------------------------

SUBROUTINE spline (break, coef, nldim, sptau, spgtau, ng)

!------------------------------------------------------------------------------
!
! Description:
!  spline calculates the values of the interpolation function for given
!  interpolation coefficients and arguments
!
! Method:
!  Spline interpolation.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  nldim, ng    ! Dimensions of the variables
!
! Array arguments with intent(in):
REAL (KIND=ireals),       INTENT(IN)    ::  &
  break(nldim),      & ! arguments for which the function has to be calculated
  coef(4,nldim),     & ! coefficients of the interpolation function
  sptau(ng)            !

!
! Array arguments with intent(out):
REAL (KIND=ireals),       INTENT(OUT)   ::  &
  spgtau(ng)

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers)    ::  i, j
REAL (KIND=ireals)          ::  dx, prod

!
!- End of header
!==============================================================================


! Calculate the product (sptau(i)-break(j)) * (sptau(i)-break(j+1))
! and look for the interval where sptau is
  DO i  = 1, ng
    DO j  = 1, nldim-1
      prod   = (sptau(i) - break(j)) * (sptau(i) - break(j+1))

      ! calculate the splines
      IF( prod <= 0.0_ireals ) THEN
        dx = sptau(i) - break(j)
        spgtau(i) = coef(1,j) + dx * (coef(2,j) + dx*0.5*(coef(3,j)     &
                                    + dx/3.0* coef(4,j)))
        EXIT       
      ENDIF
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE SPLINE

!==============================================================================

!option! -pvctl _on_adb
SUBROUTINE spline2D (break_vec, coef_vec, jline, nldim, sptau, ng,    &
                     results, yname)

!------------------------------------------------------------------------------
!
! Description:
!  spline calculates the values of the interpolation function for given
!  interpolation coefficients and arguments
!
! Method:
!  Spline interpolation.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  ng, jline            ! Dimensions of variables and line to be processed

INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  nldim(ie,1)

!
! Array arguments with intent(in):
REAL (KIND=ireals),       INTENT(IN)    ::  &
  break_vec(ie,*), & ! arguments for which the function is calculated
  coef_vec (ie,4,*),&! coefficients of the interpolation function
  sptau(ng)            !

CHARACTER (LEN=10),       INTENT(IN)    :: yname

! Array arguments with intent(out):
REAL(KIND=ireals),        INTENT(OUT)   :: results(ie,je,ng)

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers)    :: i, j, k, maxng, nzl(ie), nzcount, nzelements
REAL (KIND=ireals)          :: dx, prod
LOGICAL                     :: ldone(ie), levelready

!
!- End of header
!==============================================================================

  j = jline
  nzl(:) = 1
  nzelements = iendpar - istartpar + 1

  maxng = MAXVAL (nldim(istartpar:iendpar,1))

  DO k = 1, maxng  !ke+4*3
    nzcount = 0
    ldone(:) = .FALSE.

    DO WHILE (nzcount < nzelements)
      DO i = istartpar, iendpar

        IF (.NOT. ldone(i)) THEN
          ! sptau has to be between break_vec(i,k) and break_vec(i,k+1)

          levelready=(nzl(i) > ng)
          IF (.NOT. levelready) THEN
            IF (sptau(nzl(i)) > break_vec(i,k+1)) THEN
              levelready=.TRUE.
            END IF
          END IF

          IF (levelready) THEN
            ! the level nzl(i) is ready
            ldone(i) = .TRUE.
            nzcount    = nzcount + 1
          ELSE
            dx = sptau(nzl(i)) - break_vec(i,k)
            results(i,j,nzl(i)) = coef_vec(i,1,k) + dx     *(coef_vec(i,2,k) +   &
                  dx*0.5*(coef_vec(i,3,k) + dx/3.0 * coef_vec(i,4,k)))
            nzl(i) = nzl(i) + 1
          END IF

! if this if is splitted in 2 ifs (as below), it will not vectorize
!          IF (nzl(i) <= ng)  THEN
!            IF (sptau(nzl(i)) <= break_vec(i,k+1)) THEN
!
!              dx = sptau(nzl(i)) - break_vec(i,k)
!              results(i,j,nzl(i)) = coef_vec(i,1,k) + dx     *(coef_vec(i,2,k) +   &
!                    dx*0.5*(coef_vec(i,3,k) + dx/3.0 * coef_vec(i,4,k)))
!              nzl(i) = nzl(i) + 1
!            ELSE
!              ! the level nzl(i) is ready
!              ldone(i) = .TRUE.
!              nzcount    = nzcount + 1
!            ENDIF
!          ELSE
!            ! the level nzl(i) is ready
!            ldone(i) = .TRUE.
!            nzcount    = nzcount + 1
!          ENDIF
        ENDIF

      ENDDO
    ENDDO

  ENDDO

  SELECT CASE (yname)
  CASE ('RELHUM    ')
    DO i = istartpar, iendpar
      results(i,j,1:ng) = MAX (0.0_ireals, MIN(100.0_ireals,results(i,j,1:ng)))
    ENDDO
  CASE ('QV        ','QC        ','QI        ',    &
        'QR        ','QS        ','QG        ')
    DO i = istartpar, iendpar
      results(i,j,1:ng) = MAX (0.0_ireals, results(i,j,1:ng) )
    ENDDO
  END SELECT

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE spline2D

#ifdef NETCDF
!==============================================================================
!+ Module procedure in "src_output" for writing a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE write_nc_gdefs (ncid, outblock, icomm, npes, yextension,      &
                           yerrmsg, ierror)
 
!------------------------------------------------------------------------------
!
! Description:
!   This routine initializes global definitions for a NetCDF output.
!   It writes the latitude and longitudes values and the vertical coordinates.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
  TYPE(pp_nl),              INTENT(IN)     ::    &
    outblock         ! pointer to the namelist group

! Scalar arguments with intent(in):
  INTEGER (KIND=iintegers),   INTENT(IN)   ::    &
    ncid            ! NetCDF file IDr

  INTEGER (KIND=iintegers),   INTENT(IN)   ::    &
   icomm,    & ! MPI communicator
   npes        ! number of PEs

  CHARACTER (LEN=1),          INTENT(IN)   ::    &
    yextension    ! indicates model variables (''), p-('p') or z-levels ('z')

  CHARACTER (LEN=*),          INTENT(OUT)  ::    &
    yerrmsg     ! string for error messages

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers),   INTENT(OUT)  ::    &
    ierror          ! error index

!------------------------------------------------------------------------------
!
! Local scalars:

INTEGER (KIND=iintegers), PARAMETER  ::          &
  nbnds=2,   & !
  ntime=1,   & !
  nheight=1, & !
  nsynmsg=32

INTEGER (KIND=iintegers)  :: i, j, timestep, isect

INTEGER (KIND=iintegers)  :: &
    jgridVarID,     & ! NetCDF ID for rotated_pole
    jlonVarID,      & ! NetCDF ID for longitude
    jlatVarID,      & ! NetCDF ID for latitude
    jrlonVarID,     & ! NetCDF ID for rotated longitude
    jrlatVarID,     & ! NetCDF ID for rotated latitude
    jslonuVarID,    & ! NetCDF ID for U-direction shifted longitude
    jslatuVarID,    & ! NetCDF ID for U-direction shifted latititude
    jslonvVarID,    & ! NetCDF ID for V-direction shifted longitude
    jslatvVarID,    & ! NetCDF ID for V-direction shifted latitude
    jsrlonVarID,    & ! NetCDF ID for shifted rotated longitude
    jsrlatVarID,    & ! NetCDF ID for shifted rotated latitude
    jvcVarID,       & ! NetCDF ID for the vertical component
    jsectVarID,     & ! NetCDF ID for nhori (number of sectors for the horizon)
    jh2VarID,       & ! NetCDF ID for the 2m height
    jh10VarID,      & ! NetCDF ID for the 10m height
    jhtoaVarID,     & ! NetCDF ID for the TOA height
    jwbt13VarID,    & ! NetCDF ID for the 1.3 Celsius wet bulb temperature
    jsynmsgVarID,   & ! NetCDF ID for the synthetic satellite data (only: SYNMSG)
    jmsgchanVarID,  & ! NetCDF ID for the synthetic satellite data (grouped by channel)
    jsoilVarID,     & ! NetCDF ID for the multi soil layer component
    jsoilbdsID,     & ! NetCDF ID for the multi soil layer bounds
    jtimeID,        & ! NetCDF ID for the time
    jtbdsID           ! NetCDF ID for the time bounds
    
CHARACTER (LEN=40) :: ydate

CHARACTER (LEN= 1) :: &
  cvctype           ! character version of ivctype

CHARACTER (LEN=19) :: &
  creation_date     ! actual time when the data file is written

!Local arrays:
REAL (KIND=ireals) :: &
  zlongitude(outblock%ie_out_tot,outblock%je_out_tot),  & ! geographic longitudes
  zlatitude (outblock%ie_out_tot,outblock%je_out_tot),  & ! geographic latitudes
  zrotlon   (outblock%ie_out_tot),                      & ! rotated longitudes
  zrotlat   (outblock%je_out_tot),                      & ! rotated latitudes
  zslonu    (outblock%ie_out_tot,outblock%je_out_tot),  & ! U-direction shifted longitudes
  zslatu    (outblock%ie_out_tot,outblock%je_out_tot),  & ! U-direction shifted latitudes
  zslonv    (outblock%ie_out_tot,outblock%je_out_tot),  & ! V-direction shifted longitudes
  zslatv    (outblock%ie_out_tot,outblock%je_out_tot),  & ! V-direction shifted latitudes
  zsrlon    (outblock%ie_out_tot),                      & ! shifted rotated longitudes
  zsrlat    (outblock%je_out_tot)                         ! shifted rotated latitudes

REAL (KIND=ireals) :: &
  time(ntime), time_bnds(nbnds, ntime)
    
REAL (KIND=irealgrib) :: &
  zsoil_bnds(nbnds, ke_soil+1), &
  zczhls(0:ke_soil+1)
                      
REAL (KIND=ireals), ALLOCATABLE   :: &
  zvcoord(:)        ! vertical coordinate
    
REAL (KIND=ireals) :: &
     zheight_2m(nheight),     & !  2m height
     zheight_10m(nheight),    & ! 10m height
     zheight_toa(nheight),    & ! TOA height (actually model TOA)
     zwbtemp_13c(nheight),    & ! 1.3 Celsius wet bulb temperature
     zmsgchan_wave(nmsgchan)    ! synthetic satellite channel wavelenghts

INTEGER (KIND=iintegers) :: &
  zsect(nhori),   & !
  ztime_values(8)   ! holds the date and time information taken from 
                    ! internal subroutine date_and_time

INTEGER (KIND=iintegers) :: my_comm_id, implcode
 
!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF
IF (npes > 1) THEN
 ! Get id in communicator comm
  CALL MPI_COMM_RANK(icomm,my_comm_id,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_COMM_RANK failed'
    ierror  = 2
    RETURN
  ENDIF
ELSE
  my_comm_id = 0
ENDIF

! processor 0 does the job
IF (my_comm_id == 0) THEN

! some pre-settings
  ierror = 0
  timestep = ntstep
 
  time = 0._ireals
  time_bnds = 0._ireals
  

! determine the creation date
  CALL DATE_AND_TIME(values=ztime_values)
  WRITE (creation_date,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') &
            ztime_values(1),'-',ztime_values(2),'-',ztime_values(3),' ', &
            ztime_values(5),':',ztime_values(6),':',ztime_values(7)
  
! write the global attributes
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "title",    TRIM(yncglob_title))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "institution",               &
                                          TRIM(yncglob_institution))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "source",                    &
                                          TRIM(yncglob_source))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "project_id",                &
                                          TRIM(yncglob_project_id))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "experiment_id",             &
                                         TRIM(yncglob_experiment_id))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "realization",ncglob_realization)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "Conventions","CF-1.0")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL,                             &
                 "conventionsURL",                                   &
        "http://www.unidata.ucar.edu/packages/netcdf/conventions.html")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "contact",                  &
                                         TRIM(yncglob_contact))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL,"references",                &
                                         TRIM(yncglob_references))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "nccf_version","2.0")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, NF90_GLOBAL, "creation_date", creation_date)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

! define dimensions of time and time_bounds

  ierror=nf90_def_dim(ncid,"time", NF90_UNLIMITED, idims_id_out(5))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_var(ncid, "time", NF90_DOUBLE, (/ idims_id_out(5) /), jtimeID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jtimeID, "standard_name", "time")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jtimeID, "long_name", "time")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_dim(ncid,"bnds", nbnds, idims_id_out(6))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_def_var(ncid, "time_bnds", NF90_DOUBLE,       &
                (/ idims_id_out(6), idims_id_out(5) /), jtbdsID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jtbdsID, "long_name", "time bounds")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

! determine the actual forecast time

  IF (outblock%lanalysis) THEN
    ! When nudging is active the output fields are treated as analyses
    ! and the reference time is the actual forecast time
    ydate = 'seconds since '//yakdat1(1:4)//'-'//yakdat1(5:6)//'-'//     &
             yakdat1(7:8)//' '//yakdat1(9:10)//':00:00'
  ELSE
    ydate = 'seconds since '// ydate_ini(1:4)//'-'//ydate_ini(5:6)//'-'  &
            //ydate_ini(7:8)//' '//ydate_ini(9:10)//':00:00'
  ENDIF
  ierror = nf90_put_att (ncid, jtimeID, "units", ydate)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror = nf90_put_att (ncid, jtbdsID, "units", ydate)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF     (itype_calendar == 0) THEN
    ierror = nf90_put_att (ncid, jtimeID, "calendar", "proleptic_gregorian")
  ELSEIF (itype_calendar == 1) THEN
    ierror = nf90_put_att (ncid, jtimeID, "calendar", "360_day")
  ENDIF
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror = nf90_put_att (ncid, jtimeID, "bounds", "time_bnds")

! the time will be in seconds

  time(1)        = timestep*dt

!
!HJP Begin Change relevant for restart
!
! IF(lbdclim .AND. (outblock%nextstep /= 1)) THEN
!    time_bnds(1,1) = time(1) - &
!      (outblock%ngrib(outblock%nextstep) - outblock%ngrib(outblock%nextstep-1))*dt
! ELSE
!    time_bnds(1,1) = 0._ireals
! ENDIF
!

  time_bnds(1,1) = 0._ireals

  IF (lbdclim) THEN

     IF (outblock%nextstep /= 1) THEN
       time_bnds(1,1) = time(1) - &
           (outblock%ngrib(outblock%nextstep) - outblock%ngrib(outblock%nextstep-1))*dt
     ELSE
      IF (outblock%ngrib( outblock%nextstep) .GT. 0 ) THEN
         time_bnds(1,1) = time(1) -                             &
           (outblock%ngrib( outblock%nextstep+1 )-outblock%ngrib( outblock%nextstep ))*dt
      ENDIF
     ENDIF

  ENDIF
!
!HJP End Change relevant for restart
!
  time_bnds(2,1) = time(1)

! set the values of rotated North Pole in the grid_mapping attribute
  grid_mapping = 'rotated_pole'

  ierror=nf90_def_var(ncid, "rotated_pole", NF90_CHAR, jgridVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jgridVarID, "long_name",                 &
                            "coordinates of the rotated North Pole")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  
  ierror=nf90_put_att(ncid, jgridVarID, "grid_mapping_name",                &
                            "rotated_latitude_longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  
  ierror=nf90_put_att(ncid, jgridVarID, "grid_north_pole_latitude",         &
                                        REAL(pollat))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  
  ierror=nf90_put_att(ncid, jgridVarID, "grid_north_pole_longitude",        &
                                        REAL(pollon))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF (polgam /= 0._ireals) THEN

   ierror=nf90_put_att(ncid, jgridVarID, "north_pole_grid_longitude",      &
                                        REAL(polgam))

    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

   ENDIF

! determine rotated lat/lon values

 DO i = 1, outblock%ie_out_tot
    zrotlon(i) = startlon_tot + (i-1)*dlon
  ENDDO
  DO j = 1, outblock%je_out_tot
    zrotlat(j) = startlat_tot + (j-1)*dlat
  ENDDO

! define the attributes of longitude and latitude
  ierror=nf90_def_dim(ncid,"rlon",outblock%ie_out_tot, idims_id_out(1))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_dim(ncid,"rlat",outblock%je_out_tot, idims_id_out(2))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! rotated longitude
  ierror=nf90_def_var(ncid, "rlon", NF90_FLOAT, (/ idims_id_out(1) /),     &
                                                jrlonVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jrlonVarID, "standard_name", "grid_longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlonVarID, "long_name", "rotated longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlonVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! rotated latitude
  ierror=nf90_def_var(ncid, "rlat", NF90_FLOAT, (/ idims_id_out(2) /),      &
                                                jrlatVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "standard_name", "grid_latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "long_name", "rotated latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jrlatVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF (.NOT. outblock%luvmasspoint .AND. yextension /= 'c') THEN  
  zsrlon = zrotlon + dlon*0.5
  zsrlat = zrotlat + dlat*0.5

  ierror=nf90_def_dim(ncid,"srlon",outblock%ie_out_tot, idims_id_out(9))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_def_dim(ncid,"srlat",outblock%je_out_tot, idims_id_out(10))
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted rotated longitude
  ierror=nf90_def_var(ncid, "srlon", NF90_FLOAT, (/ idims_id_out(9) /),     &
                                                 jsrlonVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jsrlonVarID, "standard_name", "grid_longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlonVarID, "long_name",                       &
                                         "staggered rotated longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlonVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted rotated latitude
  ierror=nf90_def_var(ncid, "srlat", NF90_FLOAT, (/ idims_id_out(10) /),    &
                                                 jsrlatVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlatVarID, "standard_name", "grid_latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlatVarID, "long_name",                       &
                                         "staggered rotated latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jsrlatVarID, "units", "degrees")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF


  DO j = 1, outblock%je_out_tot
    DO i = 1, outblock%ie_out_tot
      zslonu(i,j) = rlarot2rla (zrotlat(j), zsrlon(i), pollat, pollon, polgam)
      zslatu(i,j) = phirot2phi (zrotlat(j), zsrlon(i), pollat, pollon, polgam)
      zslonv(i,j) = rlarot2rla (zsrlat(j), zrotlon(i), pollat, pollon, polgam)
      zslatv(i,j) = phirot2phi (zsrlat(j), zrotlon(i), pollat, pollon, polgam)
    ENDDO
  ENDDO

  ENDIF  ! luvmasspoint

  DO j = 1, outblock%je_out_tot
    DO i = 1, outblock%ie_out_tot
      zlongitude(i,j) = rlarot2rla (zrotlat(j), zrotlon(i), pollat, pollon, polgam)
      zlatitude(i,j)  = phirot2phi (zrotlat(j), zrotlon(i), pollat, pollon, polgam)
    ENDDO
  ENDDO

! set the values for rotated and true geographic longitudes and latitudes

! geographic longitude
  ierror=nf90_def_var(ncid, "lon", NF90_FLOAT, (/ idims_id_out(1),        &
                                               idims_id_out(2) /), jlonVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jlonVarID, "standard_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlonVarID, "long_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlonVarID, "units", "degrees_east")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! geographic latitude
  ierror=nf90_def_var(ncid, "lat", NF90_FLOAT, (/ idims_id_out(1),         &
                                               idims_id_out(2) /), jlatVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlatVarID, "standard_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlatVarID, "long_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jlatVarID, "units", "degrees_north")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF


  IF (.NOT. outblock%luvmasspoint .AND. yextension /= 'c') THEN  

 ! shifted geographic longitude or U wind component
  ierror=nf90_def_var(ncid, "slonu", NF90_FLOAT,                         &
                     (/ idims_id_out(9), idims_id_out(2) /), jslonuVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jslonuVarID, "standard_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonuVarID, "long_name",                    &
                                         "staggered U-wind longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonuVarID, "units", "degrees_east")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted geographic latitude in U wind component
  ierror=nf90_def_var(ncid, "slatu", NF90_FLOAT,                         &
                     (/ idims_id_out(9), idims_id_out(2) /), jslatuVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatuVarID, "standard_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatuVarID, "long_name",                    &
                                         "staggered U-wind latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatuVarID, "units", "degrees_north")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted geographic longitude or V wind component
  ierror=nf90_def_var(ncid, "slonv", NF90_FLOAT,                           &
                     (/ idims_id_out(1), idims_id_out(10) /), jslonvVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_att(ncid, jslonvVarID, "standard_name", "longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonvVarID, "long_name",                      &
                                         "staggered V-wind longitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslonvVarID, "units", "degrees_east")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! shifted geographic latitude in V wind component
  ierror=nf90_def_var(ncid, "slatv", NF90_FLOAT,                           &
                     (/ idims_id_out(1), idims_id_out(10) /), jslatvVarID)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatvVarID, "standard_name", "latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatvVarID, "long_name",                      &
                                         "staggered V-wind latitude")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_att(ncid, jslatvVarID, "units", "degrees_north")
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ENDIF ! luvmasspoint

! take care of different height co-ordinates and define the vertical 
! axis accordingly
  IF (yextension == 'p') THEN   ! pressure co-ordinate

    ALLOCATE (zvcoord(outblock%kepin))
    zvcoord(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
    
    ierror=nf90_def_dim(ncid,"pressure",outblock%kepin,idims_id_out(3))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_def_var(ncid, "pressure", NF90_FLOAT,                    &
                       (/ idims_id_out(3) /), jvcVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "long_name", "pressure")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "units", "Pa")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "positive", "down")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

  ELSEIF (yextension == 'z') THEN  ! geometric height co-ordinate
  
    ALLOCATE (zvcoord(outblock%kezin))
    zvcoord(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)

    ierror=nf90_def_dim(ncid,"altitude",outblock%kezin,idims_id_out(3))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_def_var(ncid, "altitude", NF90_FLOAT,                    &
                       (/ idims_id_out(3) /), jvcVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                                        "height above mean sea level")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "units", "m")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jvcVarID, "positive", "up")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

!_br 25.11.08
  ELSE ! model coordinate

      ALLOCATE (zvcoord(ke1))
      zvcoord(1:ke1) = vcoord(1:ke1)
      
      ierror=nf90_def_dim(ncid,"level",ke, idims_id_out(3))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_def_dim(ncid,"level1",ke1, idims_id_out(4))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_def_var(ncid, "vcoord", NF90_FLOAT,                     &
                         (/ idims_id_out(4) /), jvcVarID)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      IF (ivctype == 1) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                                           "Pressure based hybrid coordinate")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF (ivctype == 2) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                                    "Height-based hybrid Gal-Chen coordinate")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF (ivctype == 3) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "long_name",                     &
                                                           "SLEVE coordinate")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF

      ierror=nf90_put_att(ncid, jvcVarID, "units", "Pa")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "ivctype", ivctype)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "irefatm", irefatm)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "p0sl", p0sl)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "t0sl", t0sl)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "dt0lp", dt0lp)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror=nf90_put_att(ncid, jvcVarID, "vcflat", vcflat)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF

      IF (ivctype == 3) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "svc1", svc1)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jvcVarID, "svc2", svc2)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jvcVarID, "nfltvc", nfltvc)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF
      IF (irefatm == 2) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "delta_t", delta_t)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jvcVarID, "h_scal", h_scal)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF (irefatm == 3) THEN
        ierror=nf90_put_att(ncid, jvcVarID, "bvref", bvref)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF

      IF (lradtopo) THEN
         ierror=nf90_def_dim(ncid,"nhori", nhori, idims_id_out(15))
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
         ierror=nf90_def_var(ncid, "nhori", NF90_FLOAT,                 &
              (/ idims_id_out(15) /), jsectVarID)
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
         ierror=nf90_put_att(ncid, jsectVarID, "standard_name", "sectors")
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
         ierror=nf90_put_att(ncid, jsectVarID, "long_name", "sectors of horizontal angles around grid cells")
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
         ierror=nf90_put_att(ncid, jsectVarID, "units", "sector number starting from the North clockwise")
         IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
         ENDIF
      ENDIF

      IF (yextension /= 'c') THEN
        ! height_2m, height_10m, height_toa, and wbtemp_13c
        ! not needed in files holding the constant fields

        zheight_2m  = 2.
        zheight_10m = 10.
        zheight_toa = hhl(1,1,1)
        zwbtemp_13c = 1.3
        ierror=nf90_def_dim(ncid,"height_2m", nheight, idims_id_out(11))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_dim(ncid,"height_10m", nheight, idims_id_out(12))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_dim(ncid,"height_toa", nheight, idims_id_out(13))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_dim(ncid,"wbt_13c", nheight, idims_id_out(14))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_var(ncid, "height_2m", NF90_FLOAT,                 &
                           (/ idims_id_out(11) /), jh2VarID)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh2VarID, "standard_name", "height")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh2VarID, "long_name",                   &
                                            "height above the surface")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh2VarID, "units", "m")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh2VarID, "positive", "up")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_var(ncid, "height_10m", NF90_FLOAT,                &
                           (/ idims_id_out(12) /), jh10VarID)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh10VarID, "standard_name", "height")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh10VarID, "long_name",                  &
                                             "height above the surface")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh10VarID, "units", "m")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jh10VarID, "positive", "up")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_var(ncid, "height_toa", NF90_FLOAT, (/ idims_id_out(13) /), jhtoaVarID)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jhtoaVarID, "standard_name", "height")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jhtoaVarID, "long_name", "height of top of model")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jhtoaVarID, "units", "m")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jhtoaVarID, "positive", "up")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_def_var(ncid, "wbt_13c", NF90_FLOAT, (/ idims_id_out(14) /), jwbt13VarID)
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jwbt13VarID, "long_name", "wet bulb temperature")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        ierror=nf90_put_att(ncid, jwbt13VarID, "units", "Celsius")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF

      ENDIF
!_br 25.11.08 end

  ENDIF

! define fields of variable multi layer soil model
  IF (lmulti_layer .AND. &
     yextension /= 'c' .AND. yextension /= 'p' .AND. yextension /= 'z') THEN
  

    ierror=nf90_def_dim(ncid, "soil", ke_soil, idims_id_out(7))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_def_dim(ncid, "soil1", ke_soil+1, idims_id_out(8))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF


    ierror=nf90_def_var(ncid, "soil1", NF90_FLOAT, (/ idims_id_out(8) /),  &
                                                   jsoilVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
 
     ierror=nf90_put_att(ncid, jsoilVarID, "standard_name", "depth")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilVarID, "long_name", "depth of soil layers")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilVarID, "units", "m")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilVarID, "positive", "down")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilVarID, "bounds", "soil1_bnds")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_def_var(ncid, "soil1_bnds", NF90_FLOAT,                    &
                       (/ idims_id_out(6), idims_id_out(8) /), jsoilbdsID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jsoilbdsID, "long_name",                     &
                                          "boundaries of soil layers")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

  ENDIF  ! end defining fields of variable multi layer soil model

  IF (luse_rttov .AND. yextension == 's') THEN
    ierror=nf90_def_dim(ncid, "nsynmsg", nsynmsg, idims_id_out(15))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_def_dim(ncid, "msgchan", nmsgchan, idims_id_out(16))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    zmsgchan_wave = (/ 1.6, 6.2, 7.3, 3.8, 8.7, 10.8, 12.0, 9.7 /)
    ierror=nf90_def_var(ncid, "msgchan_wave", NF90_FLOAT, (/ idims_id_out(16) /),  &
                                                   jmsgchanVarID)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    ierror=nf90_put_att(ncid, jmsgchanVarID, "standard_name", "channel_wave_length")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jmsgchanVarID, "long_name", "channel wave length")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_att(ncid, jmsgchanVarID, "units", "micrometer")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF

! End of definition mode
!!! no more attribute definitions beyond this line !!!

  ierror=nf90_enddef(ncid)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

! write the values of longitude, latitude and vertical axis to file
  ierror=nf90_put_var(ncid, jlonVarID, zlongitude)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jlatVarID, zlatitude)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jrlonVarID, zrotlon)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF
  ierror=nf90_put_var(ncid, jrlatVarID, zrotlat)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF (.NOT. outblock%luvmasspoint .AND. yextension /= 'c') THEN

    ierror=nf90_put_var(ncid, jsrlonVarID, zsrlon)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jsrlatVarID, zsrlat)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jslonuVarID, zslonu)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jslatuVarID, zslatu)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jslonvVarID, zslonv)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jslatvVarID, zslatv)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

  ENDIF


  ierror=nf90_put_var(ncid, jvcVarID, zvcoord)   
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  IF ( yextension /= 'c' .AND. yextension /= 'p' .AND. yextension /= 'z') THEN 
    ierror=nf90_put_var(ncid, jh2VarID, zheight_2m)   
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jh10VarID, zheight_10m)   
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jhtoaVarID, zheight_toa)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror=nf90_put_var(ncid, jwbt13VarID, zwbtemp_13c)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF

    IF (lmulti_layer .AND. &
     yextension /= 'c' .AND. yextension /= 'p' .AND. yextension /= 'z') THEN
      ierror=nf90_put_var(ncid, jsoilVarID, czmls(1:ke_soil+1))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      zsoil_bnds(1, 1) = 0.
      zsoil_bnds(1, 2:ke_soil+1) = REAL(czhls(1:ke_soil))
      zsoil_bnds(2, 1:ke_soil+1) = REAL(czhls(1:ke_soil+1))
      ierror=nf90_put_var(ncid, jsoilbdsID, zsoil_bnds)
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF
  ENDIF

  IF (lradtopo .AND. yextension /= 'p' .AND. yextension /= 'z') THEN
     DO isect=1,nhori
        zsect(isect)=isect
     ENDDO
     ierror=nf90_put_var(ncid, jsectVarID, zsect)
     IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
     ENDIF
  ENDIF

  ! for the synthetic satellite images
  IF (luse_rttov .AND. yextension == 's') THEN
    ierror=nf90_put_var(ncid, jmsgchanVarID, zmsgchan_wave)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF

! write the values of time and time_bnds to file
  ierror=nf90_put_var(ncid, jtimeID, time)   
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ierror=nf90_put_var(ncid, jtbdsID, time_bnds)   
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  DEALLOCATE (zvcoord)

ENDIF   ! end of IF statement for processor 0

IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE write_nc_gdefs

!==============================================================================
!+ Module procedure in "io_utilities" for writing a NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE write_nc_vdefs (ncid, numlist, ilist, var_id, luvmasspoint,      &
                           icomm, npes, yextension, yerrmsg, ierror)

!------------------------------------------------------------------------------
!
! Description:
!   write NetCDF attributes for each output parameter
!
!------------------------------------------------------------------------------

! Scalar arguments with intent(in):
  INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
    ncid             ! NetCDF file ID

  INTEGER (KIND=iintegers), INTENT(IN)     ::  &
    numlist          ! number of variables for output

  INTEGER  (KIND=iintegers), INTENT(IN)    ::  &
    icomm,    & ! MPI communicator
    npes        ! number of PEs

   LOGICAL                   , INTENT(IN) :: &
   luvmasspoint     ! interpolate horizontal winds to mass grid points

  CHARACTER (LEN=1),          INTENT(IN)   ::  &
    yextension    ! indicates model variables (''), p-('p') or z-levels ('z')

! Array arguments with intent(in):
  INTEGER (KIND=iintegers), INTENT(IN)     ::  &
    ilist(3,numlist) ! 

! Array arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)    ::  &
    var_id(nzmxid)      ! NetCDF-ID of each variable 

  CHARACTER (LEN= *)       , INTENT(OUT)   ::  &
    yerrmsg     ! string for error messages

! Scalar arguments with intent(out):
  INTEGER (KIND=iintegers), INTENT(OUT)     ::  &
    ierror          ! error index

!------------------------------------------------------------------------------

! Local scalars:
  INTEGER (KIND=iintegers)  :: i1,i2,i3, iid, jid, n, timestep
  CHARACTER (LEN= 60)       :: ymethod

  INTEGER (kind=iintegers)  :: my_comm_id, implcode

  REAL (kind=ireals)        :: zdtime

!
! Local arrays:
!
!- End of header
!==============================================================================

IF (npes > 1) THEN
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF
IF (npes > 1) THEN
 ! Get id in communicator comm
  CALL MPI_COMM_RANK(icomm,my_comm_id,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_COMM_RANK failed'
    ierror  = 2
    RETURN
  ENDIF
ELSE
  my_comm_id = 0
ENDIF

var_id(:) = 0

IF (my_comm_id == 0) THEN
  ierror=0
  timestep = ntstep

! Re-Enter definition mode
  ierror = NF90_redef (ncid)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

  ! loop over all variables that should be written and loop over all variables
  ! in the LM variable table until equal elements are found
  DO n = 1, numlist
    ! indices of field in variable table
    i1 = ilist(1,n)
    i2 = ilist(2,n)
    i3 = ilist(3,n)

    ! select netCDF dimension IDs
    IF( ( (TRIM(var(i1,i2,i3)%name) == 'U')            .OR.      &
          (TRIM(var(i1,i2,i3)%name) == 'AUMFL_S') )    .AND.     &
        (.NOT. luvmasspoint) ) THEN
       iid = 9
       jid = 2
    ELSEIF( ( (TRIM(var(i1,i2,i3)%name) == 'V')        .OR.      &
          (TRIM(var(i1,i2,i3)%name) == 'AVMFL_S') )    .AND.     &
        (.NOT. luvmasspoint) ) THEN
       iid = 1
       jid = 10
    ELSE
       iid = 1
       jid = 2
    ENDIF


! set the dimensions of the variable regarding to its rank and other criterea
    SELECT CASE (var(i1,i2,i3)%rank)

    CASE(4)

      IF ( (yextension == 'p') .OR. (yextension == 'z') .OR.               &
                            (UBOUND(var(i1,i2,i3)%p4,3) == ke) ) THEN
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(3), idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSEIF (UBOUND(var(i1,i2,i3)%p4,3) == ke1) THEN
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(4), idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSEIF (UBOUND(var(i1,i2,i3)%p4,3) == ke_soil+1) THEN
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(8), idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF

    CASE(3)

      IF (yextension == 'p' .OR. yextension == 'z') THEN
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(3), idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE
        IF (.NOT. ASSOCIATED(var(i1,i2,i3)%p3)) THEN
          IF (var(i1,i2,i3)%levtyp == 109) THEN
            ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,&
                            (/ idims_id_out(iid), idims_id_out(jid),         &
                               idims_id_out(4), idims_id_out(5) /), var_id(n))
           IF (ierror /= NF90_NOERR) THEN
              yerrmsg = TRIM(NF90_strerror(ierror))
              RETURN
            ENDIF
          ELSEIF (var(i1,i2,i3)%levtyp == 110) THEN
            ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,&
                            (/ idims_id_out(iid), idims_id_out(jid),         &
                               idims_id_out(3), idims_id_out(5) /), var_id(n))
            IF (ierror /= NF90_NOERR) THEN
              yerrmsg = TRIM(NF90_strerror(ierror))
              RETURN
            ENDIF
          ELSE
            ierror = 1
            WRITE(yerrmsg,'(A,I4,A)') 'levtyp = ',var(i1,i2,i3)%levtyp,     &
                                      ' not implemented'
          ENDIF
        ELSEIF (UBOUND(var(i1,i2,i3)%p3,3) <= 3) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSEIF (UBOUND(var(i1,i2,i3)%p3,3) == ke) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(3), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSEIF (UBOUND(var(i1,i2,i3)%p3,3) == ke1) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(4), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSEIF (UBOUND(var(i1,i2,i3)%p3,3) == nhori) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
               (/ idims_id_out(iid), idims_id_out(jid),        &
               idims_id_out(15), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
             yerrmsg = TRIM(NF90_strerror(ierror))
             RETURN
          ENDIF
        ELSE
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(3), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ENDIF
      ENDIF
    
      ! for the synthetic satellite images
      IF (var(i1,i2,i3)%levtyp == 222) THEN
        IF (TRIM(var(i1,i2,i3)%name) == 'SYNMSG') THEN
            ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,&
                            (/ idims_id_out(iid), idims_id_out(jid),         &
                               idims_id_out(15), idims_id_out(5) /), var_id(n))
            IF (ierror /= NF90_NOERR) THEN
              yerrmsg = TRIM(NF90_strerror(ierror))
              RETURN
            ENDIF
        ELSE
            ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,&
                            (/ idims_id_out(iid), idims_id_out(jid),         &
                               idims_id_out(16), idims_id_out(5) /), var_id(n))
            IF (ierror /= NF90_NOERR) THEN
              yerrmsg = TRIM(NF90_strerror(ierror))
              RETURN
            ENDIF
        ENDIF
      ENDIF

    CASE(2)

      IF (var(i1,i2,i3)%levtyp == 105) THEN
      
        ! variables on a specific altitude
        IF (var(i1,i2,i3)%levbot == 2) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(11), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSE IF (var(i1,i2,i3)%levbot == 10) THEN
          ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                            (/ idims_id_out(iid), idims_id_out(jid),        &
                               idims_id_out(12), idims_id_out(5) /), var_id(n))
          IF (ierror /= NF90_NOERR) THEN
            yerrmsg = TRIM(NF90_strerror(ierror))
            RETURN
          ENDIF
        ELSE
          ierror = 1
          yerrmsg = 'invalid "levbot" '
          RETURN
        ENDIF
       
      ELSE IF (var(i1,i2,i3)%levtyp == 8) THEN

        ! variables at TOA (actually top of model)
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                       (/ idims_id_out(iid), idims_id_out(jid),           &
                          idims_id_out(13),  idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF

      ELSE IF (TRIM(var(i1,i2,i3)%name) == 'SNOWLMT') THEN

        ! Snow limit has coordinate variable wet bulb temp.
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT, &
                       (/ idims_id_out(iid), idims_id_out(jid),           &
                          idims_id_out(14),  idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF

      ELSE
      
        ierror = nf90_def_var(ncid, TRIM(var(i1,i2,i3)%name), NF90_FLOAT,  &
                          (/ idims_id_out(iid), idims_id_out(jid),         &
                             idims_id_out(5) /), var_id(n))
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
        
      ENDIF
    
    END SELECT
   
! set the attributes of the variable
    IF (var(i1,i2,i3)%standard_name(1:1) /= '-') THEN
      ierror = nf90_put_att (ncid, var_id(n), "standard_name",               &
                                              TRIM(var(i1,i2,i3)%standard_name))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF
    ierror = nf90_put_att (ncid, var_id(n), "long_name",                   &
                                            TRIM(var(i1,i2,i3)%long_name))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror = nf90_put_att (ncid, var_id(n), "units", TRIM(var(i1,i2,i3)%units))
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror = nf90_put_att (ncid, var_id(n), "grid_mapping", grid_mapping)
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    ierror = nf90_put_att (ncid, var_id(n), "coordinates", "lon lat")
    IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
    ENDIF
    IF(  (TRIM(var(i1,i2,i3)%name) == 'HTOP_CON')     .OR.      &
         (TRIM(var(i1,i2,i3)%name) == 'HBAS_CON')     .OR.      &
         (TRIM(var(i1,i2,i3)%name) == 'HZEROCL')      .OR.      &
         (TRIM(var(i1,i2,i3)%name) == 'SNOWLMT')      .OR.      &
         (var(i1,i2,i3)%lsm /= ' ') ) THEN
      ierror = nf90_put_att (ncid, var_id(n), "_FillValue", undefncdf)
      IF (ierror /= NF90_NOERR) THEN
      yerrmsg = TRIM(NF90_strerror(ierror))
      RETURN
      ENDIF
    ENDIF

    IF  (TRIM(var(i1,i2,i3)%name) == 'HMO3') THEN
      ierror = nf90_put_att (ncid, var_id(n), "cell_methods", &
                          "mole_fraction_of_ozone_in_air: maximum")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

    IF (TRIM(var(i1,i2,i3)%name) == 'SOILTYP') THEN
      ierror = nf90_put_att (ncid, var_id(n), "flag_values", &
                            (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 0 /))
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
      ierror = nf90_put_att (ncid, var_id(n), "flag_meanings", &
         "ice rock sand sandy_loam loam clay_loam clay peat sea_water sea_ice")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

! set attribute for vertical coordinate
    IF (TRIM(var(i1,i2,i3)%name) == 'HHL') THEN
      ierror=nf90_put_att(ncid, var_id(n), "positive", "up")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

    IF( ((TRIM(var(i1,i2,i3)%name) == 'U')            .OR.      &
        (TRIM(var(i1,i2,i3)%name) == 'AUMFL_S'))      .AND.     &
        (.NOT. luvmasspoint) ) THEN
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "slonu slatu")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ELSEIF( ((TRIM(var(i1,i2,i3)%name) == 'V')        .OR.      &
        (TRIM(var(i1,i2,i3)%name) == 'AVMFL_S'))      .AND.     &
        (.NOT. luvmasspoint) ) THEN
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "slonv slatv")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ELSE
      ierror = nf90_put_att (ncid, var_id(n), "coordinates", "lon lat")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF

! Time Range Indicator determines the cell_methods attribute
! time will be in seconds
    SELECT CASE(var(i1,i2,i3)%ntri)
    CASE(2)
      IF(var(i1,i2,i3)%name == 'TMIN_2M') THEN 
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: minimum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF(var(i1,i2,i3)%name == 'TMAX_2M') THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: maximum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF(var(i1,i2,i3)%name == 'VMAX_10M') THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: maximum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF(var(i1,i2,i3)%name == 'VGUST_CON') THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: maximum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ELSE IF(var(i1,i2,i3)%name == 'VGUST_DYN') THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: maximum")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF
    CASE(3)
      ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: mean")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    CASE(4)
      ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "time: sum")
      IF (ierror /= NF90_NOERR) THEN
        yerrmsg = TRIM(NF90_strerror(ierror))
        RETURN
      ENDIF
    END SELECT 

! Set attribute for layer means of the soil layers
    IF (ASSOCIATED(var(i1,i2,i3)%p4)) THEN
      IF ( (var(i1,i2,i3)%rank == 4) .AND.                                 &
           (UBOUND(var(i1,i2,i3)%p4,3) == ke_soil+1) ) THEN
        ierror = nf90_put_att (ncid, var_id(n), "cell_methods", "soil1: mean")
        IF (ierror /= NF90_NOERR) THEN
          yerrmsg = TRIM(NF90_strerror(ierror))
          RETURN
        ENDIF
      ENDIF
    ENDIF

  ENDDO     ! End of loop over all variables


! End of definition mode
!!! no more attribute definitions beyond this line !!!

  ierror=nf90_enddef(ncid)
  IF (ierror /= NF90_NOERR) THEN
    yerrmsg = TRIM(NF90_strerror(ierror))
    RETURN
  ENDIF

ENDIF

IF (npes > 1) THEN
  CALL distribute_values  (var_id, nzmxid, 0, imp_integers,  icomm, implcode)
  CALL MPI_BARRIER(icomm,implcode)
  IF ( implcode /= 0 ) THEN
    yerrmsg = 'MPI_BARRIER failed'
    ierror  = 1
    RETURN
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine 
!------------------------------------------------------------------------------

END SUBROUTINE write_nc_vdefs
#endif

!==============================================================================
!+ Module procedure in src_output for the horizontal smoothing
!------------------------------------------------------------------------------

SUBROUTINE smooth_pmsl ( pmsl, hsurf, ie, je )

!------------------------------------------------------------------------------
! special smoothing of pmsl in mountainous terrain
!------------------------------------------------------------------------------

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)    :: ie, je
REAL    (KIND=ireals),    INTENT (IN)    :: hsurf(ie,je)

REAL    (KIND=irealgrib), INTENT (INOUT) :: pmsl(ie,je)

! Local Variables
REAL    (KIND=ireals)                    :: zp, hsurf_max, wgt, pmsl_sm(ie,je)
INTEGER (KIND=iintegers)                 :: i, j, ii, jj, n, nmean

!------------------------------------------------------------------------------

! Begin subroutine smooth_pmsl

  PRINT *,' smoothing pmsl over mountainous terrain '

  pmsl_sm(:,:) = REAL (pmsl(:,:), ireals)
  hsurf_max = 750.0_ireals

  DO n = 1,2
    DO j = 1, je
      DO i = 1, ie
        IF (hsurf(i,j) >= hsurf_max) THEN
          IF (hsurf(i,j) >= 1000.0_ireals) THEN
            nmean = 10
          ELSE
            nmean = 5
          ENDIF
          pmsl_sm(i,j) = 0.0_ireals
          zp = 0.0_ireals
          DO jj = j-nmean,j+nmean
            DO ii = i-nmean,i+nmean
              IF ( (ii >= 1) .AND. (ii <= ie) .AND.     &
                   (jj >= 1) .AND. (jj <= je) ) THEN
                IF (hsurf(ii,jj) < hsurf_max) THEN
                  wgt = 2.0_ireals
                ELSE
                  wgt = 1.0_ireals
                ENDIF
                pmsl_sm(i,j) = pmsl_sm(i,j) + wgt * REAL (pmsl(ii,jj), ireals)
                zp = zp + wgt
              ENDIF
            ENDDO
          ENDDO
          pmsl_sm(i,j) = pmsl_sm(i,j) / zp
        ENDIF
      ENDDO
    ENDDO
    pmsl(:,:) = REAL (pmsl_sm(:,:), irealgrib)
  ENDDO

END SUBROUTINE smooth_pmsl

!==============================================================================
!+ Module procedure in src_output for the horizontal smoothing
!------------------------------------------------------------------------------

SUBROUTINE smooth_geopot ( geopot, hsurf, ie, je )

!------------------------------------------------------------------------------
! special smoothing of geopot in mountainous terrain
!------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    :: ie, je
REAL    (KIND=ireals),    INTENT (IN)    :: hsurf(ie,je)
REAL    (KIND=irealgrib), INTENT (INOUT) :: geopot(ie,je)

! Local Variables
REAL    (KIND=ireals)                    :: zp, hsurf_max, wgt, geopot_diff, &
                                            geopot_sm(ie,je)
INTEGER (KIND=iintegers)                 :: i, j, ii, jj, n, nmean

!------------------------------------------------------------------------------

! Begin subroutine smooth_geopot

  PRINT*,' smoothing of geopotential height over mountainous terrain '

  geopot_sm(:,:) = REAL (geopot(:,:), ireals)

  DO n = 1,2
    DO j = 1,je
      DO i = 1,ie
        geopot_diff = g * hsurf(i,j) - REAL(geopot(i,j), ireals)
        IF (geopot_diff > 1000.0_ireals) THEN
          nmean = 10
        ELSE
          nmean = 5
        ENDIF
        geopot_sm(i,j) = 0.0_ireals
        zp = 0.0_ireals
        DO jj = j-nmean,j+nmean
          DO ii = i-nmean,i+nmean
            IF ( (ii >= 1) .AND. (ii <= ie) .AND.     &
                 (jj >= 1) .AND. (jj <= je) ) THEN
              geopot_diff = g * hsurf(ii,jj) - REAL(geopot(ii,jj), ireals)
              wgt = 1.0_ireals
              geopot_sm(i,j) = geopot_sm(i,j) + wgt*REAL(geopot(ii,jj), ireals)
              zp = zp + wgt
            ENDIF
          ENDDO
        ENDDO
        geopot_sm(i,j) = geopot_sm(i,j) / zp
      ENDDO
    ENDDO
    geopot(:,:) = REAL (geopot_sm(:,:), irealgrib)
  ENDDO

END SUBROUTINE smooth_geopot

!==============================================================================
!==============================================================================

SUBROUTINE calc_sdi( sdi_1, sdi_2 )

!------------------------------------------------------------------------------
!
! Description:
!   calculation of the 2 supercell detection indices (SDI)
!
! Method:
!   defined in:
!   Wicker L, J. Kain, S. Weiss and D. Bright, A Brief Description of the
!          Supercell Detection Index, (available from
!   http://www.spc.noaa.gov/exper/Spring_2005/SDI-docs.pdf)
!
!------------------------------------------------------------------------------

REAL (kind=ireals), DIMENSION(1:ie,1:je), INTENT(OUT) :: sdi_1, sdi_2

INTEGER (kind=iintegers) :: k_center
INTEGER (kind=iintegers) :: d_idx_x, d_idx_y, d_idx_z
INTEGER (kind=iintegers) :: i, j, k
REAL    (kind=ireals)    :: dx, dx_aequ, dy

REAL (kind=ireals), ALLOCATABLE, DIMENSION(:,:) :: &
  w_mean,       & ! <w>          mean over box in 'height' k_center
  w2_mean,      & ! <w'w'>       mean over box in 'height' k_center
  zeta_mean,    & ! <zeta>       mean over box in 'height' k_center
  zeta2_mean,   & ! <zeta'zeta'> mean over box in 'height' k_center
  helic_mean      ! <w'zeta'>    mean over box in 'height' k_center

REAL (kind=ireals), ALLOCATABLE, DIMENSION(:,:,:) :: &
  zeta,            & ! vorticity
  w_s                ! vertical velocity at scalar position

CHARACTER (LEN=80)       :: yzerrmsg
INTEGER (KIND=iintegers) :: izerror

INTEGER (KIND=iintegers) :: kzdims(24)

REAL (KIND=ireals) :: helic_w_corr          ! Correlation coefficient
REAL (KIND=ireals) :: zeta_vert_mean(ie,je)
REAL (KIND=ireals) :: w_crit_SDI2   (ie,je) ! Criterion, if SDI2 = 0 or not

REAL (KIND=ireals) :: EPS = 1.0E-20_ireals

!------------------------------------------------------------------------------

  izerror = 0_iintegers

  ! definition of the integration box:
  !     [ i_center-d_idx_x .. i_center+d_idx_x ]
  !   * [ j_center-d_idx_y .. j_center+d_idx_y ]
  !   * [ k_center-d_idx_z .. k_center+d_idx_z ]
  d_idx_x  = 3
  d_idx_y  = 3
  d_idx_z  = 6

  k_center = 30


  ! Allocations:
  ALLOCATE( zeta       ( 1:ie, 1:je, 1:ke) )
  ALLOCATE( w_s        ( 1:ie, 1:je, 1:ke) )
  ALLOCATE( w_mean     ( 1:ie, 1:je ) )
  ALLOCATE( w2_mean    ( 1:ie, 1:je ) )
  ALLOCATE( zeta_mean  ( 1:ie, 1:je ) )
  ALLOCATE( zeta2_mean ( 1:ie, 1:je ) )
  ALLOCATE( helic_mean ( 1:ie, 1:je ) )

  ! to prevent errors at the boundaries, set some fields to 0:
  sdi_1(:,:)  = 0.0_ireals
  sdi_2(:,:)  = 0.0_ireals
  zeta(:,:,:) = 0.0_ireals
  w_s(:,:,:)  = 0.0_ireals

  ! consistency checks:

  IF ( ( d_idx_x > nboundlines ) .OR. ( d_idx_y > nboundlines ) ) THEN
    yzerrmsg="integration box is too big in horizontal direction!"
    ! if such a big value for d_idx_x or d_idx_y is really needed, then you must increase
    ! nboundlines
    CALL model_abort (my_cart_id, 100, yzerrmsg, 'calc_sdi')
  END IF

  IF ( ( k_center - d_idx_z < 2 ) .OR. ( k_center + d_idx_z > ke ) ) THEN
    yzerrmsg="integration box is too big in vertical direction!"
    CALL model_abort (my_cart_id, 100, yzerrmsg, 'calc_sdi')
  END IF

  ! calculate vorticity and vertical velocity at scalar points:

  dx_aequ = r_earth * (pi/180.0_ireals) * dlon
  dy      = r_earth * (pi/180.0_ireals) * dlat

  DO k = k_center - d_idx_z, k_center + d_idx_z
    DO j = jstart, jend
      dx = dx_aequ * crlat(j,1)
      DO i = istart, iend
        zeta(i,j,k ) = (   ( v(i+1,j,  k,nnew) + v(i+1,j-1,k,nnew) )         &
          &              - ( v(i-1,j,  k,nnew) + v(i-1,j-1,k,nnew) ) )       &
          &                        * 0.5_ireals / dx                         &
          &          - (   ( u(i,  j+1,k,nnew) + u(i-1,j+1,k,nnew) )         &
          &              - ( u(i,  j-1,k,nnew) + u(i-1,j-1,k,nnew) ) )       &
          &                        * 0.5_ireals / dy

        w_s(i,j,k) = 0.5_ireals * ( w(i,j,k,nnew) + w(i,j,k-1,nnew) )

      END DO
    END DO
  END DO

  kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
  CALL exchg_boundaries                                                &
    (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute, ie, je,     &
    kzdims, jstartpar, jendpar, nbl_exchg, nboundlines, my_cart_neigh, &
    lperi_x, lperi_y, l2dim, &
    20000+ntstep, .FALSE.   , ncomm_type, izerror, yzerrmsg,           &
    zeta(:,:,:), w_s(:,:,:) )

  ! (exchange of w_s would not be necessary, if it is calculated also at the boundary lines)

  ! --- calculate mean values over the integration box: ------

  CALL mean_over_box( w_s,  w_mean,    k_center, d_idx_x, d_idx_y, d_idx_z,   &
                      crlat, sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  CALL mean_over_box( zeta, zeta_mean, k_center, d_idx_x, d_idx_y, d_idx_z,   &
                      crlat, sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  ! (no exchange needed for w_mean, zeta_mean)

  ! --- calculate covariances over the integration box: -----------

  CALL mean_cov_over_box( w_s,  w_mean,    zeta, zeta_mean, helic_mean,   &
                          k_center, d_idx_x, d_idx_y, d_idx_z, crlat,     &
                          sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  CALL mean_cov_over_box( w_s,  w_mean,    w_s,  w_mean,    w2_mean,      &
                          k_center, d_idx_x, d_idx_y, d_idx_z, crlat,     &
                          sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  CALL mean_cov_over_box( zeta, zeta_mean, zeta, zeta_mean, zeta2_mean,   &
                          k_center, d_idx_x, d_idx_y, d_idx_z, crlat,     &
                          sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  ! calculate SDI_1, SDI_2:
  ! call to vectorized version of vert_avg
  CALL vert_avg( zeta_vert_mean,zeta, sqrtg_r_s, ie, je, ke, istart, iend,  &
                 jstart, jend, k_center, d_idx_z)

  ! The meaning of 'w>0' in Wicker et al. is not completely clear, I assume
  ! the following:
  CALL vert_avg( w_crit_SDI2, w_s, sqrtg_r_s, ie, je, ke, istart, iend,     &
                 jstart, jend, k_center, d_idx_z)

  DO j = jstart, jend
    dx = dx_aequ * crlat(j,1)
    DO i = istart, iend

      IF ( ( w2_mean(i,j) > EPS ) .AND. ( zeta2_mean(i,j) > EPS ) ) THEN

        helic_w_corr = helic_mean(i,j) / SQRT(w2_mean(i,j) * zeta2_mean(i,j))

        sdi_1(i,j) = helic_w_corr * zeta_vert_mean(i,j)

        IF ( w_crit_SDI2(i,j) > 0 ) THEN
          sdi_2(i,j) = helic_w_corr * ABS( zeta_vert_mean(i,j) )
        ELSE
          sdi_2(i,j) = 0.0
        END IF

      ELSE
        sdi_1(i,j) = 0.0
        sdi_2(i,j) = 0.0
      END IF

    END DO
  END DO

  DEALLOCATE( zeta, w_s, w_mean, zeta_mean, w2_mean, zeta2_mean, helic_mean )

END SUBROUTINE calc_sdi

!==============================================================================

!==============================================================================
END MODULE src_output
