
#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>
!!  !MODULE:  mo_nwp_phy_state\\
!!
!! Description:  Contains the data structures
!!  to store the physical model state and other auxiliary variables
!!  in order to run the NWP land physics.
!!  Constructors and destructors for these data structures as well as
!!  initialization of fields are also defined here.
!!  This module should be an analogon to 'mo_hydro_state.f90'

!!  TODO/To think about:
!     - should physics be called before or after dynamics?
!     - allocate fluxes at edges instead at the centers?
!     - horizontal/vertical tracer flux (reconstruct q'v_n' into q'u' and q'v') ?
!     - provide the "virt_inc" with meaning
!     - where to provide the lat/lon info for radiation?
!     - how to implement the echam-modules - rewriting them or "capsulate"?
!     - revision of fields if there are needed or tp be replaced
!     - fill the physics tendency construction/destruction subroutine
!     - later implement already calculated icon gradients for echam physics
!     - think about variables for flexible time steps
!!
!! @author Kristina Froehlich, DWD
!! @author Marco Giorgetta, MPI-M
!!
!!
!! @par Revision History
!! Initial  by Kristina Froehlich (2010-11-09)
!! Modification by Daniel Reinert, DWD (2012-04-03)
!! - encapsulated type definitions (mo_nwp_lnd_types)
!! Modifications by Dmitrii Mironov, DWD (2016-08-04)
!! - Changes related to the use of a rate equation 
!!   for the sea-ice albedo.
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nwp_lnd_state

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH, HINTP_TYPE_LONLAT_NNB, &
    &                                TLEV_NNOW_RCF, ALB_SI_MISSVAL
  USE mo_parallel_config,      ONLY: nproma
  USE mo_nwp_lnd_types,        ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_exception,            ONLY: message, finish
  USE mo_model_domain,         ONLY: t_patch, p_patch !(SBr, CHa: ad p_patch)
  USE mo_math_utilities,       ONLY: t_geographical_coordinates ! SBr, CHa
  USE mo_ext_data_state,       ONLY: ext_data ! SBr, CHa
  USE mo_grid_config,          ONLY: n_dom
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,       ONLY: nlev_soil, nlev_snow, ntiles_total, &
    &                                lmulti_snow, ntiles_water, lseaice, llake, &
    &                                itype_interception, l2lay_rho_snow
  USE mo_io_config,            ONLY: lnetcdf_flt64_output
  USE mo_gribout_config,       ONLY: gribout_config
  USE mo_linked_list,          ONLY: t_var_list
  USE mo_var_list,             ONLY: default_var_list_settings,  &
    &                                add_var, add_ref,           &
    &                                new_var_list,               &
    &                                delete_var_list
  USE mo_var_metadata_types,   ONLY: POST_OP_SCALE, CLASS_TILE, CLASS_TILE_LAND
  USE mo_var_metadata,         ONLY: create_hor_interp_metadata, &
    &                                groups, post_op
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var, t_grib2_int_key, OPERATOR(+)
  USE mo_cdi,                  ONLY: DATATYPE_PACK16, DATATYPE_PACK24, DATATYPE_FLT32, &
    &                                TSTEP_ACCUM, GRID_UNSTRUCTURED, DATATYPE_FLT64
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL,                     & 
    &                                GRID_CELL, ZA_SURFACE, ZA_SNOW,             &
    &                                ZA_SNOW_HALF, ZA_DEPTH_BELOW_LAND,          &
    &                                ZA_DEPTH_BELOW_LAND_P1,                     &
    &                                ZA_DEPTH_RUNOFF_S, ZA_DEPTH_RUNOFF_G,       &
    &                                ZA_SEDIMENT_BOTTOM_TW_HALF, ZA_LAKE_BOTTOM, &
    &                                ZA_LAKE_BOTTOM_HALF, ZA_MIX_LAYER
  ! SBr, CHa
  USE mo_les_config,           ONLY: les_config
  USE mo_sync,                 ONLY: global_max, global_min
  USE mo_mpi,                  ONLY: p_pe




  IMPLICIT NONE
  PRIVATE


  !public interface
  !
  ! subroutines
  PUBLIC :: construct_nwp_lnd_state
  PUBLIC :: destruct_nwp_lnd_state

  !
  !variables
  PUBLIC :: p_lnd_state  !> state vector (variable) for land scheme


! complete state vector
!
  TYPE(t_lnd_state), TARGET, ALLOCATABLE :: p_lnd_state(:) 

  CONTAINS


!-------------------------------------------------------------------------
!!            SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
!-------------------------------------------------------------------------
!
  !>
  !! Constructor for prognostic and diagnostic states.
  !!
  !! It calls constructors to
  !! single time level prognostic states, and diagnostic states.
  !! Initialization of all components with zero.
  !!
  !! @par Revision History
  !! Initial release by Kristina Froehlich (2010-11-09)
  !!
  SUBROUTINE construct_nwp_lnd_state(p_patch, p_lnd_state, n_timelevels)
!
    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch(n_dom) ! patch
    INTEGER, OPTIONAL, INTENT(IN)       :: n_timelevels   ! number of timelevels

    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(n_dom)
                                           ! nh state at different grid levels
    INTEGER :: ntl, &! local number of timelevels
      &        ist, &! status
      &        jg,  &! grid level counter
      &        jt,  &! time level counter
      &      nblks_c ! number of blocks of cells

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname, varname_prefix

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nwp_lnd_state:construct_nwp_lnd_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'land state construction started')

    DO jg = 1, n_dom


      IF(PRESENT(n_timelevels))THEN
        ntl = n_timelevels
      ELSE
        ntl = 1
      ENDIF
   
      !
      !determine size of arrays
      nblks_c = p_patch(jg)%nblks_c


      !
      !create state arrays and lists
      !
      ALLOCATE(p_lnd_state(jg)%prog_lnd(1:ntl), &
           &   p_lnd_state(jg)%lnd_prog_nwp_list(1:ntl),STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_nwp_lnd_state:construct_lnd_state', &
             'allocation of land prognostic state array failed')
      ENDIF

      ALLOCATE(p_lnd_state(jg)%prog_wtr(1:ntl), &
               p_lnd_state(jg)%wtr_prog_nwp_list(1:ntl),STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish ('mo_nwp_lnd_state:construct_lnd_state', &
           'allocation of water prognostic state array failed')
      ENDIF

      !
      !construct prognostic state
      !
      DO jt = 1, ntl

        WRITE(listname,'(a,i2.2,a,i2.2)') 'lnd_prog_of_domain_',jg, &
          &                               '_and_timelev_',jt

        !varname_prefix = 'lnd_prog_'
        varname_prefix = ''
        CALL  new_nwp_lnd_prog_list(jg, nblks_c, TRIM(listname),             &
          &     TRIM(varname_prefix), p_lnd_state(jg)%lnd_prog_nwp_list(jt), &
          &     p_lnd_state(jg)%prog_lnd(jt), jt)


        WRITE(listname,'(a,i2.2,a,i2.2)') 'wtr_prog_of_domain_',jg, &
          &                               '_and_timelev_',jt

        varname_prefix = ''
        CALL new_nwp_wtr_prog_list(jg, nblks_c, TRIM(listname),             &
          &     TRIM(varname_prefix), p_lnd_state(jg)%wtr_prog_nwp_list(jt), &
          &     p_lnd_state(jg)%prog_wtr(jt), jt)

      ENDDO

      !
      !construct diagnostic state
      !
      WRITE(listname,'(a,i2.2)') 'lnd_diag_of_domain_',jg

      !varname_prefix = 'lnd_diag_'
      varname_prefix = ''
      CALL new_nwp_lnd_diag_list(jg, nblks_c, TRIM(listname),        &
        &   TRIM(varname_prefix), p_lnd_state(jg)%lnd_diag_nwp_list, &
        &   p_lnd_state(jg)%diag_lnd)

    ENDDO !ndom

    CALL message (TRIM(routine), 'land state construction completed')

  END SUBROUTINE construct_nwp_lnd_state

!-------------------------------------------------------------------------
!
!
  !>
  !! Destructor for prognostic and diagnostic states.
  !!
  !! It calls destructors to
  !! single time level prognostic states, and diagnostic states.
  !!
  !! @par Revision History
  !! Initial release by Kristina Froehlich (2010-11-09)
  !!
  SUBROUTINE destruct_nwp_lnd_state(p_lnd_state)
   !
    TYPE(t_lnd_state),  INTENT(INOUT) :: & 
      &   p_lnd_state(n_dom)             ! land state at different grid levels

    INTEGER :: ntl, &! local number of timelevels
      &        ist, &! status
      &         jg, &! grid level counter
      &         jt   ! time level counter

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nwp_lnd_state:destruct_nwp_lnd_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of nwp land state started')

    DO jg = 1, n_dom

      ntl = SIZE(p_lnd_state(jg)%prog_lnd(:))
      IF(ntl==0)THEN
        CALL finish(TRIM(routine), 'prognostic array has no timelevels')
      ENDIF


      DO jt = 1, ntl
        ! delete prognostic state list elements
        CALL delete_var_list(p_lnd_state(jg)%lnd_prog_nwp_list(jt) )
      ENDDO

      IF (lseaice .OR. llake) THEN
        DO jt = 1, ntl
          ! delete prognostic state list elements
          CALL delete_var_list(p_lnd_state(jg)%wtr_prog_nwp_list(jt) )
        ENDDO
      ENDIF

      ! delete diagnostic state list elements
      CALL delete_var_list( p_lnd_state(jg)%lnd_diag_nwp_list )


      ! destruct state lists and arrays
      DEALLOCATE(p_lnd_state(jg)%prog_lnd,&
           &     p_lnd_state(jg)%lnd_prog_nwp_list, STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine),  &
             'deallocation of land prognostic state array failed')
      ENDIF

      IF (lseaice .OR. llake) THEN
        ! destruct state lists and arrays
        DEALLOCATE(p_lnd_state(jg)%prog_wtr,&
             &     p_lnd_state(jg)%wtr_prog_nwp_list, STAT=ist)
        IF(ist/=SUCCESS)THEN
          CALL finish (TRIM(routine),  &
               'deallocation of water prognostic state array failed')
        ENDIF
      ENDIF

    ENDDO

    CALL message (TRIM(routine), 'Destruction of nwp land state completed')

  END SUBROUTINE destruct_nwp_lnd_state

!-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of prognostic land state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by , MPI-M (2011-07-01)
  !!
  !!
  SUBROUTINE new_nwp_lnd_prog_list( p_jg, kblks, listname, vname_prefix, &
    &                               prog_list, p_prog_lnd, timelev )

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes
    INTEGER,INTENT(IN) ::  p_jg  !< patch id

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix
    CHARACTER(LEN=2)            :: csfc

    TYPE(t_var_list),INTENT(INOUT) :: prog_list
    TYPE(t_lnd_prog),INTENT(INOUT) :: p_prog_lnd

    INTEGER, INTENT(IN) :: timelev

    ! Local variables
    !
    TYPE(t_cf_var)    :: cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d_subs(3), shape3d_subsw(3)
    INTEGER :: shape4d_snow_subs(4)
    INTEGER :: ibits
    INTEGER :: jsfc          !< tile counter
    INTEGER :: datatype_flt

    CHARACTER(len=4) suffix

    ! SBr, CHa
    INTEGER :: jc, jb, jt, jk, ic
    TYPE(t_geographical_coordinates) :: point
    REAL(wp) :: minlat, maxlat, minlon, maxlon, &
                minlat_loc, maxlat_loc, mean_wso(1:nlev_soil), &
                coeff_A(1:nlev_soil), tot_area, wso_mean(1:nlev_soil), &
                wso_min(1:nlev_soil), &
                wso_max(1:nlev_soil), wso_spread, r(1:nlev_soil), l(1:nlev_soil), &
                wso_grad(1:nlev_soil), wso_grad_fun4(1:nlev_soil), l_fun4(1:nlev_soil),&
                wso_grad_fun5(1:nlev_soil), l_fun5(1:nlev_soil),&
                wso_grad_fun5cyc2(1:nlev_soil), l_fun5cyc2(1:nlev_soil)
    INTEGER :: wso_fun, wso_cyc, mproc_id
    REAL(wp), PARAMETER :: pi = 4._wp*atan(1._wp)
    CHARACTER(256) :: sval

!-----------------------------------------------------------------------

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d              = (/nproma,            kblks            /)
    shape3d_subs         = (/nproma,            kblks, ntiles_total /)
    shape3d_subsw        = (/nproma,            kblks, ntiles_total+ntiles_water /)
    shape4d_snow_subs    = (/nproma, nlev_snow, kblks, ntiles_total /)


    ! Suffix (mandatory for time level dependent variables)

    WRITE(suffix,'(".TL",i1)') timelev


    !------------------------------
    ! Ensure that all pointers have a defined association status.
    !------------------------------
    NULLIFY(p_prog_lnd%t_s_t, &
    &       p_prog_lnd%t_g, &
    &       p_prog_lnd%t_g_t, &
    &       p_prog_lnd%w_i_t, &
    &       p_prog_lnd%w_p_t, &
    &       p_prog_lnd%w_s_t, &
    &       p_prog_lnd%t_so_t, &
    &       p_prog_lnd%w_so_t, &
    &       p_prog_lnd%w_so_ice_t, &
    &       p_prog_lnd%t_snow_t, &
    &       p_prog_lnd%w_snow_t, &
    &       p_prog_lnd%rho_snow_t, &
    &       p_prog_lnd%t_snow_mult_t, &
    &       p_prog_lnd%wtot_snow_t, &
    &       p_prog_lnd%wliq_snow_t, &
    &       p_prog_lnd%rho_snow_mult_t, &
    &       p_prog_lnd%dzh_snow_t)


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( prog_list, TRIM(listname), patch_id=p_jg )
    CALL default_var_list_settings( prog_list,                 &
                                  & lrestart=.TRUE.  )

    !------------------------------

    ! & p_prog_lnd%t_g(nproma,nblks_c), STAT = ist)
    cf_desc    = t_cf_var('t_g', 'K', 'weighted surface temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_g'//suffix, p_prog_lnd%t_g,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & ldims=shape2d,                                                      &
         & tlev_source=TLEV_NNOW_RCF,                      &! for output take field from nnow_rcf slice
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in",     &
         &                 "mode_iau_fg_in","mode_iau_old_fg_in",              &
         &                 "mode_combined_in","mode_cosmode_in") ) 


    ! & p_prog_lnd%t_g_t(nproma,nblks_c,ntiles_total+ntiles_water), STAT = ist)
    cf_desc    = t_cf_var('t_g_t', 'K', 'weighted surface temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_g_t'//suffix, p_prog_lnd%t_g_t,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,           &
         & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  


    ! fill the separate variables belonging to the container t_gt
    ALLOCATE(p_prog_lnd%t_gt_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        NULLIFY(p_prog_lnd%t_gt_ptr(jsfc)%p_2d, p_prog_lnd%t_gt_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'t_g_t'//suffix,                &
               & vname_prefix//'t_g_t_'//TRIM(ADJUSTL(csfc))//suffix,          &
               & p_prog_lnd%t_gt_ptr(jsfc)%p_2d,                               &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
               & t_cf_var('t_g_t_'//TRIM(csfc), '', '', datatype_flt),       &
               & grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
               & ldims=shape2d,                                                &
               & var_class=CLASS_TILE,                                         &
               & tlev_source=TLEV_NNOW_RCF,                                    & ! for output take field from nnow_rcf slice
               & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
      ENDDO

    IF ( atm_phy_nwp_config(p_jg)%inwp_surface > 0 ) THEN

    ! & p_prog_lnd%t_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('t_s_t', 'K', 'temperature of ground surface', datatype_flt)
    grib2_desc = grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_s_t'//suffix, p_prog_lnd%t_s_t,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

    ! fill the separate variables belonging to the container t_s
    ALLOCATE(p_prog_lnd%t_s_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      NULLIFY(p_prog_lnd%t_s_ptr(jsfc)%p_2d, p_prog_lnd%t_s_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'t_s_t'//suffix,              &
           & vname_prefix//'t_s_t_'//TRIM(ADJUSTL(csfc))//suffix,          &
           & p_prog_lnd%t_s_ptr(jsfc)%p_2d,                                &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('t_s_t_'//csfc, '', '', datatype_flt),             &
           & grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
           & ldims=shape2d,                                                &
           & var_class=CLASS_TILE,                                         &
           & tlev_source=TLEV_NNOW_RCF, in_group=groups("land_tile_vars") ) ! for output take field from nnow_rcf slice
    ENDDO




    ! & p_prog_lnd%w_i_t(nproma,nblks_c,ntiles_total)
    cf_desc     = t_cf_var('w_i_t', 'm H2O', 'water content of interception water', datatype_flt)
    new_cf_desc = t_cf_var('w_i_t', 'kg/m**2', 'weighted water content of interception water', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 0, 13, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_i_t'//suffix, p_prog_lnd%w_i_t,     &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,               &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

    ! fill the separate variables belonging to the container w_i
    ALLOCATE(p_prog_lnd%w_i_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%w_i_ptr(jsfc)%p_2d, p_prog_lnd%w_i_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_i_t'//suffix,                &
           & vname_prefix//'w_i_t_'//TRIM(ADJUSTL(csfc))//suffix,            &
           & p_prog_lnd%w_i_ptr(jsfc)%p_2d,                                  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
           & t_cf_var('w_i_t_'//csfc, '', '', datatype_flt),               &
           & grib2_var(2, 0, 13, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
           & ldims=shape2d,                                                  &
           & var_class=CLASS_TILE_LAND,                                      &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t"),          &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )
    ENDDO


    IF (itype_interception == 2) THEN

      ! & p_prog_lnd%w_p_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('w_p_t', 'm H2O', 'water content of interception water', datatype_flt)
      grib2_desc = grib2_var(2, 0, 14, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'w_p_t'//suffix, p_prog_lnd%w_p_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

      ! fill the separate variables belonging to the container w_p
      ALLOCATE(p_prog_lnd%w_p_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%w_p_ptr(jsfc)%p_2d, p_prog_lnd%w_p_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'w_p_t'//suffix,                &
             & vname_prefix//'w_p_t_'//TRIM(ADJUSTL(csfc))//suffix,            &
             & p_prog_lnd%w_p_ptr(jsfc)%p_2d,                                  &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
             & t_cf_var('w_p_t_'//csfc, '', '', datatype_flt),               &
             & grib2_var(2, 0, 14, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
             & ldims=shape2d,                                                  &
             & var_class=CLASS_TILE_LAND,                                      &
             & tlev_source=TLEV_NNOW_RCF,                                      &
             & in_group=groups("land_tile_vars")) ! for output take field from nnow_rcf slice
      ENDDO

      ! & p_prog_lnd%w_s_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('w_s_t', 'm H2O', 'water content of interception water', datatype_flt)
      grib2_desc = grib2_var(2, 0, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'w_s_t'//suffix, p_prog_lnd%w_s_t,     &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,grib2_desc,                 &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

      ! fill the separate variables belonging to the container w_s
      ALLOCATE(p_prog_lnd%w_s_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%w_s_ptr(jsfc)%p_2d, p_prog_lnd%w_s_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'w_s_t'//suffix,                &
             & vname_prefix//'w_s_t_'//TRIM(ADJUSTL(csfc))//suffix,            &
             & p_prog_lnd%w_s_ptr(jsfc)%p_2d,                                  &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
             & t_cf_var('w_s_t_'//csfc, '', '', datatype_flt),               &
             & grib2_var(2, 0, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
             & ldims=shape2d,                                                  &
             & var_class=CLASS_TILE_LAND,                                      &
             & tlev_source=TLEV_NNOW_RCF,                                      &
             & in_group=groups("land_tile_vars")) ! for output take field from nnow_rcf slice
      ENDDO
    END IF  ! itype_interception == 2


    ! & p_prog_lnd%t_so_t(nproma,nlev_soil+1,nblks_c,ntiles_total) 
    cf_desc    = t_cf_var('t_so_t', 'K', 'soil temperature (main level)', datatype_flt)
    grib2_desc = grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_so_t'//suffix, p_prog_lnd%t_so_t,  &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND_P1, cf_desc, grib2_desc,  &
         & ldims=(/nproma,nlev_soil+1,kblks,ntiles_total/),                      &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

    ! fill the separate variables belonging to the container t_so
    !
    ALLOCATE(p_prog_lnd%t_so_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%t_so_ptr(jsfc)%p_2d, p_prog_lnd%t_so_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'t_so_t'//suffix,               &
           & vname_prefix//'t_so_t_'//TRIM(ADJUSTL(csfc))//suffix,           &
           & p_prog_lnd%t_so_ptr(jsfc)%p_3d,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND_P1,                 &
           & t_cf_var('t_so_t_'//csfc, '', '', datatype_flt),              &
           & grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
           & ldims=(/nproma,nlev_soil+1,kblks/),                             &
           & var_class=CLASS_TILE_LAND,                                      &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
    ENDDO



   ! & p_prog_lnd%w_so_t(nproma,nlev_soil,nblks_c,ntiles_total)
    cf_desc     = t_cf_var('w_so_t', 'm H20', 'total water content (ice + liquid water)', &
         &                datatype_flt)

    new_cf_desc = t_cf_var('w_so_t', 'kg/m**2', 'total water content (ice + liquid water)', datatype_flt)
    grib2_desc  = grib2_var(2, 3, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_so_t'//suffix, p_prog_lnd%w_so_t,  &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc,     &
         & ldims=(/nproma,nlev_soil,kblks,ntiles_total/),                        &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container w_so
    ALLOCATE(p_prog_lnd%w_so_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%w_so_ptr(jsfc)%p_2d, p_prog_lnd%w_so_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_so_t'//suffix,               &
           & vname_prefix//'w_so_t_'//TRIM(ADJUSTL(csfc))//suffix,           &
           & p_prog_lnd%w_so_ptr(jsfc)%p_3d,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                    &
           & t_cf_var('w_so_t_'//csfc, '', '', datatype_flt),              &
           & grib2_var(2, 3, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
           & ldims=(/nproma,nlev_soil,kblks/),                               &
           & var_class=CLASS_TILE_LAND,                                      &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ), & 
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t"),          &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )
    ENDDO



    ! & p_prog_lnd%w_so_ice_t(nproma,nlev_soil,nblks_c,ntiles_total)
    cf_desc     = t_cf_var('w_so_ice_t', 'm H20', 'ice content', datatype_flt)
    new_cf_desc = t_cf_var('w_so_ice_t', 'kg/m**2', 'ice content', datatype_flt)
    grib2_desc = grib2_var(2, 3, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_so_ice_t'//suffix,                 &
         & p_prog_lnd%w_so_ice_t, GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,   &
         & cf_desc, grib2_desc, ldims=(/nproma,nlev_soil,kblks,ntiles_total/),   &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container w_so_ice
    ALLOCATE(p_prog_lnd%w_so_ice_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%w_so_ice_ptr(jsfc)%p_2d, p_prog_lnd%w_so_ice_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_so_ice_t'//suffix,           &
           & vname_prefix//'w_so_ice_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
           & p_prog_lnd%w_so_ice_ptr(jsfc)%p_3d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                    &
           & t_cf_var('w_so_ice_t_'//csfc, '', '', datatype_flt),          &
           & grib2_var(2, 3, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
           & ldims=(/nproma,nlev_soil,kblks/),                               &
           & var_class=CLASS_TILE_LAND,                                      &
           & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ), &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t"),          &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )
    ENDDO




    ! & p_prog_lnd%t_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('t_snow_t', 'K', 'temperature of the snow-surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_snow_t'//suffix, p_prog_lnd%t_snow_t,&
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,               &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. ) 

    ! fill the separate variables belonging to the container t_snow
    ALLOCATE(p_prog_lnd%t_snow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%t_snow_ptr(jsfc)%p_2d, p_prog_lnd%t_snow_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'t_snow_t'//suffix,             &
             & vname_prefix//'t_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
             & p_prog_lnd%t_snow_ptr(jsfc)%p_2d,                             &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
             & t_cf_var('t_snow_t_'//csfc, '', '', datatype_flt),          &
             & grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
             & ldims=shape2d,                                                &
             & var_class=CLASS_TILE_LAND,                                    &
             & tlev_source=TLEV_NNOW_RCF,                                    & ! for output take field from nnow_rcf slice
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
    ENDDO


    ! & p_prog_lnd%w_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc     = t_cf_var('w_snow_t', 'm H2O', 'water equivalent of snow', datatype_flt)
    new_cf_desc = t_cf_var('w_snow_t', 'kg m-2', 'water equivalent of snow', datatype_flt)
    grib2_desc = grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_snow_t'//suffix, p_prog_lnd%w_snow_t,&
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,               &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container w_snow
    ALLOCATE(p_prog_lnd%w_snow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%w_snow_ptr(jsfc)%p_2d, p_prog_lnd%w_snow_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_snow_t'//suffix,           &
           & vname_prefix//'w_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
           & p_prog_lnd%w_snow_ptr(jsfc)%p_2d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('w_snow_t_'//csfc, '', '', datatype_flt),          &
           & grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
           & ldims=shape2d,                                                &
           & var_class=CLASS_TILE_LAND,                                    &
           & tlev_source=TLEV_NNOW_RCF,                                    & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t"),        &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )
    ENDDO


    ! & p_prog_lnd%rho_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('rho_snow_t', 'kg/m**3', 'snow density', datatype_flt)
    grib2_desc = grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'rho_snow_t'//suffix, p_prog_lnd%rho_snow_t, &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,                    &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )
    

    ! fill the separate variables belonging to the container rho_snow
    ALLOCATE(p_prog_lnd%rho_snow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%rho_snow_ptr(jsfc)%p_2d, p_prog_lnd%rho_snow_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( prog_list, vname_prefix//'rho_snow_t'//suffix,           &
           & vname_prefix//'rho_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
           & p_prog_lnd%rho_snow_ptr(jsfc)%p_2d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
           & t_cf_var('rho_snow_t_'//csfc, '', '', datatype_flt),          &
           & grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
           & ldims=shape2d,                                                  &
           & var_class=CLASS_TILE_LAND,                                      &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
    END DO


    IF (lmulti_snow .OR. l2lay_rho_snow) THEN

      ! & p_prog_lnd%rho_snow_mult_t(nproma,nlev_snow,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('rho_snow_mult_t', 'kg/m**3', 'snow density', datatype_flt)
      grib2_desc = grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'rho_snow_mult_t'//suffix,              &
           & p_prog_lnd%rho_snow_mult_t, GRID_UNSTRUCTURED_CELL, ZA_SNOW,            &
           & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                           &
           & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ! fill the separate variables belonging to the container rho_snow_mult
      !
      ALLOCATE(p_prog_lnd%rho_snow_mult_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%rho_snow_mult_ptr(jsfc)%p_2d, p_prog_lnd%rho_snow_mult_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'rho_snow_mult_t'//suffix,       &
             & vname_prefix//'rho_snow_mult_t_'//TRIM(ADJUSTL(csfc))//suffix,   &
             & p_prog_lnd%rho_snow_mult_ptr(jsfc)%p_3d,                         &
             & GRID_UNSTRUCTURED_CELL, ZA_SNOW,                                 &
             & t_cf_var('rho_snow_mult_t_'//csfc, '', '', datatype_flt),      &
             & grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
             & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,               &
             & var_class=CLASS_TILE_LAND,                                       &
             & tlev_source=TLEV_NNOW_RCF,                                       &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice 
      ENDDO

    ENDIF

    IF (lmulti_snow) THEN

      ! & p_prog_lnd%t_snow_mult_t(nproma,nlev_snow+1,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('t_snow_mult_t', 'K', 'temperature of the snow-surface', &
           & datatype_flt)
      grib2_desc = grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'t_snow_mult_t'//suffix, p_prog_lnd%t_snow_mult_t, &
       & GRID_UNSTRUCTURED_CELL, ZA_SNOW_HALF, cf_desc, grib2_desc,               &
       & ldims=(/nproma,nlev_snow+1,kblks,ntiles_total/),                         &
       & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. ) 

      ! fill the separate variables belonging to the container t_snow_mult
      !
      ALLOCATE(p_prog_lnd%t_snow_mult_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%t_snow_mult_ptr(jsfc)%p_2d, p_prog_lnd%t_snow_mult_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'t_snow_mult_t'//suffix,      &
             & vname_prefix//'t_snow_mult_t_'//TRIM(ADJUSTL(csfc))//suffix,  &
             & p_prog_lnd%t_snow_mult_ptr(jsfc)%p_3d,                        &
             & GRID_UNSTRUCTURED_CELL, ZA_SNOW_HALF,                         &
             & t_cf_var('t_snow_mult_t_'//csfc, '', '', datatype_flt),     &
             & grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
             & ldims=(/nproma,nlev_snow+1,kblks/),                           &
             & var_class=CLASS_TILE_LAND,                                    &
             & lrestart=.TRUE.,                                              &
             & tlev_source=TLEV_NNOW_RCF,                                    &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice 
      ENDDO



      ! & p_prog_lnd%wtot_snow_t(nproma,nlev_snow,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('wtot_snow_t', 'm H2O', 'total water content in snow', datatype_flt)
      grib2_desc = grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'wtot_snow_t'//suffix,                &
           & p_prog_lnd%wtot_snow_t, GRID_UNSTRUCTURED_CELL, ZA_SNOW,              &
           & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                         &
           & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ! fill the separate variables belonging to the container wtot_snow
      !
      ALLOCATE(p_prog_lnd%wtot_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%wtot_snow_ptr(jsfc)%p_2d, p_prog_lnd%wtot_snow_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'wtot_snow_t'//suffix,          &
             & vname_prefix//'wtot_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,      &
             & p_prog_lnd%wtot_snow_ptr(jsfc)%p_3d,                            &
             & GRID_UNSTRUCTURED_CELL, ZA_SNOW,                                &
             & t_cf_var('wtot_snow_t_'//csfc, '', '', datatype_flt),           &
             & grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
             & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,              &
             & var_class=CLASS_TILE_LAND, tlev_source=TLEV_NNOW_RCF,           &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice
      ENDDO



      ! & p_prog_lnd%wliq_snow_t(nproma,nlev_snow,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('wliq_snow_t', 'm H2O', 'liquid water content in snow', datatype_flt)
      grib2_desc = grib2_var(0, 1, 210, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'wliq_snow_t'//suffix,                &
           & p_prog_lnd%wliq_snow_t, GRID_UNSTRUCTURED_CELL, ZA_SNOW,              & 
           & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                         &
           & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ! fill the separate variables belonging to the container wliq_snow
      !
      ALLOCATE(p_prog_lnd%wliq_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%wliq_snow_ptr(jsfc)%p_2d, p_prog_lnd%wliq_snow_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'wliq_snow_t'//suffix,          &
             & vname_prefix//'wliq_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,      &
             & p_prog_lnd%wliq_snow_ptr(jsfc)%p_3d,                            &
             & GRID_UNSTRUCTURED_CELL, ZA_SNOW,                                &
             & t_cf_var('wliq_snow_t_'//csfc, '', '', datatype_flt),           &
             & grib2_var(0, 1, 210, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
             & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,              &
             & var_class=CLASS_TILE_LAND, tlev_source=TLEV_NNOW_RCF,           &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice
      ENDDO


      ! & p_prog_lnd%dzh_snow_t(nproma,nlev_snow,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('dzh_snow_t', 'm', 'layer thickness between half levels in snow', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'dzh_snow_t'//suffix,                 &
           & p_prog_lnd%dzh_snow_t, GRID_UNSTRUCTURED_CELL, ZA_SNOW,               &
           & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                         &
           & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ! fill the separate variables belonging to the container dzh_snow
      !
      ALLOCATE(p_prog_lnd%dzh_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%dzh_snow_ptr(jsfc)%p_2d, p_prog_lnd%dzh_snow_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'dzh_snow_t'//suffix,             &
               & vname_prefix//'dzh_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
               & p_prog_lnd%dzh_snow_ptr(jsfc)%p_3d,                             &
               & GRID_UNSTRUCTURED_CELL, ZA_SNOW,                                &
               & t_cf_var('dzh_snow_t_'//csfc, '', '', datatype_flt),            &
               & grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
               & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,              &
               & var_class=CLASS_TILE_LAND, tlev_source=TLEV_NNOW_RCF,           &
               & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice
      ENDDO

    ENDIF  ! lmulti_snow


    p_prog_lnd%t_so_t(:,1,:,:) = 290.4_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,2,:,:) = 290.4_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,3,:,:) = 290.7_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,4,:,:) = 291.4_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,5,:,:) = 292.7_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,6,:,:) = 293.2_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,7,:,:) = 291.1_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,8,:,:) = 283.1_wp !!HW: be careful about the indices!!
    p_prog_lnd%t_so_t(:,9,:,:) = 282.1_wp !!HW: be careful about the indices!!

    !p_prog_lnd%w_so_t(:,1,:,:) = 1.2E-3_wp*1._wp           !! JH [kg/m**2]
    p_prog_lnd%w_so_t(:,1,:,:) = les_config(p_jg)%wso_mean(1)         !! JH [kg/m**2]
    p_prog_lnd%w_so_t(:,2,:,:) = p_prog_lnd%w_so_t(:,1,:,:)*2._wp     !! JH
    p_prog_lnd%w_so_t(:,3,:,:) = p_prog_lnd%w_so_t(:,1,:,:)*6._wp     !! JH
    p_prog_lnd%w_so_t(:,4,:,:) = p_prog_lnd%w_so_t(:,1,:,:)*18._wp    !! JH
    p_prog_lnd%w_so_t(:,5,:,:) = p_prog_lnd%w_so_t(:,1,:,:)*54._wp    !! JH
    p_prog_lnd%w_so_t(:,6,:,:) = p_prog_lnd%w_so_t(:,1,:,:)*162._wp   !! JH
    p_prog_lnd%w_so_t(:,7,:,:) = p_prog_lnd%w_so_t(:,1,:,:)*486._wp   !! JH
    p_prog_lnd%w_so_t(:,8,:,:) = p_prog_lnd%w_so_t(:,1,:,:)*1458._wp  !! JH

    maxlat_loc = MAXVAL(p_patch(p_jg)%cells%center(:,:)%lat)
    minlat_loc = MINVAL(p_patch(p_jg)%cells%center(:,:)%lat)

    ! comparing zdist_cmp over all domains
    ! local variable mproc_id must be pre-set to p_pe in all domains
    mproc_id = p_pe
    maxlat = global_max(maxlat_loc,mproc_id)
    minlat = global_max(minlat_loc,mproc_id)

    maxlon = MAXVAL(p_patch(p_jg)%cells%center(:,:)%lon)
    minlon = MINVAL(p_patch(p_jg)%cells%center(:,:)%lon)

    coeff_A    = p_prog_lnd%w_so_t(1,:,1,1)
    !coeff_B    = les_config(p_jg)%wso_gradient
    wso_fun    = les_config(p_jg)%wso_fun
    wso_cyc    = les_config(p_jg)%wso_cyc
    wso_min    = les_config(p_jg)%wso_min
    wso_max    = les_config(p_jg)%wso_max
    wso_spread = les_config(p_jg)%wso_spread
    !wso_mean   = les_config(p_jg)%wso_mean
    wso_max(1)   = les_config(p_jg)%wso_max(1) * 1._wp
    wso_max(2)   = les_config(p_jg)%wso_max(1) * 2._wp
    wso_max(3)   = les_config(p_jg)%wso_max(1) * 6._wp
    wso_max(4)   = les_config(p_jg)%wso_max(1) * 18._wp
    wso_max(5)   = les_config(p_jg)%wso_max(1) * 54._wp
    wso_max(6)   = les_config(p_jg)%wso_max(1) * 162._wp
    wso_max(7)   = les_config(p_jg)%wso_max(1) * 486._wp
    wso_max(8)   = les_config(p_jg)%wso_max(1) * 1458._wp
    wso_mean(1)   = les_config(p_jg)%wso_mean(1) * 1._wp
    wso_mean(2)   = les_config(p_jg)%wso_mean(1) * 2._wp
    wso_mean(3)   = les_config(p_jg)%wso_mean(1) * 6._wp
    wso_mean(4)   = les_config(p_jg)%wso_mean(1) * 18._wp
    wso_mean(5)   = les_config(p_jg)%wso_mean(1) * 54._wp
    wso_mean(6)   = les_config(p_jg)%wso_mean(1) * 162._wp
    wso_mean(7)   = les_config(p_jg)%wso_mean(1) * 486._wp
    wso_mean(8)   = les_config(p_jg)%wso_mean(1) * 1458._wp
    wso_grad(1)   = les_config(p_jg)%wso_gradient * 1._wp
    wso_grad(2)   = les_config(p_jg)%wso_gradient * 2._wp
    wso_grad(3)   = les_config(p_jg)%wso_gradient * 6._wp
    wso_grad(4)   = les_config(p_jg)%wso_gradient * 18._wp
    wso_grad(5)   = les_config(p_jg)%wso_gradient * 54._wp
    wso_grad(6)   = les_config(p_jg)%wso_gradient * 162._wp
    wso_grad(7)   = les_config(p_jg)%wso_gradient * 486._wp
    wso_grad(8)   = les_config(p_jg)%wso_gradient * 1458._wp
  
    wso_grad_fun4(:) = (wso_max(:) - wso_min(:))/maxlat * 1.0001_wp
    l(:) = (wso_max(:) - wso_min(:)) / wso_grad(:)
    l_fun4(:) = (wso_max(:) - wso_min(:)) / wso_grad_fun4(:)
    r(:) = 0.5_wp*((wso_mean(:)-wso_min(:))*(maxlat-minlat) / (wso_max(:)-wso_min(:)) - l)
    !WRITE(sval,*) wso_grad_fun4(:)
    !CALL message('Cunbo, out loop, wso_grad_fun4: ', sval)
    !WRITE(sval,*) l_fun4(:)
    !CALL message('Cunbo, out loop, l_fun4: ', sval)


    !wso_grad_fun5(:) = 2.0_wp*((wso_max(:) - wso_min(:))/maxlat) * 1.0001_wp
    !l_fun5(:) = (wso_max(:) - wso_min(:)) / wso_grad_fun5(:)
    !WRITE(sval,*) wso_grad_fun5(:)
    !CALL message('Cunbo, out loop, wso_grad_fun5: ', sval)
    !WRITE(sval,*) l_fun5(:)
    !CALL message('Cunbo, out loop, l_fun5: ', sval)


    !SBr, CHa
    !DO jk = 1, nlev_soil
    !  IF (2_wp*(r(jk)+l(jk)) > wso_spread) THEN
    !   WRITE(sval,*) '2_wp*(r+l) > wso_spread at layer', jk, r(jk), ' ', l(jk)
    !    CALL finish('SBr, CHa: init moisture', sval)
    !  END IF
    !  IF (r(jk)<0_wp .OR. l(jk)<0._wp) THEN
    !    WRITE(sval,*) 'r<0_wp .OR. l<0._wp at layer', jk, r(jk), ' ', l(jk)
    !    CALL finish('SBr, CHa: init moisture', sval)
    !  END IF
    !END DO

    ! SBr, CHa
    ! compute mean
    !mean_wso = 0._wp
    !tot_area = 0._wp
    !DO jb = 1, p_patch(p_jg)%nblks_c
    !  DO jt = 1, ntiles_total
    !    DO jk = 1, nlev_soil
    !      DO ic = 1, ext_data(p_jg)%atm%lp_count_t(jb,jt)
    !        jc = ext_data(p_jg)%atm%idx_lst_lp_t(ic,jb,jt)
    !        point = p_patch(p_jg)%cells%center(jc,jb)
    !        IF (abs(point%lat)> r(jk)+l(jk)) THEN
    !          mean_wso(jk) = mean_wso(jk) + wso_min(jk) * &
    !                         p_patch(p_jg)%cells%area(jc,jb)
    !        ELSE IF (abs(point%lat)< r(jk)) THEN
    !          mean_wso(jk) = mean_wso(jk) + wso_max(jk) * &
    !                         p_patch(p_jg)%cells%area(jc,jb)
    !        ELSE IF (point%lat > 0) THEN
    !          mean_wso(jk) = mean_wso(jk) + &
    !                         ((-wso_max(jk)*(point%lat-r(jk)-l(jk)) + &
    !                          wso_min(jk)*(point%lat-r(jk)))/l(jk)) *p_patch(p_jg)%cells%area(jc,jb)
    !        ELSE
    !          mean_wso(jk) = mean_wso(jk) + &
    !                         ((wso_max(jk)*(point%lat+r(jk)+l(jk)) - &
    !                         wso_min(jk)*(point%lat+r(jk)))/l(jk)) *p_patch(p_jg)%cells%area(jc,jb)
    !        END IF
    !        tot_area = tot_area + p_patch(p_jg)%cells%area(jc,jb)
    !      END DO
    !    END DO
    !  END DO
    !END DO
    !
    !mean_wso = mean_wso / tot_area

    ! add heterogeneous soil moisture here
    DO jb = 1, p_patch(p_jg)%nblks_c
      DO jt = 1, ntiles_total
        DO jk = 1, nlev_soil
          DO ic = 1, ext_data(p_jg)%atm%lp_count_t(jb,jt)
            jc = ext_data(p_jg)%atm%idx_lst_lp_t(ic,jb,jt)
            point = p_patch(p_jg)%cells%center(jc,jb)

            IF (wso_fun == 1) THEN !1=homo.
              p_prog_lnd%w_so_t(jc,jk,jb,jt) = p_prog_lnd%w_so_t(jc,jk,jb,jt)

            ELSEIF(wso_fun == 2) THEN !2=function2
              !Function2
              IF (abs(point%lat)> r(jk)+l(jk)) THEN
                p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_min(jk)
              ELSE IF (abs(point%lat)< r(jk)) THEN
                p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_max(jk)
              ELSE IF (0 < point%lat ) THEN
                p_prog_lnd%w_so_t(jc,jk,jb,jt) = (-wso_max(jk)*(point%lat-r(jk)-l(jk)) + wso_min(jk)*(point%lat-r(jk)))/l(jk)
              ELSE
                p_prog_lnd%w_so_t(jc,jk,jb,jt) = (wso_max(jk)*(point%lat+r(jk)+l(jk)) - wso_min(jk)*(point%lat+r(jk)))/l(jk)
              END IF
            
            ELSEIF(wso_fun == 3) THEN !3=function3
              !Function3
              IF (abs(point%lat)>= l(jk)) THEN
                p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_min(jk)
              ELSE IF (0 < point%lat ) THEN
                p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_min(jk)
              ELSE
                p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_min(jk)
              END IF

            ELSEIF(wso_fun == 4) THEN !4=function4
              !Function4
              !wso_grad(:) = (wso_max(:) - wso_min(:))/maxlat * 1.001_wp
              !PRINT *, "wso_max: ", wso_max
              !PRINT *, "wso_min: ", wso_min
              !PRINT *, "maxlat: ", maxlat
              !WRITE(sval,*) wso_grad_fun4(:)
              !CALL message('Cunbo, in loop, wso_grad_fun4: ', sval)

              IF (abs(point%lat)>= l_fun4(jk)) THEN
                p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_min(jk)
              ELSE IF (0 < point%lat ) THEN
                p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun4(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_min(jk)
              ELSE
                p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun4(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_min(jk)
              END IF

            ELSEIF(wso_fun == 5) THEN !5=function5
              IF (wso_cyc == 1) THEN
                wso_grad_fun5(:) = wso_cyc*((wso_max(:) - wso_min(:))/maxlat) * 1.0001_wp
                l_fun5(:) = (wso_max(:) - wso_min(:)) / wso_grad_fun5(:)
                IF (abs(point%lat)>= wso_cyc*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_min(jk)
                ELSE IF (0 < point%lat ) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSE
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_min(jk)
                END IF

              ELSEIF (wso_cyc == 2) THEN
                wso_grad_fun5(:) = wso_cyc*((wso_max(:) - wso_min(:))/maxlat) * 1.0001_wp
                l_fun5(:) = (wso_max(:) - wso_min(:)) / wso_grad_fun5(:)
                !WRITE(sval,*) wso_grad_fun5(:)
                !CALL message('Cunbo, in loop, wso_grad_fun5: ', sval)
                !WRITE(sval,*) l_fun5(:)
                !CALL message('Cunbo, in loop, l_fun5: ', sval)

                IF (abs(point%lat) >= wso_cyc*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_min(jk)
                ELSE IF (point%lat >= l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSE IF (point%lat >= 0 .AND. point%lat < l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + wso_min(jk)
                ELSE IF (point%lat >= -l_fun5(jk) .AND. point%lat < 0) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + wso_min(jk)
                ELSE
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_max(jk)
                END IF


              ELSEIF (wso_cyc == 3) THEN
                wso_grad_fun5(:) = wso_cyc*((wso_max(:) - wso_min(:))/maxlat) * 1.0001_wp
                l_fun5(:) = (wso_max(:) - wso_min(:)) / wso_grad_fun5(:)
                IF (abs(point%lat) >= wso_cyc*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_min(jk)
                ELSEIF (point%lat >= 2.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= l_fun5(jk) .AND. point%lat < 2.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat - (wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= 0 .AND. point%lat < l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -l_fun5(jk) .AND. point%lat < 0) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -2.0_wp*l_fun5(jk) .AND. point%lat < -l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat - 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSE
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 3.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ENDIF

              ELSEIF (wso_cyc == 4) THEN
                wso_grad_fun5(:) = wso_cyc*((wso_max(:) - wso_min(:))/maxlat) * 1.0001_wp
                l_fun5(:) = (wso_max(:) - wso_min(:)) / wso_grad_fun5(:)
                IF (abs(point%lat)>= wso_cyc*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_min(jk)
                ELSEIF (point%lat >= 3.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + 3.0_wp*(wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= 2.0_wp*l_fun5(jk) .AND. point%lat < 3.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat - 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= l_fun5(jk) .AND. point%lat < 2.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= 0 .AND. point%lat < l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + wso_min(jk)
                ELSEIF (point%lat >= -l_fun5(jk) .AND. point%lat < 0) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + wso_min(jk)
                ELSEIF (point%lat >= -2.0_wp*l_fun5(jk) .AND. point%lat < -l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -3.0_wp*l_fun5(jk) .AND. point%lat < -2.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat - 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSE
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 4.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ENDIF

              ELSEIF (wso_cyc == 6) THEN
                wso_grad_fun5(:) = wso_cyc*((wso_max(:) - wso_min(:))/maxlat) * 1.0001_wp
                l_fun5(:) = (wso_max(:) - wso_min(:)) / wso_grad_fun5(:)

                IF (abs(point%lat)>= wso_cyc*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_min(jk)
                ELSEIF (point%lat >= 5.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + 5.0_wp*(wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= 4.0_wp*l_fun5(jk) .AND. point%lat < 5.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat - 4.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= 3.0_wp*l_fun5(jk) .AND. point%lat < 4.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + 3.0_wp*(wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= 2.0_wp*l_fun5(jk) .AND. point%lat < 3.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat - 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= 1.0_wp*l_fun5(jk) .AND. point%lat < 2.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= 0 .AND. point%lat < 1.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + wso_min(jk)
                ELSEIF (point%lat >= -1.0_wp*l_fun5(jk) .AND. point%lat < 0) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + wso_min(jk)
                ELSEIF (point%lat >= -2.0_wp*l_fun5(jk) .AND. point%lat < -1.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -3.0_wp*l_fun5(jk) .AND. point%lat < -2.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat - 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -4.0_wp*l_fun5(jk) .AND. point%lat < -3.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 4.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -5.0_wp*l_fun5(jk) .AND. point%lat < -4.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat - 4.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSE
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 6.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ENDIF

              ELSEIF (wso_cyc == 8) THEN
                wso_grad_fun5(:) = wso_cyc*((wso_max(:) - wso_min(:))/maxlat) * 1.0001_wp
                l_fun5(:) = (wso_max(:) - wso_min(:)) / wso_grad_fun5(:)

                IF (abs(point%lat)>= wso_cyc*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = wso_min(jk)
                ELSEIF (point%lat >= 7.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + 7.0_wp*(wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= 6.0_wp*l_fun5(jk) .AND. point%lat < 7.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat - 6.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= 5.0_wp*l_fun5(jk) .AND. point%lat < 6.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + 5.0_wp*(wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= 4.0_wp*l_fun5(jk) .AND. point%lat < 5.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat - 4.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= 3.0_wp*l_fun5(jk) .AND. point%lat < 4.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + 3.0_wp*(wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= 2.0_wp*l_fun5(jk) .AND. point%lat < 3.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat - 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= 1.0_wp*l_fun5(jk) .AND. point%lat < 2.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + (wso_max(jk) - wso_min(jk)) + wso_max(jk)
                ELSEIF (point%lat >= 0 .AND. point%lat < 1.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + wso_min(jk)
                ELSEIF (point%lat >= -1.0_wp*l_fun5(jk) .AND. point%lat < 0) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat + wso_min(jk)
                ELSEIF (point%lat >= -2.0_wp*l_fun5(jk) .AND. point%lat < -1.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -3.0_wp*l_fun5(jk) .AND. point%lat < -2.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat - 2.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -4.0_wp*l_fun5(jk) .AND. point%lat < -3.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 4.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -5.0_wp*l_fun5(jk) .AND. point%lat < -4.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat - 4.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -6.0_wp*l_fun5(jk) .AND. point%lat < -5.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 6.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSEIF (point%lat >= -7.0_wp*l_fun5(jk) .AND. point%lat < -6.0_wp*l_fun5(jk)) THEN
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) = -wso_grad_fun5(jk)*point%lat - 6.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ELSE
                  p_prog_lnd%w_so_t(jc,jk,jb,jt) =  wso_grad_fun5(jk)*point%lat + 8.0_wp*(wso_max(jk) - wso_min(jk)) + wso_min(jk)
                ENDIF

              ENDIF !wso_cyc

            ELSEIF(wso_grad(1) <= 0) THEN  !negetative grad.
              p_prog_lnd%w_so_t(jc,jk,jb,jt) = p_prog_lnd%w_so_t(jc,jk,jb,jt)
            
            END IF !wso_fun
            
            ! SBr, CHa: check that the value does not get limited from below in
            ! the later initialization of w_so
            IF( p_prog_lnd%w_so_t(jc,jk,jb,jt) < 0_wp ) &
              CALL finish('ERROR assigning p_prog_lnd%w_so_t:: w_so-t is negative! Please increase wso_mean in les_nml! Exiting...')

            !p_prog_lnd%w_so_t(jc,jk,jb,jt) = coeff_A(jk) * EXP(-(coeff_B * point%lat)**2) &
            !  - mean_wso(jk) + les_config(p_jg)%wso_mean(jk)
            !p_prog_lnd%w_so_t(jc,jk,jb,jt) = p_prog_lnd%w_so_t(jc,jk,jb,jt) * 2.0_wp * EXP(-2.0_wp * (point%lat*les_config(p_jg)%wso_gradient)**2)
            !p_prog_lnd%w_so_t(jc,jk,jb,jt) = p_prog_lnd%w_so_t(jc,jk,jb,jt) * 2.0_wp * EXP(-2.0_wp * (point%lat*les_config(p_jg)%wso_gradient)**2) - mean_wso + les_config(p_jg)%wso_mean
            !p_prog_lnd%w_so_t(jc,jk,jb,jt) = p_prog_lnd%w_so_t(jc,jk,jb,jt) *1.0_wp * EXP(-2.0_wp * ((point%lon-0.5*(maxlon+minlon))/3.1415925)**2)
            !p_prog_lnd%w_so_t(jc,jk,jb,jt) = (0.1E-3_wp*(maxlat-point%lat) + 1E-2_wp*(point%lat-minlat))/(maxlat-minlat)
          END DO
        END DO
      END DO
    END Do


    END IF !inwp_surface > 0

  END SUBROUTINE new_nwp_lnd_prog_list



  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of prognostic lake/seaice state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl, DWD (2012-08-08)
  !!
  !! Modifications by Dmitrii Mironov, DWD (2016-08-04)
  !! - Prognostic sea-ice albedo is added.
  !!
  SUBROUTINE new_nwp_wtr_prog_list( p_jg, kblks, listname, vname_prefix, &
    &                               prog_list, p_prog_wtr, timelev )

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes
    INTEGER,INTENT(IN) ::  p_jg  !< patch id

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix

    TYPE(t_var_list),INTENT(INOUT) :: prog_list
    TYPE(t_wtr_prog),INTENT(INOUT) :: p_prog_wtr

    INTEGER, INTENT(IN) :: timelev

    ! Local variables
    !
    TYPE(t_cf_var)    :: cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2)
    INTEGER :: ibits

    CHARACTER(len=4) :: suffix
    INTEGER          :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

!-----------------------------------------------------------------------

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d              = (/nproma, kblks/)

    ! Suffix (mandatory for time level dependent variables)

    WRITE(suffix,'(".TL",i1)') timelev


    !------------------------------
    ! Ensure that all pointers have a defined association status.
    !------------------------------
    NULLIFY(p_prog_wtr%t_ice, &
    &       p_prog_wtr%h_ice, &
    &       p_prog_wtr%t_snow_si, &
    &       p_prog_wtr%h_snow_si, &
    &       p_prog_wtr%alb_si,    &
    &       p_prog_wtr%t_snow_lk,  &
    &       p_prog_wtr%h_snow_lk,  &
    &       p_prog_wtr%t_mnw_lk,   &
    &       p_prog_wtr%t_wml_lk,   &
    &       p_prog_wtr%h_ml_lk,    &
    &       p_prog_wtr%t_bot_lk,   &
    &       p_prog_wtr%c_t_lk,     &
    &       p_prog_wtr%t_b1_lk,    &
    &       p_prog_wtr%h_b1_lk)


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( prog_list, TRIM(listname), patch_id=p_jg )
    CALL default_var_list_settings( prog_list,                 &
                                  & lrestart=.TRUE.  )

    !------------------------------

    ! h_ice-field is also needed by turbdiff irrespective of whether 
    ! the seaice model is used or not.

    !
    ! sea-ice and lake-model
    !
    ! since it is currently not envisaged to have mixed sea-lake gridpoints, t_ice 
    ! is used for both sea- and lake-ice temperatures. This is in accordance with 
    ! the COSMO implementation 
    ! & p_prog_wtr%t_ice(nproma,nblks_c)
    cf_desc    = t_cf_var('t_ice', 'K', 'sea/lake-ice temperature', datatype_flt)
    grib2_desc = grib2_var(10, 2, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_ice'//suffix, p_prog_wtr%t_ice,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & ldims=shape2d, tlev_source=TLEV_NNOW_RCF,                           &
         & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_fg_in", &
         &                 "mode_iau_old_fg_in","mode_combined_in","mode_cosmode_in") )   


    ! since it is currently not envisaged to have mixed sea-lake gridpoints, h_ice 
    ! is used for both sea- and lake-ice thickness. This is in accordance with 
    ! the COSMO implementation 
    ! & p_prog_wtr%h_ice(nproma,nblks_c)
    cf_desc    = t_cf_var('h_ice', 'm', 'sea/lake-ice depth', datatype_flt)
    grib2_desc = grib2_var(10, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'h_ice'//suffix, p_prog_wtr%h_ice,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & ldims=shape2d, tlev_source=TLEV_NNOW_RCF,                           &
         & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_fg_in", &
         &                 "mode_iau_old_fg_in","mode_combined_in","mode_cosmode_in") )   


    !
    ! sea-ice model specific
    !

    ! & p_prog_wtr%t_snow_si(nproma,nblks_c)
    cf_desc    = t_cf_var('t_snow_si', 'K', 'temperature of snow on sea ice', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_snow_si'//suffix, p_prog_wtr%t_snow_si,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,     &
         & tlev_source=TLEV_NNOW_RCF)


    ! & p_prog_wtr%h_snow_si(nproma,nblks_c)
    cf_desc    = t_cf_var('h_snow_si', 'm', 'depth of snow on sea ice', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'h_snow_si'//suffix, p_prog_wtr%h_snow_si,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,     &
         & tlev_source=TLEV_NNOW_RCF)

    ! & p_prog_wtr%alb_si(nproma,nblks_c)
    cf_desc     = t_cf_var('sea_ice_albedo', '-', 'sea ice albedo (diffuse)', datatype_flt)
    new_cf_desc = t_cf_var('sea_ice_albedo', '%', 'sea ice albedo (diffuse)', datatype_flt)
    grib2_desc = grib2_var(0, 19, 234, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'alb_si'//suffix, p_prog_wtr%alb_si,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
         & ldims=shape2d, tlev_source=TLEV_NNOW_RCF,                             &
         & initval = ALB_SI_MISSVAL,                                             &
         & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ), & 
         & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in", "mode_iau_fg_in", &
         &                 "mode_iau_old_fg_in","mode_cosmode_in"),              &
         & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc) )   


    !
    !FLAKE
    !

    IF (llake) THEN   ! lake model switched on

      ! p_prog_wtr%t_snow_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_snow_lk', 'K', 'temperature of snow on lake ice', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'t_snow_lk'//suffix, p_prog_wtr%t_snow_lk,  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,     &
           & tlev_source=TLEV_NNOW_RCF)


      ! p_prog_wtr%h_snow_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('h_snow_lk', 'm', 'depth of snow on lake ice', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'h_snow_lk'//suffix, p_prog_wtr%h_snow_lk,  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,     &
           & tlev_source=TLEV_NNOW_RCF)


      ! p_prog_wtr%t_mnw_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_mnw_lk', 'K', 'mean temperature of the water column', datatype_flt)
      grib2_desc = grib2_var(1, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
      CALL add_var( prog_list, vname_prefix//'t_mnw_lk'//suffix, p_prog_wtr%t_mnw_lk,                               &
           & GRID_UNSTRUCTURED_CELL, ZA_LAKE_BOTTOM, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF, &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in") )


      ! p_prog_wtr%t_wml_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_wml_lk', 'K', 'mixed-layer temperature', datatype_flt)
      grib2_desc = grib2_var(1, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
      CALL add_var( prog_list, vname_prefix//'t_wml_lk'//suffix, p_prog_wtr%t_wml_lk,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_MIX_LAYER, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF, &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in") )


      ! p_prog_wtr%h_ml_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('h_ml_lk', 'm', 'mixed-layer thickness', datatype_flt)
      grib2_desc = grib2_var(1, 2, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
      CALL add_var( prog_list, vname_prefix//'h_ml_lk'//suffix, p_prog_wtr%h_ml_lk,                               &
           & GRID_UNSTRUCTURED_CELL, ZA_MIX_LAYER, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF, &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in") )


      ! p_prog_wtr%t_bot_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_bot_lk', 'K', 'temperature at the water-bottom sediment interface', &
        &          datatype_flt)
      grib2_desc = grib2_var(1, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'t_bot_lk'//suffix, p_prog_wtr%t_bot_lk,                                    &
           & GRID_UNSTRUCTURED_CELL, ZA_LAKE_BOTTOM_HALF, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF, &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in") )


      ! p_prog_wtr%c_t_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('c_t_lk', '-', 'shape factor (temp. profile in lake thermocline)', &
        &          datatype_flt)
      grib2_desc = grib2_var(1, 2, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfSecondFixedSurface", 162)
      CALL add_var( prog_list, vname_prefix//'c_t_lk'//suffix, p_prog_wtr%c_t_lk,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_MIX_LAYER, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF, &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in") )


      ! p_prog_wtr%t_b1_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_b1_lk', 'K',                                         &
        &          'temperature at the bottom of the upper layer of the sediments', &
        &          datatype_flt)
      grib2_desc = grib2_var(1, 2, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'t_b1_lk'//suffix, p_prog_wtr%t_b1_lk, &
           & GRID_UNSTRUCTURED_CELL, ZA_SEDIMENT_BOTTOM_TW_HALF, cf_desc,           &
           & grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF,                  &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in") )


      ! p_prog_wtr%h_b1_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('h_b1_lk', 'm',                                         &
        &          'thickness of the upper layer of the sediments', datatype_flt)
      grib2_desc = grib2_var(1, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 165)
      CALL add_var( prog_list, vname_prefix//'h_b1_lk'//suffix, p_prog_wtr%h_b1_lk,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_LAKE_BOTTOM, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF, &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in") )

    ENDIF  ! llake

  END SUBROUTINE new_nwp_wtr_prog_list


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of diagnostic land state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by , MPI-M (2011-07-01)
  !!
  !!
  SUBROUTINE new_nwp_lnd_diag_list( p_jg, kblks, listname, vname_prefix, &
    &                               diag_list, p_diag_lnd)

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes
    INTEGER,INTENT(IN) ::  p_jg !< patch id

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix
    CHARACTER(LEN=2)            :: csfc

    TYPE(t_var_list),INTENT(INOUT) :: diag_list
    TYPE(t_lnd_diag),INTENT(INOUT) :: p_diag_lnd


    ! Local variables
    TYPE(t_cf_var)    :: cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d_subs(3), shape3d_subsw(3)
    INTEGER :: jsfc          !< tile counter
    INTEGER :: ibits
    INTEGER :: DATATYPE_PACK_VAR  !< variable "entropy" for some thermodynamic fields
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

!--------------------------------------------------------------

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    IF (gribout_config(p_jg)%lgribout_24bit) THEN  ! analysis
      ! higher accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK24
    ELSE
      ! standard accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK16
    ENDIF

    ! predefined array shapes
    shape2d       = (/nproma, kblks            /)
    shape3d_subs  = (/nproma, kblks, ntiles_total /)
    shape3d_subsw = (/nproma, kblks, ntiles_total+ntiles_water /)


    !------------------------------
    ! Ensure that all pointers have a defined association status.
    !------------------------------
    NULLIFY(p_diag_lnd%qv_s, &
    &       p_diag_lnd%t_s, &
    &       p_diag_lnd%t_seasfc, &
    &       p_diag_lnd%w_i, &
    &       p_diag_lnd%w_p, &
    &       p_diag_lnd%w_s, &
    &       p_diag_lnd%t_so, &
    &       p_diag_lnd%w_so, &
    &       p_diag_lnd%w_so_ice, &
    &       p_diag_lnd%runoff_s, &
    &       p_diag_lnd%runoff_g, &
    &       p_diag_lnd%fr_seaice, &
    &       p_diag_lnd%qv_s_t, &
    &       p_diag_lnd%runoff_s_t, &
    &       p_diag_lnd%runoff_g_t, &
    &       p_diag_lnd%rstom, &
    &       p_diag_lnd%rstom_t, &
    &       p_diag_lnd%t_snow, &
    &       p_diag_lnd%rho_snow, &
    &       p_diag_lnd%w_snow, &
    &       p_diag_lnd%h_snow, &
    &       p_diag_lnd%h_snow_t, &
    &       p_diag_lnd%freshsnow, &
    &       p_diag_lnd%freshsnow_t, &
    &       p_diag_lnd%snowfrac, &
    &       p_diag_lnd%snowfrac_t, &
    &       p_diag_lnd%snowfrac_lc_t, &
    &       p_diag_lnd%t_snow_mult, &
    &       p_diag_lnd%rho_snow_mult, &
    &       p_diag_lnd%wliq_snow, &
    &       p_diag_lnd%wtot_snow, &
    &       p_diag_lnd%dzh_snow)


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( diag_list, TRIM(listname), patch_id=p_jg )
    CALL default_var_list_settings( diag_list,                 &
                                  & lrestart=.TRUE.  )

    !------------------------------

    ! & p_diag_lnd%qv_s(nproma,nblks_c)
    cf_desc    = t_cf_var('qv_s', 'kg kg-1', 'specific humidity at the surface', datatype_flt)
    grib2_desc = grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      ! Quick hack: shortName.def should be revised, instead
      &           + t_grib2_int_key("scaleFactorOfFirstFixedSurface", 0)
    CALL add_var( diag_list, vname_prefix//'qv_s', p_diag_lnd%qv_s,          &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,        &
           & ldims=shape2d,                                                  &
           & initval=0.001_wp,                                               &
           & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in", &
           &                 "mode_iau_fg_in","mode_iau_old_fg_in",          &
           &                 "mode_combined_in","mode_cosmode_in") )      


    ! & p_diag_lnd%fr_seaice(nproma,nblks_c)
    cf_desc    = t_cf_var('fr_seaice', '1', 'fraction of sea ice', datatype_flt)
    grib2_desc = grib2_var(10, 2, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'fr_seaice', p_diag_lnd%fr_seaice,  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,                   &
           & initval=0._wp,                                                    &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_ana_in",&
           &                 "mode_iau_old_ana_in","mode_combined_in") )



    ! & p_diag_lnd%qv_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('qv_s_t', 'kg kg-1', 'specific humidity at the surface', datatype_flt)
    grib2_desc = grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'qv_s_t', p_diag_lnd%qv_s_t,        &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE.,         &
           & loutput=.FALSE.,                                                  &
           & initval=0.001_wp )

    ! fill the separate variables belonging to the container qv_s_t
    ALLOCATE(p_diag_lnd%qv_st_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      NULLIFY(p_diag_lnd%qv_st_ptr(jsfc)%p_2d, p_diag_lnd%qv_st_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'qv_s_t',                       &
             & vname_prefix//'qv_s_t_'//ADJUSTL(TRIM(csfc)),                 &
             & p_diag_lnd%qv_st_ptr(jsfc)%p_2d,                              &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
             & t_cf_var('qv_s_t_'//csfc, '', '', datatype_flt),            &
             & grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
             & var_class=CLASS_TILE,                                         &
             & ldims=shape2d,                                                &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
    ENDDO

    IF ( atm_phy_nwp_config(p_jg)%inwp_surface > 0) THEN

    ! & p_diag_lnd%t_s(nproma,nblks_c)
    cf_desc    = t_cf_var('t_s', 'K', 'weighted temperature of ground surface', datatype_flt)
    grib2_desc = grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_s', p_diag_lnd%t_s,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE. )


    ! & p_diag_lnd%t_seasfc(nproma,nblks_c)
    cf_desc    = t_cf_var('t_seasfc', 'K', 'sea surface temperature', datatype_flt)
    grib2_desc = grib2_var(10, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_seasfc', p_diag_lnd%t_seasfc,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,                       &
         & lmiss=.TRUE., missval=0.0_wp,                                         &
         & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ), & 
         & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in",                  &
         &                 "mode_iau_ana_in","mode_iau_old_ana_in") )


    ! & p_diag_lnd%w_i(nproma,nblks_c)
    cf_desc     = t_cf_var('w_i', 'm H2O', 'weighted water content of interception water', &
         &                datatype_flt)
    new_cf_desc = t_cf_var('w_i', 'kg m-2', 'weighted water content of interception water', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 0, 13, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_i', p_diag_lnd%w_i,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in",       &
         &                 "mode_iau_fg_in","mode_iau_old_fg_in","mode_combined_in","mode_cosmode_in"), &
         & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )


    IF (itype_interception == 2) THEN

      ! & p_diag_lnd%w_p(nproma,nblks_c)
      cf_desc     = t_cf_var('w_p', 'm H2O', 'water content of pond interception water', &
           &                datatype_flt)
      new_cf_desc = t_cf_var('w_p', 'kg m-2', 'water content of pond interception water', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'w_p', p_diag_lnd%w_p,                &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
           & in_group=groups("land_vars"),                                         &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )

      ! & p_diag_lnd%w_s(nproma,nblks_c)
      cf_desc     = t_cf_var('w_s', 'm H2O', 'water content of interception snow', &
           &                datatype_flt)
      new_cf_desc = t_cf_var('w_s', 'kg m-2', 'water content of interception snow', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'w_s', p_diag_lnd%w_s,                &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
           & in_group=groups("land_vars"),                                         &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )

    END IF  ! itype_interception == 2



    ! & p_diag_lnd%t_so(nproma,nlev_soil+1,nblks_c)
    cf_desc    = t_cf_var('t_so', 'K', 'weighted soil temperature (main level)', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_so', p_diag_lnd%t_so,              &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND_P1, cf_desc, grib2_desc,  &
         & ldims=(/nproma,nlev_soil+1,kblks/),                                   &
         & lrestart=.FALSE., loutput=.TRUE.,                                     &
         & in_group=groups("land_vars","dwd_fg_sfc_vars",                        &
         &                 "mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in",&
         &                 "mode_dwd_ana_in","mode_iau_ana_in",                  &
         &                 "mode_iau_old_ana_in","mode_combined_in","mode_cosmode_in") )


    ! & p_diag_lnd%w_so(nproma,nlev_soil,nblks_c)
    cf_desc     = t_cf_var('w_so', 'm H20', 'total water content (ice + liquid water)', datatype_flt)
    new_cf_desc = t_cf_var('w_so', 'kg m-2', 'total water content (ice + liquid water)', datatype_flt)
    grib2_desc  = grib2_var(2, 3, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_so', p_diag_lnd%w_so,              &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc,     &
         & ldims=(/nproma,nlev_soil,kblks/),                                     &
         & lrestart=.FALSE., loutput=.TRUE.,                                     &
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_ana_in",      &
         &                 "mode_iau_fg_in","mode_iau_old_fg_in","mode_iau_ana_in", &
         &                 "mode_iau_old_ana_in","mode_combined_in",             &
         &                 "mode_cosmode_in"),                                   &
         & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ),&
         & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )

    ! & p_diag_lnd%w_so_ice(nproma,nlev_soil,nblks_c)
    cf_desc      = t_cf_var('w_so_ice', 'm H20',   'ice content', datatype_flt)
    new_cf_desc  = t_cf_var('w_so_ice', 'kg m-2', 'ice content', datatype_flt)
    grib2_desc   = grib2_var(2, 3, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)      
    CALL add_var( diag_list, vname_prefix//'w_so_ice',                           &
         & p_diag_lnd%w_so_ice, GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,     &
         & cf_desc, grib2_desc, ldims=(/nproma,nlev_soil,kblks/),                &
         & lrestart=.FALSE., loutput=.TRUE.,                                     &
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in",       &
         &                 "mode_iau_fg_in","mode_iau_old_fg_in"),               &
         & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ),&
         & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )


    ! & p_diag_lnd%runoff_s(nproma,nblks_c)
    cf_desc    = t_cf_var('runoff_s', 'kg m-2', &
         &                'weighted surface water runoff; sum over forecast', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_s', p_diag_lnd%runoff_s,        &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_S, cf_desc, grib2_desc,       &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
           & isteptype=TSTEP_ACCUM )


    ! & p_diag_lnd%runoff_g(nproma,nblks_c)
    cf_desc    = t_cf_var('runoff_g', 'kg m-2', &
         &                'weighted soil water runoff; sum over forecast', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_g', p_diag_lnd%runoff_g,        &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G, cf_desc, grib2_desc,       &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
           & isteptype=TSTEP_ACCUM )


    ! & p_diag_lnd%runoff_s_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('runoff_s_t', 'kg m-2', &
         &                'surface water runoff; sum over forecast', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_s_t', p_diag_lnd%runoff_s_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_S, cf_desc, grib2_desc,       &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container runoff_s
    ALLOCATE(p_diag_lnd%runoff_s_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_diag_lnd%runoff_s_ptr(jsfc)%p_2d, p_diag_lnd%runoff_s_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'runoff_s_t',                       &
               & vname_prefix//'runoff_s_t_'//ADJUSTL(TRIM(csfc)),               &
               & p_diag_lnd%runoff_s_ptr(jsfc)%p_2d,                             &
               & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_S,                      &
               & t_cf_var('runoff_s_t_'//csfc, '', '', datatype_flt),          &
               & grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
               & var_class=CLASS_TILE_LAND,                                      &
               & ldims=shape2d )
    END DO


    ! & p_diag_lnd%runoff_g_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('runoff_g_t', 'kg m-2', &
         &                'soil water runoff; sum over forecast', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_g_t', p_diag_lnd%runoff_g_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G, cf_desc, grib2_desc,       &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container runoff_g
    ALLOCATE(p_diag_lnd%runoff_g_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_diag_lnd%runoff_g_ptr(jsfc)%p_2d, p_diag_lnd%runoff_g_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'runoff_g_t',                       &
               & vname_prefix//'runoff_g_t_'//ADJUSTL(TRIM(csfc)),               &
               & p_diag_lnd%runoff_g_ptr(jsfc)%p_2d,                             &
               & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G,                      &
               & t_cf_var('runoff_g_t_'//csfc, '', '', datatype_flt),          &
               & grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
               & var_class=CLASS_TILE_LAND,                                      &
               & ldims=shape2d )
    END DO



    ! & p_diag_lnd%rstom(nproma,nblks_c)
    cf_desc    = t_cf_var('rstom', 's m-1','stomatal resistance', datatype_flt)
    grib2_desc = grib2_var(2, 0, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'rstom', p_diag_lnd%rstom,             &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE.,                                     &
           & loutput=.TRUE. )


    ! & p_diag_lnd%rstom_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('rstom_t', 's m-1','tile based stomatal resistance', &
      &          datatype_flt)
    grib2_desc = grib2_var(2, 0, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'rstom_t', p_diag_lnd%rstom_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,        &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,        &
           & loutput=.FALSE. )



    ! & p_diag_lnd%t_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('t_snow', 'K', 'weighted temperature of the snow-surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 18, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_snow', p_diag_lnd%t_snow,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,           &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                    &
         & in_group=groups("land_vars", "snow_vars","dwd_fg_sfc_vars",         &
         &                 "mode_dwd_ana_in","mode_iau_fg_in",                 &
         &                 "mode_iau_old_fg_in","mode_combined_in",            &
         &                 "mode_cosmode_in") )


    ! & p_diag_lnd%w_snow(nproma,nblks_c)
    cf_desc     = t_cf_var('w_snow', 'm H2O',   'weighted water eqivalent of snow', datatype_flt)
    new_cf_desc = t_cf_var('w_snow', 'kg m-2', 'weighted water equivalent of snow', datatype_flt)
    grib2_desc  = grib2_var(0, 1, 60, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_snow', p_diag_lnd%w_snow,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,           &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                    &
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in",     &
         &                 "mode_iau_old_ana_in",                              &
         &                 "mode_combined_in","mode_cosmode_in"),              &
         & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )


    ! & p_diag_lnd%rho_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('rho_snow', 'kg m-3', 'weighted snow density', datatype_flt)
    grib2_desc = grib2_var(0, 1, 61, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'rho_snow', p_diag_lnd%rho_snow,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,             &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
         & in_group=groups("land_vars", "snow_vars","dwd_fg_sfc_vars",           &
         &                 "mode_dwd_fg_in","mode_iau_fg_in",                    &
         &                 "mode_iau_old_ana_in","mode_combined_in",             &
         &                 "mode_cosmode_in") )


    ! & p_diag_lnd%h_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('h_snow', 'm', 'weighted snow depth', datatype_flt)
    grib2_desc = grib2_var(0, 1, 11, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'h_snow', p_diag_lnd%h_snow,        &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                  &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_fg_in", &
           &                 "mode_iau_ana_in","mode_iau_old_ana_in",              &
           &                 "mode_combined_in","mode_cosmode_in") )    


    ! & p_diag_lnd%h_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('h_snow_t', 'm', 'snow height', datatype_flt)
    grib2_desc = grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'h_snow_t', p_diag_lnd%h_snow_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )  

    ! fill the separate variables belonging to the container h_snow
    ALLOCATE(p_diag_lnd%h_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_diag_lnd%h_snow_ptr(jsfc)%p_2d, p_diag_lnd%h_snow_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'h_snow_t',                     &
               & vname_prefix//'h_snow_t_'//ADJUSTL(TRIM(csfc)),               &
               & p_diag_lnd%h_snow_ptr(jsfc)%p_2d,                             &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
               & t_cf_var('h_snow_t_'//csfc, '', '', datatype_flt),          &
               & grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
               & var_class=CLASS_TILE_LAND,                                    &
               & ldims=shape2d,                                                &
               & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
      ENDDO


    ! & p_diag_lnd%freshsnow(nproma,nblks_c)
    cf_desc    = t_cf_var('freshsnow', '1', &
           & 'weighted indicator for age of snow in top of snow layer', datatype_flt)
    grib2_desc = grib2_var(0, 1, 203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'freshsnow', p_diag_lnd%freshsnow,     &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                     &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_fg_in",&
           &                 "mode_iau_ana_in","mode_iau_old_ana_in",             &
           &                 "mode_combined_in","mode_cosmode_in") )


    ! & p_diag_lnd%freshsnow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('freshsnow_t', '1', &
         &                'indicator for age of snow in top of snow layer', datatype_flt)
    grib2_desc = grib2_var(0, 1, 203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'freshsnow_t', p_diag_lnd%freshsnow_t, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container freshsnow
    ALLOCATE(p_diag_lnd%freshsnow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_diag_lnd%freshsnow_ptr(jsfc)%p_2d, p_diag_lnd%freshsnow_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'freshsnow_t',                   &
               & vname_prefix//'freshsnow_t_'//ADJUSTL(TRIM(csfc)),           &
               & p_diag_lnd%freshsnow_ptr(jsfc)%p_2d,                         &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
               & t_cf_var('freshsnow_t_'//csfc, '', '', datatype_flt),      &
               & grib2_var(0, 1, 203, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
               & var_class=CLASS_TILE_LAND,                                   &
               & ldims=shape2d,                                               &
               & in_group=groups("dwd_fg_sfc_vars_t") )
    END DO


    ! & p_diag_lnd%snowfrac(nproma,nblks_c)
    cf_desc    = t_cf_var('snowfrac', '1' , 'snow-cover fraction', datatype_flt)
    new_cf_desc= t_cf_var('snowfrac', '% ', 'snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac', p_diag_lnd%snowfrac,       &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                     &
           & in_group=groups("land_vars"),                                        &
           & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc) )


    ! & p_diag_lnd%snowfrac_lc(nproma,nblks_c)
    cf_desc    = t_cf_var('snowfrac_lc', '1 ', 'snow-cover fraction', DATATYPE_FLT32)
    new_cf_desc= t_cf_var('snowfrac_lc', '% ', 'snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = grib2_var(0, 1, 42, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac_lc', p_diag_lnd%snowfrac_lc, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                     &
           & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_iau_fg_in"),     &
           & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc) )


    ! & p_diag_lnd%snowfrac_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('snowfrac_t', '1 ', 'local tile-based snow-cover fraction', datatype_flt)
    new_cf_desc= t_cf_var('snowfrac_t', '% ', 'local tile-based snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac_t', p_diag_lnd%snowfrac_t, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container snowfrac
    ALLOCATE(p_diag_lnd%snowfrac_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_diag_lnd%snowfrac_ptr(jsfc)%p_2d, p_diag_lnd%snowfrac_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'snowfrac_t',                     &
               & vname_prefix//'snowfrac_t_'//ADJUSTL(TRIM(csfc)),             &
               & p_diag_lnd%snowfrac_ptr(jsfc)%p_2d,                           &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
               & t_cf_var('snowfrac_t_'//csfc, '', '', datatype_flt),        &
               & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
               & var_class=CLASS_TILE_LAND,                                    &
               & ldims=shape2d,                                                &
               & in_group=groups("land_tile_vars"),                            &
               & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc) )
    END DO


    ! & p_diag_lnd%snowfrac_lc_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('snowfrac_lc_t', '1 ', 'tile-based snow-cover fraction', datatype_flt)
    new_cf_desc= t_cf_var('snowfrac_lc_t', '% ', 'tile-based snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = grib2_var(0, 1, 42, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac_lc_t', p_diag_lnd%snowfrac_lc_t, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ! fill the separate variables belonging to the container snowfrac
    ALLOCATE(p_diag_lnd%snowfrac_lc_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_diag_lnd%snowfrac_lc_ptr(jsfc)%p_2d, p_diag_lnd%snowfrac_lc_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'snowfrac_lc_t',                  &
               & vname_prefix//'snowfrac_lc_t_'//ADJUSTL(TRIM(csfc)),          &
               & p_diag_lnd%snowfrac_lc_ptr(jsfc)%p_2d,                        &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
               & t_cf_var('snowfrac_lc_t_'//csfc, '', '', datatype_flt),     &
               & grib2_var(0, 1, 42, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
               & var_class=CLASS_TILE_LAND, ldims=shape2d,                     &
               & in_group=groups("land_tile_vars", "dwd_fg_sfc_vars_t"),       &
               & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc) )
    END DO


    IF (lmulti_snow .OR. l2lay_rho_snow) THEN

      ! & p_diag_lnd%rho_snow_mult(nproma,nlev_snow,nblks_c)
      cf_desc    = t_cf_var('rho_snow_mult', 'kg m-3', 'weighted snow density', datatype_flt)
      grib2_desc = grib2_var(0, 1, 61, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'rho_snow_mult',                      &
           & p_diag_lnd%rho_snow_mult, GRID_UNSTRUCTURED_CELL, ZA_SNOW,            &
           & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
           & lrestart=.FALSE., loutput=.TRUE.,                                     &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",  &
           &                 "mode_iau_old_fg_in","multisnow_vars"))

    ENDIF

    IF (lmulti_snow) THEN

      ! & p_diag_lnd%t_snow_mult(nproma,nlev_snow+1,nblks_c)
      cf_desc    = t_cf_var('t_snow_mult', 'K', 'weighted temperature of the snow', datatype_flt)
      grib2_desc = grib2_var(0, 0, 18, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'t_snow_mult', p_diag_lnd%t_snow_mult, &
       & GRID_UNSTRUCTURED_CELL, ZA_SNOW_HALF, cf_desc, grib2_desc,                 &
       & ldims=(/nproma,nlev_snow+1,kblks/),                                        &
       & lrestart=.FALSE., loutput=.TRUE.,                                          &
       & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",       &
       &                 "mode_iau_old_fg_in","multisnow_vars") )


      ! & p_diag_lnd%wliq_snow(nproma,nlev_snow,nblks_c)
      cf_desc    = t_cf_var('wliq_snow', 'm H2O', 'weighted liquid water content in snow', &
        &                   datatype_flt)
      new_cf_desc= t_cf_var('wliq_snow', 'kg m-2', 'weighted liquid water content in snow', &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 210, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'wliq_snow',                          &
           & p_diag_lnd%wliq_snow, GRID_UNSTRUCTURED_CELL, ZA_SNOW,                & 
           & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
           & lrestart=.FALSE., loutput=.TRUE.,                                     &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),    &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",  &
           &                 "mode_iau_old_fg_in","multisnow_vars", "snow_vars") )


      ! & p_diag_lnd%wtot_snow(nproma,nlev_snow,nblks_c)
      cf_desc    = t_cf_var('wtot_snow', 'm H2O', 'weighted total water content in snow', &
           &                datatype_flt)
      new_cf_desc= t_cf_var('wtot_snow', 'kg m-2', 'weighted total water content in snow', &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 60, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'wtot_snow',                          &
           & p_diag_lnd%wtot_snow, GRID_UNSTRUCTURED_CELL, ZA_SNOW,                &
           & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
           & lrestart=.FALSE., loutput=.TRUE.,                                     &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),    & 
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",  &
           &                 "mode_iau_old_fg_in","multisnow_vars","snow_vars") )


      ! & p_diag_lnd%dzh_snow(nproma,nlev_snow,nblks_c)
      cf_desc    = t_cf_var('dzh_snow', 'm', &
           &                'weighted layer thickness between half levels in snow', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 1, 11, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'dzh_snow',                           &
           & p_diag_lnd%dzh_snow, GRID_UNSTRUCTURED_CELL, ZA_SNOW,                 &
           & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
           & lrestart=.FALSE., loutput=.TRUE.,                                     &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",  &
           &                 "mode_iau_old_fg_in","multisnow_vars","snow_vars"))

    ENDIF  ! lmulti_snow

    ENDIF  ! inwp_surface > 0

    

  END SUBROUTINE  new_nwp_lnd_diag_list


END MODULE mo_nwp_lnd_state

