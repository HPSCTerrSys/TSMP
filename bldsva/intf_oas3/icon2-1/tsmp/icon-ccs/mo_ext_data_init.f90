!>
!! Initialization/reading reading of external datasets
!!
!! This module contains read and initialization routines for the external data state.
!!
!! @author Daniel Reinert, DWD
!! @author Hermann Asensio, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-07-12)
!! Modification by Hermann Asensio, DWD (2010-07-16)
!!  - add miscellaneous variables for external parameters
!! Modification by Daniel Reinert, DWD (2011-05-03)
!! - Memory allocation method changed from explicit allocation to Luis'
!!   infrastructure
!! Modification by Daniel Reinert, DWD (2012-02-23)
!! - Routine smooth_topography moved to a new module named mo_smooth_topo
!! Modification by Daniel Reinert, DWD (2012-03-22)
!! - Type declaration moved to new module mo_ext_data_types
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ext_data_init

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_impl_constants,     ONLY: inwp, iecham, ildf_echam, io3_clim, io3_ape,                     &
    &                              ihs_atm_temp, ihs_atm_theta, inh_atmosphere,                     &
    &                              max_char_length, min_rlcell_int, min_rlcell,                     &
    &                              MODIS, GLOBCOVER2009, GLC2000, SUCCESS, SSTICE_ANA_CLINC,        &
    &                              SSTICE_CLIM
  USE mo_math_constants,     ONLY: dbl_eps, rad2deg
  USE mo_physical_constants, ONLY: ppmv2gg, zemiss_def, tmelt
  USE mo_run_config,         ONLY: msg_level, iforcing, check_uuid_gracefully
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_lnd_nwp_config,     ONLY: ntiles_total, ntiles_lnd, ntiles_water, lsnowtile, frlnd_thrhld, &
                                   frlndtile_thrhld, frlake_thrhld, frsea_thrhld, isub_water,       &
                                   isub_seaice, isub_lake, sstice_mode, sst_td_filename,            &
                                   ci_td_filename, itype_lndtbl, c_soil, c_soil_urb
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_extpar_config,      ONLY: itopo, l_emiss, extpar_filename, generate_filename, &
    &                              generate_td_filename, extpar_varnames_map_file, &
    &                              n_iter_smooth_topo, i_lctype, nclass_lu, nmonths_ext
  USE mo_dynamics_config,    ONLY: iequations
  USE mo_radiation_config,   ONLY: irad_o3, irad_aero, albedo_type
  USE mo_echam_phy_config,   ONLY: echam_phy_config
  USE mo_smooth_topo,        ONLY: smooth_topo_real_data
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_grid_config,        ONLY: n_dom, nroot
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work, my_process_is_mpi_workroot
  USE mo_sync,               ONLY: global_sum_array
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_ext_data_state,     ONLY: construct_ext_data, levelname, cellname, o3name, o3unit, &
    &                              nlev_o3, nmonths
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mo_time_config,        ONLY: time_config
  USE mo_io_config,          ONLY: default_read_method
  USE mo_read_interface,     ONLY: nf, openInputFile, closeFile, on_cells, &
    &                              t_stream_id, read_2D, read_2D_int, &
    &                              read_3D_extdim, read_2D_extdim
  USE mo_phyparam_soil,      ONLY: c_lnd, c_sea
  USE mo_util_cdi,           ONLY: get_cdi_varID, test_cdi_varID, read_cdi_2d,     &
    &                              read_cdi_3d, t_inputParameters,                 &
    &                              makeInputParameters, deleteInputParameters,     &
    &                              has_filetype_netcdf
  USE mo_util_uuid_types,    ONLY: t_uuid, uuid_string_length
  USE mo_util_uuid,          ONLY: OPERATOR(==), uuid_unparse
  USE mo_dictionary,         ONLY: t_dictionary, dict_init, dict_finalize,         &
    &                              dict_loadfile
  USE mo_initicon_config,    ONLY: timeshift
  USE mo_nwp_tuning_config,  ONLY: itune_albedo
  USE mo_master_config,      ONLY: isRestart
  USE mo_cdi,                ONLY: FILETYPE_GRB2, streamOpenRead, streamInqFileType, &
    &                              streamInqVlist, vlistInqVarZaxis, zaxisInqSize,   &
    &                              vlistNtsteps, vlistInqVarGrid, vlistInqAttTxt,    &
    &                              vlistInqVarIntKey, CDI_GLOBAL, gridInqUUID, &
    &                              streamClose, cdiStringError
  USE mo_math_gradients,     ONLY: grad_fe_cell
  USE mo_fortran_tools,      ONLY: var_scale, var_add
  USE mtime,                 ONLY: datetime, newDatetime, deallocateDatetime,        &
    &                              MAX_DATETIME_STR_LEN, datetimetostring,           &
    &                              OPERATOR(+)
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights,         &
    &                                  calculate_time_interpolation_weights
  USE mo_coupling_config,    ONLY: is_coupled_run


  IMPLICIT NONE

  ! required for reading external data
  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_ext_data_init'

  ! Number of landcover classes provided by external parameter data
  ! Needs to be changed into a variable if landcover classifications
  ! with a different number of classes become available
  INTEGER, PARAMETER :: num_lcc = 23, n_param_lcc = 7

  LOGICAL, ALLOCATABLE :: is_frglac_in(:) !< checks whether the extpar file contains fr_glac
  LOGICAL :: read_netcdf_data             !< control variable if extpar data are in GRIB2 for NetCDF format


  PUBLIC :: init_ext_data
  PUBLIC :: init_index_lists
  PUBLIC :: interpol_monthly_mean
  PUBLIC :: diagnose_ext_aggr


!-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Init external data for atmosphere
  !!
  !! 1. Build data structure, including field lists and
  !!    memory allocation.
  !! 2. External data are read in from netCDF file or set analytically
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-16)
  !!
  SUBROUTINE init_ext_data (p_patch, p_int_state, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_int_state), INTENT(IN)        :: p_int_state(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)


    INTEGER              :: jg, ist
    INTEGER, ALLOCATABLE :: cdi_extpar_id(:)  !< CDI stream ID (for each domain)
    INTEGER, ALLOCATABLE :: cdi_filetype(:)   !< CDI filetype (for each domain)
    ! dictionary which maps internal variable names onto
    ! GRIB2 shortnames or NetCDF var names.
    TYPE (t_dictionary) :: extpar_varnames_dict

    TYPE(datetime), POINTER :: this_datetime
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':init_ext_data'

    !SBr, CHa
    CHARACTER(len = 500)  :: sval
    !-------------------------------------------------------------------------
    CALL message (TRIM(routine), 'Start')

    !-------------------------------------------------------------------------
    !  1.  inquire external files for their data structure
    !-------------------------------------------------------------------------

    ALLOCATE(is_frglac_in(n_dom))
    ! Set default value for is_frglac_in. Will be overwritten, if external data
    ! contain fr_glac
    is_frglac_in(1:n_dom) = .FALSE.

    ! Allocate and open CDI stream (files):
    ALLOCATE (cdi_extpar_id(n_dom), cdi_filetype(n_dom), stat=ist)
    IF (ist /= SUCCESS)  CALL finish(TRIM(routine),'ALLOCATE failed!')
    CALL inquire_external_files(p_patch, cdi_extpar_id, cdi_filetype)

    ! read the map file (internal -> GRIB2) into dictionary data structure:
    CALL dict_init(extpar_varnames_dict, lcase_sensitive=.FALSE.)
    IF (ANY(cdi_filetype(:) == FILETYPE_GRB2)) THEN
      IF(extpar_varnames_map_file /= ' ') THEN
        CALL dict_loadfile(extpar_varnames_dict, TRIM(extpar_varnames_map_file))
      END IF
      read_netcdf_data = .FALSE.
    ELSE
      read_netcdf_data = .TRUE.
    END IF

    !------------------------------------------------------------------
    !  2.  construct external fields for the model
    !------------------------------------------------------------------

    ! SBr, CHa: nclass_lu is zero otherwise
    DO jg = 1, n_dom
       nclass_lu(jg) = 20
    END DO

    ! top-level procedure for building data structures for
    ! external data.
    CALL construct_ext_data(p_patch, ext_data)

    !-------------------------------------------------------------------------
    !  3.  read the data into the fields
    !-------------------------------------------------------------------------

    ! Check, whether external data should be read from file

    SELECT CASE(itopo)

    CASE(0) ! itopo, do not read external data except in some cases (see below)
      !
      ! initalize external data with meaningful data, in the case that they
      ! are not read in from file.
      IF ( iforcing == inwp ) THEN
        DO jg = 1, n_dom
          !SBr, CHa: comment  out
          !ext_data(jg)%atm%fr_land(:,:)     = 0._wp      ! land fraction
          !ext_data(jg)%atm%fr_land_smt(:,:) = 0._wp      ! land fraction (smoothed)
          !ext_data(jg)%atm%fr_glac_smt(:,:) = 0._wp      ! glacier fraction (smoothed)
          !ext_data(jg)%atm%llsm_atm_c(:,:)  = .FALSE.    ! land-sea mask
          !ext_data(jg)%atm%llake_c(:,:)     = .FALSE.    ! lake mask
          !ext_data(jg)%atm%plcov_mx(:,:)    = 0.5_wp     ! plant cover
          !ext_data(jg)%atm%lai_mx(:,:)      = 3._wp      ! max Leaf area index
          !ext_data(jg)%atm%rootdp(:,:)      = 1._wp      ! root depth
          !ext_data(jg)%atm%rsmin(:,:)       = 150._wp    ! minimal stomata resistence
          !ext_data(jg)%atm%soiltyp(:,:)     = 8          ! soil type
          !ext_data(jg)%atm%z0(:,:)          = 0.001_wp   ! roughness length

          !SBr, CHa : change fr_land to 1 from 0 also fr_land_smt  
          ext_data(jg)%atm%fr_land(:,:)     = 1._wp      ! land fraction
          ext_data(jg)%atm%fr_land_smt(:,:) = 1._wp      ! land fraction (smoothed)
          ext_data(jg)%atm%fr_glac_smt(:,:) = 0._wp      ! glacier fraction (smoothed)
          ext_data(jg)%atm%llsm_atm_c(:,:)  = .FALSE.    ! land-sea mask
          ext_data(jg)%atm%llake_c(:,:)     = .FALSE.    ! lake mask
          ext_data(jg)%atm%plcov_mx(:,:)    = 0.0_wp     ! plant cover
          ext_data(jg)%atm%lai_mx(:,:)      = 3._wp      ! max Leaf area index
          ext_data(jg)%atm%rootdp(:,:)      = 1._wp      ! root depth
          ext_data(jg)%atm%rsmin(:,:)       = 150._wp    ! minimal stomata resistence
          ext_data(jg)%atm%soiltyp(:,:)     = 5          ! soil type 5=loam, 6=loam clay
          ext_data(jg)%atm%z0(:,:)          = 0.05_wp    ! roughness length


          !Special setup for EDMF
          ext_data(jg)%atm%soiltyp_t(:,:,:) = 8           ! soil type
          ext_data(jg)%atm%frac_t(:,:,:)    = 0._wp       ! set all tiles to 0
          ext_data(jg)%atm%frac_t(:,:,isub_water) = 1._wp ! set only ocean to 1
          ext_data(jg)%atm%lc_class_t(:,:,:) = 1          ! land cover class

          ! SBr, CHa
          !ext_data(jg)%atm%fr_sand(:,:)       = 41._wp
          !ext_data(jg)%atm%fr_silt            = 21._wp
          !ext_data(jg)%atm%fr_clay            = 38._wp
          !ext_data(jg)%atm%fr_oc              = 1.7_wp
          !ext_data(jg)%atm%bulk_dens          = 1.3_wp
          ext_data(jg)%atm_td%fr_ice_m        = 0._wp
          ext_data(jg)%atm%i_lc_urban         = 0._wp
          ext_data(jg)%atm%for_d              = 0.3_wp
          ext_data(jg)%atm%for_e              = 0.09_wp
          ext_data(jg)%atm%emis_rad           = 0.96_wp  !SBr, CHa: change from 0.99 to 0.96
          ext_data(jg)%atm%ndvi_max           = 0.9_wp
          ext_data(jg)%atm%lu_class_fraction  = 0  
          ext_data(jg)%atm%t_cl               = 282._wp
          ext_data(jg)%atm%sso_stdh           = 15._wp
          ext_data(jg)%atm%sso_theta          = 0.05_wp
          ext_data(jg)%atm%sso_gamma          = 0.42_wp
          ext_data(jg)%atm%sso_sigma          = 0.1_wp
          ext_data(jg)%atm%depth_lk           = -0.9_wp

          ext_data(jg)%atm%topography_c     = 10._wp
          ext_data(jg)%atm%fis              = 10._wp * 9.81_wp
          ext_data(jg)%atm%grad_topo        = 0._wp
          ext_data(jg)%atm%llsm_atm_c       = .TRUE.
          ext_data(jg)%atm%llake_c          = .FALSE.
          ext_data(jg)%atm%fr_lake          = 0._wp
          ext_data(jg)%atm%sso_stdh         = 0._wp
          ext_data(jg)%atm%z0_lcc(:)        = 0.05_wp    ! land cove related roughness length
          ext_data(jg)%atm%z0_lcc_min(:)    = 0.0002_wp    ! land cove related roughness length
        END DO
      END IF
      ! SBr, CHa: comment out the default emis_rad value: zemiss_def
      !IF ( iforcing == inwp .OR. iforcing == iecham .OR. iforcing == ildf_echam ) THEN
      !  DO jg = 1, n_dom
      !    ext_data(jg)%atm%emis_rad(:,:)    = zemiss_def ! longwave surface emissivity
      !  END DO
      !END IF

      ! call read_ext_data_atm to read O3
      ! topography is used from analytical functions, except for ljsbach=.TRUE. in which case
      ! elevation of cell centers is read in and the topography is "grown" gradually to this elevation
      IF ( irad_o3 == io3_clim .OR. irad_o3 == io3_ape .OR. sstice_mode == SSTICE_CLIM .OR. &
         & echam_phy_config%ljsbach) THEN
        IF ( echam_phy_config%ljsbach .AND. (iequations /= inh_atmosphere) ) THEN
          CALL message( TRIM(routine),'topography is grown to elevation' )
        ELSE
          CALL message( TRIM(routine),'Running with analytical topography' )
        END IF
        CALL read_ext_data_atm (p_patch, ext_data, nlev_o3, cdi_extpar_id, &
          &                     extpar_varnames_dict)
        CALL message( TRIM(routine),'read_ext_data_atm completed' )
      END IF

    CASE(1) ! itopo, read external data from file

      CALL message( TRIM(routine),'Start reading external data from file' )

      CALL read_ext_data_atm (p_patch, ext_data, nlev_o3, cdi_extpar_id, &
        &                     extpar_varnames_dict)

      CALL message( TRIM(routine),'Finished reading external data' )

      IF ( iforcing == inwp ) THEN

        DO jg = 1, n_dom
         ! topography smoothing
         IF (n_iter_smooth_topo(jg) > 0) THEN
            CALL smooth_topo_real_data ( p_patch(jg)                   ,&
              &                          p_int_state(jg)               ,&
              &                          ext_data(jg)%atm%fr_land      ,&
              &                          ext_data(jg)%atm%fr_lake      ,&
              &                          ext_data(jg)%atm%topography_c ,&
              &                          ext_data(jg)%atm%sso_stdh     )
          ENDIF

          ! calculate gradient of orography for resolved surface drag
          !
          call grad_fe_cell  ( ext_data(jg)%atm%topography_c, &
            &                  p_patch(jg),                   &
            &                  p_int_state(jg),               &
            &                  ext_data(jg)%atm%grad_topo )
        END DO



        ! Get interpolated ndviratio, alb_dif, albuv_dif and albni_dif. Interpolation
        ! is done in time, based on ini_datetime (midnight). Fields are updated on a
        ! daily basis.

        ! When initializing the model we set the target hour to 0 (midnight) as well.
        ! When restarting, the target interpolation time must be set to cur_datetime
        ! midnight.
        !
        IF (.NOT. isRestart()) THEN
          this_datetime => newDatetime(time_config%tc_startdate)
          IF (timeshift%dt_shift < 0._wp) THEN
            this_datetime = this_datetime + timeshift%mtime_shift
          END IF
        ELSE
          this_datetime => newDatetime(time_config%tc_current_date)
        END IF  ! isRestart
        
        ! always assume midnight
        this_datetime%time%hour   = 0   
        this_datetime%time%minute = 0
        this_datetime%time%second = 0
        this_datetime%time%ms     = 0

        DO jg = 1, n_dom
          CALL interpol_monthly_mean(p_patch(jg), this_datetime,         &! in
            &                        ext_data(jg)%atm_td%ndvi_mrat,      &! in
            &                        ext_data(jg)%atm%ndviratio          )! out
        ENDDO

        IF ( albedo_type == MODIS) THEN
          DO jg = 1, n_dom
            CALL interpol_monthly_mean(p_patch(jg), this_datetime,       &! in
              &                        ext_data(jg)%atm_td%alb_dif,      &! in
              &                        ext_data(jg)%atm%alb_dif          )! out

            CALL interpol_monthly_mean(p_patch(jg), this_datetime,       &! in
              &                        ext_data(jg)%atm_td%albuv_dif,    &! in
              &                        ext_data(jg)%atm%albuv_dif        )! out

            CALL interpol_monthly_mean(p_patch(jg), this_datetime,       &! in
              &                        ext_data(jg)%atm_td%albni_dif,    &! in
              &                        ext_data(jg)%atm%albni_dif        )! out
          ENDDO
        ENDIF  ! albedo_type

        ! clean up
        CALL deallocateDatetime(this_datetime)

      END IF

    CASE DEFAULT ! itopo

      CALL finish( TRIM(routine), 'topography selection not supported' )

    END SELECT ! itopo

    ! close CDI stream (file):
    DO jg=1,n_dom
      IF (cdi_extpar_id(jg) == -1) CYCLE
      IF (my_process_is_mpi_workroot())  CALL streamClose(cdi_extpar_id(jg))
    END DO
    DEALLOCATE (cdi_extpar_id, cdi_filetype, stat=ist)
    IF (ist /= SUCCESS)  CALL finish(TRIM(routine),'DEALLOCATE failed!')

    ! destroy variable name dictionary:
    CALL dict_finalize(extpar_varnames_dict)

  END SUBROUTINE init_ext_data


  !-------------------------------------------------------------------------
  ! Open ExtPar file and investigate the data structure of the
  ! external parameters.
  !
  ! Note: This subroutine opens the file and returns a CDI file ID.
  !
  ! @author F. Prill, DWD (2014-01-07)
  !-------------------------------------------------------------------------
  SUBROUTINE inquire_extpar_file(p_patch, jg, cdi_extpar_id, cdi_filetype, &
    &                            is_frglac_in)
    TYPE(t_patch), INTENT(IN)      :: p_patch(:)
    INTEGER,       INTENT(IN)      :: jg
    INTEGER,       INTENT(INOUT)   :: cdi_extpar_id     !< CDI stream ID
    INTEGER,       INTENT(INOUT)   :: cdi_filetype      !< CDI filetype
    LOGICAL,       INTENT(OUT)     :: is_frglac_in      !< check for fr_glac in Extpar file

    ! local variables
    CHARACTER(len=max_char_length), PARAMETER :: routine = modname//'::inquire_extpar_file'
    INTEGER                 :: mpi_comm, vlist_id, lu_class_fraction_id, zaxis_id, var_id
    LOGICAL                 :: l_exist
    CHARACTER(filename_max) :: extpar_file !< file name for reading in
    INTEGER :: extpar_file_namelen

    TYPE(t_uuid)            :: extpar_uuidOfHGrid             ! uuidOfHGrid contained in the
                                                              ! extpar file

    CHARACTER(len=uuid_string_length) :: grid_uuid_unparsed   ! unparsed grid uuid (human readable)
    CHARACTER(len=uuid_string_length) :: extpar_uuid_unparsed ! same for extpar-file uuid

    LOGICAL                 :: lmatch                         ! for comparing UUIDs

    INTEGER                 :: cdiGridID

    INTEGER :: lu_var_id, localInformationNumber
    INTEGER :: ret
    CHARACTER(len=max_char_length) :: rawdata_attr

    !---------------------------------------------!
    ! Check validity of external parameter file   !
    !---------------------------------------------!
    IF (my_process_is_mpi_workroot()) THEN
      ! generate file name
      extpar_file = generate_filename(extpar_filename,                   &
        &                             getModelBaseDir(),                 &
        &                             TRIM(p_patch(jg)%grid_filename),   &
        &                             nroot,                             &
        &                             p_patch(jg)%level, p_patch(jg)%id)
      extpar_file_namelen = LEN_TRIM(extpar_file)
      CALL message(routine, "extpar_file = "//extpar_file(1:extpar_file_namelen))

      INQUIRE (FILE=extpar_file, EXIST=l_exist)
      IF (.NOT.l_exist)  CALL finish(routine,'external data file is not found.')

      ! open file
      cdi_extpar_id = streamOpenRead(extpar_file(1:extpar_file_namelen))
      IF (cdi_extpar_id < 0) THEN
        WRITE (message_text, '(133a)') "Cannot open external parameter file ", &
             cdiStringError(cdi_extpar_id)
        CALL finish(routine, TRIM(message_text))
      END IF
      cdi_filetype  = streamInqFileType(cdi_extpar_id)


      ! get the number of landuse classes
      lu_class_fraction_id = get_cdi_varID(cdi_extpar_id, "LU_CLASS_FRACTION")
      vlist_id             = streamInqVlist(cdi_extpar_id)
      zaxis_id             = vlistInqVarZaxis(vlist_id, lu_class_fraction_id)
      nclass_lu(jg)        = zaxisInqSize(zaxis_id)

      ! get time dimension from external data file
      nmonths_ext(jg)      = vlistNtsteps(vlist_id)

      ! make sure that num_lcc is equal to nclass_lu. If not, then the internal
      ! land-use lookup tables and the external land-use class field are inconsistent.

      IF (nclass_lu(jg) /= num_lcc) THEN
        WRITE(message_text,'(A,I3,A,I3)')  &
          & 'Number of land-use classes in external file ', nclass_lu(jg), &
          & ' does not match ICON-internal value num_lcc ', num_lcc
        CALL finish(routine,TRIM(message_text))
      ENDIF

      IF ( msg_level>10 ) THEN
        WRITE(message_text,'(A,I4)')  &
          & 'Number of land_use classes in external data file = ', nclass_lu(jg)
        CALL message(routine,message_text)

        WRITE(message_text,'(A,I4)')  &
          & 'Number of months in external data file = ', nmonths_ext(jg)
        CALL message(routine,message_text)
      ENDIF


      ! Compare UUID of external parameter file with UUID of grid.
      !
      ! get horizontal grid UUID contained in extpar file
      ! use lu_class_fraction as sample field
      cdiGridID = vlistInqVarGrid(vlist_id, lu_class_fraction_id)
      CALL gridInqUUID(cdiGridID, extpar_uuidOfHGrid%DATA)
      !
      ! --- compare UUID of horizontal grid file with UUID from extpar file
      lmatch = (p_patch(jg)%grid_uuid == extpar_uuidOfHGrid)

      IF (.NOT. lmatch) THEN
        CALL uuid_unparse(p_patch(jg)%grid_uuid, grid_uuid_unparsed)
        CALL uuid_unparse(extpar_uuidOfHGrid   , extpar_uuid_unparsed)
        WRITE(message_text,'(a,a)') 'uuidOfHgrid from gridfile: ', TRIM(grid_uuid_unparsed)
        CALL message(routine,message_text)
        WRITE(message_text,'(a,a)') 'uuidOfHgrid from extpar file: ', TRIM(extpar_uuid_unparsed)
        CALL message(routine,message_text)

        WRITE(message_text,'(a)') 'Extpar file and horizontal grid file do not match!'
        IF (check_uuid_gracefully) THEN
          CALL message(routine, TRIM(message_text))
        ELSE
          CALL finish(routine, TRIM(message_text))
        END IF
      ENDIF



      ! Determine which data source has been used to generate the
      ! external perameters: For NetCDF format, we check the
      ! global attribute "rawdata". For GRIB2 format we check the
      ! key "localInformationNumber".
      IF (has_filetype_netcdf(cdi_filetype)) THEN
        ret      = vlistInqAttTxt(vlist_id, CDI_GLOBAL, 'rawdata', max_char_length, rawdata_attr)
        IF (INDEX(rawdata_attr,'GLC2000') /= 0) THEN
          i_lctype(jg) = GLC2000
        ELSE IF (INDEX(rawdata_attr,'GLOBCOVER2009') /= 0) THEN
          i_lctype(jg) = GLOBCOVER2009
        ELSE
          CALL finish(routine,'Unknown landcover data source')
        ENDIF
      ELSE IF (cdi_filetype == FILETYPE_GRB2) THEN
        lu_var_id              = get_cdi_varID(cdi_extpar_id, 'LU_CLASS_FRACTION')
        localInformationNumber = vlistInqVarIntKey(vlist_id, lu_var_id, "localInformationNumber")
        SELECT CASE (localInformationNumber)
        CASE (2)  ! 2 = GLC2000
          i_lctype(jg) = GLC2000
        CASE (1)  ! 1 = ESA GLOBCOVER
          i_lctype(jg) = GLOBCOVER2009
        CASE DEFAULT
          CALL finish(routine,'Unknown landcover data source')
        END SELECT
      END IF

      ! Check whether external parameter file contains MODIS albedo-data
      IF ( albedo_type == MODIS ) THEN
        IF ( (test_cdi_varID(cdi_extpar_id, 'ALB')   == -1) .OR.    &
          &  (test_cdi_varID(cdi_extpar_id, 'ALNID') == -1) .OR.    &
          &  (test_cdi_varID(cdi_extpar_id, 'ALUVD') == -1) ) THEN
          CALL finish(routine,'MODIS albedo fields missing in '//TRIM(extpar_filename))
        ENDIF
      ENDIF

      ! Check whether external parameter file contains SST climatology
      IF ( sstice_mode == SSTICE_ANA_CLINC ) THEN
        IF ( test_cdi_varID(cdi_extpar_id, 'T_SEA')  == -1 ) THEN
          CALL finish(routine,'SST climatology missing in '//TRIM(extpar_filename))
        ENDIF
      ENDIF

      ! Search for glacier fraction in Extpar file
      !
      IF (has_filetype_netcdf(cdi_filetype)) THEN
        var_id = test_cdi_varID(cdi_extpar_id,'ICE')
      ELSE IF (cdi_filetype == FILETYPE_GRB2) THEN
        var_id = test_cdi_varID(cdi_extpar_id,'FR_ICE')
      ENDIF
      IF (var_id == -1) THEN
        is_frglac_in = .FALSE.
      ELSE
        is_frglac_in = .TRUE.
      ENDIF

    ENDIF ! my_process_is_mpi_workroot()

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF
    ! broadcast nclass_lu from I-Pe to WORK Pes
    CALL p_bcast(nclass_lu(jg), p_io, mpi_comm)
    ! broadcast nmonths from I-Pe to WORK Pes
    CALL p_bcast(nmonths_ext(jg), p_io, mpi_comm)
    ! broadcast is_frglac_in from I-Pe to WORK Pes
    CALL p_bcast(is_frglac_in, p_io, mpi_comm)
    ! broadcast i_lctype from I-Pe to WORK Pes
    CALL p_bcast(i_lctype(jg), p_io, mpi_comm)
    ! broadcast cdi filetype
    CALL p_bcast(cdi_filetype, p_io, mpi_comm)

  END SUBROUTINE inquire_extpar_file

  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE inquire_external_files(p_patch, cdi_extpar_id, cdi_filetype)

    !-------------------------------------------------------
    !
    ! open netcdf files and investigate the data structure
    ! of the external parameters
    !
    !-------------------------------------------------------

    TYPE(t_patch), INTENT(IN)      :: p_patch(:)
    INTEGER,       INTENT(INOUT)   :: cdi_extpar_id(:)  !< CDI stream ID
    INTEGER,       INTENT(INOUT)   :: cdi_filetype(:)   !< CDI filetype

    INTEGER :: jg, mpi_comm
    INTEGER :: no_cells
    INTEGER :: ncid, dimid

    LOGICAL :: l_exist

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':inquire_external_files'

    CHARACTER(filename_max) :: ozone_file  !< file name for reading in

!--------------------------------------------------------------------------

    ! set stream IDs to "uninitialized":
    IF(my_process_is_mpi_workroot()) THEN
      cdi_extpar_id(:) = -1
    END IF

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    DO jg= 1,n_dom

      !------------------------------------------------!
      ! 1. Check validity of external parameter file   !
      !------------------------------------------------!

      IF ( itopo == 1 .AND. iforcing == inwp) THEN
        CALL inquire_extpar_file(p_patch, jg, cdi_extpar_id(jg), cdi_filetype(jg), &
          &                      is_frglac_in(jg))
      END IF

      !------------------------------------------------!
      ! 2. Check validity of ozone file                !
      !------------------------------------------------!

      ! default values for nlev_o3 and nmonths
      nlev_o3 = 1
      nmonths   = 1

      O3 : IF ((irad_o3 == io3_clim) .OR. (irad_o3 == io3_ape )) THEN

        IF(irad_o3 == io3_ape ) THEN
          levelname = 'level'
          cellname  = 'ncells'
          o3name    = 'O3'
          o3unit    = 'g/g'
        ELSE ! o3_clim
          levelname = 'plev'
          cellname  = 'ncells'
          o3name    = 'O3'
          o3unit    = 'g/g' !this unit ozon will have after being read out and converted from ppmv
        ENDIF

        IF_IO : IF(my_process_is_stdio()) THEN

          WRITE(ozone_file,'(a,i2.2,a)') 'o3_icon_DOM',jg,'.nc'

          ! Note resolution assignment is done per script by symbolic links

          INQUIRE (FILE=ozone_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            WRITE(0,*) 'DOMAIN=',jg
            CALL finish(TRIM(routine),'ozone file of domain is not found.')
          ENDIF

          !
          ! open file
          !
          CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid), routine)
          WRITE(0,*)'open ozone file'

          !
          ! get number of cells
          !
          CALL nf(nf_inq_dimid (ncid, TRIM(cellname), dimid), routine)
          CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)
          WRITE(0,*)'number of cells are', no_cells

          !
          ! check the number of cells and verts
          !
          IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
            CALL finish(TRIM(ROUTINE),&
              & 'Number of patch cells and cells in ozone file do not match.')
          ENDIF

          !
          ! check the time structure
          !
          CALL nf(nf_inq_dimid (ncid, 'time', dimid), routine)
          CALL nf(nf_inq_dimlen(ncid, dimid, nmonths), routine)
          WRITE(message_text,'(A,I4)')  &
            & 'Number of months in ozone file = ', nmonths
          CALL message(TRIM(ROUTINE),message_text)

          !
          ! check the vertical structure
          !
          CALL nf(nf_inq_dimid (ncid,TRIM(levelname), dimid), routine)
          CALL nf(nf_inq_dimlen(ncid, dimid, nlev_o3), routine)

          WRITE(message_text,'(A,I4)')  &
            & 'Number of pressure levels in ozone file = ', nlev_o3
          CALL message(TRIM(ROUTINE),message_text)

          !
          ! close file
          !
          CALL nf(nf_close(ncid), routine)

        END IF IF_IO ! pe

        CALL p_bcast(nlev_o3, p_io, mpi_comm)
        CALL p_bcast(nmonths,   p_io, mpi_comm)

      END IF O3 !o3

    ENDDO ! ndom

  END SUBROUTINE inquire_external_files

  !-------------------------------------------------------------------------
  !>
  !! Read atmospheric external data
  !!
  !! Read atmospheric external data from netcdf
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !!
  SUBROUTINE read_ext_data_atm (p_patch, ext_data, nlev_o3, cdi_extpar_id, &
    &                           extpar_varnames_dict)

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    INTEGER,               INTENT(IN)    :: nlev_o3

    INTEGER,               INTENT(IN)    :: cdi_extpar_id(:)      !< CDI stream ID
    TYPE (t_dictionary),   INTENT(IN)    :: extpar_varnames_dict  !< variable names dictionary (for GRIB2)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':read_ext_data_atm'

    CHARACTER(filename_max) :: ozone_file  !< file name for reading in
    CHARACTER(filename_max) :: sst_td_file !< file name for reading in
    CHARACTER(filename_max) :: ci_td_file  !< file name for reading in

    INTEGER :: jg, jc, jb, i, mpi_comm, ilu,im
    INTEGER :: jk
    INTEGER :: ncid, varid
    TYPE(t_stream_id) :: stream_id

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk   !> blocks
    INTEGER :: i_startidx, i_endidx   !< slices
    INTEGER :: i_nchdom               !< domain index

    REAL(wp):: zdummy_o3lev(nlev_o3) ! will be used for pressure and height levels
    REAL(wp):: albfac, albthresh             ! for MODIS albedo tuning

    REAL(wp), DIMENSION(num_lcc*n_param_lcc)         :: lu_glc2000   ! < lookup table landuse class GLC2000
    REAL(wp), DIMENSION(num_lcc*n_param_lcc), TARGET :: lu_gcv2009   ! < lookup table landuse class GlobCover2009
    REAL(wp), DIMENSION(num_lcc*n_param_lcc), TARGET :: lu_gcv2009_v2 ! < modified lookup table landuse class GlobCover2009
    REAL(wp), DIMENSION(num_lcc*n_param_lcc), TARGET :: lu_gcv2009_v3 ! < even less evaporating lookup table landuse class GlobCover2009
    REAL(wp), DIMENSION(num_lcc*n_param_lcc), TARGET :: lu_gcv2009_v4 ! < retuned lookup table landuse class GlobCover2009
    REAL(wp), POINTER :: lu_gcv(:)

    LOGICAL :: l_exist
    CHARACTER(filename_max) :: extpar_file

    TYPE(t_inputParameters) :: parameters

!                    z0         pcmx      laimx rd      rsmin      snowalb snowtile
!
 DATA lu_glc2000 /   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp,-1._wp, & ! evergreen broadleaf forest
                 &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 150.0_wp,  0.31_wp,-1._wp, & ! deciduous broadleaf closed forest
                 &   0.15_wp,  0.8_wp,  4.0_wp, 2.0_wp, 150.0_wp,  0.31_wp,-1._wp, & ! deciduous broadleaf open   forest
                 &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.27_wp,-1._wp, & ! evergreen needleleaf forest
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.33_wp,-1._wp, & ! deciduous needleleaf forest
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 150.0_wp,  0.29_wp,-1._wp, & ! mixed leaf trees
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  -1.0_wp,-1._wp, & ! fresh water flooded trees
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  -1.0_wp,-1._wp, & ! saline water flooded trees
                 &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  -1.0_wp, 1._wp, & ! mosaic tree / natural vegetation
                 &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp, 150.0_wp,  -1.0_wp, 1._wp, & ! burnt tree cover
                 &   0.20_wp,  0.8_wp,  3.0_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! evergreen shrubs closed-open
                 &   0.15_wp,  0.8_wp,  1.5_wp, 2.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! decidous shrubs closed-open
                 &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp,  40.0_wp,  -1.0_wp, 1._wp, & ! herbaceous vegetation closed-open
                 &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp,  40.0_wp,  -1.0_wp, 1._wp, & ! sparse herbaceous or grass
                 &   0.05_wp,  0.8_wp,  2.0_wp, 0.4_wp,  40.0_wp,  -1.0_wp,-1._wp, & ! flooded shrubs or herbaceous
                 &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! cultivated & managed areas
                 &   0.25_wp,  0.8_wp,  3.0_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! mosaic crop / tree / natural vegetation
                 &   0.07_wp,  0.9_wp,  3.5_wp, 1.0_wp, 100.0_wp,  -1.0_wp, 1._wp, & ! mosaic crop / shrub / grass
                 &   0.05_wp,  0.05_wp, 0.6_wp, 0.3_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! bare areas
                 &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! water
                 &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! snow & ice
                 &   1.00_wp,  0.2_wp,  1.0_wp, 0.6_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! artificial surface
                 &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp,  40.0_wp,  -1.0_wp,-1._wp / ! undefined

 DATA lu_gcv2009 /   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 120.0_wp,  0.72_wp, 1._wp, & ! irrigated croplands
                 &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 120.0_wp,  0.72_wp, 1._wp, & ! rainfed croplands
                 &   0.25_wp,  0.8_wp,  3.0_wp, 1.0_wp, 120.0_wp,  0.55_wp, 1._wp, & ! mosaic cropland (50-70%) - vegetation (20-50%)
                 &   0.07_wp,  0.9_wp,  3.5_wp, 1.0_wp, 100.0_wp,  0.72_wp, 1._wp, & ! mosaic vegetation (50-70%) - cropland (20-50%)
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp, 1._wp, & ! closed broadleaved evergreen forest
                 &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 150.0_wp,  0.31_wp, 1._wp, & ! closed broadleaved deciduous forest
                 &   0.15_wp,  0.8_wp,  4.0_wp, 2.0_wp, 150.0_wp,  0.31_wp, 1._wp, & ! open broadleaved deciduous forest
                 &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.27_wp, 1._wp, & ! closed needleleaved evergreen forest
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.33_wp, 1._wp, & ! open needleleaved deciduous forest
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 150.0_wp,  0.29_wp, 1._wp, & ! mixed broadleaved and needleleaved forest
                 &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  0.60_wp, 1._wp, & ! mosaic shrubland (50-70%) - grassland (20-50%)
                 &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  0.65_wp, 1._wp, & ! mosaic grassland (50-70%) - shrubland (20-50%)
                 &   0.15_wp,  0.8_wp,  2.5_wp, 1.5_wp, 120.0_wp,  0.65_wp, 1._wp, & ! closed to open shrubland
                 &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp,  40.0_wp,  0.82_wp, 1._wp, & ! closed to open herbaceous vegetation
                 &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp,  40.0_wp,  0.76_wp, 1._wp, & ! sparse vegetation
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  0.30_wp, 1._wp, & ! closed to open forest regulary flooded
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  0.30_wp, 1._wp, & ! closed forest or shrubland permanently flooded
                 &   0.05_wp,  0.8_wp,  2.0_wp, 1.0_wp,  40.0_wp,  0.76_wp, 1._wp, & ! closed to open grassland regularly flooded
                 &   1.00_wp,  0.2_wp,  1.6_wp, 0.6_wp, 120.0_wp,  0.50_wp, 1._wp, & ! artificial surfaces
                 &   0.05_wp,  0.05_wp, 0.6_wp, 0.3_wp, 120.0_wp,  0.82_wp, 1._wp, & ! bare areas
                 &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp,-1._wp, & ! water bodies
                 &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! permanent snow and ice
                 &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp, 250.0_wp,  -1.0_wp,-1._wp  / !undefined

! Tuned version of gcv2009 based on IFS values (Juergen Helmert und Martin Koehler)
 DATA lu_gcv2009_v2 /  0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 180.0_wp,  0.72_wp, 1._wp, & ! irrigated croplands
                   &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 140.0_wp,  0.72_wp, 1._wp, & ! rainfed croplands
                   &   0.25_wp,  0.8_wp,  3.0_wp, 1.0_wp, 130.0_wp,  0.55_wp, 1._wp, & ! mosaic cropland (50-70%) - vegetation (20-50%)
                   &   0.07_wp,  0.9_wp,  3.5_wp, 1.0_wp, 120.0_wp,  0.72_wp, 1._wp, & ! mosaic vegetation (50-70%) - cropland (20-50%)
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp, 1._wp, & ! closed broadleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 175.0_wp,  0.31_wp, 1._wp, & ! closed broadleaved deciduous forest
                   &   0.15_wp,  0.8_wp,  4.0_wp, 1.5_wp, 175.0_wp,  0.31_wp, 1._wp, & ! open broadleaved deciduous forest
                   &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 250.0_wp,  0.27_wp, 1._wp, & ! closed needleleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 250.0_wp,  0.33_wp, 1._wp, & ! open needleleaved deciduous forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 210.0_wp,  0.29_wp, 1._wp, & ! mixed broadleaved and needleleaved forest
                   &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  0.60_wp, 1._wp, & ! mosaic shrubland (50-70%) - grassland (20-50%)
                   &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  0.65_wp, 1._wp, & ! mosaic grassland (50-70%) - shrubland (20-50%)
                   &   0.15_wp,  0.8_wp,  2.5_wp, 1.5_wp, 225.0_wp,  0.65_wp, 1._wp, & ! closed to open shrubland
                   &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp, 100.0_wp,  0.82_wp, 1._wp, & ! closed to open herbaceous vegetation
                   &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp,  80.0_wp,  0.76_wp, 1._wp, & ! sparse vegetation
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  0.30_wp, 1._wp, & ! closed to open forest regulary flooded
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  0.30_wp, 1._wp, & ! closed forest or shrubland permanently flooded
                   &   0.05_wp,  0.8_wp,  2.0_wp, 1.0_wp,  80.0_wp,  0.76_wp, 1._wp, & ! closed to open grassland regularly flooded
                   &   1.00_wp,  0.2_wp,  1.6_wp, 0.6_wp, 180.0_wp,  0.50_wp, 1._wp, & ! artificial surfaces
                   &   0.05_wp,  0.05_wp, 0.6_wp, 0.3_wp, 200.0_wp,  0.82_wp, 1._wp, & ! bare areas
                   &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 150.0_wp,  -1.0_wp,-1._wp, & ! water bodies
                   &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! permanent snow and ice
                   &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp, 250.0_wp,  -1.0_wp,-1._wp  / ! undefined

! Even more tuned version of gcv2009 by Guenther Zaengl (appears to produce the smallest temperature biases)
 DATA lu_gcv2009_v3 /  0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 190.0_wp,  0.72_wp, 1._wp, & ! irrigated croplands
                   &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 170.0_wp,  0.72_wp, 1._wp, & ! rainfed croplands
                   &   0.25_wp,  0.8_wp,  3.0_wp, 0.5_wp, 160.0_wp,  0.55_wp, 1._wp, & ! mosaic cropland (50-70%) - vegetation (20-50%)
                   &   0.07_wp,  0.9_wp,  3.5_wp, 0.7_wp, 150.0_wp,  0.72_wp, 1._wp, & ! mosaic vegetation (50-70%) - cropland (20-50%)
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 280.0_wp,  0.38_wp, 1._wp, & ! closed broadleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 225.0_wp,  0.31_wp, 1._wp, & ! closed broadleaved deciduous forest
                   &   0.15_wp,  0.8_wp,  4.0_wp, 1.5_wp, 225.0_wp,  0.31_wp, 1._wp, & ! open broadleaved deciduous forest
                   &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 300.0_wp,  0.27_wp, 1._wp, & ! closed needleleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 300.0_wp,  0.33_wp, 1._wp, & ! open needleleaved deciduous forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 270.0_wp,  0.29_wp, 1._wp, & ! mixed broadleaved and needleleaved forest
                   &   0.20_wp,  0.8_wp,  2.5_wp, 0.8_wp, 200.0_wp,  0.60_wp, 1._wp, & ! mosaic shrubland (50-70%) - grassland (20-50%)
                   &   0.20_wp,  0.8_wp,  2.5_wp, 0.6_wp, 200.0_wp,  0.65_wp, 1._wp, & ! mosaic grassland (50-70%) - shrubland (20-50%)
                   &   0.15_wp,  0.8_wp,  2.5_wp, 0.9_wp, 265.0_wp,  0.65_wp, 1._wp, & ! closed to open shrubland
                   &   0.03_wp,  0.9_wp,  3.1_wp, 0.4_wp, 140.0_wp,  0.82_wp, 1._wp, & ! closed to open herbaceous vegetation
                   &   0.05_wp,  0.5_wp,  0.6_wp, 0.2_wp, 120.0_wp,  0.76_wp, 1._wp, & ! sparse vegetation
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 190.0_wp,  0.30_wp, 1._wp, & ! closed to open forest regulary flooded
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 190.0_wp,  0.30_wp, 1._wp, & ! closed forest or shrubland permanently flooded
                   &   0.05_wp,  0.8_wp,  2.0_wp, 0.7_wp, 120.0_wp,  0.76_wp, 1._wp, & ! closed to open grassland regularly flooded
                   &   1.00_wp,  0.2_wp,  1.6_wp, 0.2_wp, 300.0_wp,  0.50_wp, 1._wp, & ! artificial surfaces
                   &   0.05_wp,  0.05_wp, 0.6_wp,0.05_wp, 300.0_wp,  0.82_wp, 1._wp, & ! bare areas
                   &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 150.0_wp,  -1.0_wp,-1._wp, & ! water bodies
                   &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! permanent snow and ice
                   &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp, 250.0_wp,  -1.0_wp,-1._wp  / ! undefined

! Yet another tuned version by Guenther Zaengl (adjusted to resistance-based bare soil evaporation scheme)
 DATA lu_gcv2009_v4 /  0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 190.0_wp,  0.72_wp, 1._wp, & ! irrigated croplands
                   &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 140.0_wp,  0.72_wp, 1._wp, & ! rainfed croplands
                   &   0.25_wp,  0.8_wp,  3.0_wp, 0.8_wp, 130.0_wp,  0.55_wp, 1._wp, & ! mosaic cropland (50-70%) - vegetation (20-50%)
                   &   0.07_wp,  0.9_wp,  3.5_wp, 1.0_wp, 120.0_wp,  0.72_wp, 1._wp, & ! mosaic vegetation (50-70%) - cropland (20-50%)
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp, 1._wp, & ! closed broadleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 1.0_wp, 300.0_wp,  0.31_wp, 1._wp, & ! closed broadleaved deciduous forest
                   &   0.50_wp,  0.8_wp,  4.0_wp, 1.5_wp, 225.0_wp,  0.31_wp, 1._wp, & ! open broadleaved deciduous forest
                   &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 300.0_wp,  0.27_wp, 1._wp, & ! closed needleleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 300.0_wp,  0.33_wp, 1._wp, & ! open needleleaved deciduous forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 270.0_wp,  0.29_wp, 1._wp, & ! mixed broadleaved and needleleaved forest
                   &   0.20_wp,  0.8_wp,  2.5_wp, 0.8_wp, 200.0_wp,  0.60_wp, 1._wp, & ! mosaic shrubland (50-70%) - grassland (20-50%)
                   &   0.20_wp,  0.8_wp,  2.5_wp, 0.6_wp, 180.0_wp,  0.65_wp, 1._wp, & ! mosaic grassland (50-70%) - shrubland (20-50%)
                   &   0.15_wp,  0.8_wp,  2.5_wp, 0.9_wp, 265.0_wp,  0.65_wp, 1._wp, & ! closed to open shrubland
                   &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp, 100.0_wp,  0.82_wp, 1._wp, & ! closed to open herbaceous vegetation
                   &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp, 140.0_wp,  0.76_wp, 1._wp, & ! sparse vegetation
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 190.0_wp,  0.30_wp, 1._wp, & ! closed to open forest regulary flooded
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 190.0_wp,  0.30_wp, 1._wp, & ! closed forest or shrubland permanently flooded
                   &   0.05_wp,  0.8_wp,  2.0_wp, 1.0_wp,  80.0_wp,  0.76_wp, 1._wp, & ! closed to open grassland regularly flooded
                   &   1.00_wp,  0.2_wp,  1.6_wp, 0.6_wp, 300.0_wp,  0.50_wp, 1._wp, & ! artificial surfaces
                   &   0.05_wp,  0.01_wp, 0.2_wp, 0.3_wp, 300.0_wp,  0.82_wp, 1._wp, & ! bare areas
                   &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 150.0_wp,  -1.0_wp,-1._wp, & ! water bodies
                   &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp, 1._wp, & ! permanent snow and ice
                   &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp, 250.0_wp,  -1.0_wp,-1._wp  / ! undefined


    !----------------------------------------------------------------------

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    IF ( itopo == 1 .AND. ( iforcing == iecham .OR. iforcing == ildf_echam ) ) THEN

      ! Read elevation of grid cells centers from grid file; this is then used to dynamically "grow" a topography for
      ! the hydrostatic model (in mo_ha_diag_util). This should be removed once the echam atmosphere is realistically
      ! initialized and uses a real topography.

      DO jg = 1,n_dom

        stream_id = openInputFile(p_patch(jg)%grid_filename, p_patch(jg), &
          &                       default_read_method)

        ! get land-sea-mask on cells, integer marks are:
        ! inner sea (-2), boundary sea (-1, cells and vertices), boundary (0, edges),
        ! boundary land (1, cells and vertices), inner land (2)
        CALL read_2D_int(stream_id, on_cells, 'cell_sea_land_mask', &
          &              ext_data(jg)%atm%lsm_ctr_c)

        ! get topography [m]
        ! - The hydrostatic AMIP setup grows the topography form zero to the elevation
        !   read from the grid file. Therefore the read in topography is stored in
        !   'elevation_c' and the actual 'topography_c' is computed later.
        ! - The non-hydrostatic AMIP setup starts directly from the topography read from
        !   the grid file. Hence the read in topography is stored in 'topography_c'.
        SELECT CASE (iequations)
        CASE (ihs_atm_temp,ihs_atm_theta) ! iequations
          ! Read topography
          CALL read_2D(stream_id, on_cells, 'cell_elevation', &
            &          ext_data(jg)%atm%elevation_c)
          ! Mask out ocean
          ext_data(jg)%atm%elevation_c(:,:) = MERGE(ext_data(jg)%atm%elevation_c(:,:), 0._wp, &
            &                                       ext_data(jg)%atm%lsm_ctr_c(:,:)  > 0     )
        CASE (inh_atmosphere) ! iequations
          ! Read topography
          CALL read_2D(stream_id, on_cells, 'cell_elevation', &
            &          ext_data(jg)%atm%topography_c)
          ! Mask out ocean
          ext_data(jg)%atm%topography_c(:,:) = MERGE(ext_data(jg)%atm%topography_c(:,:), 0._wp, &
            &                                        ext_data(jg)%atm%lsm_ctr_c(:,:)   > 0     )
        END SELECT ! iequations

        CALL closeFile(stream_id)

        ! LW surface emissivity
        !
        ! (eventually emis_rad should be read from file)
        !
        ext_data(jg)%atm%emis_rad(:,:)= zemiss_def


      END DO

    END IF

    ! Open/Read slm for HDmodel in configuration with jsbach, used in yac-coupler
    IF ( is_coupled_run() .AND. echam_phy_config%ljsbach) THEN

      DO jg = 1,n_dom

        stream_id = openInputFile('hd_mask.nc', p_patch(jg), default_read_method)
     
        ! get land-sea-mask on cells, integer marks are:
        ! inner sea (-2), boundary sea (-1, cells and vertices), boundary (0, edges),
        ! boundary land (1, cells and vertices), inner land (2)
        CALL read_2D_int(stream_id, on_cells, 'cell_sea_land_mask', &
          &              ext_data(jg)%atm%lsm_hd_c)

        CALL closeFile(stream_id)

      END DO

    END IF

    !------------------------------------------------!
    ! Read data from ExtPar file                     !
    !------------------------------------------------!

    IF (itopo == 1 .AND. iforcing == inwp) THEN
      DO jg = 1,n_dom

        ! Preset parameter fields with the correct table values
        ilu = 0
        IF (i_lctype(jg) == GLC2000) THEN
          ext_data(jg)%atm%i_lc_snow_ice = 21
          ext_data(jg)%atm%i_lc_water    = 20
          ext_data(jg)%atm%i_lc_urban    = 22
          ext_data(jg)%atm%i_lc_shrub_eg = 11
          ext_data(jg)%atm%i_lc_shrub    = 12
          ext_data(jg)%atm%i_lc_grass    = 13
          ext_data(jg)%atm%i_lc_bare_soil= 19
          ext_data(jg)%atm%i_lc_sparse   = 14
          DO i = 1, num_lcc*n_param_lcc, n_param_lcc
            ilu=ilu+1
            ext_data(jg)%atm%z0_lcc(ilu)          = lu_glc2000(i  )  ! Land-cover related roughness length
            ext_data(jg)%atm%plcovmax_lcc(ilu)    = lu_glc2000(i+1)  ! Maximum plant cover fraction for each land-cover class
            ext_data(jg)%atm%laimax_lcc(ilu)      = lu_glc2000(i+2)  ! Maximum leaf area index for each land-cover class
            ext_data(jg)%atm%rootdmax_lcc(ilu)    = lu_glc2000(i+3)  ! Maximum root depth for each land-cover class
            ext_data(jg)%atm%stomresmin_lcc(ilu)  = lu_glc2000(i+4)  ! Minimum stomata resistance for each land-cover class
            ext_data(jg)%atm%snowalb_lcc(ilu)     = lu_glc2000(i+5)  ! Albedo in case of snow cover for each land-cover class
            ext_data(jg)%atm%snowtile_lcc(ilu)    = &
              &          MERGE(.TRUE.,.FALSE.,lu_glc2000(i+6)>0._wp) ! Existence of snow tiles for land-cover class
          ENDDO
        ELSE IF (i_lctype(jg) == GLOBCOVER2009) THEN
          SELECT CASE (itype_lndtbl)
          CASE (1)
            lu_gcv => lu_gcv2009
          CASE (2)
            lu_gcv => lu_gcv2009_v2
          CASE (3)
            lu_gcv => lu_gcv2009_v3
          CASE (4)
            lu_gcv => lu_gcv2009_v4
          END SELECT

          ext_data(jg)%atm%i_lc_snow_ice = 22
          ext_data(jg)%atm%i_lc_water    = 21
          ext_data(jg)%atm%i_lc_urban    = 19
          ext_data(jg)%atm%i_lc_shrub_eg = 12
          ext_data(jg)%atm%i_lc_shrub    = 13
          ext_data(jg)%atm%i_lc_grass    = 14
          ext_data(jg)%atm%i_lc_bare_soil= 20
          ext_data(jg)%atm%i_lc_sparse   = 15
          DO i = 1, num_lcc*n_param_lcc, n_param_lcc
            ilu=ilu+1
            ext_data(jg)%atm%z0_lcc(ilu)          = lu_gcv(i  )  ! Land-cover related roughness length
            ext_data(jg)%atm%plcovmax_lcc(ilu)    = lu_gcv(i+1)  ! Maximum plant cover fraction for each land-cover class
            ext_data(jg)%atm%laimax_lcc(ilu)      = lu_gcv(i+2)  ! Maximum leaf area index for each land-cover class
            ext_data(jg)%atm%rootdmax_lcc(ilu)    = lu_gcv(i+3)  ! Maximum root depth for each land-cover class
            ext_data(jg)%atm%stomresmin_lcc(ilu)  = lu_gcv(i+4)  ! Minimum stomata resistance for each land-cover class
            ext_data(jg)%atm%snowalb_lcc(ilu)     = lu_gcv(i+5)  ! Albedo in case of snow cover for each land-cover class
            ext_data(jg)%atm%snowtile_lcc(ilu)    = &
              &          MERGE(.TRUE.,.FALSE.,lu_gcv(i+6)>0._wp) ! Existence of snow tiles for land-cover class
          ENDDO
        ENDIF

        ! Derived parameter: minimum allowed land-cover related roughness length in the
        ! presence of low ndvi and/or snow cover
        DO ilu = 1, num_lcc
          IF (ilu == ext_data(jg)%atm%i_lc_urban .OR. ilu == ext_data(jg)%atm%i_lc_water) THEN
            ext_data(jg)%atm%z0_lcc_min(ilu) = ext_data(jg)%atm%z0_lcc(ilu) ! no reduction in urban regions and over water
          ELSE IF (ext_data(jg)%atm%z0_lcc(ilu) >= 0.1) THEN
            ext_data(jg)%atm%z0_lcc_min(ilu) = 0.3_wp*ext_data(jg)%atm%z0_lcc(ilu) ! 30% for nominal roughness lengths > 10 cm
          ELSE
            ext_data(jg)%atm%z0_lcc_min(ilu) = 0.1_wp*ext_data(jg)%atm%z0_lcc(ilu) ! 10% otherwise
          ENDIF
        ENDDO

        ! Start reading external parameter data
        ! The cdi-based read routines are used for GRIB2 input data only due to performance problems
        IF (read_netcdf_data) THEN
          extpar_file = generate_filename(extpar_filename, getModelBaseDir(), &
            &                             TRIM(p_patch(jg)%grid_filename),    &
            &                              nroot,                             &
            &                             p_patch(jg)%level, p_patch(jg)%id)
          stream_id   = openInputFile(extpar_file, p_patch(jg), default_read_method)
        ELSE
          parameters = makeInputParameters(cdi_extpar_id(jg), p_patch(jg)%n_patch_cells_g, p_patch(jg)%comm_pat_scatter_c, &
          &                                opt_dict=extpar_varnames_dict)
        ENDIF

        !--------------------------------------------------------------------
        !
        ! Read topography for triangle centers (triangular grid)
        !
        !--------------------------------------------------------------------
        CALL read_extdata('topography_c', ext_data(jg)%atm%topography_c)

        !
        ! other external parameters on triangular grid
        !
        CALL read_extdata('FR_LAND', ext_data(jg)%atm%fr_land)

        SELECT CASE ( iforcing )
        CASE ( inwp )
        CALL read_extdata('NDVI_MAX',  ext_data(jg)%atm%ndvi_max)
        CALL read_extdata('SOILTYP',   arr2di=ext_data(jg)%atm%soiltyp)
        CALL read_extdata('T_CL',      ext_data(jg)%atm%t_cl)
        CALL read_extdata('SSO_STDH',  ext_data(jg)%atm%sso_stdh)
        CALL read_extdata('SSO_THETA', ext_data(jg)%atm%sso_theta)
        CALL read_extdata('SSO_GAMMA', ext_data(jg)%atm%sso_gamma)
        CALL read_extdata('SSO_SIGMA', ext_data(jg)%atm%sso_sigma)
        CALL read_extdata('FR_LAKE',   ext_data(jg)%atm%fr_lake)
        CALL read_extdata('DEPTH_LK',  ext_data(jg)%atm%depth_lk)

        CALL read_extdata('LU_CLASS_FRACTION', arr3d=ext_data(jg)%atm%lu_class_fraction,ltime=.FALSE.) 

          ! The following fields are only required without surface tiles
          IF (ntiles_lnd == 1) THEN
            CALL read_extdata('PLCOV_MX', ext_data(jg)%atm%plcov_mx)
            CALL read_extdata('LAI_MX',   ext_data(jg)%atm%lai_mx)
            CALL read_extdata('ROOTDP',   ext_data(jg)%atm%rootdp)
            CALL read_extdata('RSMIN',    ext_data(jg)%atm%rsmin)
            CALL read_extdata('FOR_D',    ext_data(jg)%atm%for_d)
            CALL read_extdata('FOR_E',    ext_data(jg)%atm%for_e)
          ENDIF

          IF (atm_phy_nwp_config(jg)%itype_z0 == 1) THEN
            ! only read, if contribution from sub-scale orography should be included in z0
            CALL read_extdata('Z0', ext_data(jg)%atm%z0)
          ENDIF

          IF (is_frglac_in(jg)) THEN
            ! for backward compatibility with extpar files generated prior to 2014-01-31
             CALL read_extdata('ICE', ext_data(jg)%atm%fr_glac)
          ELSE
            ! for new extpar files (generated after 2014-01-31)
            ! take it from lu_class_fraction
            ext_data(jg)%atm%fr_glac(:,:) = ext_data(jg)%atm%lu_class_fraction(:,:,ext_data(jg)%atm%i_lc_snow_ice)
          ENDIF

          IF ( l_emiss ) THEN
            CALL read_extdata('EMIS_RAD', ext_data(jg)%atm%emis_rad)
          ELSE
            ext_data(jg)%atm%emis_rad(:,:)= zemiss_def
          ENDIF

          ! Copy sso_stdh to sso_stdh_raw before applying correction for orography filtering
          ext_data(jg)%atm%sso_stdh_raw(:,:) = ext_data(jg)%atm%sso_stdh(:,:)


          ! Read time dependent data
          IF ( irad_aero == 6 .OR. irad_aero == 9) THEN
            CALL read_extdata('AER_SS',   arr3d=ext_data(jg)%atm_td%aer_ss)
            CALL read_extdata('AER_DUST', arr3d=ext_data(jg)%atm_td%aer_dust)
            CALL read_extdata('AER_ORG',  arr3d=ext_data(jg)%atm_td%aer_org)
            CALL read_extdata('AER_SO4',  arr3d=ext_data(jg)%atm_td%aer_so4)
            CALL read_extdata('AER_BC',   arr3d=ext_data(jg)%atm_td%aer_bc)
          ENDIF  ! irad_aero
          CALL read_extdata('NDVI_MRAT', arr3d=ext_data(jg)%atm_td%ndvi_mrat)

          IF (sstice_mode == SSTICE_ANA_CLINC) THEN
            CALL read_extdata('T_SEA', arr3d=ext_data(jg)%atm_td%sst_m)
          ENDIF

          !--------------------------------
          ! If MODIS albedo is used
          !--------------------------------
          IF ( albedo_type == MODIS) THEN
            CALL read_extdata('ALB',   arr3d=ext_data(jg)%atm_td%alb_dif)
            CALL read_extdata('ALUVD', arr3d=ext_data(jg)%atm_td%albuv_dif)
            CALL read_extdata('ALNID', arr3d=ext_data(jg)%atm_td%albni_dif)

            rl_start = 1
            rl_end   = min_rlcell

            i_startblk = p_patch(jg)%cells%start_block(rl_start)
            i_endblk   = p_patch(jg)%cells%end_block(rl_end)

            albthresh = 0.3_wp ! threshold value for albedo modification

!$OMP PARALLEL
            ! Scale from [%] to [1]
            CALL var_scale(ext_data(jg)%atm_td%alb_dif(:,:,:), 1._wp/100._wp)
            CALL var_scale(ext_data(jg)%atm_td%albuv_dif(:,:,:), 1._wp/100._wp)
            CALL var_scale(ext_data(jg)%atm_td%albni_dif(:,:,:), 1._wp/100._wp)
!$OMP BARRIER


            IF (itune_albedo >= 1) THEN
              ! Test: reduce albedo over land where modis albedo is higher than 0.3 (variable albthresh)
!$OMP DO PRIVATE(jb,jc,im,i_startidx,i_endidx,albfac)
              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

                DO im = 1, 12
                  DO jc = i_startidx,i_endidx
                    IF (ext_data(jg)%atm%soiltyp(jc,jb) >= 2 .AND. ext_data(jg)%atm%soiltyp(jc,jb) <= 8) THEN
                      IF (ext_data(jg)%atm_td%alb_dif(jc,jb,im) > albthresh) THEN
                        albfac = (albthresh+2._wp*ext_data(jg)%atm_td%alb_dif(jc,jb,im))/ &
                          (3._wp*ext_data(jg)%atm_td%alb_dif(jc,jb,im))
                        ext_data(jg)%atm_td%alb_dif(jc,jb,im)   = albfac*ext_data(jg)%atm_td%alb_dif(jc,jb,im)
                        ext_data(jg)%atm_td%albuv_dif(jc,jb,im) = albfac*ext_data(jg)%atm_td%albuv_dif(jc,jb,im)
                        ext_data(jg)%atm_td%albni_dif(jc,jb,im) = albfac*ext_data(jg)%atm_td%albni_dif(jc,jb,im)
                      ENDIF
                    ENDIF
                  ENDDO
               ENDDO
              ENDDO
!$OMP END DO
            ENDIF  ! Sahara albedo tuning

            IF (itune_albedo >= 2) THEN
              ! Increase albedo over the Antarctic plateau by 5% (from 70% to 75%) in order to get rid of summertime warm bias
!$OMP DO PRIVATE(jb,jc,im,i_startidx,i_endidx,albfac)
              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

                DO im = 1, 12
                  DO jc = i_startidx,i_endidx
                    IF (ext_data(jg)%atm%soiltyp(jc,jb) == 1 .AND. p_patch(jg)%cells%center(jc,jb)%lat*rad2deg < -65._wp ) THEN
                      IF (ext_data(jg)%atm%topography_c(jc,jb) > 1000._wp) THEN
                        albfac = MIN(1._wp,1.e-3_wp*(ext_data(jg)%atm%topography_c(jc,jb)-1000._wp))
                        ext_data(jg)%atm_td%alb_dif(jc,jb,im)   = 0.05_wp*albfac + ext_data(jg)%atm_td%alb_dif(jc,jb,im)
                        ext_data(jg)%atm_td%albuv_dif(jc,jb,im) = 0.05_wp*albfac + ext_data(jg)%atm_td%albuv_dif(jc,jb,im)
                        ext_data(jg)%atm_td%albni_dif(jc,jb,im) = 0.05_wp*albfac + ext_data(jg)%atm_td%albni_dif(jc,jb,im)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
!$OMP END DO
            ENDIF  ! Antarctic albedo tuning
!$OMP END PARALLEL

          END IF  !  albedo_type

        END SELECT ! iforcing

        IF (read_netcdf_data) THEN
          CALL closeFile(stream_id)
        ELSE
          CALL deleteInputParameters(parameters)
        ENDIF

        !
        ! derived external parameter fields
        !

        ! land sea mask at cell centers (LOGICAL)
        !
        i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

        rl_start = 1
        rl_end   = min_rlcell

        i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
        i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

        DO jb = i_startblk, i_endblk
          CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
            &                i_startidx, i_endidx, rl_start, rl_end)

          ! Loop starts with 1 instead of i_startidx
          ! because the start index is missing in RRTM
          DO jc = 1,i_endidx
            IF (ext_data(jg)%atm%fr_land(jc,jb) > 0.5_wp) THEN
              ext_data(jg)%atm%llsm_atm_c(jc,jb) = .TRUE.  ! land point
            ELSE
              ext_data(jg)%atm%llsm_atm_c(jc,jb) = .FALSE.  ! water point
            ENDIF
            IF (ext_data(jg)%atm%fr_lake(jc,jb) >= 0.5_wp) THEN
              ext_data(jg)%atm%llake_c(jc,jb) = .TRUE.   ! lake point
            ELSE
              ext_data(jg)%atm%llake_c(jc,jb) = .FALSE.  ! no lake point
            ENDIF
          ENDDO
        ENDDO

        ! As long as a routine for computing smoothed external
        ! parameter fields is missing. Just copy.
        !
        DO jb = i_startblk, i_endblk
          CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
            &                i_startidx, i_endidx, rl_start, rl_end)
          ! Loop starts with 1 instead of i_startidx
          ! because the start index is missing in RRTM
          DO jc = 1,i_endidx
            ext_data(jg)%atm%fr_land_smt(jc,jb) = ext_data(jg)%atm%fr_land(jc,jb)
            ext_data(jg)%atm%fr_glac_smt(jc,jb) = ext_data(jg)%atm%fr_glac(jc,jb)
          ENDDO

        ENDDO

      ENDDO  ! jg

    ENDIF ! (itopo == 1)

    !-------------------------------------------------------
    ! Read ozone
    !-------------------------------------------------------

    IF((irad_o3 == io3_clim) .OR. (irad_o3 == io3_ape)) THEN

      DO jg = 1,n_dom

        WRITE(ozone_file,'(a,I2.2,a)') 'o3_icon_DOM',jg,'.nc'

        IF(my_process_is_stdio()) THEN
          ! open file
          !
          CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid), routine)
          WRITE(0,*)'read ozone levels'
          CALL nf(nf_inq_varid(ncid, TRIM(levelname), varid), routine)
          CALL nf(nf_get_var_double(ncid, varid, zdummy_o3lev(:)), routine)
          CALL nf(nf_close(ncid), routine)
          !
        ENDIF ! pe

        CALL p_bcast(zdummy_o3lev(:), p_io, mpi_comm)

        !         SELECT CASE (iequations)
        !         CASE(ihs_atm_temp,ihs_atm_theta)

        DO jk=1,nlev_o3
          ext_data(jg)%atm_td%pfoz(jk)=zdummy_o3lev(jk)
        ENDDO

        ! define half levels of ozone pressure grid
        ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
        ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
        ext_data(jg)%atm_td%phoz(1)           = 0._wp
        ext_data(jg)%atm_td%phoz(2:nlev_o3) = (ext_data(jg)%atm_td%pfoz(1:nlev_o3-1) &
          &                                   +  ext_data(jg)%atm_td%pfoz(2:nlev_o3))*.5_wp
        ext_data(jg)%atm_td%phoz(nlev_o3+1) = 125000._wp

        DO i=1,nlev_o3

          WRITE(message_text,'(a,i4,f12.4,f12.4)')'full/half level press ozone ', &
                              i, ext_data(jg)%atm_td%pfoz(i), ext_data(jg)%atm_td%phoz(i+1)
          CALL message(routine, TRIM(message_text))
        ENDDO

        stream_id = openInputFile(ozone_file, p_patch(jg), default_read_method)

        CALL read_3D_extdim(stream_id, on_cells, TRIM(o3name), &
          &                 ext_data(jg)%atm_td%O3)

        WRITE(message_text,'(a,f12.4,f12.4)')'MAX/MIN o3 ppmv', &
           MAXVAL(ext_data(jg)%atm_td%O3(:,:,:,:)), MINVAL(ext_data(jg)%atm_td%O3(:,:,:,:))
        CALL message(routine, TRIM(message_text))

        ! convert from ppmv to g/g only in case of APE ozone
        ! whether o3mr2gg or ppmv2gg is used to convert O3 to gg depends on the units of
        ! the incoming ozone file.  Often, the incoming units are not ppmv.
        IF(irad_o3 == io3_ape) &
         ! &         ext_data(jg)%atm_td%O3(:,:,:,:)= ext_data(jg)%atm_td%O3(:,:,:,:)*o3mr2gg
          &         ext_data(jg)%atm_td%O3(:,:,:,:)= ext_data(jg)%atm_td%O3(:,:,:,:)*ppmv2gg


        WRITE(message_text,'(a,f12.4,f12.4)')'MAX/MIN o3 g/g', &
           MAXVAL(ext_data(jg)%atm_td%O3(:,:,:,:)), MINVAL(ext_data(jg)%atm_td%O3(:,:,:,:))
        CALL message(routine, TRIM(message_text))

        ! close file
        CALL closeFile(stream_id)

      ENDDO ! ndom
    END IF ! irad_o3


    !------------------------------------------
    ! Read time dependent SST and ICE Fraction
    !------------------------------------------
    IF (sstice_mode == SSTICE_CLIM .AND. iforcing == inwp) THEN

      DO jg = 1,n_dom
       !Read the climatological values for SST and ice cover

        DO im=1,12

         sst_td_file= generate_td_filename(sst_td_filename,                &
           &                             getModelBaseDir(),                &
           &                             TRIM(p_patch(jg)%grid_filename),  &
           &                             im,clim=.TRUE.                   )

         IF(my_process_is_mpi_workroot()) THEN

          CALL message  (routine, TRIM(sst_td_file))

          INQUIRE (FILE=sst_td_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            CALL finish(routine,'td sst external data file is not found.')
          ENDIF

         ENDIF

         stream_id = openInputFile(sst_td_file, p_patch(jg), &
          &                        default_read_method)
         CALL read_2D (stream_id, on_cells, 'SST', &
          &            ext_data(jg)%atm_td%sst_m(:,:,im))
         CALL closeFile(stream_id)

         ci_td_file= generate_td_filename(ci_td_filename,                  &
           &                             getModelBaseDir(),                &
           &                             TRIM(p_patch(jg)%grid_filename),  &
           &                             im,clim=.TRUE.                   )

         IF(my_process_is_stdio()) THEN

          CALL message  (routine, TRIM(ci_td_file))

          INQUIRE (FILE=ci_td_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            CALL finish(routine,'td ci external data file is not found.')
          ENDIF

         ENDIF

         stream_id = openInputFile(ci_td_file, p_patch(jg), default_read_method)
         CALL read_2D(stream_id, on_cells, 'CI', &
          &           ext_data(jg)%atm_td%fr_ice_m(:,:,im))
         CALL closeFile(stream_id)

        END DO

      END DO ! ndom

   END IF ! sstice_mode

   CONTAINS

     ! Wrapper routine for reading input data via cdilib for GRIB2 or via optimized netcdf routines
     !
     SUBROUTINE read_extdata(varname,arr2d,arr2di,arr3d,ltime)
       CHARACTER(LEN=*), INTENT(IN) :: varname  ! name of input variable
       REAL(wp), OPTIONAL, INTENT(INOUT) :: arr2d(:,:), arr3d(:,:,:) ! alternative I/O arrays
       INTEGER, OPTIONAL, INTENT(INOUT) :: arr2di(:,:)
       LOGICAL, OPTIONAL, INTENT(IN) :: ltime ! .true. if third dimension is time

       LOGICAL :: dim3_is_time

       IF (PRESENT(arr3d) .AND. PRESENT(ltime)) THEN
         dim3_is_time = ltime
       ELSE
         dim3_is_time = .TRUE.
       ENDIF

       IF (PRESENT(arr2d)) THEN
         IF (read_netcdf_data) THEN
           CALL read_2D(stream_id, on_cells, TRIM(varname), arr2d)
         ELSE
           CALL read_cdi_2d(parameters, TRIM(varname), arr2d)
         ENDIF
       ELSE IF (PRESENT(arr2di)) THEN
         IF (read_netcdf_data) THEN
           CALL read_2D_int(stream_id, on_cells, TRIM(varname), arr2di)
         ELSE
           CALL read_cdi_2d(parameters, TRIM(varname), arr2di)
         ENDIF
       ELSE IF (PRESENT(arr3d)) THEN
         IF (read_netcdf_data) THEN
           CALL read_2D_extdim(stream_id, on_cells, TRIM(varname), arr3d)
         ELSE IF (dim3_is_time) THEN
           CALL read_cdi_2d(parameters, SIZE(arr3d,3), TRIM(varname), arr3d)
         ELSE
           CALL read_cdi_3d(parameters, TRIM(varname), SIZE(arr3d,3), arr3d, opt_lev_dim=3 )
         ENDIF
       ENDIF

     END SUBROUTINE read_extdata

  END SUBROUTINE read_ext_data_atm
  !-------------------------------------------------------------------------



  SUBROUTINE init_index_lists (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    INTEGER :: i_lu, jb,jc, jg, i_count, i_count_flk, ic, jt, jt_in
    INTEGER :: i_count_sea             ! number of sea points
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    LOGICAL  :: tile_mask(num_lcc)
    REAL(wp) :: tile_frac(num_lcc), sum_frac
    INTEGER  :: lu_subs, it_count(ntiles_total)
    INTEGER  :: npoints, npoints_sea, npoints_lake
    INTEGER  :: i_lc_water

    REAL(wp), POINTER  ::  &  !< pointer to proportion of actual value/maximum
      &  ptr_ndviratio(:,:)   !< NDVI (for starting time of model integration)

    REAL(wp) :: scalfac       ! scaling factor for inflating dominating land fractions
                              ! to fr_land (or fr_land + fr_lake)
    REAL(wp) :: zfr_land      ! fr_land derived from land tile fractions

!!$    CHARACTER(len=max_char_length), PARAMETER :: &
!!$      routine = modname//':init_index_lists'

    !-------------------------------------------------------------------------

    WRITE(message_text,'(a,i4)') &
      &  'Index list generation - number of land tiles: ', ntiles_lnd
    CALL message('', TRIM(message_text))
    WRITE(message_text,'(a,i4)')  'Total number of tiles: ', ntiles_total
    CALL message('', TRIM(message_text))


    DO jg = 1, n_dom

       ptr_ndviratio => ext_data(jg)%atm%ndviratio(:,:)

       i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

       i_lc_water = ext_data(jg)%atm%i_lc_water

       ! Initialization of index list counts - moved here in order to avoid uninitialized elements
       ! along nest boundaries
       ext_data(jg)%atm%lp_count(:) = 0
       ext_data(jg)%atm%sp_count(:) = 0
       ext_data(jg)%atm%fp_count(:) = 0

       ext_data(jg)%atm%spi_count(:) = 0
       ext_data(jg)%atm%spw_count(:) = 0

       ext_data(jg)%atm%gp_count_t(:,:) = 0
       ext_data(jg)%atm%lp_count_t(:,:) = 0


!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
       !
       ! exclude the boundary interpolation zone of nested domains
       rl_start = grf_bdywidth_c+1
       rl_end   = min_rlcell_int

       i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
       i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_lu,i_startidx,i_endidx,i_count,i_count_sea,i_count_flk,tile_frac,&
!$OMP            tile_mask,lu_subs,sum_frac,scalfac,zfr_land,it_count,ic,jt,jt_in ) ICON_OMP_DEFAULT_SCHEDULE
       DO jb=i_startblk, i_endblk

         CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

         i_count                       = 0   ! counter for land points
         i_count_sea                   = 0   ! counter for sea points
         i_count_flk                   = 0   ! counter for lake points

         it_count(:)                   = 0 ! counter for tiles

         DO jc = i_startidx, i_endidx
           ext_data(jg)%atm%lc_class_t(jc,jb,:) = -1    ! dummy value for undefined points

           IF (ext_data(jg)%atm%fr_land(jc,jb)> frlnd_thrhld) THEN ! searching for land-points
             i_count=i_count+1
             ext_data(jg)%atm%idx_lst_lp(i_count,jb) = jc  ! write index of land-points

             ! SBr, CHa: set tile_frac to 1
             !tile_frac(:)= ext_data(jg)%atm%lu_class_fraction(jc,jb,:)
             tile_frac(:)= 1._wp
             tile_mask(:)=.true.
             tile_mask(i_lc_water)=.false. ! exclude water points

             ext_data(jg)%atm%lp_count(jb) = i_count

             IF (ntiles_lnd == 1) THEN

               ! i_lu=1 contains grid-box mean values from EXTPAR!
               !
               ext_data(jg)%atm%lc_frac_t(jc,jb,1)  = 1._wp
               ext_data(jg)%atm%lc_class_t(jc,jb,1) = MAXLOC(tile_frac,1,tile_mask)
               !
               ! root depth
               ext_data(jg)%atm%rootdp_t (jc,jb,1)  = ext_data(jg)%atm%rootdp(jc,jb)

               ! plant cover
               ext_data(jg)%atm%plcov_t  (jc,jb,1)  = ptr_ndviratio(jc,jb)  &
                 &     * MIN(ext_data(jg)%atm%ndvi_max(jc,jb),ext_data(jg)%atm%plcov_mx(jc,jb))
               ! transpiration area index
               ext_data(jg)%atm%tai_t    (jc,jb,1)  = ext_data(jg)%atm%plcov_t  (jc,jb,1)  &
                 &                                  * ext_data(jg)%atm%lai_mx(jc,jb)
               ! surface area index
               ext_data(jg)%atm%sai_t    (jc,jb,1)  = c_lnd+ext_data(jg)%atm%tai_t(jc,jb,1)
               ! evaporative soil area index
               ext_data(jg)%atm%eai_t(jc,jb,1) = &
                 MERGE(c_soil_urb,c_soil,ext_data(jg)%atm%lc_class_t(jc,jb,1) == ext_data(jg)%atm%i_lc_urban)
               ! minimal stomata resistance
               ext_data(jg)%atm%rsmin2d_t(jc,jb,1)  = ext_data(jg)%atm%rsmin(jc,jb)
               ! soil type
               ext_data(jg)%atm%soiltyp_t(jc,jb,1)  = ext_data(jg)%atm%soiltyp(jc,jb)

               ! Workaround for GLC2000 hole below 60 deg S
               ! (only necesary for old extpar files generated prior to 2014-01-31)
               IF (is_frglac_in(jg)) THEN
                 IF (tile_frac(ext_data(jg)%atm%lc_class_t(jc,jb,1))<=0._wp) &
                   ext_data(jg)%atm%lc_class_t(jc,jb,1) = ext_data(jg)%atm%i_lc_snow_ice
               ENDIF

               ! static index list and corresponding counter
               ext_data(jg)%atm%idx_lst_lp_t(i_count,jb,1)  = jc
               ext_data(jg)%atm%lp_count_t(jb,1)            = i_count

               ! initialize dynamic index list (in case of lsnowtile=true) with the same values
               ext_data(jg)%atm%idx_lst_t(i_count,jb,1) = jc
               ext_data(jg)%atm%gp_count_t(jb,1)        = i_count

               ! initialize snowtile flag with 1 if the tile is eligible for separate treatment of
               ! a snow-covered and a snow-free part, otherwise with -1
               IF (ext_data(jg)%atm%snowtile_lcc(ext_data(jg)%atm%lc_class_t(jc,jb,1))) THEN
                 ext_data(jg)%atm%snowtile_flag_t(jc,jb,1)  = 1
               ELSE
                 ext_data(jg)%atm%snowtile_flag_t(jc,jb,1)  = -1
               ENDIF

             ELSE
               ext_data(jg)%atm%lc_frac_t(jc,jb,:)  = 0._wp ! to be really safe

               DO i_lu = 1, ntiles_lnd
                 lu_subs = MAXLOC(tile_frac,1,tile_mask)
                 ! Note that we have to take into account points with fr_land > frlnd_thrhld but
                 ! maximum_tile_frac <frlndtile_thrhld (for tile=1).
                 ! This e.g. may be the case at non-dominant land points with very inhomogeneous land class coverage.
                 ! That's why checking for (tile_frac(lu_subs) >= frlndtile_thrhld) in the next If-statement is not enough.
                 ! We accept class fractions for tile=1 even if tile_frac << frlndtile_thrhld.
                 !
                 ! The additional check tile_frac(lu_subs) >= 1.e-03_wp is only added for backward compatibility and is
                 ! required by all extpar-files generated prior to 2014-01-30. In these files it is not assured that
                 ! SUM(tile_frac(:))=1. I.e. glacier below 60degS are missing, such that SUM(tile_frac(:))=0 in these cases.
                 !
                 IF ( (i_lu==1 .AND. tile_frac(lu_subs) >= 1.e-03_wp) .OR. (tile_frac(lu_subs) >= frlndtile_thrhld) ) THEN
                   it_count(i_lu)    = it_count(i_lu) + 1
                   tile_mask(lu_subs)= .FALSE.

                   ! static index list and corresponding counter
                   ext_data(jg)%atm%idx_lst_lp_t(it_count(i_lu),jb,i_lu) = jc
                   ext_data(jg)%atm%lp_count_t(jb,i_lu)                  = it_count(i_lu)

                   ! initialize dynamic index list (in case of lsnowtile=true) with the same values
                   ext_data(jg)%atm%idx_lst_t(it_count(i_lu),jb,i_lu) = jc
                   ext_data(jg)%atm%gp_count_t(jb,i_lu)               = it_count(i_lu)

                   ! initialize snowtile flag with 1 if the tile is eligible for separate treatment of
                   ! a snow-covered and a snow-free part, otherwise with -1
                   IF (ext_data(jg)%atm%snowtile_lcc(lu_subs)) THEN
                     ext_data(jg)%atm%snowtile_flag_t(jc,jb,i_lu)  = 1
                   ELSE
                     ext_data(jg)%atm%snowtile_flag_t(jc,jb,i_lu)  = -1
                   ENDIF

                   ext_data(jg)%atm%lc_frac_t(jc,jb,i_lu)  = tile_frac(lu_subs)
                   ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = lu_subs
                 ELSE
                   EXIT ! no more land cover classes exceeding the threshold
                 ENDIF

               END DO

               ! fix for non-dominant land points
               !!! only active for 'old' extpar datasets (20131009 and earlier) !!!
               IF (ext_data(jg)%atm%fr_land(jc,jb) < 0.5_wp) THEN
                 IF (ext_data(jg)%atm%soiltyp(jc,jb) == 9) THEN  ! sea water
                   ! reset soil type to sandy loam ...
                   ext_data(jg)%atm%soiltyp(jc,jb) = 4
                 ENDIF
                 IF (ptr_ndviratio(jc,jb) <= 0.0_wp) THEN  ! here: 0=extpar_missval
                   ! ... and reset ndviratio
                   ptr_ndviratio(jc,jb) = 0.5_wp
                 ENDIF
                 IF (ext_data(jg)%atm%ndvi_max(jc,jb) <= 0.0_wp ) THEN  ! here: 0=extpar_missval
                   ! ... and reset ndvi_max to meaningful value (needed for plant cover)
                   ext_data(jg)%atm%ndvi_max(jc,jb) = 0.8_wp
                 ENDIF
               ENDIF
!!$               IF (ext_data(jg)%atm%fr_land(jc,jb) < 0.5_wp) THEN
!!$                 ! fix for non-dominant land points: reset soil type to sandy loam ...
!!$                 ext_data(jg)%atm%soiltyp(jc,jb) = 4
!!$                 ! ... and reset ndviratio to 0.5
!!$                 ptr_ndviratio(jc,jb) = 0.5_wp
!!$               ENDIF

               sum_frac = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd))

               DO i_lu = 1, ntiles_lnd

                 !  Workaround for GLC2000 hole below 60 deg S
                 ! (only necesary for old extpar files generated prior to 2014-01-31)
                 IF (is_frglac_in(jg)) THEN
                   IF ( sum_frac < 1.e-10_wp) THEN
                     IF (i_lu == 1) THEN
                       it_count(i_lu)    = it_count(i_lu) + 1
                       ! static index list and corresponding counter
                       ext_data(jg)%atm%idx_lst_lp_t(it_count(i_lu),jb,i_lu) = jc
                       ext_data(jg)%atm%lp_count_t(jb,i_lu)               = it_count(i_lu)

                       ! initialize dynamic index list (in case of lsnowtile=true) with the same values
                       ext_data(jg)%atm%idx_lst_t(it_count(i_lu),jb,i_lu) = jc
                       ext_data(jg)%atm%gp_count_t(jb,i_lu)               = it_count(i_lu)

                       ! the snowtile flag is initialized with -1 here because the snow/ice class is not
                       ! eligible for separate consideration of a snow-free and a snow-covered part
                       ext_data(jg)%atm%snowtile_flag_t(jc,jb,i_lu)         = -1

                       ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = ext_data(jg)%atm%i_lc_snow_ice
                       ext_data(jg)%atm%lc_frac_t(jc,jb,i_lu)  = ext_data(jg)%atm%fr_land(jc,jb)
                     ELSE
                       ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = -1
                       ext_data(jg)%atm%lc_frac_t(jc,jb,i_lu)  = 0._wp
                     ENDIF
                   END IF  ! sum_frac < 1.e-10_wp
                 ENDIF  ! is_frglac_in(jg)

                 lu_subs = ext_data(jg)%atm%lc_class_t(jc,jb,i_lu)
                 IF (lu_subs < 0) CYCLE

                 ! root depth
                 ext_data(jg)%atm%rootdp_t (jc,jb,i_lu)  = ext_data(jg)%atm%rootdmax_lcc(lu_subs)
                 ! plant cover
                 ext_data(jg)%atm%plcov_t  (jc,jb,i_lu)  = ptr_ndviratio(jc,jb) &
                   & * MIN(ext_data(jg)%atm%ndvi_max(jc,jb),ext_data(jg)%atm%plcovmax_lcc(lu_subs))
                 ! transpiration area index
                 ext_data(jg)%atm%tai_t    (jc,jb,i_lu)  = ext_data(jg)%atm%plcov_t(jc,jb,i_lu) &
                   & * ext_data(jg)%atm%laimax_lcc(lu_subs)

                 ! surface area index
                 ext_data(jg)%atm%sai_t    (jc,jb,i_lu)  = c_lnd+ ext_data(jg)%atm%tai_t (jc,jb,i_lu)

                 ! evaporative soil area index
                 ext_data(jg)%atm%eai_t (jc,jb,i_lu)  = MERGE(c_soil_urb,c_soil,lu_subs == ext_data(jg)%atm%i_lc_urban)

                 ! minimal stomata resistance
                 ext_data(jg)%atm%rsmin2d_t(jc,jb,i_lu)  = ext_data(jg)%atm%stomresmin_lcc(lu_subs)

                 ! soil type
                 ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu)  = ext_data(jg)%atm%soiltyp(jc,jb)

                 ! consistency corrections for partly glaciered points
                 ! a) set soiltype to ice if landuse = ice (already done in extpar for dominant glacier points)
                 IF (ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) == ext_data(jg)%atm%i_lc_snow_ice) &
                   & ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu) = 1
                 ! b) set soiltype to rock or sandy loam if landuse /= ice and soiltype = ice
                 IF (ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu) == 1 .AND. &
                   & ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) /= ext_data(jg)%atm%i_lc_snow_ice) THEN
                   IF (ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) == ext_data(jg)%atm%i_lc_bare_soil) THEN
                     ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu) = 2 ! assume rock in case of bare soil
                   ELSE
                     ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu) = 4 ! otherwise assume sandy loam
                   ENDIF
                 ENDIF

               END DO
             END IF ! ntiles
           ELSE  ! fr_land(jc,jb)<= frlnd_thrhld
             ! measures for land-specific fields that are also defined on non-dominating land points:
             !
             ! for_d, for_e: only accessed via land point index list
             !               -> non-dominant land points do not matter when running without tiles
             ! rootdp, rsmin, lai_mx, plcov_mx, ndvi_max : only accessed via land point index list
             ! ndvi_mrat -> ndviratio : only accessed via land point index list
             ! soiltyp :
             !
             ! glacier fraction
             ext_data(jg)%atm%fr_glac(jc,jb) = 0._wp  ! for frlnd_thrhld=0.5 (i.e. without tiles) this is
                                                      ! identical to what has previously been done within
                                                      ! EXTPAR crosschecks.
             ext_data(jg)%atm%fr_glac_smt(jc,jb) = 0._wp  ! note that this one is used rather than fr_glac !!
           ENDIF



           !
           ! searching for lake-points
           !
           IF (ext_data(jg)%atm%fr_lake(jc,jb) >= frlake_thrhld) THEN
             i_count_flk=i_count_flk+1
             ext_data(jg)%atm%idx_lst_fp(i_count_flk,jb) = jc  ! write index of lake-points
             ext_data(jg)%atm%fp_count(jb) = i_count_flk
             ! set land-cover class
             ext_data(jg)%atm%lc_class_t(jc,jb,isub_lake) = ext_data(jg)%atm%i_lc_water
             ! set also area fractions
             ext_data(jg)%atm%lc_frac_t(jc,jb,isub_lake)  = ext_data(jg)%atm%fr_lake(jc,jb)

             ! set surface area index (needed by turbtran)
             ext_data(jg)%atm%sai_t    (jc,jb,isub_lake)  = c_sea
           ENDIF

           !
           ! searching for sea points
           !
           IF (1._wp-ext_data(jg)%atm%fr_land(jc,jb)-ext_data(jg)%atm%fr_lake(jc,jb) &
             &   >= frsea_thrhld) THEN
             i_count_sea=i_count_sea + 1
             ext_data(jg)%atm%idx_lst_sp(i_count_sea,jb) = jc  ! write index of sea-points
             ext_data(jg)%atm%sp_count(jb) = i_count_sea
             ! set land-cover class
             ext_data(jg)%atm%lc_class_t(jc,jb,isub_water) = ext_data(jg)%atm%i_lc_water
             ! set also area fractions
             ext_data(jg)%atm%lc_frac_t(jc,jb,isub_water)  = 1._wp                        &
               &         -ext_data(jg)%atm%fr_land(jc,jb) - ext_data(jg)%atm%fr_lake(jc,jb)

             ! set surface area index (needed by turbtran)
             ext_data(jg)%atm%sai_t    (jc,jb,isub_water)  = c_sea

             ! set land-cover class for seaice tile
             ! sea-ice and sea have the same land cover class. This is consistent with the
             ! applied GRIB2 tile template, where sea-ice and sea are treated as two
             ! attributes of the same tile. Per definition, different attributes of the
             ! same tile have the same land-cover class
             ext_data(jg)%atm%lc_class_t(jc,jb,isub_seaice) = ext_data(jg)%atm%lc_class_t(jc,jb,isub_water)
           ENDIF

           !
           ! index list for seaice points is generated in mo_nwp_sfc_utils/init_sea_lists
           !
           ! note that in principle, sai_t for seaice should be different from sai_t for
           ! open water points. However, for the time being, sai_t=c_sea is also used
           ! for seaice points.

         END DO ! jc



         ! Inflate dominating land-tile fractions to fr_land or fr_land + fr_lake, depending
         ! on whether a lake tile is present (fr_lake >= frlake_thrhld), or not
         ! (fr_lake < frlake_thrhld).
         IF (ntiles_lnd > 1) THEN
           ! Inflate fractions for land points
           DO ic = 1, ext_data(jg)%atm%lp_count(jb)

             jc = ext_data(jg)%atm%idx_lst_lp(ic,jb)

             ! sum up fractions of dominating land tiles
             sum_frac = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd))

             IF (ext_data(jg)%atm%fr_lake(jc,jb) >= frlake_thrhld) THEN
               ! cell with lake point
               ! inflate dominating land fractions to fr_land
               scalfac = ext_data(jg)%atm%fr_land(jc,jb)/sum_frac
             ELSE
               ! cell without lake point
               ! inflate dominating land fractions to (fr_land + fr_lake)
               scalfac = (ext_data(jg)%atm%fr_land(jc,jb) + ext_data(jg)%atm%fr_lake(jc,jb))/sum_frac
             ENDIF

             ! inflate land fractions
             DO jt = 1, ntiles_total
               ext_data(jg)%atm%lc_frac_t(jc,jb,jt) = ext_data(jg)%atm%lc_frac_t(jc,jb,jt) * scalfac
             ENDDO
           ENDDO  ! ic


           ! Inflate fractions for
           ! - sea-water only points
           ! - lake only points.
           ! As a side effect, fractions for land-only points (with 0<fr_sea<frsea_thrhld)
           ! are also corrected.
           ! Note that, for simplicity, we loop over all points. At mixed land/water points this
           ! should do no harm, since these have already been inflated in the loop above.
           DO jc = i_startidx, i_endidx
             sum_frac = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd)) + &
                        SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,isub_water:isub_lake))

             DO jt = 1, ntiles_total + MIN(2,ntiles_water)
               ext_data(jg)%atm%lc_frac_t(jc,jb,jt) = ext_data(jg)%atm%lc_frac_t(jc,jb,jt) / sum_frac
             ENDDO
           ENDDO  ! jc

           ! Ensure consistency between fr_land and the adjusted sum of the tile fractions
           DO jc = i_startidx, i_endidx
             ext_data(jg)%atm%fr_land(jc,jb)     = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd))
             ext_data(jg)%atm%fr_land_smt(jc,jb) = ext_data(jg)%atm%fr_land(jc,jb)
           ENDDO  ! jc

         ELSE ! overwrite fractional settings over water points if tile approach is turned off
           DO jc = i_startidx, i_endidx
             ext_data(jg)%atm%lc_frac_t(jc,jb,1) = 1._wp
           ENDDO
         ENDIF


         ! Compute inverse of fr_land based on land tile fractions.
         ! Required for proper aggregation of land-only variables
         DO jc = i_startidx, i_endidx
           ext_data(jg)%atm%inv_frland_from_tiles(jc,jb) = 0._wp
           zfr_land = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd))

           IF (zfr_land > 0._wp) THEN
             ext_data(jg)%atm%inv_frland_from_tiles(jc,jb) = 1._wp/zfr_land
           ENDIF
         ENDDO  ! jc


         IF (lsnowtile) THEN ! copy external data fields to snow tile grid points
           DO jt = ntiles_lnd+1, ntiles_total

             jt_in = jt - ntiles_lnd
             ext_data(jg)%atm%lp_count_t(jb,jt)     = ext_data(jg)%atm%lp_count_t(jb,jt_in)
             ext_data(jg)%atm%idx_lst_lp_t(:,jb,jt) = ext_data(jg)%atm%idx_lst_lp_t(:,jb,jt_in)
             !
             ! the following two fields are reset in init_snowtile_lists, but presetting them here
             ! avoids complications in initicon
             ext_data(jg)%atm%gp_count_t(jb,jt)     = ext_data(jg)%atm%gp_count_t(jb,jt_in)
             ext_data(jg)%atm%idx_lst_t(:,jb,jt)    = ext_data(jg)%atm%idx_lst_t(:,jb,jt_in)

!CDIR NODEP
             DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
               jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)
               ext_data(jg)%atm%rootdp_t(jc,jb,jt)   = ext_data(jg)%atm%rootdp_t(jc,jb,jt_in)
               ext_data(jg)%atm%plcov_t(jc,jb,jt)    = ext_data(jg)%atm%plcov_t(jc,jb,jt_in)
               ext_data(jg)%atm%tai_t(jc,jb,jt)      = ext_data(jg)%atm%tai_t(jc,jb,jt_in)
               ext_data(jg)%atm%sai_t(jc,jb,jt)      = ext_data(jg)%atm%sai_t(jc,jb,jt_in)
               ext_data(jg)%atm%eai_t(jc,jb,jt)      = ext_data(jg)%atm%eai_t(jc,jb,jt_in)
               ext_data(jg)%atm%rsmin2d_t(jc,jb,jt)  = ext_data(jg)%atm%rsmin2d_t(jc,jb,jt_in)
               ext_data(jg)%atm%soiltyp_t(jc,jb,jt)  = ext_data(jg)%atm%soiltyp_t(jc,jb,jt_in)
               ext_data(jg)%atm%lc_class_t(jc,jb,jt) = ext_data(jg)%atm%lc_class_t(jc,jb,jt_in)
               ext_data(jg)%atm%lc_frac_t(jc,jb,jt)  = ext_data(jg)%atm%lc_frac_t(jc,jb,jt_in)
             ENDDO

           ENDDO
         ENDIF



         ! Initialize frac_t with lc_frac_t on all static grid points
         ! Recall: frac_t differs from lc_frac_t in the presence of time-dependent sub-lists
         !         (snow tiles or sea ice)
         ! In this case, frac_t refers to the time-dependent sub-tiles.
         ! ** Aggregation operations always must use frac_t **
         DO jt = 1, ntiles_lnd
           DO jc = i_startidx, i_endidx
             ext_data(jg)%atm%frac_t(jc,jb,jt)  = ext_data(jg)%atm%lc_frac_t(jc,jb,jt)
           ENDDO
         ENDDO
         DO jc = i_startidx, i_endidx
           ext_data(jg)%atm%frac_t(jc,jb,isub_water) = ext_data(jg)%atm%lc_frac_t(jc,jb,isub_water)
         ENDDO
         DO jc = i_startidx, i_endidx
           ext_data(jg)%atm%frac_t(jc,jb,isub_lake)  = ext_data(jg)%atm%lc_frac_t(jc,jb,isub_lake)
           !
           ! If tiles are active, ensure consistency between fr_lake and the rescaled tile fraction
           IF (ntiles_lnd > 1) ext_data(jg)%atm%fr_lake(jc,jb) = ext_data(jg)%atm%frac_t(jc,jb,isub_lake)
         ENDDO
         ! frac_t(jc,jb,isub_seaice) is set in init_sea_lists

       END DO !jb
!$OMP END DO


!$OMP SINGLE
       ! Some useful diagnostics
       npoints = SUM(ext_data(jg)%atm%lp_count(i_startblk:i_endblk))
       npoints = global_sum_array(npoints)
       WRITE(message_text,'(a,i3,a,i10)') 'Number of land points in domain',jg,':', npoints
       CALL message('', TRIM(message_text))
       npoints_sea = SUM(ext_data(jg)%atm%sp_count(i_startblk:i_endblk))
       npoints_sea = global_sum_array(npoints_sea)
       WRITE(message_text,'(a,i3,a,i10)') 'Number of sea points in domain',jg,':', npoints_sea
       CALL message('', TRIM(message_text))
       npoints_lake = SUM(ext_data(jg)%atm%fp_count(i_startblk:i_endblk))
       npoints_lake = global_sum_array(npoints_lake)
       WRITE(message_text,'(a,i3,a,i10)') 'Number of lake points in domain',jg,':', npoints_lake
       CALL message('', TRIM(message_text))
       !
       !
       DO i_lu = 1, ntiles_lnd
         npoints = SUM(ext_data(jg)%atm%gp_count_t(i_startblk:i_endblk,i_lu))
         npoints = global_sum_array(npoints)
         WRITE(message_text,'(a,i2,a,i10)') 'Number of points in tile',i_lu,':',npoints
         CALL message('', TRIM(message_text))
       ENDDO
!$OMP END SINGLE NOWAIT


       !
       ! For consistency: remove depth_lk information, where fr_lake < frlake_thrhld.
       ! Boundary interpolation zone of nested domains is explicitly included.
       !
       ! In case of ntiles_lnd > 1, fr_lake ranges from 0<=fr_lake<=1 at nest 
       ! boundaries, whereas at prognostic points fr_lake ranges from 
       ! frlake_thrhld<=fr_lake<=1. I.e. fr_lake consistency adjustment 
       ! was not performed for nest boundary.
       rl_start = 1
       rl_end   = min_rlcell_int

       i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
       i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
       DO jb=i_startblk, i_endblk

         CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

         DO jc = i_startidx, i_endidx
           !
           ! For consistency: remove depth_lk information, where fr_lake=0
           ext_data(jg)%atm%depth_lk(jc,jb) = MERGE(ext_data(jg)%atm%depth_lk(jc,jb), &
             &                                      -1._wp,                           &
             &                                      ext_data(jg)%atm%fr_lake(jc,jb) >= frlake_thrhld)
         ENDDO

       ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END DO  !jg



    ! Diagnose aggregated external parameter fields
    ! (mainly for output purposes)
    !
    CALL diagnose_ext_aggr (p_patch, ext_data)


  END SUBROUTINE init_index_lists



  !-------------------------------------------------------------------------
  !>
  !! Diagnose aggregated external fields
  !!
  !! Aggregated external fields are diagnosed based on tile based external
  !! fields. In addition, fr_land, fr_lake and depth_lk are re-diagnosed,
  !! in order to be consistent with tile-information. Note that the latter 
  !! re-diagnosis has been moved to init_index_lists in order not to 
  !! compromise restart reproducibility.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2013-01-23)
  !!
  SUBROUTINE diagnose_ext_aggr (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    INTEGER  :: jg,jb,jt,ic,jc
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk    !> blocks
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: i_nchdom                !< domain index
    INTEGER  :: i_count
    REAL(wp) :: area_frac

!!$    CHARACTER(len=max_char_length), PARAMETER :: &
!!$      routine = modname//':diagnose_ext_aggr'

    !-------------------------------------------------------------------------

    DO jg = 1, n_dom

      i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

      ! exclude the boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
      i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

      ! Fill nest boundary points of sai with c_sea because the initial call of turbtran
      ! may produce invalid operations otherwise
      ext_data(jg)%atm%sai(:,1:i_startblk) = c_sea

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,ic,i_startidx,i_endidx,i_count,jc,area_frac)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)


        ext_data(jg)%atm%plcov (i_startidx:i_endidx,jb) = 0._wp
        ext_data(jg)%atm%rootdp(i_startidx:i_endidx,jb) = 0._wp
        ext_data(jg)%atm%lai   (i_startidx:i_endidx,jb) = 0._wp
        ext_data(jg)%atm%rsmin (i_startidx:i_endidx,jb) = 0._wp
        ext_data(jg)%atm%tai   (i_startidx:i_endidx,jb) = 0._wp
        ext_data(jg)%atm%eai   (i_startidx:i_endidx,jb) = 0._wp
        ext_data(jg)%atm%sai   (i_startidx:i_endidx,jb) = 0._wp



        DO jt = 1, ntiles_total
          i_count = ext_data(jg)%atm%gp_count_t(jb,jt)
          IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

          DO ic = 1, i_count
            jc = ext_data(jg)%atm%idx_lst_t(ic,jb,jt)

            ! note that frac_t must be re-scaled such that sum(frac_t(1:ntiles_lnd)) = 1
            ! therefore we multiply by inv_frland_from_tiles
            area_frac = ext_data(jg)%atm%frac_t(jc,jb,jt)           &
              &       * ext_data(jg)%atm%inv_frland_from_tiles(jc,jb)

            ! plant cover (aggregated)
            ext_data(jg)%atm%plcov(jc,jb) = ext_data(jg)%atm%plcov(jc,jb)       &
              &              + ext_data(jg)%atm%plcov_t(jc,jb,jt) * area_frac

            ! root depth (aggregated)
            ext_data(jg)%atm%rootdp(jc,jb) = ext_data(jg)%atm%rootdp(jc,jb)     &
              &              + ext_data(jg)%atm%rootdp_t(jc,jb,jt) * area_frac

            ! surface area index (aggregated)
            ext_data(jg)%atm%lai(jc,jb) = ext_data(jg)%atm%lai(jc,jb)           &
              &              + ( ext_data(jg)%atm%tai_t(jc,jb,jt)                &
              &              /(ext_data(jg)%atm%plcov_t(jc,jb,jt)+dbl_eps) * area_frac )

            ! evaporative soil area index (aggregated)
            ext_data(jg)%atm%eai(jc,jb) = ext_data(jg)%atm%eai(jc,jb)           &
              &              +  ext_data(jg)%atm%eai_t(jc,jb,jt) * area_frac

            ! transpiration area index (aggregated)
            ext_data(jg)%atm%tai(jc,jb) = ext_data(jg)%atm%tai(jc,jb)           &
              &              +  ext_data(jg)%atm%tai_t(jc,jb,jt) * area_frac

            ! minimal stomata resistance (aggregated)
            ext_data(jg)%atm%rsmin(jc,jb) = ext_data(jg)%atm%rsmin(jc,jb)       &
              &              + ext_data(jg)%atm%rsmin2d_t(jc,jb,jt) * area_frac

          ENDDO  !ic
        ENDDO  !jt


        ! aggregate fields with water tiles
        DO jt = 1,ntiles_total + ntiles_water
          DO jc = i_startidx, i_endidx

            area_frac = ext_data(jg)%atm%frac_t(jc,jb,jt)

            ! surface area index (aggregated)
            ext_data(jg)%atm%sai(jc,jb) = ext_data(jg)%atm%sai(jc,jb)           &
              &             +  ext_data(jg)%atm%sai_t(jc,jb,jt) * area_frac
          ENDDO  ! jc
        ENDDO  !jt

      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    ENDDO  !jg

  END SUBROUTINE diagnose_ext_aggr



  !-------------------------------------------------------------------------
  !>
  !! Get interpolated field from monthly mean climatology
  !!
  !! Get interpolated field from monthly mean climatology. A linear interpolation
  !! in time between successive months is performed, assuming that the monthly field
  !! applies to the 15th of the month.
  !!
  !! @par Revision History
  !! Initial revision by Juergen Helmert, DWD (2012-04-17)
  !! Modification by Daniel Reinert, DWD (2013-05-03)
  !! Generalization to arbitrary monthly mean climatologies
  !!
  SUBROUTINE interpol_monthly_mean(p_patch, mtime_date, monthly_means, out_field)

    TYPE(datetime),   POINTER      :: mtime_date
    TYPE(t_patch),     INTENT(IN)  :: p_patch
    REAL(wp),          INTENT(IN)  :: monthly_means(:,:,:)  ! monthly mean climatology
    REAL(wp),          INTENT(OUT) :: out_field(:,:)        ! interpolated output field

    INTEGER                             :: jc, jb               !< loop index
    INTEGER                             :: i_startblk, i_endblk
    INTEGER                             :: rl_start, rl_end
    INTEGER                             :: i_startidx, i_endidx
    INTEGER                             :: mo1, mo2             !< nearest months
    REAL(wp)                            :: zw1, zw2
    TYPE(datetime), POINTER             :: mtime_hour

    TYPE(t_time_interpolation_weights)  :: current_time_interpolation_weights
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dtime_string
    
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//': interpol_monthly_mean'

    !---------------------------------------------------------------
    ! Find the 2 nearest months mo1, mo2 and the weights zw1, zw2
    ! to the actual date and time

    mtime_hour => newDatetime(mtime_date)
    mtime_hour%time%minute = 0
    mtime_hour%time%second = 0
    mtime_hour%time%ms     = 0     
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_hour)
    call deallocateDatetime(mtime_hour)    
    mo1 = current_time_interpolation_weights%month1
    mo2 = current_time_interpolation_weights%month2
    zw1 = current_time_interpolation_weights%weight1
    zw2 = current_time_interpolation_weights%weight2

    ! consistency check
    IF ((MIN(mo1,mo2) < 1) .OR. (MAX(mo1,mo2) > SIZE(monthly_means,3))) THEN
      CALL datetimeToString(mtime_date, dtime_string)
      CALL message('','Result of call to calculate_time_interpolation_weights:')
      WRITE (message_text,'(a,a)')      '   mtime_date = ', TRIM(dtime_string)
      CALL message('', message_text)
      WRITE (message_text,'(a,i2.2)')   '   mo1        = ', mo1
      CALL message('', message_text)
      WRITE (message_text,'(a,i2.2)')   '   mo2        = ', mo2
      CALL message('', message_text)
      WRITE (message_text,'(a,f25.15)') '   weight1    = ', zw1
      CALL message('', message_text)
      WRITE (message_text,'(a,f25.15)') '   weight2    = ', zw2
      CALL message('', message_text)      
      CALL finish(routine, "Error!")
    END IF

    ! exclude the boundary interpolation zone of nested domains
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! Get interpolated field
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb=i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
         & i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx
        out_field(jc,jb) = zw1*monthly_means(jc,jb,mo1) &
          &              + zw2*monthly_means(jc,jb,mo2)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE interpol_monthly_mean

END MODULE mo_ext_data_init

