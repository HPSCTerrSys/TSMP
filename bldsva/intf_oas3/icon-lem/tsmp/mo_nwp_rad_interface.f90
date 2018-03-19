!>
!! This module is the interface between nwp_nh_interface to the radiation schemes
!! (RRTM or Ritter-Geleyn).
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach
!!
!! @par Revision History
!! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_nwp_rad_interface

  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma, parallel_radiation_mode
  USE mo_impl_constants,       ONLY: max_char_length, MODIS 
  USE mo_kind,                 ONLY: wp
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_radiation_config,     ONLY: albedo_type
  USE mo_radiation,            ONLY: pre_radiation_nwp_steps
  USE mo_nwp_rrtm_interface,   ONLY: nwp_rrtm_radiation,             &
    &                                nwp_rrtm_radiation_reduced,     &
    &                                nwp_rrtm_radiation_repartition, &
    &                                nwp_ozon_aerosol
  USE mo_nwp_rg_interface,     ONLY: nwp_rg_radiation,               &
    &                                nwp_rg_radiation_reduced
  USE mo_albedo,               ONLY: sfc_albedo, sfc_albedo_modis
  USE mtime,                   ONLY: datetime
  
  IMPLICIT NONE

  PRIVATE



  PUBLIC :: nwp_radiation
  

 CONTAINS
  
  !---------------------------------------------------------------------------------------
  !>
  !! This subroutine is the interface between nwp_nh_interface to the radiation schemes.
  !! Depending on inwp_radiation, it can call RRTM (1) or Ritter-Geleyn (2).
  !!
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_radiation ( lredgrid, p_sim_time, mtime_datetime, pt_patch,pt_par_patch, &
    & ext_data, lnd_diag, pt_prog, pt_diag, prm_diag, lnd_prog, wtr_prog )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
      &  routine = 'mo_nwp_rad_interface:nwp_radiation'
    
    LOGICAL,                 INTENT(in)    :: lredgrid        !< use reduced grid for radiation

    REAL(wp),                INTENT(in)    :: p_sim_time

    TYPE(datetime), POINTER, INTENT(in)    :: mtime_datetime
    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_patch     !<grid/patch info.
    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_external_data),   INTENT(inout) :: ext_data
    TYPE(t_lnd_diag),        INTENT(in)    :: lnd_diag   !<diag vars for sfc
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog    !<the prognostic variables
    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag    !<the diagnostic variables
    TYPE(t_nwp_phy_diag),    INTENT(inout) :: prm_diag
    TYPE(t_lnd_prog),        INTENT(inout) :: lnd_prog   ! time level new
    TYPE(t_wtr_prog),        INTENT(in)    :: wtr_prog   ! time level new

    REAL(wp) :: &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    
    INTEGER :: jg, irad

    REAL(wp):: zsct        ! solar constant (at time of year)
    REAL(wp):: cosmu0_dark ! minimum cosmu0, for smaller values no shortwave calculations



    ! patch ID
    jg = pt_patch%id



    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    SELECT CASE (atm_phy_nwp_config(jg)%inwp_radiation )
    CASE (1, 3)
      ! RRTM
      ! In radiative transfer routine RRTM skips all points with cosmu0<=0. That's why 
      ! points to be skipped need to be marked with a value <=0
      cosmu0_dark = -1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    CASE (2)
      ! Ritter-Geleyn
      ! Skipping of points is performed on block- rather than cell-level. I.e. if a block 
      ! contains at least 1 point with cosmu0>1.E-8, radiatve transfer is computed for 
      ! the entire block. Therefore cosmu0_dark = -1.e-9_wp does not work here (crashes).
      ! For all points cosmu0 must be <0.
      cosmu0_dark =  1.e-9_wp   ! minimum cosmu0, for smaller values no shortwave calculations
    END SELECT


    ! Calculation of zenith angle optimal during dt_rad.
    ! (For radheat, actual zenith angle is calculated separately.)
    CALL pre_radiation_nwp_steps (                        &
      & kbdim        = nproma,                            & !in
      & cosmu0_dark  = cosmu0_dark,                       & !in
      & p_inc_rad    = atm_phy_nwp_config(jg)%dt_rad,     & !in
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_fastphy, & !in
      & p_sim_time   = p_sim_time,                        & !in
      & pt_patch     = pt_patch,                          & !in
      & zsmu0        = prm_diag%cosmu0(:,:),              & !out
      & zsct         = zsct                               ) !out, optional



#ifndef COUP_OAS_ICON
    ! Compute tile-based and aggregated surface-albedo
    !
    IF ( albedo_type == MODIS ) THEN
      ! MODIS albedo
      CALL sfc_albedo_modis(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)
    ELSE
      ! albedo based on tabulated bare soil values
      CALL sfc_albedo(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)
    ENDIF
#else
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    i_startblk = p_patch(1)%cells%start_blk(rl_start, 1)
    i_endblk   = p_patch(1)%cells%end_blk(rl_end,MAX(1,p_patch(1)%n_childdom))
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch(1), jb, i_startblk, i_endblk, i_startidx, &
        i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        prm_diag%albdif(jc,jb) = oas_rcv_field_icon(2,jc,jb)
        prm_diag%albvisdif(jc,jb) =  prm_diag%albdif(jc,jb)
        prm_diag%albnirdif(jc,jb) = prm_diag%albdif(jc,jb)
        prm_diag%albnirdir(jc,jb) = oas_rcv_field_icon(3,jc,jb)
        prm_diag%albvisdir(jc,jb) = prm_diag%albnirdir(jc,jb)
      END DO
    END DO
#endif



    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------
    !
    SELECT CASE (atm_phy_nwp_config(jg)%inwp_radiation)
    CASE (1, 3) ! RRTM / PSRAD

      irad = atm_phy_nwp_config(jg)%inwp_radiation

      CALL nwp_ozon_aerosol ( p_sim_time, mtime_datetime, pt_patch, ext_data, &
        & pt_diag, prm_diag, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5 )
    
      IF ( .NOT. lredgrid ) THEN

        SELECT CASE(parallel_radiation_mode(jg))
        CASE(1) 
          CALL nwp_rrtm_radiation_repartition ( mtime_datetime, pt_patch, ext_data, &
            & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                    &
            & pt_diag, prm_diag, lnd_prog )
          
        CASE default
          CALL nwp_rrtm_radiation ( mtime_datetime, pt_patch, ext_data, &
            & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,        &
            & pt_diag, prm_diag, lnd_prog, irad )

        END SELECT
       
      ELSE 

        CALL nwp_rrtm_radiation_reduced ( mtime_datetime, pt_patch,pt_par_patch, ext_data, &
          & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                             &
          & pt_diag, prm_diag, lnd_prog, irad )
          
      ENDIF

    CASE (2) ! Ritter-Geleyn

      IF (.NOT. lredgrid) THEN
        CALL nwp_rg_radiation ( p_sim_time, mtime_datetime, pt_patch, &
          & ext_data,pt_prog,pt_diag,prm_diag, lnd_prog, zsct )
      ELSE
        CALL nwp_rg_radiation_reduced ( p_sim_time, mtime_datetime, pt_patch,pt_par_patch, &
          & ext_data, pt_prog, pt_diag, prm_diag, lnd_prog, zsct )
      ENDIF

    END SELECT ! inwp_radiation

  END SUBROUTINE nwp_radiation


END MODULE mo_nwp_rad_interface

