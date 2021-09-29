!>
!! This module contains various interfaces to the Ritter-Geleyn radiation scheme.
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

MODULE mo_nwp_rg_interface

  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_exception,            ONLY: message
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma, p_test_run

  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi
  USE mo_impl_constants,       ONLY: min_rlcell_int  ! io3_ape 
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c, grf_ovlparea_start_c
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog
  USE mo_model_domain,         ONLY: t_patch, p_patch_local_parent
  USE mo_phys_nest_utilities,  ONLY: upscale_rad_input_rg, downscale_rad_output_rg
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_radiation_config,     ONLY: vmr_co2
  USE mo_radiation_rg,         ONLY: fesft
  USE mo_nwp_rrtm_interface,   ONLY: nwp_ozon_aerosol
  USE mo_satad,                ONLY: qsat_rho
  USE mtime,                   ONLY: datetime

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: nwp_rg_radiation
  PUBLIC :: nwp_rg_radiation_reduced

 CONTAINS
  


  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_rg_radiation ( p_sim_time, mtime_datetime, pt_patch, &
    & ext_data,pt_prog,pt_diag,prm_diag,lnd_prog,zsct )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rg_interface:'

!    REAL(wp), PARAMETER::  &
!      & zqco2 = 0.5014E-03_wp*353.9_wp/330._wp ! CO2 (mixing ratio 353.9 ppm (like vmr_co2))

    REAL(wp),                INTENT(in)    :: p_sim_time
    
    TYPE(datetime), POINTER, INTENT(in)    :: mtime_datetime
    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_patch     !<grid/patch info.
    TYPE(t_external_data),   INTENT(inout) :: ext_data
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog     !<the prognostic variables
    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),    INTENT(inout) :: prm_diag
    TYPE(t_lnd_prog),        INTENT(inout) :: lnd_prog
    REAL(wp),                INTENT(in)    :: zsct        ! solar constant (at time of year)

    REAL(wp):: zi0        (nproma)  !< solar incoming radiation at TOA   [W/m2]
    ! for Ritter-Geleyn radiation:
    REAL(wp):: zqco2
    REAL(wp):: zsqv (nproma,pt_patch%nlev) !< saturation water vapor mixing ratio
    REAL(wp):: zfls (nproma,pt_patch%nlevp1)
    REAL(wp):: zduo3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zduco2(nproma,pt_patch%nlev)
    REAL(wp):: alb_ther(nproma)
    LOGICAL :: lo_sol (nproma)
    LOGICAL :: losol

    ! Local scalars:
    INTEGER:: jc,jk,jb
    INTEGER:: nlev, nlevp1      !< number of full and half levels

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index

    i_nchdom  = MAX(1,pt_patch%n_childdom)

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    ! CO2 help variable for Ritter-Geleyn scheme
    zqco2 = 0.5014E-03_wp*vmr_co2/330.e-6_wp

    CALL nwp_ozon_aerosol ( p_sim_time, mtime_datetime, pt_patch, ext_data, &
      & pt_diag, prm_diag, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, zduo3 )
    


    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

    IF (msg_level >= 15) &
      & CALL message('mo_nwp_rg_interface:nwp_rg_radiation', 'RG radiation on full grid')

    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO &
#ifndef __xlC__
!$OMP& SCHEDULE(guided), &
#endif
!$OMP& PRIVATE(jb,i_startidx,i_endidx,jc,jk,zi0,losol,lo_sol,&
!$OMP& alb_ther,zfls,zduco2,zsqv)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1,nlev
        DO jc = i_startidx,i_endidx
          zsqv (jc,jk) = qsat_rho(pt_diag%temp(jc,jk,jb),pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDDO

      ! geographical dependent thermal albedo
      alb_ther(i_startidx:i_endidx) = 1._wp-ext_data%atm%emis_rad(i_startidx:i_endidx,jb)

      prm_diag%tsfctrad(i_startidx:i_endidx,jb) = lnd_prog%t_g(i_startidx:i_endidx,jb)

      ! CO2 (mixing ratio 353.9 ppm as vmr_co2)
      DO jk = 1,nlev
        DO jc = i_startidx,i_endidx
          zduco2(jc,jk) = zqco2 * pt_diag%dpres_mc (jc,jk,jb)
        ENDDO
      ENDDO

      ! Switch off solar radiation calculations where sun is below horizon:
      lo_sol(i_startidx:i_endidx) = prm_diag%cosmu0(i_startidx:i_endidx,jb) > 1.e-8_wp
      losol = ANY(lo_sol(i_startidx:i_endidx))

      CALL fesft ( &
                                !  Input:
        & pti = pt_diag%temp_ifc (:,:,jb) , &! Temperature at layer boundaries
        & pdp = pt_diag%dpres_mc (:,:,jb), &! pressure thickness
        & pclc_in= prm_diag%clc  (:,:,jb)   , &
        & pqv = prm_diag%tot_cld(:,:,jb,iqv), &
        & pqvs = zsqv, &!saturation water vapor
        & pqcwc = prm_diag%tot_cld    (:,:,jb,iqc) ,&
        & pqiwc = prm_diag%tot_cld    (:,:,jb,iqi) ,&
        & pduco2 = zduco2, &! layer CO2 content
        & pduo3  = zduo3(:,:,jb),&! layer O3 content
        & paeq1 = zaeq1(:,:,jb), &
        & paeq2 = zaeq2(:,:,jb),&
        & paeq3 = zaeq3(:,:,jb),&
        & paeq4 = zaeq4(:,:,jb),&
        & paeq5 = zaeq5(:,:,jb),&
        & papre_in =  pt_diag%pres_sfc (:,jb), & ! Surface pressure
        & psmu0 = prm_diag%cosmu0 (:,jb) , & ! Cosine of zenith angle
        & palso = prm_diag%albdif(:,jb), & ! solar surface albedo
        & palth = alb_ther, & ! thermal surface albedo
        & psct = zsct, &! solar constant (at time of year)
        & kig1s = 1 ,&
        & kig1e = nproma , &
        & ki3s = 1, &
        & ki3e = nlev,&
        & ki1sc= i_startidx, &
        & ki1ec= i_endidx, &
        & lsolar = losol, &! control switch for solar calculations
        !          & lsolar = .TRUE., &! control switch for solar calculations
        & lthermal =.TRUE., &
        & lcrf = .FALSE., &! control switch for cloud-free calcul.
                                ! Output:
        & pflt  = prm_diag%lwflxall(:,:,jb),& !Thermal radiative fluxes at each layer boundary
        & pfls  = zfls &! solar radiative fluxes at each layer boundary
        & )

      zi0 (i_startidx:i_endidx) = 1._wp / (prm_diag%cosmu0(i_startidx:i_endidx,jb) * zsct)
      ! compute sw transmissivity trsolall from sw fluxes
      DO jk = 1,nlevp1
        DO jc = i_startidx,i_endidx
          prm_diag%trsolall(jc,jk,jb) &
            = MERGE(0.0_wp, zfls(jc,jk) * zi0(jc), prm_diag%cosmu0(jc,jb) < 1.e-8_wp)
        ENDDO
      ENDDO

    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE nwp_rg_radiation
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_rg_radiation_reduced ( p_sim_time, mtime_datetime, pt_patch,pt_par_patch, &
    &                                   ext_data,pt_prog,pt_diag,prm_diag, &
    &                                   lnd_prog,zsct )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rg_interface:'

!    REAL(wp), PARAMETER::  &
!      & zqco2 = 0.5014E-03_wp*353.9_wp/330._wp ! CO2 (mixing ratio 353.9 ppm (like vmr_co2))

    
    REAL(wp),                INTENT(in) :: p_sim_time

    TYPE(datetime), POINTER, INTENT(in) :: mtime_datetime 
    TYPE(t_patch), TARGET,   INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_patch), TARGET,   INTENT(in) :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_external_data),   INTENT(inout):: ext_data
    TYPE(t_nh_prog), TARGET, INTENT(inout)  :: pt_prog     !<the prognostic variables
    TYPE(t_nh_diag), TARGET, INTENT(inout)  :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag),    INTENT(inout):: prm_diag
    TYPE(t_lnd_prog),        INTENT(inout):: lnd_prog
    REAL(wp),                INTENT(in)  :: zsct        ! solar constant (at time of year)

    ! For radiation on reduced grid
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_cosmu0   (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albdif(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albeff(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_albefffac(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_tsfc     (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_o3       (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_tot_cld  (:,:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_clc      (:,:,:)
    ! Output fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_lwflxall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_trsolall (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_fls (:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_flsp (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_flsd (:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_flsu (:,:)
    ! Pointer to parent patach or local parent patch for reduced grid
    TYPE(t_patch), POINTER       :: ptr_pp

    REAL(wp):: zi0        (nproma)  !< solar incoming radiation at TOA   [W/m2]
    ! for Ritter-Geleyn radiation:
    REAL(wp):: zqco2
    REAL(wp):: zsqv (nproma,pt_patch%nlev,pt_patch%nblks_c) !< saturation water vapor mixing ratio
    REAL(wp):: zduo3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: zduco2(nproma,pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp):: alb_ther    (nproma,pt_patch%nblks_c) !!
    LOGICAL :: lo_sol (nproma)
    LOGICAL :: losol
    ! For Ritter-Geleyn radiation on reduced grid additionally
    ! These fields need to be allocatable because they have different dimensions for
    ! the global grid and nested grids, and for runs with/without MPI parallelization
    ! Input fields
    REAL(wp), ALLOCATABLE, TARGET:: zrg_alb_ther(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_pres_sfc(:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_temp_ifc(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_dpres_mc(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_sqv(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_duco2(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq1(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq2(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq3(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq4(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET:: zrg_aeq5(:,:,:)

    ! Local scalars:
    INTEGER:: jc,jk,jb
    INTEGER:: jg                                !domain id
    INTEGER:: nlev, nlevp1, nlev_rg, nlevp1_rg  !< number of full and half levels
    INTEGER:: nblks_par_c                       !nblks for reduced grid

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index
    INTEGER:: i_chidx

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    jg        = pt_patch%id

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    ! CO2 help variable for Ritter-Geleyn scheme
    zqco2 = 0.5014E-03_wp*vmr_co2/330.e-6_wp


    CALL nwp_ozon_aerosol ( p_sim_time, mtime_datetime, pt_patch, ext_data, &
      & pt_diag, prm_diag, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, zduo3 )


    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------

      ! section for computing radiation on reduced grid

    IF (p_test_run) THEN
      prm_diag%lwflxall(:,:,:) = 0._wp
      prm_diag%trsolall(:,:,:) = 0._wp
    ENDIF

    IF (msg_level >= 15) &
      &  CALL message('mo_nwp_rg_interface:nwp_rg_radiation_reduced', &
      &               'Ritter-Geleyn radiation on reduced grid')

    i_chidx     =  pt_patch%parent_child_index

    IF (jg == 1) THEN
      ptr_pp => pt_par_patch
      nblks_par_c = pt_par_patch%nblks_c

      ! number of vertical levels
      ! ** for the time being, the radiation grid is assumed to have the same
      !    levels as the main grid **
      ! nlev   = ptr_pp%nlev
      ! nlevp1 = ptr_pp%nlevp1
    ELSE ! Nested domain with MPI parallelization
      ptr_pp      => p_patch_local_parent(jg)
      nblks_par_c =  ptr_pp%nblks_c

      ! number of vertical levels
      ! ** for the time being, the radiation grid is assumed to have the same
      !    levels as the main grid **
      ! nlev   = ptr_pp%nlev
      ! nlevp1 = ptr_pp%nlevp1
    ENDIF

    ! Add extra layer for atmosphere above model top if requested
    IF (atm_phy_nwp_config(jg)%latm_above_top) THEN
      nlev_rg   = nlev + 1
      nlevp1_rg = nlevp1 + 1
    ELSE
      nlev_rg   = nlev
      nlevp1_rg = nlevp1
    ENDIF

    ALLOCATE (                                       &
      zrg_cosmu0   (nproma,          nblks_par_c),   &
      zrg_albdif   (nproma,          nblks_par_c),   &
      zrg_albeff   (nproma,          nblks_par_c),   &
      zrg_albefffac(nproma,          nblks_par_c),   &
      zrg_alb_ther (nproma,          nblks_par_c),   &
      zrg_pres_sfc (nproma,          nblks_par_c),   &
      zrg_tsfc     (nproma,          nblks_par_c),   &
      zrg_temp_ifc (nproma,nlevp1_rg,nblks_par_c),   &
      zrg_dpres_mc (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_sqv      (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_duco2    (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_o3       (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq1     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq2     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq3     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq4     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_aeq5     (nproma,nlev_rg  ,nblks_par_c),   &
      zrg_tot_cld  (nproma,nlev_rg  ,nblks_par_c,3), &
      zrg_clc      (nproma,nlev_rg  ,nblks_par_c)  , &
      zrg_fls      (nproma,nlevp1_rg,nblks_par_c),   &
      zrg_flsp     (nproma,          nblks_par_c),   &
      zrg_flsd     (nproma,          nblks_par_c),   &
      zrg_flsu     (nproma,          nblks_par_c),   &
      zrg_lwflxall (nproma,nlevp1_rg,nblks_par_c),   &
      zrg_trsolall (nproma,nlevp1_rg,nblks_par_c)  )


    rl_start = 1  !DR 3 seems to be sufficient 
                  ! SR radiation is not set up to handle boundaries of nested domains
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
    
    ! *** this parallel section will be removed later once real data are
    !     are available as input for radiation ***
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
    !
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1,nlev
        DO jc = i_startidx,i_endidx
          zsqv (jc,jk,jb) = qsat_rho(pt_diag%temp(jc,jk,jb),pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDDO

      ! geographical dependent thermal albedo
      alb_ther(i_startidx:i_endidx,jb) = 1._wp-ext_data%atm%emis_rad(i_startidx:i_endidx,jb)

      prm_diag%tsfctrad(i_startidx:i_endidx,jb) = lnd_prog%t_g(i_startidx:i_endidx,jb)

      ! CO2 (mixing ratio 353.9 ppm as vmr_co2)
      DO jk = 1,nlev
        DO jc = i_startidx,i_endidx
          zduco2(jc,jk,jb) = zqco2 * pt_diag%dpres_mc (jc,jk,jb)
        ENDDO
      ENDDO

    ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

    CALL upscale_rad_input_rg( pt_patch%id, pt_par_patch%id,  nlev_rg, nlevp1_rg,            &
      & prm_diag%cosmu0, prm_diag%albdif, alb_ther, pt_diag%temp_ifc,                        &
      & pt_diag%dpres_mc, prm_diag%tot_cld, prm_diag%clc, zsqv ,zduco2, zduo3,               &
      & zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,pt_diag%pres_sfc,pt_diag%pres_ifc,                     &
      & zrg_cosmu0, zrg_albdif, zrg_alb_ther, zrg_temp_ifc, zrg_dpres_mc,                    &
      & zrg_tot_cld, zrg_clc,zrg_sqv ,zrg_duco2, zrg_o3,                                     &
      & zrg_aeq1,zrg_aeq2,zrg_aeq3,zrg_aeq4,zrg_aeq5,zrg_pres_sfc     )

    rl_start = grf_ovlparea_start_c
    rl_end   = min_rlcell_int

    i_startblk = ptr_pp%cells%start_blk(rl_start,i_chidx)
    i_endblk   = ptr_pp%cells%end_blk(rl_end,i_chidx)


!$OMP PARALLEL
#ifdef __xlC__
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,lo_sol,losol,zi0)
#else
!$OMP DO PRIVATE(jb,jk,i_startidx,i_endidx,lo_sol,losol,zi0) ,SCHEDULE(guided)
#endif
    !
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
        &                         i_startidx, i_endidx, rl_start, rl_end)

      ! Switch off solar radiation calculations where sun is below horizon:
      WHERE ( zrg_cosmu0(i_startidx:i_endidx,jb) > 1.e-8_wp ) !zepmu0 )
        lo_sol(i_startidx:i_endidx) = .TRUE.
      ELSEWHERE
        lo_sol(i_startidx:i_endidx) = .FALSE.
      END WHERE
      losol = ANY(lo_sol(i_startidx:i_endidx))

      CALL fesft ( &
                                         !  Input:
        & pti = zrg_temp_ifc (:,:,jb) , &! Temperature at layer boundaries
        & pdp = zrg_dpres_mc (:,:,jb), &! pressure thickness
        & pclc_in= zrg_clc   (:,:,jb) , &
        & pqv = zrg_tot_cld(:,:,jb,iqv), &
        & pqvs = zrg_sqv(:,:,jb), &!saturation water vapor
        & pqcwc = zrg_tot_cld    (:,:,jb,iqc) ,&
        & pqiwc = zrg_tot_cld    (:,:,jb,iqi) ,&
        & pduco2 = zrg_duco2 (:,:,jb), &! layer CO2 content
        & pduo3  = zrg_o3 (:,:,jb),&! layer O3 content
        & paeq1 = zrg_aeq1(:,:,jb), &
        & paeq2 = zrg_aeq2(:,:,jb),&
        & paeq3 = zrg_aeq3(:,:,jb),&
        & paeq4 = zrg_aeq4(:,:,jb),&
        & paeq5 = zrg_aeq5(:,:,jb),&
        & papre_in = zrg_pres_sfc (:,jb), & ! Surface pressure
        & psmu0 = zrg_cosmu0 (:,jb) , & ! Cosine of zenith angle
        & palso = zrg_albdif(:,jb),   & ! solar surface albedo
        & palth = zrg_alb_ther(:,jb), & ! thermal surface albedo
        & psct = zsct, &! solar constant (at time of year)
        & kig1s = 1 ,&
        & kig1e = nproma , &
        & ki3s = 1, &
        & ki3e = nlev_rg,&
        & ki1sc= i_startidx, &
        & ki1ec= i_endidx, &
        & lsolar = losol, &! control switch for solar calculations
        !          & lsolar = .TRUE., &! control switch for solar calculations
        & lthermal =.TRUE., &
        & lcrf = .FALSE., &! control switch for cloud-free calcul.
                                ! Output:
        & pflt  = zrg_lwflxall(:,:,jb) ,& ! Thermal radiative fluxes at each layer boundary
        & pfls  = zrg_fls  (:,:,jb),  &! solar radiative fluxes at each layer boundary
        & pflsp = zrg_flsp (:,jb), &
        & pflsd = zrg_flsd (:,jb), &
        & pflsu = zrg_flsu (:,jb) &
        & )

      zi0 (i_startidx:i_endidx) = zrg_cosmu0(i_startidx:i_endidx,jb) * zsct
      ! compute sw transmissivity trsolall from sw fluxes
      DO jk = 1,nlevp1_rg
        DO jc = i_startidx,i_endidx
          ! This is needed to avoid false synchronization errors
          IF (zrg_cosmu0(jc,jb) < 1.e-8_wp) THEN
            zrg_trsolall(jc,jk,jb) = 0._wp
          ELSE
            zrg_trsolall(jc,jk,jb) = zrg_fls(jc,jk,jb) / zi0(jc)
          ENDIF
        ENDDO
      ENDDO

      DO jc = i_startidx,i_endidx
        
        IF (zrg_cosmu0(jc,jb) < 1.e-8_wp) THEN
          zrg_albefffac(jc,jb) = 1._wp
          zrg_albeff(jc,jb)    = 0.5_wp
        ELSE
          zrg_albeff(jc,jb) = zrg_flsu(jc,jb) / ( zrg_flsp(jc,jb) + zrg_flsd(jc,jb) )
          zrg_albefffac(jc,jb) = zrg_albeff(jc,jb) / zrg_albdif(jc,jb)
        ENDIF
        
        zrg_tsfc(jc,jb) = zrg_temp_ifc(jc,nlevp1_rg,jb)
        
      ENDDO

      ! to avoid division by zero inside downscale_rad_output_rg
      zrg_flsp(i_startidx:i_endidx,jb)=MAX(zrg_flsp(i_startidx:i_endidx,jb),1.e-9_wp)
      zrg_flsd(i_startidx:i_endidx,jb)=MAX(zrg_flsd(i_startidx:i_endidx,jb),1.e-9_wp)

    ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL


    CALL downscale_rad_output_rg(          &
      & jg           = pt_patch%id,        &
      & jgp          = pt_par_patch%id,    &
      & nlev_rg      = nlev_rg,            &
      & rg_lwflxall  = zrg_lwflxall,       &
      & rg_trsolall  = zrg_trsolall,       &
      & tsfc_rg      = zrg_tsfc,           &
      & albeff_rg    = zrg_albeff,         &
      & albefffac_rg = zrg_albefffac,      &
      & flsp_rg      = zrg_flsp,           &
      & flsd_rg      = zrg_flsd,           & 
      & alb_ther_rg  = zrg_alb_ther,       &
      & cosmu0_rg    = zrg_cosmu0,         &
      & tot_cld_rg   = zrg_tot_cld,        &
      & dpres_mc_rg  = zrg_dpres_mc,       &
      & pres_sfc_rg  = zrg_pres_sfc,       &
      & tsfc         = prm_diag%tsfctrad,  &
      & albdif       = prm_diag%albdif,    &
      & zsct         = zsct,               &
      & lwflxall     = prm_diag%lwflxall,  &
      & trsolall     = prm_diag%trsolall )

    DEALLOCATE (zrg_cosmu0,zrg_tsfc,zrg_albdif,zrg_alb_ther,                &
      & zrg_pres_sfc,zrg_temp_ifc,zrg_dpres_mc,zrg_sqv,zrg_duco2,zrg_o3,    &
      & zrg_aeq1,zrg_aeq2,zrg_aeq3,zrg_aeq4,zrg_aeq5,                       &
      & zrg_tot_cld, zrg_clc, zrg_fls, zrg_flsp, zrg_flsd, zrg_flsu,        &
      & zrg_lwflxall,zrg_trsolall)

  END SUBROUTINE nwp_rg_radiation_reduced
  

END MODULE mo_nwp_rg_interface

