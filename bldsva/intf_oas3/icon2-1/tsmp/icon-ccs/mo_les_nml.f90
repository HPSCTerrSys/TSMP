!>
!! Contains the setup of variables related to large eddy simulation setup
!!        
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_les_nml

  USE mo_les_config,          ONLY: les_config
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_les_namelist
  PUBLIC :: turb_profile_list, turb_tseries_list

  REAL(wp) :: sst        ! prescribed SST
  REAL(wp) :: psfc       ! prescribed surface pressure
  REAL(wp) :: shflx      ! prescribed sensible heat flux (Km/s)
  REAL(wp) :: lhflx      ! prescribed latent heat flux   (Km/s)
  INTEGER  :: isrfc_type ! 1=fixed sst, 2=fixed flux, 3=fixed buyancy flux

  REAL(wp) :: ufric      ! friction velocity
 
  LOGICAL  :: is_dry_cbl  !special case for CBL testcase
 
  !For isrf_type==3
  REAL(wp) :: bflux      !Buoyancy flux
  REAL(wp) :: tran_coeff !Surface transfer coefficient in units of velocity (m/s)

  !Some parameters
  REAL(wp) :: smag_constant
  REAL(wp) :: turb_prandtl 
  REAL(wp) :: km_min         !min mass weighted turbulent viscosity 
  REAL(wp) :: max_turb_scale !max turbulence length scale
  REAL(wp) :: min_sfc_wind  !min sfc wind in free convection limit

  !Scheme for vertical discretization
  INTEGER :: vert_scheme_type !1=explicit, 2=implicit

  !Parameters for additional diagnostic output
  LOGICAL  :: ldiag_les_out                    !.TRUE. to turn it on
  REAL(wp) :: avg_interval_sec, sampl_freq_sec !averaging and sampling time 
  CHARACTER(LEN=7) :: turb_tseries_list(19), turb_profile_list(44) !list of variables  
  CHARACTER(MAX_CHAR_LENGTH) :: expname        !name of experiment for naming the file
  LOGICAL  :: les_metric

  ! SBr, CHa: vars for idealized profile of soil moisture
  REAL(wp) :: wso_gradient, wso_mean(8), wso_max(8), wso_min(8), wso_spread
  INTEGER  :: wso_fun !1=homo, 2=function2, 3=function3
  INTEGER  :: wso_cyc !wso_cyc=1,2,3,4,6

  NAMELIST/les_nml/ sst, psfc, shflx, lhflx, isrfc_type, ufric, is_dry_cbl, &
                    smag_constant, turb_prandtl, bflux, tran_coeff,   &
                    vert_scheme_type, avg_interval_sec, sampl_freq_sec,  &
                    expname, ldiag_les_out, km_min, min_sfc_wind, les_metric, &
                    max_turb_scale, wso_gradient, wso_mean, wso_max, wso_min, &
                    wso_spread, wso_fun, wso_cyc

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read Namelist for LES
  !!
  !! This subroutine 
  !! - reads the Namelist 
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state 
  !!
  !! @par Revision History
  !!  by Anurag Dipankar, MPIM (2013-04)
  !!
  SUBROUTINE read_les_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename 
    INTEGER :: istat, funit, jg
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_les_nml: read_les_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    sst          = 300._wp
    psfc         = -999._wp
    shflx        = 0.1_wp 
    lhflx        = 0._wp 
    isrfc_type   = 1 
    ufric        = -999._wp 

    is_dry_cbl   = .FALSE.

    !parameters
    smag_constant    = 0.23_wp
    turb_prandtl     = 0.33333333333_wp
    km_min           = 0.0_wp  
    max_turb_scale   = 300._wp
    min_sfc_wind     = 1._wp !Default from Holstag and Boville 1991

    bflux       = 0.0007_wp
    tran_coeff  = 0.02_wp

    vert_scheme_type = 2 !implicit

    !output parameters
    ldiag_les_out = .FALSE. 
    expname  = 'ICOLES'
    avg_interval_sec = 900._wp
    sampl_freq_sec   = 60._wp

    turb_profile_list = (/                                                     &
      'u      ','v      ','w      ','th     ','exner  ','rho    ','qv     ',   & !1-7
      'qc     ','wu     ','wv     ','wth    ','wqv    ','wqc    ','ww     ',   & !8-14
      'thth   ','qvqv   ','qcqc   ','uu     ','vv     ','kh     ','km     ',   & !15-21
      'thv    ','wthv   ','wqvd   ','wthd   ','wqcd   ','bruvais','mechprd',   & !22-28
      'wud    ','wvd    ','wthsfs ','rh     ','clc    ','qi     ','qs     ',   & !29-35
      'qr     ','qg     ','qh     ','lwf    ','swf    ','dt_t_sw','dt_t_lw',   & !36-42
      'dt_t_tb','dt_t_mc' /)    !43-44 

    turb_tseries_list = (/                                          &
      'ccover ','shflx  ','lhflx  ','ustress','vstress','tsfc   ',  & !1-6
      'qsfc   ','hbl    ','psfc   ','swf_tom','lwf_tom','swf_sfc',  & !7-12
      'lwf_sfc','precp_t','precp_r','precp_s','precp_g','precp_h',  & !13-18
      'precp_i' /)                                                    !19

    ! grid metric terms in the les diffusion
    les_metric = .FALSE.

    ! SBr, CHa
    wso_gradient = 0.0_wp
    wso_fun      = 4
    wso_cyc      = 1

    !wso_mean(1) = 1.8E-3_wp*2._wp
    !wso_mean(2) = 3.7E-3_wp*2._wp
    !wso_mean(3) = 11.3E-3_wp*2._wp
    !wso_mean(4) = 35.1E-3_wp*2._wp
    !wso_mean(5) = 56.7E-3_wp*2._wp
    !wso_mean(6) = 254.6E-3_wp*2._wp
    !wso_mean(7) = 763.9E-3_wp*2._wp
    !wso_mean(8) = 2291.5E-3_wp*2._wp

    wso_mean(1) = 1.2E-3_wp*1._wp
    wso_mean(2) = wso_mean(1)*2._wp
    wso_mean(3) = wso_mean(1)*6._wp
    wso_mean(4) = wso_mean(1)*18._wp
    wso_mean(5) = wso_mean(1)*54._wp
    wso_mean(6) = wso_mean(1)*162._wp
    wso_mean(7) = wso_mean(1)*486._wp
    wso_mean(8) = wso_mean(1)*1458._wp

    wso_max(1) = 3.40E-3_wp*1._wp
    wso_max(2) = wso_max(1)*2._wp
    wso_max(3) = wso_max(1)*6._wp
    wso_max(4) = wso_max(1)*18._wp
    wso_max(5) = wso_max(1)*54._wp
    wso_max(6) = wso_max(1)*162._wp
    wso_max(7) = wso_max(1)*486._wp
    wso_max(8) = wso_max(1)*1458._wp

    wso_min(1) = 0.35E-3_wp*1._wp
    wso_min(2) = wso_min(1)*2._wp
    wso_min(3) = wso_min(1)*6._wp
    wso_min(4) = wso_min(1)*18._wp
    wso_min(5) = wso_min(1)*54._wp
    wso_min(6) = wso_min(1)*162._wp
    wso_min(7) = wso_min(1)*486._wp
    wso_min(8) = wso_min(1)*1458._wp

    wso_mean(1) = 1.2E-3_wp*1._wp
    wso_mean(2) = wso_mean(1)*2._wp
    wso_mean(3) = wso_mean(1)*6._wp
    wso_mean(4) = wso_mean(1)*18._wp
    wso_mean(5) = wso_mean(1)*54._wp
    wso_mean(6) = wso_mean(1)*162._wp
    wso_mean(7) = wso_mean(1)*486._wp
    wso_mean(8) = wso_mean(1)*1458._wp

    wso_spread = 0.1_wp
    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('les_nml')
      READ(funit,NML=les_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('les_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, les_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, les_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, les_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    DO jg = 1 , max_dom
      les_config(jg)% sst          =  sst
      les_config(jg)% psfc         =  psfc
      les_config(jg)% shflx        =  shflx
      les_config(jg)% lhflx        =  lhflx
      les_config(jg)% isrfc_type   =  isrfc_type
      les_config(jg)% ufric        =  ufric
      les_config(jg)% is_dry_cbl   =  is_dry_cbl
      les_config(jg)% smag_constant     =  smag_constant
      les_config(jg)% turb_prandtl      =  turb_prandtl
      les_config(jg)% rturb_prandtl     =  1._wp/turb_prandtl
      les_config(jg)% bflux             =  bflux
      les_config(jg)% tran_coeff        =  tran_coeff
      les_config(jg)% vert_scheme_type  =  vert_scheme_type
      les_config(jg)% ldiag_les_out     =  ldiag_les_out
      les_config(jg)% expname           =  expname
      les_config(jg)% avg_interval_sec  =  avg_interval_sec
      les_config(jg)% sampl_freq_sec    =  sampl_freq_sec
      les_config(jg)% km_min            =  km_min
      les_config(jg)% max_turb_scale    =  max_turb_scale
      les_config(jg)% min_sfc_wind      =  min_sfc_wind
      les_config(jg)% les_metric        =  les_metric
      ! SBr, CHa
      les_config(jg)% wso_gradient      =  wso_gradient
      les_config(jg)% wso_mean          =  wso_mean
      les_config(jg)% wso_max           =  wso_max
      les_config(jg)% wso_min           =  wso_min
      les_config(jg)% wso_spread        =  wso_spread
      les_config(jg)% wso_fun           =  wso_fun
      les_config(jg)% wso_cyc           =  wso_cyc
    END DO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=les_nml)                    
      CALL store_and_close_namelist(funit,'les_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=les_nml)

  END SUBROUTINE read_les_namelist

END MODULE mo_les_nml
